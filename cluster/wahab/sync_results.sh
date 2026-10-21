#!/bin/bash
# Mac-side loop: pull every finished well from Wahab as soon as its DONE marker appears, and
# keep outputs/wells/progress.md up to date.
#
#   caffeinate -i cluster/wahab/sync_results.sh            # foreground
#   nohup caffeinate -i cluster/wahab/sync_results.sh > outputs/wells/sync.log 2>&1 &
#
# Needs the ssh control-master session (`ssh -fN wahab`). Pulls flat.zarr + masks + detections
# + tracks + analysis + qc (~4.5 GB per well; raw.zarr and the chunk cache stay on the cluster),
# so each well can be opened locally with
#   python -m celltrack.view napari --config outputs/wells/<well>/config.yaml --name hoct
set -uo pipefail
cd "$(dirname "$0")/../.."
REMOTE=wahab
RPROJ=/scratch/dwadduwa/celltrack
LOCAL=outputs/wells
INTERVAL=${INTERVAL:-600}
mkdir -p "$LOCAL"

pull_well() {
    local w=$1
    echo "$(date '+%F %T')  syncing $w"
    rsync -a --partial \
        --include='/analysis/***' --include='/qc/***' --include='/masks.zarr/***' --include='/flat.zarr/***' \
        --include='/tracks/' --include='/tracks/*.parquet' --include='/tracks/*.json' --include='/tracks/*.geff/***' \
        --include='/*.parquet' --include='/*.json' --include='/*.npy' --include='/DONE' \
        --exclude='*' \
        "$REMOTE:$RPROJ/outputs/$w/" "$LOCAL/$w/" || return 1
    # local config so the viewer / re-analysis work from the Mac
    sed -e "s#^  raw_tif: .*#  raw_tif: /dev/null   \# raw TIFF not synced#" \
        -e "s#^\(  [a-z_]*: \)\.\./outputs#\1.#" -e 's/device: cuda/device: mps/' configs/wahab.yaml \
        | sed '1s#.*#\# Local copy of the Wahab run of '"$w"' (synced by cluster/wahab/sync_results.sh)#' > "$LOCAL/$w/config.yaml"
    date '+%F %T' > "$LOCAL/$w/SYNCED"
    echo "$(date '+%F %T')  $w synced ($(du -sh "$LOCAL/$w" | cut -f1))"
}

write_progress() {
    local remote_state=$1
    {
        echo "# Plate progress ($(date '+%F %T'))"
        echo
        echo "| well | cluster | job | synced to Mac |"
        echo "|---|---|---|---|"
        while read -r w state info; do
            local synced="—"
            [ -f "$LOCAL/$w/SYNCED" ] && synced=$(cat "$LOCAL/$w/SYNCED")
            echo "| $w | $state | $info | $synced |"
        done <<< "$remote_state"
        echo
        echo "Queue:"; echo '```'; ssh -o BatchMode=yes $REMOTE 'squeue -u $USER -o "%.12i %.3t %.10M %R" | grep -E "JOBID|celltrack"' 2>/dev/null; echo '```'
    } > "$LOCAL/progress.md"
}

while true; do
    state=$(ssh -o BatchMode=yes $REMOTE bash -s <<'EOF' 2>/dev/null
cd /scratch/$USER/celltrack
for w in $(cat cluster/wahab/wells.txt) well_A1; do
  d=outputs/$w
  if   [ -f $d/DONE ];    then echo "$w done $(cut -d' ' -f1-2 $d/DONE | tr ' ' 'T')"
  elif [ -f $d/FAILED ];  then echo "$w FAILED $(cut -c1-60 $d/FAILED | tr ' ' '_')"
  elif [ -f $d/RUNNING ]; then
       n=$(python3 -c "import json,os;p='$d/masks_done.json';print(len(json.load(open(p))) if os.path.exists(p) else 0)" 2>/dev/null)
       st=$(grep "=====" logs/celltrack_*_*.out 2>/dev/null | grep " $w " | tail -1 | sed 's/.*  '"$w"'  //' | cut -c1-30)
       echo "$w running ${st:-starting}_(seg_${n}_frames)"
  else echo "$w queued -"; fi
done
EOF
)
    [ -z "$state" ] && { echo "$(date '+%F %T')  ssh to $REMOTE failed (session expired? run: ssh -fN wahab)"; sleep "$INTERVAL"; continue; }
    while read -r w st _; do
        if [ "$st" = "done" ] && [ ! -f "$LOCAL/$w/SYNCED" ]; then pull_well "$w" || echo "$(date '+%F %T')  rsync of $w failed, will retry"; fi
    done <<< "$state"
    write_progress "$state"
    n_done=$(grep -c ' done ' <<< "$state"); n_syn=$(ls "$LOCAL"/*/SYNCED 2>/dev/null | wc -l | tr -d ' ')
    echo "$(date '+%F %T')  $n_done done on cluster, $n_syn synced, $(grep -c ' running ' <<< "$state") running, $(grep -c ' FAILED ' <<< "$state") failed"
    [ "$n_syn" -ge "$(( $(wc -l < cluster/wahab/wells.txt) + 1 ))" ] && { echo "all wells synced"; break; }
    sleep "$INTERVAL"
done
