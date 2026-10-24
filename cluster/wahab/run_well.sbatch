#!/bin/bash
#SBATCH --job-name=celltrack
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%A_%a.out
#
# One movie per array task on a 4-GPU Wahab node. Wells come from cluster/wahab/wells.txt
# (one per line, e.g. well_B3). Submit the plate, at most 4 nodes at a time:
#
#   cd /scratch/$USER/celltrack && mkdir -p logs
#   sbatch --array=0-$(( $(wc -l < cluster/wahab/wells.txt) - 1 ))%4 cluster/wahab/run_well.sbatch
#
# A single well: WELL=well_B3 sbatch cluster/wahab/run_well.sbatch
# Every stage is resumable; re-submitting a failed task continues where it stopped.
# On success the task writes outputs/<well>/DONE, which the Mac-side sync loop
# (cluster/wahab/sync_results.sh) watches.
set -uo pipefail

PROJ=${PROJ:-/scratch/$USER/celltrack}
cd "$PROJ"
if [ -z "${WELL:-}" ]; then
    WELL=$(sed -n "$(( ${SLURM_ARRAY_TASK_ID:-0} + 1 ))p" cluster/wahab/wells.txt)
fi
[ -n "$WELL" ] || { echo "no well for task ${SLURM_ARRAY_TASK_ID:-?}"; exit 1; }
OUT=outputs/$WELL
mkdir -p "$OUT" logs

module load python3/2026.0-py312          # only this container module (others shadow `crun`)
PY="crun -p $PROJ/env python"
export MPLCONFIGDIR=$PROJ/.cache/mpl NUMBA_CACHE_DIR=$PROJ/.cache/numba
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
NGPU=$(nvidia-smi -L | wc -l)
echo "$(date)  $WELL on $(hostname), $NGPU GPUs, $(nproc) cores, job ${SLURM_JOB_ID:-?}"
stage() { echo; echo "===== $(date)  $WELL  $*"; }
echo "$(date) running $(hostname) ${SLURM_JOB_ID:-}" > "$OUT/RUNNING"; rm -f "$OUT/FAILED"
fail() { echo "$(date) $*" > "$OUT/FAILED"; rm -f "$OUT/RUNNING"; echo "FAILED: $*"; exit 1; }

# 0. per-well config (paths + n_frames from the TIFF)
CFG=configs/wells/$WELL.yaml
$PY cluster/wahab/make_well_config.py "$WELL" || fail "config"
T=$($PY -c "import yaml;print(yaml.safe_load(open('$CFG'))['acquisition']['n_frames'])")
CHUNK=$($PY -c "import yaml;print(yaml.safe_load(open('$CFG'))['tracking']['chunk_frames'])")

# 1. ingest
if [ ! -d "$OUT/flat.zarr" ]; then
    stage ingest
    $PY -m celltrack.ingest --config "$CFG" || fail "ingest"
fi

# 2. segmentation: one Cellpose-SAM process per GPU on a contiguous frame range
NDONE=$($PY -c "import json,os;p='$OUT/masks_done.json';print(len(json.load(open(p))) if os.path.exists(p) else 0)")
if [ "$NDONE" -lt "$T" ]; then
    stage "segment ($NDONE / $T done) on $NGPU GPUs"
    PER=$(( (T + NGPU - 1) / NGPU ))
    for i in $(seq 0 $((NGPU - 1))); do
        A=$((i * PER)); B=$(( (i + 1) * PER )); [ $B -gt $T ] && B=$T
        CUDA_VISIBLE_DEVICES=$i OMP_NUM_THREADS=6 $PY -m celltrack.segment --config "$CFG" run --start $A --stop $B \
            > "logs/segment_${WELL}_gpu${i}.log" 2>&1 &
        sleep 20
    done
    wait
    NDONE=$($PY -c "import json;print(len(json.load(open('$OUT/masks_done.json'))))")
    [ "$NDONE" -lt "$T" ] && fail "segmentation incomplete ($NDONE / $T), see logs/segment_${WELL}_gpu*.log"
fi

# 3. features
stage features
$PY -m celltrack.segment --config "$CFG" features --workers 24 || fail "features"

# 4. HOCT forward (GPUs 0,1) and time-reversed (GPUs 2,3), several chunks in flight each.
#    The ILP (SCIP) memory grows with cell density: 12 concurrent solvers fit in 180 GB up to
#    ~1.2 M detections (~2000 cells/frame); denser wells (rows G-H: 1.7-1.8 M) ran out of
#    memory, so the concurrency is scaled by the detection count (override with HOCT_JOBS).
NDET=$($PY -c "import pyarrow.parquet as pq;print(pq.ParquetFile('$OUT/detections.parquet').metadata.num_rows)")
if [ -n "${HOCT_JOBS:-}" ]; then JOBS=$HOCT_JOBS
elif [ "$NDET" -gt 1500000 ]; then JOBS=2
elif [ "$NDET" -gt 1100000 ]; then JOBS=3
else JOBS=6; fi
stage "tracking ($NDET detections, $JOBS chunks in flight per direction)"
$PY -m celltrack.track --config "$CFG" hoct --stop $T --chunk $CHUNK --name hoct --jobs $JOBS --gpus 0,1 \
    > "logs/track_${WELL}_hoct.log" 2>&1 &
P_FWD=$!
$PY -m celltrack.track --config "$CFG" hoct --reverse --stop $T --chunk $CHUNK --name hoct --jobs $JOBS --gpus 2,3 \
    > "logs/track_${WELL}_hoct_rev.log" 2>&1 &
P_REV=$!
wait $P_FWD || fail "forward HOCT, see logs/track_${WELL}_hoct.log"
wait $P_REV || echo "WARNING: reversed HOCT failed (high-confidence matrices skipped)"

# 5. validation -> states -> transition -> motility -> movies (no comparison linkers for the plate)
stage analysis
export CUDA_VISIBLE_DEVICES=0
$PY -m celltrack.pipeline --config "$CFG" --skip features track compare --compare-frames 0 || fail "analysis"

rm -rf "$OUT"/tracks/*_work
echo "$(date) $(hostname) ${SLURM_JOB_ID:-}" > "$OUT/DONE"; rm -f "$OUT/RUNNING"
echo; echo "===== $(date)  $WELL done. Report: $OUT/analysis/transition_hoct/report.html"
