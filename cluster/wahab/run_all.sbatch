#!/bin/bash
#SBATCH --job-name=celltrack
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=logs/celltrack_%j.out
#
# Full pipeline on one Wahab GPU node (4x V100-16GB, 32 cores).
#   cd /scratch/$USER/celltrack && mkdir -p logs && sbatch cluster/wahab/run_all.sbatch
# Every stage is resumable: re-submitting the job continues where it stopped.
set -uo pipefail

PROJ=${PROJ:-/scratch/$USER/celltrack}
CFG=$PROJ/configs/wahab.yaml
cd "$PROJ"

# only this module: other container modules (e.g. ffmpeg) put their own `crun` first on PATH.
# ffmpeg for the overlay movies comes bundled with the imageio-ffmpeg wheel.
module load python3/2026.0-py312
PY="crun -p $PROJ/env python"

export MPLCONFIGDIR=$PROJ/.cache/mpl NUMBA_CACHE_DIR=$PROJ/.cache/numba
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" logs
NGPU=$(nvidia-smi -L | wc -l)
echo "$(date)  node $(hostname), $NGPU GPUs, $(nproc) cores"
nvidia-smi -L

stage() { echo; echo "===== $(date)  $*"; }

# 1. ingest: TIFF -> zarr, background flattening, bubble mask (CPU)
if [ ! -d outputs/flat.zarr ]; then
    stage ingest
    $PY -m celltrack.ingest --config "$CFG" || exit 1
fi

# 2. segmentation: one Cellpose-SAM process per GPU, each on a contiguous frame range.
#    Workers share masks.zarr / masks_done.json (disjoint ranges; concurrent-safe).
T=$($PY -c "import yaml;print(yaml.safe_load(open('$CFG'))['acquisition']['n_frames'])")
NDONE=$($PY -c "import json,os;p='outputs/masks_done.json';print(len(json.load(open(p))) if os.path.exists(p) else 0)")
if [ "$NDONE" -lt "$T" ]; then
    stage "segment ($NDONE / $T frames done) on $NGPU GPUs"
    PER=$(( (T + NGPU - 1) / NGPU ))
    for i in $(seq 0 $((NGPU - 1))); do
        A=$((i * PER)); B=$(( (i + 1) * PER )); [ $B -gt $T ] && B=$T
        CUDA_VISIBLE_DEVICES=$i $PY -m celltrack.segment --config "$CFG" run --start $A --stop $B \
            > logs/segment_gpu${i}_${SLURM_JOB_ID}.log 2>&1 &
        sleep 20   # let the first worker create masks.zarr before the others open it
    done
    wait
    NDONE=$($PY -c "import json;print(len(json.load(open('outputs/masks_done.json'))))")
    echo "segmented frames: $NDONE / $T"
    [ "$NDONE" -lt "$T" ] && { echo "segmentation incomplete, see logs/segment_gpu*"; exit 1; }
fi

# 3. features (mask cleaning + per-object measurements; CPU process pool)
stage features
$PY -m celltrack.segment --config "$CFG" features --workers 24 || exit 1

# 4. tracking: HOCT forward (GPU0) and time-reversed (GPU1) over the whole movie, plus the
#    Trackastra / LapTrack comparison window (GPU2); all three concurrently.
stage tracking
CHUNK=$($PY -c "import yaml;print(yaml.safe_load(open('$CFG'))['tracking']['chunk_frames'])")
CUDA_VISIBLE_DEVICES=0 $PY -m celltrack.track --config "$CFG" hoct --stop $T --chunk $CHUNK --name hoct \
    > logs/track_hoct_${SLURM_JOB_ID}.log 2>&1 &
P_FWD=$!
CUDA_VISIBLE_DEVICES=$(( NGPU > 1 ? 1 : 0 )) $PY -m celltrack.track --config "$CFG" hoct --reverse --stop $T --chunk $CHUNK --name hoct \
    > logs/track_hoct_rev_${SLURM_JOB_ID}.log 2>&1 &
P_REV=$!
(
    export CUDA_VISIBLE_DEVICES=$(( NGPU > 2 ? 2 : 0 ))
    CS=113
    $PY -m celltrack.track --config "$CFG" trackastra --mode ilp_nodiv --stop $CS --chunk 12 --name ta_nodiv
    $PY -m celltrack.track --config "$CFG" trackastra --mode greedy --stop $CS --name ta_greedy
    $PY -m celltrack.track --config "$CFG" laptrack --stop $CS --name laptrack
) > logs/track_compare_${SLURM_JOB_ID}.log 2>&1 &
P_CMP=$!
wait $P_FWD || { echo "forward HOCT failed, see logs/track_hoct_${SLURM_JOB_ID}.log"; exit 1; }
wait $P_REV || echo "WARNING: reversed HOCT failed (high-confidence matrices will be skipped)"
wait $P_CMP || echo "WARNING: comparison linkers failed (linker comparison will be incomplete)"

# 5. validation -> states -> transition matrices -> motility -> overlay movies
stage "analysis"
export CUDA_VISIBLE_DEVICES=0
$PY -m celltrack.pipeline --config "$CFG" --skip features track compare
echo
echo "===== $(date)  done. Report: outputs/analysis/transition_hoct/report.html"
