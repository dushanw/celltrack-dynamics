#!/bin/bash
# One-time environment setup on the ODU Wahab cluster (run on the login node).
#
#   bash cluster/wahab/setup_env.sh
#
# Wahab ships software as Apptainer containers launched through `crun`; `crun -c -p DIR`
# creates a persistent venv stacked on the container's Python 3.12, and every later
# `crun -p DIR python ...` runs inside it. Compute nodes have no internet access, so the
# model weights (Cellpose-SAM, Trackastra, HOCT) are also downloaded here.
set -euo pipefail

PROJ=${PROJ:-/scratch/$USER/celltrack}
ENV=$PROJ/env
mkdir -p "$PROJ"
cd "$PROJ"

module load python3/2026.0-py312

if [ ! -d "$ENV" ]; then
    crun -c -p "$ENV"
fi
PY="crun -p $ENV python"

$PY -m pip install --upgrade pip wheel

# --- torch: V100 GPUs are compute capability 7.0 (sm_70). Recent CUDA-12.8+ wheels dropped
#     Volta, so fall back to the cu126 wheel index if the default wheel lacks sm_70.
$PY -m pip install torch torchvision
if ! $PY -c "import torch,sys; sys.exit(0 if 'sm_70' in torch.cuda.get_arch_list() else 1)"; then
    echo "default torch wheel has no sm_70 support -> installing cu126 wheels"
    $PY -m pip install --force-reinstall --index-url https://download.pytorch.org/whl/cu126 torch torchvision
    $PY -c "import torch; assert 'sm_70' in torch.cuda.get_arch_list(), torch.cuda.get_arch_list()"
fi

# --- the pipeline's dependencies (same packages as the local .venv, minus napari / Qt)
$PY -m pip install \
    "cellpose>=4.2" "trackastra>=0.5" "hoct>=0.3" laptrack motile traccuracy geff \
    "zarr>=3" dask scikit-image scikit-learn numpy pandas pyarrow tifffile \
    imageio imageio-ffmpeg networkx scipy pyyaml matplotlib numba tqdm

$PY -c "import torch, cellpose, trackastra, hoct, laptrack, zarr; print('torch', torch.__version__, torch.version.cuda, torch.cuda.get_arch_list())"

# --- model weights (cached under $HOME, which compute nodes can read)
$PY -c "from cellpose import models; models.CellposeModel(gpu=False, pretrained_model='cpsam_v2'); print('cellpose ok')"
$PY -c "from trackastra.model import Trackastra; Trackastra.from_pretrained('general_2d', device='cpu'); print('trackastra ok')"
$PY -c "from hoct._models import resolve_model; print('hoct', resolve_model())"

echo "environment ready: $ENV"
