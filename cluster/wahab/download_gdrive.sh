#!/bin/bash
# Download the shared Google Drive folder (48 well stacks + plate-layout spreadsheet) to scratch.
#
#   bash cluster/wahab/download_gdrive.sh [FOLDER_URL]
#
# Run on the login node (compute nodes have no internet). Resumable: gdown skips files that
# already exist with the right size. Requires `pip install gdown` in the project env.
set -euo pipefail

FOLDER=${1:-"https://drive.google.com/drive/folders/1pCGU_lYkf4yXD0daXI6n6f7ZUcpaw0R5"}
PROJ=${PROJ:-/scratch/$USER/celltrack}

module load python3/2026.0-py312
cd "$PROJ"
mkdir -p data/gdrive
crun -p "$PROJ/env" python -m gdown --folder --continue "$FOLDER" -O data/gdrive
