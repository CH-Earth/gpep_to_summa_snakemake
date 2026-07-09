#!/bin/bash
#SBATCH --job-name=download_modis
#SBATCH --output=download_modis_tuolumne.out
#SBATCH --error=download_modis_tuolumne.err
#SBATCH --time=18:00:00
#SBATCH --account=ees250064
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail
# Earthdata: for sbatch, use ~/.netrc for urs.earthdata.nasa.gov if the job is
# non-interactive. With credentials in place, long runs are usually CMR paging
# (many granules) or thousands of sequential downloads — not a "hang".
# Unbuffered Python so .out updates live under Slurm.
export PYTHONUNBUFFERED=1

# Load required modules
module use "$HOME/privatemodules"
module load conda-env/gpep_to_summa_snakemake-py3.12.8

# Initialize conda for non-login shells
source /apps/anvil/external/apps/conda/2025.02/etc/profile.d/conda.sh
conda activate gpep_to_summa_snakemake

PYTHON_SCRIPT="/home/x-dcasson/GitRepos/gpep_to_summa_snakemake/workflow/scripts/download_modis.py"
OUTPUT_DIR="/anvil/projects/x-ees240082/users/dcasson/modis/bow/"
# Tuolumne bounding box
#BOUNDING_BOX="-121,37,-118,39"
#Chena bounding box
#BOUNDING_BOX="-147,64.5,-144,65.5"
#Bow bounding box
BOUNDING_BOX="-119,49,-111,53"

mkdir -p "$OUTPUT_DIR"

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: download_aso.py not found at ${PYTHON_SCRIPT}" >&2
    exit 1
fi

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Error: output directory not found at ${OUTPUT_DIR}" >&2
    exit 1
fi

cd "$OUTPUT_DIR"
BOUNDING_BOX="$BOUNDING_BOX" python -u "$PYTHON_SCRIPT" "$@"
