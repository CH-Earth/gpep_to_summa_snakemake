#!/bin/bash
#SBATCH --job-name=download_aso
#SBATCH --output=download_aso_new.out
#SBATCH --error=download_aso_new.err
#SBATCH --time=04:00:00
#SBATCH --account=ees250064
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

# Load required modules
module use "$HOME/privatemodules"
module load conda-env/gpep_to_summa_snakemake-py3.12.8

# Initialize conda for non-login shells
source /apps/anvil/external/apps/conda/2025.02/etc/profile.d/conda.sh
conda activate gpep_to_summa_snakemake

PYTHON_SCRIPT="/home/x-dcasson/GitRepos/gpep_to_summa_snakemake/workflow/scripts/download_aso_tuolumne.py"
OUTPUT_DIR="/anvil/projects/x-ees240082/users/dcasson/lidar/tuolumne/"

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
python "$PYTHON_SCRIPT" "$@"
