#!/bin/bash
#SBATCH --job-name=process_modis
#SBATCH --output=process_modis_tuolumne.out
#SBATCH --error=process_modis_tuolumne.err
#SBATCH --time=01:00:00
#SBATCH --account=ees250064
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail
export PYTHONUNBUFFERED=1

# Conda stack (match run_download_modis.sh)
module use "$HOME/privatemodules"
module load conda-env/gpep_to_summa_snakemake-py3.12.8
source /apps/anvil/external/apps/conda/2025.02/etc/profile.d/conda.sh
conda activate gpep_to_summa_snakemake

SCRIPT="/home/x-dcasson/GitRepos/gpep_to_summa_snakemake/workflow/scripts/process_modis_hru.py"
MODIS_DIR="/anvil/projects/x-ees240082/users/dcasson/modis/chena/"
HRUS="/anvil/projects/x-ees240082/users/dcasson/gpep/chena/gis/chena_tdx.gpkg"
OUT_NC="/anvil/projects/x-ees240082/users/dcasson/modis/chena_mod10a1_ndsi_snow_cover.nc"

# Chena / Bow: edit paths and sbatch from submit directory if log names should differ

if [ ! -f "$SCRIPT" ]; then
    echo "Error: process_modis_hru.py not found at ${SCRIPT}" >&2
    exit 1
fi

if [ ! -d "$MODIS_DIR" ]; then
    echo "Error: MODIS directory not found: ${MODIS_DIR}" >&2
    exit 1
fi

if [ ! -f "$HRUS" ]; then
    echo "Error: HRU file not found: ${HRUS}" >&2
    exit 1
fi

python -u "$SCRIPT" \
    --modis-dir "$MODIS_DIR" \
    --hrus "$HRUS" \
    --output "$OUT_NC" \
    "$@"
