#!/bin/bash
#SBATCH --job-name=process_mixed_precip
#SBATCH --output=process_mixed_precip_tuolumne.out
#SBATCH --error=process_mixed_precip_tuolumne.err
#SBATCH --time=03:00:00
#SBATCH --account=ees250064
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB

set -euo pipefail
export PYTHONUNBUFFERED=1

# Conda stack (match run_process_modis.sh)
module use "$HOME/privatemodules"
module load conda-env/gpep_to_summa_snakemake-py3.12.8
source /apps/anvil/external/apps/conda/2025.02/etc/profile.d/conda.sh
conda activate gpep_to_summa_snakemake

SCRIPT="/home/x-dcasson/GitRepos/gpep_to_summa_snakemake/workflow/scripts/process_mixed_precip_netcdf.py"
GPEP_ROOT="/anvil/projects/x-ees240082/users/dcasson/gpep"
CATCHMENT="bow"

# Edit CATCHMENT, --sources, and sbatch log names before submitting.
# Examples: --sources rf   or   --sources rf,casr,era5
# Add --overwrite to replace existing NetCDFs.

if [ ! -f "$SCRIPT" ]; then
    echo "Error: process_mixed_precip_netcdf.py not found at ${SCRIPT}" >&2
    exit 1
fi

python -u "$SCRIPT" \
    --catchment "$CATCHMENT" \
    --gpep-root "$GPEP_ROOT" \
    --sources casr,era5 \
    --time-join inner \
    --compress \
    "$@"
