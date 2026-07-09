#!/bin/bash

module --force purge
ml load gcc/14.2.0
ml load openblas/0.3.17
ml load netcdf-fortran/4.5.3
ml load  openmpi/5.0.5

# Load required modules
module use "$HOME/privatemodules"
module load conda-env/gpep_to_summa_snakemake-py3.12.8

# Initialize conda for non-login shells
source /apps/anvil/external/apps/conda/2025.02/etc/profile.d/conda.sh
conda activate gpep_to_summa_snakemake


# List of configuration files
config_files=(
    "../config/gpep_to_summa_bow.yaml"
    #"../config/gpep_to_summa_chena.yaml"
    #"../config/gpep_to_summa_tuolumne.yaml"
)

# Directory for logs
log_dir="logs"

# Where to put individual working directories (adjust as needed)
base_workdir="../workdir"

# Check if log directory exists, if not create it
if [ ! -d "$log_dir" ]; then
    mkdir -p "$log_dir"
fi

# Loop over each configuration file
for config_file in "${config_files[@]}"
do
    # Create a log file name based on the configuration file
    log_file="${log_dir}/${config_file%.yaml}_log.txt"
    mkdir -p "$(dirname "$log_file")"

    # A simple per‐config working directory based on the file name (minus .yaml)
    # e.g. ../workdir/run_gpep_rf_gamma_predictor_screening
    work_dir="${base_workdir}/$(basename "${config_file%.yaml}")"
    mkdir -p "$work_dir"

    echo "Running Snakemake for configuration file: $config_file" | tee -a "$log_file"

    # First, unlock (in case a previous run left a lock)
    snakemake \
        -s ../rules/run_pysumma_single_all_gru.smk \
        --configfile "$config_file" \
        -c 6 \
        --profile snakemake_default \
        --unlock \
        --directory "$work_dir" \
        &> /dev/null

    # Then run Snakemake in the background, using the separate working directory
    snakemake \
        -s ../rules/run_pysumma_single_all_gru.smk \
        --configfile "$config_file" \
        -c 6 \
        --profile snakemake_default \
        --directory "$work_dir" \
        2>&1 | tee -a "$log_file" &

done

# Wait for all background Snakemake processes
wait
echo "All workflows completed."