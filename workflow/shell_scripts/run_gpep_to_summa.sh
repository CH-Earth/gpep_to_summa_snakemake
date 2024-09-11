#!/bin/bash

module load StdEnv/2020
module load mpi4py netcdf nco libspatialindex

source $HOME/summa_snakemake/bin/activate


# Directory for logs
log_dir="logs"

# Check if log directory exists, if not create it
if [ ! -d "$log_dir" ]; then
    mkdir -p "$log_dir"
fi

case_names=(
    "rf_dynamic_boxcox_high_predictors"
    # Add additional case names here, each on a new line
)

# Define the result directory
result_dir='/scratch/dcasson/gpep/ensemble/bow/'

# Loop over each case name and run Snakemake
for case_name in "${case_names[@]}"; do
    # Create a log file name based on the configuration file
    log_file="${log_dir}/${case_name%.yaml}_log.txt"
    mkdir -p $(dirname $log_file)
    
    snakemake -s ../gpep_to_summa_multiple.smk -c 16 --config case_name="$case_name" result_dir="$result_dir" --profile summa_profile --unlock &> /dev/null
    snakemake -s ../gpep_to_summa_multiple.smk -c 16 --config case_name="$case_name" result_dir="$result_dir" --profile summa_profile 2>&1 | tee -a $log_file &
done