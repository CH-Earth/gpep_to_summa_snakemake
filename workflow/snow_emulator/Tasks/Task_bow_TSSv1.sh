#!/bin/bash
########## SLURM Parameters - ANVIL #############
#SBATCH -A ees250064
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=61
#SBATCH --mem-per-cpu=1800
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH -J bow_TSSv1
#SBATCH -o bow_TSSv1_%j.out
#SBATCH -e bow_TSSv1_%j.err
#SBATCH --mail-user=nicolas.vasquez@ucalgary.ca
#SBATCH --mail-type=ALL

################# Load Modules #########################
module --force purge
module load r/4.4.1
module load netcdf-c/4.7.4
module load netcdf-fortran/4.5.3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

################# Run script ###########################
dirScript="/anvil/projects/x-ees240082/users/nicolas-vasquez/SnowStations/Rfiles"
dirYamlFile="/anvil/projects/x-ees240082/users/nicolas-vasquez/SnowStations/Rfiles"
yamlFileName="config_LSE.yaml" 
nCores=60
basin="bow"            # Options: "bow", "chena", "tuolumne"
targetMetric="TSSv1"     # Options: "rmse", "spearman", "TSSv1"

cd ${dirScript}

Rscript 02_ObtainParams.R ${dirYamlFile} ${yamlFileName} ${nCores} ${basin} ${targetMetric}
