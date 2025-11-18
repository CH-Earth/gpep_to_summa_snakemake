''' 
gmet to summa snakemake master snakemake file

This snakemake file runs all the steps required to convert GPEP forcings to SUMMA forcings.

Original process code: Andy Wood
Adapted to Snakemake: Dave Casson

'''
from pathlib import Path
import sys
#Update the current working directory to the directory of the file
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))
print(f"Scripts path: {scripts_path}")
import gpep_to_summa_utils as utils

# Resolve paths from the configuration file
config = utils.resolve_paths(config, log_config = True)

# Read in all local snakemake files and rules
include: './gpep_file_prep.smk'
include: './remap_gpep_to_shp.smk'
include: './metsim_file_prep.smk'
include: './run_metsim.smk'
include: './metsim_to_summa.smk'

ruleorder: 
    gpep_to_summa >
    gpep_file_prep >
    remap_gpep_to_shp >
    metsim_file_prep >
    run_metsim >
    metsim_to_summa
    
# Read all forcing files and create a list based on the output directory (i.e. ens/filename.nc)
_, file_path_list = utils.build_ensemble_list(config['gpep_forcing_dir'])

# Run the snakemake file, so that that it produces a summa input file for each of the gpep forcing files
rule gpep_to_summa:
    input:
        expand(Path(config['summa_forcing_dir'],'{forcing_file}.nc'), forcing_file = file_path_list)
        #Path(config['summa_forcing_dir'], "summa_mean_forcing.nc")
