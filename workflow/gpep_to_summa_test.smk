''' 
gmet to summa snakemake master snakemake file

This snakemake file runs all the steps required to convert GPEP forcings to SUMMA forcings.

Original process code: Andy Wood
Adapted to Snakemake: Dave Casson

'''
from pathlib import Path
from scripts import gpep_to_summa_utils as utils

# Resolve paths from the configuration file
config = utils.resolve_paths(config, log_config = True)

# Read in all local snakemake files and rules
include: './rules/gpep_file_prep.smk'
include: './rules/remap_gpep_to_shp.smk'
include: './rules/metsim_file_prep.smk'
include: './rules/run_metsim.smk'
include: './rules/metsim_to_summa.smk'

ruleorder: 
    gpep_to_summa >
    gpep_file_prep >
    remap_gpep_to_shp >
    metsim_file_prep >
    run_metsim >
    metsim_to_summa
    
case_name = 'lwr_static_gamma_low_predictors'
result_dir = '/scratch/dcasson/gpep/transform_tests/'

config['case_name'] = case_name
config['gpep_forcing_dir'] =  Path(result_dir, case_name, 'ensembles')
config['working_dir'] = Path(result_dir, case_name, 'summa_forcing')
config['summa_dir'] = Path(result_dir, case_name, 'summa_forcing',summa)
config['summa_forcing_dir'] = Path(result_dir, case_name, 'summa_forcing',summa)

# Read all forcing files and create a list based on the output directory (i.e. ens/filename.nc)
_, file_path_list = utils.build_ensemble_list(config['gpep_forcing_dir'])

# Run the snakemake file, so that that it produces a summa input file for each of the gpep forcing files
rule gpep_to_summa:
    input:
        expand(Path(config['summa_forcing_dir'],'{forcing_file}.nc'), forcing_file = file_path_list)
