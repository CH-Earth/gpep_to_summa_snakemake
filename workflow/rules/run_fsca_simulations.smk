"""

Snakemake file to calculate fSCA based on SUMMA model simulations.

"""

from pathlib import Path
import sys
# Import local packages
sys.path.append(str(Path('../').resolve()))
sys.path.append(str(Path('../../').resolve()))
sys.path.append(str(Path('../../snow_dist').resolve()))
from scripts import ss_utils
from snow_dist import summa_simulation

# Resolve all file paths and directories in the config file
config['summa_output_dir'] = '/Users/drc858/Data/gpep/RF_ens/summa/output'
config['working_dir'] = '/Users/drc858/Data/gpep/RF_ens/summa/output'
config['summa_interim_dir'] = '/Users/drc858/Data/gpep/RF_ens/summa/interim'
config['fsca_output_dir'] =  '/Users/drc858/Data/gpep/RF_ens/summa/output/fsca'
config_file_path = '/Users/drc858/GitHub/snow_dist/settings/config_summa_model_tuolumne.yaml'

#Read first raw forcing file for easymore remapping
summa_result_filepaths = list(Path(config['summa_output_dir']).glob('*.nc'))
summa_filenames = [Path(filepath).name for filepath in summa_result_filepaths]

fsca_simulation = summa_simulation.SummaSimulation(config_file_path)

# Prepare summa model results for fsca
# This includes summing the input and output snowpack variables
rule run_fsca_model:
    input:
        expand(Path(config['summa_interim_dir'], "{filenames}"))

rule prepare_fsca_model_input:
    input:
        summa_result = Path(config['summa_output_dir'],"{filenames}")
    output:
        fsca_input = Path(config['summa_interim_dir'],"{filenames}")
    run:
        summa_ds = xr.open_dataset(input.input_summa_result)
        fsca_ds = summa_fsca_model.prepare_summa_fsca_input(summa_ds)
        fsca_ds.to_netcdf(output.fsca_input)

# Run fsca settings
rule run_fsca_simulations:
    input:
        prepped_files = expand(Path(config['summa_interim_dir'], "{filenames}"), filenames=summa_filenames)
    output:
        fsca_result = Path(config['fsca_output_dir'], "{filenames}")
    run:
        sim = summa_simulation.SummaSimulation(config)
        summa_fsca_model.run_summa_fsca_model(input.prepped_files,output.fsca_result)

