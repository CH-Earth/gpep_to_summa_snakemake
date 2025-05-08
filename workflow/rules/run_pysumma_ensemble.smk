"""

Snakemake file to run the base SUMMA simulations.

The model simulation is chunked by GRU to allow for parallelization on a cluster.

The chunks of GRUs are defined by the user

"""

from pathlib import Path
import sys
import pysumma as ps
import xarray as xr
# Import local packages
sys.path.append(str(Path('../').resolve()))

# UPDATE LOCAL SUMMA PATH
config['summa_exe'] = '/Users/dcasson/GitHub/summa/bin/summa.exe'
config['summa_forcing_dir'] = Path('/Users/dcasson/Data/gpep/chena/forcing/summa/')

# Resolve all file paths and directories in the config file
config['file_manager'] = '/Users/dcasson/Data/yukon_esp/summa/settings/fileManagerSWE.txt'
config['summa_output_dir'] = '/Users/dcasson/Data/yukon_esp/summa_output/'
config['case_name'] = 'chena'
config['run_suffix'] = 'base'

def build_ensemble_list(directory):
    ''' Build a list of the ensemble name and the file name for each file in the directory
        e.g. for each file in the directory: ens_forc.tuolumne.01d.2020.001.nc' --> 001/ens_forc.tuolumne.01d.2020.001.nc
    '''
    
    path = Path(directory)
    files = list(path.rglob('*.nc'))
    file_path_list = set()
    ensemble_list = set()  

    for file in files:
        if file.exists():
            filename = file.stem  # Get the filename without the extension
            ens = filename[-3:]  # Get the last three characters
            ensemble_list.add(ens)
            file_path_list.add(Path(ens, filename)) # Create and add the ens/filename.nc path

    return ensemble_list, file_path_list

def list_files_in_subdirectory(directory, suffix_to_remove='.nc'):
    path = Path(directory)
    file_paths = [file.as_posix().replace(suffix_to_remove, "") for file in path.glob('**/*nc') if file.is_file()]

    return file_paths

ens, _ = build_ensemble_list(config['summa_forcing_dir'])
ens = list(ens)
file_paths = list_files_in_subdirectory(config['summa_forcing_dir'],'')

rule run_summa_base_simulations:
    input:
        expand(Path(config['summa_output_dir'],f"{config['case_name']}_{{ens_member}}_timestep.nc"),ens_member=ens)

rule run_summa_ensemble_simulations:
    input:
        file_manager = Path(config['file_manager']),
        forcing_file = lambda wildcards: file_paths[ens.index(wildcards.ens_member)]
    output:
        summa_chunked_output = Path(config['summa_output_dir'],f"{config['case_name']}_{{ens_member}}_timestep.nc")
    params:
        summa_exe = config['summa_exe'],
        run_suffix = lambda wildcards: wildcards.ens_member,
    run:
        sim = ps.Simulation(params.summa_exe, input.file_manager)
        forcing_ds = xr.open_dataset(input.forcing_file)
        sim.assign_forcing_file(params.run_suffix, forcing_ds)
        sim.run(run_suffix=str(params.run_suffix))


        
