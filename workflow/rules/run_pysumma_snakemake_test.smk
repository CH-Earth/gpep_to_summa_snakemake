"""

Snakemake file to run the base SUMMA simulations.

The model simulation is chunked by GRU to allow for parallelization on a cluster.

The chunks of GRUs are defined by the user

"""

from pathlib import Path
import sys
import pysumma as ps
# Import local packages
sys.path.append(str(Path('../').resolve()))

sys.path.append('../')
from scripts import gpep_to_summa_utils as gts_utils

# UPDATE LOCAL SUMMA PATH
config['summa_exe'] = '/Users/dcasson/GitHub/summa/bin/summa.exe'
config['summa_forcing_dir'] = Path('/Users/dcasson/Data/gpep/chena/forcing/summa/')

# Resolve all file paths and directories in the config file
config['file_manager'] = '/Users/dcasson/Data/yukon_esp/summa/settings/fileManagerSWE.txt'
config['summa_output_dir'] = '/Users/dcasson/Data/yukon_esp/summa_output/'
config['attributes_nc'] = '/Users/dcasson/Data/yukon_esp/summa/settings/attributes.nc'
config['case_name'] = 'chena'
config['run_suffix'] = 'base'


ens_set, _ = gts_utils.build_ensemble_list(config['summa_forcing_dir'])
file_paths = gts_utils.list_files_in_subdirectory(config['summa_forcing_dir'])
ens = list(ens_set)

#Filter files remove those that end in .txt
file_paths = [file for file in file_paths if not file.endswith('.txt')]

forcing_file_list = []
for ens_member in ens:
    # write a new file containing the forcing files which correspond to the ensemble member
    forcing_list_file = f'{config["summa_forcing_dir"]}/forcing_files_{ens_member}.txt'
    with open(forcing_list_file, 'w') as f:
        for file in file_paths:
            #if ensemble member in first three characters of file name
            file_path = Path(file)
            if str(file_path.parent) == ens_member:
                f.write(f'{file}.nc\n')
    forcing_file_list.append(forcing_list_file)

rule run_summa_base_simulations:
    input:
        expand(Path(config['summa_output_dir'],f"{config['case_name']}_{{ens_member}}_timestep.nc"),ens_member=ens)

rule run_summa_ensemble_simulations:
    input:
        file_manager = Path(config['file_manager']),
        forcing_file = lambda wildcards: forcing_file_list[ens.index(wildcards.ens_member)]
    output:
        summa_chunked_output = Path(config['summa_output_dir'],f"{config['case_name']}_{{ens_member}}_timestep.nc")
    params:
        summa_exe = config['summa_exe'],
        run_suffix = lambda wildcards: wildcards.ens_member,
    run:
        print(input.forcing_file)
        sim = ps.Simulation(params.summa_exe, input.file_manager)
        sim._force_file_list = Path(input.forcing_file)
        print(sim.force_file_list)
        print(str(params.run_suffix))
        sim.run(run_suffix=str(params.run_suffix))


        
