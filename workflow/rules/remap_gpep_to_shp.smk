
# Import needed packages
from pathlib import Path
import sys
sys.path.append('../')
from scripts import remap_forcing_to_shp
from scripts import gpep_to_summa_utils as gts_utils

# Resolve paths from the configuration file
config = gts_utils.resolve_paths(config)

# Read all forcing files and create a list based on the output directory (i.e. ens/filename.nc)
ensemble_list, file_path_list = gts_utils.build_ensemble_list(config['gpep_forcing_dir'])

rule remap_gpep_to_shp:
    input:
        expand(Path(config['easymore_output_dir'],"{file}.nc"), file=file_path_list)

# Define rule to run file remapping when remap file exists
rule create_remap_file:
    input:
        input_forcing_files = expand(Path(config['gpep_tmp_forcing_dir'],"{forcing_file}.nc"), forcing_file=file_path_list),
        input_shp = config['catchment_shp']
    output:
        remap_nc = config['remap_file']
    resources:
        runtime= 2,
        mem_mb= 500
    run:
        remap_forcing_to_shp.remap_with_easymore(config, input.input_forcing_files[0],input.input_shp, output.remap_nc, only_create_remap_nc=True)

# Define rule to run file remapping when remap file exists
rule remap_with_easymore:
    input:
        input_forcing = Path(config['gpep_tmp_forcing_dir'],"{id}.nc"),
        input_shp = config['catchment_shp'],
        remap_nc = config['remap_file'],
    output:
        output_forcing = Path(config['easymore_output_dir'],"{id}.nc")
    params:
        file_path = "{id}"
    resources:
        runtime= 120,
        mem_mb= 3000
    run:
        try:
            remap_forcing_to_shp.remap_with_easymore(config, input.input_forcing, input.input_shp, input.remap_nc, only_create_remap_nc=False, file_path=params.file_path)
        except Exception as e:
            # Log the error or take some other action
            print(f"Error occurred: {e}")
            # Optionally, raise the exception to stop the workflow
            raise e
