# Import needed packages
from pathlib import Path
import sys
#Update the current working directory to the directory of the file
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))
print(f"Scripts path: {scripts_path}")
import remap_forcing_to_shp
import gpep_to_summa_utils as gts_utils

# Resolve paths from the configuration file
config = gts_utils.resolve_paths(config)

# Read all forcing files and create a list based on the output directory (i.e. ens/filename.nc)
ensemble_list, file_path_list = gts_utils.build_ensemble_list(config['gpep_forcing_dir'])
print(f"File path list: {file_path_list}")
# Convert the set to a sorted list so you can index into it
file_path_list = sorted(file_path_list)

# Now you can grab the "first" entry
first_file = (Path(config['gpep_tmp_forcing_dir']) / file_path_list[0]).with_suffix(".nc")

print(f"First file: {first_file}")
rule remap_gpep_to_shp:
    input:
        expand(Path(config['gpep_tmp_forcing_dir'],"{file}.nc"), file=file_path_list)

# Define rule to run file remapping when remap file exists
rule create_remap_file:
    input:
        input_file = first_file,
        input_shp = config['catchment_shp']
    output:
        remap_nc = config['remap_file']
    resources:
        mem_mb = 10000,
        runtime = 10
    run:
        remap_forcing_to_shp.remap_with_easymore(config, input.input_file,input.input_shp, output.remap_nc, only_create_remap_nc=True)

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
        mem_mb = 10000,
        runtime = 10
    run:
        remap_forcing_to_shp.remap_with_easymore(config, input.input_forcing,input.input_shp,input.remap_nc,only_create_remap_nc=False,file_path=params.file_path)
