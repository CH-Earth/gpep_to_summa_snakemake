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
from pprint import pprint


# Resolve paths from the configuration file
config = gts_utils.resolve_paths(config)

config['easymore_input_dir'] = config['wind_input_dir']
config['easymore_output_dir'] = config['wind_output_dir']
config['easymore_temp_dir'] = config['wind_temp_dir']
remapping_file_name = config['case_name'] + "_remapping.nc"
config['remap_file'] = Path(config['wind_temp_dir'],remapping_file_name)
config['easymore_input_var'] = config['wind_input_var']

file_path_list = sorted([f for f in os.listdir(config['wind_input_dir']) if f.endswith('.nc')])
print(f"File path list: {file_path_list}")

# Now you can grab the "first" entry
first_file = (Path(config['wind_input_dir']) / file_path_list[0]).with_suffix(".nc")
print(f"First file: {first_file}")

rule all:
    input:
       Path(config['wind_output_dir'],"final", "wind.nc")

# Define rule to run file remapping when remap file exists
rule create_remap_file:
    input:
        input_file = first_file,
        input_shp = config['catchment_shp']
    output:
        remap_nc = config['remap_file']
    params:
        var_lat = 'lat',
        var_lon = 'lon'
    resources:
        mem_mb = 10000,
        runtime = 5
    run:
        remap_forcing_to_shp.remap_with_easymore(config, input.input_file,input.input_shp, output.remap_nc, only_create_remap_nc=True, var_lat=params.var_lat, var_lon=params.var_lon)

# Define rule to run file remapping when remap file exists
rule remap_with_easymore:
    input:
        input_forcing = Path(config['easymore_input_dir'],"{id}"),
        input_shp = config['catchment_shp'],
        remap_nc = config['remap_file'],
    output:
        output_forcing = Path(config['easymore_output_dir'],"{id}")
    params:
        var_lat = 'lat',
        var_lon = 'lon'
    resources:
        mem_mb = 10000,
        runtime = 5
    run:
        remap_forcing_to_shp.remap_with_easymore(config, input.input_forcing,input.input_shp,input.remap_nc,only_create_remap_nc=False,var_lat=params.var_lat, var_lon=params.var_lon)

rule merge_wind_files_with_cdo:
    input:
        expand(Path(config['easymore_output_dir'],"{id}"), id=file_path_list)
    output:
        output_file = Path(config['wind_output_dir'],"final", "wind.nc")
    resources:
        mem_mb = 10000,
        runtime = 10
    shell:
        "cdo mergetime {input} {output}"