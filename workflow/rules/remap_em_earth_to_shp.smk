from pathlib import Path
import sys, os

# 1) Set up scripts path
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))

import remap_forcing_to_shp
import gpep_to_summa_utils as gts_utils

config = gts_utils.resolve_paths(config)

def clean_str(x):
    return x.strip() if isinstance(x, str) else x


case       = clean_str(config["case_name"])
grid_in    = Path(clean_str("/anvil/projects/x-ees240082/data/CONFLUENCE_data/domain_Bow_at_Banff_TDX/forcing/raw_data_em_earth/"))
grid_out   = Path(clean_str("/anvil/projects/x-ees240082/data/CONFLUENCE_data/domain_Bow_at_Banff_TDX/forcing/easymore_output/"))
grid_tmp   = Path(clean_str("/anvil/projects/x-ees240082/data/CONFLUENCE_data/domain_Bow_at_Banff_TDX/forcing/temp/"))
remapping_final_dir = Path(clean_str("/anvil/projects/x-ees240082/data/CONFLUENCE_data/domain_Bow_at_Banff_TDX/forcing/em_earth_remapped/"))
shp_file   = Path(clean_str(config["catchment_shp"]))

# 3) Populate easymore settings in config\ńconfig["easymore_input_dir"]  = grid_in
config["easymore_output_dir"] = grid_out
config["easymore_temp_dir"]   = grid_tmp
config["easymore_input_dir"]  = grid_in

remapping_fn = 'bow_remapping.nc'
config["remap_file"]          = Path(grid_tmp, 'bow_remapping.nc')
config["easymore_input_var"]  = ['prcp','prcp_corrected','tmean']

# List input files and pick first one
file_path_list = sorted(p.name for p in grid_in.iterdir() if p.suffix == ".nc")
first_file     = grid_in / file_path_list[0]

# 6) Snakemake rules
rule all:
    input:
       expand(Path(remapping_final_dir, "remapped_{id}"), id=file_path_list)

rule create_remap_file:
    input:
        input_file = first_file,
        input_shp  = shp_file
    output:
        remap_nc = config["remap_file"]
    params:
        var_lat = 'lat',
        var_lon = 'lon'
    resources:
        mem_mb = 1000,
        runtime = 1
    run:
        remap_forcing_to_shp.remap_with_easymore(
            config,
            input.input_file,
            input.input_shp,
            output.remap_nc,
            only_create_remap_nc=True,
            var_lat=params.var_lat,
            var_lon=params.var_lon
        )

rule remap_with_easymore:
    input:
        input_forcing = config["easymore_input_dir"] / "{id}",
        input_shp     = shp_file,
        remap_nc      = config["remap_file"]
    output:
        output_forcing = config["easymore_output_dir"] / "{id}"
    params:
        var_lat = 'lat',
        var_lon = 'lon'
    resources:
        mem_mb = 1000,
        runtime = 1
    run:
        remap_forcing_to_shp.remap_with_easymore(
            config,
            input.input_forcing,
            input.input_shp,
            input.remap_nc,
            only_create_remap_nc=False,
            var_lat=params.var_lat,
            var_lon=params.var_lon
        )

rule update_nc_id_with_hru:
    input:
        input_forcing = config["easymore_output_dir"] / "{id}"
    output:
        temp_forcing = temp(config["easymore_output_dir"] / "temp_{id}"),
        output_forcing = Path(remapping_final_dir, "remapped_{id}")
    resources:
        mem_mb = 1000,
        runtime = 1
    shell:
        """
        ncrename -O -v .ID,hruId {input.input_forcing}
        ncrename -d .ID,hru {input.input_forcing}
        ncap2 -O -s "hru=array(0,1,hruId)" {input.input_forcing} {output.temp_forcing};
        ncatted -O -a long_name,hru,a,c,"hru coordinate index" {output.temp_forcing} {output.output_forcing}
        """

