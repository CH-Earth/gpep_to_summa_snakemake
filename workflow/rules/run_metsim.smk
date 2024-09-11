from pathlib import Path

sys.path.append('../')
from scripts import gpep_to_summa_utils as gts_utils
from scripts import metsim_utils as ms_utils

# Resolve paths from the configuration file
config = gts_utils.resolve_paths(config)


# Generate list of forcing files to create wildcards
input_forcing_list = gts_utils.list_files_in_subdirectory(Path(config['easymore_output_dir']))

rule run_metsim:
    input:
        expand(Path(config['metsim_output_dir'],"{id}.nc"), id=input_forcing_list)

rule update_metsim_base_time:
    input:
        metsim_input_forcing = Path(config['metsim_input_dir'],"{id}.nc")
    output:
        metsim_temp_input_forcing = temp(Path(config['metsim_input_dir'],'temp',"{id}_input.nc"))
    group:
        "run_metsim"
    resources:
        runtime= 2,
        mem_mb= 500
    run:
        gts_utils.update_time_units(input.metsim_input_forcing, output.metsim_temp_input_forcing)

rule generate_metsim_output:
    input:
        metsim_input_forcing = Path(config['metsim_input_dir'],'temp',"{id}_input.nc"),
        metsim_input_state = Path(config['metsim_input_dir'],"{id}_state.nc"),
        metsim_input_domain = Path(config["metsim_dir"], config["metsim_domain_nc"])
    output:
        metsim_output_forcing = Path(config['metsim_output_dir'],"{id}.nc")
    group:
        "run_metsim"
    resources:
        runtime= 30,
        mem_mb=10000
    run:
        ms = ms_utils.create_metsim_config(config, input.metsim_input_forcing,input.metsim_input_state,output.metsim_output_forcing)
        ms.run()
        ms_utils.rename_metsim_output(ms,output.metsim_output_forcing)

