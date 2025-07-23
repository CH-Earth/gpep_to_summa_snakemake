from pathlib import Path
import sys

#Update the current working directory to the directory of the file
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))
print(f"Scripts path: {scripts_path}")
import gpep_to_summa_utils as gts_utils
import metsim_utils as ms_utils

# Resolve paths from the configuration file
config = gts_utils.resolve_paths(config)


# Generate list of forcing files to create wildcards
input_forcing_list = gts_utils.list_files_in_subdirectory(Path(config['metsim_input_dir']))

# Remove any files that contain the string "state", note that it is a list of PosixPaths
input_forcing_list = [file for file in input_forcing_list if "state" not in str(file)]

print(f"Input forcing list: {input_forcing_list}")
rule run_metsim:
    input:
        expand(Path(config['metsim_output_dir'],"{id}.nc"), id=input_forcing_list)

rule update_metsim_base_time:
    input:
        metsim_input_forcing = Path(config['metsim_input_dir'],"{id}.nc")
    output:
        metsim_temp_input_forcing = temp(Path(config['metsim_input_dir'],'temp',"{id}_input.nc"))
    resources:
        runtime=1
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
        "gpep_to_summa"
    resources:
        runtime=60,
        mem_mb=10000
    run:
        ms = ms_utils.create_metsim_config(config, input.metsim_input_forcing,input.metsim_input_state,output.metsim_output_forcing)
        ms.run()
        ms_utils.rename_metsim_output(ms,output.metsim_output_forcing)

