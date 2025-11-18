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

# Create a list of the forcing files produced in the last workflow
input_forcing_list = gts_utils.list_files_in_subdirectory(config['easymore_output_dir'])
ensemble_list, file_path_list = gts_utils.build_ensemble_list(config['gpep_forcing_dir'])

easymore_output = Path(config['easymore_output_dir']) 
metsim_input = Path(config['metsim_input_dir'])

# Main rule to define the files produced by this workflow
rule metsim_file_prep:
    input:
        Path(config["metsim_dir"], config["metsim_domain_nc"]),
        expand(Path(metsim_input,"{forcing}.nc"), forcing = input_forcing_list),
        expand(Path(metsim_input,"{forcing}_state.nc"), forcing = input_forcing_list)
         
# Create metsim domain file from an existing summa attribute file
rule create_metsim_domain_summa_attr:
    input:
        attr_nc = Path(config["attribute_nc"])
    output:
        domain_nc = Path(config["metsim_dir"], config["metsim_domain_nc"])
    resources:
        runtime=2,
        mem_mb=1000
    shell:
        'ncap2 -O -s "mask=elevation*0+1" {input.attr_nc} {output.domain_nc}'

# Define rule to run file remapping when remap file exists
rule prep_forcing_files_with_hru_id:
    input:
        input_forcing = Path(easymore_output,"{forcing}.nc")
    output:
        hru_id_temp = temp(Path(easymore_output,"{forcing}_hruId.nc")),
        hru_id = Path(metsim_input,"{forcing}.nc")
    group:
        "metsim_file_prep"
    resources:
        runtime=2,
        mem_mb=1000
    shell:
        """
        ncrename -O -v .ID,hruId {input.input_forcing}
        ncrename -d .ID,hru {input.input_forcing}
        ncap2 -O -s "hru=array(0,1,hruId)" {input.input_forcing} {output.hru_id_temp};
        ncatted -O -a long_name,hru,a,c,"hru coordinate index" {output.hru_id_temp} {output.hru_id}
        """
        
rule create_state_file:
    input:
        input_forcing_file = Path(metsim_input,"{forcing}.nc")
    output:
        output_state_file = temp(Path(metsim_input,"{forcing}_state_temp.nc"))
    group:
        "metsim_file_prep"
    resources:
        runtime=2,
        mem_mb=1000
    run:
        ms_utils.create_state_file(input.input_forcing_file, output.output_state_file)

rule update_state_file_time:
    input:
        input_state_file = Path(metsim_input,"{forcing}_state_temp.nc")
    output:
        output_state_file = Path(metsim_input,"{forcing}_state.nc")
    group:
        "metsim_file_prep"
    resources:
        runtime=2,
        mem_mb=1000
    run:
        gts_utils.update_time_units(input.input_state_file, output.output_state_file)

