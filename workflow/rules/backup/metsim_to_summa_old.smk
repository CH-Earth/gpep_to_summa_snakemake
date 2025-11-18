from pathlib import Path
import sys
# Import custom functions
#Update the current working directory to the directory of the file
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))
print(f"Scripts path: {scripts_path}")
import gpep_to_summa_utils as gts_utils

# Resolve paths from the configuration file
config = gts_utils.resolve_paths(config)

input_forcing_list = gts_utils.list_files_in_subdirectory(config['metsim_output_dir'], '.nc')

rule metsim_to_summa:
    input:
        expand(Path(config['summa_forcing_dir'],"{id}.nc"), id=input_forcing_list),
        output_ens_file = Path(config['summa_forcing_dir'],"summa_mean_input.nc")

rule create_hru_id_file:
    input:
        subset_domain_file = Path(config["metsim_dir"], config["metsim_domain_nc"])
    output:
        hru_id_file = temp(Path(config["metsim_dir"], 'hruId.nc'))
    group:
        "gpep_to_summa"
    resources:
        runtime=5,
        mem_mb=5000
    shell:
        'ncks -v hruId {input.subset_domain_file} {output.hru_id_file}'

rule append_hru_id_and_datastep_to_metsim_output:
    input:
        input_metsim_file = Path(config['metsim_output_dir'],"{id}.nc"),
        hru_id_file = Path(config["metsim_dir"], 'hruId.nc')
    output:
        output_metsim_file_temp = temp(Path(config['summa_forcing_dir'],"{id}_temp.nc")),
        output_metsim_file = Path(config['summa_forcing_dir'],"{id}.nc")
    group:
        "gpep_to_summa"
    resources:
        runtime=5,
        mem_mb=5000
    params:
        timestep = int(config["metsim_timestep_minutes"]) * 60
    shell:
        """
        if ! ncdump -h {input.input_metsim_file} | grep -q "hruId"; then
            ncks -h -A {input.hru_id_file} {input.input_metsim_file}
        fi
        ncks -O -C -x -v hru {input.input_metsim_file} {output.output_metsim_file_temp}
        ncrename -O -v .SWradAtm,SWRadAtm {output.output_metsim_file_temp}
        ncrename -O -v .LWradAtm,LWRadAtm {output.output_metsim_file_temp}
        ncap2 -s "data_step={params.timestep}" {output.output_metsim_file_temp} --append {output.output_metsim_file}
        """

rule calculate_mean:
    input:
        input_ens_files = expand(Path(config['summa_forcing_dir'],"{id}.nc"), id=input_forcing_list)
    output:
        output_ens_file = Path(config['summa_forcing_dir'],"summa_mean_input.nc")
    shell:
        """
        cdo -O ensmean {input.input_ens_files} {output.output_ens_file}
        """
        
