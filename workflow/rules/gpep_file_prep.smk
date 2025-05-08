
# This Snakemake file prepares the GPEP  data for use in 
from pathlib import Path
import sys

sys.path.append('../')
from scripts import gpep_to_summa_utils as gts_utils

# Resolve paths from the configuration file
config = gts_utils.resolve_paths(config)

#Set the list of forcing files to process
# Read all forcing files and create a list based on the output directory (i.e. ens/filename.nc)
ensemble_list, file_path_list = gts_utils.build_ensemble_list(config['gpep_forcing_dir'])

#This first rule establishes the output files that will be created
rule gpep_file_prep:
    input:
        expand(Path(config['gpep_tmp_forcing_dir'],"{forcing_file}.nc"), forcing_file=file_path_list)

#Add greogrian calendar to the time variable, needed for easymore  
rule add_gregorian_to_nc:
    input:  
        input_forcing = Path(config['gpep_forcing_dir'],"{id}.nc")
    output: 
        output_forcing = temp(Path(config['gpep_tmp_forcing_dir'],"{id}_greg.nc"))
    group:
        "gpep_to_summa"
    shell: 
        'ncatted -a "calendar,time,o,c,"gregorian"" {input.input_forcing} {output.output_forcing}'
#Process temperature data to create t_max and t_min
rule add_t_max_and_t_min:
    input: 
        input_file = Path(config['gpep_tmp_forcing_dir'],"{id}_greg.nc")
    output:
        temp = temp(Path(config['gpep_tmp_forcing_dir'],"{ens}","{id}_temp.nc")),
        output_file = Path(config['gpep_tmp_forcing_dir'],"{ens}","{id}.nc")
    params:
        t_mean_var = 'tmean',
        t_range_var = 'trange'
    group:
        "gpep_to_summa"
    shell:
        """
        ncap2 -s 'where({params.t_range_var} < 0) {params.t_range_var}=0;' {input.input_file} {output.temp};
        ncap2 -s "t_max = {params.t_mean_var}+0.5*{params.t_range_var}" -A {output.temp};
        ncap2 -s "t_min = {params.t_mean_var}-0.5*{params.t_range_var}" -A {output.temp};
        ncatted -O -a long_name,t_max,o,c,"estimated daily maximum temperature" {output.temp};
        ncatted -O -a long_name,t_min,o,c,"estimated daily minimum temperature" {output.temp};
        ncrename -O -v .lat,latitude {output.temp};
        ncrename -O -v .lon,longitude {output.temp};
        cp {output.temp} {output.output_file}
        """
