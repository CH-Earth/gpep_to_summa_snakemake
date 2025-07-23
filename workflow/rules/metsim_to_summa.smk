from pathlib import Path
import sys

# ensure your custom utils are on PYTHONPATH
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))

import gpep_to_summa_utils as gts_utils

# resolve all your paths
config = gts_utils.resolve_paths(config)
input_forcing_list = gts_utils.list_files_in_subdirectory(
    config['metsim_output_dir'], '.nc'
)
wind_file = Path(config['wind_output_dir'], "final", "wind.nc")

rule metsim_to_summa:
    input:
        Path(config['summa_forcing_dir'], "summa_mean_forcing.nc")

rule create_hru_id_file:
    """
    Extract hruId from your domain file once.
    """
    input:
        domain = Path(config["metsim_dir"], config["metsim_domain_nc"])
    output:
        hru_id = Path(config["metsim_dir"], "hruId.nc")
    resources:
        mem_mb = 10000,
        runtime = 1
    shell:
        "ncks -v hruId {input.domain} {output.hru_id}"


rule build_summa_forcing:
    """
    For each {id}.nc:
      1) copy to tmp
      2) append hruId (if missing)
      3) rename radiation vars
      4) add data_step
      5) pull in wind (knots→m/s), aligning timesteps and fixing attributes
      6) mv tmp → final
    """
    input:
        mets   = Path(config['metsim_output_dir'], "{id}.nc"),
        hru_id = rules.create_hru_id_file.output.hru_id,
        wind   = wind_file
    output:
        out           = Path(config['summa_forcing_dir'], "{id}.nc"),
        out_tmp       = temp(Path(config['summa_forcing_dir'], "{id}.nc.tmp")),
        out_tmp_wind  = temp(Path(config['summa_forcing_dir'], "{id}.nc.wind.tmp"))
    params:
        dt = int(config["metsim_timestep_minutes"]) * 60
    resources:
        mem_mb  = 10000,
        runtime = 5
    shell:
        r"""
        # 1) copy to tmp
        cp {input.mets} {output.out_tmp}

        # 2) add hruId if missing
        if ! ncdump -h {output.out_tmp} | grep -q 'hruId'; then
            ncks -h -A {input.hru_id} {output.out_tmp}
        fi

        # 3) rename radiation vars
        ncrename -O -v '.SWradAtm,SWRadAtm' {output.out_tmp}
        ncrename -O -v '.LWradAtm,LWRadAtm' {output.out_tmp}



        # 4) replace windspd from wind file (time-aligned, new metadata)
        if [ -f "{input.wind}" ]; then

            # a) grab only time & wind
            ncks -O -C -v time,CaSR_v3.1_P_UVC_10m \
                {input.wind} {output.out_tmp_wind}

            # b) rename dims & coords
            ncrename -O -d ID,hru {output.out_tmp_wind}
            ncrename -O -v .ID,hru  {output.out_tmp_wind}

            # c) rename var & convert to m/s
            ncrename -O -v 'CaSR_v3.1_P_UVC_10m,windspd' {output.out_tmp_wind}
            ncap2 -O -s 'windspd=windspd*0.51444444' \
                {output.out_tmp_wind} {output.out_tmp_wind}

            # d) fix windspd metadata
            ncatted -O -a units,windspd,o,c,"m s-1" \
                    -a long_name,windspd,o,c,"wind speed" \
                    -a standard_name,windspd,o,c,"wind_speed" \
                    {output.out_tmp_wind}

            # e) convert time to minutes & update its units
            ncap2 -O -s 'time=time*60' \
                {output.out_tmp_wind} {output.out_tmp_wind}
            ncatted -O -a units,time,o,c,"minutes since 2000-01-01 00:00:00" \
                {output.out_tmp_wind}

            # f) strip the old windspd from your tmp forcing file
            ncks -O -C -x -v windspd {output.out_tmp} {output.out_tmp}.nowspd
            mv {output.out_tmp}.nowspd {output.out_tmp}

            # g) merge in the new windspd
            cdo -O merge {output.out_tmp} {output.out_tmp_wind} {output.out_tmp}.mrg
            mv {output.out_tmp}.mrg {output.out_tmp}

        else
            echo "[build_summa_forcing] no wind file at {input.wind}, skipping"
            touch {output.out_tmp_wind}
        fi

        # 5) add data_step (append only)
        # create and copy to temp file
        tmp_file="{output.out_tmp}.bak"
        cp {output.out_tmp} "$tmp_file"
        ncap2 -O -s "data_step={params.dt}" "$tmp_file" {output.out_tmp}
        rm "$tmp_file"

        # 6) finalize
        cp {output.out_tmp} {output.out}
    """
rule calculate_mean:
    """
    Compute the ensemble mean over all {id}.nc
    """
    input:
        expand(Path(config['summa_forcing_dir'], "{id}.nc"), id=input_forcing_list)
    output:
        mean = Path(config['summa_forcing_dir'], "summa_mean_forcing.nc")
    resources:
        mem_mb = 50000,
        runtime = 40
    shell:
        "cdo -O ensmean {input} {output.mean}"