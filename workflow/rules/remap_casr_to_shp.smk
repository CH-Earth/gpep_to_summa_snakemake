from pathlib import Path
import sys, os
import re
import unicodedata
from pathlib import Path

# 1) Set up scripts path
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))

import remap_forcing_to_shp
import gpep_to_summa_utils as gts_utils

config = gts_utils.resolve_paths(config)

_WS_RE = re.compile(r'\s+', re.UNICODE)

# remove whitespace touching common punctuation/separators (left or right)
_ADJ_PUNCT_RE_L = re.compile(r'\s+([_\.\-\\/])')
_ADJ_PUNCT_RE_R = re.compile(r'([_\.\-\\/])\s+')

def _normalize_text(s: str) -> str:
    # Unicode normalize (NBSP, compatibility forms, etc.)
    s = unicodedata.normalize("NFKC", s)
    # Map any Unicode whitespace to a simple space
    s = ''.join(' ' if (c.isspace() or unicodedata.category(c).startswith('Z')) else c for c in s)
    # Collapse runs and trim ends
    s = _WS_RE.sub(' ', s).strip()
    # Critically: remove whitespace stuck to punctuation (prevents "name _x", "name .nc")
    s = _ADJ_PUNCT_RE_L.sub(r'\1', s)
    s = _ADJ_PUNCT_RE_R.sub(r'\1', s)
    return s

def clean_str(x):
    """Robustly clean strings, Paths, bytes, and basic containers."""
    if x is None:
        return None
    if isinstance(x, bytes):
        x = x.decode('utf-8', errors='replace')
    if isinstance(x, str):
        return _normalize_text(x)
    if isinstance(x, Path):
        anchor = x.anchor
        parts = [_normalize_text(p) for p in x.parts]
        if anchor:
            cleaned = [p for p in parts if p != anchor]
            return Path(anchor, *cleaned)
        return Path(*parts)
    if isinstance(x, (list, tuple)):
        return type(x)(clean_str(v) for v in x)
    if isinstance(x, dict):
        return {k: clean_str(v) for k, v in x.items()}
    return x

case       = clean_str(config["case_name"])
grid_in    = Path(clean_str(config["grid_input_dir"]))
grid_out   = Path(clean_str(config["grid_output_dir"]))
grid_tmp   = Path(clean_str(config["grid_temp_dir"]))
shp_file   = Path(clean_str(config["catchment_shp"]))
metsim_dir = Path(clean_str(config.get("metsim_dir", "")))

remap_fn = clean_str(f'{case}_remapped.nc')
proc_fn  = clean_str(f'{case}_processed.nc')

config["easymore_output_dir"] = grid_out
config["easymore_temp_dir"]   = grid_tmp
config["easymore_input_dir"]  = grid_in

remapping_fn = clean_str(f'{case}_remapping.nc')
config["remap_file"]          = grid_tmp / remapping_fn
config["easymore_input_var"]  = clean_str(config["grid_input_var"])

# 4) Define all derived Path variables
tmp_out_dir       = clean_str(grid_out / "final")
remapping_file    = Path(clean_str(config["remap_file"]))

final_nc_path     = clean_str(tmp_out_dir / remap_fn)
processed_nc_path = clean_str(tmp_out_dir / proc_fn)
hru_id_nc_path    = clean_str(metsim_dir / "hruId.nc")

# 5) List input files and pick first one
file_path_list = sorted(p.name for p in grid_in.iterdir() if p.suffix == ".nc")
first_file     = clean_str(grid_in / file_path_list[0])

# 5) List input files and pick first one
file_path_list = sorted(p.name for p in grid_in.iterdir() if p.suffix == ".nc")
first_file     = clean_str(grid_in / file_path_list[0])

# 6) Snakemake rules
rule all:
    input:
        processed_nc_path

rule create_hru_id_file:
    """
    Extract hruId from your domain file once.
    """
    input:
        domain = metsim_dir / clean_str(config.get("attribute_nc", ""))
    output:
        hru_id = hru_id_nc_path
    resources:
        mem_mb = 10000,
        runtime = 1
    shell:
        "ncks -v hruId {input.domain} {output.hru_id}"

rule create_remap_file:
    input:
        input_file = first_file,
        input_shp  = shp_file
    output:
        remap_nc = remapping_file
    params:
        var_lat = 'lat',
        var_lon = 'lon'
    resources:
        mem_mb = 10000,
        runtime = 5
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
        remap_nc      = remapping_file
    output:
        output_forcing = config["easymore_output_dir"] / "{id}"
    params:
        var_lat = 'lat',
        var_lon = 'lon'
    resources:
        mem_mb = 10000,
        runtime = 5
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

rule merge_grid_files_with_cdo:
    input:
        expand(config["easymore_output_dir"] / "{id}", id=file_path_list)
    output:
        output_file = temp(final_nc_path)
    resources:
        mem_mb = 10000,
        runtime = 10
    shell:
        "cdo mergetime {input} {output}"

rule process_netcdf:
    input:
        merged = final_nc_path,
        hru_id = hru_id_nc_path
    output:
        processed = processed_nc_path
    params:
        dt = 3600
    resources:
        mem_mb = 10000,
        runtime = 20
    shell:
        r"""
        tmp=$(mktemp -u {output.processed}.XXXX.nc)
        cp {input.merged} "$tmp"

        # 1) Rename variables
        ncrename -O -v CaSR_v3.1_A_PR0_SFC,pptrate "$tmp"
        ncrename -O -v CaSR_v3.1_P_TT_1.5m,airtemp "$tmp"
        ncrename -O -v CaSR_v3.1_P_P0_SFC,airpres "$tmp"
        ncrename -O -v CaSR_v3.1_P_FI_SFC,LWRadAtm "$tmp"
        ncrename -O -v CaSR_v3.1_P_FB_SFC,SWRadAtm "$tmp"
        ncrename -O -v CaSR_v3.1_P_HU_1.5m,spechum "$tmp"
        ncrename -O -v CaSR_v3.1_P_UVC_10m,windspd "$tmp"

        # 2) Update units
        ncatted -O -a units,pptrate,o,c,'millimeter / second' "$tmp"
        ncatted -O -a units,airtemp,o,c,'kelvin' "$tmp"
        ncatted -O -a units,airpres,o,c,'pascal' "$tmp"
        ncatted -O -a units,LWRadAtm,o,c,'watt / meter ** 2' "$tmp"
        ncatted -O -a units,SWRadAtm,o,c,'watt / meter ** 2' "$tmp"
        ncatted -O -a units,spechum,o,c,'dimensionless' "$tmp"
        ncatted -O -a units,windspd,o,c,'meter / second' "$tmp"

        # 3) Convert values
        ncap2 -O -s "pptrate=pptrate*1000/3600" "$tmp" "$tmp"
        ncap2 -O -s "airtemp=airtemp+273.15" "$tmp" "$tmp"
        ncap2 -O -s "airpres=airpres*100" "$tmp" "$tmp"
        ncap2 -O -s "windspd=windspd*0.514444" "$tmp" "$tmp"

        # 4) Shift time to UTC
        tmp_shift=$(mktemp -u {output.processed}.shift.XXXX.nc)
        tmp_ref=$(mktemp -u {output.processed}.ref.XXXX.nc)
        cdo -L -shifttime,-12hours "$tmp" "$tmp_shift"
        cdo setreftime,1970-01-01,00:00:00,1hour "$tmp_shift" "$tmp_ref"
        cp "$tmp_ref" "$tmp"
        rm -f "$tmp_shift" "$tmp_ref"

        # 5) Rename variables
        ncrename -O -d .ID,hru "$tmp"
        ncrename -O -v .ID,hruId "$tmp"

        # 6) Add hruId as a coordinate variable
        ncap2 -O -s "hru=array(0,1,hruId)" "$tmp" "$tmp"
        ncatted -O -a long_name,hru,a,c,"hru coordinate index" "$tmp" "$tmp"

        # 7) Add data_step
        bak=$(mktemp -u {output.processed}.bak.XXXX.nc)
        cp "$tmp" "$bak"
        ncap2 -O -s "data_step={params.dt}" "$bak" {output.processed}

        rm -f "$bak" "$tmp"
        """
