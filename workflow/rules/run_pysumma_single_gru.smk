from pathlib import Path
import sys
import pysumma as ps
import xarray as xr
import shutil
import time
import pandas as pd

# Add scripts to path
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))
import gpep_to_summa_utils as gts_utils

def clean_path(p):
    return Path(str(p).replace(" ", "").replace("\n", "").replace("\t", ""))

# ---------------------------------------
# Load configuration and resolve paths
# ---------------------------------------
config = gts_utils.resolve_paths(config)
CASE_NAME         = config['case_name']
# Add _deterministic to the output directory with is a path
SUMMA_OUTPUT_DIR  = clean_path(Path(config['summa_det_output_dir']))
FILE_MANAGER      = clean_path(config['file_manager'])
SUMMA_EXE         = clean_path(config['summa_exe'])
SUMMA_FORCING_DIR = clean_path(config['summa_det_forcing_dir'])
SUMMA_SETTINGS_DIR = clean_path(config['summa_settings_dir'])
SUMMA_GRU_CSV = clean_path(config['summa_gru_csv'])
z_fill_width = config['base_10_grus']

print(f"SUMMA_OUTPUT_DIR: {SUMMA_OUTPUT_DIR}")
print(f"CASE_NAME: {CASE_NAME}")
print(f"FILE_MANAGER: {FILE_MANAGER}")
print(f"SUMMA_EXE: {SUMMA_EXE}")
print(f"SUMMA_FORCING_DIR: {SUMMA_FORCING_DIR}")
print(f"SUMMA_SETTINGS_DIR: {SUMMA_SETTINGS_DIR}")
print(f"SUMMA_GRU_CSV: {SUMMA_GRU_CSV}")

gru_df = pd.read_csv(SUMMA_GRU_CSV)
list_of_grus_full = gru_df['GRU_ID'].apply(lambda x: int(float(x))).tolist()
list_of_grus = sorted(set(list_of_grus_full)) 

# Get max GRU to determine padding width
max_gru = max(list_of_grus)
zfill_width = len(str(max_gru))

# Build zero-padded GRU strings
grus = [str(gru).zfill(z_fill_width) for gru in list_of_grus]

# ---------------------------------------
# List forcing files (cleaned, Path objects)
# ---------------------------------------
file_paths = gts_utils.list_files_in_subdirectory(SUMMA_FORCING_DIR, filenames_only=True)
file_paths = [clean_path(fp) for fp in file_paths]
file_paths.sort(key=lambda p: p.name)
print(f"File paths: {[str(fp) for fp in file_paths]}")

# ---------------------------------------
# Create forcing_files_<member>.txt for SUMMA
# ---------------------------------------
forcing_file_dict = {}

forcing_list_path = clean_path(SUMMA_FORCING_DIR / f"forcing_files_{CASE_NAME}.txt")
print(f"Forcing list path: {forcing_list_path}")
forcing_list_path.parent.mkdir(parents=True, exist_ok=True)

with forcing_list_path.open('w') as fh:
    for fp in file_paths:
        full_path = Path(*[part.strip() for part in fp.parts]).with_suffix(".nc")
        fh.write(full_path.as_posix().strip() + "\n")
    print(f'forcing_list_path: {forcing_list_path}')
    forcing_file_dict[CASE_NAME] = forcing_list_path

# ---------------------------------------
# Snakemake Rules
# ---------------------------------------
rule all:
    input: 
        expand(str(SUMMA_OUTPUT_DIR / (CASE_NAME + "_" + CASE_NAME + "_G{gru}-{gru}_timestep.nc")), gru=grus)

rule write_summa_configuration:
    input:
        file_manager = FILE_MANAGER,
        forcing_file = lambda wildcards: str(forcing_file_dict[CASE_NAME])
    output:
        attributes = str(SUMMA_SETTINGS_DIR/ ".pysumma" / CASE_NAME / "attributes.nc"),
        file_manager = str(SUMMA_SETTINGS_DIR / ".pysumma" / CASE_NAME / "fileManager.txt")
    params:
        summa_exe = SUMMA_EXE,
        forcing_path = SUMMA_FORCING_DIR,
        output_path = SUMMA_OUTPUT_DIR,
    resources:
        runtime=1,
        mem_mb=1000
    run:
        forcing_list_name = Path(input.forcing_file).name
        print(f"Writing SUMMA config for{CASE_NAME}")
        print(f"Forcing file: {forcing_list_name}")
        sim = ps.Simulation(params.summa_exe, input.file_manager)

        sim_config = {
            'manager': {
                'forcingPath': str(params.forcing_path) + '/',
                'outputPath': str(params.output_path) + '/',
                'forcingListFile': str(forcing_list_name),
                'tmZoneInfo': 'localTime'
            }
        }

        sim.apply_config(sim_config)
       
        sim.run_suffix = CASE_NAME
        sim._write_configuration(name=CASE_NAME)

        forcing_list_output = Path(Path(input.file_manager).parent,'.pysumma',CASE_NAME,forcing_list_name)
        print(f"Copying forcing list from {input.forcing_file} to {forcing_list_output}")
        if forcing_list_output.exists():
            forcing_list_output.unlink()
        shutil.copy(input.forcing_file,forcing_list_output)

rule run_summa_single_simulations:
    input:
        file_manager = str(SUMMA_SETTINGS_DIR / ".pysumma" / CASE_NAME / "fileManager.txt"),
        attributes = str(SUMMA_SETTINGS_DIR / ".pysumma" / CASE_NAME / "attributes.nc")
    output:
        summa_chunked_output = str(SUMMA_OUTPUT_DIR / (CASE_NAME + "_" + CASE_NAME + "_G{gru}-{gru}_timestep.nc"))
    params:
        summa_exe = SUMMA_EXE,
        forcing_path = SUMMA_FORCING_DIR,
        output_path = SUMMA_OUTPUT_DIR,
        run_suffix = CASE_NAME,
        gru = lambda wildcards: str(wildcards.gru)
    resources:
        runtime=10,
        mem_mb=20000
    run:

        sim = ps.Simulation(params.summa_exe, input.file_manager)

        sim.run(run_suffix=str(params.run_suffix), startGRU=int(params.gru), countGRU=1, write_config=False)
        print(sim.status)
