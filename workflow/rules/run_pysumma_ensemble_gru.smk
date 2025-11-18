from pathlib import Path
import sys
import pandas as pd
import pysumma as ps
import xarray as xr
import shutil
import time

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
SUMMA_OUTPUT_DIR  = clean_path(config['summa_output_dir'])
FILE_MANAGER      = clean_path(config['file_manager'])
SUMMA_EXE         = clean_path(config['summa_exe'])
SUMMA_FORCING_DIR = clean_path(config['summa_forcing_dir'])
SUMMA_SETTINGS_DIR = clean_path(config['summa_settings_dir'])
BASE_10_GRUS = int(config['base_10_grus'])

SUMMA_GRU_CSV = clean_path(config['summa_gru_csv'])


gru_df = pd.read_csv(SUMMA_GRU_CSV)
list_of_grus_full = gru_df['GRU_ID'].apply(lambda x: int(float(x))).tolist()
list_of_grus = sorted(set(list_of_grus_full)) 

# Get max GRU to determine padding width
max_gru = max(list_of_grus)
zfill_width = len(str(max_gru))

# Build zero-padded GRU strings

grus = [str(gru).zfill(BASE_10_GRUS) for gru in list_of_grus]

print(f"SUMMA_OUTPUT_DIR: {SUMMA_OUTPUT_DIR}")
print(f"CASE_NAME: {CASE_NAME}")
print(f"FILE_MANAGER: {FILE_MANAGER}")
print(f"SUMMA_EXE: {SUMMA_EXE}")
print(f"SUMMA_FORCING_DIR: {SUMMA_FORCING_DIR}")
print(f"SUMMA_SETTINGS_DIR: {SUMMA_SETTINGS_DIR}")
print(f"SUMMA_GRU_CSV: {SUMMA_GRU_CSV}")
print(f"grus: {grus}")


# ---------------------------------------
# Ensemble members (subdirectories)
# ---------------------------------------
ens, _ = gts_utils.build_ensemble_list(SUMMA_FORCING_DIR)
ens = [member.replace(" ", "").replace("\n", "").replace("\t", "") for member in ens if member]

# ---------------------------------------
# List forcing files (cleaned, Path objects)
# ---------------------------------------
file_paths = gts_utils.list_files_in_subdirectory(SUMMA_FORCING_DIR, filenames_only=True)
file_paths = [clean_path(fp) for fp in file_paths]
print(f"File paths: {[str(fp) for fp in file_paths]}")

# ---------------------------------------
# Create forcing_files_<member>.txt for SUMMA
# ---------------------------------------
forcing_file_dict = {}
for member_id in ens:
    member_id_clean = member_id.replace(" ", "").replace("\n", "").replace("\t", "")
    forcing_list_path = clean_path(SUMMA_FORCING_DIR / f"forcing_files_{member_id_clean}.txt")
    print(f"Forcing list path: {forcing_list_path}")
    forcing_list_path.parent.mkdir(parents=True, exist_ok=True)

    with forcing_list_path.open('w') as fh:
        for fp in file_paths:
            if member_id_clean in fp.parts[0].strip():
                full_path = Path(*[part.strip() for part in fp.parts]).with_suffix(".nc")
                fh.write(full_path.as_posix().strip() + "\n")
    print(f'member_id_clean: {member_id_clean}')
    print(f'forcing_list_path: {forcing_list_path}')
    forcing_file_dict[member_id_clean] = forcing_list_path

# ---------------------------------------
# Snakemake Rules
# ---------------------------------------
rule all:
    input:
        expand(str(SUMMA_SETTINGS_DIR / ".pysumma" / "{ens_member}" / "attributes.nc"), ens_member=ens), 
        expand(str(SUMMA_OUTPUT_DIR / (CASE_NAME + "_{ens_member}_G{gru}-{gru}_timestep.nc")), ens_member=ens, gru=grus)

rule write_summa_configuration:
    input:
        file_manager = FILE_MANAGER,
        forcing_file = lambda wildcards: str(forcing_file_dict[wildcards.ens_member])
    output:
        attributes = str(SUMMA_SETTINGS_DIR/ ".pysumma" / "{ens_member}" / "attributes.nc"),
        file_manager = str(SUMMA_SETTINGS_DIR / ".pysumma" / "{ens_member}" / "fileManager.txt")
    params:
        summa_exe = SUMMA_EXE,
        forcing_path = SUMMA_FORCING_DIR,
        output_path = SUMMA_OUTPUT_DIR,
        run_suffix = lambda wildcards: str(wildcards.ens_member),
        ens_member = lambda wildcards: str(wildcards.ens_member)
    resources:
        runtime=1,
        mem_mb=1000
    run:
        forcing_list_name = Path(input.forcing_file).name
        print(f"Writing SUMMA config for ensemble {params.run_suffix}")
        print(f"Forcing file: {forcing_list_name}")
        sim = ps.Simulation(params.summa_exe, input.file_manager)

        sim_config = {
            'manager': {
                'forcingPath': str(params.forcing_path / params.run_suffix) + '/',
                'outputPath': str(params.output_path) + '/',
                'forcingListFile': str(forcing_list_name)

            }
        }

        sim.apply_config(sim_config)
       
        sim.run_suffix = params.run_suffix
        sim._write_configuration(name=params.run_suffix)

        forcing_list_output = Path(Path(input.file_manager).parent,'.pysumma',params.ens_member,forcing_list_name)
        print(f"Copying forcing list from {input.forcing_file} to {forcing_list_output}")
        if forcing_list_output.exists():
            forcing_list_output.unlink()
        shutil.copy(input.forcing_file,forcing_list_output)

rule run_summa_ensemble_simulations:
    input:
        file_manager = str(SUMMA_SETTINGS_DIR / ".pysumma" / "{ens_member}" / "fileManager.txt"),
        attributes = str(SUMMA_SETTINGS_DIR / ".pysumma" / "{ens_member}" / "attributes.nc")
    output:
        summa_chunked_output = str(SUMMA_OUTPUT_DIR / (CASE_NAME + "_{ens_member}_G{gru}-{gru}_timestep.nc"))
    params:
        summa_exe = SUMMA_EXE,
        forcing_path = SUMMA_FORCING_DIR,
        output_path = SUMMA_OUTPUT_DIR,
        run_suffix = lambda wildcards: str(wildcards.ens_member),
        gru = lambda wildcards: str(wildcards.gru),
        ens_member = lambda wildcards: str(wildcards.ens_member)
    resources:
        runtime=10,
        mem_mb=20000
    run:
        print(f"File manager: {input.file_manager}")
        print(f"Summa exe: {params.summa_exe}")
        print(f"Run suffix: {params.run_suffix}")
        print(f"GRU: {params.gru}")
        print(f"Ensemble member: {params.ens_member}")
        print(f"Summa output dir: {params.output_path}")
        print(f"Summa forcing dir: {params.forcing_path}")
        sim = ps.Simulation(params.summa_exe, input.file_manager)
        sim.run(run_suffix=str(params.run_suffix), startGRU=int(params.gru), countGRU=1, write_config=False)
        print(sim.stderr)
        print(sim.stdout)
        print(sim.status)

