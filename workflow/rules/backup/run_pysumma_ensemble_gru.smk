from pathlib import Path
import sys
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

print(f"SUMMA_OUTPUT_DIR: {SUMMA_OUTPUT_DIR}")
print(f"CASE_NAME: {CASE_NAME}")
print(f"FILE_MANAGER: {FILE_MANAGER}")
print(f"SUMMA_EXE: {SUMMA_EXE}")
print(f"SUMMA_FORCING_DIR: {SUMMA_FORCING_DIR}")

list_of_grus = [89,141,184,245]
grus = []
for gru in list_of_grus:
    str_num = str(gru).zfill(4)
    grus.append(str_num)
print(f"grus: {grus}")
# ---------------------------------------
# Ensemble members (subdirectories)
# ---------------------------------------
ens, _ = gts_utils.build_ensemble_list(SUMMA_FORCING_DIR)
ens = [member.replace(" ", "").replace("\n", "").replace("\t", "") for member in ens if member]

print(f"Ensemble members: {ens}")

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
        expand(str(SUMMA_OUTPUT_DIR / (CASE_NAME + "_{ens_member}_G{gru}-{gru}_timestep.nc")), ens_member=ens, gru=grus)

rule run_summa_ensemble_simulations:
    input:
        file_manager = FILE_MANAGER,
        forcing_file = lambda wildcards: str(forcing_file_dict[wildcards.ens_member])
    output:
        summa_chunked_output = str(SUMMA_OUTPUT_DIR / (CASE_NAME + "_{ens_member}_G{gru}-{gru}_timestep.nc"))
    params:
        summa_exe = SUMMA_EXE,
        forcing_path = SUMMA_FORCING_DIR,
        output_path = SUMMA_OUTPUT_DIR,
        run_suffix = lambda wildcards: str(wildcards.ens_member),
        gru = lambda wildcards: str(wildcards.gru)
    resources:
        runtime=5,
        mem_mb=20000
    run:
        print(f"Running member: {wildcards.ens_member}")
        print(f"Forcing file: {input.forcing_file}")
        forcing_list_name = Path(input.forcing_file).name
        print(f"Forcing list name: {forcing_list_name}")

        sim = ps.Simulation(params.summa_exe, input.file_manager)

        sim_config = {
            'manager': {
                'forcingPath': str(params.forcing_path / params.run_suffix) + '/',
                #'forcingPath': str(params.forcing_path) + '/',
                'outputPath': str(params.output_path) + '/',
                'forcingListFile': str(forcing_list_name)
            }
        }

        sim.apply_config(sim_config)

        forcing_list_output = Path(Path(input.file_manager).parent,'.pysumma',params.run_suffix,forcing_list_name)
        print(f"Copying forcing list from {input.forcing_file} to {forcing_list_output}")
        forcing_list_output.parent.mkdir(parents=True, exist_ok=True)
        time.sleep(2)  # Small delay to help file system sync
        shutil.copy(input.forcing_file,forcing_list_output)
        
        sim.run(run_suffix=str(params.run_suffix), startGRU=int(params.gru), countGRU=1)
        print(sim.status)
