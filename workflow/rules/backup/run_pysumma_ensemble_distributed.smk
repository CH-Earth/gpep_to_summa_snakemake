from pathlib import Path
import sys
import pysumma as ps
import shutil
import time

# Add scripts to path
workflow_dir = Path(workflow.basedir).parent
scripts_path = workflow_dir / "scripts"
sys.path.append(str(scripts_path))
import gpep_to_summa_utils as gts_utils
import ss_utils

def clean_path(p):
    return Path(str(p).replace(" ", "").replace("\n", "").replace("\t", ""))

def parse_gru_chunk(chunk_str):
    chunk_str = chunk_str.strip('G')
    start_str, end_str = chunk_str.split('-')
    gru_start = int(start_str)
    gru_end = int(end_str)
    gru_count = gru_end - gru_start + 1
    return gru_start, gru_count
# --- Load configuration ---
config = gts_utils.resolve_paths(config)
CASE_NAME         = config['case_name']
SUMMA_OUTPUT_DIR  = clean_path(config['summa_output_dir'])
FILE_MANAGER      = clean_path(config['file_manager'])
SUMMA_EXE         = clean_path(config['summa_exe'])
SUMMA_FORCING_DIR = clean_path(config['summa_forcing_dir'])
ATTRIBUTES_NC     = clean_path(config['attribute_nc'])
GRU_CHUNK_SIZE    = int(config.get('gru_chunk_size', 50))

# --- Ensemble members ---
ens, _ = gts_utils.build_ensemble_list(SUMMA_FORCING_DIR)
ens = [member.strip() for member in ens if member]

# --- GRU chunks ---
num_grus = ss_utils.calc_num_grus(ATTRIBUTES_NC)
gru_chunks = ss_utils.generate_gru_start_and_count(num_grus, chunk_size=GRU_CHUNK_SIZE)

# --- Forcing file lists per ensemble member ---
file_paths = gts_utils.list_files_in_subdirectory(SUMMA_FORCING_DIR, filenames_only=True)
file_paths = [clean_path(fp) for fp in file_paths]

forcing_file_dict = {}
for member_id in ens:
    member_id_clean = member_id.strip()
    forcing_list_path = clean_path(SUMMA_FORCING_DIR / f"forcing_files_{member_id_clean}.txt")
    forcing_list_path.parent.mkdir(parents=True, exist_ok=True)
    with forcing_list_path.open('w') as fh:
        for fp in file_paths:
            if member_id_clean in fp.parts[0]:
                full_path = Path(*[part.strip() for part in fp.parts]).with_suffix(".nc")
                fh.write(full_path.as_posix() + "\n")
    forcing_file_dict[member_id_clean] = forcing_list_path

# --- Rules ---

rule all:
    input:
        expand(
            str(SUMMA_OUTPUT_DIR / "{case_name}_{ens_member}_timestep.nc"),
            case_name=CASE_NAME,
            ens_member=ens
        )

# Run SUMMA in GRU chunks for each ensemble member
rule run_summa_ensemble_gru_chunk:
    input:
        file_manager = FILE_MANAGER,
        forcing_file = lambda wildcards: str(forcing_file_dict[wildcards.ens_member])
    output:
        chunked_output = temp(str(SUMMA_OUTPUT_DIR / "{case_name}_{ens_member}_{gru_chunk}_timestep.nc"))
    params:
        summa_exe = SUMMA_EXE,
        forcing_path = SUMMA_FORCING_DIR,
        output_path = SUMMA_OUTPUT_DIR,
        run_suffix = lambda wildcards: wildcards.ens_member,
        gru_start = lambda wildcards: parse_gru_chunk(wildcards.gru_chunk)[0],
        gru_count = 50
        #gru_count = lambda wildcards: parse_gru_chunk(wildcards.gru_chunk)[1]
    resources:
        runtime=1,
        mem_mb=1000
    run:
        print(f"Running member: {wildcards.ens_member}, GRU chunk: {wildcards.gru_chunk}")
        sim = ps.Simulation(params.summa_exe, input.file_manager)
        sim_config = {
            'manager': {
                'forcingPath': str(params.forcing_path / params.run_suffix) + '/',
                'outputPath': str(params.output_path) + '/',
                'forcingListFile': Path(input.forcing_file).name
            }
        }
        sim.apply_config(sim_config)

        # Copy forcing list for SUMMA
        forcing_list_name = Path(input.forcing_file).name
        forcing_list_output = Path(Path(input.file_manager).parent, '.pysumma', params.run_suffix, forcing_list_name)
        forcing_list_output.parent.mkdir(parents=True, exist_ok=True)
        time.sleep(2)
        if not forcing_list_output.exists():
            shutil.copy(input.forcing_file, forcing_list_output)
            print(f"Copied to {forcing_list_output}")
        else:
            print(f"{forcing_list_output} already exists, skipping copy.")

        attribute_file_check = Path(Path(input.file_manager).parent,'attributes.nc')

        if attribute_file_check.exists():
            print(f"Attribute file found at: {attribute_file_check}")
            # Actually run SUMMA for this GRU chunk!
            sim.run(
                run_suffix=str(params.run_suffix),
                startGRU=params.gru_start,
                countGRU=params.gru_count,
                write_config=False
            )
            print(sim.status)
            print(sim.stdout)
            print(sim.stderr)
        else:
            print(f"Warning: Attribute file not found at: {attribute_file_check}")
            print("Will set write_config=True to generate default configuration")
            # Actually run SUMMA for this GRU chunk!
            sim.run(
                run_suffix=str(params.run_suffix),
                startGRU=params.gru_start,
                countGRU=params.gru_count,
                write_config=True
            )
            print(sim.status)
            print(sim.stdout)
            print(sim.stderr)

# Merge SUMMA chunk outputs for each ensemble member
rule merge_summa_output_files:
    input:
        chunked_outputs = expand(
            str(SUMMA_OUTPUT_DIR / "{case_name}_{ens_member}_{gru_chunk}_timestep.nc"),
            case_name=CASE_NAME,
            ens_member="{ens_member}",
            gru_chunk=gru_chunks
        )
    output:
        merged_output = str(SUMMA_OUTPUT_DIR / "{case_name}_{ens_member}_timestep.nc")
    run:
        print(f"Merging SUMMA outputs for {wildcards.ens_member} into {output.merged_output}.")
        ss_utils.merge_netcdf_output(input.chunked_outputs, output.merged_output)