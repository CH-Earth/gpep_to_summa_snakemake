"""

Remapping of forcing files to catchment hrus using easymore

"""

from pathlib import Path
import json
import sys

# Print the path to the Python interpreter
print("Python interpreter path:", sys.executable)
from easymore import Easymore


def rename_easymore_output(esmr):
    """Rename easymore output to original filename"""

    # Read data from easymore object.
    # Strict output formatting is done by easymore
    forcing_input = Path(esmr.source_nc)

    easymore_output = Path(
        esmr.output_dir,
        f'{esmr.case_name}_remapped_{forcing_input.name}'
    )

    easymore_output.rename(Path(easymore_output.parent,forcing_input.name))


def remap_with_easymore(
    config, input_forcing, input_shp, remap_file, only_create_remap_nc=False, file_path=None):


    # initializing EASYMORE object
    esmr = Easymore()

    if file_path is None:
        ens_member_str = None
        output_dir = str(config["easymore_output_dir"]) + "/"
    else:
        #file_path_str = file_path.pop()
        ens_member_str = file_path[:3]
        output_dir = str(config["easymore_output_dir"]) + "/" + ens_member_str + "/"

    temp_dir = str(config["easymore_temp_dir"]) + "/"
    print(temp_dir)
    easymore_input_var = config["easymore_input_var"]

    json_dict = {
        "case_name": config["case_name"],
        "temp_dir": temp_dir,
        "target_shp": input_shp,
        "source_nc": input_forcing,
        "var_names": easymore_input_var,
        "var_lon": "longitude",
        "var_lat": "latitude",
        "var_time": "time",
        "output_dir": output_dir,
        "target_shp_ID": config["catchment_shp_hru_id_field"],
        "target_shp_lat": config["catchment_shp_lat_id_field"],
        "target_shp_lon": config["catchment_shp_lon_id_field"],
        "format_list": ["f4"],
        "fill_value_list": ["-9999.00"],
        #"save_csv": True   
    }

    # Convert the dictionary to a json string
    json_str = json.dumps(json_dict)

    esmr = esmr.from_json(json_str)

    if only_create_remap_nc:
        # update the status of easymore, so the GIS tasks will be skipped in following calculation
        #esmr.remap_csv = ""
        # Name of column id in shp file
        #esmr.target_shp_ID = config["easymore"]["catchment_shp_hru_id_field"]
        esmr.only_create_remap_nc = True
        # execute EASYMORE
        print('Creating remap nc')
        esmr.nc_remapper()
        print('Creating remap nc complete')
        # execute EASYMORE
    else:
        #esmr.remap_csv = remap_file
        print('Starting remapping')
        esmr.remap_nc = remap_file
        esmr.nc_remapper()
        rename_easymore_output(esmr)
