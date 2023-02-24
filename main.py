import tomllib
import sys
import numpy as np
import time
import tomli_w
import utils

default_conf_path = "config.toml"


def main(config_path):
    with open(config_path, "rb") as f:
        config = tomllib.load(f)

    config_lower = utils.dict_2_lower(config)  # sanitize the config dict
    meas_name = config_lower["measurement"]["type"]

    measurement_init = identify_measurement_type(meas_name)
    result_dict, used_config = measurement_init(config_lower)  # Begin the measurement!

    # Extract results from dict and put in list
    resultarray = []
    keys = result_dict.keys()
    no_of_points = len(result_dict[keys[0]])  # How many rows there are

    for i in range(no_of_points):
        row = []
        for key in keys:
            item = result_dict[key][i]

            if type(item) is list:
                row += item
            else:
                row.append(item)
        resultarray.append(row)

    result_headers = []
    for key in keys:
        data = result_dict[key][0]
        if type(data) is list:
            length = f"( {len(data)})"

        result_headers.append(key + f"({length})")

    # Logic on where to save file
    save_folder = config_lower["measurement"]["save_folder"]
    if save_folder[-1:] != "/":
        save_folder = save_folder + "/"  # make sure it ends with /

    timestamp = time.strftime(rf"%Y%m%d-%H%M%S")  # get the current time
    save_file = save_folder + meas_name + "-" + timestamp

    # Save the data
    np.savetxt(save_file + ".txt", resultarray, header=" ".join(result_headers))

    # Save the config
    with open(save_file + ".toml", "wb") as f:
        tomli_w.dump(used_config, f)


def identify_measurement_type(measurement: str):
    # Matches measurement name with correct module
    match measurement:
        case "ipv":
            import measurement_type.ipv

            return measurement_type.ipv.init

        case "spectrum":
            import measurement_type.spectrum

            return measurement_type.spectrum.init

        case _:
            # TODO Change this
            raise Exception(f"No measurement of type {measurement} found.")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        config_path = default_conf_path
    else:
        # for optional system path
        config_path = sys.argv[1]

    main(config_path)
