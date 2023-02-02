def argument_checker(config: dict, expected_keys: list):
    # Behavior:
    # If extra parameter is found, warn user but continue program
    # If parameter missing, raise exception

    config_values = set(config.keys())
    expected_set = set(expected_keys)

    Extra_parameters = config_values - expected_set
    Missing_parameters = expected_set - config_values

    if Extra_parameters != set() and Missing_parameters != set():
        raise Exception(
            f"Extra parameter {Extra_parameters} and missing parameter {Missing_parameters} detected."
        )

    elif Extra_parameters != set():
        print(f"Warning: Unused parameter {Extra_parameters}")

    elif Missing_parameters != set():
        raise Exception(f"Missing parameters {Missing_parameters}")


if __name__ == "__main__":

    test_dict = {
        "ABC": "EFG",
        "HI": {"JK": "HK"},
        "volTage": 10,
        "A": {"B": {"C": {"D": 7}}},
        "78": 12,
        "12A": "abc",
    }

    argument_checker(test_dict, ["ABC"])