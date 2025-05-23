"""General-purpose functions used in other modules.

# Functions
* `confirm_save` -- checks for user confirmation before overwriting
existing files.
* `save_plot` -- saves current plot, making parent folders if needed.
* `unique` -- if input list or DataFrame series only contains one unique
value, returns it. Otherwise, raises AssertionError.

# Variables
* `current_version` -- module version stored as `VersionNumber` object.
"""

import json
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import re
import typing

from . import classes
from . import constants


current_version = classes.VersionNumber(constants.VERSION)


def confirm_save(
        file_path: str,
        filename: str,
        output_dict: dict[str, typing.Any],
        force_overwrite: bool = False
) -> None:
    """Asks for confirmation before overwriting existing JSON file.

    Call this function as a final step to save JSON files.
    """
    # Construct target path for JSON file:
    save_file = file_path + filename
    if force_overwrite:
        # Do not ask for confirmation before overwriting if a file
        # exists at the target path:
        json.dump(output_dict, open(save_file, 'w'), indent=4)
        saved = True
    else:
        try:
            # Try to create a new file to save to at the target path:
            json.dump(output_dict, open(save_file, 'x'), indent=4)
            saved = True
        except FileExistsError:
            # If file already exists at target path, ask user for manual
            # confirmation before overwriting:
            confirm_overwrite = input(
                "The file you are trying to save to already exists. Do you "
                "want to overwrite it? \n"
                "Enter [Y/y] to overwrite, enter anything else or press "
                "[Escape] to skip."
            )
            # Regex checks if user input contains only `Y` or `y`
            if re.search(r"^(?i:y)$", confirm_overwrite):
                json.dump(output_dict, open(save_file, 'w'), indent=4)
                saved = True
            else:
                saved = False
    # Notify user whether file was saved or not:
    if saved:
        print(f"SUCCESS: `{filename}` saved to `{file_path}`")
    else:
        print(f"FAILURE: `{filename}` skipped")
    return


def save_plot(
        plot_type: str,
        target_name: str
) -> None:
    """Saves current plot.
    
    Identifies correct directory for plot type and creates parent
    directories if necessary.
    """
    # Define target directory:
    try:
        target_path = (constants.SAVE_PATHS['plot_root'] +
                       constants.SAVE_PATHS[plot_type])
    except KeyError:
        # If plot_type does not have entry in `constants.SAVE_PATHS`,
        # use its name as the directory name:
        target_path = constants.SAVE_PATHS['plot_root'] + f'{plot_type}\\'
    # Save figure:
    try:
        plt.savefig(f'{target_path}{target_name}.pdf')
    except FileNotFoundError:
        # If target directory does not exist, create it and its parents
        # before saving figure:
        pathlib.Path(target_path).mkdir(
            parents=True,
            exist_ok=True
        )
        plt.savefig(f'{target_path}{target_name}.pdf')
    return


def unique(
        values: list | pd.Series,
) -> typing.Any:
    """Returns unique value in list or raises error if there isn't one.
    
    If `values` is a list or DataFrame series which contains only one
    unique value, that value is returned. Otherwise, AssertionError is
    raised with a message displaying list of unique values.
    """
    # Get unique values:
    if isinstance(values, pd.Series):
        unique = values.unique()
    else:
        unique = set(values)
    # Check number of unique values:
    assert len(unique) == 1, (
        f"Supplied list contains multiple values: {unique}"
    )
    # Assign first element of `unique` to `value` variable and return:
    value = None
    for value in unique:
        break
    return value
