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

from collections import abc
import json
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import numpy as np
import re
import typing

from . import classes
from . import constants


small_float = np.nextafter(0, 1)
current_version = classes.VersionNumber(constants.core.VERSION)

# Type hints
T = typing.TypeVar('T')
S = typing.TypeVar('S')
CanDoMathsT = float | np.typing.NDArray[np.floating]


# Functions
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
    try:
        # Try to create a new file to save to at the target path:
        json.dump(output_dict, open(save_file, 'x'), indent=4)
        saved = 1
    except FileNotFoundError:
        # Make directories if they don't exist
        pathlib.Path(file_path).mkdir(parents=True, exist_ok=True)
        confirm_save(file_path, filename, output_dict, force_overwrite)
        return
    except FileExistsError:
        if not force_overwrite:
            # If file already exists at target path, ask user for manual
            # confirmation before overwriting:
            confirm_overwrite = input(
                "The file you are trying to save to already exists. Do you "
                "want to overwrite it? \n"
                "Enter [Y/y] to overwrite, enter anything else or press "
                "[Escape] to skip."
            )
            manual_overwrite = bool(re.search(r"^(?i:y)$", confirm_overwrite))
        # Regex checks if user input contains only `Y` or `y`
        if force_overwrite or manual_overwrite: # type: ignore
            json.dump(output_dict, open(save_file, 'w'), indent=4)
            saved = 2
        else:
            saved = 0
    # Notify user whether file was saved or not:
    if saved == 1:
        print(f"SUCCESS: `{filename}` saved to `{file_path}`")
    elif saved == 2:
        print(f"SUCCESS: `{filename}` overwritten at `{file_path}`")
    else:
        print(f"FAILURE: `{filename}` skipped")
    return


def convert_to_json_dict(obj: typing.Any) -> typing.Any:
    """Converts `obj` to a format compatible with the JSON decoder.
    
    Searches recursively through `dicts`, `lists`, and `tuples` to
    perform the following conversions:
    * Items with the `__dict__` attribute are converted to `dict` using
    `vars()` (once converted, any such items will also be searched).
    * `range` objects are converted to `dict` with keys `start`, `stop`,
    and `step`.
    * 1-dimensional `ndarray` objects are converted to `list` and
    multi-dimensional `ndarray` objects are converted to nested `list`s.
    Other NumPy types are not converted - note that NumPy `float` is an
    instance of `float` and is compatible with the JSON decoder, but
    NumPy `int` is not an instance of `int` and is therefore
    incompatible with the JSON decoder.

    Other types which are not compatible with the JSON decoder are not
    converted, nor will they raise an error.
    """
    try:
        if isinstance(obj, range):
            return {
                'start': obj.start,
                'stop': obj.stop,
                'step': obj.step
            }
        elif isinstance(obj, (list, np.ndarray)):
            return [convert_to_json_dict(x) for x in obj]
        elif isinstance(obj, abc.Mapping):
            return {k:convert_to_json_dict(v) for k,v in obj.items()}
        elif hasattr(obj, '__dict__'):
            return convert_to_json_dict(vars(obj))
        else:
            return obj
    except RecursionError:
        print(obj)
        print(type(obj))
        raise


def detect_edges(
        bools: abc.Sequence[bool],
        direction: typing.Literal['falling', 'rising']
) -> list[int]:
    try:
        assert direction=='falling' or direction=='rising'
    except AssertionError:
        raise ValueError(f"{direction} is not a valid value for `direction`.")
    # Convolve binary list to detect left edges:
    # Note that `np.convolve()` flips smaller array before performing
    # convolution.
    convolution = np.convolve(bools, [1,-1], 'same')
    edge_value = -1 if direction == 'falling' else 1
    # Correct index of left edge (i.e. time in samples) for known delay
    # between trigger and recording channels:
    return [i for i,x in enumerate(convolution) if x==edge_value]


def remove_empty_lists(
        input_dict: abc.Mapping[T, list[S]]
) -> dict[T, list[S]]:
    return {k: v for k, v in input_dict.items() if v}


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
        target_path = (constants.core.SAVE_PATHS['plot_root'] +
                       constants.core.SAVE_PATHS[plot_type])
    except KeyError:
        # If plot_type does not have entry in `constants.SAVE_PATHS`,
        # use its name as the directory name:
        target_path = constants.core.SAVE_PATHS['plot_root'] + f'{plot_type}\\'
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


def unique(values: abc.Sequence[T] | pd.Series,) -> T:
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
    return next(iter(unique))
