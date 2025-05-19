"""General-purpose functions used in other modules.

# Functions
* `typed_dataframe` -- applies types to DataFrame columns (handles both
categorical and non-categorical types).
* `output_filename` -- builds filename for single-trial outputs.
* `save_plot` -- saves current plot, making parent folders if needed.
* `unique` -- if input list or DataFrame series only contains one unique
value, returns it. Otherwise, raises AssertionError.

## Importing LabChart recordings
* `read_record` -- extracts data from a LabChart recording segment.
* `trim_1darray` -- removes some number of values from start and end of
1D NumPy array.

## Separating trials and epochs
* `trigger_value` -- determines value of trigger peaks.
* `triggers` -- finds onset index of triggers in input.
* `separate_sweep_phases` -- splits trigger indexes by experimental
phase.

## Detecting and sorting spike responses
* `detect_spikes` -- splits input into epochs and finds peaks above
threshold.

## Saving and loading JSON files
* `convert_to_json_dict` -- converts input to JSON-compatible format.
* `confirm_save` -- checks for user confirmation before overwriting
existing files.
* `get_json_filename` -- gets list of saved JSON files matching
specified type.
* `load_spikes_trials` -- loads retrieved list of saved spike files.
* `tabulate_spikes` -- creates DataFrame summary of saved spike data.
"""

import adi
import json
import math
import numpy as np
from numpy.typing import NDArray
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
from os import walk
from scipy import signal
import re
import typing

from . import classes
from . import constants


current_version = classes.VersionNumber(constants.VERSION)


def typed_dataframe(
        df: pd.DataFrame,
        noncategorical_types_dict: dict[str, str],
        categorical_types_dict: dict[str, list]
) -> pd.DataFrame:
    """Applies dtypes to DataFrame columns.
    
    This function is only useful if at least one dtype is categorical.
    Otherwise, simply use `df.astype(noncategorical_types_dict)`.
    """
    # Define categorical dtypes using a dictionary:
    categorical_types = {
        column: pd.CategoricalDtype(categories, False) for column, categories
        in categorical_types_dict.items()
    }
    # Concatenate categorical and non-categorical dtype dictionaries:
    column_types = (noncategorical_types_dict | categorical_types)
    # Return the output DataFrame and its dtypes:
    return df.astype(column_types)


def output_filename(
        save_type: str,
        animal_id: str,
        position: int,
        repetition: int,
        recording_segment: int,
        test_code: str
) -> str:
    """Builds standardised filename for single-trial outputs.
    
    Currently used for all `main_part1` outputs (clusters, frs, epochs,
    and spikes).
    """
    unit_id = f'{animal_id}-{position}'
    recording_id = f'[{repetition}-{recording_segment}]'
    return '_'.join([
        save_type,
        unit_id,
        recording_id,
        test_code,
        str(current_version)
    ])


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


################### *IMPORTING LABCHART RECORDINGS* ####################

def read_record(
        data: adi.read.File,
        record_number: int
) -> dict[str, typing.Any]:
    """Reads a recording segment.
    
    A recording segment is created whenever recording is started within
    a `.adicht` file. Returns a dictionary containing some attributes
    for a `Recording` object. The `animal_id`, `position`, and `test`
    attributes cannot be determined from the recording.

    # Error Handling
    Raises `AssertionError` if either of two conditions occurs:
    1. Comments have not been stored in sequential order within the
    recording segment. Handling of this case has not been implemented,
    but a solution is suggested within the function body should this
    issue arise.
    2. Artefact trimming causes a comment to fall outside the usable
    portion of the recording. Again, handling of this case has not been
    implemented, but more information is provided within the error
    message.

    The error message will identify which of these conditions has been
    encountered.
    """
    # Prepare variables for timing:
    # - `tick_dt` -- duration of a sample in seconds.
    # - `b` and `a` -- used for filtering.
    # - `artefact_width_samples` -- used when trimming data to remove
    #   edge artefacts from filtering and adjusting comment timings to
    #   match.
    # TODO Stretch goal -- replace `constants.NOTCH_FILTER_F0` with a
    # TODO dynamically calculated frequency, using a Fourier transform
    # TODO to determine actual peak frequency.
    record = data.records[record_number]
    tick_dt = record.tick_dt
    b, a = signal.iirnotch(
        constants.NOTCHFILT_F0,
        constants.NOTCHFILT_Q,
        fs=1/tick_dt
    )
    artefact_width_samples = int(constants.NOTCHFILT_ARTEFACT_WIDTH_S/tick_dt)

    # Extract signal and trigger data:
    # The signal is inverted and filtered during this step, and all
    # traces and marker timings are trimmed to remove edge artefacts or
    # remain in sync with the trimmed signal.
    signal_raw = data.channels[0].get_data(record_number+1)*-1
    signal_filtered = signal.filtfilt(b, a, signal_raw)
    signal_trimmed = trim_1darray(signal_filtered, artefact_width_samples)

    mechstim = data.channels[1].get_data(record_number+1)
    mechstim_trimmed = trim_1darray(mechstim, artefact_width_samples)

    elecstim = data.channels[2].get_data(record_number+1)
    elecstim_trimmed = trim_1darray(elecstim, artefact_width_samples)

    # Extract comments and adjust timings to match trimmed signals:
    #!This will fail if comments are not stored sequentially
    # Solution: make a list of comments from `record.comments` which has
    # explicitly sorted by `record.comments.time`
    markers = []
    for i, comment in enumerate(record.comments):
        start = int(comment.time/tick_dt - artefact_width_samples)
        if i < len(record.comments)-1:
            end = int((record.comments[i+1].time)/tick_dt -
                      artefact_width_samples)
            assert start < end, (
                "Comments have not been stored in sequential order. This is "
                "known to break detection of recording trials. If this error "
                "is ever raised, `src.read_record()` must be modified to be "
                "robust to out-of-order comment storage. A suggested solution "
                "is provided within the function.\n"
                f"    Comment: {comment.text} ({i})\n"
                f"    Start: {start}\n"
                f"    End: {end}"
            )
        else:
            end =  int(len(signal_trimmed)/tick_dt - artefact_width_samples)
            assert start < end, (
                "Artefact trimming has resulted in this comment falling "
                "outside the usable portion of the recording.\n"
                f"    Comment: {comment.text} ({i})\n"
                f"    Start: {start}\n"
                f"    End: {end}\n"
                "Consult the file to decide if this recording segment should "
                "be used, and skip it if it is not useful. However, if it is "
                "useful, either reduce `artefact_width_samples` or modify "
                "this function to allow manual skipping of specific comments. "
                "Another solution would be to modify the function such that "
                "the two problems raise different error types, and handle "
                "this problem by automatically skipping any recording "
                "segments which raise it."
            )

        # Add marker to list:
        markers.append(classes.recordings.Marker(comment.text, start, end))

    # Print the duration of the record for inspection purposes:
    signal_length_s = len(signal_trimmed)*tick_dt
    print(f"[{record_number}]: {signal_length_s} s")

    # Return a `Recording` object, populating fields which are known and
    # initialising others as `None`:
    return {
        'tick_dt': tick_dt,
        'signal_data': signal_trimmed,
        'mech_triggers': mechstim_trimmed,
        'elec_triggers': elecstim_trimmed,
        'markers': markers
    }


def trim_1darray(trace: np.ndarray, trim_width: int) -> np.ndarray:
    """Removes `trim_width` values from each end of 1D NumPy array."""
    trimmed_start = trim_width
    trimmed_end = len(trace)-trim_width
    return trace[trimmed_start:trimmed_end]


#################### *SEPARATING TRIALS AND EPOCHS* ####################

def trigger_value(
        trigger_data: np.ndarray,
        correction_factor: float = 1.0
) -> float:
    """Returns value of trigger peaks in `trigger_data`."""
    # Previously, this looked for the minimum value within the trigger
    # data (i.e. the peak discharge amplitude) and corrected it to peak
    # amplitude. I believe that this was done to solve an imaginary
    # issue so have changed it to simply look for the peak amplitude
    # directly. However, if issues arise, this is a potential cause.
    return round(max(trigger_data)/correction_factor, 1)


def triggers(
        trigger_data: np.ndarray,
        trigger_threshold: int | float
) -> list[int]:
    """Returns onset index of all triggers in `trigger_data` as list."""
    # Convert `trigger_data` into list of binary values according to
    # whether each item exceeds `trigger_threshold`:
    triggers_binary = [x>trigger_threshold for x in trigger_data]
    # Convolve binary list to detect left edges:
    # Note that `np.convolve()` flips smaller array before performing
    # convolution.
    triggers_edges = np.convolve(triggers_binary, [1,-1], 'same')
    # Correct index of left edge (i.e. time in samples) for known delay
    # between trigger and recording channels:
    return [i-constants.TRIGGER_DELAY_SAMPLES for i,x in
            enumerate(triggers_edges) if x==1]


def separate_sweep_phases(
        test: str,
        mech_val: float,
        elec_val: float,
        triggers_mech: list[int],
        triggers_elec: list[int]
) -> classes.epochs.TriggersTrial:
    """Generates `TriggersTrial` object from paired trigger data.
    
    # Error Handling
    * Raises `KeyError` if test type is neither `'frequency'` nor
    `'amplitude'`.
    * Raises `ValueError` if neither stimulation value is 0, as that
    should not occur during either sweep.
    """
    # Determine the separation points between experimental phases:
    if mech_val == 0:
        start_itlv = triggers_mech[0]
        start_rcvr = triggers_mech[
            constants.INTERLEAVED_EPOCHS_FREQUENCY *
            constants.INTERLEAVED_DURATION_SECONDS
        ]
        if elec_val == 0:
            test_stim = 'control'
            stim_value = 0
        else:
            test_stim = 'electrical'
            stim_value = elec_val
    elif elec_val == 0:
        test_stim = 'mechanical'
        stim_value = mech_val
        if test == 'frequency':
            first_mech_interleaved = math.floor(
                mech_val*constants.SHORT_CONDITIONING_DURATION_SECONDS
            )
        elif test == 'amplitude':
            first_mech_interleaved = math.floor(
                constants.AMPLITUDE_SWEEP_CONDITIONING_FREQUENCY *
                constants.SHORT_CONDITIONING_DURATION_SECONDS
            )
        else:
            raise KeyError(f"Test type not recognised ({test}).")
        start_itlv = triggers_mech[first_mech_interleaved]
        start_rcvr = triggers_mech[
            first_mech_interleaved +
            constants.INTERLEAVED_EPOCHS_FREQUENCY *
            constants.INTERLEAVED_DURATION_SECONDS
        ]
    else:
        # Raise error if neither stimulation value is 0:
        raise ValueError(
            "Neither stimulation value is 0, suggesting that"
            "stimulation was performed incorrectly."
        )
    # Collect triggers into lists according to phase:
    triggers_mech_cond = [x for x in triggers_mech if x<start_itlv]
    triggers_elec_cond = [x for x in triggers_elec if x<start_itlv]
    triggers_mech_itlv = [x for x in triggers_mech if start_itlv<=x<start_rcvr]
    triggers_elec_itlv = [x for x in triggers_elec if start_itlv<=x<start_rcvr]
    triggers_mech_rcvr = [x for x in triggers_mech if x>=start_rcvr]
    triggers_elec_rcvr = [x for x in triggers_elec if x>=start_rcvr]
    # Create and return a `TriggerSet` object:
    return classes.epochs.TriggersTrial(
        test,
        test_stim,
        stim_value,
        triggers_mech_cond,
        triggers_elec_cond,
        triggers_mech_itlv,
        triggers_elec_itlv,
        triggers_mech_rcvr,
        triggers_elec_rcvr
    )


############### *DETECTING AND SORTING SPIKE RESPONSES* ################

def detect_spikes(
        signal_data: NDArray[np.floating],
        triggers: list[int],
        epoch_timing_ms: tuple[int | float, int | float],
        tick_dt: float,
        threshold: int | float,
        phase: str,
        stim_type: str
) -> tuple[
    list[classes.spikes.Spike],
    list[classes.epochs.DataEpoch]
]:
    """Detects peaks above `threshold` in `signal_data`.
    
    Returns a list of `Spike` objects containing data about each
    detected peak, and a list of `DataEpoch` objects which split the
    signal into epochs for plotting.

    # Arguments
    * `signal_data` -- signal in which to detect peaks
    * `triggers` -- list of indices indicating when each stimulation
    begins
    * `epoch_timing_ms` -- tuple of two numeric values describing the
    timing window during which spikes may occur. The first value of the
    tuple represents start time and the second value stop time (both in
    milliseconds after each stimulus onset).
    * `tick_dt` -- number of samples per second in `signal_data`.
    * `threshold` -- threshold for peak detection.
    * `phase` -- experimental phase during which `signal_data` was
    recorded.
    * `stim_type` -- type of stimulation being delivered for each epoch.
    Only one type of stimulation should be provided at a time; this
    argument does not accept tuples or lists.
    """
    # Convert `tick_dt` into milliseconds:
    tick_dt_ms = tick_dt*constants.MILLISECONDS_PER_SECOND
    # Define the epoch start and end time in samples, relative to each trigger:
    epoch_start = int(epoch_timing_ms[0]/tick_dt_ms)
    epoch_end = int(epoch_timing_ms[1]/tick_dt_ms)
    # Create empty output lists:
    spikes = []
    data_epochs = []
    for epoch_number, trigger in enumerate(triggers):
        # Create `DataEpoch` object describing current epoch and add to
        # `data_epochs` list:
        trace = signal_data[trigger+epoch_start:trigger+epoch_end]
        epoch = classes.epochs.DataEpoch(
            trace,
            epoch_timing_ms[0],
            tick_dt_ms,
            phase,
            stim_type
        )
        data_epochs.append(epoch)
        # Detect peaks greater than `threshold` within current epoch and
        # add to `spikes` list:
        peaks, properties = signal.find_peaks(trace, threshold)
        for i, peak in enumerate(peaks):
            spikes.append(classes.spikes.Spike(
                epoch_number,
                peak*tick_dt_ms + epoch_timing_ms[0],
                properties['peak_heights'][i]
            ))
    # Return populated output lists:
    return (spikes, data_epochs)


################### *SAVING AND LOADING JSON FILES* ####################

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
    if isinstance(obj, dict):
        return {key: convert_to_json_dict(value) for key, value in obj.items()}
    elif isinstance(obj, (tuple, list, np.ndarray)):
        return [convert_to_json_dict(x) for x in obj]
    elif isinstance(obj, range):
        return {
            'start': obj.start,
            'stop': obj.stop,
            'step': obj.step
        }
    elif hasattr(obj, '__dict__'):
        return convert_to_json_dict(vars(obj))
    else:
        return obj


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


def get_json_filenames(json_type: str) -> tuple[list[str], list[str]]:
    """Gets list of compatible JSON files matching `json_type`.
    
    **Currently does not check if multiple files relating to the same
    trial are present.** User should manually verify this is not true.
    """
    try:
        json_path = (constants.SAVE_PATHS['json_root'] +
                     constants.SAVE_PATHS[json_type])
    except KeyError:
        json_path = (constants.SAVE_PATHS['json_root'] +
                     f'{json_type}\\')
    json_filenames = [
        f for f in next(walk(json_path), (None, None, []))[2]
    ]
    incompatible_files = []
    for filename in json_filenames:
        file_version = classes.VersionNumber(filename)
        (compatible, message) = file_version.iscompatible(current_version)
        print(f"{filename}: {message.capitalize()}")
        if not compatible:
            incompatible_files.append(filename)
    for filename in incompatible_files:
        json_filenames.remove(filename)
    return (json_filenames, incompatible_files)


def load_spikes_trials(
        json_spike_filenames: list[str]
) -> list[list[classes.spikes.SpikesTrial]]:
    """Loads list of JSON spike files.

    Converts each loaded file to `SpikesTrial` object and returns all
    `SpikesTrial` objects as list."""
    all_trials = []
    spike_filepath = (constants.SAVE_PATHS['json_root'] +
                      constants.SAVE_PATHS['spikes'])
    for spike_filename in json_spike_filenames:
        # Construct path from which JSON file is to be read:
        save_file = spike_filepath + spike_filename
        # Read JSON file and convert to list of `SpikesTrial` objects:
        spikes_dict = json.load(open(save_file))
        assert spikes_dict['json_type'] == 'spikes'
        spikes = [classes.spikes.SpikesTrial.from_dict(dictionary)
                  for dictionary in spikes_dict['spikes']]
        all_trials.append(spikes)
    return all_trials


def tabulate_spikes(
        spikes: list[list[classes.spikes.SpikesTrial]]
) -> pd.DataFrame:
    """Loads or builds DataFrame summary of saved spike data."""
    spikes_df = pd.DataFrame()
    for list in spikes:
        for trial in list:
            spikes_df = pd.concat(
                [spikes_df, trial.to_df()],
                ignore_index=True
            )
    categories = {
        'Test': constants.TEST_CODES.values(),
        'Test Stimulus': constants.STIMULATION_TYPES,
        'Phase': constants.EXPERIMENTAL_PHASES,
        'Epoch Stimulus': constants.STIMULATION_TYPES
    }
    return typed_dataframe(
        spikes_df,
        constants.ALL_SPIKES_TYPES,
        categories
    )
