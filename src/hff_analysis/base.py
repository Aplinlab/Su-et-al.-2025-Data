"""#TODO Write module docstring.
"""

import adi
import json
import math
import numpy as np
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
from scipy import signal
import re

from . import classes
from . import constants


# Regex pattern for parsing filenames:
# It is compiled here once, rather than within a function, to reduce
# unnecessary computations.
filename_regex_pattern = re.compile(constants.FILENAME_REGEX)


##################### *AUXILIARY HELPER FUNCTIONS* #####################

def typed_dataframe(
        columns: list[str],
        noncategorical_types_dict: dict[str, str],
        categorical_types_dict: dict[str, list],
        index: list | None = None
) -> pd.DataFrame:
    """Returns empty DataFrame with typed columns.
    
    Note that this function is only useful if at least one column dtype
    is categorical. If no column dtypes are categorical, the desired
    outcome can be achieved in a single line using:

        df = pd.DataFrame(index=index, columns=columns).astype(column_types)
    """
    # Define categorical dtypes using a dictionary:
    categorical_types = {
        column: pd.CategoricalDtype(categories, False) for column, categories
        in categorical_types_dict.items()
    }
    # Concatenate categorical and non-categorical dtype dictionaries:
    column_types = (noncategorical_types_dict | categorical_types)
    # Return the output DataFrame and its dtypes:
    if index:
        return pd.DataFrame(
            index=index,
            columns=columns
        ).astype(column_types)
    else:
        return pd.DataFrame(columns=columns).astype(column_types)
    

################### *IMPORTING LABCHART RECORDINGS* ####################

def trim_1darray(trace: np.ndarray, trim_width: int) -> np.ndarray:
    """Removes `trim_width` values from each end of a 1D NumPy array."""
    trimmed_start = trim_width
    trimmed_end = len(trace)-trim_width
    return trace[trimmed_start:trimmed_end]


def read_adicht(
        filename: str,
        data_segments: list | None = None
) -> list[classes.recordings.Recording]:
    """Reads a `.adicht` file and returns a list of `Recording` objects.
    
    # Arguments
    * `filename` -- name of the file to be read. If an extension is
    included, it will be used. If no extension is provided, `.adicht`
    will be used. Note that extensions other than `.adicht` are untested
    and may result in an error.
    * `data_segments` -- list of recording segments to read. If a list
    is not provided, all recording segments in the file will be read.

    # Error Handling
    Raises `ValueError` if any part of the expected filename pattern
    other than the extension is missing.
    Raises `KeyError` if test type is not recognised.
    """
    # Parse the input filename:
    m = filename_regex_pattern.search(filename)
    try:
        name = m.group('name')
        extension = m.group('extension') if m.group('extension') else '.adicht'
        animal_id = m.group('id')
        position = m.group('position')
        test_code = m.group('testcode')
        test = constants.TEST_CODE_CONVERSION_TABLE[test_code]
    except AttributeError as e:
        # AttributeError is raised if m is None (i.e. if regex pattern
        # didn't match)
        raise ValueError(
            f"Filename does not match the expected format ({filename})."
        ) from e
    except KeyError as e:
        # KeyError raised if `[test_code]` not found in
        # `constants.TEST_CODE_CONVERSION_TABLE`
        raise KeyError(f"Unable to match test type ({filename}).") from e
    
    # Read the specified file and break it into recording segments (each
    # time recording was started/stopped within the file), extracting
    # data from each segment separately:
    # Note: Previously, all segments were stitched together and treated
    # as a single recording. However, I believe that separating them
    # allows more flexibility in handling data (for example, if a file
    # contains multiple thresholding sweeps and only the final one is of
    # interest). It should also improve the notch filter results.
    data = adi.read_file(
        rf'{constants.RAW_DATA_FOLDER}{animal_id}\{name}{extension}'
    )
    # If a list has been specified using the `data_segments` argument,
    # it will be used; otherwise, all available segments will be read.
    if data_segments:
        records = [read_record(data, i) for i in data_segments]
    else:
        records = [read_record(data, i) for i, _x in enumerate(data.records)]
    # Populate `Recording` attributes which `read_record()` cannot:
    for record in records:
        record.animal_id = animal_id
        record.position = position
        record.test = test
    return records


def read_record(
        data: adi.read.File,
        record_number: int
) -> classes.recordings.Recording:
    """Reads a recording segment and returns a `Recording` object.
    
    A recording segment is created whenever recording is started within
    a `.adicht` file. Note that the `animal_id`, `position`, and `test`
    properties of the returned `Recording` object are initialised as
    `None` and must be populated afterwards.

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
        constants.NOTCH_FILTER_F0,
        constants.NOTCH_FILTER_Q,
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
    return classes.recordings.Recording(
        None,
        None,
        None,
        tick_dt,
        signal_trimmed,
        mechstim_trimmed,
        elecstim_trimmed,
        markers
    )


#################### *SEPARATING TRIALS AND EPOCHS* ####################

def trigger_value(
        trigger_data: np.ndarray,
        correction_factor: float = 1.0
) -> float:
    """Returns the value of trigger peaks in `trigger_data`."""
    # Previously, this looked for the minimum value within the trigger
    # data (i.e. the peak discharge amplitude) and corrected it to peak
    # amplitude. I believe that this was done to solve an imaginary
    # issue so have changed it to simply look for the peak amplitude
    # directly. However, if issues arise, this is a potential cause.
    return round(max(trigger_data)*correction_factor, 1)


def triggers(
        trigger_data: np.ndarray,
        trigger_threshold: int | float
) -> list[int]:
    """Returns onset index of all triggers in `trigger_data` as a list.

    `trigger_data` is first converted into a list of binary values
    according to whether it exceeds `trigger_threshold`. The binary list
    is then convolved with the array `[1,-1]` to detect left edges (note
    that `np.convolve` flips the smaller array before performing the
    convolution). The index (i.e. time in samples) of each left edge is
    corrected for the known latency delay between the trigger and
    recording channels, and the adjusted times are collected into a list
    which is returned.
    """
    triggers_binary = [x>trigger_threshold for x in trigger_data]
    triggers_edges = np.convolve(triggers_binary, [1,-1], 'same')
    return [i-constants.TRIGGER_DELAY_SAMPLES for i,x in
            enumerate(triggers_edges) if x==1]


def separate_sweep_phases(
        test: str,
        mech_val: float,
        elec_val: float,
        triggers_mech: np.ndarray,
        triggers_elec: np.ndarray
) -> classes.epochs.TriggersTrial:
    """Generates `TriggersTrial` object from paired trigger data.
    
    # Error Handling
    Raises `ValueError` if neither stimulation value is 0, as that
    should not occur during either sweep.
    """
    # Determine the separation points between experimental phases:
    if mech_val == 0:
        start_itlv = triggers_mech[0]
        start_rcvr = triggers_mech[
            constants.SWEEPS_INTERLEAVED_EPOCHS_EACH_STIMULUS
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
            interleaved_triggers_count = math.floor(
                mech_val*constants.SWEEPS_CONDITIONING_PHASE_SECONDS
            )
        elif test == 'amplitude':
            interleaved_triggers_count = math.floor(
                constants.AMPLITUDE_SWEEP_CONDITIONING_FREQUENCY*
                constants.SWEEPS_CONDITIONING_PHASE_SECONDS
            )
        start_itlv = triggers_mech[interleaved_triggers_count]
        start_rcvr = triggers_mech[
            interleaved_triggers_count+
            constants.SWEEPS_INTERLEAVED_EPOCHS_EACH_STIMULUS
        ]
    else:
        # Raise error if neither stimulation value is 0:
        raise ValueError(
            "Neither stimulation value is 0, suggesting that stimulation was "
            "performed incorrectly."
        )
    # Collect triggers into lists according to phase:
    triggers_mech_cond = [x for x in triggers_mech if x<start_itlv]
    triggers_elec_cond = [x for x in triggers_elec if x<start_itlv]
    triggers_mech_itlv = [x for x in triggers_mech if x>=start_itlv and
                        x<start_rcvr]
    triggers_elec_itlv = [x for x in triggers_elec if x>=start_itlv and
                        x<start_rcvr]
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
        signal_data: np.ndarray,
        triggers: list[int],
        epoch_timing_ms: tuple[int | float, int | float],
        tick_dt: float,
        threshold: int | float,
        spike_criteria: dict[str, classes.recordings.SpikeCriteria],
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
    data_epochs = []
    spikes = []
    for epoch_id, trigger in enumerate(triggers):
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
            # Check if all criteria match, skipping if criteria is None
            criteria = spike_criteria[stim_type]
            latency_ms = peak*tick_dt_ms + epoch_timing_ms[0]
            size_uV = properties['peak_heights'][i]
            latency_min = (
                latency_ms >= criteria.latency_min_ms if
                criteria.latency_min_ms is not None else True
            )
            latency_max = (
                latency_ms < criteria.latency_max_ms if
                criteria.latency_max_ms is not None else True
            )
            size_min = (
                size_uV >= criteria.size_min_uV if
                criteria.size_min_uV is not None else True
            )
            size_max = (
                size_uV < criteria.size_max_uV if
                criteria.size_max_uV is not None else True
            )
            if latency_min and latency_max and size_min and size_max:
                spikes.append(classes.spikes.Spike(
                    epoch_id,
                    latency_ms,
                    size_uV
                ))
    # Return populated output lists:
    return spikes, data_epochs


def plot_spikes(
    spikes: list[classes.spikes.SpikesTrial],
    epochs: list[classes.epochs.EpochsTrial],
    save_figures: bool = False,
    filename: str | None = None,
    recording_segment: int | None = None
) -> None:
    """Plots peaks and traces by epoch to visualise spike detection.
    
    # Error Handling
    * Prints an error message if the stimulus type of any given epoch is
    neither `'mechanical'` nor `'electrical'` but continues plotting
    other epochs. The error message contains information which may help
    to identify the probematic epoch.
    * Raises `AssertionError` if `save_figures` set to `True` but one or
    both of `filename` and `recording_segment` are not provided.
    """
    for spikes_trial in spikes:
        # Plot detected peaks as scatterplots:
        # `figure(0)` displays mechanical stimulation epochs
        # `figure(1)` displays electrical stimulation epochs
        # zorder has been set so that 
        plt.figure(0)
        plt.scatter(
            [x.time_ms for x in spikes_trial.conditioning.mechanical],
            [x.size_uV for x in spikes_trial.conditioning.mechanical],
            constants.SCATTERPLOT_POINT_SIZE,
            constants.PLOT_COLOURS['peaks'],
            zorder=1
        )
        plt.scatter(
            [x.time_ms for x in spikes_trial.interleaved.mechanical],
            [x.size_uV for x in spikes_trial.interleaved.mechanical],
            constants.SCATTERPLOT_POINT_SIZE,
            constants.PLOT_COLOURS['peaks'],
            zorder=1
        )
        plt.scatter(
            [x.time_ms for x in spikes_trial.recovery.mechanical],
            [x.size_uV for x in spikes_trial.recovery.mechanical],
            constants.SCATTERPLOT_POINT_SIZE,
            constants.PLOT_COLOURS['peaks'],
            zorder=1
        )
        plt.figure(1)
        plt.scatter(
            [x.time_ms for x in spikes_trial.conditioning.electrical],
            [x.size_uV for x in spikes_trial.conditioning.electrical],
            constants.SCATTERPLOT_POINT_SIZE,
            constants.PLOT_COLOURS['peaks'],
            zorder=1
        )
        plt.scatter(
            [x.time_ms for x in spikes_trial.interleaved.electrical],
            [x.size_uV for x in spikes_trial.interleaved.electrical],
            constants.SCATTERPLOT_POINT_SIZE,
            constants.PLOT_COLOURS['peaks'],
            zorder=1
        )
        plt.scatter(
            [x.time_ms for x in spikes_trial.recovery.electrical],
            [x.size_uV for x in spikes_trial.recovery.electrical],
            constants.SCATTERPLOT_POINT_SIZE,
            constants.PLOT_COLOURS['peaks'],
            zorder=1
        )
        # Plot traces by epoch:
        for epochs_trial in epochs:
            for epoch in epochs_trial.epochs:
                # If an epoch has an invalid stimulus type, try-except block
                # prints an error message and skips it:
                try:
                    plot_id = constants.STIMULATION_TYPES.index(epoch.stimulus)
                except ValueError:
                    print(
                        "Epoch stimulus type not recognised.\n"
                        f"    Animal ID: {spikes_trial.animal_id}\n"
                        f"    Position: {spikes_trial.position}\n"
                        f"    Test: {spikes_trial.test}\n"
                        f"    Test Stimulus: {spikes_trial.test_stim}\n"
                        f"    Test Frequency: {spikes_trial.test_freqency} Hz"
                        "\n"
                        f"    Test Amplitude: {spikes_trial.test_amplitude} μA"
                    )
                    continue
                # Plot figure:
                plt.figure(plot_id)
                plt.plot(
                    [i * epoch.tick_dt_ms + epoch.start_ms
                    for i, _x in enumerate(epoch.trace)],
                    epoch.trace,
                    c=constants.PLOT_COLOURS[f'{epoch.phase}_traces'],
                    linewidth=constants.EPOCH_LINE_WIDTH,
                    zorder=0
                )
    # Save figures if specified:
    if save_figures:
        assert filename is not None and recording_segment is not None, (
            "When saving figures, `filename` and `recording_segment` "
            "arguments must both be provided."
        )
        common_file_path = (
            f'{constants.CLUSTER_PLOTS_FOLDER}{filename}'
            f'-[{recording_segment}]-'
        )
        plt.figure(0)
        plt.xlabel("Time (ms)")
        plt.ylabel("Signal voltage (μV)")
        plt.savefig(common_file_path + 'mechanical.pdf')
        plt.figure(1)
        plt.xlabel("Time (ms)")
        plt.ylabel("Signal voltage (μV)")
        plt.savefig(common_file_path + 'electrical.pdf')
    return


################### *SAVING AND LOADING JSON FILES* ####################

def convert_to_json_dict(obj) -> any:
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


def save_to_json(
        file_path: str,
        input_filename: str,
        recording_index: int,
        output_suffix: str,
        force_overwrite: bool = False,
        **kwargs
) -> None:
    """Saves a list of objects as a JSON file.
    
    # Arguments
    * `file_path` -- path to directory where JSON file is to be saved.
    * `input_filename` -- name of LabChart file from which data was
    obtained.
    * `recording_segment` -- index of recording segment from which data
    was obtained, relative to the entire file at `input_filename`.
    entire file.
    * `force_overwrite` -- whether to overwrite an existing file at the
    target path without asking. If set to `True`, an existing file will
    be overwritten without manual confirmation. If set to `False`, user
    will be asked to confirm or reject overwrite if a file exists at the
    target path.
    * **kwargs -- data to include in the saved JSON file.
        * Any objects which have __dict__ attribute will be expanded

    # Error Handling
    Raises `TypeError` if any values in `kwargs` contain items which are
    not JSON serialisable after applying `convert_to_json_dict()`.
    """
    # Construct target filename for JSON file:
    filename = f'{input_filename}-[{recording_index}]-{output_suffix}.json'
    # Construct dict to save as JSON:
    output_dict = {}
    for key, value in kwargs.items():
        output_dict[key] = convert_to_json_dict(value)
    # Save to JSON and print success/failure message:
    try:
        confirm_save(file_path, filename, output_dict, force_overwrite)
    except FileNotFoundError:
        # Make directories if they don't exist
        pathlib.Path(file_path).mkdir(parents=True, exist_ok=True)
        confirm_save(file_path, filename, output_dict, force_overwrite)
    except TypeError:
        # Raise error if `output_dict` contains non-JSON serialisable
        # objects:
        raise
    return


def confirm_save(
        file_path: str,
        filename: str,
        output_dict: dict[str, any],
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