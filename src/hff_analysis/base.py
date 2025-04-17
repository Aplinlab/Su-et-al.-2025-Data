"""#TODO Write module docstring.
"""

import adi
import json
import math
import numpy as np
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
from os import walk
from scipy import signal
import re

from . import classes
from . import constants


# Regex pattern for parsing filenames:
# It is compiled here once, rather than within a function, to reduce
# unnecessary computations.

adicht_regex_pattern = re.compile(constants.ADICHT_FILENAME_REGEX)


##################### *AUXILIARY HELPER FUNCTIONS* #####################

def typed_dataframe(
        df: pd.DataFrame,
        noncategorical_types_dict: dict[str, str],
        categorical_types_dict: dict[str, list]
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
    return df.astype(column_types)
    

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
    m = adicht_regex_pattern.search(filename)
    try:
        name = m.group('name')
        extension = m.group('extension') if m.group('extension') else '.adicht'
        animal_id = m.group('id')
        position = int(m.group('position'))
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
        rf'{constants.RAW_DATA_PATH}{animal_id}\{name}{extension}'
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
                constants.AMPLITUDE_SWEEP_CONDITIONING_FREQUENCY*
                constants.SHORT_CONDITIONING_DURATION_SECONDS
            )
        start_itlv = triggers_mech[first_mech_interleaved]
        start_rcvr = triggers_mech[
            first_mech_interleaved+
            constants.INTERLEAVED_EPOCHS_FREQUENCY *
            constants.INTERLEAVED_DURATION_SECONDS
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
    return spikes, data_epochs


def filter_spikes(
        trials: list[classes.spikes.SpikesTrial],
        criteria: dict[str, dict[str, int | float | None]],
        exclude_frequencies: list[int | float],
        exclude_amplitudes: list[int | float]
) -> list[classes.spikes.SpikesTrial]:
    try:
        spike_criteria_mech = classes.recordings.SpikeCriteria(
            **criteria['mechanical']
        )
        spike_criteria_elec = classes.recordings.SpikeCriteria(
            **criteria['electrical']
        )
    # TODO Write error messages for both exceptions:
    except KeyError:
        # `criteria` missing `'mechanical'` and/or `'electrical'` key(s)
        raise
    except TypeError:
        # `'mechanical'` and/or `'electrical'` don't match SpikeCriteria
        raise
    filtered_trials = []
    for trial in trials:
        if (trial.test_frequency in exclude_frequencies or
            trial.test_amplitude in exclude_amplitudes):
            continue
        filtered_phases = []
        for phase in [trial.conditioning, trial.interleaved, trial.recovery]:
            filtered_spikes_mech = []
            filtered_spikes_elec = []
            for spike in phase.mechanical:
                latency_min = (
                    spike.time_ms >= spike_criteria_mech.latency_min_ms if
                    spike_criteria_mech.latency_min_ms is not None else True
                )
                latency_max = (
                    spike.time_ms < spike_criteria_mech.latency_max_ms if
                    spike_criteria_mech.latency_max_ms is not None else True
                )
                size_min = (
                    spike.size_uV >= spike_criteria_mech.size_min_uV if
                    spike_criteria_mech.size_min_uV is not None else True
                )
                size_max = (
                    spike.size_uV < spike_criteria_mech.size_max_uV if
                    spike_criteria_mech.size_max_uV is not None else True
                )
                if latency_min and latency_max and size_min and size_max:
                    filtered_spikes_mech.append(spike)
            for spike in phase.electrical:
                latency_min = (
                    spike.time_ms >= spike_criteria_elec.latency_min_ms if
                    spike_criteria_elec.latency_min_ms is not None else True
                )
                latency_max = (
                    spike.time_ms < spike_criteria_elec.latency_max_ms if
                    spike_criteria_elec.latency_max_ms is not None else True
                )
                size_min = (
                    spike.size_uV >= spike_criteria_elec.size_min_uV if
                    spike_criteria_elec.size_min_uV is not None else True
                )
                size_max = (
                    spike.size_uV < spike_criteria_elec.size_max_uV if
                    spike_criteria_elec.size_max_uV is not None else True
                )
                if latency_min and latency_max and size_min and size_max:
                    filtered_spikes_elec.append(spike)
            filtered_phases.append(classes.spikes.SpikesPhase(
                phase.epochs_mech,
                phase.epochs_elec,
                filtered_spikes_mech,
                filtered_spikes_elec
            ))
        filtered_trials.append(classes.spikes.SpikesTrial(
            trial.animal_id,
            trial.position,
            trial.test,
            trial.test_stim,
            trial.test_frequency,
            trial.test_amplitude,
            *filtered_phases
        ))
    return filtered_trials


def save_plot(
        plot_type: str,
        target_name: str
) -> None:
    # TODO write docstring
    target_path = (f'{constants.PLOT_PATH}{plot_type}\\')
    try:
        plt.savefig(f'{target_path}{target_name}.pdf')
    except FileNotFoundError:
        pathlib.Path(target_path).mkdir(
            parents=True,
            exist_ok=True
        )
        plt.savefig(f'{target_path}{target_name}.pdf')
    

def plot_clusters(
        spikes: list[classes.spikes.SpikesTrial],
        epochs: list[classes.epochs.EpochsTrial],
        repetition: int,
        recording_segment: int,
        exclude_frequencies: list[int | float] = [],
        exclude_amplitudes: list[int | float] = [],
        save_figures: bool = False,
        filename: str | None = None
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
    # `figure(0)` displays mechanical stimulation epochs
    # `figure(1)` displays electrical stimulation epochs
    # `zorder` has been set so that points appear above traces

    # Plot detected peaks as scatterplots:
    plt.figure(figsize=constants.FIGSIZE)
    max_voltage = 0
    for spikes_trial in spikes:
        phase_stims_dict = {
            1: [
                spikes_trial.conditioning.mechanical,
                spikes_trial.interleaved.mechanical,
                spikes_trial.recovery.mechanical
            ],
            2: [
                spikes_trial.conditioning.electrical,
                spikes_trial.interleaved.electrical,
                spikes_trial.recovery.electrical
            ]
        }
        for (plot_id, phase_stims) in phase_stims_dict.items():
            plt.subplot(2, 1, plot_id)
            for phase_stim in phase_stims:
                plt.scatter(
                    [x.time_ms for x in phase_stim],
                    [x.size_uV for x in phase_stim],
                    constants.CLUSTER_POINT_SIZE,
                    constants.CLUSTER_COLOURS['peaks'],
                    zorder=1
                )
                try:
                    max_voltage = max(max_voltage, max(x.size_uV for x in phase_stim))
                except ValueError:
                    pass
    # Plot traces by epoch:
    for epochs_trial in epochs:
        if (epochs_trial.test_frequency in exclude_frequencies or
            epochs_trial.test_amplitude in exclude_amplitudes):
            continue
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
                    f"    Test Frequency: {spikes_trial.test_frequency} Hz"
                    "\n"
                    f"    Test Amplitude: {spikes_trial.test_amplitude} μA"
                )
                continue
            # Plot figure:
            plt.subplot(2, 1, plot_id+1)
            plt.plot(
                [i * epoch.tick_dt_ms + epoch.start_ms
                for i, _x in enumerate(epoch.trace)],
                epoch.trace,
                c=constants.CLUSTER_COLOURS[f'{epoch.phase}_traces'],
                linewidth=constants.CLUSTER_LINE_WIDTH,
                zorder=0
                )
    plot_title = (
        f'{spikes_trial.animal_id.upper()}-{spikes_trial.position} '
        f'[{repetition}-{recording_segment}] ('
        f'{constants.METADATA[spikes_trial.animal_id.upper()][spikes_trial.position]})'
    )
    plt.suptitle(plot_title)
    for plot_id in range(1,3):
        plt.subplot(2, 1, plot_id)
        plt.xlabel("Time (ms)")
        plt.ylabel("Signal voltage (μV)")
        ax = plt.gca()
        ax.set_ylim([
            constants.CLUSTER_YMIN,
            max_voltage*constants.CLUSTER_YMAX_SCALE
        ])
    plt.subplot(2, 1, 1)
    plt.title("Mechanical")
    plt.subplot(2, 1, 2)
    plt.title("Electrical")
    # Save figures if specified:
    if save_figures:
        assert filename is not None, (
            "When saving figures, `filename`, `recording_segment`, and "
            "`repetition` must all be provided."
        )
        save_plot(
            'clusters',
            f'{filename}-[{repetition}-{recording_segment}]'
        )
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
        input_filename: str,
        repetition: int,
        recording_index: int,
        save_type: str,
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
    # Construct target name and path for JSON file:
    file_path = (
        constants.JSON_PATHS['rootpath'] +
        constants.JSON_PATHS[save_type]['path']
    )
    output_suffix = constants.JSON_PATHS[save_type]['suffix']
    filename = (
        f'{input_filename}-[{repetition}-{recording_index}]-'
        f'{output_suffix}-{constants.VERSION}.json'
    )
    # Construct dict to save as JSON:
    output_dict = {'version': constants.VERSION, 'repetition': repetition}
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


def get_json_filenames(json_type: str) -> tuple[list[str], list[str]]:
    # TODO write docstring
    current_version = classes.VersionNumber(constants.VERSION)
    json_path = (
        constants.JSON_PATHS['rootpath'] +
        constants.JSON_PATHS[json_type]['path']
    )
    json_filenames = [
        f for f in next(walk(json_path), (None, None, []))[2]
    ]
    incompatible_files = []
    for filename in json_filenames:
        file_version = classes.VersionNumber(filename)
        if file_version.major != current_version.major:
            incompatible_files.append(filename)
            print(f"{filename}: EXCLUDE (incompatible version)")
        elif file_version.minor != current_version.minor:
            print(f"{filename}: INCLUDE (warning: minor version mismatch)")
        else:
            print(f"{filename}: INCLUDE")
    for filename in incompatible_files:
        json_filenames.remove(filename)
    return (json_filenames, incompatible_files)


def load_spikes_trials(
        json_spike_filenames: list[str]
) -> list[tuple[list[classes.spikes.SpikesTrial], int]]:
    # TODO write docstring
    all_trials = []
    for spike_filename in json_spike_filenames:
        # Construct path from which JSON file is to be read:
        save_file = (
            constants.JSON_PATHS['rootpath'] +
            constants.JSON_PATHS['spikes']['path'] +
            spike_filename
        )
        # Read JSON file and convert to list of `SpikesTrial` objects:
        spikes_dict = json.load(open(save_file))
        spikes = [classes.spikes.SpikesTrial.from_dict(dictionary)
                  for dictionary in spikes_dict['spikes']]
        all_trials.append((spikes, spikes_dict['repetition']))
    return all_trials


def tabulate_spikes(
        spikes: list[list[classes.spikes.SpikesTrial]]
) -> pd.DataFrame:
    # TODO write docstring
    spikes_df = pd.DataFrame()
    for (list, repetition) in spikes:
        for trial in list:
            spikes_df = pd.concat(
                [spikes_df, trial.to_df(repetition)],
                ignore_index=True
            )
    return typed_dataframe(
        spikes_df,
        constants.ALL_SPIKES_TYPES,
        constants.ALL_SPIKES_CATEGORIES
    )


def spikes_table(load_df_json: bool = True) -> pd.DataFrame:
    # TODO write docstring
    df_json_path = (
        constants.JSON_PATHS['rootpath'] +
        constants.JSON_PATHS['spikes_df']['path']
    )
    df_json_name = (
        constants.JSON_PATHS['spikes_df']['filename'] +
        f'-{constants.VERSION}.json'
    )

    if load_df_json:
        try:
            df_json = json.load(open(df_json_path+df_json_name))
            spikes_df = pd.DataFrame.from_dict(df_json)
            print(f"`spikes_df` loaded from file: {df_json_name}")
            print('\n'.join([unit_id for unit_id in spikes_df['Unit ID'].unique()]))
            return spikes_df
        except OSError:
            pass

    print("Building `spikes_df` from saved spikes...")
    (spikes_files, _incompatible_files) = get_json_filenames('spikes')
    all_trials = load_spikes_trials(spikes_files)
    spikes_df = tabulate_spikes(all_trials)
    try:
        spikes_df.to_json(df_json_path+df_json_name,orient='records')
    except OSError:
        pathlib.Path(df_json_path).mkdir(
            parents=True,
            exist_ok=True
        )
        spikes_df.to_json(df_json_path+df_json_name,orient='records')
    return spikes_df
