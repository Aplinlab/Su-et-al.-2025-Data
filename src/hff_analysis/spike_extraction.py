"""Functions used for initial processing of LabChart recordings.

Includes functions for reading LabChart recordings, signal processing,
and peak detection.

# Functions
## User-exposed functions
See package docstring for function summaries.
* `load_filereadsettings`
* `read_adicht`
* `spikes_info`

## Backend functions
* `detect_spikes` -- splits input into epochs and finds peaks above
threshold.
* `max_amp` -- calculates maximum amplitude used in amplitude sweep.
* `read_record` -- extracts data from a LabChart recording segment.
* `separate_sweep_phases` -- splits trigger indexes by experimental
phase.
* `trigger_value` -- determines value of trigger peaks.
* `triggers` -- finds onset index of triggers in input.
* `trim_1darray` -- removes some number of values from start and end of
1D NumPy array.

## Variables
* `long_pattern` - regex pattern for parsing 5 minute trial comments.
* `nineone_pattern` - regex pattern for parsing comments in 9:1
interleaved trials.
* `sweep_pattern` - regex pattern for parsing comments in frequency and
amplitude sweeps.
"""
"""Functions used for refining peak detection and isolating responses.

Includes functions for filtering spikes as well as plotting and saving
isolated responses.

# Functions
## User-exposed functions
See package docstring for function summaries.
* `filter_spikes`
* `plot_clusters`
* `save_to_json`

## Backend functions
* `convert_to_json_dict` -- converts input to JSON-compatible format.
* `output_filename` -- builds filename for single-trial outputs.
"""

import adi
from collections import abc
import itertools
import json
import math
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy import signal
import statistics
from typing import Any

from . import classes
from . import constants
from . import utils


# Regex patterns for parsing marker comments:
# Compiled once, rather than within a function, to reduce computations.
long_pattern = re.compile(constants.regex.LONGDURATION_REGEX)
nineone_pattern = re.compile(constants.regex.NINEONE_REGEX)
sweep_pattern = re.compile(constants.regex.SWEEP_REGEX)


####################### *USER-EXPOSED FUNCTIONS* #######################

def filter_trials(
        trials: abc.Iterable[classes.SpikesTrial],
        criteria: abc.Mapping[str, abc.Mapping[str, float | None]],
        exclude_frequencies: abc.Container[float],
        exclude_amplitudes: abc.Container[float],
        enforce_max_failrate: bool = False
) -> tuple[
    list[classes.SpikesTrial],
    classes.ISIResult,
    classes.SpikeCriteria,
    classes.SpikeCriteria
]:
    """Applies exclusion criteria to detected peaks.
    
    The purpose of this function is to isolate responses from a single
    unit.

    # Arguments
    * `trials` -- list of `SpikesTrial` objects containing peaks to
    filter.
    * `criteria` -- dictionary with keys `'mechanical'` and
    `'electrical'`, whose values are dictionaries themselves having keys
    for 'latency_min_ms', 'latency_max_ms', 'size_min_uV', and
    'size_max_uV'. Each of these keys has a numeric value specifying an
    exclusion criteria, or None for no exclusion.
    * `exclude_frequencies` and `exclude_amplitudes` -- if a trial has
    test frequency or amplitude (respectively) whose value is in one or
    both of these lists, the entire trial will be excluded.
    """
    try:
        spike_criteria_mech = classes.SpikeCriteria(**criteria['mechanical'])
        spike_criteria_elec = classes.SpikeCriteria(**criteria['electrical'])
    except KeyError:
        print("`criteria` missing `'mechanical'` and/or `'electrical'` key(s)")
        raise
    except TypeError:
        print("`'mechanical'` and/or `'electrical'` don't match SpikeCriteria")
        raise
    # Build new list of `SpikesTrial` objects:
    filtered_trials = []
    for trial in trials:
        # Test for trial exclusion criteria:
        if (trial.test_frequency in exclude_frequencies or
            trial.test_amplitude in exclude_amplitudes):
            continue
        # Test each experimental phase separately:
        # `filtered_phases` will contain one `SpikesPhase` object for
        # each experimental phase
        filtered_phases = {}
        rejected_spikes = {stim_type: {} for stim_type in
                           constants.experiments.STIMULATION_TYPES}
        failed_spikes = {stim_type: {} for stim_type in
                         constants.experiments.STIMULATION_TYPES}
        failed_spikes_count = {stim_type: 0 for stim_type in
                               constants.experiments.STIMULATION_TYPES}
        for attribute, spikes_phase in (
            ('spikes_cond', trial.conditioning),
            ('spikes_itlv', trial.interleaved),
            ('spikes_rcvr', trial.recovery)
        ):
            filtered_spikes = {}
            for stim_type, phase_stim, spike_criteria in (
                ('mechanical', spikes_phase.mechanical, spike_criteria_mech),
                ('electrical', spikes_phase.electrical, spike_criteria_elec)
            ):
                # Apply spike exclusion criteria:
                (
                    filtered_spikes_stim,
                    rejected_spikes_stim,
                    failed_spikes_stim,
                    failed_count_stim
                ) = filter_spikes(
                    phase_stim, # type: ignore
                    spike_criteria,
                    attribute[-4:]
                )
                filtered_spikes[stim_type] = filtered_spikes_stim
                rejected_spikes[stim_type].update(rejected_spikes_stim)
                failed_spikes[stim_type].update(failed_spikes_stim)
                failed_spikes_count[stim_type] += failed_count_stim
            filtered_phases[attribute] = classes.SpikesPhase(
                spikes_phase.epochs_mech,
                spikes_phase.epochs_elec,
                utils.remove_empty_lists(filtered_spikes['mechanical']),
                utils.remove_empty_lists(filtered_spikes['electrical'])
            )
        rejected_spikes_noempty = {
            k: utils.remove_empty_lists(v) for k,v in rejected_spikes.items()
        }
        filtered_phases['spikes_rejected'] = classes.SpikesPhase(
            len(rejected_spikes_noempty['mechanical']),
            len(rejected_spikes_noempty['electrical']),
            rejected_spikes_noempty['mechanical'],
            rejected_spikes_noempty['electrical']
        )
        filtered_phases['spikes_failed'] = classes.SpikesPhase(
            failed_spikes_count['mechanical'],
            failed_spikes_count['electrical'],
            failed_spikes['mechanical'],
            failed_spikes['electrical']
        )
        # Generate `SpikesTrial` object from filtered spikes:
        filtered_trials.append(classes.SpikesTrial(
            trial.animal_id,
            trial.position,
            trial.test,
            trial.test_stim,
            trial.test_frequency,
            trial.test_frequency_minor,
            trial.test_amplitude,
            trial.repetition,
            trial.stddev,
            trial.threshold,
            spike_criteria_mech,
            spike_criteria_elec,
            **filtered_phases
        ))
    failed_epochs_total = sum((
        sum(trial.failed.epochs_mech for trial in filtered_trials),
        sum(trial.failed.epochs_elec for trial in filtered_trials)
    ))
    epochs_total = sum((
        sum(len(trial.conditioning.mechanical) for trial in filtered_trials),
        sum(len(trial.interleaved.mechanical) for trial in filtered_trials),
        sum(len(trial.recovery.mechanical) for trial in filtered_trials),
        sum(len(trial.conditioning.electrical) for trial in filtered_trials),
        sum(len(trial.interleaved.electrical) for trial in filtered_trials),
        sum(len(trial.recovery.electrical) for trial in filtered_trials)
    ))
    failure_rate = failed_epochs_total/epochs_total if epochs_total else 0
    if (
        enforce_max_failrate and
        (failure_rate>constants.core.MAXIMUM_ISI_FAILRATE)
    ):
        spike_criteria_mech.increment_minimum_size(
            constants.core.THRESHOLD_INCREMENT
        )
        spike_criteria_elec.increment_minimum_size(
            constants.core.THRESHOLD_INCREMENT
        )
        spike_criteria = {
            'mechanical': vars(spike_criteria_mech),
            'electrical': vars(spike_criteria_elec)
        }
        return filter_trials(
            trials,
            spike_criteria,
            exclude_frequencies,
            exclude_amplitudes,
            True
        )
    isi_result = classes.ISIResult(
        failed_epochs_total,
        epochs_total,
        failure_rate
    )
    print(
        'ISI failure rate: %d / %d = %.4f' % (
            failed_epochs_total,
            epochs_total,
            failure_rate
        )
    )
    # Return complete list of filtered peaks:
    return (
        filtered_trials,
        isi_result,
        spike_criteria_mech,
        spike_criteria_elec
    )


def load_filereadsettings(
        frs_filename: str
) -> classes.FileReadSettings:
    """Loads `FileReadSettings` object from JSON file.
    
    The returned object can be used to reproduce an earlier run of
    `main_part1`.
    """
    # Construct file path:
    frs_filepath = (constants.core.SAVE_PATHS['json_root'] +
                    constants.core.SAVE_PATHS['frs'])
    save_file = frs_filepath + frs_filename
    # Add file extension if not already present:
    if save_file[-5:].lower() != '.json':
        save_file += '.json'
    # Load file as `FileReadSettings` object:
    frs = classes.FileReadSettings.from_dict(json.load(open(
        save_file
    )))
    return frs
    

def plot_clusters(
        spikes: abc.Iterable[classes.SpikesTrial],
        epochs: abc.Iterable[classes.EpochsTrial],
        repetition: int,
        recording_number: int,
        xlim: tuple[float, float],
        spike_criteria_mech: classes.SpikeCriteria,
        spike_criteria_elec: classes.SpikeCriteria,
        isi_result: classes.ISIResult,
        exclude_frequencies: abc.Container[float] = [],
        exclude_amplitudes: abc.Container[float] = [],
        save_plot: bool = False
) -> None:
    """Plots peaks and traces by epoch to visualise spike detection.
    
    **Spikes from different recordings cannot be plotted together.**
    This is primarily a conceit of the naming system for saved plots,
    and is possible to implement. However, while this functionality
    could be useful (e.g. to verify that units across tests are the
    same), it is not useful for the current analysis and therefore is
    not planned.
    
    # Error Handling
    * Prints an error message if the stimulus type of any given epoch is
    neither `'mechanical'` nor `'electrical'` but continues plotting
    other epochs. The error message contains information which may help
    to identify the probematic epoch.
    * Raises `AssertionError` if there is more than one unique value for
    any of `animal_id`, `position`, or `test` present between `spikes`
    and `epochs`.
    """
    # Top subplot displays electrical stimulation epochs and lower epoch
    # displays mechanical stimulation epochs. `zorder` is set in the
    # following order (from most obscured layer to most visible layer):
    # 0: traces
    # 1: threshold
    # 2: manual spike boundaries
    # 3: excluded spikes
    # 4: correct spikes
    # 5: incorrect spikes (ISI too small)

    # Initialise variables for later:
    plt.figure(figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.CLUSTER_FIG_HEIGHT
    ))
    max_voltage = 0
    animal_id = utils.unique(
        [spikes_trial.animal_id for spikes_trial in spikes] +
        [epochs_trial.animal_id for epochs_trial in epochs]
    )
    position = utils.unique(
        [spikes_trial.position for spikes_trial in spikes] +
        [epochs_trial.position for epochs_trial in epochs]
    )
    test = utils.unique(
        [spikes_trial.test for spikes_trial in spikes] +
        [epochs_trial.test for epochs_trial in epochs]
    )
    threshold = utils.unique(
        [spikes_trial.threshold for spikes_trial in spikes] +
        [epochs_trial.threshold for epochs_trial in epochs]
    )
    # Plot detected peaks as scatterplots:
    for spikes_trial in spikes:
        max_voltage_elec = plot_peaks(
            1,
            spikes_trial.rejected.electrical.values(),
            (
                spikes_trial.conditioning.electrical,
                spikes_trial.interleaved.electrical,
                spikes_trial.recovery.electrical
            ),
            spikes_trial.failed.electrical.values()
        )
        max_voltage_mech = plot_peaks(
            2,
            spikes_trial.rejected.mechanical.values(),
            (
                spikes_trial.conditioning.mechanical,
                spikes_trial.interleaved.mechanical,
                spikes_trial.recovery.mechanical
            ),
            spikes_trial.failed.mechanical.values()
        )
        max_voltage = max(max_voltage, max_voltage_mech, max_voltage_elec)
    # Plot traces by epoch:
    for epochs_trial in epochs:
        # Test for trial exclusion criteria:
        # `filter_spikes()` only applies these to spikes, not epochs
        if (epochs_trial.test_frequency in exclude_frequencies or
            epochs_trial.test_amplitude in exclude_amplitudes):
            continue
        overlay_epochs(epochs_trial, animal_id, position, test)
    # Plot threshold and manual criteria:
    y_max = max_voltage * constants.plots.CLUSTER_YMAX_SCALE
    plot_filter_bounds(
        1,
        threshold,
        spike_criteria_elec,
        xlim,
        (constants.plots.CLUSTER_YMIN, y_max)
    )
    plot_filter_bounds(
        2,
        threshold,
        spike_criteria_mech,
        xlim,
        (constants.plots.CLUSTER_YMIN, y_max)
    )
    # Prettify figure:
    plot_title = (
        f'{animal_id.upper()}-{position} [{repetition}-{recording_number}] (' +
        constants.experiments.ANIMAL_DATA[animal_id.upper()][position]["type"]+
        ') | %s | Failure rate: %.2f' % (test.capitalize(), isi_result.result)
    )
    for plot_id in (1, 2):
        plt.subplot(2, 1, plot_id)
        plt.xlabel("Time (ms)")
        plt.ylabel("Signal voltage (μV)")
        ax = plt.gca()
        ax.set_ylim((
            constants.plots.CLUSTER_YMIN,
            y_max
        ))
        ax.set_xlim(xlim)
    plt.subplot(2, 1, 1)
    plt.title("Electrical")
    plt.subplot(2, 1, 2)
    plt.title("Mechanical")
    plt.suptitle(plot_title)
    plt.tight_layout()
    # Save figures if specified:
    if save_plot:
        plot_type = 'clusters'
        target_name = output_filename(
            plot_type,
            animal_id,
            position,
            repetition,
            recording_number,
            test[0:4]
        )
        utils.save_plot(
            plot_type,
            target_name
        )
    return


def read_adicht(
        filename: str,
        show_feedback: bool = True
) -> list[classes.Recording]:
    """Reads a `.adicht` file and returns a list of `Recording` objects.
    
    # Arguments
    * `filename` -- name of the file to be read. If an extension is
    included, it will be used. If no extension is provided, `.adicht`
    will be used. **Extensions other than `.adicht` are untested
    and may result in an error.**
    * `data_segments` -- list of recording segments to read. If a list
    is not provided, all recording segments in the file will be read.

    # Error Handling
    Raises `ValueError` if any part of the expected filename pattern
    other than the extension is missing.
    Raises `KeyError` if test type is not recognised.
    """
    # Parse the input filename:
    filename_info = classes.FilenameInfo.from_filename(filename)
    # Read the specified file and break it into recording segments (each
    # time recording was started/stopped within the file), extracting
    # data from each segment separately:
    # Note: Previously, all segments were stitched together and treated
    # as a single recording. However, I believe that separating them
    # allows more flexibility in handling data (for example, if a file
    # contains multiple thresholding sweeps and only the final one is of
    # interest). It should also improve the notch filter results.
    data = adi.read_file(
        constants.core.RAW_DATA_PATH + filename_info.animal_id + '\\' +
        filename_info.name + filename_info.extension
    )
    # Read all available data segments:
    # The first segment must be read to calculate stimulation threshold,
    # so all segments are processed to reduce code complexity. If
    # `data_segments` has been provided, only the specified segments
    # will be returned.
    records = [read_record(
        data,
        i,
        show_feedback
    ) for i, _ in enumerate(data.records)]
    # Calculate threshold from portion of recording before first comment:
    stddev = threshold_stddev(
        records[0]['signal_data'],
        records[0]['bessel'],
        records[0]['notch'],
        records[0]['tick_dt'],
        records[0]['markers'],
        show_feedback
    )
    # Return list of `Recording` objects, filling in attributes which
    # `read_record()` cannot:
    output = [
        classes.Recording(
            filename_info.animal_id,
            filename_info.position,
            filename_info.test,
            record['tick_dt'],
            record['signal_data'],
            record['mech_triggers'],
            record['elec_triggers'],
            record['markers'],
            stddev,
            record['bessel'],
            record['notch']
        ) for record in records
    ]
    return output


def save_to_json(
        input_filename: str,
        repetition: int,
        recording_number: int,
        json_type: str,
        force_overwrite: bool = False,
        **kwargs
) -> None:
    """Saves a list of objects as a JSON file.
    
    # Arguments
    * `input_filename` -- name of LabChart file from which data was
    obtained.
    * `repetition` - which repetition of its test the file is.
    * `recording_segment` -- index of recording segment from which data
    was obtained, relative to the entire file at `input_filename`.
    * `json_type` -- category of JSON file to save (determines save path
    and included in output filename).
    * `force_overwrite` -- whether to overwrite an existing file at the
    target path without asking. If set to `True`, an existing file will
    be overwritten without manual confirmation. If set to `False`, user
    will be asked to confirm or reject overwrite if a file exists at the
    target path.
    * `**kwargs` -- data to include in the saved JSON file.

    # Error Handling
    Raises `TypeError` if any values in `kwargs` contain items which are
    not JSON serialisable after applying `convert_to_json_dict()`.
    """
    # Construct target name and path for JSON file:
    try:
        file_path = (constants.core.SAVE_PATHS['json_root'] +
                     constants.core.SAVE_PATHS[json_type])
    except KeyError:
        file_path = (constants.core.SAVE_PATHS['json_root'] +
                     f'{json_type}\\')
    input_info = classes.FilenameInfo.from_filename(input_filename)
    filename = output_filename(
        json_type,
        input_info.animal_id,
        input_info.position,
        repetition,
        recording_number,
        input_info.test[0:4]
    ) + '.json'
    # Construct dict to save as JSON:
    output_dict = {
        'json_type': json_type,
        'version': str(utils.current_version),
        'repetition': repetition,
        'recording_number': recording_number
    }
    for key, value in kwargs.items():
        output_dict[key] = utils.convert_to_json_dict(value)
    # Save to JSON and print success/failure message:
    try:
        utils.confirm_save(file_path, filename, output_dict, force_overwrite)
    except TypeError:
        # Raise error if `output_dict` contains non-JSON serialisable
        # objects:
        raise
    return


def spikes_info(
        recording: classes.Recording,
        repetition: int,
        epoch_timing_ms: tuple[float, float],
        skip_superfast: bool = True
) -> tuple[list[classes.SpikesTrial], list[classes.EpochsTrial]]:
    # Extract triggers as list of all triggers and individual lists with marker metadata
    # Filter using list of all triggers
    # Extract spikes using lists with metadata
    """Returns spikes and traces by epoch from extracted data.
    
    Applies signal processing (frequency filters) to extracted data,
    separates it into epochs, and detects peaks above `threshold_uV`.
    
    # Arguments
    * `recording` -- `Recording` object to be analysed.
    * `epoch_timing_ms` -- tuple of two numeric values describing timing
    window during which spikes may occur. The first value represents
    start time and the second value stop time in milliseconds after
    stimulus onset.
    * `threshold_uV` -- threshold above which spikes should be detected,
    in microvolts.

    # Error Handling
    * Upon encountering a comment which does not match the expected
    format for frequency sweeps, prints an error message containing the
    comment text and moves onto the next comment.
    * Raises `AssertionError` if the number of epochs extracted from any
    trial does not match the number of triggers present in that trial.
    """
    # Convert `tick_dt` into milliseconds:
    tick_dt_ms = recording.tick_dt * constants.core.MILLISECONDS_PER_SECOND
    # Calculate spike and trigger thresholds:
    threshold = recording.stddev * constants.core.SPIKE_NOISE_FLOOR
    trigger_threshold_values = trigger_thresholds(
        recording
    )
    # Create empty lists to populate and return later:
    triggers_by_stim = {'mechanical': [], 'electrical': []}
    triggers_by_marker = []
    output_spikes = []
    output_epochs = []
    for marker in recording.markers:
        try:
            triggers_info = trigger_info(
                marker,
                recording,
                *trigger_threshold_values
            )
        except (AssertionError, AttributeError):
            continue
        triggers_by_stim['mechanical'] += [x+marker.start_sample for x in
                                           triggers_info.mech_triggers]
        triggers_by_stim['electrical'] += [x+marker.start_sample for x in
                                           triggers_info.elec_triggers]
        # Sort triggers by experimental phase using frequencies:
        triggers_by_phase = separate_sweep_phases(
            recording.test,
            triggers_info
        )
        # test_frequency, test_frequency_minor, test_amplitude = test_variables(
        #     recording.test,
        #     triggers_by_phase
        # )
        # Collect variables used when defining both `SpikesTrial`
        # and `EpochsTrial` objects to reduce repetition:
        effective_frequency = (triggers_by_phase.frequency if
                               recording.test!='nine-one' else
                               triggers_by_phase.frequency +
                               triggers_by_phase.minor_frequency)
        # Skip if frequency is too fast
        if (skip_superfast and (1/effective_frequency)<
            (epoch_timing_ms[1]/constants.core.MILLISECONDS_PER_SECOND)):
            continue
        else:
            triggers_by_marker.append({
                'animal_id': recording.animal_id,
                'position': recording.position,
                'test': recording.test,
                'test_stim': triggers_by_phase.test_stim,
                'test_frequency': triggers_by_phase.frequency,
                'test_frequency_minor': triggers_by_phase.minor_frequency,
                'test_amplitude': triggers_by_phase.amplitude,
                'repetition': repetition,
                'marker_text': marker.comment,
                'start_sample': marker.start_sample,
                'end_sample': marker.end_sample,
                'triggers': triggers_by_phase
            })
    signal_data = recording.signal_data
    # for stim_type in constants.experiments.STIMULATION_TYPES:
    #     signal_data = blank_artefacts(
    #         signal_data,
    #         triggers_by_stim[stim_type],
    #         constants.experiments.PULSEWIDTH_MS[stim_type],
    #         tick_dt_ms
    #     )
    # signal_filtered = filter_signal(
    #     signal_data,
    #     recording.bessel,
    #     recording.notch
    # )
    signal_filtered = process_signal(
        signal_data,
        triggers_by_stim['mechanical'],
        triggers_by_stim['electrical'],
        tick_dt_ms,
        recording.bessel,
        recording.notch
    )
    for trigger_dict in triggers_by_marker:
        # Constrain signal_data to the current trial:
        signal_data_trial = signal_filtered[
            trigger_dict['start_sample']:trigger_dict['end_sample']
        ]
        # Detect spikes within each epoch and extract traces:
        spike_detection_results = {
            phase: {} for phase in constants.experiments.PHASES
        }
        epochs = []
        for phase in spike_detection_results.keys():
            for stim_type in constants.experiments.STIMULATION_TYPES:
                detected_spikes, new_epochs = detect_spikes(
                    signal_data_trial,
                    trigger_dict['triggers'].triggers[phase][stim_type],
                    epoch_timing_ms,
                    tick_dt_ms,
                    threshold,
                    phase,
                    stim_type
                )
                spike_detection_results[phase][stim_type] = detected_spikes
                epochs += new_epochs
        # Define `SpikesTrial` object:
        # It is important to count the number of triggers rather
        # supplying a precalculated value, since trimming of the
        # recording in an earlier step (to remove artefacts caused
        # by filtering) results in some markers being lost from the
        # recovery phase of the final trial in each recording.
        spikes_phases = {
            phase: classes.SpikesPhase(
                len(trigger_dict['triggers'].triggers[phase]['mechanical']),
                len(trigger_dict['triggers'].triggers[phase]['electrical']),
                spike_detection_results[phase]['mechanical'],
                spike_detection_results[phase]['electrical']
            ) for phase in constants.experiments.PHASES
        }
        empty_criteria = classes.SpikeCriteria(
            None,
            None,
            None,
            None,
        )
        empty_phase = classes.SpikesPhase(
            0,
            0,
            {},
            {}
        )
        output_spikes.append(classes.SpikesTrial(
            trigger_dict['animal_id'],
            trigger_dict['position'],
            trigger_dict['test'],
            trigger_dict['test_stim'],
            trigger_dict['test_frequency'],
            trigger_dict['test_frequency_minor'],
            trigger_dict['test_amplitude'],
            trigger_dict['repetition'],
            recording.stddev,
            threshold,
            empty_criteria,
            empty_criteria,
            spikes_phases['conditioning'],
            spikes_phases['interleaved'],
            spikes_phases['recovery'],
            empty_phase,
            empty_phase
        ))
        # Check that number of epochs matches number of triggers:
        assert len(epochs) == sum(
            len(trigger_dict['triggers'].triggers[phase]['mechanical'])+
            len(trigger_dict['triggers'].triggers[phase]['electrical'])
            for phase in constants.experiments.PHASES
        ), (
            "The number of epochs does not match the number of triggers "
            f"({trigger_dict['marker_text'].comment}). There may be an "
            "issue with the code."
        )
        output_epochs.append(classes.EpochsTrial(
            trigger_dict['animal_id'],
            trigger_dict['position'],
            trigger_dict['test'],
            trigger_dict['test_stim'],
            trigger_dict['test_frequency'],
            trigger_dict['test_frequency_minor'],
            trigger_dict['test_amplitude'],
            trigger_dict['repetition'],
            recording.stddev,
            threshold,
            tick_dt_ms,
            epoch_timing_ms[0],
            epochs
        ))
    return output_spikes, output_epochs


########################## *BACKEND FUNCTIONS* #########################

def blank_artefacts(
        signal_data: np.typing.NDArray[np.floating],
        triggers: abc.Iterable[int],
        pulsewidth_ms: float,
        tick_dt_ms: float
) -> np.typing.NDArray[np.floating]:
    blank_width = int(
        (pulsewidth_ms+constants.core.BLANKING_PADDING_MS)/tick_dt_ms
    )
    for trigger in triggers:
        signal_data[trigger:trigger+blank_width] = np.zeros(blank_width)
    return signal_data


def detect_spikes(
        signal_data: np.typing.NDArray[np.floating],
        triggers: abc.Iterable[int],
        epoch_timing_ms: tuple[float, float],
        tick_dt_ms: float,
        threshold: float,
        phase: str,
        stim_type: str
        # epoch_adjustment: abc.Callable[[int], int]
) -> tuple[
    dict[int, list[classes.Spike]],
    list[classes.DataEpoch]
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
    # Define epoch bounds in samples relative to the trigger:
    epoch_start = int(epoch_timing_ms[0]/tick_dt_ms)
    epoch_end = int(epoch_timing_ms[1]/tick_dt_ms) + 1
    # Create empty output lists:
    spikes = {}
    data_epochs = []
    for i, trigger in enumerate(triggers):
        trigger_spikes = []
        # Create `DataEpoch` object describing current epoch and add to
        # `data_epochs` list:
        trace = np.array(signal_data[trigger+epoch_start:trigger+epoch_end])
        data_epochs.append(classes.DataEpoch(
            trace,
            i,
            phase,
            stim_type
        ))
        # Detect peaks greater than `threshold` within current epoch and
        # add to `spikes` list:
        peaks, properties = signal.find_peaks(
            trace,
            threshold,
            distance=constants.core.MINIMUM_SPIKE_DISTANCE_MS/tick_dt_ms
        )
        for j, peak in enumerate(peaks):
            trigger_spikes.append(classes.Spike(
                i,
                peak*tick_dt_ms + epoch_timing_ms[0],
                properties['peak_heights'][j]
            ))
        spikes[i] = trigger_spikes
    # Return populated output lists:
    return utils.remove_empty_lists(spikes), data_epochs


def filter_signal(
        signal_data: np.typing.NDArray[np.floating],
        bessel_coefficients: tuple[np.typing.ArrayLike, np.typing.ArrayLike],
        notch_coefficients: tuple[np.typing.ArrayLike, np.typing.ArrayLike]
) -> np.typing.NDArray[np.floating]:
    notch_filtered = signal.filtfilt(*notch_coefficients, signal_data)
    bessel_filtered = signal.filtfilt(*bessel_coefficients, notch_filtered)
    return bessel_filtered

def filter_spikes(
        spikes_dict: abc.Mapping[int, abc.Iterable[classes.Spike]],
        spike_criteria: classes.SpikeCriteria,
        phase_code: str,
) -> tuple[
    dict[int, list[classes.Spike]],
    dict[str, list[classes.Spike]],
    dict[str, list[classes.Spike]],
    int
]:
    spike_list = list(itertools.chain.from_iterable(spikes_dict.values()))
    filtered_spikes = {i: [] for i in spikes_dict.keys()}
    rejected_spikes = {
        f'{phase_code}-{i}': [] for i in spikes_dict.keys()
    }
    failed_spikes_list = []
    failed_counter = {i: False for i in spikes_dict.keys()}
    for spike in spike_list:
        if passes_spikecriteria(spike, spike_criteria):
            previous_spike = (
                filtered_spikes[spike.epoch_number][-1] if
                filtered_spikes[spike.epoch_number] else
                classes.Spike(
                    spike.epoch_number,
                    -2*constants.core.MINIMUM_ISI_MS,
                    0.0
                )
            )
            if ((spike.time_ms-previous_spike.time_ms)<
                constants.core.MINIMUM_ISI_MS):
                failed_spikes_list.append([previous_spike, spike])
                failed_counter[spike.epoch_number] = True
            # New spike must be appended AFTER checking last entry!
            filtered_spikes[spike.epoch_number].append(spike)
        else:
            rejected_spikes[f'{phase_code}-{spike.epoch_number}'].append(spike)
    failed_spikes = {
        f'{phase_code}-{x[0].epoch_number}.{i}': x for i, x in
        enumerate(failed_spikes_list)
    }
    failed_count = sum(failed_counter.values())
    return filtered_spikes, rejected_spikes, failed_spikes, failed_count


def max_amp_multiplier(
        markers: abc.Iterable[classes.Marker],
        stim_type: str
) -> float:
    """Determines greatest amplitude multiplier used in amplitude sweep.
    
    Returns largest multiple greater than 0, or 0 if there aren't any,
    for a single stimulation type over amplitude sweep (i.e. multiple
    trials). Since minimum returned value is 0 and absolute value is not
    used, **should only be used if all amplitude multipliers are
    positive**.

    # Arguments
    * `markers` -- list of `Marker` objects relating to amplitude trial
    * `stim_type` -- `'mech'` or `'elec'` indicating stimulation type to
    calculate maximum amplitude for.

    # Error Handling
    If the comment for a given trial does not match the known pattern or
    indicates that the trial is not part of an amplitude sweep, prints
    an error message and skips that trial.
    """
    max_amp = 0
    for marker in markers:
        m = sweep_pattern.search(marker.comment)
        try:
            assert m.group('testvar') == 'amplitude' # type: ignore
            max_amp = max(
                max_amp,
                float(m.group(f'{stim_type}val')) # type: ignore
            )
        except (AssertionError, AttributeError):
            # AssertionError is raised if test type is not amplitude
            # AttributeError is raised if comment does not match pattern
            print(
                "Trial which is not part of an amplitude sweep has been"
                "passed to amplitude sweep pipeline"
                f"({marker.comment})."
            )
    return max_amp


def output_filename(
        save_type: str,
        animal_id: str,
        position: int,
        repetition: int,
        recording_number: int,
        test_code: str
) -> str:
    """Builds standardised filename for single-trial outputs.
    
    Currently used for all `main_part1` outputs (clusters, frs, epochs,
    and spikes).
    """
    unit_id = f'{animal_id}-{position}'
    recording_file_id = f'[{repetition}-{recording_number}]'
    return '_'.join([
        save_type,
        unit_id,
        recording_file_id,
        test_code,
        str(utils.current_version)
    ])


def overlay_epochs(
        epochs_trial: classes.EpochsTrial,
        animal_id: str,
        position: int,
        test: str
) -> None:
    # Initialise variables to identify trial if error is raised:
    test_stim = epochs_trial.test_stim
    test_frequency = epochs_trial.test_frequency
    test_amplitude = epochs_trial.test_amplitude
    for epoch in epochs_trial.epochs:
        # If an epoch has an invalid stimulus type, try-except block
        # prints an error message and skips it:
        try:
            plot_id = (constants.experiments.STIMULATION_TYPES.
                       index(epoch.stimulus) + 1)
            pulsewidth_ms = constants.experiments.PULSEWIDTH_MS[epoch.stimulus]
        except ValueError:
            print(
                "Epoch stimulus type not recognised.\n"
                f"    Animal ID: {animal_id}\n"
                f"    Position: {position}\n"
                f"    Test: {test}\n"
                f"    Test Stimulus: {test_stim}\n"
                f"    Test Frequency: {test_frequency} Hz"
                "\n"
                f"    Test Amplitude: {test_amplitude} μA"
            )
            continue
        # Plot figure:
        phase_index = constants.experiments.PHASES.index(epoch.phase)
        colour = list(constants.plots.PALETTE.values())[phase_index+3] + '50'
        blanking_end = int(
            (pulsewidth_ms+constants.core.BLANKING_PADDING_MS)/
            epochs_trial.tick_dt_ms
        )
        x_ms = [i * epochs_trial.tick_dt_ms + epochs_trial.start_ms
                for i,_ in enumerate(epoch.trace)]
        plt.subplot(2, 1, plot_id)
        plt.plot(
            x_ms[:blanking_end],
            epoch.trace[:blanking_end],
            c=constants.plots.PALETTE['light_grey']+'50',
            linewidth=constants.plots.TRACE_LINEWIDTH,
            zorder=0
        )
        plt.plot(
            x_ms[blanking_end-1:],
            epoch.trace[blanking_end-1:],
            c=colour,
            linewidth=constants.plots.TRACE_LINEWIDTH,
            zorder=0
        )
    return


def passes_spikecriteria(
        spike: classes.Spike,
        spike_criteria: classes.SpikeCriteria
) -> bool | np.bool:
    latency_min = (
        spike.time_ms >= spike_criteria.latency_min_ms if
        spike_criteria.latency_min_ms is not None else True
    )
    latency_max = (
        spike.time_ms < spike_criteria.latency_max_ms if
        spike_criteria.latency_max_ms is not None else True
    )
    size_min = (
        spike.size_uV >= spike_criteria.size_min_uV if
        spike_criteria.size_min_uV is not None else True
    )
    size_max = (
        spike.size_uV < spike_criteria.size_max_uV if
        spike_criteria.size_max_uV is not None else True
    )
    return latency_min and latency_max and size_min and size_max


def plot_filter_bounds(
        plot_id: int,
        threshold: float,
        spike_criteria: classes.SpikeCriteria,
        xlim: tuple[float, float],
        ylim: tuple[float, float]
) -> None:
    plt.subplot(2, 1, plot_id)
    plt.plot(
        xlim,
        [threshold, threshold],
        c=constants.plots.PALETTE['black'],
        linewidth=constants.plots.CLUSTER_THRESHOLD_LINEWIDTH/2,
        zorder=1
    )
    for x, y in (
        (xlim, [spike_criteria.size_min_uV, spike_criteria.size_min_uV]),
        (xlim, [spike_criteria.size_max_uV, spike_criteria.size_max_uV]),
        ([spike_criteria.latency_min_ms, spike_criteria.latency_min_ms], ylim),
        ([spike_criteria.latency_max_ms, spike_criteria.latency_max_ms], ylim)
    ):
        plt.plot(
            x, # type: ignore
            y, # type: ignore
            c=constants.plots.PALETTE['green'],
            linewidth=constants.plots.CLUSTER_THRESHOLD_LINEWIDTH,
            zorder=1
        )
    return


def plot_peaks(
        plot_id: int,
        rejected_peaks: abc.Iterable[abc.Iterable[classes.Spike]],
        filtered_phases: abc.Iterable[
            abc.Mapping[Any, abc.Iterable[classes.Spike]]
        ],
        failed_peaks: abc.Iterable[abc.Iterable[classes.Spike]]
) -> float:
    # Set subplot:
    plt.subplot(2, 1, plot_id)
    # Plot rejected peaks:
    plt.scatter(
        [x.time_ms for spike_list in rejected_peaks for x in spike_list],
        [x.size_uV for spike_list in rejected_peaks for x in spike_list],
        constants.plots.CLUSTER_POINT_SIZE,
        constants.plots.PALETTE['dark_grey'],
        zorder=3
    )
    # Concatenate filtered peaks:
    filtered_peaks = list(itertools.chain.from_iterable(
        list(itertools.chain.from_iterable(
            phase_dict.values() for phase_dict in filtered_phases
        ))
    ))
    # Plot filtered peaks and calculate largest voltage:
    plt.scatter(
        [x.time_ms for x in filtered_peaks],
        [x.size_uV for x in filtered_peaks],
        constants.plots.CLUSTER_POINT_SIZE,
        constants.plots.PALETTE['black'],
        zorder=4
    )
    try:
        max_voltage = max(x.size_uV for x in filtered_peaks)
    except ValueError:
        max_voltage = 0
    # Plot failed peaks:
    for failed_pair in failed_peaks:
        x_points = [x.time_ms for x in failed_pair]
        y_points = [x.size_uV for x in failed_pair]
        plt.scatter(
            x_points,
            y_points,
            constants.plots.CLUSTER_POINT_SIZE,
            constants.plots.PALETTE['vermillion'],
            zorder=5
        )
        plt.plot(
            x_points,
            y_points,
            c=constants.plots.PALETTE['vermillion'],
            linewidth=constants.plots.TRACE_LINEWIDTH,
            zorder=5
        )
    return max_voltage


def process_signal(
        signal_data: np.typing.NDArray[np.floating],
        mech_triggers: abc.Iterable[int],
        elec_triggers: abc.Iterable[int],
        tick_dt_ms: float,
        bessel_coefficients: tuple[np.typing.ArrayLike, np.typing.ArrayLike],
        notch_coefficients: tuple[np.typing.ArrayLike, np.typing.ArrayLike]
) -> np.typing.NDArray[np.floating]:
    signal_data = blank_artefacts(
        signal_data,
        mech_triggers,
        constants.experiments.PULSEWIDTH_MS['mechanical'],
        tick_dt_ms
    )
    signal_data = blank_artefacts(
        signal_data,
        elec_triggers,
        constants.experiments.PULSEWIDTH_MS['electrical'],
        tick_dt_ms
    )
    return filter_signal(
        signal_data,
        bessel_coefficients,
        notch_coefficients
    )


def read_record(
        data: adi.read.File,
        record_number: int,
        print_info: bool = True
) -> dict[str, Any]:
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
    record = data.records[record_number]
    tick_dt = record.tick_dt
    # Extract signal and trigger data:
    # The signal is inverted during this step.
    signal_data = data.channels[0].get_data(record_number+1)*-1
    mechstim = data.channels[1].get_data(record_number+1)
    elecstim = data.channels[2].get_data(record_number+1)
    # Extract comments and adjust timings to match trimmed signals:
    #!This will fail if comments are not stored sequentially
    # Solution: make a list of comments from `record.comments` which has
    # explicitly sorted by `record.comments.time`
    markers = []
    for i, comment in enumerate(record.comments):
        start = int(comment.time/tick_dt)
        if i < len(record.comments)-1:
            end = int(record.comments[i+1].time/tick_dt)
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
            end = int(len(signal_data)/tick_dt)
            # assert start < end, (
            #     "Artefact trimming has resulted in this comment falling "
            #     "outside the usable portion of the recording.\n"
            #     f"    Comment: {comment.text} ({i})\n"
            #     f"    Start: {start}\n"
            #     f"    End: {end}\n"
            #     "Consult the file to decide if this recording segment should "
            #     "be used, and skip it if it is not useful. However, if it is "
            #     "useful, either reduce `artefact_width_samples` or modify "
            #     "this function to allow manual skipping of specific comments. "
            #     "Another solution would be to modify the function such that "
            #     "the two problems raise different error types, and handle "
            #     "this problem by automatically skipping any recording "
            #     "segments which raise it."
            # )
        # Add marker to list:
        markers.append(classes.Marker(comment.text, start, end))
    # Print the duration of the record for inspection purposes:
    if print_info:
        print("[%s]: %.1f s" % (record_number, len(signal_data)*tick_dt))
    # Calculate filter arrays:
    bessel = signal.bessel(
        constants.core.BESSELHIGHPASS_N,
        constants.core.BESSELHIGHPASS_FREQUENCY,
        'highpass',
        fs=1/tick_dt
    )
    notch = signal.iirnotch(
        constants.core.NOTCHFILT_FREQUENCY,
        constants.core.NOTCHFILT_Q,
        fs=1/tick_dt
    )
    # Return a dictionary of partial data for a `Recording` object:
    return {
        'tick_dt': tick_dt,
        'signal_data': signal_data,
        'mech_triggers': mechstim,
        'elec_triggers': elecstim,
        'markers': markers,
        'bessel': bessel,
        'notch': notch
    }


def separate_sweep_phases(
        test: str,
        triggers_info: classes.TriggersInfo
) -> classes.TriggersTrial:
    """Generates `TriggersTrial` object from paired trigger data.
    
    # Error Handling
    * Raises `AssertionError` if neither stimulation value is 0 or both
    are 0, as that should not occur.
    * May raise `KeyError` if test type is not `'frequency'`,
    `'amplitude'`, `'long-duration'`, or `'nine-one'` (test type is not
    checked if the test is not `'nine-one'` and the conditioning
    stimulus is electrical).
    """
    assert bool(triggers_info.mech_value) != bool(triggers_info.elec_value), (
        "Neither stimulation value is 0 or both stimulation values are "
        "0 , suggesting that stimulation was performed incorrectly."
    )
    frequency = constants.experiments.DEFAULT_CONDITIONING_FREQUENCY
    minor_frequency = (constants.experiments.
                       DEFAULT_CONDITIONING_FREQUENCY_MINOR)
    amplitude = constants.experiments.DEFAULT_CONDITIONING_AMPLITUDE
    duration = constants.experiments.DEFAULT_CONDITIONING_DURATION_SECONDS
    # Determine the separation points between experimental phases:
    if test == 'nine-one':
        minor_frequency = (constants.experiments.
                           NINEONE_CONDITIONING_FREQUENCY_MINOR)
        minor_stimulation_count = math.floor(minor_frequency * duration)
        if triggers_info.elec_value:
            test_stim = 'electrical'
            start_itlv_index = minor_stimulation_count
            stim_value = triggers_info.elec_value
        else:
            test_stim = 'mechanical'
            start_itlv_index = (math.floor(frequency * duration) -
                                minor_stimulation_count)
            stim_value = triggers_info.mech_value
        frequency -= minor_frequency
    else:
        # if triggers_info.elec_value:
        #     frequency = 0
        #     test_stim = 'electrical'
        #     stim_value = triggers_info.elec_value
        # else:
        #     test_stim = 'mechanical'
        #     stim_value = triggers_info.mech_value
        #     if test == 'frequency':
        #         frequency = stim_value
        #     elif test == 'amplitude':
        #         pass
        #     elif test == 'long-duration':
        #         frequency = constants.experiments.LONG_CONDITIONING_FREQUENCY
        #         duration = (constants.experiments.
        #                     LONG_CONDITIONING_DURATION_SECONDS)
        #     else:
        #         raise KeyError(f"Test type not recognised ({test}).")
        # start_itlv_index = math.floor(frequency*duration)
        if triggers_info.elec_value:
            test_stim = 'electrical'
            stim_value = triggers_info.elec_value
        else:
            test_stim = 'mechanical'
            stim_value = triggers_info.mech_value
        if test == 'frequency':
            frequency = stim_value
        elif test == 'amplitude':
            amplitude = stim_value
        elif test == 'long-duration':
            frequency = constants.experiments.LONG_CONDITIONING_FREQUENCY
            duration = (constants.experiments.
                        LONG_CONDITIONING_DURATION_SECONDS)
        else:
            raise KeyError(f"Test type not recognised ({test}).")
        start_itlv_index = (0 if triggers_info.elec_value else
                            math.floor(frequency * duration))
    start_rcvr_index = (
        start_itlv_index +
        constants.experiments.INTERLEAVED_FREQUENCY *
        constants.experiments.INTERLEAVED_DURATION_SECONDS
    )
    start_itlv = triggers_info.mech_triggers[start_itlv_index]
    start_rcvr = triggers_info.mech_triggers[start_rcvr_index]
    # Collect triggers into lists according to phase:
    triggers_mech_cond = [x for x in triggers_info.mech_triggers if
                          x<start_itlv]
    triggers_elec_cond = [x for x in triggers_info.elec_triggers if
                          x<start_itlv]
    triggers_mech_itlv = [x for x in triggers_info.mech_triggers if
                          start_itlv<=x<start_rcvr]
    triggers_elec_itlv = [x for x in triggers_info.elec_triggers if
                          start_itlv<=x<start_rcvr]
    triggers_mech_rcvr = [x for x in triggers_info.mech_triggers if
                          x>=start_rcvr]
    triggers_elec_rcvr = [x for x in triggers_info.elec_triggers if
                          x>=start_rcvr]
    # Create and return a `TriggerSet` object:
    return classes.TriggersTrial(
        test,
        test_stim,
        frequency,
        minor_frequency,
        amplitude,
        triggers_mech_cond,
        triggers_elec_cond,
        triggers_mech_itlv,
        triggers_elec_itlv,
        triggers_mech_rcvr,
        triggers_elec_rcvr
    )


# def test_variables(
#         test: str,
#         triggers_by_phase: classes.TriggersTrial
# ) -> tuple[float, float, float]:
#     test_frequency = constants.experiments.DEFAULT_CONDITIONING_FREQUENCY
#     test_frequency_minor = (constants.experiments.
#                             DEFAULT_CONDITIONING_FREQUENCY_MINOR)
#     test_amplitude = constants.experiments.DEFAULT_CONDITIONING_AMPLITUDE
#     if test == 'frequency':
#         test_frequency = triggers_by_phase.stim_value
#     elif test == 'amplitude':
#         test_amplitude = triggers_by_phase.stim_value
#     elif test == 'long-duration':
#         test_frequency = constants.experiments.LONG_CONDITIONING_FREQUENCY
#     elif test == 'nine-one':
#         test_frequency_minor = (constants.experiments.
#                                 NINEONE_CONDITIONING_FREQUENCY_MINOR)
#         test_frequency -= test_frequency_minor
#     return test_frequency, test_frequency_minor, test_amplitude


def threshold_stddev(
        signal_data: np.typing.NDArray[np.floating],
        bessel: tuple[np.typing.ArrayLike, np.typing.ArrayLike],
        notch: tuple[np.typing.ArrayLike, np.typing.ArrayLike],
        tick_dt: float,
        markers: abc.Sequence[classes.Marker],
        show_threshold: bool = True
) -> float:
    artefact_width_samples = int(constants.core.NOTCHFILT_ARTEFACT_WIDTH_S/
                                 tick_dt)
    signal_filtered = filter_signal(signal_data, bessel, notch)
    assert artefact_width_samples < markers[0].start_sample
    premarker = signal_filtered[artefact_width_samples:markers[0].start_sample]
    stddev_premarker = statistics.stdev(premarker)
    if show_threshold:
        threshold = stddev_premarker * constants.core.SPIKE_NOISE_FLOOR
        plot_x = [i*tick_dt for i, _ in enumerate(premarker)]
        print(
            "Threshold: %.2f mV (%d × %.2f mV)" %
            (threshold, constants.core.SPIKE_NOISE_FLOOR, stddev_premarker)
        )
        plt.plot(plot_x, premarker)
        plt.plot([min(plot_x), max(plot_x)], [stddev_premarker, stddev_premarker])
        plt.plot([min(plot_x), max(plot_x)], [threshold, threshold])
    return stddev_premarker


def triggers(
        trigger_data: np.typing.NDArray[np.floating],
        trigger_threshold: float
) -> list[int]:
    """Returns onset index of all triggers in `trigger_data` as list."""
    # Convert `trigger_data` into list of binary values according to
    # whether each item exceeds `trigger_threshold`:
    triggers_binary = [x>trigger_threshold for x in trigger_data]
    # Get indices of left edges:
    trigger_indices = utils.detect_edges(triggers_binary, 'rising')
    # Correct index of left edge (i.e. time in samples) for known delay
    # between trigger and recording channels:
    return [x-constants.core.TRIGGER_DELAY_SAMPLES for x in trigger_indices]


def trigger_info(
        marker: classes.Marker,
        recording: classes.Recording,
        trigger_threshold_mech: float,
        trigger_threshold_elec: float
) -> classes.TriggersInfo:
    if recording.test == 'frequency' or recording.test == 'amplitude':
        p = sweep_pattern
    elif recording.test == "long-duration":
        p = long_pattern
    elif recording.test == 'nine-one':
        p = nineone_pattern
    else:
        raise KeyError(f"Test type not recognised ({recording.test}).")
    # Parse marker text and check that it fits expected format:
    m = p.search(marker.comment)
    adjustment = 1
    try:
        if recording.test == 'long-duration':
            mech_value = float(m.group('stimtype')=='mech') # type: ignore
            elec_value = float(m.group('stimtype')=='elec') # type: ignore
        else:
            try:
                assert (recording.test=='nine-one' or
                        recording.test==m.group('testvar')) # type: ignore
            except AssertionError:
                print(
                    "Comment indicates a different test type to filename "
                    f"({marker.comment})."
                )
                raise
            # If test was frequency sweep, values are test frequencies
            # If test was amplitude sweep, values are test amplitudes
            # If test was 9:1, values indicate primary stimulus
            mech_value = float(m.group('mechval')) # type: ignore
            elec_value = float(m.group('elecval')) # type: ignore
            # Extract values from marker text:
            if recording.test == 'amplitude':
                max_amplitude = max(mech_value, elec_value)
                if (
                    0.5 <= max_amplitude <= 2 and
                    min(mech_value, elec_value) == 0
                ):
                    adjustment = min(1, max_amplitude)
                else:
                    raise IndexError(
                        "Amplitude values are invalid and will interfere "
                        "with trigger detection. One amplitude should be "
                        "0 and the other between 0.5 and 2 (inclusive).\n"
                        f"Mechanical amplitude: {mech_value}\n"
                        f"Electrical amplitude: {elec_value}"
                    )
    # Print error message and skip marker if no match is found:
    except AttributeError:
        # AttributeError is raised if m is None (i.e. if regex
        # pattern didn't match)
        print(
            "Comment does not match the expected format for "
            f"{recording.test} trials ({marker.comment})."
        )
        raise
    # Detect triggers:
    return classes.TriggersInfo(
        mech_value,
        triggers(
            recording.mech_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_mech * adjustment
        ),
        elec_value,
        triggers(
            recording.elec_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_elec * adjustment
        )
    )


def trigger_thresholds(
        recording: classes.Recording
        # trigger_data: np.typing.NDArray[np.floating],
        # correction_factor: float = 1.0
) -> tuple[float, float]:
    """Returns value of trigger peaks in `trigger_data`."""
    # Define trigger thresholds:
    if recording.test == 'amplitude':
        max_amp_mech = max_amp_multiplier(recording.markers, "mech")
        max_amp_elec = max_amp_multiplier(recording.markers, "elec")
        trigger_correction_mech = max_amp_mech if max_amp_mech!=0 else 1
        trigger_correction_elec = max_amp_elec if max_amp_elec!=0 else 1
    else:
        trigger_correction_mech = 1
        trigger_correction_elec = 1
    # Previously, this looked for the minimum value within the trigger
    # data (i.e. the peak discharge amplitude) and corrected it to peak
    # amplitude. I believe that this was done to solve an imaginary
    # issue so have changed it to simply look for the peak amplitude
    # directly. However, if issues arise, this is a potential cause.
    trigger_value_mech = round(
        max(recording.mech_triggers)/trigger_correction_mech,
        1
    )
    trigger_value_elec = round(
        max(recording.elec_triggers)/trigger_correction_elec,
        1
    )
    return (
        trigger_value_mech * constants.core.TRIGGER_DETECTION_NOISE_WINDOW,
        trigger_value_elec * constants.core.TRIGGER_DETECTION_NOISE_WINDOW
    )


def trim_1darray(trace: np.ndarray, trim_width: int) -> np.ndarray:
    """Removes `trim_width` values from each end of 1D NumPy array."""
    trimmed_start = trim_width
    trimmed_end = len(trace)-trim_width
    return trace[trimmed_start:trimmed_end]
