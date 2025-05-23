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
* `long_pattern` - regex pattern for parsing comments in frequency and
amplitude sweeps.
* `sweep_pattern` - regex pattern for parsing comments in frequency and
amplitude sweeps.
"""

import adi
import json
import math
import numpy as np
from numpy.typing import NDArray
import re
from scipy import signal
import typing

from . import classes
from . import constants


# Regex patterns for parsing marker comments:
# Compiled once, rather than within a function, to reduce computations.
long_pattern = re.compile(constants.LONGDURATION_REGEX)
sweep_pattern = re.compile(constants.SWEEP_REGEX)


####################### *USER-EXPOSED FUNCTIONS* #######################

def load_filereadsettings(
        frs_filename: str
) -> classes.recordings.FileReadSettings:
    """Loads `FileReadSettings` object from JSON file.
    
    The returned object can be used to reproduce an earlier run of
    `main_part1`.
    """
    # Construct file path:
    frs_filepath = (constants.SAVE_PATHS['json_root'] +
                    constants.SAVE_PATHS['frs'])
    save_file = frs_filepath + frs_filename
    # Add file extension if not already present:
    if save_file[-5:].lower() != '.json':
        save_file += '.json'
    # Load file as `FileReadSettings` object:
    frs = classes.recordings.FileReadSettings.from_dict(json.load(open(
        save_file
    )))
    return frs


def read_adicht(
        filename: str,
        data_segments: list | None = None
) -> list[classes.recordings.Recording]:
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
    filename_info = classes.recordings.FilenameInfo.from_filename(filename)
    # Read the specified file and break it into recording segments (each
    # time recording was started/stopped within the file), extracting
    # data from each segment separately:
    # Note: Previously, all segments were stitched together and treated
    # as a single recording. However, I believe that separating them
    # allows more flexibility in handling data (for example, if a file
    # contains multiple thresholding sweeps and only the final one is of
    # interest). It should also improve the notch filter results.
    data = adi.read_file(
        constants.RAW_DATA_PATH + filename_info.animal_id + '\\' +
        filename_info.name + filename_info.extension
    )
    # If a list has been specified using the `data_segments` argument,
    # it will be used; otherwise, all available segments will be read.
    if data_segments:
        records = [read_record(data, i) for i in data_segments]
    else:
        records = [read_record(data, i) for i, _x in
                   enumerate(data.records)]
    # Return list of `Recording` objects, filling in attributes which
    # `read_record()` cannot:
    return [
        classes.recordings.Recording(
            filename_info.animal_id,
            filename_info.position,
            filename_info.test,
            record['tick_dt'],
            record['signal_data'],
            record['mech_triggers'],
            record['elec_triggers'],
            record['markers']
        ) for record in records
    ]

def spikes_info(
        recording: classes.recordings.Recording,
        repetition: int,
        epoch_timing_ms: tuple[int | float, int | float],
        threshold_uV: int | float
) -> tuple[list[classes.spikes.SpikesTrial], list[classes.epochs.EpochsTrial]]:
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
    if recording.test == 'frequency' or recording.test == 'amplitude':
        p = sweep_pattern
    elif recording.test == "long-duration":
        p = long_pattern
    else:
        raise KeyError(f"Test type not recognised ({recording.test}).")
    # Define trigger thresholds:
    if recording.test == 'amplitude':
        trigger_correction_mech = max_amp(recording.markers, "mech")
        trigger_correction_elec = max_amp(recording.markers, "elec")
    else:
        trigger_correction_mech = 1
        trigger_correction_elec = 1
    trigger_threshold_mech = (
        trigger_value(recording.mech_triggers, trigger_correction_mech)*
        constants.TRIGGER_DETECTION_NOISE_WINDOW
    )
    trigger_threshold_elec = (
        trigger_value(recording.elec_triggers, trigger_correction_elec)*
        constants.TRIGGER_DETECTION_NOISE_WINDOW
    )
    # Create empty lists to populate and return later:
    output_spikes = []
    output_epochs = []
    for marker in recording.markers:
        # Parse marker text and check that it fits expected format:
        m = p.search(marker.comment)
        adjustment = 1
        if recording.test == 'frequency' or recording.test == 'amplitude':
            try:
                assert m.group('testvar') == recording.test # type: ignore
                mech_value = float(m.group('mechval')) # type: ignore
                elec_value = float(m.group('elecval')) # type: ignore
            # Print error message and skip marker if no match is found:
            except AttributeError:
                # AttributeError is raised if m is None (i.e. if regex
                # pattern didn't match)
                print(
                    "Comment does not match the expected format for "
                    f"{recording.test} trials ({marker.comment})."
                )
                continue
            except AssertionError:
                print(
                    "Comment indicates a different test type to filename "
                    f"({marker.comment})."
                )
                continue
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
                        "with trigger detection. One amplitude should be 0 "
                        "and the other between 0.5 and 2 (inclusive).\n"
                        f"Mechanical amplitude: {mech_value}\n"
                        f"Electrical amplitude: {elec_value}"
                    )
            else:
                # Cannot be reached due to assertion at start of function,
                # but is included to fix `reportPossiblyUnboundError` later
                continue
        elif recording.test == 'long-duration':
            try:
                mech_value = int(m.group('stimtype') == 'mech') # type: ignore
                elec_value = int(m.group('stimtype') == 'elec') # type: ignore
            # Print error message and skip marker if no match is found:
            except AttributeError:
                # AttributeError is raised if m is None (i.e. if regex
                # pattern didn't match)
                print(
                    "Comment does not match the expected format for "
                    f"{recording.test} trials ({marker.comment})."
                )
                continue
        else:
            # Cannot be reached due to assertion at start of function,
            # but is included to fix `reportPossiblyUnboundError` later
            continue
        # Detect triggers:
        triggers_mech = triggers(
            recording.mech_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_mech * adjustment
        )
        triggers_elec = triggers(
            recording.elec_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_elec * adjustment
        )
        # Sort triggers by experimental phase using frequencies:
        triggers_by_phase = separate_sweep_phases(
            recording.test,
            mech_value,
            elec_value,
            triggers_mech,
            triggers_elec
        )
        test_frequency = constants.DEFAULT_CONDITIONING_FREQUENCY
        test_amplitude = constants.DEFAULT_CONDITIONING_AMPLITUDE
        if recording.test == 'frequency':
            test_frequency = triggers_by_phase.stim_value
        elif recording.test == 'amplitude':
            test_amplitude = triggers_by_phase.stim_value
        # Constrain signal_data to the current trial:
        signal_data_trial = recording.signal_data[
            marker.start_sample:marker.end_sample
        ]
        # Detect spikes within each epoch and extract traces:
        spike_detection_results = {
            phase: {
                stim_type: {} for stim_type in constants.STIMULATION_TYPES
            } for phase in constants.EXPERIMENTAL_PHASES
        }
        for phase in spike_detection_results.keys():
            for stim_type in spike_detection_results[phase]:
                (
                    spike_detection_results[phase][stim_type]['spikes'],
                    spike_detection_results[phase][stim_type]['epochs']
                ) = detect_spikes(
                    signal_data_trial,
                    triggers_by_phase.triggers[phase][stim_type],
                    epoch_timing_ms,
                    recording.tick_dt,
                    threshold_uV,
                    phase,
                    stim_type
                )
        # Collect variables used when defining both `SpikesTrial`
        # and `EpochsTrial` objects to reduce repetition:
        common_attributes = {
            'animal_id': recording.animal_id,
            'position': recording.position,
            'test': recording.test,
            'test_stim': triggers_by_phase.test_stim,
            'test_frequency': test_frequency,
            'test_amplitude': test_amplitude,
            'repetition': repetition
        }
        # Define `SpikesTrial` object:
        # It is important to count the number of triggers rather
        # supplying a precalculated value, since trimming of the
        # recording in an earlier step (to remove artefacts caused
        # by filtering) results in some markers being lost from the
        # recovery phase of the final trial in each recording.
        spikes_phases = {
            phase: classes.spikes.SpikesPhase(
                len(triggers_by_phase.triggers[phase]['mechanical']),
                len(triggers_by_phase.triggers[phase]['electrical']),
                spike_detection_results[phase]['mechanical']['spikes'],
                spike_detection_results[phase]['electrical']['spikes']
            ) for phase in constants.EXPERIMENTAL_PHASES
        }
        output_spikes.append(classes.spikes.SpikesTrial(
            **common_attributes,
            spikes_cond=spikes_phases['conditioning'],
            spikes_itlv=spikes_phases['interleaved'],
            spikes_rcvr=spikes_phases['recovery']
        ))
        # Define `EpochsTrial` object:
        data_epochs = (
            spike_detection_results['conditioning']['mechanical']['epochs']+
            spike_detection_results['conditioning']['electrical']['epochs']+
            spike_detection_results['interleaved']['mechanical']['epochs']+
            spike_detection_results['interleaved']['electrical']['epochs']+
            spike_detection_results['recovery']['mechanical']['epochs']+
            spike_detection_results['recovery']['electrical']['epochs']
        )
        # Check that number of epochs matches number of triggers:
        assert len(data_epochs) == sum(
            len(triggers_by_phase.triggers[phase]['mechanical'])+
            len(triggers_by_phase.triggers[phase]['electrical'])
            for phase in constants.EXPERIMENTAL_PHASES
        ), (
            "The number of epochs does not match the number of triggers "
            f"({marker.comment}). There may be an issue with the code."
        )
        output_epochs.append(classes.epochs.EpochsTrial(
            **common_attributes,
            epochs=data_epochs
        ))
    return (output_spikes, output_epochs)


########################## *BACKEND FUNCTIONS* #########################


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


def max_amp(markers: list[classes.recordings.Marker], stim_type: str):
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
                round(float(m.group(f'{stim_type}val'))) # type: ignore
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
        start_itlv_index = 0
        if elec_val == 0:
            test_stim = 'control'
            stim_value = 0
        else:
            test_stim = 'electrical'
            stim_value = elec_val
    elif elec_val == 0:
        test_stim = 'mechanical'
        stim_value = mech_val
        cond_frequency = constants.DEFAULT_CONDITIONING_FREQUENCY
        cond_duration = constants.SHORT_CONDITIONING_DURATION_SECONDS
        if test == 'frequency':
            cond_frequency = stim_value
        elif test == 'amplitude':
            pass
        elif test == 'long-duration':
            cond_duration = constants.LONG_CONDITIONING_DURATION_SECONDS
        else:
            raise KeyError(f"Test type not recognised ({test}).")
        start_itlv_index = math.floor(cond_frequency * cond_duration)
    else:
        # Raise error if neither stimulation value is 0:
        raise ValueError(
            "Neither stimulation value is 0, suggesting that"
            "stimulation was performed incorrectly."
        )
    start_rcvr_index = (
        start_itlv_index +
        constants.INTERLEAVED_FREQUENCY *
        constants.INTERLEAVED_DURATION_SECONDS
    )
    start_itlv = triggers_mech[start_itlv_index]
    start_rcvr = triggers_mech[start_rcvr_index]
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


def trim_1darray(trace: np.ndarray, trim_width: int) -> np.ndarray:
    """Removes `trim_width` values from each end of 1D NumPy array."""
    trimmed_start = trim_width
    trimmed_end = len(trace)-trim_width
    return trace[trimmed_start:trimmed_end]
