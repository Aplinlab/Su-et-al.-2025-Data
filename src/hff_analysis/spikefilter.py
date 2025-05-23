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

import matplotlib.pyplot as plt
import numpy as np
import pathlib
import typing

from . import classes
from . import constants
from . import utils


####################### *USER-EXPOSED FUNCTIONS* #######################


def filter_spikes(
        trials: list[classes.spikes.SpikesTrial],
        criteria: dict[str, dict[str, int | float | None]],
        exclude_frequencies: list[int | float],
        exclude_amplitudes: list[int | float]
) -> list[classes.spikes.SpikesTrial]:
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
    # Convert `criteria` from `dict` to `SpikeCriteria`:
    try:
        spike_criteria_mech = classes.recordings.SpikeCriteria(
            **criteria['mechanical']
        )
        spike_criteria_elec = classes.recordings.SpikeCriteria(
            **criteria['electrical']
        )
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
        for (argument, phase) in {
            'spikes_cond': trial.conditioning,
            'spikes_itlv': trial.interleaved,
            'spikes_rcvr': trial.recovery
        }.items():
            # Apply spike exclusion criteria:
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
            filtered_phases[argument] = classes.spikes.SpikesPhase(
                phase.epochs_mech,
                phase.epochs_elec,
                filtered_spikes_mech,
                filtered_spikes_elec
            )
        # Generate `SpikesTrial` object from filtered spikes:
        filtered_trials.append(classes.spikes.SpikesTrial(
            trial.animal_id,
            trial.position,
            trial.test,
            trial.test_stim,
            trial.test_frequency,
            trial.test_amplitude,
            trial.repetition,
            **filtered_phases
        ))
    # Return complete list of filtered peaks:
    return filtered_trials
    

def plot_clusters(
        spikes: list[classes.spikes.SpikesTrial],
        epochs: list[classes.epochs.EpochsTrial],
        repetition: int,
        recording_segment: int,
        exclude_frequencies: list[int | float] = [],
        exclude_amplitudes: list[int | float] = [],
        save_figures: bool = False
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
    # Top subplot displays mechanical stimulation epochs and lower epoch
    # displays electrical stimulation epochs. `zorder` is set so points
    # appear above traces.

    # Initialise variables for later:
    plt.figure(figsize=constants.FIGSIZE)
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
    colour_list = list(constants.PALETTE.values())
    # Plot detected peaks as scatterplots:
    for spikes_trial in spikes:
        # Map points from each stim type to correct subplot:
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
                    constants.PALETTE['black'],
                    zorder=1
                )
                try:
                    max_voltage = max(
                        max_voltage,
                        max(x.size_uV for x in phase_stim)
                    )
                except ValueError:
                    # ValueError is raised if no spikes in a given phase
                    pass
    # Plot traces by epoch:
    for epochs_trial in epochs:
        # Test for trial exclusion criteria:
        # `filter_spikes()` only applies these to spikes, not epochs
        if (epochs_trial.test_frequency in exclude_frequencies or
            epochs_trial.test_amplitude in exclude_amplitudes):
            continue
        # Initialise variables to identify trial if error is raised:
        test_stim = epochs_trial.test_stim
        test_frequency = epochs_trial.test_frequency
        test_amplitude = epochs_trial.test_amplitude
        for epoch in epochs_trial.epochs:
            # If an epoch has an invalid stimulus type, try-except block
            # prints an error message and skips it:
            try:
                plot_id = constants.STIMULATION_TYPES.index(epoch.stimulus)
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
            phase_index = constants.EXPERIMENTAL_PHASES.index(epoch.phase)
            colour = colour_list[phase_index+1] + '50'
            plt.subplot(2, 1, plot_id+1)
            plt.plot(
                [i * epoch.tick_dt_ms + epoch.start_ms
                for i, _x in enumerate(epoch.trace)],
                epoch.trace,
                c=colour,
                linewidth=constants.CLUSTER_LINE_WIDTH,
                zorder=0
                )
    # Prettify figure:
    plot_title = (
        f'{animal_id.upper()}-{position} [{repetition}-{recording_segment}] ('
        f'{constants.METADATA[animal_id.upper()][position]})'
    )
    for plot_id in range(1,3):
        plt.subplot(2, 1, plot_id)
        plt.xlabel("Time (ms)")
        plt.ylabel("Signal voltage (μV)")
        ax = plt.gca()
        ax.set_ylim((
            constants.CLUSTER_YMIN,
            max_voltage*constants.CLUSTER_YMAX_SCALE
        ))
    plt.subplot(2, 1, 1)
    plt.title("Mechanical")
    plt.subplot(2, 1, 2)
    plt.title("Electrical")
    plt.suptitle(plot_title)
    plt.tight_layout()
    # Save figures if specified:
    if save_figures:
        plot_type = 'clusters'
        target_name = output_filename(
            plot_type,
            animal_id,
            position,
            repetition,
            recording_segment,
            test[0:4]
        )
        utils.save_plot(
            plot_type,
            target_name
        )
    return


def save_to_json(
        input_filename: str,
        repetition: int,
        recording_segment: int,
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
        file_path = (constants.SAVE_PATHS['json_root'] +
                    constants.SAVE_PATHS[json_type])
    except KeyError:
        file_path = (constants.SAVE_PATHS['json_root'] +
                    f'{json_type}\\')
    input_info = classes.recordings.FilenameInfo.from_filename(input_filename)
    filename = output_filename(
        json_type,
        input_info.animal_id,
        input_info.position,
        repetition,
        recording_segment,
        input_info.test[0:4]
    ) + '.json'
    # Construct dict to save as JSON:
    output_dict = {
        'json_type': json_type,
        'version': str(utils.current_version),
        'repetition': repetition,
        'recording_segment': recording_segment
    }
    for key, value in kwargs.items():
        output_dict[key] = convert_to_json_dict(value)
    # Save to JSON and print success/failure message:
    try:
        utils.confirm_save(file_path, filename, output_dict, force_overwrite)
    except FileNotFoundError:
        # Make directories if they don't exist
        pathlib.Path(file_path).mkdir(parents=True, exist_ok=True)
        utils.confirm_save(file_path, filename, output_dict, force_overwrite)
    except TypeError:
        # Raise error if `output_dict` contains non-JSON serialisable
        # objects:
        raise
    return


########################## *BACKEND FUNCTIONS* #########################

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
        str(utils.current_version)
    ])
