"""Functions exposed to the user.

See package docstring for function summaries.
"""

import adi
import json
import pathlib
import matplotlib.pyplot as plt
import pandas as pd

from . import classes
from . import constants
from . import utils
from . import updater
from . import freq


############################# *main_part1* #############################

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
        records = [utils.read_record(data, i) for i in data_segments]
    else:
        records = [utils.read_record(data, i) for i, _x in
                   enumerate(data.records)]
    # Populate `Recording` attributes which `read_record()` cannot:
    for record in records:
        record.animal_id = filename_info.animal_id
        record.position = filename_info.position
        record.test = filename_info.test
    return records


def update_outputs(
        save_outputs: str | list[str],
        input_filenames: list[str],
        input_folder: str | None = None,
        show_plots = False
) -> None:
    """Reads outdated FRS files and saves updated outputs.
    
    # Arguments
    * `save_outputs` -- string or list of strings indicating clusters to
    be saved:
      * `'clusters'` -- cluster plots as PDF image.
      * `'frs'` -- file_read_settings as JSON.
      * `'epochs'` -- traces for each epoch as JSON.
      * `'spikes'` -- spike data as JSON.
      * `'all'` -- all of the above.
    * `input_filenames` -- list of filenames which can identify FRS
    files to read. They do not have to the FRS files themselves, but
    can be other output files corresponding to the desired FRS files.
    **However, they must be of the same version as the FRS files to be
    read.**
    * `input_folder` -- path to location where relevant outputs are
    saved. For example, to read FRS files in the folder
    `.\\outputs\\archive\\v1.1.0\\JSON\\file_read_settings\\`, set
    `input_folder` to `'outputs\\archive\\v1.1.0'`.
    * `show_plots` -- whether cluster plots should be displayed.

    # Error Handling:
    * If an updater does not exist matching the versino of any filename
    in `input_filenames`, prints an error and skips that filename.
    """
    # Count number of updated files:
    updated_count = 0
    for filename in input_filenames:
        # Determine version of input file and pass to correct reader:
        file_version = classes.VersionNumber(filename)
        if file_version < classes.VersionNumber('2.0.0'):
            frs = updater.read_frs_v1(filename, input_folder)
        elif file_version >= utils.current_version:
            # Skip loop if file is up-to-date (does not update count):
            print(f"Up to date: {filename}")
            continue
        else:
            # Skip loop if matching updater not implemented:
            print(f"No updater for {file_version}: {filename}")
            continue
        # Set variables to be passed to other functions:
        if save_outputs == 'all' or save_outputs == ['all']:
            save_outputs = constants.OUTPUT_TYPES
        try:
            save_clusters = 'clusters' in save_outputs
        except TypeError:
            save_clusters = False
        force_overwrite = False
        # Read and process LabChart data:
        recording = read_adicht(frs.filename, [frs.recording_segment])[0]
        (spikes, epochs) = spikes_info(
            recording,
            frs.repetition,
            frs.epoch_timing_ms,
            frs.threshold_uV
        )
        filtered_spikes = filter_spikes(
            spikes,
            frs.spike_criteria,
            frs.exclude_frequencies,
            frs.exclude_amplitudes
        )
        # Draw cluster plots:
        if show_plots or save_clusters:
            plot_clusters(
                filtered_spikes,
                epochs,
                frs.repetition,
                frs.recording_segment,
                frs.exclude_frequencies,
                frs.exclude_amplitudes,
                save_clusters
            )
            if not show_plots:
                plt.close()
        # Save specified JSON outputs:
        if save_outputs:
            common_inputs = [
                frs.filename,
                frs.repetition,
                frs.recording_segment
            ]
            if 'frs' in save_outputs:
                save_to_json(
                    *common_inputs,
                    'frs',
                    force_overwrite,
                    filename=frs.filename,
                    epoch_timing_ms=frs.epoch_timing_ms,
                    threshold_uV=frs.threshold_uV,
                    spike_criteria=frs.spike_criteria,
                    exclude_frequencies=frs.exclude_frequencies,
                    exclude_amplitudes=frs.exclude_amplitudes
                )
            if 'epochs' in save_outputs:
                save_to_json(
                    *common_inputs,
                    'epochs',
                    force_overwrite,
                    epochs=epochs,
                    exclude_frequencies=frs.exclude_frequencies,
                    exclude_amplitudes=frs.exclude_amplitudes
                )
            if 'spikes' in save_outputs:
                save_to_json(
                    *common_inputs,
                    'spikes',
                    force_overwrite,
                    spikes=filtered_spikes
                )
        # Print completion message and update total updated files count:
        print(f"UPDATED from {file_version}: {frs.filename}")
        updated_count += 1
    # Print total updated files count:
    print(f"\nCOMPLETE: {updated_count} files updated.")
    return


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
    if recording.test == 'frequency':
        return freq.spikes_info(
            recording,
            repetition,
            epoch_timing_ms,
            threshold_uV
        )
    # elif recording.test == "amplitude":
    #     return ampl.spikes_info(
    #         recording,
    #         repetition,
    #         epoch_timing_ms,
    #         threshold_uV
    #     )
    # elif recording.test == "nine-one":
    #     return nine.analyse_nineone()
    # elif recording.test == "long-duration":
    #     return long.analyse_longduration()


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
        # `filtered_phases` will contain one list of `SpikesPhase`
        # objects for each experimental phase
        filtered_phases = []
        for phase in [trial.conditioning, trial.interleaved, trial.recovery]:
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
            filtered_phases.append(classes.spikes.SpikesPhase(
                phase.epochs_mech,
                phase.epochs_elec,
                filtered_spikes_mech,
                filtered_spikes_elec
            ))
        # Generate `SpikesTrial` object from filtered spikes:
        # Each list in `filtered_phases` is passed as different argument
        filtered_trials.append(classes.spikes.SpikesTrial(
            trial.animal_id,
            trial.position,
            trial.test,
            trial.test_stim,
            trial.test_frequency,
            trial.test_amplitude,
            trial.repetition,
            *filtered_phases
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
        ax.set_ylim([
            constants.CLUSTER_YMIN,
            max_voltage*constants.CLUSTER_YMAX_SCALE
        ])
    plt.subplot(2, 1, 1)
    plt.title("Mechanical")
    plt.subplot(2, 1, 2)
    plt.title("Electrical")
    plt.suptitle(plot_title)
    plt.tight_layout()
    # Save figures if specified:
    if save_figures:
        plot_type = 'clusters'
        target_name = utils.output_filename(
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
    filename = utils.output_filename(
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
        output_dict[key] = utils.convert_to_json_dict(value)
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


############################# *main_part2* #############################

def spikes_table(load_df_json: bool | str = False) -> pd.DataFrame:
    """Loads or builds DataFrame summary of saved spikes data.

    If a save file is found, a new DataFrame will not be generated (i.e.
    the loaded DataFrame will not be updated with new values).

    # Arguments
    * `load_df_json` -- bool or string indicating whether the function
    should check for an existing save.
      * `True` or `'current'` -- look for save file matching current
      version.
      * `'latest'` -- if any save files exist, load save with latest
      compatible version.
      * `'recent'` -- if any save files exist, load save with latest
      compatible version not exceeding current version.
      * A specific version number can be supplied as `str` and the
      function will only check for a save matching that version. If no
      save is found but a save matching the current version exists, the
      current save will be overwritten by a newly generated spikes_df.
      * `False` -- do not look for an existing save. As above, if a save
      matching the current version exists, it will be overwritten by a
      newly generated spikes_df.

    # Error Handling
    * If no file is matching `load_df_json` is found, or if
    `load_df_json` does not match any of the listed formats, a new
    DataFrame will be generated (overwriting a saved file if one exists
    for the current version).
    """
    # Build path and name for current version:
    # If a new DataFrame must be generated, it will be saved here. In
    # additional, if `load_df_json` is `True` or `'current'`, only this
    # location will be checked for an existing file.
    df_json_path = (constants.SAVE_PATHS['json_root'] +
                    constants.SAVE_PATHS['spikes_df'])
    # Define `spikes_df` so whether it is loaded can be checked later:
    spikes_df = None
    # Try to load `spikes_df` from file:
    if load_df_json:
        print("Trying to load `spikes_df` from file...")
        try:
            if load_df_json is True or load_df_json == 'current':
                # Build filename for current version:
                df_json_name = (constants.SPIKES_DF_JSON_NAME +
                                f'_{constants.VERSION}.json')
            elif load_df_json == 'latest' or load_df_json == 'recent':
                # Get list of compatible filenames:
                (compatible, _incompatible) = utils.get_json_filenames(
                    'spikes_df'
                )
                # Get list of versions, preserving index in `compatible`:
                versions = [
                    (i, classes.VersionNumber(filename)) for (i, filename) in
                    enumerate(compatible)
                ]
                # Get index and version of latest version number:
                latest = max(versions, key=lambda x: x[1])
                if load_df_json == 'recent':
                    # Cull versions later than current:
                    while latest[1] > utils.current_version:
                        versions.remove(latest)
                        latest = max(versions, key=lambda x: x[1])
                # Find filename at index of latest compatible version:
                df_json_name = compatible[latest[0]]
            else:
                # Build filename for specified version:
                df_json_name = (
                    constants.SPIKES_DF_JSON_NAME +
                    f'_{classes.VersionNumber(load_df_json)}.json'
                )
            # Attempt to load JSON at target filename:
            df_json = json.load(open(df_json_path+df_json_name))
            spikes_df = pd.DataFrame.from_dict(df_json)
            print(f"\n`spikes_df` loaded from file: {df_json_name}")
        except (IndexError, OSError):
            print("Unable to load `spikes_df`: no compatible file.\n")
        except ValueError:
            print("Unable to load `spikes_df`: invalid `load_df_json`.\n")
    if spikes_df is None:
        # If spikes_df was not loaded, build new DataFrame:
        print("Building `spikes_df` from saved spikes...")
        (spikes_files, _incompatible_files) = utils.get_json_filenames('spikes')
        all_trials = utils.load_spikes_trials(spikes_files)
        spikes_df = utils.tabulate_spikes(all_trials)
        print("\n`spikes_df` built.")
        # Save new DataFrame:
        df_json_name = (constants.SPIKES_DF_JSON_NAME +
                        f'_{constants.VERSION}.json')
        try:
            spikes_df.to_json(df_json_path+df_json_name,orient='records')
            print("\n`spikes_df` saved.")
        except OSError:
            pathlib.Path(df_json_path).mkdir(
                parents=True,
                exist_ok=True
            )
            spikes_df.to_json(df_json_path+df_json_name,orient='records')
    # Print completion message and return DataFrame:
    units = spikes_df['Unit ID'].unique()
    print(f'Contains {len(units)} units:')
    print('\n'.join([unit_id for unit_id in units]))
    return spikes_df


def calculate_ssr(
        spikes_df: pd.DataFrame
) -> pd.DataFrame:
    """Calculates mean spike rate per stimulation from saved spike data.
    
    For each trial in `spikes_df`, calculates mean spike rate per
    stimulation separated by experimental phase and stimulation type.
    Other factors, such as if spikes are uniformly distributed, are not
    considered.

    # Arguments
    * `spikes_df` -- DataFrame containing spike data.
    """
    freq_spikes = spikes_df.loc[spikes_df['Test'] == 'frequency']
    # ampl_spikes = spikes_df.loc[spikes_df['Test'] == 'amplitude']
    # nine_spikes = spikes_df.loc[spikes_df['Test'] == 'nine-one']
    # long_spikes = spikes_df.loc[spikes_df['Test'] == 'long-duration']
    freq_ssr = (freq.simple_spikerate_df(freq_spikes) if not freq_spikes.empty
                else pd.DataFrame())
    # ampl_ssr = (ampl.simple_spikerate_df(freq_spikes) if not ampl_spikes.empty
    #             else pd.DataFrame())
    # nine_ssr = (nine.simple_spikerate_df(freq_spikes) if not nine_spikes.empty
    #             else pd.DataFrame())
    # long_ssr = (long.simple_spikerate_df(freq_spikes) if not long_spikes.empty
    #             else pd.DataFrame())
    return freq_ssr
    # return pd.concat([freq_ssr, ampl_ssr, nine_ssr, long_ssr])


def plot_combined_ssr(
        ssr_df: pd.DataFrame,
        plot_title: str,
        save_figure: bool = False,
        filename: str | None = None
) -> None:
    """Plots calculated mean spike rates by experimental phase."""
    try:
        # Determine correct module to handle plotting:
        test = utils.unique(ssr_df['Test'])
        if test == 'frequency':
            ssr_plot_df = freq.simple_spikerate_plotdf(ssr_df)
            plot_ssr = freq.plot_single_ssr
        # elif test == 'amplitude':
        #     ssr_plot_df = ampl.simple_spikerate_plotdf(ssr_df)
        #     plot_ssr = ampl.plot_single_ssr
        # elif test == 'nine-one':
        #     ssr_plot_df = nine.simple_spikerate_plotdf(ssr_df)
        #     plot_ssr = nine.plot_single_ssr
        # elif test == 'long-duration':
        #     ssr_plot_df = long.simple_spikerate_plotdf(ssr_df)
        #     plot_ssr = long.plot_single_ssr
        else:
            raise KeyError(f"Invalid `Test` value: {test}.")
        # Plot subplots:
        f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=constants.FIGSIZE)
        plot_ssr(
            ssr_plot_df,
            'Conditioning',
            ax1
        )
        plot_ssr(
            ssr_plot_df,
            'Interleaved',
            ax2
        )
        plot_ssr(
            ssr_plot_df,
            'Recovery',
            ax3
        )
        # Prettify figure:
        plt.suptitle(plot_title)
        plt.tight_layout()
        # Save figure if specified:
        if save_figure:
            if not filename:
                filename = plot_title
            filename = filename.lower().replace(' ', '-').replace('_', '-')
            utils.save_plot(
                'ssr',
                f"ssr_{filename}_{test[0:4]}_{utils.current_version}"
            )
        return
    except AssertionError:
        # TODO generate four SSR plots together, one for each test
        raise NotImplementedError("Plotting multiple test types together.")


def plot_combined_raster(
        spikes_df: pd.DataFrame,
        save_figure: bool = False
 ) -> None:
    """Plots rasters for all recordings.
    
    A separate plot is generated for each test type, separated into
    three subplots by experimental phase, and with each trial on a
    separate row.
    """
    for test in spikes_df['Test'].unique():
        # Separate DataFrame by test type:
        test_df = spikes_df.loc[spikes_df['Test'] == test]
        plt.figure(figsize=constants.FIGSIZE)
        # Determine correct module:
        if test == 'frequency':
            plot_split = 'Test Frequency'
            plot_raster = freq.plot_single_raster
            ylabel_suffix = ' Hz'
        # elif test == 'amplitude':
        #     plot_split = 'Test Amplitude'
        #     plot_raster = ampl.plot_single_raster
        #     ylabel_suffix = ' μA'
        # else:
        #     plot_split = 'Test Stimulus'
        #     plot_raster = freq.plot_single_raster
        #     ylabel_suffix = ''
        # Separate, order, and plot subplots:
        plot_names = sorted(test_df[plot_split].unique())
        plots_count = len(plot_names)
        for plot_id, plot_name in enumerate(plot_names):
            plot_df = test_df.loc[test_df[plot_split] == plot_name]
            plt.subplot(plots_count, 1, plot_id+1)
            plot_raster(
                plot_df,
                'Trial ID'
            )
            plt.ylabel(f'{plot_name}{ylabel_suffix}')
        # Prettify figure:
        plt.suptitle(f"Rasters by {plot_split.lower()} and trial")
        plt.tight_layout()
        # Save figure if specified:
        if save_figure:
            utils.save_plot(
                'rasters',
                f"rasters_{test[0:4]}_[{plot_split}]-[Trial ID]_{utils.current_version}"
            )
    return
