"""Functions for quantifying and plotting spike response rates.

# Functions
## User-exposed functions
See package docstring for function summaries.
* `calculate_ssr`
* `plot_combined_raster`
* `plot_combined_ssr`
* `spikes_table`

## Backend functions
* `get_json_filename` -- gets list of saved JSON files matching
specified type.
* `load_spikes_trials` -- loads retrieved list of saved spike files.
* `plot_single_raster` -- plots raster with input variable on each row.
* `plot_single_ssr` -- plots mean spike rate per stimulation.
* `simple_spikerate_plotdf` -- organises SSR DataFrame for plotting.
* `tabulate_spikes` -- creates DataFrame summary of saved spike data.
* `typed_dataframe` -- applies types to DataFrame columns (handles both
categorical and non-categorical types).
"""

import json
import math
import matplotlib.axes as axes
import matplotlib.pyplot as plt
from os import walk
import pandas as pd
import pathlib
import statistics

from . import classes
from . import constants
from . import utils


####################### *USER-EXPOSED FUNCTIONS* #######################

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
                (compatible, _incompatible) = get_json_filenames(
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
        (spikes_files, _incompatible_files) = get_json_filenames('spikes')
        all_trials = load_spikes_trials(spikes_files)
        spikes_df = tabulate_spikes(all_trials)
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
    #!Effect of missing values (e.g. if specific frequencies are removed
    #!from analysis of a particular unit) is currently untested. Will
    #!likely raise `ZeroDivisionError` which can be handled to skip that
    #!frequency. However, must ensure that missing values do not end up
    #!in the output DataFrame as `0`!
    # Define output DataFrame:
    simple_spikerate_df = pd.DataFrame()
    for trial_id in spikes_df['Trial ID'].unique():
        test_spikes_df = spikes_df.loc[spikes_df['Trial ID'] == trial_id]
        # Get columns to be retained:
        animal_id = utils.unique(test_spikes_df['Animal ID'])
        sex = utils.unique(test_spikes_df['Sex'])
        position = utils.unique(test_spikes_df['Position'])
        unit_id = utils.unique(test_spikes_df['Unit ID'])
        unit_type = utils.unique(test_spikes_df['Unit Type'])
        test = utils.unique(test_spikes_df['Test'])
        stimulus = utils.unique(test_spikes_df['Test Stimulus'])
        frequency = utils.unique(test_spikes_df['Test Frequency'])
        amplitude = utils.unique(test_spikes_df['Test Amplitude'])
        test_id = utils.unique(test_spikes_df['Test ID'])
        repetition = utils.unique(test_spikes_df['Repetition'])
        # Prepare spike rates dict:
        simple_spikerates = {
            'Conditioning': 0.0,
            'Interleaved Mechanical': 0.0,
            'Interleaved Electrical': 0.0,
            'Recovery Mechanical': 0.0,
            'Recovery Electrical': 0.0
        }
        for phase_id in test_spikes_df['Phase ID'].unique():
            phase_spikes_df = test_spikes_df.loc[
                spikes_df['Phase ID'] == phase_id
            ]
            # Get columns to be retained:
            epochs = utils.unique(phase_spikes_df['Total Epochs'])
            phase = utils.unique(phase_spikes_df['Phase'])
            stim = utils.unique(phase_spikes_df['Epoch Stimulus'])
            # Identify and write relevant dict entry for phase and stim:
            if epochs > 0:
                if phase == 'conditioning':
                    #!This works for all test types except nine-one, as
                    #!that uniquely contains both stim types during the
                    #!conditioning phase
                    column = 'Conditioning'
                else:
                    column = f'{phase.capitalize()} {stim.capitalize()}'
                simple_spikerates[column] = len(phase_spikes_df) / epochs
            else:
                # Only conditioning phase may have no epochs:
                # This occurs for stim type other than test stim.
                if phase == 'conditioning':
                    pass
                else:
                    raise AssertionError(
                        f"Unexpected empty phase ({phase_id})."
                    )
        # Append results to output dict:
        simple_spikerate_df = pd.concat([
            simple_spikerate_df,
            pd.DataFrame([[
                animal_id,
                sex,
                position,
                unit_id,
                unit_type,
                test,
                stimulus,
                frequency,
                amplitude,
                test_id,
                repetition,
                trial_id,
                simple_spikerates['Conditioning'],
                simple_spikerates['Interleaved Mechanical'],
                simple_spikerates['Interleaved Electrical'],
                simple_spikerates['Recovery Mechanical'],
                simple_spikerates['Recovery Electrical']
            ]], columns=constants.SIMPLE_SPIKERATE_COLUMNS)
        ], ignore_index=True)
    # Apply dtypes to output dict:
    categories = {
        'Test': constants.TEST_CODES.values(),
        'Test Stimulus': constants.STIMULATION_TYPES
    }
    return typed_dataframe(
        simple_spikerate_df,
        constants.SIMPLE_SPIKERATE_TYPES,
        categories
    )


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
        if test == 'frequency' or test == 'amplitude':
            plot_split = f'Test {test.capitalize()}'
            plot_raster = plot_single_raster
        else:
            print(
                "Rasters may only be plotted for frequency sweeps and "
                "amplitude sweeps."
            )
            continue
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
            plt.ylabel(f'{plot_name} {constants.TEST_UNITS[test]}')
        # Prettify figure:
        plt.suptitle(f"Rasters by {plot_split.lower()} and trial")
        plt.tight_layout()
        # Save figure if specified:
        if save_figure:
            utils.save_plot(
                'rasters',
                (f"rasters_{test[0:4]}_[{plot_split}]-[Trial ID]_" +
                 str(utils.current_version))
            )
    return


def plot_combined_rvt(
        spikes_df: pd.DataFrame,
        save_figure: bool = False
) -> None:
    return


def plot_combined_ssr(
        ssr_df: pd.DataFrame,
        plot_title: str,
        save_figure: bool = False,
        filename: str | None = None
) -> None:
    """Plots calculated mean spike rates by experimental phase.
    
    # Error Handling
    If `ssr_df` contains test types which cannot be plotted by this
    function, prints an error message and skips those test types.
    """
    for test in ssr_df['Test'].unique():
        if test == 'frequency' or test == 'amplitude':
            ssr_plot_df = simple_spikerate_plotdf(
                ssr_df.loc[ssr_df['Test'] == test],
                f'Test {test.capitalize()}'
            )
            plot_ssr = plot_single_ssr
        else:
            print(
                "Simple spikerate should only be plotted for frequency "
                f"sweeps or amplitude sweeps (input test type: {test})."
            )
            continue
        # Plot subplots:
        _f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=constants.FIGSIZE)
        for (phase, ax) in {
            'Conditioning': ax1,
            'Interleaved': ax2,
            'Recovery': ax3
        }.items():
            plot_ssr(
                ssr_plot_df,
                test,
                phase,
                ax
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


########################## *BACKEND FUNCTIONS* #########################

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
        (compatible, message) = file_version.iscompatible(utils.current_version)
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


def plot_single_raster(
        spikes_df: pd.DataFrame,
        row_split: str,
) -> None:
    """Plots raster with rows separated by `row_split`."""
    # Ensure only one valid test type is present in input DataFrame:
    try:
        test = utils.unique(spikes_df['Test'])
        assert test == 'frequency' or test == 'amplitude'
    except AssertionError:
        print(
            "Before running this function, filter input DataFrame to "
            "contain only frequency sweep trials or only amplitude "
            "sweep trials."
        )
        raise
    # Get list of rows to be drawn:
    row_names = spikes_df[row_split].unique()
    # Define initial y-coordinate for mechanical and electrical rows:
    y_mech = 0
    y_elec = len(row_names) / 2
    for row_name in row_names:
        row_df = spikes_df.loc[spikes_df[row_split] == row_name]
        test_stim = utils.unique(row_df['Test Stimulus'])
        for (_i, spike) in row_df.iterrows():
            # Get variables for plotting:
            phase = spike['Phase']
            epoch_number = spike['Epoch Number']
            point_colour = (
                'orange' if spike['Epoch Stimulus'] == 'mechanical' else (
                    'sky' if spike['Epoch Stimulus'] == 'electrical' else
                    'black'
                )
            )
            # Get x-value (latency) and make necessary adjustments:
            if phase == 'conditioning':
                x = epoch_number/spike['Test Frequency']
            else:
                if phase == 'interleaved':
                    frequency = constants.INTERLEAVED_FREQUENCY * 2
                elif phase == 'recovery':
                    frequency = constants.RECOVERY_FREQUENCY * 20
                else:
                    raise KeyError(
                        f"Experimental phase not recognised ({phase})."
                    )
                x_diff = int(spike['Epoch Stimulus'] == 'mechanical')
                #!`RASTER_PHASE_SPACING` is only valid for frequency
                #!sweeps and amplitude sweeps
                x = (constants.RASTER_PHASE_SPACING[phase] +
                     (epoch_number * 2 - x_diff) / frequency)
            x += (spike['Latency (ms)'] / 1000)
            # Get y-value (no point will be plotted if y is None)
            y = (y_mech if test_stim == 'mechanical' else
                 (y_elec if test_stim == 'electrical' else None))
            # Plot point:
            plt.scatter(
                x,
                y, # type: ignore
                constants.RASTER_POINT_SIZE,
                constants.PALETTE[point_colour]
            )
        # Iterate relevant y-coordinate:
        if test_stim == 'mechanical':
            y_mech += 1
        elif test_stim == 'electrical':
            y_elec += 1
    # Remove axis ticks:
    plt.xticks([])
    plt.yticks([])
    return


def plot_single_rvt(
        spikes_df: pd.DataFrame,
        test: str,
        stim_type: str,
        phase: str,
        ax: axes.Axes | None = None
) -> None:
    # TODO write docstring and comments
    # ! Completely untested
    if ax is None:
        ax = plt.gca()
    spikes_df_stim = spikes_df.loc[spikes_df['Test Stimulus'] == stim_type]
    spikes_df_phase = spikes_df_stim.loc[spikes_df_stim['Phase'] == phase]
    if test == 'frequency' or test == 'amplitude':
        filter_col = f'Test {test.capitalize()}'
        unit = f' {constants.TEST_UNITS[test]}'
    else:
        filter_col = 'Test Stimulus'
        unit = ''
    filter_values = spikes_df_phase[filter_col].unique()
    for value in filter_values:
        spikes_df_filtered = spikes_df_phase.loc[
            spikes_df_phase[filter_col] == value
        ]
        frequency = utils.unique(spikes_df_filtered['Test Frequency'])
        max_epoch_count = max(spikes_df_filtered['Total Epochs'])
        bin_width = min(math.floor(max_epoch_count/10), 100)
        bins = [x+bin_width/2 for x in range(
            0,
            math.floor(max_epoch_count/bin_width)*bin_width,
            bin_width
        )]
        response_rates = [[] for _x in bins]
        for trial_id in spikes_df_filtered['Trial ID'].unique():
            spikes_df_trial = spikes_df_filtered.loc[
                spikes_df_filtered['Trial ID'] == trial_id
            ]
            epoch_count = utils.unique(spikes_df_trial['Total Epochs'])
            for (i, x) in enumerate(range(
                0,
                math.floor(epoch_count/bin_width)*bin_width,
                bin_width
            )):
                response_rates[i].append(len(spikes_df_trial.loc[
                    x<=spikes_df_trial['Epoch Number']<x+bin_width
                ])/bin_width)
        mean_rates = [statistics.mean(x) for x in response_rates]
        ax.plot(
            [x/frequency for x in bins],
            mean_rates,
            label=f'{value}{unit}'
        )
    plt.setp(
        ax,
        xlabel="Time (s)",
        ylabel="Mean spike rate per bin",
        title=phase
    )
    ax.legend(
        prop={'size': constants.QUANTIFICATION_LEGEND_SIZE}
    )
    ax.set_ylim(
        -constants.QUANTIFICATION_YLIM_BORDER,
        1+constants.QUANTIFICATION_YLIM_BORDER
    )
    return


def plot_single_ssr(
        ssr_plot_df: pd.DataFrame,
        test: str,
        phase: str | None = None,
        ax: axes.Axes | None = None
) -> None:
    """Plots the specified columns (phases) of `simple_plot_df`.
    
    If no `phase` argument is provided, or there is no data matching the
    provided `phase`, every column will be plotted.
    """
    # Plot specified columns:
    if ax is None:
        ax = plt.gca()
    try:
        ssr_plot_df = ssr_plot_df.loc[
            :,
            ssr_plot_df.columns.str.contains(phase) # type: ignore
        ]
    except TypeError:
        pass
    max_point = 0
    for (series_name, series) in ssr_plot_df.items():
        ax.plot(series, label=series_name)
        max_point = max(max_point, series.max())
    # Prettify figure:
    plt.setp(
        ax,
        xticks=list(ssr_plot_df.index),
        xlabel=f"Conditioning phase {test} ({constants.TEST_UNITS[test]})",
        ylabel="Mean spike rate per stimulation",
        title=phase
    )
    ax.legend(
        prop={'size': constants.QUANTIFICATION_LEGEND_SIZE}
    )
    ax.set_ylim(
        -constants.QUANTIFICATION_YLIM_BORDER,
        max_point+constants.QUANTIFICATION_YLIM_BORDER
    )
    return


def simple_spikerate_plotdf(
        ssr_df: pd.DataFrame,
        index_col: str
) -> pd.DataFrame:
    """Returns a DataFrame for plotting simple spike rate analysis.
    
    Calculates mean of relevant values in `ssr_df` and arranges them for
    easy plotting.

    # Arguments
    * `ssr_df` -- input DataFrame following format output by
    `base.calculate_ssr()`.
    """
    #!As with previous function, effect of missing values (e.g. if
    #!specific frequencies are removed from analysis of a particular
    #!unit) is currently untested. Must ensure that whatever value ends
    #!up in the table can be differentiated from `0` and excluded from
    #!mean calculation!
    # Check that only one test type is present in input DataFrame:
    try:
        utils.unique(ssr_df['Test'])
    except AssertionError:
        print(
            "Before running this function, filter input DataFrame to "
            "contain only one test type."
        )
        raise
    # Identify test variables present in the analysis:
    index = ssr_df[index_col].unique()
    # Define DataFrame to be plotted:
    simple_plot_df = pd.DataFrame(
            index=index,
            columns=constants.SSR_PLOT_COLUMNS
        ).astype(constants.SSR_PLOT_TYPES)
    for test_id in ssr_df['Test ID'].unique():
        # Filter analysis DataFrame to relevant rows:
        filtered_df = ssr_df.loc[ (ssr_df['Test ID'] == test_id)]
        stimulus = utils.unique(filtered_df['Test Stimulus'])
        index = utils.unique(filtered_df[index_col])
        # Filter matching rows and calculate mean:
        # Conditioning phase is handled first, separately.
        simple_plot_df.loc[
            index,
            f'SSR Conditioning - {stimulus.capitalize()}'
        ] = statistics.mean(filtered_df['SSR Conditioning'])
        for input_column in filtered_df.columns[-4:]:
            simple_plot_df.loc[
                index,
                f'{input_column} Post-{stimulus.capitalize()}'
            ] = statistics.mean(filtered_df[input_column])
    # Sort and return DataFrame
    simple_plot_df.sort_index(inplace=True)
    return simple_plot_df


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
