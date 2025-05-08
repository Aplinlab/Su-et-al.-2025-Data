"""Functions specific to the analysis of frequency sweeps.

# Functions
## Spike detection
* `spikes_info` -- extracts spikes and epochs from imported data.

## Spike- and following-rate analysis
* `simple_spikerate_df` -- calculates mean spike rate per stimulation.
* `simple_spikerate_plotdf` -- organises SSR DataFrame for plotting.
* `plot_single_ssr` -- plots mean spike rate per stimulation.
* `plot_single_raster` -- plots raster with input variable on each row.
"""

import matplotlib.axes as axes
import matplotlib.pyplot as plt
import pandas as pd
import re
import statistics

from . import classes
from . import constants
from . import utils


# Regex pattern for parsing marker comments:
# It is compiled here once, rather than within a function, to reduce
# unnecessary computations.
p = re.compile(constants.SWEEP_REGEX)


########################## *SPIKE DETECTION* ###########################

def spikes_info(
        recording: classes.recordings.Recording,
        repetition: int,
        epoch_timing_ms: tuple[int | float, int | float],
        threshold_uV: int | float
) -> tuple[list[classes.spikes.SpikesTrial], list[classes.epochs.EpochsTrial]]:
    """Returns spikes and traces by epoch from extracted data.

    See `base.spikes_info()` docstring for more information.
    """
    # Define trigger thresholds:
    trigger_threshold_mech = (
        utils.trigger_value(recording.mech_triggers)*
        constants.TRIGGER_DETECTION_NOISE_WINDOW
    )
    trigger_threshold_elec = (
        utils.trigger_value(recording.elec_triggers)*
        constants.TRIGGER_DETECTION_NOISE_WINDOW
    )
    # Create empty lists to populate and return later:
    output_spikes = []
    output_epochs = []
    for marker in recording.markers:
        # Parse marker text and check that it fits expected format:
        m = p.search(marker.comment)
        try:
            assert m.group('testvar')=='frequency'
        # Print error message and skip marker if no match is found:
        except AttributeError:
            # AttributeError is raised if m is None (i.e. if regex
            # pattern didn't match)
            print(
                "Comment does not match the expected format for frequency "
                f"sweeps ({marker.comment})."
            )
            continue
        except AssertionError:
            print(
                "Comment indicates that this was not a frequency sweep, but "
                f"the file was labelled `freqsweep` ({marker.comment})."
            )
            continue
        # Extract frequencies from marker text:
        mech_cond_frequency = float(m.group('mechval'))
        elec_cond_frequency = float(m.group('elecval'))
        # Detect triggers and sort by experimental phase:
        triggers_mech = utils.triggers(
            recording.mech_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_mech
        )
        triggers_elec = utils.triggers(
            recording.elec_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_elec
        )
        triggers_by_phase = utils.separate_sweep_phases(
            'frequency',
            mech_cond_frequency,
            elec_cond_frequency,
            triggers_mech,
            triggers_elec
        )
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
                ) = utils.detect_spikes(
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
            'test_frequency': triggers_by_phase.stim_value,
            'test_amplitude': constants.FREQUENCY_SWEEP_CONDITIONING_AMPLITUDE,
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


################# *SPIKE- AND FOLLOWING-RATE ANALYSIS* #################

def simple_spikerate_df(
    spikes_df: pd.DataFrame
) -> pd.DataFrame:
    """Calculates mean spike rate per stimulation from saved spike data.

    See `base.calculate_ssr()` docstring for more information.
    """
    #!Effect of missing values (e.g. if specific frequencies are removed
    #!from analysis of a particular unit) is currently untested. Will
    #!likely raise `ZeroDivisionError` which can be handled to skip that
    #!frequency. However, must ensure that missing values do not end up
    #!in the output DataFrame as `0`!
    # Filter input to frequency tests if others are present:
    try:
        test = utils.unique(spikes_df['Test'])
        assert test == 'frequency'
    except AssertionError:
        spikes_df = spikes_df.loc[spikes_df['Test'] == 'frequency']
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
            # Identify relevant dict entry for this phase and stim:
            if epochs > 0:
                if phase == 'conditioning':
                    column = 'Conditioning'
                else:
                    column = f'{phase.capitalize()} {stim.capitalize()}'
            else:
                # Only conditioning phase may have no epochs:
                # This occurs for stim type other than test stim.
                if phase == 'conditioning':
                    pass
                else:
                    raise AssertionError(f"Unexpected empty phase ({phase_id}).")
            # Write relevant dict entry:
            simple_spikerates[column] = len(phase_spikes_df) / epochs
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
    return utils.typed_dataframe(
        simple_spikerate_df,
        constants.SIMPLE_SPIKERATE_TYPES,
        categories
    )


def simple_spikerate_plotdf(
        ssr_df: pd.DataFrame
) -> pd.DataFrame:
    """Returns a DataFrame for plotting simple spike rate analysis.
    
    Calculates mean of relevant values across rows in
    `simple_frequency_table` and arranges them so that they can be
    easily plotted using `simple_plot_df.plot()`.
    """
    #!As with previous function, effect of missing values (e.g. if
    #!specific frequencies are removed from analysis of a particular
    #!unit) is currently untested. Must ensure that whatever value ends
    #!up in the table can be differentiated from `0` and excluded from
    #!mean calculation!
    # Filter input to frequency tests if others are present:
    try:
        test = utils.unique(ssr_df['Test'])
        assert test == 'frequency'
    except AssertionError:
        ssr_df = ssr_df.loc[
            ssr_df['Test'] == 'frequency'
        ]
    # Identify frequencies present in the analysis:
    frequencies = ssr_df['Test Frequency'].unique()
    # Define DataFrame to be plotted:
    simple_plot_df = pd.DataFrame(
            index=frequencies,
            columns=constants.SSR_PLOT_COLUMNS
        ).astype(constants.SSR_PLOT_TYPES)
    for test_id in ssr_df['Test ID'].unique():
        # Filter analysis DataFrame to relevant rows:
        filtered_df = ssr_df.loc[ (ssr_df['Test ID'] == test_id)]
        stimulus = utils.unique(filtered_df['Test Stimulus'])
        frequency = utils.unique(filtered_df['Test Frequency'])
        # Filter matching rows and calculate mean:
        # Conditioning phase is handled first, separately.
        simple_plot_df.loc[
            frequency,
            f'SSR Conditioning - {stimulus.capitalize()}'
        ] = statistics.mean(filtered_df['SSR Conditioning'])
        for input_column in filtered_df.columns[-4:]:
            simple_plot_df.loc[
                frequency,
                f'{input_column} Post-{stimulus.capitalize()}'
            ] = statistics.mean(filtered_df[input_column])
    # Sort and return DataFrame
    simple_plot_df.sort_index(inplace=True)
    return simple_plot_df


def plot_single_ssr(
        ssr_plot_df: pd.DataFrame,
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
        ssr_plot_df = ssr_plot_df.loc[:, ssr_plot_df.columns.str.contains(phase)]
    except TypeError:
        pass
    for (series_name, series) in ssr_plot_df.items():
        ax.plot(series, label=series_name)
    # Prettify figure:
    plt.setp(
        ax,
        xticks=list(ssr_plot_df.index),
        xlabel="Conditioning phase frequency (Hz)",
        ylabel="Mean spike rate per stimulation",
        title=phase
    )
    ax.legend(
        prop={'size': constants.SSR_PLOT_LEGEND_SIZE}
    )
    ax.set_ylim(constants.SSR_PLOT_YLIM)
    return


def plot_single_raster(
        spikes_df: pd.DataFrame,
        row_split: str,
) -> None:
    """Plots raster with rows separated by `row_split`."""
    # Filter input to frequency tests if others are present:
    try:
        test = utils.unique(spikes_df['Test'])
        assert test == 'frequency'
    except AssertionError:
        spikes_df = spikes_df.loc[spikes_df['Test'] == 'frequency']
    # Get list of rows to be drawn:
    row_names = spikes_df[row_split].unique()
    # Define initial y-coordinate for mechanical and electrical rows:
    y_mech = 0
    y_elec = len(row_names) / 2
    for row_name in row_names:
        row_df = spikes_df.loc[spikes_df[row_split] == row_name]
        for (_i, spike) in row_df.iterrows():
            # Get variables for plotting:
            phase = spike['Phase']
            test_stim = spike['Test Stimulus']
            epoch_number = spike['Epoch Number']
            point_colour = ('orange' if spike['Epoch Stimulus'] == 'mechanical' else
                      ('sky' if spike['Epoch Stimulus'] == 'electrical' else
                       'black'))
            # Get x-value (latency) and make necessary adjustments:
            if phase == 'conditioning':
                x = epoch_number/spike['Test Frequency']
            else:
                if phase == 'interleaved':
                    frequency = constants.INTERLEAVED_EPOCHS_FREQUENCY * 2
                elif phase == 'recovery':
                    frequency = constants.RECOVERY_EPOCHS_FREQUENCY * 20
                x_diff = int(spike['Epoch Stimulus'] == 'mechanical')
                x = (constants.RASTER_PHASE_SPACING[phase] +
                     (epoch_number * 2 - x_diff) / frequency)
            x += (spike['Latency (ms)'] / 1000)
            # Get y-value (no point will be plotted if y is None)
            y = (y_mech if test_stim == 'mechanical' else
                 (y_elec if test_stim == 'electrical' else None))
            # Plot point:
            plt.scatter(
                x,
                y,
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
