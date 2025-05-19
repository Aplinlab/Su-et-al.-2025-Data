"""Functions specific to the analysis of frequency and amplitude sweeps.

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
    assert recording.test == 'frequency' or recording.test == 'amplitude'
    # Define trigger thresholds:
    if recording.test == 'amplitude':
        trigger_correction_mech = max_amp("mech", recording.markers)
        trigger_correction_elec = max_amp("elec", recording.markers)
    else:
        trigger_correction_mech = 1
        trigger_correction_elec = 1
    trigger_threshold_mech = (
        utils.trigger_value(recording.mech_triggers, trigger_correction_mech)*
        constants.TRIGGER_DETECTION_NOISE_WINDOW
    )
    trigger_threshold_elec = (
        utils.trigger_value(recording.elec_triggers, trigger_correction_elec)*
        constants.TRIGGER_DETECTION_NOISE_WINDOW
    )
    # Create empty lists to populate and return later:
    output_spikes = []
    output_epochs = []
    for marker in recording.markers:
        # Parse marker text and check that it fits expected format:
        m = p.search(marker.comment)
        try:
            assert m.group('testvar') == recording.test # type: ignore
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
        mech_value = float(m.group('mechval')) # type: ignore
        elec_value = float(m.group('elecval')) # type: ignore
        if recording.test == 'frequency':
            adjustment = 1
        elif recording.test == 'amplitude':
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
        # Detect triggers:
        triggers_mech = utils.triggers(
            recording.mech_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_mech * adjustment
        )
        triggers_elec = utils.triggers(
            recording.elec_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_elec * adjustment
        )
        # Sort triggers by experimental phase using frequencies:
        triggers_by_phase = utils.separate_sweep_phases(
            recording.test,
            mech_value,
            elec_value,
            triggers_mech,
            triggers_elec
        )
        if recording.test == 'frequency':
            test_frequency = triggers_by_phase.stim_value
            test_amplitude = constants.FREQUENCY_SWEEP_CONDITIONING_AMPLITUDE
        elif recording.test == 'amplitude':
            test_frequency = constants.AMPLITUDE_SWEEP_CONDITIONING_FREQUENCY
            test_amplitude = triggers_by_phase.stim_value
        else:
            # Cannot be reached due to assertion at start of function,
            # but is included to fix `reportPossiblyUnboundError` later
            continue
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


def max_amp(stim_type: str, markers: list[classes.recordings.Marker]):
    max_amp = 1
    for marker in markers:
        m = p.search(marker.comment)
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
                "Trial which is not amplitude sweep has been passed to "
                f"amplitude sweep pipeline ({marker.comment})."
            )
    return max_amp


################# *SPIKE- AND FOLLOWING-RATE ANALYSIS* #################

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
    for (series_name, series) in ssr_plot_df.items():
        ax.plot(series, label=series_name)
    # Prettify figure:
    plt.setp(
        ax,
        xticks=list(ssr_plot_df.index),
        xlabel=f"Conditioning phase {test} ({constants.TEST_UNITS[test]})",
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
                    frequency = constants.INTERLEAVED_EPOCHS_FREQUENCY * 2
                elif phase == 'recovery':
                    frequency = constants.RECOVERY_EPOCHS_FREQUENCY * 20
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
