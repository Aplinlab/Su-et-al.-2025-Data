"""Functions specific to the analysis of frequency sweeps."""

# TODO Allow user to select which phases are plotted using
# TODO `phase` argument (line 289)

import json
import matplotlib.pyplot as plt
import pandas as pd
import re
import statistics

from . import base
from . import classes
from . import constants


# Regex pattern for parsing marker comments:
# It is compiled here once, rather than within a function, to reduce
# unnecessary computations.
p = re.compile(constants.SWEEP_REGEX)


########################## *SPIKE DETECTION* ###########################

def freqsweep_spikes(
        recording: classes.recordings.Recording,
        epoch_timing_ms: tuple[int | float, int | float],
        threshold_uV: int | float
) -> tuple[list[classes.spikes.SpikesTrial], list[classes.epochs.EpochsTrial]]:
    """Returns information about spikes in a frequency sweep.
    
    # Arguments
    * `recording` -- `Recording` object to be analysed.
    * `threshold_uV` -- threshold above which spikes should be detected,
    in microvolts.
    * `epoch_timing_ms_mech` and `epoch_timing_ms_elec` -- tuples of two
    numeric values describing the timing window during which spikes may
    occur. The first value of each tuple represents start time and the
    second value stop time, in milliseconds after each stimulus onset).
    Separate tuples are provided for mechanical and electrical
    stimulation, respectively.

    # Error Handling
    * Upon encountering a comment which does not match the expected
    format for frequency sweeps, prints an error message containing the
    comment text and moves onto the next comment.
    * Raises AssertionError if the number of epochs extracted from any
    trial does not match the number of triggers present in that trial.
    """
    # Define trigger thresholds:
    trigger_threshold_mech = (
        base.trigger_value(recording.mech_triggers)*
        constants.TRIGGER_DETECTION_NOISE_WINDOW
    )
    trigger_threshold_elec = (
        base.trigger_value(recording.elec_triggers)*
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
        triggers_mech = base.triggers(
            recording.mech_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_mech
        )
        triggers_elec = base.triggers(
            recording.elec_triggers[marker.start_sample:marker.end_sample],
            trigger_threshold_elec
        )
        triggers_by_phase = base.separate_sweep_phases(
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
                ) = base.detect_spikes(
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
            'test_freqency': triggers_by_phase.stim_value,
            'test_amplitude': constants.FREQUENCY_SWEEP_CONDITIONING_AMPLITUDE
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
    return output_spikes, output_epochs


################## *SPIKE- & FOLLOWING-RATE ANALYSIS* ##################

def simple_spikerate_df(json_spike_filenames: list[str]) -> pd.DataFrame:
    """Calculates mean spike rate per stimulation from saved spike data.
    
    For each JSON file in `json_spike_filenames`, calculates mean spike
    rate per stimulation separated by trial, experimental phase, and
    stimulation type. Other factors, such as if spikes are uniformly
    distributed in time, are not calculated.

    # Arguments
    * `json_spike_filenames` -- list of filenames pointing to JSON files
    containing spike data.
    """
    #!Effect of missing values (e.g. if specific frequencies are removed
    #!from analysis of a particular unit) is currently untested. Will
    #!likely raise `ZeroDivisionError` which can be handled to skip that
    #!frequency. However, must ensure that missing values do not end up
    #!in the output DataFrame as `0`!
    # Define output DataFrame:
    simple_frequency_table = base.typed_dataframe(
        constants.SIMPLE_FREQUENCY_COLUMNS,
        constants.SIMPLE_FREQUENCY_TYPES,
        constants.SIMPLE_FREQUENCY_CATEGORIES
    )
    for spike_filename in json_spike_filenames:
        # Construct path from which JSON file is to be read:
        save_file = constants.SPIKES_JSON_FOLDER + spike_filename
        # Read JSON file and convert to list of `SpikesTrial` objects:
        spikes_dict = json.load(open(save_file))
        # Check that version is correct:
        # TODO Handle incorrect version numbers
        assert spikes_dict['version'] == constants.VERSION
        spikes = [classes.spikes.SpikesTrial.from_dict(dictionary)
                  for dictionary in spikes_dict['spikes']]
        for trial in spikes:
            # Calculate mean spike rate per trial in conditioning phase:
            if trial.conditioning.epochs_mech:
                cond_rate = (len(trial.conditioning.mechanical)/
                             trial.conditioning.epochs_mech)
            elif trial.conditioning.epochs_elec:
                cond_rate = (len(trial.conditioning.electrical)/
                             trial.conditioning.epochs_elec)
            # Calculate mean spike rate per trial in other phases:
            mech_itlv_rate = (len(trial.interleaved.mechanical)/
                              trial.interleaved.epochs_mech)
            elec_itlv_rate = (len(trial.interleaved.electrical)/
                              trial.interleaved.epochs_elec)
            mech_rcvr_rate = (len(trial.recovery.mechanical)/
                              trial.recovery.epochs_mech)
            elec_rcvr_rate = (len(trial.recovery.electrical)/
                              trial.recovery.epochs_elec)
            # Append results to output DataFrame:
            simple_frequency_table = pd.concat([
                simple_frequency_table,
                pd.DataFrame([[
                    trial.animal_id,
                    trial.position,
                    trial.test_freqency,
                    trial.test_stim.capitalize(),
                    cond_rate,
                    mech_itlv_rate,
                    elec_itlv_rate,
                    mech_rcvr_rate,
                    elec_rcvr_rate
                ]], columns=constants.SIMPLE_FREQUENCY_COLUMNS)
            ], ignore_index=True)
    # Return output DataFrame:
    return simple_frequency_table


def simple_spikerate_plotdf(
        simple_frequency_table: pd.DataFrame
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
    # Identify frequencies present in the analysis:
    frequencies = simple_frequency_table[
        'Conditioning Frequency (Hz)'
    ].unique()
    # Define DataFrame to be plotted:
    simple_plot_df = pd.DataFrame(
            index=frequencies,
            columns=constants.SIMPLE_PLOT_COLUMNS
        ).astype(constants.SIMPLE_PLOT_TYPES)
    
    for stimulus in constants.SIMPLE_FREQUENCY_CATEGORIES[
        'Conditioning Stimulus'
    ]:
        for frequency in frequencies:
            # Filter analysis DataFrame to relevant rows:
            filtered_df = simple_frequency_table.loc[
                (simple_frequency_table['Conditioning Stimulus'] == stimulus) &
                (simple_frequency_table['Conditioning Frequency (Hz)'] ==
                 frequency)
            ]
            # Filter matching rows and calculate mean:
            # Conditioning phase is handled first, separately.
            simple_plot_df.loc[
                frequency,
                f'FR Conditioning - {stimulus}'
            ] = statistics.mean(filtered_df['FR Conditioning'])
            for input_column in constants.SIMPLE_FREQUENCY_COLUMNS[5:]:
                simple_plot_df.loc[
                    frequency,
                    f'{input_column} Post-{stimulus}'
                ] = statistics.mean(filtered_df[input_column])
    # Sort and return DataFrame
    simple_plot_df.sort_index(inplace=True)
    return simple_plot_df


def plot_simple_spikerate(
        simple_plot_df: pd.DataFrame,
        phase: str | None = None,
        save_figure = False
) -> None:
    """Plots the specified columns (phases) of `simple_plot_df`.
    
    If no `phase` argument is provided, or there is no data matching the
    provided `phase`, every column will be plotted.
    """
    # Plot specified columns:
    try:
        simple_plot_df.loc[:, simple_plot_df.columns.str.contains(phase)].plot()
    except TypeError:
        simple_plot_df.plot()
    # Prettify figure:
    plt.xticks(list(simple_plot_df.index))
    plt.xlabel("Conditioning phase frequency (Hz)")
    plt.ylabel("Mean spike rate per stimulation")
    plt.legend(
        prop={'size': constants.SIMPLE_PLOT_LEGEND_SIZE},
        loc=constants.SIMPLE_PLOT_LEGEND_LOC
    )
    ax = plt.gca()
    ax.set_ylim(0)
    # Save figure if specified:
    if save_figure:
        filename = (f'simple_spikerate-{phase}.pdf' if phase else
                    'simple_spikerate-all.pdf')
        plt.savefig(constants.PLOTS_FOLDER + filename)
    return
