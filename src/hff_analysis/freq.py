"""Functions specific to the analysis of frequency sweeps."""

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
            'test_frequency': triggers_by_phase.stim_value,
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

def simple_spikerate_df(
    spikes_df: pd.DataFrame
) -> pd.DataFrame:
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
    
    test, *test_empty = spikes_df['Test'].unique()
    try:
        assert len(test_empty) == 0 and test == 'Frequency'
    except AssertionError:
        spikes_df = spikes_df.loc[spikes_df['Test'] == 'Frequency']

    # Define output DataFrame:
    simple_spikerate_df = pd.DataFrame()
    for trial_id in spikes_df['Trial ID'].unique():
        test_spikes_df = spikes_df.loc[spikes_df['Trial ID'] == trial_id]

        animal_id, *animal_empty = test_spikes_df['Animal ID'].unique()
        sex, *sex_empty = test_spikes_df['Sex'].unique()
        position, *position_empty = test_spikes_df['Position'].unique()
        unit_id, *unit_id_empty = test_spikes_df['Unit ID'].unique()
        unit_type, *unit_type_empty = test_spikes_df['Unit Type'].unique()
        stimulus, *stimulus_empty = test_spikes_df['Test Stimulus'].unique()
        frequency, *frequency_empty = test_spikes_df['Test Frequency'].unique()
        amplitude, *amplitude_empty = test_spikes_df['Test Amplitude'].unique()
        test_id, *test_id_empty = test_spikes_df['Test ID'].unique()
        repetition, *repetition_empty = test_spikes_df['Repetition'].unique()
        assert (
            len(animal_empty) == 0 and
            len(sex_empty) == 0 and
            len(position_empty) == 0 and
            len(unit_id_empty) == 0 and
            len(unit_type_empty) == 0 and
            len(stimulus_empty) == 0 and
            len(frequency_empty) == 0 and
            len(amplitude_empty) == 0 and
            len(test_id_empty) == 0 and
            len(repetition_empty) == 0
        )
        
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

            epochs, *epochs_empty = phase_spikes_df['Total Epochs'].unique()
            phase, *phase_empty = phase_spikes_df['Phase'].unique()
            stim, *stim_empty = phase_spikes_df['Epoch Stimulus'].unique()
            assert (
                len(epochs_empty) == 0 and
                len(phase_empty) == 0 and
                len(stim_empty) == 0
            )

            if epochs > 0:
                if phase == 'Conditioning':
                    column = 'Conditioning'
                else:
                    column = f'{phase} {stim}'
            else:
                if phase == 'Conditioning':
                    pass
                else:
                    raise AssertionError(f"Unexpected empty phase ({phase_id}).")
            
            simple_spikerates[column] = len(phase_spikes_df) / epochs

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
    
    return base.typed_dataframe(
        simple_spikerate_df,
        constants.SIMPLE_SPIKERATE_TYPES,
        constants.SIMPLE_SPIKERATE_CATEGORIES
    )


def simple_spikerate_plotdf(
        simple_spikerate_df: pd.DataFrame
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
    
    test, *test_empty = simple_spikerate_df['Test'].unique()
    try:
        assert len(test_empty) == 0 and test == 'Frequency'
    except AssertionError:
        simple_spikerate_df = simple_spikerate_df.loc[
            simple_spikerate_df['Test'] == 'Frequency'
        ]

    # Identify frequencies present in the analysis:
    frequencies = simple_spikerate_df[
        'Test Frequency'
    ].unique()
    # Define DataFrame to be plotted:
    simple_plot_df = pd.DataFrame(
            index=frequencies,
            columns=constants.SSR_PLOT_COLUMNS
        ).astype(constants.SSR_PLOT_TYPES)

    for test_id in simple_spikerate_df['Test ID'].unique():
        # Filter analysis DataFrame to relevant rows:
        filtered_df = simple_spikerate_df.loc[
            (simple_spikerate_df['Test ID'] == test_id)
        ]
        stimulus, *stimulus_empty = filtered_df['Test Stimulus'].unique()
        frequency, *frequency_empty = filtered_df['Test Frequency'].unique()
        assert (
            len(stimulus_empty) == 0 and
            len(frequency_empty) == 0
        )
        # Filter matching rows and calculate mean:
        # Conditioning phase is handled first, separately.
        simple_plot_df.loc[
            frequency,
            f'SSR Conditioning - {stimulus}'
        ] = statistics.mean(filtered_df['SSR Conditioning'])
        for input_column in filtered_df.columns[-4:]:
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
        plot_title = None,
        save_figure = False,
        filename_suffix = None,
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
    plt.title(plot_title)
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
        if filename_suffix:
            filename = f'simple_spikerate-{filename_suffix}.pdf'
        elif phase:
            filename = f'simple_spikerate-{phase}.pdf'
        else:
            filename = 'simple_spikerate-all.pdf'
        base.save_plot('ssr', filename)
    return


def plot_raster(
        spikes_df: pd.DataFrame,
        plot_id: int,
        plot_count: int,
        row_criteria: str,
        save_figure: bool = False
) -> None:
    plt.subplot(plot_count, 1, plot_id)
    row_ids = spikes_df[row_criteria].unique()
    y_mech = 0
    y_elec = len(row_ids) / 2
    for row_id in row_ids:
        row_df = spikes_df.loc[spikes_df[row_criteria] == row_id]
        for epoch_id in row_df['Epoch ID'].unique():
            epoch_df = row_df.loc[row_df['Epoch ID'] == epoch_id]

            phase, *phase_empty = epoch_df['Phase'].unique()
            test_stim, *test_stim_empty = epoch_df['Test Stimulus'].unique()
            test_frequency, *test_frequency_empty = epoch_df['Test Frequency'].unique()
            epoch_stim, *epoch_stim_empty = epoch_df['Epoch Stimulus'].unique()
            epoch_number, *epoch_number_empty = epoch_df['Epoch Number'].unique()
            assert(
                len(phase_empty) == 0 and
                len(test_stim_empty) == 0 and
                len(test_frequency_empty) == 0 and
                len(epoch_stim_empty) == 0 and
                len(epoch_number_empty) == 0
            )

            epoch_number = float(epoch_number)

            if phase == 'Conditioning':
                x = epoch_number/test_frequency
            else:
                if phase == 'Interleaved':
                    frequency = constants.INTERLEAVED_EPOCHS_FREQUENCY * 2
                elif phase == 'Recovery':
                    frequency = constants.RECOVERY_EPOCHS_FREQUENCY * 20
                x_diff = int(epoch_stim == 'Mechanical')
                x = (constants.RASTER_PHASE_SPACING[phase] +
                     (epoch_number * 2 - x_diff) / frequency)

            if test_stim == 'Mechanical':
                plt.scatter(
                    x,
                    y_mech,
                    len(epoch_df)*constants.RASTER_POINT_SCALE,
                    constants.RASTER_COLOURS[epoch_stim]
                )
            elif test_stim == 'Electrical':
                plt.scatter(
                    x,
                    y_elec,
                    len(epoch_df)*constants.RASTER_POINT_SCALE,
                    constants.RASTER_COLOURS[epoch_stim]
                )
        if test_stim == 'Mechanical':
            y_mech += 1
        elif test_stim == 'Electrical':
            y_elec += 1

    if save_figure:
        base.save_plot('raster',constants.VERSION)
