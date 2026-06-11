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

from collections import abc
from cycler import cycler
import json
import math
import matplotlib as mpl
import matplotlib.axes as axes
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import ticker
from os import walk
import numpy as np
import pandas as pd
import re
from scipy import signal
from scipy.optimize import curve_fit
import statistics
from statistics import StatisticsError
import typing

from . import classes
from . import constants
from . import spike_extraction
from . import utils

cmap = mpl.colormaps[constants.plots.COLOURMAP]

# Regex patterns for parsing marker comments:
# Compiled once, rather than within a function, to reduce computations.
saved_pattern = re.compile(constants.regex.SAVED_FILENAME_REGEX)


####################### *USER-EXPOSED FUNCTIONS* #######################

def longduration_probability(
        spikes_df: pd.DataFrame,
        conduction_velocities: abc.Mapping[str, float],
        cv_norm: colors.Normalize,
        trace_units: abc.Collection[str],
        slow_trace: tuple[str, int],
        slow_spikes: str,
        slow_yspacing: float,
        slow_yscale: float,
        fast_trace: tuple[str, int],
        fast_spikes: str,
        fast_yspacing: float,
        fast_yscale: float,
        trace_times: tuple[float, float],
        constant_units: abc.Collection[str] = [],
        slowest_on_top: bool = False,
        rolling_bins: bool = False,
        save_plot: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df_categories = {'sex': constants.experiments.SEXES}
    fitted_df = typed_dataframe(
        pd.DataFrame(columns=constants.dataframes.FIG3_FITTEDDF_COLUMNS),
        constants.dataframes.FIG3_FITTEDDF_TYPES,
        df_categories
    )
    binned_df = typed_dataframe(
        pd.DataFrame(columns=constants.dataframes.FIG3_BINNEDDF_COLUMNS),
        constants.dataframes.FIG3_BINNEDDF_TYPES,
        df_categories
    )

    single_exponential_fig = plt.figure(0, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    single_exponential_gs = single_exponential_fig.add_gridspec(9, 7)
    double_exponential_fig = plt.figure(1, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    double_exponential_gs = double_exponential_fig.add_gridspec(9, 7)
    axis = (
        single_exponential_fig.add_subplot(single_exponential_gs[:2, :]),
        double_exponential_fig.add_subplot(double_exponential_gs[:2, :]),
        single_exponential_fig.add_subplot(single_exponential_gs[2:4, :]),
        double_exponential_fig.add_subplot(double_exponential_gs[2:4, :]),
        single_exponential_fig.add_subplot(single_exponential_gs[4:7, :]),
        double_exponential_fig.add_subplot(double_exponential_gs[4:7, :]),
        single_exponential_fig.add_subplot(single_exponential_gs[7:, :3]),
        single_exponential_fig.add_subplot(single_exponential_gs[7:, 4:]),
        double_exponential_fig.add_subplot(double_exponential_gs[7:, :2]),
        double_exponential_fig.add_subplot(double_exponential_gs[7:, 2:4]),
        double_exponential_fig.add_subplot(double_exponential_gs[7:, 4:6])
    )
    for ax in axis[:2]:
        trace_grid(
            ax,
            {
                slow_trace[0]: {
                    slow_trace[1]: (slow_spikes,(
                        (
                            'long_elec',
                            0,
                            0,
                            trace_times[0],
                            0,
                            f'{trace_times[0]} s'
                        ),
                        (
                            'long_elec',
                            1,
                            0,
                            trace_times[1],
                            0,
                            f'{trace_times[1]} s'
                        )
                    ))
                }
            },
            0.2,
            None,
            slow_yspacing,
            slow_yscale
        )
        # ax.set_title(f"{trace_times[0]} s | {trace_times[1]} s")
        ax.scatter(
            (0.21, 0.215, 0.22),
            (0, 0, 0),
            5,
            constants.plots.PALETTE['black']
        )
    for ax in axis[2:4]:
        trace_grid(
            ax,
            {
                fast_trace[0]: {
                    fast_trace[1]: (fast_spikes, (
                        (
                            'long_elec',
                            0,
                            0,
                            trace_times[0],
                            0,
                            f'{trace_times[0]} s'
                        ),
                        (
                            'long_elec',
                            1,
                            0,
                            trace_times[1],
                            0,
                            f'{trace_times[1]} s'
                        )
                    ))
                }
            },
            0.2,
            0.02,
            fast_yspacing,
            fast_yscale
        )
        ax.scatter(
            (0.21, 0.215, 0.22),
            (0, 0, 0),
            5,
            constants.plots.PALETTE['black']
        )
    long_df = spikes_df.loc[
        (spikes_df['Test'] == 'long-duration') &
        (spikes_df['Test Stimulus'] == 'electrical')
    ]
    for unit_id, conduction_velocity in sort_by_conduction_velocity(
        conduction_velocities,
        long_df['Unit ID'].unique(),
        slowest_on_top
    ):
        unit_df = long_df.loc[long_df['Unit ID'] == unit_id]
        cond_df = unit_df.loc[unit_df['Phase'] == 'conditioning']
        binned_data = calculate_binned_probability(
            cond_df,
            (0, constants.experiments.LONG_CONDITIONING_DURATION_SECONDS),
            constants.plots.LONG_BIN_WIDTH_S,
            rolling_bins
        )

        for ax in axis[4:6]:
            ax.plot(
                binned_data.x,
                binned_data.y,
                color=(
                    cmap(cv_norm(conduction_velocity)),
                    constants.plots.QUANTIFICATION_RAW_OPACITY
                ),
                linewidth=constants.plots.QUANTIFICATION_LINEWIDTH/2,
                clip_on=False,
                zorder=100
            )
            # for trace_time in trace_times:
            #     ax.plot(
            #         (trace_time, trace_time),
            #         (0, 1.1),
            #         color=constants.plots.PALETTE['light_grey'],
            #         linestyle='dashed',
            #         linewidth=constants.plots.QUANTIFICATION_LINEWIDTH,
            #         clip_on=False,
            #         zorder=102
            #     )
        if unit_id in constant_units:
            popt_single = ()
            popt_double = ()
            a_single = classes.FitParameter(float('nan'), float('nan'))
            k_single = classes.FitParameter(
                float('nan'),
                float('nan')
            )
            c_single = classes.FitParameter(
                statistics.mean(binned_data.y),
                np.std(binned_data.y) # type: ignore
            )
            fit_single = lambda t: np.array([c_single.value for _ in t])
            # a_double = classes.FitParameter(float('nan'), float('nan'))
            k1_double = classes.FitParameter(
                float('nan'),
                float('nan')
            )
            # b_double = classes.FitParameter(float('nan'), float('nan'))
            k2_double = classes.FitParameter(
                float('nan'),
                float('nan')
            )
            c_double = classes.FitParameter(
                statistics.mean(binned_data.y),
                np.std(binned_data.y) # type: ignore
            )
            h_double = classes.FitParameter(float('nan'), float('nan'))
            fit_double = lambda t: np.array([c_double.value for _ in t])
        else:
            popt_single, pcov_single = curve_fit(
                exponential_decay,
                binned_data.x,
                binned_data.y
            )
            perr_single = np.sqrt(np.diag(pcov_single))
            a_single = classes.FitParameter(popt_single[0],perr_single[0])
            k_single = classes.FitParameter(popt_single[1],perr_single[1])
            c_single = classes.FitParameter(
                1 - popt_single[0],
                float('nan')
            )
            fit_single = exponential_decay
            popt_double_unshifted, _pcov = curve_fit(
                double_exponential_decay,
                binned_data.x,
                binned_data.y
            )
            popt_double, pcov_double = curve_fit(
                double_exponential_decay_shifted,
                binned_data.x,
                binned_data.y,
                [*popt_double_unshifted, 0]
            )
            perr_double = np.sqrt(np.diag(pcov_double))
            if popt_double[1] >= popt_double[3]:
                # a_double = classes.FitParameter(popt_double[0],perr_double[0])
                k1_double = classes.FitParameter(popt_double[1],perr_double[1])
                # b_double = classes.FitParameter(popt_double[2],perr_double[2])
                k2_double = classes.FitParameter(popt_double[3],perr_double[3])
                h_double = classes.FitParameter(popt_double[4],perr_double[4])
            else:
                # a_double = classes.FitParameter(popt_double[2],perr_double[2])
                k1_double = classes.FitParameter(popt_double[3],perr_double[3])
                # b_double = classes.FitParameter(popt_double[0],perr_double[0])
                k2_double = classes.FitParameter(popt_double[1],perr_double[1])
                h_double = classes.FitParameter(popt_double[4],perr_double[4])
            c_double = classes.FitParameter(
                1 - popt_double[0] - popt_double[2],
                float('nan')
            )
            fit_double = double_exponential_decay_shifted
        error_single = rmse(
            fit_single(binned_data.x, *popt_single), # type: ignore
            binned_data.y
        )
        inspect_fit(
            fit_single, # type: ignore
            popt_single,
            binned_data.x,
            ax=axis[4],
            # error=error_single,
            color=cmap(cv_norm(conduction_velocity)),
            linewidth=constants.plots.QUANTIFICATION_LINEWIDTH
        )
        # error_double = rmse(
        #     fit_double(binned_data.x, *popt_double), # type: ignore
        #     binned_data.y
        # )
        inspect_fit(
            fit_double, # type: ignore
            popt_double,
            binned_data.x,
            max(0, h_double.value),
            ax=axis[5],
            # error=error_double,
            color=cmap(cv_norm(conduction_velocity)),
            linewidth=constants.plots.QUANTIFICATION_LINEWIDTH
        )
        for i, y in enumerate((
            c_single,
            k_single,
            c_double,
            k1_double,
            k2_double
        )):
            if math.isnan(y.value):
                y_value = 0
                size = constants.plots.QUANTIFICATION_POINT_SIZE*2
                colour = constants.plots.PALETTE['vermillion']
                marker = 'x'
            else:
                y_value = y.value
                size = constants.plots.QUANTIFICATION_POINT_SIZE
                colour = constants.plots.PALETTE['black']
                marker = 'o'
            axis[i+6].scatter(
                conduction_velocity,
                y_value,
                s=size,
                c=colour,
                marker=marker,
                clip_on=False,
                zorder=101
            )
            if unit_id in trace_units:
                axis[i+6].plot(
                    conduction_velocity,
                    y_value,
                    marker='s',
                    markersize=size*0.8,
                    mfc='none',
                    mec=constants.plots.PALETTE['light_grey'],
                    clip_on=False,
                    zorder=100
                )
        common = {
            'animal_id': utils.unique(unit_df['Animal ID']),
            'unit_id': unit_id,
            'test_id': utils.unique(unit_df['Test ID']),
            'sex': utils.unique(unit_df['Sex']),
            'unit_type': utils.unique(unit_df['Unit Type']),
            'conduction_velocity': conduction_velocity
        }
        fitted_df.loc[len(fitted_df)] = {
            **common,
            'a': a_single.value,
            'a_std': a_single.error,
            'k': k_single.value,
            'k_std': k_single.error,
            'c': c_single.value,
            'c_std': c_single.error,
            'rmse': error_single
        }
        # a_actual = a.value*np.exp(k1.value*h.value)
        # b_actual = b.value*np.exp(k2.value*h.value)
        # fitted_df.loc[len(fitted_df)] = {
        #     **common,
        #     'a': a.value,
        #     'a_std': a.error,
        #     'k1': k1.value,
        #     'k1_std': k1.error,
        #     'b': b.value,
        #     'b_std': b.error,
        #     'k2': k2.value,
        #     'k2_std': k2.error,
        #     'c': c.value,
        #     'c_std': c.error,
        #     'h': h.value,
        #     'h_std': h.error,
        #     'a_actual': a_actual,
        #     'b_actual': b_actual,
        #     'ab_ratio': a_actual/b_actual,
        #     'rmse': error
        # }
        for trial_id in unit_df['Trial ID'].unique():
            trial_df = unit_df.loc[
                (unit_df['Trial ID'] == trial_id) &
                (unit_df['Phase'] == 'conditioning')
            ]
            trial_bins = calculate_binned_probability(
                trial_df,
                (0, constants.experiments.LONG_CONDITIONING_DURATION_SECONDS),
                constants.plots.LONG_BIN_WIDTH_S,
                # constants.LONG_CONDITIONING_FREQUENCY,
                rolling_bins
            )
            for n, x_n in enumerate(trial_bins.x):
                y_n = trial_bins.y[n]
                log_y = np.log(y_n) if y_n!=0 else np.log(utils.small_float)
                binned_df.loc[len(binned_df)] = {
                    **common,
                    'trial_id': trial_id,
                    'time': x_n,
                    'exp_time': np.exp(x_n),
                    'response_probability': y_n,
                    'log_probability': log_y
                }
    # Finishing the plot
    for ax in axis[4:6]:
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
        colourbar(ax, cv_norm)
        ax.set_xlim(
            0,
            constants.experiments.LONG_CONDITIONING_DURATION_SECONDS
        )
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(
            constants.plots.PROBABILITY_MINORTICK_SPACING
        ))
    for ax in axis[6:]:
        ax.set_xlim(
            0,
            round_sigfigs(max(conduction_velocities.values()), 'ceil')
        )
        ax.xaxis.set_major_locator(ticker.MultipleLocator(
            constants.plots.FIG2_CONDUCTIONVELOCITY_MAJORTICK_SPACING
        ))
    axis[6].yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.PROBABILITY_MINORTICK_SPACING_SMALL
    ))
    axis[7].yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.FIG2_SINGLE_INVTAU_MINORTICK_SPACING
    ))
    axis[8].yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.PROBABILITY_MINORTICK_SPACING_SMALL
    ))
    axis[9].yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.FIG2_FAST_INVTAU_MINORTICK_SPACING
    ))
    axis[10].yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.FIG2_SLOW_INVTAU_MINORTICK_SPACING
    ))
    apply_setp(axis, 'long-duration')
    plt.figure(0)
    plt.suptitle(constants.plots.FIG2_TITLE)
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 2"
        )
    plt.figure(1)
    plt.suptitle(constants.plots.FIG2_TITLE)
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 2 (double exponential)"
        )
    plt.close(1)
    return fitted_df, binned_df


def sweeps_probability(
        spikes_df: pd.DataFrame,
        conduction_velocities: abc.Mapping[str, float],
        cv_norm: colors.Normalize,
        trace_unit: str,
        frequency_trace: tuple[str, int],
        frequency_spikes: str,
        frequency_trials: abc.Sequence[tuple[float, float]],
        frequency_yspacing: float,
        frequency_yscale: float,
        amplitude_trace: tuple[str, int],
        amplitude_spikes: str,
        amplitude_trials: abc.Sequence[tuple[float, float]],
        amplitude_yspacing: float,
        amplitude_yscale: float,
        first_n: int | None = None,
        extend_jn_boundaries: bool = False,
        slowest_on_top: bool = False,
        save_plot: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if first_n is None:
        sweeps_df = spikes_df.loc[
            ((spikes_df['Test'] == 'frequency') |
             (spikes_df['Test'] == 'amplitude')) &
            (spikes_df['Phase'] == 'conditioning')
        ]
        first_n = math.floor(
            min(sweeps_df['Test Frequency'])*
            constants.experiments.DEFAULT_CONDITIONING_DURATION_SECONDS
        )

    frequency_fig = plt.figure(0, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    frequency_gs = frequency_fig.add_gridspec(9, 7)
    amplitude_fig = plt.figure(1, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    amplitude_gs = amplitude_fig.add_gridspec(9, 7)
    axis = (
        frequency_fig.add_subplot(frequency_gs[:4, :]),
        amplitude_fig.add_subplot(amplitude_gs[:4, :]),
        frequency_fig.add_subplot(frequency_gs[4:7, :]),
        frequency_fig.add_subplot(frequency_gs[7:, 4:]),
        amplitude_fig.add_subplot(amplitude_gs[4:7, :]),
        frequency_fig.add_subplot(frequency_gs[7:, :3])
    )
    trace_grid(
        axis[0],
        {
            frequency_trace[0]: {
                frequency_trace[1]: (frequency_spikes, [
                    (
                        f'frequency_mech_0_elec_{frequency}',
                        0,
                        i,
                        trace_time,
                        0,
                        f'{frequency} {constants.core.TEST_UNITS['frequency']}'
                    ) for i,(frequency,trace_time) in enumerate(frequency_trials)
                ])
            }
        },
        0.4,
        0.02,
        frequency_yspacing,
        frequency_yscale
    )
    trace_grid(
        axis[1],
        {
            amplitude_trace[0]: {
                amplitude_trace[1]: (amplitude_spikes, [
                    (
                        f'amplitude_mech_0_elec_{amplitude}',
                        0,
                        i,
                        trace_time,
                        0,
                        f'{amplitude} {constants.core.TEST_UNITS['amplitude']}'
                    ) for i,(amplitude,trace_time) in enumerate(amplitude_trials)
                ])
            }
        },
        0.4,
        0.02,
        amplitude_yspacing,
        amplitude_yscale
    )
    outputs = {}
    for test, ax in (
        ('frequency', axis[2:4]),
        ('amplitude', [axis[4]]),
    ):
        test_df = spikes_df.loc[
            (spikes_df['Test'] == test) &
            (spikes_df['Test Stimulus'] == 'electrical')
        ]
        conduction_velocities_sorted = sort_by_conduction_velocity(
            conduction_velocities,
            test_df['Unit ID'].unique(),
            slowest_on_top
        )

        output_df = typed_dataframe(
            pd.DataFrame(columns=[
                *constants.dataframes.FIG4_OUTPUTDF_COLUMNS,
                test
            ]),
            {**constants.dataframes.FIG4_OUTPUTDF_TYPES, test: 'float64'},
            {'sex': constants.experiments.SEXES}
        )

        values = sorted(test_df[f'Test {test.capitalize()}'].unique())
        if test == 'frequency':
            values = [int(x) for x in values]

        datapoint_dicts = [{i:[] for i,_ in enumerate(values)} for _ in ax]
        for i, value in enumerate(values):
            for unit_id, conduction_velocity in conduction_velocities_sorted:
                unit_df = test_df.loc[
                    (test_df['Unit ID'] == unit_id) &
                    (test_df[f'Test {test.capitalize()}'] == value)
                ]
                if unit_df.empty:
                    continue

                animal_id = utils.unique(unit_df['Animal ID'])
                position = utils.unique(unit_df['Position'])
                amplitude = utils.unique(unit_df['Test Amplitude'])
                stim_current = (constants.experiments.ANIMAL_DATA
                                [animal_id.upper()][position]['elec_amp'] *
                                amplitude)
                
                x = np.random.normal(i+1, constants.plots.JITTER_AMOUNT_NARROW)

                unit_data = calculate_unbinned_probability(
                    unit_df,
                    first_n,
                )
                for trial_dict in unit_data:
                    absolute_charge = trial_dict.pop('epochs') * stim_current
                    probabilities = trial_dict.pop('response_probability')
                    output_df.loc[len(output_df)] = {
                        **trial_dict,
                        'conduction_velocity': conduction_velocity,
                        'absolute_charge': absolute_charge,
                        'response_probability_all': probabilities[0],
                        'response_probability_initial': probabilities[1],
                        test: value
                    }
                    trial_dict['response_probability'] = probabilities
                for j, data_dict in enumerate(datapoint_dicts):
                    probability = statistics.mean([
                        x['response_probability'][j] for x in unit_data
                    ])
                    data_dict[i].append(
                        classes.BoxPoint(
                            i,
                            value,
                            probability
                        )
                    )
                    ax[j].scatter(
                        x,
                        probability,
                        s=constants.plots.QUANTIFICATION_POINT_SIZE*2/(j+1),
                        c=conduction_velocity,
                        cmap=cmap, # type: ignore
                        norm=cv_norm,
                        alpha=constants.plots.QUANTIFICATION_RAW_OPACITY,
                        clip_on=False,
                        zorder=101
                    )
                    trace_trials = (frequency_trials if test == 'frequency' else
                                    amplitude_trials)
                    trace_values = [x[0] for x in trace_trials]
                    if unit_id == trace_unit and value in trace_values:
                        ax[j].plot(
                            x,
                            probability,
                            marker='s',
                            markersize=constants.plots.QUANTIFICATION_POINT_SIZE,
                            mfc='none',
                            mec=constants.plots.PALETTE['light_grey'],
                            clip_on=False,
                            zorder=100
                        )
        for i, data_dict in enumerate(datapoint_dicts):
            datapoints = [[point.y for point in data_dict[j]]
                          for j,_ in enumerate(values)]

            boxplot_with_mean_colour(
                ax[i],
                x=datapoints,
                tick_labels=[str(x) for x in values],
                widths=constants.plots.QUANTIFICATION_BOXPLOT_WIDTH_NARROW
            )

        outputs[test] = output_df
    johnson_neyman(
        axis[5],
        (spikes_df['Test Frequency'].min(), spikes_df['Test Frequency'].max()),
        extend_jn_boundaries
    )
    for ax in (axis[2], axis[4]):
        colourbar(ax, cv_norm)
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(
            constants.plots.PROBABILITY_MINORTICK_SPACING
        ))
    axis[3].yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.PROBABILITY_MINORTICK_SPACING_SMALL
    ))
    apply_setp(axis, 'sweeps')
    plt.figure(0)
    plt.suptitle(constants.plots.FREQUENCY_SWEEP_TITLE)
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 3"
        )
    plt.figure(1)
    plt.suptitle(constants.plots.AMPLITUDE_SWEEP_TITLE)
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 4"
        )
    return outputs['frequency'], outputs['amplitude']


def nineone_probability(
        spikes_df: pd.DataFrame,
        conduction_velocities: abc.Mapping[str, float],
        cv_norm: colors.Normalize,
        trace_unit: str,
        nineone_trace: tuple[str, int],
        nineone_spikes: str,
        nineone_time: float,
        mech_trace: tuple[str, int],
        mech_spikes: str,
        mech_time: float,
        yspacing: float,
        yscale: float,
        slowest_on_top: bool = False,
        save_plot: bool = True
) -> pd.DataFrame:
    actual_fig = plt.figure(0, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    actual_gs = actual_fig.add_gridspec(9, 7)
    differences_fig = plt.figure(1, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    differences_gs = differences_fig.add_gridspec(9, 7)
    axis = (
        actual_fig.add_subplot(actual_gs[:4, :]),
        differences_fig.add_subplot(differences_gs[:4, :]),
        actual_fig.add_subplot(actual_gs[4:7, :]),
        differences_fig.add_subplot(differences_gs[4:7, :])
    )
    for ax in axis[:2]:
        trace_grid(
            ax,
            {
                nineone_trace[0]: {
                    nineone_trace[1]: (
                        nineone_spikes,
                        [(
                            'mech_0_elec_1',
                            0,
                            0,
                            nineone_time,
                            0,
                            'Interleaved'
                        )]
                    )
                },
                mech_trace[0]: {
                    mech_trace[1]: (
                        mech_spikes,
                        [(
                            'frequency_mech_10_elec_0',
                            0,
                            1,
                            mech_time,
                            -1/constants.experiments.DEFAULT_CONDITIONING_FREQUENCY,
                            'Mechanical'
                        )]
                    )
                }
            },
            0.4,
            0.02,
            yspacing,
            yscale
        )
    output_df = typed_dataframe(
        pd.DataFrame(columns=constants.dataframes.FIG5_OUTPUTDF_COLUMNS),
        constants.dataframes.FIG5_OUTPUTDF_TYPES,
        {
            'sex': constants.experiments.SEXES,
            'test': constants.experiments.TEST_CODES.values(),
            'stimulus': constants.experiments.STIMULATION_TYPES
        }
    )
    datapoints = []
    for i, stim_type in enumerate(constants.experiments.STIMULATION_TYPES):
        stim_data = []
        interleaved_df = spikes_df.loc[
            (spikes_df['Test'] == 'nine-one') &
            (spikes_df['Test Stimulus'] == 'electrical') &
            (spikes_df['Epoch Stimulus'] == stim_type)
        ]
        x_actual = i * constants.plots.FIG5_ACTUAL_INTERSPACING
        frequency = utils.unique(interleaved_df['Minor Test Frequency'])
        if stim_type == 'electrical':
            frequency += utils.unique(interleaved_df['Test Frequency'])
        for unit_id, conduction_velocity in sort_by_conduction_velocity(
            conduction_velocities,
            interleaved_df['Unit ID'].unique(),
            slowest_on_top
        ):
            unit_interleaved_data = calculate_unbinned_probability(
                interleaved_df.loc[interleaved_df['Unit ID'] == unit_id]
            )
            unit_regular_data = calculate_unbinned_probability(
                spikes_df.loc[
                    (spikes_df['Unit ID'] == unit_id) &
                    (spikes_df['Test Stimulus'] == stim_type) &
                    (spikes_df['Epoch Stimulus'] == stim_type) &
                    (spikes_df['Test Frequency'] == frequency) &
                    (spikes_df['Test'] == 'frequency')
                ]
            )

            for trial_dict in unit_interleaved_data:
                del trial_dict['epochs']
                output_df.loc[len(output_df)] = {
                    **trial_dict,
                    'test': 'nine-one',
                    'stimulus': stim_type,
                    'conduction_velocity': conduction_velocity
                }
            for trial_dict in unit_regular_data:
                del trial_dict['epochs']
                output_df.loc[len(output_df)] = {
                    **trial_dict,
                    'test': 'frequency',
                    'stimulus': stim_type,
                    'conduction_velocity': conduction_velocity
                }

            jitter = np.random.normal(0, constants.plots.JITTER_AMOUNT_WIDE)
            unit_regular_mean = statistics.mean([
                x['response_probability'] for x in unit_regular_data
            ])
            unit_interleaved_mean = statistics.mean([
                x['response_probability'] for x in unit_interleaved_data
            ])
            for ax, x, y, in (
                (
                    axis[2],
                    x_actual + jitter - constants.plots.FIG5_ACTUAL_INTRASPACING,
                    unit_regular_mean
                ),
                (
                    axis[2],
                    x_actual + jitter + constants.plots.FIG5_ACTUAL_INTRASPACING,
                    unit_interleaved_mean
                ),
                (
                    axis[3],
                    i + jitter + 1,
                    unit_interleaved_mean - unit_regular_mean
                )
            ):
                ax.scatter(
                    x,
                    y,
                    s=constants.plots.FIG5_POINT_SIZE,
                    c=conduction_velocity,
                    cmap=cmap,
                    norm=cv_norm,
                    alpha=constants.plots.QUANTIFICATION_RAW_OPACITY,
                    clip_on=False,
                    zorder=101
                )
                if unit_id == trace_unit:
                    ax.plot(
                        x,
                        y,
                        marker='s',
                        markersize=constants.plots.FIG5_POINT_SIZE/2,
                        mfc='none',
                        mec=constants.plots.PALETTE['light_grey'],
                        clip_on=False,
                        zorder=100
                    )
            # axis[2].scatter(
            #     x_actual + jitter - constants.plots.FIG5_ACTUAL_INTRASPACING,
            #     unit_regular_mean,
            #     s=constants.plots.FIG5_POINT_SIZE,
            #     c=conduction_velocity,
            #     cmap=cmap,
            #     norm=cv_norm,
            #     alpha=constants.plots.QUANTIFICATION_RAW_OPACITY,
            #     clip_on=False,
            #     zorder=100
            # )
            # axis[2].scatter(
            #     x_actual + jitter + constants.plots.FIG5_ACTUAL_INTRASPACING,
            #     unit_interleaved_mean,
            #     s=constants.plots.FIG5_POINT_SIZE,
            #     c=conduction_velocity,
            #     cmap=cmap,
            #     norm=cv_norm,
            #     alpha=constants.plots.QUANTIFICATION_RAW_OPACITY,
            #     clip_on=False,
            #     zorder=100
            # )
            # axis[3].scatter(
            #     i + jitter + 1,
            #     unit_interleaved_mean - unit_regular_mean,
            #     s=constants.plots.FIG5_POINT_SIZE,
            #     c=conduction_velocity,
            #     cmap=cmap,
            #     norm=cv_norm,
            #     alpha=constants.plots.QUANTIFICATION_RAW_OPACITY,
            #     clip_on=False,
            #     zorder=100
            # )
            stim_data.append((unit_regular_mean, unit_interleaved_mean))
        datapoints_regular = [x[0] for x in stim_data]
        datapoints_interleaved = [x[1] for x in stim_data]
        boxplot_with_mean_colour(
            axis[2],
            x=(datapoints_regular, datapoints_interleaved),
            positions=(
                x_actual - constants.plots.FIG5_ACTUAL_INTRASPACING,
                x_actual + constants.plots.FIG5_ACTUAL_INTRASPACING
            ),
            widths = constants.plots.QUANTIFICATION_BOXPLOT_WIDTH_WIDE
        )
        datapoints.append(stim_data)
    plt.figure(0)
    plt.suptitle(constants.plots.FIG5_TITLE)
    colourbar(axis[2], cv_norm)
    axis[2].set_xticklabels(**constants.plots.FIG5_ACTUAL_PRIMARY_LABELS)
    sec = axis[2].secondary_xaxis(location=0)
    sec.set_xticks(np.arange(len(constants.experiments.STIMULATION_TYPES)) *
                   constants.plots.FIG5_ACTUAL_INTERSPACING)
    sec.set_xticklabels([f'\n\n{x.capitalize()}' for x in
                         constants.experiments.STIMULATION_TYPES])
    sec.xaxis.set_ticks_position('none')
    axis[2].yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.PROBABILITY_MINORTICK_SPACING
    ))
    apply_setp(axis, 'nine-one')
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 5 (actual)"
        )
    boxplot_with_mean_colour(
        axis[3],
        x=[[x[1]-x[0] for x in stim_data] for stim_data in datapoints],
        widths = constants.plots.QUANTIFICATION_BOXPLOT_WIDTH_NARROW
    )
    plt.figure(1)
    plt.suptitle(constants.plots.FIG5_TITLE)
    colourbar(axis[3], cv_norm)
    axis[3].set_xticklabels([x.capitalize() for x in
                        constants.experiments.STIMULATION_TYPES])
    axis[3].yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.FIG5_DIFFERENCES_MINORTICK_SPACING
    ))
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 5 (differences)"
        )
    plt.close(1)
    return output_df


def peak_property_changes(
        spikes_df: pd.DataFrame,
        initial_properties: abc.Mapping[
            str,
            classes.InitialProperties
        ],
        amplitude_sweep_outliers: abc.Collection[str],
        latency_sweep_outliers: abc.Collection[str],
        cv_norm: colors.Normalize,
        first_n: int | None = None,
        slowest_on_top: bool = False,
        save_plot: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    fitted_df = typed_dataframe(
        pd.DataFrame(columns=constants.dataframes.FIG6_FITTEDDF_COLUMNS),
        constants.dataframes.FIG6_FITTEDDF_TYPES,
        {
            'sex': constants.experiments.SEXES,
            'test': constants.experiments.TEST_CODES.values()
        }
    )
    sweeps_df = typed_dataframe(
        pd.DataFrame(columns=constants.dataframes.FIG6_SWEEPSDF_COLUMNS),
        constants.dataframes.FIG6_SWEEPSDF_TYPES,
        {
            'sex': constants.experiments.SEXES,
            'test': constants.experiments.TEST_CODES.values()
        }
    )
    latency_fig = plt.figure(0, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    latency_gs = latency_fig.add_gridspec(5, 7)
    amplitude_fig = plt.figure(1, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    amplitude_gs = latency_fig.add_gridspec(5, 7)
    axis = (
        latency_fig.add_subplot(latency_gs[:2, :]),
        # latency_fig.add_subplot(latency_gs[2:, :2]),
        latency_fig.add_subplot(latency_gs[2:, :3]),
        latency_fig.add_subplot(latency_gs[2:, 4:]),
        amplitude_fig.add_subplot(amplitude_gs[:2, :]),
        # amplitude_fig.add_subplot(amplitude_gs[2:, :2]),
        amplitude_fig.add_subplot(amplitude_gs[2:, :3]),
        amplitude_fig.add_subplot(amplitude_gs[2:, 4:])
    )
    if first_n:
        spikes_df = spikes_df.loc[spikes_df['Epoch Number'] < first_n]
    properties_sorted = sort_by_conduction_velocity(
        initial_properties,
        spikes_df['Unit ID'].unique(),
        slowest_on_top,
        lambda x: x[1].conduction_velocity
    )
    
    for unit_id, properties in properties_sorted:
        long_df = spikes_df.loc[
            (spikes_df['Unit ID'] == unit_id) &
            (spikes_df['Test'] == 'long-duration') &
            (spikes_df['Test Stimulus'] == 'electrical') &
            (spikes_df['Phase'] == 'conditioning')
        ]
        if long_df.empty:
            continue

        fitted_dict = fit_peak_properties(
            *peak_properties(long_df),
            properties,
            (axis[0], axis[3]),
            cv_norm
        )
        # axis[1].scatter(
        #     properties.conduction_velocity,
        #     fitted_dict['latency_normal_steady_cv']*100,
        #     s=constants.plots.QUANTIFICATION_POINT_SIZE,
        #     c=constants.plots.PALETTE['black'],
        #     clip_on=False,
        #     zorder=100
        # )
        # axis[5].scatter(
        #     properties.conduction_velocity,
        #     fitted_dict['amplitude_normal_slope']*100,
        #     s=constants.plots.QUANTIFICATION_POINT_SIZE,
        #     c=constants.plots.PALETTE['black'],
        #     clip_on=False,
        #     zorder=100
        # )
        fitted_df.loc[len(fitted_df)] = {
            **fitted_dict,
            'animal_id': utils.unique(long_df['Animal ID']),
            'unit_id': unit_id,
            'test_id': utils.unique(long_df['Test ID']),
            'sex': utils.unique(long_df['Sex']),
            'unit_type': utils.unique(long_df['Unit Type']),
            'test': 'long-duration',
            'test_frequency': utils.unique(long_df['Test Frequency']),
            'test_amplitude': utils.unique(long_df['Test Amplitude']),
        }

    for test, ax in (
        ('frequency', (axis[1], axis[4])),
        ('amplitude', (axis[2], axis[5]))
    ):
        test_df = spikes_df.loc[
            (spikes_df['Test'] == test) &
            (spikes_df['Test Stimulus'] == 'electrical') &
            (spikes_df['Phase'] == 'conditioning')
        ]

        values = sorted(test_df[f'Test {test.capitalize()}'].unique())
        if test == 'frequency':
            values = [int(x) for x in values]

        datapoints = ([[] for _ in values], [[] for _ in values])
        for i, value in enumerate(values):
            for unit_id, properties in properties_sorted:
                unit_df = test_df.loc[
                    (test_df['Unit ID'] == unit_id) &
                    (test_df[f'Test {test.capitalize()}'] == value)
                ]
                if unit_df.empty:
                    continue

                animal_id = utils.unique(unit_df['Animal ID'])
                position = utils.unique(unit_df['Position'])
                test_id = utils.unique(unit_df['Test ID'])
                sex = utils.unique(unit_df['Sex'])
                unit_type = utils.unique(unit_df['Unit Type'])
                test_frequency = utils.unique(unit_df['Test Frequency'])
                test_amplitude = utils.unique(unit_df['Test Amplitude'])
                stim_current = (constants.experiments.ANIMAL_DATA
                                [animal_id.upper()][position]['elec_amp'] *
                                test_amplitude)

                unit_mean_amplitudes = []
                unit_mean_latencies = []
                unit_mean_cvs = []
                for trial_id in unit_df['Trial ID'].unique():
                    trial_df = unit_df.loc[unit_df['Trial ID'] == trial_id]

                    epoch_count = utils.unique(trial_df['Total Epochs'])
                    absolute_charge = epoch_count * stim_current

                    amplitude_actual_mean = trial_df['Size (μV)'].mean()
                    amplitude_normal_mean = normalise(
                        amplitude_actual_mean,
                        properties.amplitude
                    )
                    latency_actual_mean = trial_df['Latency (ms)'].mean()
                    latency_normal_mean = normalise(
                        latency_actual_mean,
                        properties.latency
                    )
                    cv_actual_mean = properties.distance/latency_actual_mean
                    cv_normal_mean = normalise(
                        cv_actual_mean,
                        properties.conduction_velocity
                    )

                    unit_mean_amplitudes.append(
                        (amplitude_actual_mean, amplitude_normal_mean)
                    )
                    unit_mean_latencies.append(
                        (latency_actual_mean, latency_normal_mean)
                    )
                    unit_mean_cvs.append(
                        (cv_actual_mean, cv_normal_mean)
                    )

                    sweeps_df.loc[len(sweeps_df)] = {
                        'animal_id': animal_id,
                        'unit_id': unit_id,
                        'test_id': test_id,
                        'trial_id': trial_id,
                        'sex': sex,
                        'unit_type': unit_type,
                        'initial_conduction_velocity': properties.conduction_velocity,
                        'initial_amplitude': properties.amplitude,
                        'initial_latency': properties.latency,
                        'recording_distance': properties.distance,
                        'test': test,
                        'test_frequency': test_frequency,
                        'test_amplitude': test_amplitude,
                        'absolute_charge': absolute_charge,
                        'latency_actual_mean': latency_actual_mean,
                        'latency_normal_mean': latency_normal_mean,
                        'cv_actual_mean': cv_actual_mean,
                        'cv_normal_mean': cv_normal_mean,
                        'amplitude_actual_mean': amplitude_actual_mean,
                        'amplitude_normal_mean': amplitude_normal_mean
                    }

                x = np.random.normal(i+1, constants.plots.JITTER_AMOUNT_NARROW)
                cv = statistics.mean([x[1] for x in unit_mean_cvs])*100
                amplitude = statistics.mean([x[1] for x in unit_mean_amplitudes])*100

                if unit_id in latency_sweep_outliers:
                    ax[0].scatter(
                        x,
                        cv,
                        s=constants.plots.QUANTIFICATION_POINT_SIZE*2,
                        c=constants.plots.PALETTE['vermillion'],
                        marker='x',
                        clip_on=False,
                        zorder=101
                    )
                else:
                    ax[0].scatter(
                        x,
                        cv,
                        s=constants.plots.QUANTIFICATION_POINT_SIZE,
                        c=properties.conduction_velocity,
                        cmap=cmap,
                        norm=cv_norm,
                        alpha=constants.plots.QUANTIFICATION_RAW_OPACITY,
                        clip_on=False,
                        zorder=100
                    )
                    datapoints[0][i].append(cv)

                if unit_id in amplitude_sweep_outliers:
                    ax[1].scatter(
                        x,
                        amplitude,
                        s=constants.plots.QUANTIFICATION_POINT_SIZE*2,
                        c=constants.plots.PALETTE['vermillion'],
                        marker='x',
                        clip_on=False,
                        zorder=101
                    )
                else:
                    ax[1].scatter(
                        x,
                        amplitude,
                        s=constants.plots.QUANTIFICATION_POINT_SIZE,
                        c=properties.conduction_velocity,
                        cmap=cmap,
                        norm=cv_norm,
                        alpha=constants.plots.QUANTIFICATION_RAW_OPACITY,
                        clip_on=False,
                        zorder=100
                    )
                    datapoints[1][i].append(amplitude)
        for i, data_list in enumerate(datapoints):
            boxplot_with_mean_colour(
                ax[i],
                x=data_list,
                tick_labels=[str(x) for x in values],
                widths=constants.plots.QUANTIFICATION_BOXPLOT_WIDTH_NARROW
            )
    for ax in (axis[0], axis[3]):
        # legend = ax.get_legend()
        # if legend is not None:
        #     legend.remove()
        colourbar(ax, cv_norm)
        ax.set_xlim(
            0,
            constants.experiments.LONG_CONDITIONING_DURATION_SECONDS
        )
    if first_n:
        apply_setp(axis, 'peak_properties_first30')
    else:
        apply_setp(axis, 'peak_properties')
    plt.figure(0)
    plt.suptitle("Spike latency and CV changes")
    # axis[0].set_ylim((1.6, 6.4))
    axis[0].yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    for ax in axis[1:3]:
        minor_ticks = 0.4 if first_n else 1
        # ax.set_ylim(-15, 15)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d%%'))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(minor_ticks))
    plt.tight_layout()
    filename_suffix = f' (First {first_n})' if first_n else ''
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 6" + filename_suffix
        )
    plt.figure(1)
    plt.suptitle("Spike amplitude changes")
    axis[3].yaxis.set_minor_locator(ticker.MultipleLocator(40))
    for ax in axis[4:6]:
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d%%'))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(2))
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 7" + filename_suffix
        )
    return fitted_df, sweeps_df


def model_output(
        y_spacing: float,
        x_scale: float,
        y_scale: float,
        save_plot: bool = True
) -> None:
    # import csvs from filenames in constants as dataframe
    trace_dfs_7um = model_dfs(
        '7um_',
        'Hz_1T\\Abeta0.txt',
        constants.core.MODEL_FREQUENCIES
    )
    gate_dfs_7um = model_dfs(
        '7um_',
        'Hz_1T\\Abeta0_s.txt',
        constants.plots.FIG8_FREQUENCIES.keys()
    )
    trace_dfs_14um = model_dfs(
        '14um_',
        'Hz_1T\\Aalpha0.txt',
        constants.core.MODEL_FREQUENCIES
    )
    gate_dfs_14um = model_dfs(
        '14um_',
        'Hz_1T\\Aalpha0_s.txt',
        constants.plots.FIG8_FREQUENCIES.keys()
    )

    fig_7um = plt.figure(0, figsize=(
        constants.plots.FIG_WIDTH*0.8,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    gs_7um = fig_7um.add_gridspec(9, 4)
    fig_14um = plt.figure(1, figsize=(
        constants.plots.FIG_WIDTH*0.8,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    gs_14um = fig_14um.add_gridspec(9, 4)
    axis = (
        fig_7um.add_subplot(gs_7um[0, :]),
        fig_7um.add_subplot(gs_7um[1, :]),
        fig_7um.add_subplot(gs_7um[2, :]),
        fig_14um.add_subplot(gs_14um[0, :]),
        fig_14um.add_subplot(gs_14um[1, :]),
        fig_14um.add_subplot(gs_14um[2, :]),
        fig_7um.add_subplot(gs_7um[3:5, :3]),
        fig_7um.add_subplot(gs_7um[5:7, :3]),
        fig_14um.add_subplot(gs_14um[3:5, :3]),
        fig_14um.add_subplot(gs_14um[5:7, :3]),
        fig_7um.add_subplot(gs_7um[7:, :3]),
        fig_14um.add_subplot(gs_14um[7:, :3])
    )

    colour_cycle = [
        constants.plots.PALETTE['sky'],
        '#db6d00',
        '#006dff',
        '#920000',
        '#490092',
        constants.plots.PALETTE['black']
    ]
    linestyle_cycle = [':','--','-','-.',(0,(3,2,1,2,1,2)),(0,(3,1,3,1,1,1))]

    for ax, trace_dfs in (
        (axis[:3], trace_dfs_7um),
        (axis[3:6], trace_dfs_14um)
    ):
        for i, (frequency, trace) in enumerate(trace_dfs.items()):
            # trace_x = np.array(i*0.000025 for i,_x in enumerate(trace['y']))
            trace_y = np.array(trace['y']) + 80
            duration_s = 0.5
            x_spacing = duration_s * 1.15
            tick_dt = 0.000025
            tick_dt_ms = tick_dt * constants.core.MILLISECONDS_PER_SECOND
            triggers_s = np.arange(0, duration_s, 1/frequency)
            triggers = [int(x/tick_dt) for x in triggers_s]
            # trace_y_blanked = spike_extraction.blank_artefacts(
            #     trace_y,
            #     triggers,
            #     1.2,
            #     tick_dt_ms
            # )
            # Define epoch bounds in samples relative to the trigger:
            epoch_end = int(10/tick_dt_ms) + 1
            spikes = {}
            for j, trigger in enumerate(triggers):
                # Create `DataEpoch` object describing current epoch and add to
                # `data_epochs` list:
                trace = np.array(trace_y[trigger:trigger+epoch_end])
                # Detect peaks greater than `threshold` within current epoch and
                # add to `spikes` list:
                peaks, _properties = signal.find_peaks(
                    trace,
                    20,
                    distance=constants.core.MINIMUM_SPIKE_DISTANCE_MS/tick_dt_ms
                )
                spikes[str(j)] = peaks*tick_dt_ms
            plot_trace(
                ax[i],
                np.arange(stop=duration_s+tick_dt, step=tick_dt), # type: ignore
                trace_y,
                (0, int(duration_s/tick_dt)),
                classes.TraceGridCoordinates(
                    '',
                    0,
                    0,
                    0,
                    0,
                    ' '.join((str(frequency), constants.core.TEST_UNITS['frequency']))
                ),
                x_spacing,
                y_spacing,
                triggers_s, # type: ignore
                [],
                [i for i,_ in enumerate(triggers)],
                [],
                spikes,
                {},
                frequency,
                -1,
                tick_dt,
                0.1,
                linewidth=constants.plots.TRACE_LINEWIDTH*4
            )
            x_scale_pos = y_spacing * -0.5
            x_scalelabel_pos = x_scale_pos - y_spacing*0.15
            y_scale_pos = x_spacing * -0.025
            y_scalelabel_pos = y_scale_pos - x_spacing*0.025
            if i == 2:
                ax[i].plot(
                    (0, x_scale),
                    (x_scale_pos, x_scale_pos),
                    c=constants.plots.PALETTE['black'],
                    clip_on=False,
                    zorder=101
                )
                ax[i].text(
                    0,
                    x_scalelabel_pos,
                    f"{int(x_scale*constants.core.MILLISECONDS_PER_SECOND)} ms",
                    va='top'
                )
            if i == 0:
                ax[i].plot(
                    (y_scale_pos, y_scale_pos),
                    (0, y_scale),
                    c=constants.plots.PALETTE['black'],
                    clip_on=False,
                    zorder=101
                )
                ax[i].text(
                    y_scalelabel_pos,
                    0,
                    f"{int(y_scale)} mV",
                    ha='left',
                    va='bottom',
                    rotation='vertical'
                )
            ax[i].set_xlim(0, x_spacing * 0.9)
            ax[i].set_ylim(x_scalelabel_pos, 100)
            for spine in ax[i].spines.values():
                spine.set_visible(False)
            ax[i].set_xticks([])
            ax[i].set_yticks([])
    for ax, results_dict in (
        (axis[6:8], constants.plots.FIG8_FREQUENCIES),
        (axis[8:10], constants.plots.FIGS5_FREQUENCIES)
    ):
        ax[0].set_prop_cycle(
            cycler('linestyle',linestyle_cycle) +
            cycler('c', colour_cycle)
        )
        for frequency, response_probability in results_dict.items():
            ax[0].plot(
                (0, 3),
                (response_probability, response_probability),
                # c=constants.plots.PALETTE['black'],
                alpha=0.5,
                label=' '.join((str(frequency), constants.core.TEST_UNITS['frequency'])),
                clip_on=False,
                zorder=100
            )
        ax1_cmap = colors.LinearSegmentedColormap.from_list(
            'fig8_cmap',
            [constants.plots.PALETTE['sky'], '#0c0c3c']
        )
        # axis[3].set_prop_cycle(
        #     cycler(
        #         'linestyle', [':','--','-','-.']
        #     ) +
        #     cycler('marker', ['^','o','x','d']) +
        #     cycler('c', [
        #         axis3_cmap(0.25),
        #         axis3_cmap(0.5),
        #         axis3_cmap(0.75),
        #         axis3_cmap(1-utils.small_float)
        #     ])
        # )
        bar_x = np.arange(len(results_dict))
        for i, amplitude in enumerate(constants.plots.MODEL_AMPLITUDES):
            if amplitude < 1:
                y = np.zeros(len(results_dict)) + 0.015
            elif amplitude > 1:
                y = np.zeros(len(results_dict)) + 1
            else:
                y = [x for x in results_dict.values()]
            # axis[3].plot(
            #     [x for x in results_dict.keys()],
            #     y,
            #     # c=constants.plots.PALETTE['black'],
            #     alpha=0.5,
            #     label=' '.join((str(amplitude), constants.core.TEST_UNITS['amplitude'])),
            #     clip_on=False,
            #     zorder=100
            # )
            ax[1].bar(
                bar_x + i * constants.plots.FIG8_BAR_WIDTH,
                y,
                constants.plots.FIG8_BAR_WIDTH,
                color=ax1_cmap(i*0.25),
                label=f'{amplitude} {constants.core.TEST_UNITS['amplitude']}',
                # clip_on=False,
                # zorder=100
            )
        ax[0].legend(loc='center left', bbox_to_anchor=(1,0.6))
        ax[1].legend(loc='center left', bbox_to_anchor=(1,0.6))
        ax[1].set_xticks(
            bar_x + constants.plots.FIG8_BAR_WIDTH*(1.5),
            [str(x) for x in results_dict.keys()]
        )
    for ax, gate_dfs in (
        (axis[10], gate_dfs_7um),
        (axis[11], gate_dfs_14um)
    ):
        ax.set_prop_cycle(cycler('linestyle',linestyle_cycle))
        for i, (frequency, trace) in enumerate(gate_dfs.items()):
            # trace_x = np.array(i*0.000025 for i,_x in enumerate(trace['y']))
            trace_y = np.array(trace['y'])
            duration_s = 0.2
            x_spacing = duration_s * 1.15
            tick_dt = 0.000025
            tick_dt_ms = tick_dt * constants.core.MILLISECONDS_PER_SECOND
            triggers_s = np.arange(0, duration_s, 1/frequency)
            triggers = [int(x/tick_dt) for x in triggers_s]
            # Define epoch bounds in samples relative to the trigger:
            epoch_end = int(10/tick_dt_ms) + 1
            spikes = {}
            for j, trigger in enumerate(triggers):
                # Create `DataEpoch` object describing current epoch and add to
                # `data_epochs` list:
                trace = np.array(trace_y[trigger:trigger+epoch_end])
                # Detect peaks greater than `threshold` within current epoch and
                # add to `spikes` list:
                peaks, _properties = signal.find_peaks(
                    trace,
                    20,
                    distance=constants.core.MINIMUM_SPIKE_DISTANCE_MS/tick_dt_ms
                )
                spikes[str(j)] = peaks*tick_dt_ms
            plot_trace(
                ax,
                np.arange(stop=(duration_s+tick_dt)*1000, step=tick_dt*1000), # type: ignore
                trace_y,
                (0, int(duration_s/tick_dt)),
                classes.TraceGridCoordinates(
                    '',
                    0,
                    0,
                    0,
                    0,
                    None
                ),
                x_spacing,
                0,
                triggers_s, # type: ignore
                [],
                [i for i,_ in enumerate(triggers)],
                [],
                spikes,
                {},
                frequency,
                -1,
                tick_dt,
                linewidth=constants.plots.TRACE_LINEWIDTH*4,
                c=colour_cycle[i],
                alpha=0.5,
                label=' '.join((str(frequency), constants.core.TEST_UNITS['frequency'])),
                clip_on=False,
                zorder=100
            )
            ax.legend(loc='center left', bbox_to_anchor=(1,0.6))
    apply_setp(axis, 'model_output')
    plt.figure(0)
    plt.suptitle(constants.plots.FIG8_TITLE)
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Figure 8"
        )
    plt.figure(1)
    plt.suptitle(constants.plots.FIGS5_TITLE)
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Supplementary S5"
        )
    return


def unit_type_effects(save_plot: bool = True) -> None:
    df = pd.read_csv(
        constants.core.SAVE_PATHS['root'] +
        constants.core.SAVE_PATHS['plot_data'] +
        constants.core.UNIT_TYPE_NAME + '.csv'
    )

    fig = plt.figure(1, figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.QUANTIFICATION_FIG_HEIGHT
    ))
    gs = fig.add_gridspec(9, 7)
    axis = (
        fig.add_subplot(gs[:2, :2]),
        fig.add_subplot(gs[:2, 2:4]),
        fig.add_subplot(gs[:2, 4:6])
    )

    for ax_i, col in enumerate((
        'c_5 min',
        'response_probability_all_100Hz',
        'cv_normal_mean_100 Hz'
    )):
        ax = axis[ax_i]
        for i, unit_type in enumerate(constants.experiments.UNIT_TYPES_PRINT):
            unit_type_df = df.loc[df['unit_type'] == unit_type]
            y = unit_type_df[col].tolist()
            if col == 'cv_normal_mean_100 Hz':
                y = [n*100 for n in y]
            x = [np.random.normal(i+1, constants.plots.JITTER_AMOUNT_NARROW/2)
                 for n in y]
            ax.scatter(
                x,
                y,
                s = constants.plots.QUANTIFICATION_POINT_SIZE,
                c = constants.plots.PALETTE['black'],
                clip_on=False,
                zorder=100
            )
        ax.set_xticks((1,2,3,4), constants.experiments.UNIT_TYPES_PRINT)
    axis[2].yaxis.set_major_formatter(ticker.FormatStrFormatter('%d%%'))
    apply_setp(axis, 'unit_type_effects')
    plt.tight_layout()
    if save_plot:
        utils.save_plot(
            'paper',
            "Supplementary S2"
        )
    return


def initial_properties(
        load_json: bool | str,
        spikes_df: pd.DataFrame,
        length_predictor: str = 'weight',
        stim_type: str = 'electrical'
) -> dict[str, classes.InitialProperties]:
    json_path = (constants.core.SAVE_PATHS['json_root'] +
                 constants.core.SAVE_PATHS['quantification'])
    props_json = load_quantification_json(
        load_json,
        constants.core.INITIAL_PROPS_JSON_NAME,
        json_path
    )
    if props_json is None:
        print("Calculating length estimation functions...")
        length_functions = {x: length_estimation_function(x, length_predictor)
                            for x in constants.experiments.DISTANCE_COMPONENTS}
        print("Calculating unit initial properties...")
        initial_df = spikes_df.loc[
            (spikes_df['Phase'] == 'conditioning') &
            (spikes_df['Epoch Number'] == 0) &
            (spikes_df['Test Amplitude'] == 1) &
            (spikes_df['Epoch Stimulus'] == stim_type) &
            (spikes_df['Test Stimulus'] == stim_type)
        ]
        properties = {}
        for unit_id in initial_df['Unit ID'].unique():
            animal_metadata = constants.experiments.ANIMAL_DATA[
                unit_id[0:5].upper()
            ]
            x_coord = animal_metadata[length_predictor]
            distance_components = animal_metadata[int(unit_id[-1])]['distance']
            distance = sum(
                length_functions[component](x_coord) * multiplier
                for component, multiplier in distance_components.items()
            )

            unit_df = initial_df.loc[initial_df['Unit ID'] == unit_id]
            amplitudes, latencies = peak_properties(unit_df)
            amplitude = statistics.mean(amplitudes[0])
            latency = statistics.mean(latencies[0])
            properties[unit_id] = classes.InitialProperties(
                amplitude,
                latency,
                distance
            )
        print("Unit initial properties calculated.\n")
        utils.confirm_save(
            json_path,
            (f'{constants.core.INITIAL_PROPS_JSON_NAME}_'
             f'{constants.core.VERSION}.json'),
            utils.convert_to_json_dict(properties)
        )
    else:
        properties = {
            k: classes.InitialProperties(**v)
            for k,v in props_json.items()
        }
    return properties


def normalise_conduction_velocity(
    conduction_velocities: abc.Collection[float]
) -> colors.Normalize:
    return colors.Normalize(
        min(conduction_velocities),
        max(conduction_velocities)
    )


def plot_length_estimates(
        length_predictor: str = 'weight',
        save_plot: bool = True
) -> None:
    fig = plt.figure(figsize=(constants.plots.FIG_WIDTH, 6))
    for i, x in enumerate(constants.experiments.DISTANCE_COMPONENTS):
        length_estimation_function(x, length_predictor, s=5)
    experimental_weights = [x['weight'] for x in
                            constants.experiments.ANIMAL_DATA.values()]
    ax = plt.gca()
    ax.fill_between(
        [min(experimental_weights), max(experimental_weights)],
        [-1, 0-1],
        [71, 71],
        color=constants.plots.PALETTE['light_grey'],
        alpha=0.5,
        ec=constants.plots.PALETTE['dark_grey']
    )
    plt.suptitle('Recording distance estimation')
    ax.legend(loc='center left', bbox_to_anchor=(1,0.6))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(2))
    apply_setp([ax], 'length_estimation')
    fig.tight_layout()
    if save_plot:
        utils.save_plot('paper', 'Supplementary S1')
    return


def plot_unit_inspection_figs(
        unit_df: pd.DataFrame,
        title: str,
        test_stim: str = 'electrical',
        save_plot: bool = True
) -> None:
    tests = unit_df['Test'].unique()
    fig = plt.figure(figsize=(
        constants.plots.FIG_WIDTH,
        constants.plots.UNIT_INSPECTION_HEIGHT*len(tests)
    ))
    axis = fig.subplots(len(tests))
    for i, test in enumerate(tests):
        test_df = unit_df.loc[
            (unit_df['Test'] == test) &
            (unit_df['Test Stimulus'] == test_stim)
        ]
        ax = axis if type(axis) is axes.Axes else axis[i] # type: ignore
        bin_width = constants.plots.SHORT_BIN_WIDTH_S
        duration = constants.experiments.DEFAULT_CONDITIONING_DURATION_SECONDS
        label_suffix = None
        if test == 'nine-one':
            filter_col = 'Epoch Stimulus'
            filter_values = [x for x in constants.experiments.STIMULATION_TYPES
                             if x in test_df['Epoch Stimulus'].unique()]
        else:
            test_df = test_df.loc[test_df['Epoch Stimulus'] == test_stim]
            if test == 'long-duration':
                filter_col = 'Test Stimulus'
                filter_values = [test_stim]
                bin_width = constants.plots.LONG_BIN_WIDTH_S
                duration = (constants.experiments.
                            LONG_CONDITIONING_DURATION_SECONDS)
            else:
                filter_col = f'Test {test.capitalize()}'
                filter_values = test_df[filter_col].unique()
                label_suffix = constants.core.TEST_UNITS[test]
        xlim = (0, duration)
        for value in filter_values:
            filtered_df = test_df.loc[
                (test_df['Phase'] == 'conditioning') &
                (test_df[filter_col] == value)
            ]
            binned_data = calculate_binned_probability(
                filtered_df,
                xlim,
                bin_width,
                # duration - 1
            )
            label = (value if label_suffix is None else
                     ' '.join([str(value), label_suffix]))
            ax.plot(
                binned_data.x,
                binned_data.y,
                label=label,
                clip_on=False,
                zorder=100
            )
        ax.set_title(test)
        ax.legend(
            prop={'size': constants.plots.QUANTIFICATION_LEGEND_SIZE},
            loc=constants.plots.QUANTIFICATION_LEGEND_LOC
        )
        ax.set_xlim(xlim)
        apply_setp([ax], 'unit_inspection')
        # ax.set_ylim(constants.plots.PROBABILITY_YLIM)
    plt.suptitle(title)
    fig.tight_layout()
    if save_plot:
        utils.save_plot('unit_inspection', title)
    return


def spikes_table(load_df_json: bool|str = False) -> tuple[pd.DataFrame, bool]:
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
    json_path = (constants.core.SAVE_PATHS['json_root'] +
                 constants.core.SAVE_PATHS['quantification'])
    # Define `spikes_df` so whether it is loaded can be checked later:
    loaded_json = load_quantification_json(
        load_df_json,
        constants.core.SPIKESDF_JSON_NAME,
        json_path
    )
    new_df = loaded_json is None
    if new_df:
        # If spikes_df was not loaded, build new DataFrame:
        print("Building `spikes_df` from saved spikes...")
        spikes_files, _incompatible = json_filenames('spikes')
        all_trials, isi_results = load_spikes_trials(spikes_files)
        spikes_df = tabulate_spikes(all_trials).sort_values([
            'Unit ID',
            'Test',
            'Test Stimulus',
            'Test Frequency',
            'Test Amplitude',
            'Repetition',
            'Phase',
            'Epoch Stimulus',
            'Epoch Number'
        ])
        print("`spikes_df` built.\n")
        # Save new DataFrame:
        df_json = {
            'df': spikes_df.to_dict('records'),
            'isi_results': utils.convert_to_json_dict(isi_results)
        }
        utils.confirm_save(
            json_path,
            (f'{constants.core.SPIKESDF_JSON_NAME}_'
             f'{constants.core.VERSION}.json'),
            df_json,
        )
    else:
        spikes_df = pd.DataFrame.from_dict(loaded_json['df'])
        isi_results = {
            k: classes.ISIResult(**v) for k,v in
            loaded_json['isi_results'].items()
        }
    # Print completion messages and return DataFrame:
    units = spikes_df['Unit ID'].unique()
    print(f'Contains {len(units)} units:')
    print('\n'.join([*units, '']))
    isi_fails = ['%s: %.4f' % (k,v.result) for k,v in isi_results.items() if
                 v.result>constants.core.MAXIMUM_ISI_FAILRATE]
    if isi_fails:
        print(
            f"Warning: {len(isi_fails)} recordings fail ISI threshold "
            f"of {constants.core.MAXIMUM_ISI_FAILRATE}:"
        )
        print('\n'.join(isi_fails))
    highest_isi = max(isi_results.values(), key=lambda x: x.result).result
    print(f'Highest ISI result: {highest_isi}\n')
    return spikes_df, new_df


########################## *BACKEND FUNCTIONS* #########################

def apply_setp(
    axis: abc.Sequence[axes.Axes],
    figure_id: typing.Any
) -> None:
    for i, properties in enumerate(
        constants.plots.QUANTIFICATION_SETP[figure_id]
    ):
        ax = axis[i]
        plt.setp(ax, **properties)
        ax.tick_params(axis='y', which='minor', colors='gray')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    return


def boxplot_with_mean_colour(
        ax: axes.Axes,
        colour: typing.Any = 'black',
        **kwargs
) -> None:
    plot_kwargs = {'showfliers': False} | kwargs
    bp = ax.boxplot(**plot_kwargs) # type: ignore
    for median in bp['medians']:
        median.set_color(colour)
    return


def calculate_binned_probability(
        spikes_df: pd.DataFrame,
        x_lim_s: tuple[float, float],
        bin_width_s: float,
        rolling_bins: bool = False
) -> classes.PlotData:
    try:
        frequency = utils.unique(spikes_df['Epoch Frequency'])
    except AssertionError:
        if spikes_df.empty:
            bins = np.arange(*x_lim_s, bin_width_s)+bin_width_s/2
            return classes.PlotData(
                bins,
                np.zeros(len(bins)),
                None
            )
        # if spikes_df['Epoch Frequency'].unique() == []:
        #     return classes.quantification.PlotData(
        #         np.array([]),
        #         np.array([]),
        #         None
        #     )
        print(
            "Attempted to calculate binned probabilities on a DataFrame which "
            "was insufficiently filtered. Ensure that the DataFrame passed to "
            "this function contains only one epoch frequency."
        )
        raise
    xlim_lower = x_lim_s[0] * frequency
    xlim_upper = x_lim_s[1] * frequency
    bin_width = bin_width_s * frequency
    bins = np.arange(
        xlim_lower,
        xlim_upper-bin_width+1
    ) if rolling_bins else np.arange(
        xlim_lower,
        xlim_upper,
        bin_width
    )+bin_width/2
    response_rates = [[] for _ in bins]
    for unit_id in spikes_df['Unit ID'].unique():
        unit_df = spikes_df.loc[spikes_df['Unit ID'] == unit_id]
        unit_responses = [[] for _ in bins]
        for trial_id in unit_df['Trial ID'].unique():
            trial_df = unit_df.loc[unit_df['Trial ID'] == trial_id]
            try:
                epoch_count = min(
                    utils.unique(trial_df['Total Epochs']),
                    xlim_upper
                )
            except AssertionError as e:
                print(e)
                epoch_count = xlim_upper
            unique_epochs = trial_df['Epoch Number'].unique()
            for i, lower_bound in enumerate(np.arange(
                xlim_lower,
                epoch_count,
                bin_width
            )):
                upper_bound = lower_bound + bin_width
                full_bins = len(np.arange(math.ceil(lower_bound), upper_bound))
                unit_responses[i].append(
                    sum(lower_bound<=float(epoch_id)<upper_bound
                        for epoch_id in unique_epochs)/full_bins
                )
        for i, trials in enumerate(unit_responses):
            try:
                response_rates[i].append(statistics.mean(trials))
            except StatisticsError:
                if trials == []:
                    response_rates[i].append(0)
                else:
                    raise
    error = sem(response_rates)
    return classes.PlotData(
        bins/frequency,
        np.array([statistics.mean(x) for x in response_rates]),
        error
    )


def calculate_unbinned_probability(
        unit_df: pd.DataFrame,
        first_n: int | None = None
) -> list[dict[str, typing.Any]]:
    output = []
    frequency = utils.unique(unit_df['Test Frequency'])
    for trial_id in unit_df['Trial ID'].unique():
        trial_df = unit_df.loc[
            (unit_df['Trial ID'] == trial_id) &
            (unit_df['Phase'] == 'conditioning')
        ]
        if trial_df.empty:
            epoch_count = int(
                frequency * 
                constants.experiments.DEFAULT_CONDITIONING_DURATION_SECONDS
            )
        else:
            epoch_count = utils.unique(trial_df['Total Epochs'])
        unique_epochs = trial_df['Epoch Number'].unique()

        if first_n is None:
            probability = len(unique_epochs)/epoch_count
        else:
            probability = (
                len(unique_epochs)/epoch_count,
                sum(x<first_n for x in unique_epochs)/first_n
            )
        output.append({
            'animal_id': utils.unique(unit_df['Animal ID']),
            'unit_id': utils.unique(unit_df['Unit ID']),
            'test_id': utils.unique(unit_df['Test ID']),
            'trial_id': trial_id,
            'sex': utils.unique(unit_df['Sex']),
            'unit_type': utils.unique(unit_df['Unit Type']),
            'response_probability': probability,
            'epochs': epoch_count
        })
    return output


def colourbar(
        ax: axes.Axes,
        cv_norm: colors.Normalize
) -> None:
    colourbar = plt.colorbar(
        plt.cm.ScalarMappable(cv_norm, cmap),
        ax=ax
    )
    colourbar.set_label(constants.plots.COLOURBAR_LABEL)
    return


@typing.overload
def denormalise(
    y: np.typing.NDArray[np.floating],
    y0: float
) -> np.typing.NDArray[np.floating]: ...

@typing.overload
def denormalise(
    y: float,
    y0: float
) -> float: ...

def denormalise(
    y: np.typing.NDArray[np.floating] | float,
    y0: float
) -> np.typing.NDArray[np.floating] | float:
    return y0*(y+1)


@typing.overload
def double_exponential_decay(
    t: np.typing.NDArray[np.floating],
    a: float,
    k1: float,
    b: float,
    k2: float
) -> np.typing.NDArray[np.floating] | float: ...

@typing.overload
def double_exponential_decay(
    t: float,
    a: float,
    k1: float,
    b: float,
    k2: float
) -> float: ...

def double_exponential_decay(
    t: np.typing.NDArray[np.floating] | float,
    a: float,
    k1: float,
    b: float,
    k2: float
) -> np.typing.NDArray[np.floating] | float:
    for arg in (a, k1, b, k2):
        if arg < 0:
            return 1.0E10
    return a*(np.exp(-k1*t)-1)+b*(np.exp(-k2*t)-1)+1


@typing.overload
def double_exponential_decay_shifted(
    t: np.typing.NDArray[np.floating],
    a: float,
    k1: float,
    b: float,
    k2: float,
    h: float
) -> np.typing.NDArray[np.floating] | float: ...

@typing.overload
def double_exponential_decay_shifted(
    t: float,
    a: float,
    k1: float,
    b: float,
    k2: float,
    h: float
) -> float: ...

def double_exponential_decay_shifted(
    t: np.typing.NDArray[np.floating] | float,
    a: float,
    k1: float,
    b: float,
    k2: float,
    h: float
) -> np.typing.NDArray[np.floating] | float:
    for arg in (a, k1, b, k2, h):
        if arg < 0:
            return 1.0E10
    fit = a*(np.exp(-k1*(t-h))-1)+b*(np.exp(-k2*(t-h))-1)+1
    return np.array([min(x, 1) for x in fit])


@typing.overload
def exponential_decay(
    t: np.typing.NDArray[np.floating],
    a: float,
    b: float
) -> np.typing.NDArray[np.floating] | float: ...

@typing.overload
def exponential_decay(
    t: float,
    a: float,
    b: float
) -> float: ...

def exponential_decay(
    t: np.typing.NDArray[np.floating] | float,
    a: float,
    b: float
) -> np.typing.NDArray[np.floating] | float:
    for arg in (a, b):
        if arg < 0:
            return 1.0E10
    return a*(np.exp(-b*t)-1)+1


def fit_peak_properties(
        amplitudes: abc.Mapping[int, list[float]],
        latencies: abc.Mapping[int, list[float]],
        initial_properties: classes.InitialProperties,
        plot_over_time: tuple[axes.Axes, axes.Axes],
        cv_norm: colors.Normalize
) -> dict[str, float]:
    latency_data = classes.PlotData(
        np.array([x/constants.experiments.LONG_CONDITIONING_FREQUENCY
                  for x in latencies.keys()]),
        np.array([statistics.mean(y) for y in latencies.values()]),
        sem([x for x in latencies.values()])
    )
    latency_fit = lambda t, a, b, c: c*(a*(1-np.exp(-b*t))+1)
    latency_normal_y = normalise(latency_data.y, initial_properties.latency)
    popt_latency_normal, pcov_latency_normal = curve_fit(
        negative_exponential_decay,
        latency_data.x,
        latency_normal_y
    )
    perr_latency_normal = np.sqrt(np.diag(pcov_latency_normal))
    latency_k = classes.FitParameter(
        popt_latency_normal[1],
        perr_latency_normal[1]
    )
    latency_normal_a = classes.FitParameter(
        popt_latency_normal[0],
        perr_latency_normal[0]
    )
    latency_constants = (
        latency_normal_a.value,
        latency_k.value,
        initial_properties.latency
    )

    amplitude_data = classes.PlotData(
        np.array([x/constants.experiments.LONG_CONDITIONING_FREQUENCY
                  for x in amplitudes.keys()]),
        np.array([statistics.mean(y) for y in amplitudes.values()]),
        sem([x for x in amplitudes.values()])
    )
    amplitude_fit = lambda t, a, b: a*t+b
    popt_amplitude, pcov_amplitude = curve_fit(
        amplitude_fit,
        amplitude_data.x,
        amplitude_data.y,
    )
    perr_amplitude = np.sqrt(np.diag(pcov_amplitude))
    amplitude_m = classes.FitParameter(
        popt_amplitude[0],
        perr_amplitude[0]
    )
    amplitude_c = classes.FitParameter(
        popt_amplitude[1],
        perr_amplitude[1]
    )

    plot_over_time[0].scatter(
        latency_data.x,
        latency_data.y,
        s=constants.plots.FIG6_POINT_SIZE,
        c=[initial_properties.conduction_velocity for _ in latency_data.x],
        cmap=cmap,
        norm=cv_norm,
        alpha=constants.plots.QUANTIFICATION_RAW_OPACITY/10,
        clip_on=False,
        zorder=100
    )
    latency_actual_rmse = rmse(
        latency_fit(latency_data.x, *latency_constants),
        latency_data.y
    )
    inspect_fit(
        latency_fit,
        latency_constants,
        latency_data.x,
        ax=plot_over_time[0],
        # error=latency_actual_rmse,
        color=cmap(cv_norm(initial_properties.conduction_velocity)),
        linewidth=constants.plots.QUANTIFICATION_LINEWIDTH
    )
    # inspect_fit(
    #     latency_fit,
    #     latency_constants,
    #     latency_data.x,
    #     ax=plot_over_time[0],
    #     color=constants.plots.PALETTE['light_grey'],
    #     linewidth=constants.plots.QUANTIFICATION_LINEWIDTH/1.2
    # )
    plot_over_time[1].scatter(
        amplitude_data.x,
        amplitude_data.y,
        s=constants.plots.FIG6_POINT_SIZE,
        c=[initial_properties.conduction_velocity for _ in
        amplitude_data.x],
        cmap=cmap,
        norm=cv_norm,
        alpha=constants.plots.QUANTIFICATION_RAW_OPACITY/10,
        clip_on=False,
        zorder=100
    )
    amplitude_predicted = amplitude_fit(amplitude_data.x, *popt_amplitude)
    amplitude_actual_rmse = rmse(
        amplitude_predicted,
        amplitude_data.y
    )
    inspect_fit(
        amplitude_fit,
        popt_amplitude,
        amplitude_data.x,
        ax=plot_over_time[1],
        # error=amplitude_actual_rmse,
        color=cmap(cv_norm(initial_properties.conduction_velocity)),
        linewidth=constants.plots.QUANTIFICATION_LINEWIDTH
    )
    # inspect_fit(
    #     amplitude_fit,
    #     popt_amplitude,
    #     amplitude_data.x,
    #     ax=plot_over_time[1],
    #     # error=amplitude_actual_rmse,
    #     color=constants.plots.PALETTE['light_grey'],
    #     linewidth=constants.plots.QUANTIFICATION_LINEWIDTH/1.2
    # )

    latency_asymptote = denormalise(
        latency_normal_a.value,
        initial_properties.latency
    )
    latency_steady_cv = initial_properties.distance/latency_asymptote
    amplitude_mean = statistics.mean(amplitude_data.y)
    amplitude_normal_y = normalise(
        amplitude_data.y,
        initial_properties.amplitude
    )
    amplitude_normal_predicted = normalise(
        amplitude_predicted,
        initial_properties.amplitude
    )
    return {
        'initial_conduction_velocity': initial_properties.conduction_velocity,
        'initial_amplitude': initial_properties.amplitude,
        'initial_latency': initial_properties.latency,
        'recording_distance': initial_properties.distance,
        'latency_k': latency_k.value,
        'latency_k_std': latency_k.error,
        'latency_normal_a': latency_normal_a.value,
        'latency_normal_a_std': latency_normal_a.error,
        'latency_normal_rmse': rmse(
            negative_exponential_decay(
                latency_data.x,
                *popt_latency_normal
            ), # type: ignore
            latency_normal_y
        ),
        'latency_normal_steady_cv': normalise(
            latency_steady_cv,
            initial_properties.conduction_velocity
        ),
        'latency_actual_a': latency_normal_a.value*initial_properties.latency,
        'latency_actual_c': latency_asymptote,
        'latency_actual_rmse': latency_actual_rmse,
        'latency_actual_steady_cv': latency_steady_cv,
        'amplitude_actual_slope': amplitude_m.value,
        'amplitude_actual_slope_std': amplitude_m.error,
        'amplitude_actual_intercept': amplitude_c.value,
        'amplitude_actual_intercept_std': amplitude_c.error,
        'amplitude_actual_mean': amplitude_mean,
        'amplitude_actual_rmse': amplitude_actual_rmse,
        'amplitude_actual_rsquared': rsquared(
            amplitude_predicted,
            amplitude_data.y
        ),
        'amplitude_normal_slope': (amplitude_m.value/
                                   initial_properties.amplitude),
        'amplitude_normal_intercept': normalise(
            amplitude_c.value,
            initial_properties.amplitude
        ),
        'amplitude_normal_mean': normalise(
            amplitude_mean,
            initial_properties.amplitude
        ),
        'amplitude_normal_rmse': rmse(
            amplitude_normal_predicted,
            amplitude_normal_y
        ),
        'amplitude_normal_rsquared': rsquared(
            amplitude_normal_predicted,
            amplitude_normal_y
        )
    }


def inspect_fit(
        fit: abc.Callable[..., np.typing.NDArray[np.floating]],
        constants: abc.Collection[float],
        x_actual: np.typing.NDArray[np.floating],
        x_min: float | None = None,
        x_max: float | None = None,
        resolution: int = 1000,
        error: float | None = None,
        ax: axes.Axes | None = None,
        **kwargs
) -> None:
    if ax is None:
        ax = plt.gca()
    if x_min is None:
        x_min = min(0, min(x_actual))
    if x_max is None:
        x_max = max(x_actual)
    step = (x_max - x_min) / resolution # type: ignore
    x = np.arange(x_min, x_max, step)
    y = fit(x, *constants)
    plot_kwargs = {
        'clip_on': False,
        'zorder': 101
    } | kwargs
    ax.plot(x, y, **plot_kwargs)
    if error:
        ax.fill_between(x, y-error, y+error, alpha=0.2, **kwargs)
    return


def johnson_neyman(
        ax: axes.Axes,
        xlim: tuple[float | None, float | None] = (None, None),
        extend_boundaries: bool = False
) -> None:
    # import csv from filename in constants as dataframe
    jn_df = pd.read_csv(
        constants.core.SAVE_PATHS['root'] +
        constants.core.SAVE_PATHS['plot_data'] +
        constants.core.JOHNSON_NEYMAN_NAME + '.csv'
    )
    significance_list = (jn_df['Significance'] == 'Significant').tolist()
    significance_indices = (
        utils.detect_edges(significance_list, 'falling')[-1],
        utils.detect_edges(significance_list, 'rising')[-1]
    )
    significance_thresholds = (
        jn_df['frequency'][significance_indices[0]],
        jn_df['frequency'][significance_indices[1]]
    )
    xlim_lower = xlim[0] if xlim[0] is not None else jn_df['frequency'].min()
    xlim_upper = xlim[1] if xlim[1] is not None else jn_df['frequency'].max()
    # Do the following separately for each section of significance plot:
    plot_df = jn_df.loc[
        (jn_df['frequency'] >= xlim_lower) &
        (jn_df['frequency'] <= xlim_upper)
    ]
    significant_lower_df = plot_df.loc[
        plot_df['frequency'] <= significance_thresholds[0]
    ]
    significant_upper_df = plot_df.loc[
        plot_df['frequency'] >= significance_thresholds[1]
    ]
    insignificant_df = plot_df.loc[
        (plot_df['frequency'] >= significance_thresholds[0]) &
        (plot_df['frequency'] <= significance_thresholds[1])
    ]
    text_y = (
        constants.plots.QUANTIFICATION_SETP['sweeps'][-1]['ylim'][1] * 0.88 +
        constants.plots.QUANTIFICATION_SETP['sweeps'][-1]['ylim'][0] * 0.12
    )
    ax.plot(
        insignificant_df['frequency'],
        insignificant_df['Slope of conduction_velocity'],
        color=constants.plots.PALETTE['orange']
    )
    ax.fill_between(
        insignificant_df['frequency'],
        insignificant_df['Lower'],
        insignificant_df['Upper'],
        color=constants.plots.PALETTE['orange'],
        alpha=0.5,
        ec='none'
    )
    ax.text(
        insignificant_df['frequency'].mean(),
        text_y,
        "n.d.",
        ha='center'
    )
    for df in (significant_lower_df, significant_upper_df):
        ax.plot(
            df['frequency'],
            df['Slope of conduction_velocity'],
            color=constants.plots.PALETTE['green']
        )
        ax.fill_between(
            df['frequency'],
            df['Lower'],
            df['Upper'],
            color=constants.plots.PALETTE['green'],
            alpha=0.5,
            ec='none'
        )
        try:
            ax.text(
                df['frequency'].mean(),
                text_y,
                "*",
                ha='center'
            )
        except StatisticsError:
            continue
    for i in significance_indices:
        boundary_row = jn_df.iloc[i]
        y = (constants.plots.QUANTIFICATION_SETP['sweeps'][-1]['ylim'] if
             extend_boundaries else
             (boundary_row['Lower'], boundary_row['Upper']))
        ax.plot(
            (boundary_row['frequency'], boundary_row['frequency']),
            y,
            color=constants.plots.PALETTE['black'],
            linestyle='dotted',
            linewidth=constants.plots.JOHNSON_NEYMAN_LINEWIDTH
        )
    ax.set_xlim(min(0, xlim_lower), xlim_upper)
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(
        constants.plots.JOHNSON_NEYMAN_MINORTICK_SPACING
    ))
    return


def json_filenames(
    json_type: str,
    match: re.Pattern | None = None
) -> tuple[list[str], list[str]]:
    """Gets list of compatible JSON files matching `json_type`."""
    try:
        json_path = (constants.core.SAVE_PATHS['json_root'] +
                     constants.core.SAVE_PATHS[json_type])
    except KeyError:
        json_path = (constants.core.SAVE_PATHS['json_root'] +
                     f'{json_type}\\')
    json_filenames = [
        f for f in next(walk(json_path), (None, None, []))[2]
    ]
    filtered_filenames = (json_filenames if match is None else
                          [f for f in json_filenames if match.search(f)])
    incompatible_files = []
    for filename in filtered_filenames:
        file_version = classes.VersionNumber(filename)
        compatible, message = file_version.iscompatible(utils.current_version)
        print(f"{filename}: {message.capitalize()}")
        if not compatible:
            incompatible_files.append(filename)
    for filename in incompatible_files:
        filtered_filenames.remove(filename)
    return filtered_filenames, incompatible_files


def length_estimation_function(
        y_name: str,
        x_name: str = 'weight',
        inspect: bool = True,
        **kwargs
) -> abc.Callable[[utils.CanDoMathsT], utils.CanDoMathsT]:
    x_data = np.array([x[x_name] for x in
                       constants.experiments.DISTANCE_ESTIMATION_DATA])
    y_data = np.array([x[y_name] for x in
                       constants.experiments.DISTANCE_ESTIMATION_DATA])
    popt, _pcov = curve_fit(lambda t, a, b: ((t-b)**(1/3))/a, x_data, y_data)
    fit = lambda t: ((t-popt[1])**(1/3))/popt[0]
    if inspect:
        plot_kwargs = {
            # 'linewidth': constants.plots.TRACE_LINEWIDTH,
            'label': y_name.capitalize(),
            'clip_on': False,
            'zorder': 100
        } | kwargs
        plt.scatter(x_data, y_data, **plot_kwargs)
        inspect_fit(
            fit,
            (),
            x_data,
            0,
            600,
            alpha=constants.plots.QUANTIFICATION_RAW_OPACITY,
            clip_on=True
        )
        error = rmse(fit(x_data), y_data)
        print('RMSE (%s): %.4f' % (y_name, error))
    return fit


def load_quantification_json(
        keyword: bool | str,
        json_name: str,
        json_path: str
) -> typing.Any | None:
    if keyword:
        print(f"Trying to load `{json_name}` from file...")
        try:
            if keyword is True or keyword == 'current':
                # Build filename for current version:
                json_name_full = f'{json_name}_{constants.core.VERSION}.json'
            elif keyword == 'latest' or keyword == 'recent':
                # Get list of compatible filenames:
                compatible, _incompatible = json_filenames(
                    'quantification',
                    re.compile(rf'^{json_name}')
                )
                # Get list of versions, preserving index in `compatible`:
                versions = [
                    (i, classes.VersionNumber(filename)) for i, filename in
                    enumerate(compatible)
                ]
                # Get index and version of latest version number:
                latest = max(versions, key=lambda x: x[1])
                if keyword == 'recent':
                    # Cull versions later than current:
                    while latest[1] > utils.current_version:
                        versions.remove(latest)
                        latest = max(versions, key=lambda x: x[1])
                # Find filename at index of latest compatible version:
                json_name_full = compatible[latest[0]]
            else:
                # Build filename for specified version:
                json_name_full = f'{json_name}_{classes.VersionNumber(keyword)}.json'
            # Attempt to load JSON at target filename:
            json_obj = json.load(open(json_path+json_name_full))
            print(f"`{json_name}` loaded from file: {json_name_full}")
            return json_obj
        except (IndexError, OSError):
            print(f"Unable to load `{json_name}`: no compatible file.\n")
        except ValueError:
            print(
                f"Unable to load `{json_name}`: no file exists or invalid keyword.\n"
            )
    return None


def load_spikes_trials(
        json_spike_filenames: abc.Iterable[str]
) -> tuple[
    list[list[classes.SpikesTrial]],
    dict[str, classes.ISIResult]
]:
    """Loads list of JSON spike files.

    Converts each loaded file to `SpikesTrial` object and returns all
    `SpikesTrial` objects as list."""
    all_trials = []
    isi_results = {}
    spike_filepath = (constants.core.SAVE_PATHS['json_root'] +
                      constants.core.SAVE_PATHS['spikes'])
    basenames = {saved_pattern.search(x).group('basename') # type: ignore
                 for x in json_spike_filenames}
    # for spike_filename in json_spike_filenames:
    for basename in basenames:
        filenames = [
            x for x in json_spike_filenames if
            saved_pattern.search(x).group('basename')==basename # type: ignore
        ]
        try:
            filename = utils.unique(filenames)
        except AssertionError:
            versions = [classes.VersionNumber(x) for x in filenames]
            latest = max(versions)
            filename = utils.unique([
                x for x in filenames if classes.VersionNumber(x)==latest
            ])
        # Construct path from which JSON file is to be read:
        save_file = spike_filepath + filename
        # Read JSON file and convert to list of `SpikesTrial` objects:
        spikes_dict = json.load(open(save_file))
        assert spikes_dict['json_type'] == 'spikes'
        r_id = saved_pattern.search(filename).group('r_id') # type: ignore
        spikes = [classes.SpikesTrial.from_dict(dictionary)
                  for dictionary in spikes_dict['spikes']]
        isi_result = classes.ISIResult(**spikes_dict['isi_result'])
        all_trials.append(spikes)
        isi_results[r_id] = isi_result
    return all_trials, isi_results


def model_dfs(file_pre, file_post, frequencies) -> dict[int, pd.DataFrame]:
    return {x: pd.read_csv(
        constants.core.SAVE_PATHS['root'] +
        constants.core.SAVE_PATHS['plot_data'] +
        'model\\' + file_pre + str(x) + file_post,
        names=['y']
    ) for x in frequencies}


@typing.overload
def negative_exponential_decay(
    t: np.typing.NDArray[np.floating],
    a: float,
    b: float
) -> np.typing.NDArray[np.floating] | float: ...

@typing.overload
def negative_exponential_decay(
    t: float,
    a: float,
    b: float
) -> float: ...

def negative_exponential_decay(
    t: np.typing.NDArray[np.floating] | float,
    a: float,
    b: float
) -> np.typing.NDArray[np.floating] | float:
    for arg in (a, b):
        if arg < 0:
            return 1.0E10
    return a*(1-np.exp(-b*t))


@typing.overload
def normalise(
    y: np.typing.NDArray[np.floating],
    y0: float
) -> np.typing.NDArray[np.floating]: ...

@typing.overload
def normalise(
    y: float,
    y0: float
) -> float: ...

def normalise(
    y: np.typing.NDArray[np.floating] | float,
    y0: float
) -> np.typing.NDArray[np.floating] | float:
    return (y/y0)-1


def peak_properties(
        unit_df: pd.DataFrame
) -> tuple[dict[int, list[float]], dict[int, list[float]]]:
    per_trial_data = (
        {x:[] for x in unit_df['Epoch Number'].unique()},
        {x:[] for x in unit_df['Epoch Number'].unique()}
    )
    # per_trial_latencies = [[] for _ in unit_df['Epoch Number'].unique()]
    for epoch_id in unit_df['Epoch ID'].unique():
        epoch_df = unit_df.loc[unit_df['Epoch ID'] == epoch_id]
        epoch_number = utils.unique(epoch_df['Epoch Number'])
        per_trial_data[0][epoch_number].append(
            epoch_df['Size (μV)'].mean()
        )
        per_trial_data[1][epoch_number].append(
            epoch_df['Latency (ms)'].mean()
        )
    per_trial_properties = (
        {i:x for i,x in per_trial_data[0].items() if len(x)>0},
        {i:x for i,x in per_trial_data[1].items() if len(x)>0}
    )
    return per_trial_properties


def plot_recording_trace(
    ax: axes.Axes,
    recording: classes.Recording,
    marker: classes.Marker,
    coordinates: classes.TraceGridCoordinates,
    duration_s: float,
    x_spacing: float,
    y_spacing: float,
    spikes_trials: abc.Collection[classes.SpikesTrial]
) -> None:
    tick_dt_ms = recording.tick_dt*constants.core.MILLISECONDS_PER_SECOND
    trigger_threshold_values = spike_extraction.trigger_thresholds(
        recording
    )
    triggers_info = spike_extraction.trigger_info(
        marker,
        recording,
        *trigger_threshold_values
    )
    trace = recording.signal_data[
        marker.start_sample:marker.end_sample
    ]
    trace_filtered = spike_extraction.process_signal(
        trace,
        triggers_info.mech_triggers,
        triggers_info.elec_triggers,
        tick_dt_ms,
        recording.bessel,
        recording.notch
    )
    trace_filtered = spike_extraction.blank_artefacts(
        trace_filtered,
        triggers_info.mech_triggers,
        constants.experiments.PULSEWIDTH_MS['mechanical']+0.5,
        tick_dt_ms
    )
    trace_filtered = spike_extraction.blank_artefacts(
        trace_filtered,
        triggers_info.elec_triggers,
        constants.experiments.PULSEWIDTH_MS['electrical']+0.5,
        tick_dt_ms
    )
    triggers_by_phase = spike_extraction.separate_sweep_phases(
        recording.test,
        triggers_info
    )
    if recording.test == 'nine-one':
        effective_frequency = (triggers_by_phase.frequency +
                               triggers_by_phase.minor_frequency)
        if triggers_by_phase.test_stim == 'electrical':
            frequency_elec = effective_frequency
            frequency_mech = triggers_by_phase.minor_frequency
            zero_tick = triggers_info.elec_triggers[0]
            if (coordinates.start_s%(1/effective_frequency) ==
                (coordinates.start_s+1/effective_frequency)%
                (1/triggers_by_phase.minor_frequency)):
                start_epoch = int(coordinates.start_s*frequency_mech)
                start_tick = triggers_info.mech_triggers[start_epoch]
            else:
                start_epoch = int(coordinates.start_s*triggers_by_phase.frequency)
                start_tick = triggers_info.elec_triggers[start_epoch]
        elif triggers_by_phase.test_stim == 'mechanical':
            frequency_elec = triggers_by_phase.minor_frequency
            frequency_mech = effective_frequency
            zero_tick = triggers_info.mech_triggers[0]
            if (coordinates.start_s%(1/effective_frequency) ==
                (coordinates.start_s+1/effective_frequency)%
                (1/triggers_by_phase.minor_frequency)):
                start_epoch = int(coordinates.start_s*frequency_elec)
                start_tick = triggers_info.elec_triggers[start_epoch]
            else:
                start_epoch = int(coordinates.start_s*triggers_by_phase.frequency)
                start_tick = triggers_info.mech_triggers[start_epoch]
        else:
            raise KeyError(
                f"Test type not recognised ({triggers_by_phase.test_stim})."
            )
    else:
        start_epoch = int(coordinates.start_s*triggers_by_phase.frequency)
        if triggers_by_phase.test_stim == 'electrical':
            frequency_elec = triggers_by_phase.frequency
            frequency_mech = triggers_by_phase.minor_frequency
            zero_tick = triggers_info.elec_triggers[0]
            start_tick = triggers_info.elec_triggers[start_epoch]
        elif triggers_by_phase.test_stim == 'mechanical':
            frequency_elec = triggers_by_phase.minor_frequency
            frequency_mech = triggers_by_phase.frequency
            zero_tick = triggers_info.mech_triggers[0]
            start_tick = triggers_info.mech_triggers[start_epoch]
        else:
            raise KeyError(
                "Test stimulus not recognised "
                f"({triggers_by_phase.test_stim})."
            )
    start_tick -= int(coordinates.start_offset_s/recording.tick_dt)
    prestart_ticks = start_tick - zero_tick
    end_tick = start_tick + int(duration_s/recording.tick_dt)
    # trace_segment = trace_filtered[start_tick:end_tick]
    x_offset = coordinates.grid_x * x_spacing
    # y_offset = coordinates.grid_y * y_spacing
    spikes_list = [x for x in spikes_trials if
                   x.animal_id == recording.animal_id and
                   x.position == recording.position and
                   x.test == recording.test and
                   x.test_stim == triggers_by_phase.test_stim and
                   x.test_frequency == triggers_by_phase.frequency and
                   x.test_frequency_minor == triggers_by_phase.minor_frequency and
                   x.test_amplitude == triggers_by_phase.amplitude]
    assert len(spikes_list) == 1
    elec_triggers = [(x-start_tick)*recording.tick_dt+x_offset for x in
                     triggers_info.elec_triggers if start_tick<=x<end_tick]
    elec_epochs = np.array(
        [i for i,x in enumerate(triggers_info.elec_triggers)
         if start_tick<=x<end_tick]
    )
    elec_spikes = {
        str(i):np.array([x.time_ms - prestart_ticks*tick_dt_ms for x in spike_list])
        for i,spike_list in spikes_list[0].conditioning.electrical.items()
    }
    if recording.test=='nine-one':
        if triggers_by_phase.test_stim=='electrical':
            elec_epochs = np.array([
                x + math.floor(
                    x*triggers_by_phase.minor_frequency/
                    triggers_by_phase.frequency
                ) for x in elec_epochs
            ])
        elif triggers_by_phase.test_stim!='electrical':
            for key in elec_spikes.keys():
                elec_spikes[key] += int(
                    (1/frequency_elec-1/frequency_mech)*constants.core.MILLISECONDS_PER_SECOND
                )
    mech_triggers = [(x-start_tick)*recording.tick_dt+x_offset for x in
                     triggers_info.mech_triggers if start_tick<=x<end_tick]
    mech_epochs = np.array(
        [i for i,x in enumerate(triggers_info.mech_triggers)
         if start_tick<=x<end_tick]
    )
    mech_spikes = {
        str(i):np.array([x.time_ms - prestart_ticks*tick_dt_ms for x in spike_list])
        for i,spike_list in spikes_list[0].conditioning.mechanical.items()
    }
    if recording.test=='nine-one':
        if triggers_by_phase.test_stim=='mechanical':
            mech_epochs = np.array([
                x + math.floor(
                    x*triggers_by_phase.minor_frequency/
                    triggers_by_phase.frequency
                ) for x in mech_epochs
            ])
        elif triggers_by_phase.test_stim!='mechanical':
            for key in mech_spikes.keys():
                mech_spikes[key] += int(
                    (1/frequency_mech-1/frequency_elec)*constants.core.MILLISECONDS_PER_SECOND
                )
    plot_trace(
        ax,
        np.array([i*recording.tick_dt for i,_ in enumerate(trace_filtered)]),
        trace_filtered,
        (start_tick, end_tick),
        coordinates,
        x_spacing,
        y_spacing,
        elec_triggers,
        mech_triggers,
        elec_epochs,
        mech_epochs,
        elec_spikes,
        mech_spikes,
        frequency_elec,
        frequency_mech,
        recording.tick_dt
    )
    return


def plot_trace(
        ax: axes.Axes,
        trace_x: np.typing.NDArray[np.floating],
        trace_y: np.typing.NDArray[np.floating],
        plot_range: tuple[int, int],
        coordinates: classes.TraceGridCoordinates,
        x_spacing: float,
        y_spacing: float,
        elec_triggers: abc.Sequence[float],
        mech_triggers: abc.Sequence[float],
        elec_epochs: abc.Collection[int],
        mech_epochs: abc.Collection[int],
        elec_spikes: abc.Mapping[str, abc.Collection[float]],
        mech_spikes: abc.Mapping[str, abc.Collection[float]],
        elec_frequency: float,
        mech_frequency: float,
        tick_dt: float,
        vertical_label_offset: float | None = None,
        **kwargs
) -> None:
    x_offset = coordinates.grid_x * x_spacing
    y_offset = coordinates.grid_y * y_spacing
    trace_x_segment = (trace_x[plot_range[0]:plot_range[1]] + x_offset -
                       trace_x[plot_range[0]])
    trace_y_segment = trace_y[plot_range[0]:plot_range[1]] + y_offset
    plot_kwargs = {
        'c': constants.plots.PALETTE['light_grey'],
        'linewidth': constants.plots.TRACE_LINEWIDTH,
        'clip_on': False,
        'zorder': 100
    } | kwargs
    ax.plot(
        trace_x_segment,
        trace_y_segment,
        **plot_kwargs
    )
    label_position = max(trace_y_segment) + y_spacing*0.05
    trigger_line_top = min(trace_y_segment) - y_spacing*0.05
    if vertical_label_offset is not None:
        label_position += y_spacing*vertical_label_offset
        trigger_line_top -= y_spacing*vertical_label_offset
    trigger_line_bottom = trigger_line_top - y_spacing*0.1
    axis_line_bottom = statistics.mean((
        min(trace_y_segment),
        min(trace_y_segment),
        trigger_line_top
    ))
    axis_line_top = statistics.mean((
        max(trace_y_segment),
        max(trace_y_segment),
        label_position
    ))
    if coordinates.label is not None:
        ax.text(
            x_offset,
            label_position,
            coordinates.label
        )
        ax.plot(
            (x_offset, x_offset),
            (axis_line_bottom, axis_line_top),
            color=constants.plots.PALETTE['black'],
            linestyle='dotted',
            linewidth=constants.plots.TRACE_LINEWIDTH*4,
            clip_on=False
        )
    for triggers, epochs, spikes, frequency, colour in (
        (elec_triggers, elec_epochs, elec_spikes, elec_frequency, 'sky'),
        (mech_triggers, mech_epochs, mech_spikes, mech_frequency, 'pink')
    ):
        plot_triggers(
            ax,
            trace_x_segment,
            trace_y_segment,
            triggers,
            (trigger_line_bottom, trigger_line_top),
            epochs,
            spikes,
            colour,
            frequency,
            tick_dt,
            plot_kwargs
        )
    return


def plot_triggers(
        ax: axes.Axes,
        trace_x: np.typing.NDArray[np.floating],
        trace_y: np.typing.NDArray[np.floating],
        triggers: abc.Collection[float],
        trigger_y: tuple[float, float],
        epochs: abc.Collection[int],
        spike_times_dict: abc.Mapping[
            str,
            abc.Collection[float]
        ],
        colour: str,
        frequency: float,
        tick_dt: float,
        plot_kwargs: dict[str, typing.Any]
) -> None:
    for trigger in triggers:
        ax.plot(
            (trigger, trigger),
            (trigger_y),
            c=constants.plots.PALETTE[colour],
            linewidth=constants.plots.TRACE_LINEWIDTH*4,
            clip_on=False,
            zorder=101
        )
    for epoch in epochs:
        if str(epoch) not in spike_times_dict.keys():
            continue
        for spike_time_ms in spike_times_dict[str(epoch)]:
            spike_start = int(
                ((spike_time_ms - constants.plots.TRACE_SPIKEWIDTH_MS/2)/
                 constants.core.MILLISECONDS_PER_SECOND +
                 epoch/frequency)/tick_dt
            )
            spike_end = spike_start + int(
                constants.plots.TRACE_SPIKEWIDTH_MS/
                constants.core.MILLISECONDS_PER_SECOND/tick_dt
            )
            ax.plot(
                trace_x[spike_start:spike_end],
                trace_y[spike_start:spike_end],
                **plot_kwargs | {
                    'c': constants.plots.PALETTE[colour],
                    'zorder': 101
                }
            )


def rmse(
        predicted: np.typing.NDArray[np.floating],
        actual: np.typing.NDArray[np.floating]
) -> float:
    return np.sqrt(np.mean(np.square(predicted-actual)))


def rsquared(
        predicted: np.typing.NDArray[np.floating],
        actual: np.typing.NDArray[np.floating]
) -> float:
    return (np.sum(np.square(predicted-actual))/
            np.sum(np.square(actual-np.mean(actual))))


def round_sigfigs(
        x: float,
        direction: typing.Literal['ceil', 'floor'],
        sf: int = 1
) -> float:
    # ! floating point inaccuracy can cause this not to round up in
    # ! specific cases
    # e.g. 20000000000000001 rounds to 20000000000000000
    digits = math.floor(math.log10(abs(x))) - sf + 1
    significant = 10**digits
    if direction == 'ceil':
        x_rounded = round(x+(0.5*(significant-utils.small_float)), -digits)
        return math.ceil(x_rounded/significant)*significant
    elif direction == 'floor':
        x_rounded = round(x-(0.5*(significant-utils.small_float)), -digits)
        return math.floor(x_rounded/significant)*significant
    else:
        raise ValueError(
            f"Expected direction to be `ceil` or `floor` (got {direction})."
        )


def sem(data: abc.Sequence[abc.Collection[float]]) -> np.ndarray | None:
    try:
        return np.array([statistics.stdev(x)/math.sqrt(len(x)) for x in data])
    except StatisticsError:
        return None


def sort_by_conduction_velocity(
        cv_dict: abc.Mapping[str, utils.T],
        unit_ids: abc.Collection[str],
        reverse: bool = False,
        key: abc.Callable = lambda x: x[1]
) -> list[tuple[str, utils.T]]:
    return [(k,v) for k,v in sorted(
        cv_dict.items(),
        key=key,
        reverse=reverse
    ) if k in unit_ids]


def tabulate_spikes(
        spikes: list[list[classes.SpikesTrial]]
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
        'Sex': constants.experiments.SEXES,
        'Test': constants.experiments.TEST_CODES.values(),
        'Test Stimulus': constants.experiments.STIMULATION_TYPES,
        'Phase': constants.experiments.PHASES,
        'Epoch Stimulus': constants.experiments.STIMULATION_TYPES
    }
    return typed_dataframe(
        spikes_df,
        constants.dataframes.SPIKESDF_TYPES,
        categories
    )


def trace_grid(
    ax: axes.Axes,
    recording_map: abc.Mapping[
        str,
        abc.Mapping[int, tuple[
            str,
            abc.Sequence[tuple[str, int, int, float, float, str | None]]
        ]]
    ],
    duration_s: float,
    x_scale: float | None,
    y_spacing: float,
    y_scale: float | None
) -> None:
    x_spacing = duration_s * 1.15
    max_x = 0
    max_y = 0
    for recording_filename, grid_map in recording_map.items():
        recordings = spike_extraction.read_adicht(recording_filename, False)
        for (
            recording_number,
            (spikes_filename, recording_coordinates)
        ) in grid_map.items():
            recording = recordings[recording_number]
            all_trials, _isi_results = load_spikes_trials([spikes_filename])
            assert len(all_trials) == 1
            # # Calculate spike and trigger thresholds:
            comments = [marker.comment for marker in recording.markers]
            for coordinate_vars in recording_coordinates:
                coordinates = classes.TraceGridCoordinates(*coordinate_vars)
                marker = recording.markers[comments.index(coordinates.comment)]
                plot_recording_trace(
                    ax,
                    recording,
                    marker,
                    coordinates,
                    duration_s,
                    x_spacing,
                    y_spacing,
                    all_trials[0]
                )
                max_x = max(max_x, coordinates.grid_x)
                max_y = max(max_y, coordinates.grid_y)
    x_scale_pos = y_spacing * -0.4
    x_scalelabel_pos = x_scale_pos - y_spacing*0.05
    y_scale_pos = x_spacing * -0.025
    y_scalelabel_pos = y_scale_pos - x_spacing*(max_x+1)*0.025
    if x_scale is not None:
        # x_scale_pos = y_spacing * (max_y-0.4)
        ax.plot(
            (0, x_scale),
            (x_scale_pos, x_scale_pos),
            c=constants.plots.PALETTE['black'],
            clip_on=False,
            zorder=101
        )
        ax.text(
            0,
            x_scalelabel_pos,
            f"{int(x_scale*constants.core.MILLISECONDS_PER_SECOND)} ms",
            va='top'
        )
    if y_scale is not None:
        # y_scale_pos = x_spacing * (max_x+0.92)
        y_scale_start = y_spacing * max_y
        ax.plot(
            (y_scale_pos, y_scale_pos),
            (y_scale_start, y_scale_start + y_scale),
            c=constants.plots.PALETTE['black'],
            clip_on=False,
            zorder=101
        )
        ax.text(
            y_scalelabel_pos,
            y_scale_start,
            f"{int(y_scale)} μV",
            ha='left',
            va='bottom',
            rotation='vertical'
        )
    ax.set_xlim(0, x_spacing * (max_x+0.9))
    ax.set_ylim(x_scalelabel_pos, y_spacing * (max_y+0.5))
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    return


def typed_dataframe(
        df: pd.DataFrame,
        noncategorical_types_dict: dict[str, str],
        categorical_types_dict: abc.Mapping[str, abc.Collection]
) -> pd.DataFrame:
    """Applies dtypes to DataFrame columns.
    
    This function is only useful if at least one dtype is categorical.
    Otherwise, simply use `df.astype(noncategorical_types_dict)`.
    """
    # Define categorical dtypes using a dictionary:
    categorical_types = {
        column: pd.CategoricalDtype(categories, False) # type: ignore
        for column, categories in categorical_types_dict.items()
    }
    # Concatenate categorical and non-categorical dtype dictionaries:
    column_types = (noncategorical_types_dict | categorical_types)
    # Return the output DataFrame and its dtypes:
    return df.astype(column_types)
