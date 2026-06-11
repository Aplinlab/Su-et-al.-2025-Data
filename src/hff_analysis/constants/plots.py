FIG_WIDTH = 8
FIG_WIDTH_COLOURBARLESS = 6
PALETTE = {
    'black': '#000000',
    'dark_grey': '#666666',
    'light_grey': '#bbbbbb',
    'white': '#ffffff',
    'orange': '#e69f00',
    'sky': '#56b4e9',
    'pink': '#cc79a7',
    'green': '#009e73',
    'yellow': '#f0e442',
    'blue': '#0072b2',
    'vermillion': '#d55e00'
}
FIBRE_TYPE_COLOURS = {
    'FA': 'sky',
    'SA': 'green',
    'Ad': 'orange',
    'Unclassified': 'pink'
}

TRACE_LINEWIDTH = 0.3
TRACE_SPIKEWIDTH_MS = 0.8

CLUSTER_FIG_HEIGHT = 12
CLUSTER_THRESHOLD_LINEWIDTH = 0.2
CLUSTER_POINT_SIZE = 0.1
CLUSTER_YMIN = -200
CLUSTER_YMAX_SCALE = 1.2

SHORT_BIN_WIDTH_S = 0.1
LONG_BIN_WIDTH_S = 1
PROBABILITY_MINORTICK_SPACING = 0.04
PROBABILITY_MINORTICK_SPACING_SMALL = 0.05

COLOURMAP = 'viridis'
QUANTIFICATION_BOXPLOT_WIDTH_NARROW = 0.6
QUANTIFICATION_BOXPLOT_WIDTH_WIDE = 0.8
QUANTIFICATION_FIG_HEIGHT = 10
QUANTIFICATION_LEGEND_SIZE = 8
QUANTIFICATION_LEGEND_LOC = 'lower left'
QUANTIFICATION_LINEWIDTH = 2
QUANTIFICATION_POINT_SIZE = 10
QUANTIFICATION_RAW_OPACITY = 0.5
QUANTIFICATION_SETP = {
    'unit_inspection': [{'ylim': (0, 1)}],
    'long-duration': [
        {'title': "CV: 14.4 ms$^{-1}$"},
        {'title': "CV: 14.4 ms$^{-1}$"},
        {'title': "CV: 28.4 ms$^{-1}$"},
        {'title': "CV: 28.4 ms$^{-1}$"},
        {
            'xlabel': "Time (s)",
            'ylabel': "Mean response probability per bin",
            'ylim': (0, 1)
        },
        {
            'xlabel': "Time (s)",
            'ylabel': "Mean response probability per bin",
            'ylim': (0, 1)
        },
        {
            'xlabel': "Conduction velocity (ms$^{-1}$)",
            'ylabel': "Steady state probability",
            'ylim': (0, 1)
        },
        {
            'xlabel': "Conduction velocity (ms$^{-1}$)",
            'ylabel': r"Rate of decay (s$^{-1}$)",
            # 'ylabel': r"$\tau^{-1}$ (s$^{-1}$)",
            'ylim': (0, 0.4)
        },
        {
            'xlabel': "Conduction velocity (ms$^{-1}$)",
            'ylabel': "Steady state probability",
            'ylim': (0, 1)
        },
        {
            'xlabel': "Conduction velocity (ms$^{-1}$)",
            'ylabel': r"Fast rate of decay (s$^{-1}$)",
            # 'ylabel': r"Fast $\tau^{-1}$ (s$^{-1}$)",
            'ylim': (0, 0.6)
        },
        {
            'xlabel': "Conduction velocity (ms$^{-1}$)",
            'ylabel': r"Slow rate of decay (s$^{-1}$)",
            # 'ylabel': r"Slow $\tau^{-1}$ (s$^{-1}$)",
            'ylim': (0, 0.3)
        }
    ],
    'sweeps': [
        {},
        {},
        {
            'title': "Full 3 seconds",
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': "Mean response probability",
            'ylim': (0, 1)
        },
        {
            'title': "First 30 epochs",
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': "Mean response probability",
            'ylim': (0, 1)
        },
        {
            'xlabel': "Stimulation amplitude (× threshold)",
            'ylabel': "Mean response probability",
            'ylim': (0, 1)
        },
        {
            'title': "Johnson-Neyman plot",
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': "Slope (probability vs\nconduction velocity)",
            'ylim': (-0.05, 0.1)
        }
    ],
    'nine-one': [
        {},
        {},
        {
            'ylabel': "Mean response probability",
            'xlim': (-1.5, 4.5),
            'ylim': (0, 1)
        },
        {
            'ylabel': "Mean response probability",
            'xlim': (0.25, 2.75),
            'ylim': (0, 1)
        }
    ],
    'peak_properties': [
        {
            'xlabel': "Time (s)",
            'ylabel': "Spike latency (ms)",
            'ylim': (1.6, 6.4)
        },
        {
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': r"Mean Δ$_{CV}$",
            'ylim': (-15, 15)
        },
        {
            'xlabel': "Stimulation amplitude (× threshold)",
            'ylabel': r"Mean Δ$_{CV}$",
            'ylim': (-15, 15)
        },
        {
            'xlabel': "Time (s)",
            'ylabel': "Spike amplitude (μV)",
            'ylim': (0, 1200)
        },
        {
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': r"Mean Δ$_{amplitude}$",
            'ylim': (-60, 40)
        },
        {
            'xlabel': "Stimulation amplitude (× threshold)",
            'ylabel': r"Mean Δ$_{amplitude}$",
            'ylim': (-60, 40)
        },
    ],
    'peak_properties_first30': [
        {
            'xlabel': "Time (s)",
            'ylabel': "Spike latency (ms)",
            'ylim': (1.6, 6.4)
        },
        {
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': r"Mean Δ$_{CV}$",
            'ylim': (-10, 4)
        },
        {
            'xlabel': "Stimulation amplitude (× threshold)",
            'ylabel': r"Mean Δ$_{CV}$",
            'ylim': (-10, 15)
        },
        {
            'xlabel': "Time (s)",
            'ylabel': "Spike amplitude (μV)",
            'ylim': (0, 1200)
        },
        {
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': r"Mean Δ$_{amplitude}$",
            'ylim': (-20, 40)
        },
        {
            'xlabel': "Stimulation amplitude (× threshold)",
            'ylabel': r"Mean Δ$_{amplitude}$",
            'ylim': (-60, 30)
        },
    ],
    'model_output': [
        {},
        {},
        {},
        {},
        {},
        {},
        {
            'xlabel': "Time (s)",
            'ylabel': "Response probability",
            'xlim': (0, 3),
            'ylim': (0, 1)
        },
        {
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': "Mean response probability",
            # 'xlim': (0, 200),
            'ylim': (0, 1)
        },
        {
            'xlabel': "Time (s)",
            'ylabel': "Response probability",
            'xlim': (0, 3),
            'ylim': (0, 1)
        },
        {
            'xlabel': "Stimulation frequency (Hz)",
            'ylabel': "Mean response probability",
            # 'xlim': (0, 200),
            'ylim': (0, 1)
        },
        {
            'xlabel': "Time (ms)",
            'ylabel': "Slow potassium\ngating variable",
            'xlim': (0, 200),
            'ylim': (0, 0.6)
        },
        {
            'xlabel': "Time (ms)",
            'ylabel': "Slow potassium\ngating variable",
            'xlim': (0, 200),
            'ylim': (0, 0.6)
        }
    ],
    'length_estimation': [
        {
            'xlabel': "Weight (g)",
            'ylabel': "Length (mm)",
            'xlim': (0, 600),
            'ylim': (0, 70)
        }
    ],
    'unit_type_effects': [
        {
            'title': "5 min 50 Hz",
            'xlabel': "Unit type",
            'ylabel': "Steady state probability",
            'xlim': (0, 5),
            'ylim': (0, 1)
        },
        {
            'title': "3 s 100 Hz",
            'xlabel': "Unit type",
            'ylabel': "Mean response probability",
            'xlim': (0, 5),
            'ylim': (0, 1)
        },
        {
            'title': "3 s 100 Hz",
            'xlabel': "Unit type",
            'ylabel': r"Mean Δ$_{CV}$",
            'xlim': (0, 5),
            'ylim': (-5, 0)
        }
    ]
}
COLOURBAR_LABEL = "First-epoch conduction\nvelocity (ms$^{-1}$)"
JITTER_AMOUNT_NARROW = 0.1
JITTER_AMOUNT_WIDE = 0.15
JOHNSON_NEYMAN_LINEWIDTH = 1
JOHNSON_NEYMAN_MINORTICK_SPACING = 0.01
FIG2_CONDUCTIONVELOCITY_MAJORTICK_SPACING = 10
FIG2_SINGLE_INVTAU_MINORTICK_SPACING = 0.02
FIG2_FAST_INVTAU_MINORTICK_SPACING = 0.04
FIG2_SLOW_INVTAU_MINORTICK_SPACING = 0.02
FIG2_TITLE = "5-minute 50-Hz trains"
FREQUENCY_SWEEP_TITLE = "3-second varied-frequency test"
AMPLITUDE_SWEEP_TITLE = "3-second varied-amplitude test"
FIG4_SLOPEPLOT_MINORTICK_SPACING = 0.0004
FIG5_POINT_SIZE = 20
FIG5_TITLE = "9 electrical:1 mechanical interleaved-pulse trains"
FIG5_ACTUAL_INTERSPACING = 3
FIG5_ACTUAL_INTRASPACING = 0.5
FIG5_ACTUAL_PRIMARY_LABELS = {
    'labels': [
        'Regular\n(100 Hz)',
        'Interleaved\n(9/10)',
        'Regular\n(10 Hz)',
        'Interleaved\n(1/10)'
    ],
    'fontsize': 8
}
FIG5_DIFFERENCES_MINORTICK_SPACING = 0.05
FIG5_DIFFERENCES_YLIM = (-1, 1)
FIG6_POINT_SIZE = 0.5
FIG8_BAR_WIDTH = 0.1
FIG8_TITLE = "MRG model output (7.3-μm sciatic nerve)"
FIGS5_TITLE = "MRG model output (14-μm sciatic nerve)"
UNIT_INSPECTION_HEIGHT = 3

RASTER_PHASE_SPACING = {
    'interleaved': 3.5,
    'recovery': 7
}
RASTER_POINT_SIZE = 0.2

MODEL_AMPLITUDES = (0.8, 1, 1.2, 1.5)
FIG8_FREQUENCIES = {
    10: 1.0,
    25: 0.3333,
    50: 0.1667,
    100: 0.214,
    150: 0.4667,
    200: 1.0
}
FIGS5_FREQUENCIES = {
    10: 1.0,
    25: 0.3333,
    50: 0.1667,
    100: 0.1667,
    150: 0.238,
    200: 0.36
}
