"""Constants to be used in the HFF analysis pipeline.

This module contains all constants for this Python library, and should
not include anything else.
"""

VERSION = '1.1.0'


#*USER-DEFINED CONSTANTS
# These constants may be tweaked as required.

# Filenames and paths
RAW_DATA_PATH = '.\\labchart_raw_data\\'
JSON_PATHS = {
    'rootpath': '.\\outputs\\JSON\\',
    'epochs': {
        'path': 'epochs\\',
        'suffix': 'EPOCHS'
    },
    'frs': {
        'path': 'file_read_settings\\',
        'suffix': 'FILEREADSETTINGS'
    },
    'spikes': {
        'path': 'spikes\\',
        'suffix': 'SPIKES'
    },
    'spikes_df': {
        'path': 'dataframes\\',
        'filename': 'spikes_df'
    }
}
PLOT_PATH = '.\\outputs\\images\\'


# The amount of time to trim from the recording after applying a notch
# filter, in seconds. The specified duration is trimmed from both
# the beginning and end of the recording.
NOTCHFILT_ARTEFACT_WIDTH_S = 3

# Plotting variables
FIGSIZE = (8, 12)

CLUSTER_COLOURS = {
    'peaks': '#e59e00ff',
    'conditioning_traces': '#56b4e950',
    'interleaved_traces': '#cc79a750',
    'recovery_traces': '#009e7350'
}
CLUSTER_LINE_WIDTH = 0.3
CLUSTER_POINT_SIZE = 0.5
CLUSTER_YMIN = -200
CLUSTER_YMAX_SCALE = 1.2

SIMPLE_PLOT_LEGEND_LOC = 'upper right'
SIMPLE_PLOT_LEGEND_SIZE = 8

RASTER_COLOURS = {
    'Mechanical': '#56b4e9ff',
    'Electrical': '#e59e00ff'
}
RASTER_PHASE_SPACING = {
    'Interleaved': 3.5,
    'Recovery': 7
}
RASTER_POINT_SCALE = 0.1

# Metadata describing animals and units
METADATA = {
    'HFF02': {
        'sex': 'M',
        1: 'FA',
        2: 'SA'
    },
    'HFF03': {
        'sex': 'F',
        3: 'FA'
    },
    'HFF04': {
        'sex': 'M',
        2: None
    },
    'HFF05': {
        'sex': 'M',
        1: None,
        2: 'Ad'
    },
    'HFF08': {
        'sex': 'F',
        1: 'FA'
    },
    'HFF10': {
        'sex': 'F',
        1: 'FA'
    },
    'HFF11': {
        'sex': 'F',
        1: 'FA',
        2: 'FA'
    },
    'HFF12': {
        'sex': 'F',
        1: 'Ad'
    },
    'HFF13': {
        'sex': 'F',
        1: 'SA',
        2: 'FA'
    },
    'HFF15': {
        'sex': 'M',
        1: None
    },
    'HFF16': {
        'sex': 'M',
        1: 'SA'
    },
    'HFF19': {
        'sex': 'M',
        2: 'SA',
        3: 'FA',
        5: 'Ad'
    },
    'HFF20': {
        'sex': 'F',
        1: 'SA'
    }
}


### BACKEND CONSTANTS ###


TEST_CODE_CONVERSION_TABLE = {
    'freq': 'frequency',
    'ampl': 'amplitude',
    'nine': 'nine-one',
    'long': 'long-duration'
}
EXPERIMENTAL_PHASES = ['conditioning', 'interleaved', 'recovery']
STIMULATION_TYPES = ['mechanical', 'electrical']


#*UNIT CONVERSION CONSTANTS
MILLISECONDS_PER_SECOND = 1000


#*NOTCH FILTER PARAMETERS
NOTCH_FILTER_F0 = 50
NOTCH_FILTER_Q = 30


#*OUTPUT PARAMETERS
ALL_SPIKES_COLUMNS = [
    'Animal ID',
    'Sex',
    'Position',
    'Unit ID',
    'Unit Type',
    'Test',
    'Test Stimulus',
    'Test Frequency',
    'Test Amplitude',
    'Test ID',
    'Repetition',
    'Trial ID',
    'Phase',
    'Epoch Stimulus',
    'Phase ID',
    'Epoch Number',
    'Epoch ID',
    'Total Epochs',
    'Latency (ms)',
    'Size (μV)'
]
ALL_SPIKES_TYPES = {
    'Animal ID': 'string',
    'Sex': 'string',
    'Position': 'uint64',
    'Unit ID': 'string',
    'Unit Type': 'string',
    'Test Frequency': 'float64',
    'Test Amplitude': 'float64',
    'Test ID': 'string',
    'Repetition': 'uint64',
    'Trial ID': 'string',
    'Phase ID': 'string',
    'Epoch Number': 'uint64',
    'Epoch ID': 'string',
    'Total Epochs': 'uint64',
    'Latency (ms)': 'float64',
    'Size (μV)': 'float64'
}
ALL_SPIKES_CATEGORIES = {
    'Test': ['Frequency', 'Amplitude', '9:1', 'Long Duration'],
    'Test Stimulus': ['Mechanical', 'Electrical'],
    'Phase': ['Conditioning', 'Interleaved', 'Recovery'],
    'Epoch Stimulus': ['Mechanical', 'Electrical']
}

SIMPLE_SPIKERATE_COLUMNS = [
    'Animal ID',
    'Sex',
    'Position',
    'Unit ID',
    'Unit Type',
    'Test',
    'Test Stimulus',
    'Test Frequency',
    'Test Amplitude',
    'Test ID',
    'Repetition',
    'Trial ID',
    'SSR Conditioning',
    'SSR Interleaved - Mechanical',
    'SSR Interleaved - Electrical',
    'SSR Recovery - Mechanical',
    'SSR Recovery - Electrical'
]
SIMPLE_SPIKERATE_TYPES = {
    'Animal ID': 'string',
    'Sex': 'string',
    'Position': 'uint64',
    'Unit ID': 'string',
    'Unit Type': 'string',
    'Test Frequency': 'float64',
    'Test Amplitude': 'float64',
    'Test ID': 'string',
    'Repetition': 'uint64',
    'Trial ID': 'string',
    'SSR Conditioning': 'float64',
    'SSR Interleaved - Mechanical': 'float64',
    'SSR Interleaved - Electrical': 'float64',
    'SSR Recovery - Mechanical': 'float64',
    'SSR Recovery - Electrical': 'float64'
}
SIMPLE_SPIKERATE_CATEGORIES = {
    'Test': ['Frequency', 'Amplitude', '9:1', 'Long Duration'],
    'Test Stimulus': ['Mechanical', 'Electrical'],
}

SSR_PLOT_COLUMNS = [
    'SSR Conditioning - Mechanical',
    'SSR Conditioning - Electrical',
    'SSR Interleaved - Mechanical Post-Mechanical',
    'SSR Interleaved - Electrical Post-Mechanical',
    'SSR Interleaved - Mechanical Post-Electrical',
    'SSR Interleaved - Electrical Post-Electrical',
    'SSR Recovery - Mechanical Post-Mechanical',
    'SSR Recovery - Electrical Post-Mechanical',
    'SSR Recovery - Mechanical Post-Electrical',
    'SSR Recovery - Electrical Post-Electrical' 
]
SSR_PLOT_TYPES = {
    'SSR Conditioning - Mechanical': 'float64',
    'SSR Conditioning - Electrical': 'float64',
    'SSR Interleaved - Mechanical Post-Mechanical': 'float64',
    'SSR Interleaved - Electrical Post-Mechanical': 'float64',
    'SSR Interleaved - Mechanical Post-Electrical': 'float64',
    'SSR Interleaved - Electrical Post-Electrical': 'float64',
    'SSR Recovery - Mechanical Post-Mechanical': 'float64',
    'SSR Recovery - Electrical Post-Mechanical': 'float64',
    'SSR Recovery - Mechanical Post-Electrical': 'float64',
    'SSR Recovery - Electrical Post-Electrical': 'float64'
}


#*REGULAR EXPRESSIONS
# For parsing filenames and recording comments.

# Reads a `.adicht` filename with or without the extension and captures
# the following groups:
# - `name`: Filename without file extension.
# - `id`: Animal-specific identifier made up of three letters followed
#   by two digits.
# - `position`: Position from which the recording was made (digits
#   only).
# - `testcode`: Four-letter code indicating the test which was
#   performed.
# - `extension`: File extension, if the filename includes one.
ADICHT_FILENAME_REGEX = (
    r"(?:(?P<name>^(?P<id>hff\d{2})_pos(?P<position>[\d.]+)_"
    r"(?P<testcode>\w{4})[^.]*)(?P<extension>\..+)?)"
)

# The following expressions parse comments attached to different test
# trains, and capture one or more of the following groups:
# - `mechval`: The value (usually frequency in Hz) of the mechanical
#   stimulus.
# - `elecval`: The value (usually frequency in Hz) of the electrical
#   stimulus.
# - `testvar`: In case of sweeps, captures the variable which was varied
#   (i.e. 'frequency' or 'amplitude').
# - `stimtype`: In case of the long-duration train, the type of stimulus
#   used during the conditioning phase (in our data, this should always
#   be 'electrical').
SWEEP_REGEX = (
    r"(?:^(?P<testvar>\w+)_mech_(?P<mechval>\d+)_elec_(?P<elecval>\d+)$)"
)
NINEONE_REGEX = r"(?:^mech_(?P<mechval>\d)_elec_(?P<elecval>\d)$)"
LONGDURATION_REGEX = r"(?:^long_(?P<stimtype>\w{4})$)"

VERSION_REGEX = r"(?:(?P<major>\d+)\.(?P<minor>\d+)\.(?P<patch>\d+))"


#*STIMULATION PARAMETERS

AMPLITUDE_SWEEP_CONDITIONING_FREQUENCY = 100.0
FREQUENCY_SWEEP_CONDITIONING_AMPLITUDE = 1.0
SHORT_CONDITIONING_DURATION_SECONDS = 3
LONG_CONDITIONING_DURATION_SECONDS = 300
INTERLEAVED_DURATION_SECONDS = 3
SHORT_RECOVERY_DURATION_SECONDS = 20
LONG_RECOVERY_DURATION_SECONDS = 600
INTERLEAVED_EPOCHS_FREQUENCY = 25
RECOVERY_EPOCHS_FREQUENCY = 1

# The number of samples by which the trigger onset is delayed compared
# to the artefact onset (positive values indicate that the trigger onset
# occurs after the artefact onset).
TRIGGER_DELAY_SAMPLES = 4

# The fractional amount below recorded trigger value which triggers
# may be expected to dip.
# In actual testing, it seems that a much higher value (i.e. a narrower
# window for error) is fine. However, the value of 0.9 is large enough
# not to interfere with trigger detection within amplitude sweeps where
# the maximum value is no more than 4x the minimum value, so for the
# purposes of our analysis it should be completely acceptable.
TRIGGER_DETECTION_NOISE_WINDOW = 0.9
