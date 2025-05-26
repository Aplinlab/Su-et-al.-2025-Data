"""Constants to be used in the HFF analysis pipeline.

This module contains all constants for the current version of this
Python library, and should not include anything else. Legacy constants
can be found in `updater.py`.
"""

VERSION = '2.3.2'


####################### *USER-DEFINED CONSTANTS* #######################
# These constants may be tweaked as required.

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


#*NOTCH FILTER PARAMETERS
NOTCHFILT_F0 = 50
NOTCHFILT_Q = 30
# The amount of time to trim from the recording after applying a notch
# filter, in seconds. The specified duration is trimmed from both
# the beginning and end of the recording.
NOTCHFILT_ARTEFACT_WIDTH_S = 3


#*FILEPATHS
RAW_DATA_PATH = '.\\labchart_raw_data\\'
SAVE_PATHS = {
    'json_root': '.\\outputs\\JSON\\',
    'plot_root': '.\\outputs\\images\\',
    'epochs':  'epochs\\',
    'frs': 'file_read_settings\\',
    'spikes': 'spikes\\',
    'spikes_df': 'dataframes\\'
}
SPIKES_DF_JSON_NAME = 'spikesdf'


#*PLOTS
FIGSIZE = (8, 12)
PALETTE = {
    'black': '#000000',
    'orange': '#e69f00',
    'sky': '#56b4e9',
    'pink': '#cc79a7',
    'green': '#009e73',
    'yellow': '#f0e442',
    'blue': '#0072b2',
    'vermillion': '#d55e00'
}

CLUSTER_LINE_WIDTH = 0.3
CLUSTER_POINT_SIZE = 0.1
CLUSTER_YMIN = -200
CLUSTER_YMAX_SCALE = 1.2

QUANTIFICATION_LEGEND_SIZE = 8
QUANTIFICATION_YLIM_BORDER = 0.05

RASTER_PHASE_SPACING = {
    'interleaved': 3.5,
    'recovery': 7
}
RASTER_POINT_SIZE = 0.1


######################### *BACKEND CONSTANTS* ##########################

TEST_CODES = {
    'freq': 'frequency',
    'ampl': 'amplitude',
    'nine': 'nine-one',
    'long': 'long-duration'
}
EXPERIMENTAL_PHASES = ['conditioning', 'interleaved', 'recovery']
STIMULATION_TYPES = ['mechanical', 'electrical']
OUTPUT_TYPES = ['clusters', 'frs', 'epochs', 'spikes']


#*UNIT CONSTANTS
MILLISECONDS_PER_SECOND = 1000
TEST_UNITS = {
    'frequency': 'Hz',
    'amplitude': '× threshold',
}


#*STIMULATION PARAMETERS

DEFAULT_CONDITIONING_AMPLITUDE = 1.0
DEFAULT_CONDITIONING_FREQUENCY = 100.0
DEFAULT_CONDITIONING_DURATION_SECONDS = 3
LONG_CONDITIONING_FREQUENCY = 50.0
LONG_CONDITIONING_DURATION_SECONDS = 300
INTERLEAVED_DURATION_SECONDS = 3
INTERLEAVED_FREQUENCY = 25
RECOVERY_FREQUENCY = 1
SHORT_RECOVERY_DURATION_SECONDS = 20
LONG_RECOVERY_DURATION_SECONDS = 600

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

VERSION_REGEX = r"(?:(?P<major>\d+)\.(?P<minor>\d+)\.(?P<patch>\d+))"

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
    r"(?:^(?P<name>(?P<a_id>\w{3}\d{2})_pos(?P<pos>[\d]+)_"
    r"(?P<testcode>\w{4})[^.]*)(?P<extension>\..+)?$)"
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
    r"(?:^(?P<testvar>\w+)_mech_(?P<mechval>\d+(.\d+)?)_elec_(?P<elecval>\d+(.\d+)?)$)"
)
NINEONE_REGEX = r"(?:^mech_(?P<mechval>\d)_elec_(?P<elecval>\d)$)"
LONGDURATION_REGEX = r"(?:^long_(?P<stimtype>\w{4})$)"

SAVED_FILENAME_REGEX = (
    r"(?:(?P<savetype>\w+)_(?P<u_id>(?P<a_id>\w{3}\d{2})-(?P<pos>\d+))_"
    r"\[(?P<rep>\d+)-(?P<rec>\d+)\]_(?P<testcode>\w{4})_"
    r"(?P<version>v\d+\.\d+\.\d+)(?P<extension>\..+)?$)"
)
