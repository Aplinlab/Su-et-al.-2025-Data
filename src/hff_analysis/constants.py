"""Constants to be used in the HFF analysis pipeline.

This module contains all constants for this Python library, and should
not include anything else.
"""


#*USER-DEFINED CONSTANTS
# These constants may be tweaked as required.

# File paths
RAW_DATA_FOLDER = '.\\labchart_raw_data\\'
EPOCHS_JSON_FOLDER = '.\\outputs\\JSON\\epochs\\'
FILEREADSETTINGS_JSON_FOLDER = '.\\outputs\\JSON\\file_read_settings\\'
SPIKES_JSON_FOLDER = '.\\outputs\\JSON\\spikes\\'
PLOTS_FOLDER = '.\\outputs\\plots\\'
CLUSTER_PLOTS_FOLDER = '.\\outputs\\plots\\clusters\\'

# The amount of time to trim from the recording after applying a notch
# filter, in seconds. The specified duration is trimmed from both
# the beginning and end of the recording.
NOTCHFILT_ARTEFACT_WIDTH_S = 3

# Colours to use for plotting.
PLOT_COLOURS = {
    'peaks': '#e59e00ff',
    'conditioning_traces': '#56b4e920',
    'interleaved_traces': '#cc79a720',
    'recovery_traces': '#009e7320'
}

EPOCH_LINE_WIDTH = 0.1
SCATTERPLOT_POINT_SIZE = 1

SIMPLE_PLOT_LEGEND_LOC = 'upper right'
SIMPLE_PLOT_LEGEND_SIZE = 8


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
SIMPLE_FREQUENCY_COLUMNS = [
    'Animal ID',
    'Position',
    'Conditioning Frequency (Hz)',
    'Conditioning Stimulus',
    'FR Conditioning',
    'FR Interleaved - Mechanical',
    'FR Interleaved - Electrical',
    'FR Recovery - Mechanical',
    'FR Recovery - Electrical'
]
SIMPLE_FREQUENCY_TYPES = {
    'Animal ID': 'string',
    'Position': 'uint64',
    'Conditioning Frequency (Hz)': 'uint64',
    'FR Conditioning': 'float64',
    'FR Interleaved - Mechanical': 'float64',
    'FR Interleaved - Electrical': 'float64',
    'FR Recovery - Mechanical': 'float64',
    'FR Recovery - Electrical': 'float64'
}
SIMPLE_FREQUENCY_CATEGORIES = {
    'Conditioning Stimulus': ['Mechanical', 'Electrical']
}

SIMPLE_PLOT_COLUMNS = [
    'FR Conditioning - Mechanical',
    'FR Conditioning - Electrical',
    'FR Interleaved - Mechanical Post-Mechanical',
    'FR Interleaved - Electrical Post-Mechanical',
    'FR Interleaved - Mechanical Post-Electrical',
    'FR Interleaved - Electrical Post-Electrical',
    'FR Recovery - Mechanical Post-Mechanical',
    'FR Recovery - Electrical Post-Mechanical',
    'FR Recovery - Mechanical Post-Electrical',
    'FR Recovery - Electrical Post-Electrical' 
]
SIMPLE_PLOT_TYPES = {
    'FR Conditioning - Mechanical': 'float64',
    'FR Conditioning - Electrical': 'float64',
    'FR Interleaved - Mechanical Post-Mechanical': 'float64',
    'FR Interleaved - Electrical Post-Mechanical': 'float64',
    'FR Interleaved - Mechanical Post-Electrical': 'float64',
    'FR Interleaved - Electrical Post-Electrical': 'float64',
    'FR Recovery - Mechanical Post-Mechanical': 'float64',
    'FR Recovery - Electrical Post-Mechanical': 'float64',
    'FR Recovery - Mechanical Post-Electrical': 'float64',
    'FR Recovery - Electrical Post-Electrical': 'float64'
}

#?Not currently in use
# OUTPUT_COLUMNS = [
#     'Animal ID',
#     'Position',
#     'Test',
#     'Conditioning stimulus',
#     'Conditioning frequency',
#     'Conditioning period',
#     'Conditioning amplitude',
#     'Spikes per trial - conditioning mechanical',
#     'Spikes per trial - conditioning electrical',
#     'Spikes per trial - interleaved mechanical',
#     'Spikes per trial - interleaved electrical',
#     'Spikes per trial - recovery mechanical',
#     'Spikes per trial - recovery electrical'
# ]
# OUTPUT_TYPES = {
#     'Animal ID': 'string',
#     'Position': 'string',
#     'Conditioning frequency': 'uint64',
#     'Conditioning period': 'float64',
#     'Conditioning amplitude': 'float64',
#     'Spikes per trial - conditioning mechanical': 'object',
#     'Spikes per trial - conditioning electrical': 'object',
#     'Spikes per trial - interleaved mechanical': 'object',
#     'Spikes per trial - interleaved electrical': 'object',
#     'Spikes per trial - recovery mechanical': 'object',
#     'Spikes per trial - recovery electrical': 'object'
# }
# OUTPUT_CATEGORIES = {
#     'Test': ['Frequency', 'Amplitude', '9:1', 'Long Duration'],
#     'Conditioning stimulus': ['Mechanical', 'Electrical']
# }


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
FILENAME_REGEX = (
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


#*STIMULATION PARAMETERS

AMPLITUDE_SWEEP_CONDITIONING_FREQUENCY = 100
FREQUENCY_SWEEP_CONDITIONING_AMPLITUDE = 1
SWEEPS_CONDITIONING_PHASE_SECONDS = 3
SWEEPS_INTERLEAVED_EPOCHS_EACH_STIMULUS = 75
SWEEPS_RECOVERY_EPOCHS_EACH_STIMULUS = 19

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