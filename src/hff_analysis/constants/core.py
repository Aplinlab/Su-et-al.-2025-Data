"""Constants to be used in the HFF analysis pipeline.

This module contains all constants for the current version of this
Python library, and should not include anything else. Legacy constants
can be found in `updater.py`.
"""

VERSION = '3.0.1'

# Minimum threshold for spike detection (# of st.devs)
SPIKE_NOISE_FLOOR = 5
MINIMUM_SPIKE_DISTANCE_MS = 0.25
MINIMUM_ISI_MS = 1.0
MAXIMUM_ISI_FAILRATE = 0.01
THRESHOLD_INCREMENT = 5
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

#*UNIT CONSTANTS
MILLISECONDS_PER_SECOND = 1000
TEST_UNITS = {
    'frequency': 'Hz',
    'amplitude': '× threshold',
}

#*FILTER PARAMETERS
BLANKING_PADDING_MS = 0.3
BESSELHIGHPASS_N = 2
BESSELHIGHPASS_FREQUENCY = 100
NOTCHFILT_FREQUENCY = 50
NOTCHFILT_Q = 30
# The amount of time to trim from the recording after applying a notch
# filter, in seconds. The specified duration is trimmed from both
# the beginning and end of the recording.
NOTCHFILT_ARTEFACT_WIDTH_S = 1

#*FILES AND FOLDERS
RAW_DATA_PATH = '.\\labchart_raw_data\\'
SAVE_PATHS = {
    'root': '.\\',
    'json_root': '.\\outputs\\JSON\\',
    'plot_root': '.\\outputs\\images\\',
    'epochs':  'epochs\\',
    'frs': 'file_read_settings\\',
    'spikes': 'spikes\\',
    'quantification': 'quantification\\',
    'unit_inspection': 'unit_inspection\\',
    'paper': 'figures\\',
    'plot_data': 'plot_data\\'
}
INITIAL_PROPS_JSON_NAME = 'initialprops'
SPIKESDF_JSON_NAME = 'spikesdf'
JOHNSON_NEYMAN_NAME = 'johnson_neyman_output'
UNIT_TYPE_NAME = 'unit_types'
MODEL_FREQUENCIES = (10, 50, 200)
SPIKE_EXTRACTION_OUTPUT_TYPES = ('clusters', 'frs', 'epochs', 'spikes')
