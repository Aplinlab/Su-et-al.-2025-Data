"""Custom package for analysis of the Su et al. (2025) dataset.

# Functions
* `calculate_ssr` -- calculates mean spike rate over entire duration of
each phase, separated by trial and stimulation modality.
* `filter_trials` -- culls peaks to isolate single-unit responses.
* `load_filereadsettings` -- loads saved variables from earlier run.
* `plot_clusters` -- plots overlaid peaks and traces from all epochs.
* `plot_combined_raster` -- plots rasters for all trials in input.
* `plot_combined_ssr` -- plots mean spike rate values for each phase and
stimulation modality (i.e. averages all trials in input DataFrame).
* `read_adicht` -- extracts relevant data from LabChart recordings.
* `save_to_json` -- saves groups of variables to JSON format.
* `spikes_info` -- applies signal processing, epoch splitting, and peak
detection to data extracted from LabChart.
* `spikes_table` -- loads or builds summary of exported spike data.
* `unique` -- if input list or DataFrame series only contains one unique
value, returns it. Otherwise, raises AssertionError.
* `update_outputs` -- generates updated outputs from FRS files.

# Variables
* `current_version` -- module version stored as `VersionNumber` object.
* `MAXIMUM_ISI_FAILRATE` -- passing score on ISI test (if enforced).

# Compatibility
This section lists requirements to maintain compatibility with specific
versions. Changes which do not satisfy these requirements must be
accompanied by a major version iteration. Conversely, any of these
conditions which are cumbersome to maintain may be removed when
iterating the major version number.

## Version 1
### v1.1.0
* `exclude_frequencies` and `exclude_amplitudes` keys introduced for FRS
files. Default value for earlier files should be `[]`.

## Version 2
### v2.0.0
* `repetition` key introduced for FRS files. Default value for earlier
files should be `0`.
* `json_type` key introduced for all JSON outputs. Value for earlier
files should match specific file.

### v2.4.0
* `test_frequency_minor` key introduced for spikes and epochs files.
Default value for earlier files should be `-1.0`.
"""

# TODO implement check to warn of units which do not pass ISI failrate
# TODO update docstrings and comments
# TODO Identify units which don't respond well at 10 Hz (check ssr_df)
# TODO Stretch goals to get a few more datapoints
#   Concatenate two recording segments (HFF08-1)
#       Also necessary for mechanical long-duration (very few examples)
#   Exclude mechanical nine-one from a particular unit
#       Actually this can be done manually in main_part2
#   Some other stuff


from .spike_extraction import (
    filter_trials,
    load_filereadsettings,
    plot_clusters,
    read_adicht,
    save_to_json,
    spikes_info
)
from .quantification import (
    longduration_probability,
    sweeps_probability,
    nineone_probability,
    peak_property_changes,
    model_output,
    unit_type_effects,
    initial_properties,
    normalise_conduction_velocity,
    plot_length_estimates,
    plot_unit_inspection_figs,
    spikes_table
)
from .updater import update_outputs
from .utils import (current_version, unique)
from .constants.core import MAXIMUM_ISI_FAILRATE

__all__ = (
    'longduration_probability',
    'sweeps_probability',
    'nineone_probability',
    'peak_property_changes',
    'filter_trials',
    'initial_properties',
    'load_filereadsettings',
    'model_output',
    'normalise_conduction_velocity',
    'plot_clusters',
    'plot_length_estimates',
    'plot_unit_inspection_figs',
    'read_adicht',
    'save_to_json',
    'spikes_info',
    'spikes_table',
    'unique',
    'unit_type_effects',
    'update_outputs',
    'current_version',
    'MAXIMUM_ISI_FAILRATE'
)