"""Custom package for analysis of the Su et al. (2025) dataset.

# Functions
* `unique` -- if input list or DataFrame series only contains one unique
value, returns it. Otherwise, raises AssertionError.

## main_part1
* `load_filereadsettings` -- loads saved variables from earlier run.
* `read_adicht` -- extracts relevant data from LabChart recordings.
* `update_outputs` -- generates updated outputs from FRS files.
* `spikes_info` -- applies signal processing, epoch splitting, and peak
detection to data extracted from LabChart.
* `filter_spikes` -- culls peaks to isolate single-unit responses.
* `plot_clusters` -- plots overlaid peaks and traces from all epochs.
* `save_to_json` -- saves groups of variables to JSON format.

## main_part2
* `spikes_table` -- loads or builds summary of exported spike data.
* `calculate_ssr` -- calculates mean spike rate over entire duration of
each phase, separated by trial and stimulation modality.
* `plot_combined_ssr` -- plots mean spike rate values for each phase and
stimulation modality (i.e. averages all trials in input DataFrame).
* `plot_combined_raster` -- plots rasters for all trials in input.

# Compatibility
This section lists requirements to maintain compatibility with specific
versions. Changes which do not satisfy these requirements must be
accompanied by a major version iteration. Conversely, any of these
conditions which are cumbersome to maintain may be removed when
iterating the major version number.

## Version 1
### v1.1.0
* FRS files include keys for `exclude_frequencies` and
`exclude_amplitudes`, in addition to those present earlier.

### v1.0.0 - 1.0.2
* FRS files only include keys for 'version', 'filename',
'recording_segment', 'epoch_timing_ms', 'threshold', and
'spike_criteria'. Default values must be supplied for other keys.
"""

# TODO Add `repetition` to spikes and epoch objects
# TODO Try to make all *vars -> **vars
# TODO  Update read_frs_v1 if necessary
# TODO Automatically set thresholds using standard deviations
# TODO Apply 200 Hz high-pass filter to sliced epoch before spike detection
# TODO Make cluster plots smaller
# TODO Add discarded units to cluster plots
# TODO  Come up with method for deciding which ones are plotted
# TODO Rethreshold all units
# TODO  Find places where high frequency sets can be included by discarding units
# TODO Improve update function
# TODO Write read-frs-v2
# TODO Work out method for adding comment to start of recording/concatenating two recording segments
# TODO  Necessary to get an extra trial out of HFF08-1
# TODO Lots of docstrings to write and update
# TODO Add animal ages to METADATA
# TODO Improve error handling

from .utils import unique
from .base import (
    load_filereadsettings,
    read_adicht,
    update_outputs,
    spikes_info,
    filter_spikes,
    plot_clusters,
    save_to_json,
    spikes_table,
    calculate_ssr,
    plot_combined_ssr,
    plot_combined_raster
)

__all__ = [
    'unique',
    'load_filereadsettings',
    'read_adicht',
    'update_outputs',
    'spikes_info',
    'filter_spikes',
    'plot_clusters',
    'save_to_json',
    'spikes_table',
    'calculate_ssr',
    'plot_combined_ssr',
    'plot_combined_raster'
]