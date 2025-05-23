"""Custom package for analysis of the Su et al. (2025) dataset.

# Functions
* `calculate_ssr` -- calculates mean spike rate over entire duration of
each phase, separated by trial and stimulation modality.
* `filter_spikes` -- culls peaks to isolate single-unit responses.
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

# TODO Implement long
# TODO Identify units which don't respond well at 10 Hz (check ssr_df)
# TODO Check that ssr has correct average for units with multiple reps
# TODO Apply 200 Hz high-pass filter to sliced epoch before spike detection
# TODO  New pipeline?
# TODO      1a. Split epochs, apply per-epoch filters
# TODO      1b. Get baseline, calculate SD, detect spikes
# TODO      1c. Export epochs, spikes, baseline, SD
# TODO      2. Test and plot filters
# TODO      3. Apply filters and export
# TODO Automatically set thresholds using standard deviations
# TODO  Use gaps between trials to set baseline - check timing
# TODO Add discarded units to cluster plots
# TODO  Come up with method for deciding which ones are plotted
# TODO Work out method for adding comment to start of recording/concatenating two recording segments
# TODO  Necessary to get an extra trial out of HFF08-1
# TODO  Also would be much better (but not strictly necessary) if we want to include mechanical long-duration trials
# TODO Add animal ages to METADATA

from .spikedetect import (
    load_filereadsettings,
    read_adicht,
    spikes_info
)
from .spikefilter import (
    filter_spikes,
    plot_clusters,
    save_to_json
)
from .spikequantify import (
    calculate_ssr,
    plot_combined_raster,
    plot_combined_ssr,
    spikes_table
)
from .updater import update_outputs
from .utils import (current_version, unique)

__all__ = [
    'calculate_ssr',
    'filter_spikes',
    'load_filereadsettings',
    'plot_clusters',
    'plot_combined_raster',
    'plot_combined_ssr',
    'read_adicht',
    'save_to_json',
    'spikes_info',
    'spikes_table',
    'unique',
    'update_outputs',
    'current_version'
]