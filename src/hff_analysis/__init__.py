"""#TODO Write package docstring.
# Compatibility
This section lists requirements to maintain compatibility with specific
versions. Changes which do not satisfy these requirements must be
accompanied by a major version iteration. Conversely, any of these
conditions which are cumbersome to maintain may be removed when
iterating the major version number.

## v1.0.2
* Any function which reads a saved `epoch` or `frs` file file cannot
require `exclude_frequencies` or `exclude_amplitudes` keys. This can be
achieved by using the following code to access those keys after reading
the JSON file as `json_dict`:
```
exclude_frequencies = json_dict.get('exclude_frequencies', [])
exclude_amplitudes = json_dict.get('exclude_amplitudes', [])
```
"""

# TODO Lots of docstrings to write and update
# TODO Clean up `constants.py`
# TODO Add animal ages to METADATA
# TODO Improve error handling

from .constants import (
    VERSION,
    FIGSIZE
)
from .base import (
    filter_spikes,
    plot_clusters,
    read_adicht,
    save_to_json,
    spikes_table
)
from .freq import (
    freqsweep_spikes,
    simple_spikerate_df,
    simple_spikerate_plotdf,
    plot_simple_spikerate,
    plot_raster
)
# from . import ampl
# from . import nine
# from . import long

__all__ = [
    'VERSION',
    'FIGSIZE',
    'filter_spikes',
    'plot_clusters',
    'read_adicht',
    'save_to_json',
    'spikes_table',
    'freqsweep_spikes',
    'simple_spikerate_df',
    'simple_spikerate_plotdf',
    'plot_simple_spikerate',
    'plot_raster'
]