"""#TODO Write package docstring.
"""

# TODO Lots of docstrings to write and update
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