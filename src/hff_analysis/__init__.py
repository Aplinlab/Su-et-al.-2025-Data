"""#TODO Write package docstring.
"""

# TODO Classify plots by unit type (FA/SA)
# TODO Improve error handling

from .base import (
    plot_spikes,
    read_adicht,
    save_to_json
)
from .classes.recordings import (
    SpikeCriteria
)
from .constants import (
    EPOCHS_JSON_FOLDER,
    FILEREADSETTINGS_JSON_FOLDER,
    SPIKES_JSON_FOLDER,
    PLOTS_FOLDER,
    CLUSTER_PLOTS_FOLDER
)
from .freq import (
    freqsweep_spikes,
    simple_spikerate_df,
    simple_spikerate_plotdf,
    plot_simple_spikerate
)
# from . import ampl
# from . import nine
# from . import long

__all__ = [
    'plot_spikes',
    'read_adicht',
    'save_to_json',
    'SpikeCriteria',
    'EPOCHS_JSON_FOLDER',
    'FILEREADSETTINGS_JSON_FOLDER',
    'SPIKES_JSON_FOLDER',
    'PLOTS_FOLDER',
    'CLUSTER_PLOTS_FOLDER',
    'freqsweep_spikes',
    'simple_spikerate_df',
    'simple_spikerate_plotdf',
    'plot_simple_spikerate'
]