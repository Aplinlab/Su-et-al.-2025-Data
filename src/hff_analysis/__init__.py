"""#TODO Write package docstring.
"""

# TODO Classify plots by unit type (FA/SA)
# TODO Improve error handling

from .constants import VERSION
from .base import (
    plot_spikes,
    read_adicht,
    save_to_json
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
    'VERSION',
    'plot_spikes',
    'read_adicht',
    'save_to_json',
    'freqsweep_spikes',
    'simple_spikerate_df',
    'simple_spikerate_plotdf',
    'plot_simple_spikerate'
]