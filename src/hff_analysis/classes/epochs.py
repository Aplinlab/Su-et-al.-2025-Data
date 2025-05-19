"""Classes for splitting trials into epochs and storing epoch traces.

# Classes
* `DataEpoch` -- data describing a single epoch trace.
* `EpochsTrial` -- data describing all epochs in a recording trial.
* `TriggersTrial` -- list of trigger timings for a recording trial.
"""

import numpy as np
from numpy.typing import NDArray
import typing


class DataEpoch:
    """Data necessary to plot the trace for an epoch.
    
    # Attributes
    * `trace` -- recorded signal as list of voltage values.
    * `start_ms` -- start time of trace relative to epoch trigger.
    * `tick_dt_ms` -- duration of a sampling tick in milliseconds.
    * `phase` -- experimental phase to which epoch belongs.
    * `stimulus` -- stimulus modality delivered during epoch.
    """
    def __init__(
            self,
            trace: NDArray[np.floating],
            start_ms: int | float,
            tick_dt_ms: int | float,
            phase: str,
            stim_type: str
    ):
        self.trace = trace
        self.start_ms = start_ms
        self.tick_dt_ms = tick_dt_ms
        self.phase = phase
        self.stimulus = stim_type
    

class EpochsTrial:
    """Data describing epochs for a recording trial.

    # Methods
    * `from_dict` -- defines instance from dictionary entries.
    """
    def __init__(
            self,
            animal_id: str,
            position: int,
            test: str,
            test_stim: str,
            test_frequency: int | float,
            test_amplitude: int | float,
            repetition: int,
            epochs: list[DataEpoch]
    ):
        self.animal_id =  animal_id
        self.position =  position
        self.test =  test
        self.test_stim =  test_stim
        self.test_frequency =  test_frequency
        self.test_amplitude =  test_amplitude
        self.repetition = repetition
        self.epochs = epochs

    @classmethod
    def from_dict(cls, dictionary: dict[str, typing.Any]):
        return cls(
            dictionary['animal_id'],
            dictionary['position'],
            dictionary['test'],
            dictionary['test_stim'],
            dictionary['test_frequency'],
            dictionary['test_amplitude'],
            dictionary['repetition'],
            [DataEpoch(**list_item) for list_item in dictionary['epochs']]
        )


class TriggersTrial:
    """Collection of trigger timings for a recording trial."""
    def __init__(
            self,
            test: str,
            test_stim: str,
            stim_value: float,
            triggers_mech_cond: list[int],
            triggers_elec_cond: list[int],
            triggers_mech_itlv: list[int],
            triggers_elec_itlv: list[int],
            triggers_mech_rcvr: list[int],
            triggers_elec_rcvr: list[int]
    ):
        self.test = test
        self.test_stim = test_stim
        self.stim_value = stim_value
        self.triggers = {
            'conditioning': {
                'mechanical': triggers_mech_cond,
                'electrical': triggers_elec_cond
            },
            'interleaved': {
                'mechanical': triggers_mech_itlv,
                'electrical': triggers_elec_itlv
            },
            'recovery': {
                'mechanical': triggers_mech_rcvr,
                'electrical': triggers_elec_rcvr
            }
        }
        