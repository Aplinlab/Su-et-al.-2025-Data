"""Classes for separating recording trials into epochs and storing
signals which have been separated into epochs.
"""

class DataEpoch:
    """Data necessary to plot the trace for an epoch."""
    def __init__(
            self,
            trace: list[float],
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
    """Collection of `DataEpoch` objects for a recording trial."""
    def __init__(
            self,
            animal_id: str,
            position: int,
            test: str,
            test_stim: str,
            test_freqency: int | float,
            test_amplitude: int | float,
            epochs: list[DataEpoch]
    ):
        self.animal_id =  animal_id
        self.position =  position
        self.test =  test
        self.test_stim =  test_stim
        self.test_freqency =  test_freqency
        self.test_amplitude =  test_amplitude
        self.epochs = epochs

    @classmethod
    def from_dict(cls, dictionary: dict[str, any]):
        return cls(
            dictionary['animal_id'],
            dictionary['position'],
            dictionary['test'],
            dictionary['test_stim'],
            dictionary['test_freqency'],
            dictionary['test_amplitude'],
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