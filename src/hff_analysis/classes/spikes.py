"""Classes for storing data about detected spikes."""

class Spike:
    """Data to describe a single spike."""
    def __init__(
            self,
            epoch_id: int,
            time_ms: float,
            size_uV: float
    ):
        self.epoch_id = epoch_id
        self.time_ms = time_ms
        self.size_uV = size_uV


class SpikesPhase:
    """Collection of spikes making up an experimental phase, including
    mechanical and/or electrical stimulation.
    """
    def __init__(
            self,
            epochs_mech: int,
            epochs_elec: int,
            spikes_mech: list[Spike],
            spikes_elec: list[Spike]
    ):
        self.epochs_mech = epochs_mech
        self.epochs_elec = epochs_elec
        self.mechanical = spikes_mech
        self.electrical = spikes_elec

    @classmethod
    def from_dict(cls, dictionary: dict[str, any]):
        return cls(
            dictionary['epochs_mech'],
            dictionary['epochs_elec'],
            [Spike(**list_item) for list_item in dictionary['mechanical']],
            [Spike(**list_item) for list_item in dictionary['electrical']]
        )


class SpikesTrial:
    """Data describing every detected spike for a recording trial."""
    def __init__(
            self,
            animal_id: str,
            position: int,
            test: str,
            test_stim: str,
            test_freqency: float,
            test_amplitude: float,
            spikes_cond: SpikesPhase,
            spikes_itlv: SpikesPhase,
            spikes_rcvr: SpikesPhase
    ):
        self.animal_id = animal_id
        self.position = position
        self.test = test
        self.test_stim = test_stim
        self.test_freqency = test_freqency
        self.test_amplitude = test_amplitude
        self.conditioning = spikes_cond
        self.interleaved = spikes_itlv
        self.recovery = spikes_rcvr

    @classmethod
    def from_dict(cls, dictionary: dict[str, any]):
        return cls(
            dictionary['animal_id'],
            dictionary['position'],
            dictionary['test'],
            dictionary['test_stim'],
            dictionary['test_freqency'],
            dictionary['test_amplitude'],
            SpikesPhase.from_dict(dictionary['conditioning']),
            SpikesPhase.from_dict(dictionary['interleaved']),
            SpikesPhase.from_dict(dictionary['recovery'])
        )
