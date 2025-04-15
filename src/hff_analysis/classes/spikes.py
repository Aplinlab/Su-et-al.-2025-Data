"""Classes for storing data about detected spikes."""

import pandas as pd

from hff_analysis import constants


class Spike:
    """Data to describe a single spike."""
    def __init__(
            self,
            epoch_number: int,
            time_ms: float,
            size_uV: float
    ):
        self.epoch_number = epoch_number
        self.time_ms = time_ms
        self.size_uV = size_uV

    def to_df(
            self,
            animal_id: str,
            sex: str,
            position: int,
            unit_id: str,
            unit_type: str,
            test: str,
            test_stim: str,
            test_frequency: float,
            test_amplitude: float,
            test_id: str,
            repetition: int,
            trial_id: str,
            phase: str,
            epoch_stim: str,
            epochs_count: int
    ) -> pd.DataFrame:
        phase_id = '_'.join([trial_id, phase, epoch_stim])
        epoch_id = '_'.join([phase_id, str(self.epoch_number)])
        return pd.DataFrame([[
            animal_id.upper(),
            sex,
            position,
            unit_id,
            unit_type,
            test.capitalize(),
            test_stim.capitalize(),
            test_frequency,
            test_amplitude,
            test_id,
            repetition,
            trial_id,
            phase,
            epoch_stim,
            phase_id,
            self.epoch_number,
            epoch_id,
            epochs_count,
            self.time_ms,
            self.size_uV
        ]], columns=constants.ALL_SPIKES_COLUMNS)


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

    def to_df(
            self,
            animal_id: str,
            sex: str,
            position: int,
            unit_id: str,
            unit_type: str,
            test: str,
            test_stim: str,
            test_frequency: float,
            test_amplitude: float,
            test_id: str,
            repetition: int,
            trial_id: str,
            phase: str,
    ) -> pd.DataFrame:
        common_inputs = [
            animal_id,
            sex,
            position,
            unit_id,
            unit_type,
            test,
            test_stim,
            test_frequency,
            test_amplitude,
            test_id,
            repetition,
            trial_id,
            phase
        ]
        try:
            mech_df = pd.concat(
                [spike.to_df(
                    *common_inputs,
                    'Mechanical',
                    self.epochs_mech
                ) for spike in self.mechanical],
                ignore_index=True
            )
        except ValueError:
            mech_df = pd.DataFrame()

        try:
            elec_df = pd.concat(
                [spike.to_df(
                    *common_inputs,
                    'Electrical',
                    self.epochs_elec
                ) for spike in self.electrical],
                ignore_index=True
            )
        except ValueError:
            elec_df = pd.DataFrame()

        return pd.concat([mech_df, elec_df], ignore_index=True)


class SpikesTrial:
    """Data describing every detected spike for a recording trial."""
    def __init__(
            self,
            animal_id: str,
            position: int,
            test: str,
            test_stim: str,
            test_frequency: float,
            test_amplitude: float,
            spikes_cond: SpikesPhase,
            spikes_itlv: SpikesPhase,
            spikes_rcvr: SpikesPhase
    ):
        self.animal_id = animal_id
        self.position = position
        self.test = test
        self.test_stim = test_stim
        self.test_frequency = test_frequency
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
            dictionary['test_frequency'],
            dictionary['test_amplitude'],
            SpikesPhase.from_dict(dictionary['conditioning']),
            SpikesPhase.from_dict(dictionary['interleaved']),
            SpikesPhase.from_dict(dictionary['recovery'])
        )
    
    def to_df(self, repetition: int) -> pd.DataFrame:
        sex = constants.METADATA[self.animal_id.upper()]['sex']
        unit_type = constants.METADATA[self.animal_id.upper()][self.position]
        
        unit_id = '_'.join([self.animal_id, str(self.position)])
        test_id = '_'.join([
            self.test,
            str(self.test_stim),
            str(self.test_frequency),
            str(self.test_amplitude)
        ])
        trial_id = '_'.join([
            unit_id,
            test_id,
            str(repetition)
        ])
        common_inputs = [
            self.animal_id,
            sex,
            self.position,
            unit_id,
            unit_type,
            self.test,
            self.test_stim,
            self.test_frequency,
            self.test_amplitude,
            test_id,
            repetition,
            trial_id
        ]
        df = pd.concat(
            [phase_obj.to_df(
                *common_inputs,
                phase_str
            ) for phase_obj, phase_str in {
                self.conditioning: 'Conditioning',
                self.interleaved: 'Interleaved',
                self.recovery: 'Recovery'
            }.items()],
            ignore_index=True
        )
        return df
