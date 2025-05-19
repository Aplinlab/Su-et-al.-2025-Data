"""Classes for storing data about detected spikes.

# Classes
* `Spike` -- data describing a single spike.
* `SpikesPhase` -- describes spikes making up an experimental phase.
* `SpikesTrial` -- data for all spikes in a recording trial.
"""

import pandas as pd
import typing

from hff_analysis import constants


class Spike:
    """Data to describe a single spike.
    
    # Attributes
    * `epoch_number` -- index of epoch within phase and stim type.
    * `time_ms` -- latency in milliseconds.
    * `size_uV` -- amplitude in microvolts.

    # Methods
    * `to_df` -- reformats spike data as DataFrame row.
    """
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
    """Collection of spikes making up an experimental phase.

    Includes mechanical and electrical stimulation as separate lists.

    # Attributes
    * `epochs_mech` -- number of mechanical stimulation epochs.
    * `epochs_elec` -- number of electrical stimulation epochs.
    * `mechanical` -- spike responses to mechanical stimulation.
    * `electrical` -- spike responses to electrical stimulation.

    # Methods
    * `from_dict` -- defines instance from dictionary entries.
    * `to_df` -- reformats data to DataFrame.
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
    def from_dict(cls, dictionary: dict[str, typing.Any]):
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
                    epoch_stim='mechanical',
                    epochs_count=self.epochs_mech
                ) for spike in self.mechanical],
                ignore_index=True
            )
        except ValueError:
            mech_df = pd.DataFrame()

        try:
            elec_df = pd.concat(
                [spike.to_df(
                    *common_inputs,
                    epoch_stim='electrical',
                    epochs_count=self.epochs_elec
                ) for spike in self.electrical],
                ignore_index=True
            )
        except ValueError:
            elec_df = pd.DataFrame()

        return pd.concat([mech_df, elec_df], ignore_index=True)


class SpikesTrial:
    """Data describing every detected spike for a recording trial.

    # Methods
    * `from_dict` -- defines instance from dictionary entries.
    * `to_df` -- reformats data to DataFrame."""
    def __init__(
            self,
            animal_id: str,
            position: int,
            test: str,
            test_stim: str,
            test_frequency: float,
            test_amplitude: float,
            repetition: int,
            spikes_cond: SpikesPhase,
            spikes_itlv: SpikesPhase,
            spikes_rcvr: SpikesPhase
    ):
        self.animal_id = animal_id
        self.position = position
        self.repetition = repetition
        self.test = test
        self.test_stim = test_stim
        self.test_frequency = test_frequency
        self.test_amplitude = test_amplitude
        self.conditioning = spikes_cond
        self.interleaved = spikes_itlv
        self.recovery = spikes_rcvr

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
            SpikesPhase.from_dict(dictionary['conditioning']),
            SpikesPhase.from_dict(dictionary['interleaved']),
            SpikesPhase.from_dict(dictionary['recovery'])
        )
    
    def to_df(self) -> pd.DataFrame:
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
            str(self.repetition)
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
            self.repetition,
            trial_id
        ]
        df = pd.concat(
            [phase_obj.to_df(
                *common_inputs,
                phase=phase_str
            ) for phase_obj, phase_str in {
                self.conditioning: 'conditioning',
                self.interleaved: 'interleaved',
                self.recovery: 'recovery'
            }.items()],
            ignore_index=True
        )
        return df
