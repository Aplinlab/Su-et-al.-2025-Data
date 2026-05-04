"""Classes for splitting trials into epochs and storing epoch traces.

# Classes
* `DataEpoch` -- data describing a single epoch trace.
* `EpochsTrial` -- data describing all epochs in a recording trial.
* `TriggersTrial` -- list of trigger timings for a recording trial.
"""
"""Classes for reading LabChart files and storing extracted data.

# Classes
* `FileNameInfo` -- data extracted from LabChart filename.
* `SpikeCriteria` -- set of criteria for isolating responses.
* `FileReadSettings` -- stores all variables used in `main_part1`.

# Variables
* `adicht_filename_pattern` -- regex pattern for recording filenames.
"""
"""Classes for storing data about detected spikes.

# Classes
* `Spike` -- data describing a single spike.
* `SpikesPhase` -- describes spikes making up an experimental phase.
* `SpikesTrial` -- data for all spikes in a recording trial.
"""
"""Class for storing version numbers.

Implements methods for string conversion, comparison/equality, and
checking compatibility.
"""

import numpy as np
import pandas as pd
import re
import typing

from . import constants


# Regex patterns (compiled here to reduce computations)
# For parsing filenames:
adicht_filename_pattern = re.compile(constants.regex.ADICHT_FILENAME_REGEX)
# For parsing version numbers:
version_pattern = re.compile(constants.regex.VERSION_REGEX)


class VersionNumber:
    def __init__(
            self,
            version_str: str
    ):
        m = version_pattern.search(version_str)
        try:
            self.major = int(m.group('major')) # type: ignore
            self.minor = int(m.group('minor')) # type: ignore
            self.patch = int(m.group('patch')) # type: ignore
        except (AttributeError, ValueError) as e:
            # AttributeError is raised if m is None (i.e. if regex pattern
            # didn't match)
            raise ValueError(
                f"Unexpected error parsing version number ({version_str})."
            ) from e
        
    def __str__(self):
        return f'v{self.major}.{self.minor}.{self.patch}'
 
    def __eq__(self, other: object) -> bool:
        try:
            return (
                self.major == other.major and # type: ignore
                self.minor == other.minor and # type: ignore
                self.patch == other.patch # type: ignore
            )
        except AttributeError:
            return False
 
    def __ne__(self, other: object) -> bool:
        try:
            return (
                self.major != other.major or # type: ignore
                self.minor != other.minor or # type: ignore
                self.patch != other.patch # type: ignore
            )
        except AttributeError:
            return False
        
    def __lt__(self, other: typing.Self) -> bool:
        if self.major < other.major:
            return True
        elif (
            self.major == other.major and
            self.minor < other.minor
        ):
            return True
        elif (
            self.major == other.major and
            self.minor == other.minor and
            self.patch < other.patch
        ):
            return True
        else:
            return False

    def __le__(self, other: typing.Self) -> bool:
        if self.major < other.major:
            return True
        elif (
            self.major == other.major and
            self.minor < other.minor
        ):
            return True
        elif (
            self.major == other.major and
            self.minor == other.minor and
            self.patch <= other.patch
        ):
            return True
        else:
            return False

    def __gt__(self, other: typing.Self) -> bool:
        if self.major > other.major:
            return True
        elif (
            self.major == other.major and
            self.minor > other.minor
        ):
            return True
        elif (
            self.major == other.major and
            self.minor == other.minor and
            self.patch > other.patch
        ):
            return True
        else:
            return False

    def __ge__(self, other: typing.Self) -> bool:
        if self.major > other.major:
            return True
        elif (
            self.major == other.major and
            self.minor > other.minor
        ):
            return True
        elif (
            self.major == other.major and
            self.minor == other.minor and
            self.patch >= other.patch
        ):
            return True
        else:
            return False
        
    def iscompatible(
            self,
            other: typing.Self
    ) -> tuple[bool, str]:
        if self.major != other.major:
            return False, "incompatible"
        elif self.minor != other.minor:
            return True, "compatible (minor version mismatch)"
        else:
            return True, "compatible"


class BoxPoint:
    def __init__(
            self,
            column: float,
            x: float,
            y: float
    ):
        self.column = column
        self.x = x
        self.y = y


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
            trace: np.typing.NDArray[np.floating],
            epoch_number: int,
            phase: str,
            stim_type: str
    ):
        self.trace = trace
        self.epoch_number = epoch_number
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
            test_frequency: float,
            test_frequency_minor: float,
            test_amplitude: float,
            repetition: int,
            stddev: float,
            threshold: float,
            tick_dt_ms: float,
            start_ms: float,
            epochs: list[DataEpoch]
    ):
        self.animal_id =  animal_id
        self.position =  position
        self.test =  test
        self.test_stim =  test_stim
        self.test_frequency =  test_frequency
        self.test_frequency_minor =  test_frequency_minor
        self.test_amplitude =  test_amplitude
        self.repetition = repetition
        self.stddev = stddev
        self.threshold = threshold
        self.tick_dt_ms = tick_dt_ms
        self.start_ms = start_ms
        self.epochs = epochs

    @classmethod
    def from_dict(cls, dictionary: dict[str, typing.Any]):
        return cls(
            dictionary['animal_id'],
            dictionary['position'],
            dictionary['test'],
            dictionary['test_stim'],
            dictionary['test_frequency'],
            dictionary['test_frequency_minor'],
            dictionary['test_amplitude'],
            dictionary['repetition'],
            dictionary['stddev'],
            dictionary['threshold'],
            dictionary['tick_dt_ms'],
            dictionary['start_ms'],
            [DataEpoch(**list_item) for list_item in dictionary['epochs']]
        )


class FilenameInfo:
    """Information about a trial present in LabChart filename.
    
    # Attributes
    * `name` -- LabChart file name without file extension.
    * `extension` -- LabChart file extension.
    * `animal_id` -- three letters followed by two numbers identifying animal.
    * `position` -- index of recording position for animal.
    * `test` -- type of test performed during trial.

    # Methods
    * `from_filename` -- defines instance from LabChart filename.
    """
    def __init__(
            self,
            name: str,
            extension: str,
            animal_id: str,
            position: int,
            test: str
    ):
        self.name = name
        self.extension = extension
        self.animal_id = animal_id
        self.position = position
        self.test = test

    @classmethod
    def from_filename(cls, filename: str):
        # Parse the input filename:
        m = adicht_filename_pattern.search(filename)
        try:
            extension = (
                m.group('extension') if m.group('extension') # type: ignore
                else '.adicht'
            )
            testcode = m.group('testcode') # type: ignore
            return cls(
                m.group('name'), # type: ignore
                extension,
                m.group('a_id'), # type: ignore
                int(m.group('pos')), # type: ignore
                constants.experiments.TEST_CODES[testcode]
            )
        except AttributeError as e:
            # AttributeError is raised if m is None (i.e. if regex pattern
            # didn't match)
            raise ValueError(
                f"Filename does not match the expected format ({filename})."
            ) from e
        except KeyError as e:
            # KeyError raised if `[test_code]` not found in
            # `constants.TEST_CODES`
            raise KeyError(f"Unable to match test type ({filename}).") from e
    

class FileReadSettings:
    """Stores all variables used during `main_part1`.

    Can be used to reproduce an earlier run.
    
    To reproduce a saved state in `main_part1.ipynb`:
    0. Load a JSON file containing a `FileReadSettings` object using
    `frs = FileReadSettings.from_dict(json.load(open(file_path)))`.
    1. Run STEP 0 (import the `hff_analysis` module).
    2. In STEP 1, set `filename = frs.filename`, `repetition =
    frs.repetition` and `data_segments = [frs.recording_segment]`, then
    run STEP 1.
    3. In STEP 2, set `recording_id = 0`, `epoch_timing_ms =
    frs.epoch_timing_ms`, and `threshold_uV = frs.threshold_uV`, then
    run STEP 2.
    4. In STEP 3, set `spike_criteria = frs.spike_criteria`,
    `exclude_frequencies = frs.exclude_frequencies`, and
    `exclude_amplitudes = frs.exclude_amplitudes`, then run STEP 3.
    """
    # TODO Write a function which does the steps described
    def __init__(
            self,
            version: VersionNumber,
            filename: str,
            repetition: int,
            recording_number: int,
            epoch_timing_ms: tuple[float, float],
            skip_superfast: bool,
            spike_criteria: dict[str, dict[str, float | None]],
            exclude_frequencies: list[float],
            exclude_amplitudes: list[float],
            enforce_max_failrate: bool
    ):
        self.version = version
        self.filename = filename
        self.repetition = repetition
        self.recording_number = recording_number
        self.epoch_timing_ms = epoch_timing_ms
        self.skip_superfast = skip_superfast
        self.spike_criteria = spike_criteria
        self.exclude_frequencies = exclude_frequencies
        self.exclude_amplitudes = exclude_amplitudes
        self.enforce_max_failrate = enforce_max_failrate

    @classmethod
    def from_dict(cls, dictionary: dict[str, typing.Any]):
        return cls(
            VersionNumber(dictionary['version']),
            dictionary['filename'],
            dictionary['repetition'],
            dictionary['recording_id'],
            dictionary['epoch_timing_ms'],
            dictionary['skip_superfast'],
            dictionary['spike_criteria'],
            dictionary['exclude_frequencies'],
            dictionary['exclude_amplitudes'],
            dictionary['enforce_max_failrate']
        )


class FitParameter:
    def __init__(
            self,
            value: float,
            error: float
    ):
        self.value = value
        self.error = error


class InitialProperties:
    def __init__(
            self,
            amplitude: float,
            latency: float,
            distance: float,
            **_kwargs
    ):
        self.amplitude = amplitude
        self.latency = latency
        self.distance = distance
        self.conduction_velocity = distance/latency


class ISIResult:
    def __init__(
            self,
            failures: int,
            total: int,
            result: float
    ):
        self.failures = failures
        self.total = total
        self.result = result


class Marker:
    """Data from a comment for separating a recording into trials.

    # Attributes
    * `comment` -- Text of the comment.
    * `start_sample` -- Sample on which the trial starts, relative to
    beginning of trimmed signal.
    * `end_sample` -- Sample on which the trial ends, relative to
    beginning of trimmed signal.
    """
    def __init__(self, comment: str, start_sample: int, end_sample: int):
        self.comment = comment
        self.start_sample = start_sample
        self.end_sample = end_sample


class PlotData:
    def __init__(
            self,
            x: np.typing.NDArray[np.floating],
            y: np.typing.NDArray[np.floating],
            sem: np.typing.NDArray[np.floating] | None
    ):
        self.x = x
        self.y = y
        self.sem = sem


class Recording:
    """A class for variables related to a recording segment within a
    `.adicht` file:

    # Attributes
    * `animal_id` -- Animal-specific identifier made up of three letters
    followed by two digits.
    * `position` -- Integer indicating the position from which the
    recording was made.
    * `test`: String describing the test which was performed.
    * `tick_dt` -- Duration of a single sample in seconds (float).
    * `signal_data` -- NumPy array of floats representing recording from
    signal channel. Values are voltages in microvolts with one value per
    sample. Note that the signal is inverted, notch-filtered, and
    trimmed.
    * `mech_triggers` -- NumPy array of floats representing recording
    from mechanical trigger channel. Values are voltages in volts with
    one value per sample.
    * `elec_triggers` -- NumPy array of floats representing recording
    from electrical trigger channel. Values are voltages in volts with
    one value per sample.
    * `markers` -- List of Marker objects describing each comment in the
    recording.
    * `threshold` -- Spike detection threshold to use. Note that this is
    always calculated from the pre-comment section of the **first record
    of each file**.
    """

    def __init__(
            self,
            animal_id: str,
            position: int,
            test: str,
            tick_dt: float,
            signal_data: np.typing.NDArray[np.floating],
            mech_triggers: np.typing.NDArray[np.floating],
            elec_triggers: np.typing.NDArray[np.floating],
            markers: list[Marker],
            stddev: float,
            bessel: tuple[
                np.typing.NDArray[np.floating],
                np.typing.NDArray[np.floating]
            ],
            notch: tuple[
                np.typing.NDArray[np.floating],
                np.typing.NDArray[np.floating]
            ]
    ):
        self.animal_id = animal_id
        self.position = position
        self.test = test
        self.tick_dt = tick_dt
        self.signal_data = signal_data
        self.mech_triggers = mech_triggers
        self.elec_triggers = elec_triggers
        self.markers = markers
        self.stddev = stddev
        self.bessel = bessel
        self.notch = notch


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
            test_frequency_minor: float,
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
        if phase == 'conditioning':
            amplitude = test_amplitude
            if test == 'nine-one' and epoch_stim != test_stim:
                frequency = test_frequency_minor
            else:
                frequency = test_frequency
        else:
            amplitude = constants.experiments.POSTCONDITIONING_AMPLITUDE
            if phase == 'interleaved':
                frequency = constants.experiments.INTERLEAVED_FREQUENCY
            elif phase == 'recovery':
                frequency = constants.experiments.RECOVERY_FREQUENCY
            else:
                raise KeyError(f"Phase not recognised: {phase}")
        return pd.DataFrame([[
            animal_id,
            sex,
            position,
            unit_id,
            unit_type,
            test,
            test_stim,
            test_frequency,
            test_frequency_minor,
            test_amplitude,
            test_id,
            repetition,
            trial_id,
            phase,
            epoch_stim,
            phase_id,
            self.epoch_number,
            epoch_id,
            frequency,
            amplitude,
            epochs_count,
            self.time_ms,
            self.size_uV
        ]], columns=constants.dataframes.SPIKESDF_COLUMNS)


class SpikeCriteria:
    """Criteria used to isolate responses from detected peaks.
    
    These four criteria constrain peaks to a rectangular box when
    plotted as voltage vs time.
    """
    def __init__(
            self,
            latency_min_ms: float | None,
            latency_max_ms: float | None,
            size_min_uV: float | None,
            size_max_uV: float | None
    ):
        self.latency_min_ms = latency_min_ms
        self.latency_max_ms = latency_max_ms
        self.size_min_uV = size_min_uV
        self.size_max_uV = size_max_uV

    def increment_minimum_size(self, increment: float) -> None:
        if self.size_min_uV is None:
            self.size_min_uV = increment
        else:
            self.size_min_uV += increment


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
            spikes_mech: dict[int, list[Spike]] | dict[str, list[Spike]],
            spikes_elec: dict[int, list[Spike]] | dict[str, list[Spike]]
    ):
        self.epochs_mech = epochs_mech
        self.epochs_elec = epochs_elec
        self.mechanical = spikes_mech
        self.electrical = spikes_elec

    @classmethod
    def from_dict(cls, dictionary: dict[str, typing.Any]):
        # spikes_mech = {}
        # spikes_elec = {}
        # for input_dict, output_dict in (
        #     (dictionary['mechanical'], spikes_mech),
        #     (dictionary['electrical'], spikes_elec)
        # ):
        #     for epoch_id, spike_list in input_dict:
        #         output_dict[epoch_id] = [
        #             Spike(**list_item) for list_item in spike_list
        #         ]
        return cls(
            dictionary['epochs_mech'],
            dictionary['epochs_elec'],
            {k:[Spike(**x) for x in v] for k,v in
             dictionary['mechanical'].items()},
            {k:[Spike(**x) for x in v] for k,v in
             dictionary['electrical'].items()}
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
            test_frequency_minor: float,
            test_amplitude: float,
            test_id: str,
            repetition: int,
            trial_id: str,
            phase: str,
    ) -> pd.DataFrame:
        common_inputs = {
            'animal_id': animal_id,
            'sex': sex,
            'position': position,
            'unit_id': unit_id,
            'unit_type': unit_type,
            'test': test,
            'test_stim': test_stim,
            'test_frequency': test_frequency,
            'test_frequency_minor': test_frequency_minor,
            'test_amplitude': test_amplitude,
            'test_id': test_id,
            'repetition': repetition,
            'trial_id': trial_id,
            'phase': phase
        }
        try:
            mech_df = pd.concat(
                [spike.to_df(
                    **common_inputs,
                    epoch_stim='mechanical',
                    epochs_count=self.epochs_mech
                ) for _, spike_list in self.mechanical.items()
                for spike in spike_list],
                ignore_index=True
            )
        except ValueError:
            mech_df = pd.DataFrame()

        try:
            elec_df = pd.concat(
                [spike.to_df(
                    **common_inputs,
                    epoch_stim='electrical',
                    epochs_count=self.epochs_elec
                ) for _, spike_list in self.electrical.items()
                for spike in spike_list],
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
            test_frequency_minor: float,
            test_amplitude: float,
            repetition: int,
            stddev: float,
            threshold: float,
            mech_criteria: SpikeCriteria,
            elec_criteria: SpikeCriteria,
            spikes_cond: SpikesPhase,
            spikes_itlv: SpikesPhase,
            spikes_rcvr: SpikesPhase,
            spikes_rejected: SpikesPhase,
            spikes_failed: SpikesPhase
    ):
        self.animal_id = animal_id
        self.position = position
        self.test = test
        self.test_stim = test_stim
        self.test_frequency = test_frequency
        self.test_frequency_minor =  test_frequency_minor
        self.test_amplitude = test_amplitude
        self.repetition = repetition
        self.stddev = stddev
        self.threshold = threshold
        self.mech_criteria = mech_criteria
        self.elec_criteria = elec_criteria
        self.conditioning = spikes_cond
        self.interleaved = spikes_itlv
        self.recovery = spikes_rcvr
        self.rejected = spikes_rejected
        self.failed = spikes_failed

    @classmethod
    def from_dict(cls, dictionary: dict[str, typing.Any]):
        return cls(
            dictionary['animal_id'],
            dictionary['position'],
            dictionary['test'],
            dictionary['test_stim'],
            dictionary['test_frequency'],
            dictionary['test_frequency_minor'],
            dictionary['test_amplitude'],
            dictionary['repetition'],
            dictionary['stddev'],
            dictionary['threshold'],
            SpikeCriteria(
                **dictionary['mech_criteria']
            ),
            SpikeCriteria(
                **dictionary['elec_criteria']
            ),
            SpikesPhase.from_dict(dictionary['conditioning']),
            SpikesPhase.from_dict(dictionary['interleaved']),
            SpikesPhase.from_dict(dictionary['recovery']),
            SpikesPhase.from_dict(dictionary['rejected']),
            SpikesPhase.from_dict(dictionary['failed'])
        )
    
    def to_df(self) -> pd.DataFrame:
        sex = constants.experiments.ANIMAL_DATA[self.animal_id.upper()]['sex']
        unit_type = (constants.experiments.ANIMAL_DATA[self.animal_id.upper()]
                     [self.position]['type'])
        unit_id = '_'.join([self.animal_id, str(self.position)])
        test_id = '_'.join([
            self.test,
            str(self.test_stim),
            f'{self.test_frequency}~{self.test_frequency_minor}',
            str(self.test_amplitude)
        ])
        trial_id = '_'.join([
            unit_id,
            test_id,
            str(self.repetition)
        ])
        common_inputs = {
            'animal_id': self.animal_id,
            'sex': sex,
            'position': self.position,
            'unit_id': unit_id,
            'unit_type': unit_type,
            'test': self.test,
            'test_stim': self.test_stim,
            'test_frequency': self.test_frequency,
            'test_frequency_minor': self.test_frequency_minor,
            'test_amplitude': self.test_amplitude,
            'test_id': test_id,
            'repetition': self.repetition,
            'trial_id': trial_id
        }
        df = pd.concat(
            [phase_obj.to_df(
                **common_inputs,
                phase=phase
            ) for phase, phase_obj in (
                ('conditioning', self.conditioning),
                ('interleaved', self.interleaved),
                ('recovery', self.recovery)
            )],
            ignore_index=True
        )
        return df


class TraceGridCoordinates:
    def __init__(
            self,
            comment: str,
            grid_x: int,
            grid_y: int,
            start_s: float,
            start_offset_s: float,
            label: str | None
    ):
        self.comment = comment
        self.grid_x = grid_x
        self.grid_y = grid_y
        self.start_s = start_s
        self.start_offset_s = start_offset_s
        self.label = label


class TriggersInfo:
    def __init__(
            self,
            mech_value: float,
            mech_triggers: list[int],
            elec_value: float,
            elec_triggers: list[int]
    ):
        self.mech_value = mech_value
        self.mech_triggers = mech_triggers
        self.elec_value = elec_value
        self.elec_triggers = elec_triggers


class TriggersTrial:
    """Collection of trigger timings for a recording trial."""
    def __init__(
            self,
            test: str,
            test_stim: str,
            # stim_value: float,
            frequency: float,
            minor_frequency: float,
            amplitude: float,
            triggers_mech_cond: list[int],
            triggers_elec_cond: list[int],
            triggers_mech_itlv: list[int],
            triggers_elec_itlv: list[int],
            triggers_mech_rcvr: list[int],
            triggers_elec_rcvr: list[int]
    ):
        self.test = test
        self.test_stim = test_stim
        # self.stim_value = stim_value
        self.frequency = frequency
        self.minor_frequency = minor_frequency
        self.amplitude = amplitude
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
