"""Classes for reading LabChart files and storing extracted data.

# Classes
* `FileNameInfo` -- data extracted from LabChart filename.
* `SpikeCriteria` -- set of criteria for isolating responses.
* `FileReadSettings` -- stores all variables used in `main_part1`.

# Variables
* `adicht_filename_pattern` -- regex pattern for recording filenames.
"""

import numpy as np
from numpy.typing import NDArray
import re
import typing

from hff_analysis import constants
from . import version


# Regex pattern for parsing filenames:
# Compiled once, rather than within a function, to reduce computations.
adicht_filename_pattern = re.compile(constants.ADICHT_FILENAME_REGEX)


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
                constants.TEST_CODES[testcode]
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


class SpikeCriteria:
    """Criteria used to isolate responses from detected peaks.
    
    These four criteria constrain peaks to a rectangular box when
    plotted as voltage vs time.
    """
    def __init__(
            self,
            latency_min_ms: int | float | None,
            latency_max_ms: int | float | None,
            size_min_uV: int | float | None,
            size_max_uV: int | float | None
    ):
        self.latency_min_ms = latency_min_ms
        self.latency_max_ms = latency_max_ms
        self.size_min_uV = size_min_uV
        self.size_max_uV = size_max_uV
    

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
            version: version.VersionNumber,
            filename: str,
            repetition: int,
            recording_segment: int,
            epoch_timing_ms: tuple[int | float, int | float],
            threshold_uV: int | float,
            spike_criteria: dict[str, dict[str, int | float | None]],
            exclude_frequencies: list[int | float] | None = None,
            exclude_amplitudes: list[int | float] | None = None
    ):
        self.version = version
        self.filename = filename
        self.repetition = repetition
        self.recording_segment = recording_segment
        self.epoch_timing_ms = epoch_timing_ms
        self.threshold_uV = threshold_uV
        self.spike_criteria = spike_criteria
        if exclude_frequencies is None:
            self.exclude_frequencies = []
        else:
            self.exclude_frequencies = exclude_frequencies
        if exclude_amplitudes is None:
            self.exclude_amplitudes = []
        else:
            self.exclude_amplitudes = exclude_amplitudes

    @classmethod
    def from_dict(cls, dictionary: dict[str, typing.Any]):
        return cls(
            version.VersionNumber(dictionary['version']),
            dictionary['filename'],
            dictionary['repetition'],
            dictionary['recording_segment'],
            dictionary['epoch_timing_ms'],
            dictionary['threshold_uV'],
            dictionary['spike_criteria'],
            dictionary['exclude_frequencies'],
            dictionary['exclude_amplitudes']
        )


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
    """

    def __init__(
            self,
            animal_id: str,
            position: int,
            test: str,
            tick_dt: float,
            signal_data: NDArray[np.floating],
            mech_triggers: NDArray[np.floating],
            elec_triggers: NDArray[np.floating],
            markers: list[Marker]
    ):
        self.animal_id = animal_id
        self.position = position
        self.test = test
        self.tick_dt = tick_dt
        self.signal_data = signal_data
        self.mech_triggers = mech_triggers
        self.elec_triggers = elec_triggers
        self.markers = markers
        