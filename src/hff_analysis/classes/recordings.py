"""Classes used for reading `.adicht` files and storing data from them
prior to separation into epochs and peak detection.
"""

from numpy import ndarray


class SpikeCriteria:
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
    """To reproduce a saved state in `main_part1.ipynb`:

    0. Load a JSON file containing a `FileReadSettings` object using
    `frs = FileReadSettings(**json.load(open(file_path)))`.
    1. Run STEP 0 (import the `hff_analysis` module).
    2. In STEP 1, set `filename = frs.filename` and `data_segments =
    [frs.recording_segment]`, then run STEP 1.
    3. In STEP 2, set `recording_id = 0`, `threshold_uV =
    frs.threshold_uV`, and `epoch_timing_ms = frs.epoch_timing_ms`, then
    run STEP 2.
    """
    def __init__(
            self,
            filename: str,
            recording_segment: int,
            epoch_timing_ms: tuple[int | float, int | float],
            threshold_uV: int | float,
            spike_criteria: dict[str, SpikeCriteria]
    ):
        self.filename = filename
        self.recording_segment = recording_segment
        self.epoch_timing_ms = epoch_timing_ms
        self.threshold_uV = threshold_uV
        self.spike_criteria = spike_criteria

    @classmethod
    def from_dict(cls, dictionary: dict[str, any]):
        return cls(
            dictionary['filename'],
            dictionary['recording_segment'],
            dictionary['epoch_timing_ms'],
            dictionary['threshold_uV'],
            {key: SpikeCriteria(**value) for key, value in
             dictionary['spike_criteria'].items()}
        )


class Marker:
    """Data from a comment for separating a recording into trials.

    # Properties
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

    # Properties
    * `animal_id` -- Animal-specific identifier made up of three letters
    followed by two digits.
    * `position` -- Integer indicating the position from which the
    recording was made.
    * `test`: String describing the test which was performed.
    * `tick_dt` -- Duration of a single sample in seconds (float).
    * `signal_data` -- NumPy array of floats. Each value is the signal
    voltage in microvolts recorded during the sample in that
    position. Note that the signal is inverted, notch-filtered, and
    trimmed.
    * `mech_triggers` -- NumPy array of floats. Each value is the trigger
    voltage for mechanical stimulation in volts recorded during the
    sample in that position.
    * `elec_triggers` -- NumPy array of floats. Each value is the trigger
    voltage for electrical stimulation in volts recorded during the
    sample in that position.
    * `markers` -- List of Marker objects describing each comment in the
    recording.
    """

    def __init__(
            self,
            animal_id: str,
            position: int,
            test: str,
            tick_dt: float,
            signal_data: ndarray,
            mech_triggers: ndarray,
            elec_triggers: ndarray,
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