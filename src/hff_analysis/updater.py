"""Contains legacy code for updating files generated in prior versions.

Each of the variables and functions stored in this module exist for:
* v1.0.0-1.1.0

# Variables
* `filename_pattern` - regex pattern for parsing saved filename.
* `read_frs` - reads FRS file and returns `FileReadSettings` object.
"""

import json
import re

from . import classes


############################ *v1.0.0-1.1.0* ############################

filename_pattern_v1 = re.compile(
    r"(?:(?P<adicht_name>(?P<a_id>\w{3}\d{2})_pos(?P<pos>\d+)_"
    r"(?P<testcode>\w{4})[\w]*)-\[(?P<rep>\d+)-(?P<rec>\d+)\]-"
    r"(?P<savetype>\w+)-(?P<version>\d+\.\d+\.\d+)(?P<extension>\..+)?$)"
)


def read_frs_v1(
        input_filename: str,
        folder: str | None = None,
) -> classes.recordings.FileReadSettings:
    """Reads FRS files from v1.0.0-1.1.0.
    
    **Currently does not check if multiple files relating to the same
    trial are present.** User should manually verify this is not true.
    """
    # Parse the input filename:
    m = filename_pattern_v1.search(input_filename)
    try:
        adicht_name = m.group('adicht_name')
        repetition = m.group('rep')
        recording_segment = int(m.group('rec'))
        version_str = m.group('version')
    except AttributeError as e:
        # AttributeError is raised if m is None (i.e. if regex pattern
        # didn't match)
        raise ValueError(
            f"Filename does not match the expected format ({input_filename})."
        ) from e
    except KeyError as e:
        # KeyError raised if `[test_code]` not found in
        # `constants.TEST_CODES`
        raise KeyError(f"Unable to match test type ({input_filename}).") from e
    # `VersionNumber` objects up to v1.1.0 were not prefixed with 'v',
    # so the version number string found in the input filename is used
    # rather than the object.
    frs_filename = (
        f'{adicht_name}-[{repetition}-{recording_segment}]-FILEREADSETTINGS-'
        f'{version_str}.json'
    )
    if folder:
        filepath = rf'.\{folder}\JSON\file_read_settings\{frs_filename}'
    else:
        filepath = rf'.\outputs\JSON\file_read_settings\{frs_filename}'
    frs_json = json.load(open(filepath))
    # FRS files up to v1.0.2 only include keys for 'version',
    # 'filename', 'recording_segment', 'epoch_timing_ms', 'threshold',
    # and 'spike_criteria'. These keys can be accessed by indexing.
    # FRS files from v1.0.2 - 1.1.0 additionally include keys for
    # `exclude_frequencies` and `exclude_amplitudes`. These keys must be
    # accessed using `frs_json.get(key, [])` to maintain accessibility.
    # All other keys (e.g. `repetition`) are not present at all and
    # values must be supplied from elsewhere.
    exclude_frequencies = frs_json.get('exclude_frequencies', [])
    exclude_amplitudes = frs_json.get('exclude_amplitudes', [])
    frs = classes.recordings.FileReadSettings(
        classes.VersionNumber(frs_json['version']),
        frs_json['filename'],
        repetition,
        frs_json['recording_segment'],
        frs_json['epoch_timing_ms'],
        frs_json['threshold_uV'],
        frs_json['spike_criteria'],
        exclude_frequencies,
        exclude_amplitudes
    )
    # Check file name and contents match:
    assert (
        frs.version == classes.VersionNumber(version_str) and
        frs.filename == adicht_name and
        frs.recording_segment == recording_segment
    ), (
        f"Mismatch between FRS name and contents: {frs_filename}\n"
        f"VERSION: {frs.version==classes.VersionNumber(version_str)}\n"
        f"    FRS name: {classes.VersionNumber(version_str)}\n"
        f"    FRS contents: {frs.version}\n"
        f"FILENAME: {frs.filename==adicht_name}\n"
        f"    FRS name: {adicht_name}\n"
        f"    FRS contents: {frs.filename}\n"
        f"RECORDING SEGMENT: {frs.recording_segment==recording_segment}\n"
        f"    FRS name: {recording_segment}\n"
        f"    FRS contents: {frs.recording_segment}"
    )
    return frs
