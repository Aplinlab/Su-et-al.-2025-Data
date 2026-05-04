"""Contains legacy code for updating files generated in prior versions.

* `update_outputs` -- generates updated outputs from FRS files.

# Backend
Each of the following variables and functions exist for:
* v1.0.0-1.1.0

## Functions
* `read_frs` - reads FRS file and returns `FileReadSettings` object.

## Variables
* `filename_pattern` - regex pattern for parsing saved filename.
"""

from collections import abc
import json
from matplotlib import pyplot as plt
import re

from . import classes
from . import constants
from . import spike_extraction
from . import utils


def update_outputs(
        save_outputs: str | abc.Container[str],
        input_filenames: abc.Iterable[str],
        input_folder: str | None = None,
        skip_current: bool = True,
        force_overwrite: bool = False
) -> None:
    """Reads outdated FRS files and saves updated outputs.
    
    # Arguments
    * `save_outputs` -- string or list of strings indicating clusters to
    be saved:
      * `'clusters'` -- cluster plots as PDF image.
      * `'frs'` -- file_read_settings as JSON.
      * `'epochs'` -- traces for each epoch as JSON.
      * `'spikes'` -- spike data as JSON.
      * `'all'` -- all of the above.
    * `input_filenames` -- list of filenames which can identify FRS
    files to read. They do not have to the FRS files themselves, but
    can be other output files corresponding to the desired FRS files.
    **However, they must be of the same version as the FRS files to be
    read.**
    * `input_folder` -- path to location where relevant outputs are
    saved. For example, to read FRS files in the folder
    `.\\outputs\\archive\\v1.1.0\\JSON\\file_read_settings\\`, set
    `input_folder` to `'outputs\\archive\\v1.1.0'`.
    * `show_plots` -- whether cluster plots should be displayed.

    # Error Handling:
    * If an updater does not exist matching the version of any filename
    in `input_filenames`, prints an error and skips that filename.
    """
    # Count number of updated files:
    updated_count = 0
    updated_list = []
    for filename in input_filenames:
        # Start new loop
        print('')
        plt.figure()
        # Determine version of input file and pass to correct reader:
        file_version = classes.VersionNumber(filename)
        if file_version < classes.VersionNumber('2.0.0'):
            frs = read_frs_v1(filename, input_folder)
            
        elif file_version < classes.VersionNumber('3.0.0'):
            frs = read_frs_v2(filename, input_folder)
        elif file_version>=utils.current_version and skip_current:
                # Skip loop if file is up-to-date (does not update count):
                print(f"Up to date: {filename}")
                continue
        elif file_version < classes.VersionNumber('4.0.0'):
            frs = read_frs_v3(filename, input_folder)
        else:
            # Skip loop if matching updater not implemented:
            print(f"No updater for {file_version}: {filename}")
            continue
        # Set variables to be passed to other functions:
        if save_outputs == 'all' or save_outputs == ['all']:
            save_outputs = constants.core.SPIKE_EXTRACTION_OUTPUT_TYPES
        try:
            save_clusters = 'clusters' in save_outputs
        except TypeError:
            save_clusters = False
        # Read and process LabChart data:
        recording = spike_extraction.read_adicht(frs.filename)[frs.recording_id]
        spikes, epochs = spike_extraction.spikes_info(
            recording,
            frs.repetition,
            frs.epoch_timing_ms,
            frs.skip_superfast
        )
        (
            filtered_spikes,
            isi_result,
            spike_criteria_mech,
            spike_criteria_elec
        ) = spike_extraction.filter_trials(
            spikes,
            frs.spike_criteria,
            frs.exclude_frequencies,
            frs.exclude_amplitudes,
            frs.enforce_max_failrate
        )
        # Draw cluster plots:
        if save_clusters:
            spike_extraction.plot_clusters(
                filtered_spikes,
                epochs,
                frs.repetition,
                frs.recording_id,
                frs.epoch_timing_ms,
                spike_criteria_mech,
                spike_criteria_elec,
                isi_result,
                frs.exclude_frequencies,
                frs.exclude_amplitudes,
                True
            )
        # Save specified JSON outputs:
        if save_outputs:
            max_failrate = (constants.core.MAXIMUM_ISI_FAILRATE if
                            frs.enforce_max_failrate else None)
            for (file_type, file_contents) in (
                ('frs', {
                    'filename': frs.filename,
                    'epoch_timing_ms': frs.epoch_timing_ms,
                    'skip_superfast': frs.skip_superfast,
                    'spike_criteria': frs.spike_criteria,
                    'exclude_frequencies': frs.exclude_frequencies,
                    'exclude_amplitudes': frs.exclude_amplitudes,
                    'enforce_max_failrate': frs.enforce_max_failrate
                }),
                ('epochs', {
                    'epochs': epochs,
                    'exclude_frequencies': frs.exclude_frequencies,
                    'exclude_amplitudes': frs.exclude_amplitudes
                }),
                ('spikes', {
                    'isi_result': isi_result,
                    'spike_criteria': frs.spike_criteria,
                    'max_failrate': max_failrate,
                    'spikes': filtered_spikes
                })
            ):
                if file_type in save_outputs:
                    spike_extraction.save_to_json(
                        frs.filename,
                        frs.repetition,
                        frs.recording_id,
                        file_type,
                        force_overwrite,
                        **file_contents
                    )
        plt.close('all')
        # Print completion message and update total updated files count:
        print(f"UPDATED from {file_version}: {frs.filename}")
        updated_count += 1
        updated_list.append(frs.filename)
    # Print total updated files count:
    print(f"\nCOMPLETE: {updated_count} files updated.")
    for filename in updated_list:
        print(filename)
    return


############################ *v1.0.0-1.1.0* ############################

filename_pattern_v1 = re.compile(
    r"(?:(?P<adicht_name>\w{3}\d{2}_pos\d+_\w{4}[\w]*)-\[(?P<rep>\d+)-"
    r"(?P<rec>\d+)\]-\w+-(?P<version>\d+\.\d+\.\d+)(?:\..+)?$)"
)


def read_frs_v1(
        input_filename: str,
        folder: str | None = None,
) -> classes.FileReadSettings:
    """Reads FRS files from v1.0.0-1.1.0.
    
    **Currently does not check if multiple files relating to the same
    trial are present.** User should manually verify this is not true.
    """
    # Parse the input filename:
    m = filename_pattern_v1.search(input_filename)
    try:
        adicht_name = m.group('adicht_name') # type: ignore
        repetition = int(m.group('rep')) # type: ignore
        recording_id = int(m.group('rec')) # type: ignore
        version_str = m.group('version') # type: ignore
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
        f'{adicht_name}-[{repetition}-{recording_id}]-FILEREADSETTINGS-'
        f'{version_str}.json'
    )
    if folder is None:
        folder = r'outputs\archive\v1.1.0'
    filepath = rf'.\{folder}\JSON\file_read_settings\{frs_filename}'
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
    frs = classes.FileReadSettings(
        classes.VersionNumber(frs_json['version']),
        frs_json['filename'],
        repetition,
        frs_json['recording_segment'],
        frs_json['epoch_timing_ms'],
        True,
        frs_json['spike_criteria'],
        exclude_frequencies,
        exclude_amplitudes,
        False
    )
    # Check file name and contents match:
    assert (
        frs.version == classes.VersionNumber(version_str) and
        frs.filename == adicht_name and
        frs.recording_id == recording_id
    ), (
        f"Mismatch between FRS name and contents: {frs_filename}\n"
        f"VERSION: {frs.version==classes.VersionNumber(version_str)}\n"
        f"    FRS name: {classes.VersionNumber(version_str)}\n"
        f"    FRS contents: {frs.version}\n"
        f"FILENAME: {frs.filename==adicht_name}\n"
        f"    FRS name: {adicht_name}\n"
        f"    FRS contents: {frs.filename}\n"
        f"RECORDING SEGMENT: {frs.recording_id==recording_id}\n"
        f"    FRS name: {recording_id}\n"
        f"    FRS contents: {frs.recording_id}"
    )
    return frs


############################ *v2.0.0-2.4.2* ############################

filename_pattern_v2 = re.compile(
    r"(?:\w+_(?P<r_id>\w{3}\d{2}-\d+_\[\d+-\d+\]_\w{4})_"
    r"(?P<version>v\d+\.\d+\.\d+)(?:\..+)?$)"
)


def read_frs_v2(
        input_filename: str,
        folder: str | None = None,
) -> classes.FileReadSettings:
    """Reads FRS files from v1.0.0-1.1.0.
    
    **Currently does not check if multiple files relating to the same
    trial are present.** User should manually verify this is not true.
    """
    # Parse the input filename:
    m = filename_pattern_v2.search(input_filename)
    try:
        frs_id = m.group('r_id') # type: ignore
        version = classes.VersionNumber(m.group('version')) # type: ignore
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
    frs_filename = f'frs_{frs_id}_{version}.json'
    if folder is None:
        folder = r'outputs\archive\v2.4.2'
    filepath = rf'.\{folder}\JSON\file_read_settings\{frs_filename}'
    frs_json = json.load(open(filepath))
    frs = classes.FileReadSettings(
        classes.VersionNumber(frs_json['version']),
        frs_json['filename'],
        frs_json['repetition'],
        frs_json['recording_segment'],
        frs_json['epoch_timing_ms'],
        True,
        frs_json['spike_criteria'],
        frs_json['exclude_frequencies'],
        frs_json['exclude_amplitudes'],
        False
    )
    # Check file name and contents match:
    assert (
        frs.version == version
    ), (
        f"Mismatch between FRS name and contents: {frs_filename}\n"
        f"VERSION: {frs.version==version}\n"
        f"    FRS name: {version}\n"
        f"    FRS contents: {frs.version}\n"
    )
    return frs


############################### *v3.0.0+* ##############################

filename_pattern_v3 = filename_pattern_v2


def read_frs_v3(
        input_filename: str,
        folder: str | None = None,
) -> classes.FileReadSettings:
    """Reads FRS files from v1.0.0-1.1.0.
    
    **Currently does not check if multiple files relating to the same
    trial are present.** User should manually verify this is not true.
    """
    # Parse the input filename:
    m = filename_pattern_v3.search(input_filename)
    try:
        frs_id = m.group('r_id') # type: ignore
        version = classes.VersionNumber(m.group('version')) # type: ignore
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
    frs_filename = f'frs_{frs_id}_{version}.json'
    if folder is None:
        folder = r'outputs'
    filepath = rf'.\{folder}\JSON\file_read_settings\{frs_filename}'
    frs_json = json.load(open(filepath))
    frs = classes.FileReadSettings(
        classes.VersionNumber(frs_json['version']),
        frs_json['filename'],
        frs_json['repetition'],
        frs_json['recording_id'],
        frs_json['epoch_timing_ms'],
        frs_json['skip_superfast'],
        frs_json['spike_criteria'],
        frs_json['exclude_frequencies'],
        frs_json['exclude_amplitudes'],
        frs_json['enforce_max_failrate']
    )
    # Check file name and contents match:
    assert (
        frs.version == version
    ), (
        f"Mismatch between FRS name and contents: {frs_filename}\n"
        f"VERSION: {frs.version==version}\n"
        f"    FRS name: {version}\n"
        f"    FRS contents: {frs.version}\n"
    )
    return frs
