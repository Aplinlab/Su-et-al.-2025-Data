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

import json
import re

from . import classes
from . import constants
from . import spikedetect
from . import spikefilter
from . import utils


def update_outputs(
        save_outputs: str | list[str],
        input_filenames: list[str],
        input_folder: str | None = None
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
    * If an updater does not exist matching the versino of any filename
    in `input_filenames`, prints an error and skips that filename.
    """
    # Count number of updated files:
    updated_count = 0
    for filename in input_filenames:
        # Determine version of input file and pass to correct reader:
        file_version = classes.VersionNumber(filename)
        if file_version < classes.VersionNumber('2.0.0'):
            frs = read_frs_v1(filename, input_folder)
        elif file_version >= utils.current_version:
            # Skip loop if file is up-to-date (does not update count):
            print(f"Up to date: {filename}")
            continue
        else:
            # Skip loop if matching updater not implemented:
            print(f"No updater for {file_version}: {filename}")
            continue
        # Set variables to be passed to other functions:
        if save_outputs == 'all' or save_outputs == ['all']:
            save_outputs = constants.OUTPUT_TYPES
        try:
            save_clusters = 'clusters' in save_outputs
        except TypeError:
            save_clusters = False
        force_overwrite = False
        # Read and process LabChart data:
        recording = spikedetect.read_adicht(frs.filename, [frs.recording_segment])[0]
        (spikes, epochs) = spikedetect.spikes_info(
            recording,
            frs.repetition,
            frs.epoch_timing_ms,
            frs.threshold_uV
        )
        filtered_spikes = spikefilter.filter_spikes(
            spikes,
            frs.spike_criteria,
            frs.exclude_frequencies,
            frs.exclude_amplitudes
        )
        # Draw cluster plots:
        if save_clusters:
            spikefilter.plot_clusters(
                filtered_spikes,
                epochs,
                frs.repetition,
                frs.recording_segment,
                frs.exclude_frequencies,
                frs.exclude_amplitudes,
                save_clusters
            )
        # Save specified JSON outputs:
        if save_outputs:
            for (file_type, file_contents) in {
                'frs': {
                    'filename': frs.filename,
                    'epoch_timing_ms': frs.epoch_timing_ms,
                    'threshold_uV': frs.threshold_uV,
                    'spike_criteria': frs.spike_criteria,
                    'exclude_frequencies': frs.exclude_frequencies,
                    'exclude_amplitudes': frs.exclude_amplitudes
                },
                'epochs': {
                    'epochs': epochs,
                    'exclude_frequencies': frs.exclude_frequencies,
                    'exclude_amplitudes': frs.exclude_amplitudes
                },
                'spikes': {
                    'spikes': filtered_spikes
                }
            }.items():
                if file_type in save_outputs:
                    spikefilter.save_to_json(
                        frs.filename,
                        frs.repetition,
                        frs.recording_segment,
                        file_type,
                        force_overwrite,
                        **file_contents
                    )
        # Print completion message and update total updated files count:
        print(f"UPDATED from {file_version}: {frs.filename}")
        updated_count += 1
    # Print total updated files count:
    print(f"\nCOMPLETE: {updated_count} files updated.")
    return


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
        adicht_name = m.group('adicht_name') # type: ignore
        repetition = int(m.group('rep')) # type: ignore
        recording_segment = int(m.group('rec')) # type: ignore
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
