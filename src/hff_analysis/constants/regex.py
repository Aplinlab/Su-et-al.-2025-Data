#*REGULAR EXPRESSIONS
# For parsing filenames and recording comments.

VERSION_REGEX = r"(?:(?P<major>\d+)\.(?P<minor>\d+)\.(?P<patch>\d+))"

# Reads a `.adicht` filename with or without the extension and captures
# the following groups:
# - `name`: Filename without file extension.
# - `id`: Animal-specific identifier made up of three letters followed
#   by two digits.
# - `position`: Position from which the recording was made (digits
#   only).
# - `testcode`: Four-letter code indicating the test which was
#   performed.
# - `extension`: File extension, if the filename includes one.
ADICHT_FILENAME_REGEX = (
    r"(?:^(?P<name>(?P<a_id>\w{3}\d{2})_pos(?P<pos>[\d]+)_"
    r"(?P<testcode>\w{4})[^.]*)(?P<extension>\..+)?$)"
)

# The following expressions parse comments attached to different test
# trains, and capture one or more of the following groups:
# - `mechval`: The value (usually frequency in Hz) of the mechanical
#   stimulus.
# - `elecval`: The value (usually frequency in Hz) of the electrical
#   stimulus.
# - `testvar`: In case of sweeps, captures the variable which was varied
#   (i.e. 'frequency' or 'amplitude').
# - `stimtype`: In case of the long-duration train, the type of stimulus
#   used during the conditioning phase (in our data, this should always
#   be 'electrical').
SWEEP_REGEX = (
    r"(?:^(?P<testvar>\w+)_mech_(?P<mechval>\d+(.\d+)?)_elec_"
    r"(?P<elecval>\d+(.\d+)?)$)"
)
NINEONE_REGEX = r"(?:^mech_(?P<mechval>\d)_elec_(?P<elecval>\d)$)"
LONGDURATION_REGEX = r"(?:^long_(?P<stimtype>\w{4})$)"

SAVED_FILENAME_REGEX = (
    r"(?:(?P<basename>(?P<savetype>\w+)_(?P<r_id>(?P<u_id>(?P<a_id>\w{3}\d{2})-"
    r"(?P<pos>\d+))_\[(?P<rep>\d+)-(?P<rec>\d+)\]_(?P<testcode>\w{4})))_"
    r"(?P<version>v\d+\.\d+\.\d+)(?P<extension>\..+)?$)"
)
