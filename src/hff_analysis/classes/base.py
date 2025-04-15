"""#TODO Write module docstring.
"""

import re

from hff_analysis import constants


# Regex pattern for parsing version numbers:
# It is compiled here once, rather than within a function, to reduce
# unnecessary computations.
version_regex_pattern = re.compile(constants.VERSION_REGEX)


class VersionNumber:
    def __init__(
            self,
            version_str: str
    ):
        m = version_regex_pattern.search(version_str)
        try:
            self.major = int(m.group('major'))
            self.minor = int(m.group('minor'))
            self.patch = int(m.group('patch'))
        except (AttributeError, ValueError) as e:
            # AttributeError is raised if m is None (i.e. if regex pattern
            # didn't match)
            raise ValueError(
                f"Unexpected error parsing version number ({version_str})."
            ) from e