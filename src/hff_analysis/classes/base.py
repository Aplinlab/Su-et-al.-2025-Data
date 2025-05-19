"""#TODO Write module docstring.
"""

import re
import typing

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
            return (False, "incompatible")
        elif self.minor != other.minor:
            return (True, "compatible (minor version mismatch)")
        else:
            return (True, "compatible")
