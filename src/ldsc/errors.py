"""Package-specific exception hierarchy for user and internal LDSC failures.

The CLI boundary uses these classes to decide whether a failure should be
displayed as a clean user-facing message or logged with a developer traceback.
Workflow modules should raise the most specific subclass that matches the
failure source when adding new user-facing validation.
"""


class LDSCError(Exception):
    """Base class for package-specific failures."""


class LDSCUserError(LDSCError):
    """Raised when user input, CLI arguments, or configuration are invalid."""


class LDSCUsageError(LDSCUserError):
    """Raised when a public command or API call uses an unsupported option mix."""


class LDSCConfigError(LDSCUserError):
    """Raised when configuration values or provenance are incompatible."""


class LDSCInputError(LDSCUserError):
    """Raised when input paths, file schemas, or file contents are invalid."""


class LDSCDependencyError(LDSCUserError, ImportError):
    """Raised when an optional dependency is required for the requested workflow."""


class LDSCInternalError(LDSCError):
    """Raised when the package reaches an unexpected internal state."""
