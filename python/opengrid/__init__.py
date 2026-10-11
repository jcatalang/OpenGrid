"""OpenGrid: interactive structured / unstructured / hybrid mesh generator."""

from opengrid import _core

__version__ = "0.0.1"


def core_version() -> str:
    """Return the version string reported by the C++ core library."""
    return _core.version()


from opengrid.api import Session, replay  # noqa: E402  (needs _core loaded first)
from opengrid.validation import Issue, ValidationError  # noqa: E402

__all__ = ["Issue", "Session", "ValidationError", "__version__", "core_version", "replay"]
