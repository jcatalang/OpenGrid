"""OpenGrid: interactive structured / unstructured / hybrid mesh generator."""

from opengrid import _core

__version__ = "0.0.1"

__all__ = ["__version__", "core_version"]


def core_version() -> str:
    """Return the version string reported by the C++ core library."""
    return _core.version()
