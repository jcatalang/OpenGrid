"""OpenGrid: interactive structured / unstructured / hybrid mesh generator."""

from opengrid import _core

__version__ = "0.0.1"


def core_version() -> str:
    """Return the version string reported by the C++ core library."""
    return _core.version()


def threads() -> int:
    """Number of threads of the parallel loops (1 when the core was built without OpenMP)."""
    return _core.thread_count()


def set_threads(n: int) -> None:
    """Set the number of threads (default: one per core, or OMP_NUM_THREADS). Meshes do not
    depend on it: parallel loops are deterministic."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError(f"set_threads: expected a positive integer, got {n!r}")
    _core.set_thread_count(n)


from opengrid.api import Session, replay  # noqa: E402  (needs _core loaded first)
from opengrid.validation import Issue, ValidationError  # noqa: E402

__all__ = [
    "Issue",
    "Session",
    "ValidationError",
    "__version__",
    "core_version",
    "replay",
    "set_threads",
    "threads",
]
