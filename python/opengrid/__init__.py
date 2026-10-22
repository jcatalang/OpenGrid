"""OpenGrid: interactive structured / unstructured / hybrid mesh generator."""

from opengrid import _core

__version__ = "0.0.2"


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


def gpu_info() -> dict:
    """GPU compute (OpenCL): {"built", "available", "platform", "device", "version",
    "compute_units", "float64", "reason"}. The first call detects the device, builds the
    kernels and runs a self-test of their float-float arithmetic."""
    raw = dict(_core.gpu_info())
    info = {key: raw[key] for key in ("platform", "device", "version", "reason")}
    info["built"] = raw["built"] == "1"
    info["available"] = raw["available"] == "1"
    info["float64"] = raw["float64"] == "1"
    info["compute_units"] = int(raw["compute_units"])
    return info


def gpu_mode() -> str:
    """'off' (default), 'auto' or 'on'."""
    return _core.gpu_mode()


def set_gpu(mode: str) -> None:
    """Experimental: run the line sweeps of the 3D group smoothing on the GPU - 'on' (any
    size), 'auto' (groups of at least ``opengrid._core.gpu_min_nodes`` nodes) or 'off' (the
    default). GPU runs use float32 arithmetic on float-float positions and stop at a relative
    residual of at most ``opengrid._core.gpu_tolerance`` (1e-5); their grids are not
    bit-identical to the CPU's. Measured on an M1 Max (Phase 11d): 1.2x per iteration on a
    1M-cell group, slower than the CPU to the same accuracy - hence off by default."""
    if mode not in ("off", "auto", "on"):
        raise ValueError(f"set_gpu: expected 'off', 'auto' or 'on', got {mode!r}")
    _core.set_gpu_mode(mode)


from opengrid.api import Session, replay  # noqa: E402  (needs _core loaded first)
from opengrid.validation import Issue, ValidationError  # noqa: E402

__all__ = [
    "Issue",
    "Session",
    "ValidationError",
    "__version__",
    "core_version",
    "gpu_info",
    "gpu_mode",
    "replay",
    "set_gpu",
    "set_threads",
    "threads",
]
