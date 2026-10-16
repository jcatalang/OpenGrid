"""Phase 11d demo: GPU compute (OpenCL) - device detection and the experimental GPU sweeps.

  * The device found and its self-test (opengrid.gpu_info()).
  * The 1M-cell ball (10b) smoothed as a group on the CPU and with the line sweeps on the GPU
    (opengrid.set_gpu("on")): time, iterations, final residual, difference between the grids
    in units of the local cell size.
The kernel-level comparison (GPU 7-9x faster on structured sweeps) and the precision study are
in tools/gpu. Writes 11d_gpu_report.txt.

    python examples/11d_gpu.py
"""

import importlib
import sys
import time
from pathlib import Path

import numpy as np

import opengrid

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
sys.path.insert(0, str(HERE))


def local_cell(grid: np.ndarray) -> np.ndarray:
    """Smallest edge length at each node of a block grid (nk, nj, ni, 3)."""
    h = np.full(grid.shape[:3], np.inf)
    for axis in range(3):
        e = np.linalg.norm(np.diff(grid, axis=axis), axis=-1)
        before = [(0, 0)] * 3
        before[axis] = (1, 0)
        after = [(0, 0)] * 3
        after[axis] = (0, 1)
        h = np.minimum(h, np.pad(e, before, constant_values=np.inf))
        h = np.minimum(h, np.pad(e, after, constant_values=np.inf))
    return h


def main() -> int:
    OUTPUT.mkdir(exist_ok=True)
    info = opengrid.gpu_info()
    report = ["GPU compute (OpenCL)"]
    if not info["available"]:
        report.append(f"  no GPU: {info['reason']}")
        text = "\n".join(report) + "\n"
        print(text, end="")
        (OUTPUT / "11d_gpu_report.txt").write_text(text)
        return 0
    report.append(
        f"  {info['device']} ({info['platform']}, {info['version'].strip()}, "
        f"{info['compute_units']} compute units, float64: {info['float64']}); self-test passed"
    )
    ball = importlib.import_module("10b_multiblock_3d")
    og, names = ball.ball_case(n=65, radial=33)
    group = og.create_smoothing_group(blocks=names)
    runs = {}
    for mode in ("off", "on"):
        opengrid.set_gpu(mode)
        og.groups._results.clear()
        t = time.perf_counter()
        result = og.smooth_group(group)
        runs[mode] = (result, time.perf_counter() - t)
    opengrid.set_gpu("off")
    report += ["", f"1M-cell ball (10b), {opengrid.threads()} CPU threads:"]
    for mode, label in (("off", "CPU"), ("on", "GPU line sweeps")):
        result, elapsed = runs[mode]
        res = np.asarray(result.residuals)
        report.append(
            f"  {label:16s} {elapsed:5.1f} s, {result.iterations} iterations "
            f"({1e3 * elapsed / result.iterations:.0f} ms each), final relative residual "
            f"{res[-1] / res[0]:.1e}"
        )
    cpu, gpu = runs["off"][0], runs["on"][0]
    worst = max(
        float(
            (np.linalg.norm(gpu.grids[b] - cpu.grids[b], axis=-1) / local_cell(cpu.grids[b])).max()
        )
        for b in cpu.grids
    )
    report.append(
        f"  GPU grid vs CPU grid: at most {worst:.1e} of the local cell size (the GPU stops at "
        f"a relative residual of {opengrid._core.gpu_tolerance:g}, the CPU at 1e-6)"
    )
    report.append(
        f"  speed-up {runs['off'][1] / runs['on'][1]:.2f}x - the GPU path is experimental and "
        f"off by default (Anderson mixing and sliding projections stay on the CPU)"
    )
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "11d_gpu_report.txt").write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
