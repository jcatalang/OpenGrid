"""Phase 11b demo: CPU performance of the group smoothers.

  * The 10d wing stub (3D sketch on CAD faces, 159k nodes) smoothed as a group on 1, 4 and
    10 threads: time, iterations, and the grids' digest (identical for every thread count).
  * The C-grid template at a y+ ~ 1 wall with 30 wall layers (2D group), 1 and 10 threads.
Before 11b (static OpenMP schedule, heap allocation in every NURBS evaluation) the wing stub
took 39.3 s on 10 threads. Writes 11b_performance_report.txt.

    python examples/11b_performance.py
"""

import hashlib
import importlib
import sys
import time
from pathlib import Path

import numpy as np

import opengrid

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "tests" / "python"))
BEFORE_WING = 39.3  # s, 10 threads, measured at 7199fd9 (Phase 11a)


def digest(grids) -> str:
    h = hashlib.sha256()
    for key in sorted(grids):
        h.update(np.ascontiguousarray(grids[key]).tobytes())
    return h.hexdigest()[:16]


def timed(og, group, threads):
    opengrid.set_threads(threads)
    og.groups._results.clear()  # smooth again, not from the cache
    t = time.perf_counter()
    result = og.smooth_group(group)
    return time.perf_counter() - t, result


def main() -> int:
    OUTPUT.mkdir(exist_ok=True)
    default = opengrid.threads()
    report = [f"OpenGrid group smoothing, {default} threads available", ""]
    sketch = importlib.import_module("10d_sketch_topology")
    sketch.write_cad()
    og = sketch.wing_session()
    group = og.names("smoothing_group")[0]
    report.append("10d wing stub (3D sketch on the CAD model):")
    digests = set()
    times = {}
    for threads in (1, 4, default):
        elapsed, result = timed(og, group, threads)
        times[threads] = elapsed
        digests.add(digest(result.grids))
        report.append(
            f"  {threads:2d} threads: {elapsed:5.1f} s, {result.iterations} iterations, "
            f"grid digest {digest(result.grids)}"
        )
    report.append(
        f"  speed-up over 11a ({BEFORE_WING} s on {default} threads): "
        f"{BEFORE_WING / times[default]:.1f}x; threads 1 -> {default}: "
        f"{times[1] / times[default]:.1f}x; identical grids: {len(digests) == 1}"
    )
    cases = importlib.import_module("test_templates")
    og = opengrid.Session()
    og.create_template(
        "cgrid", body=cases.airfoil(og), wall=1e-5, layers=30, radial_points=(61, 31), name="c"
    )
    report += ["", "C-grid template, 1e-5 wall, 30 wall layers (2D group):"]
    digests = set()
    for threads in (1, default):
        elapsed, result = timed(og, "c-group", threads)
        digests.add(digest(result.grids))
        report.append(f"  {threads:2d} threads: {elapsed:5.2f} s, {result.iterations} iterations")
    report.append(f"  identical grids: {len(digests) == 1}")
    opengrid.set_threads(default)
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "11b_performance_report.txt").write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
