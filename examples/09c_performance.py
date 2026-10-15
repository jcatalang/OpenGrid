"""Phase 9c demo: the performance targets of SPEC.md, measured.

Targets (Release build, Apple M1):
  * 2D, 1M cells in under 10 s - structured (TFI, 1001 x 1001), unstructured (Delaunay and
    advancing front) and hybrid (boundary layers + triangles); each timed from the model to
    the exported file (mesh, quality metrics, VTKHDF export);
  * structured elliptic smoothing of a 500 x 200 block "converges interactively" - read here as
    converging to the default tolerance (1e-6 of the first residual) within 2 s, with live
    progress; an O-grid around a cylinder with wall spacing control, the hardest common case;
  * 3D, 5M tetrahedra through the backend (Netgen) in under 2 minutes.

Also measured: the same smoothing without Anderson acceleration and on one thread, and a 3D
block smoothing. Parallel loops are deterministic: the one-thread run must give the same grid
(compared by digest). Writes 09c_benchmark.txt / .json / .png.

    python examples/09c_performance.py [--quick] [--skip-3d] [--gui [--screenshot out.png]]
"""

import argparse
import contextlib
import hashlib
import json
import math
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid import _core

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
# Exported files of the timed cases (large; deleted at the end).
SCRATCH = Path(tempfile.mkdtemp(prefix="opengrid-09c-"))


class Case:
    """Timed stages of one benchmark case."""

    def __init__(self, name, target=None, target_label=""):
        self.name = name
        self.target = target
        self.target_label = target_label
        self.stages = []
        self.notes = []
        self.cells = 0

    def time(self, label, function, *args, **kwargs):
        start = time.perf_counter()
        result = function(*args, **kwargs)
        self.stages.append((label, time.perf_counter() - start))
        return result

    @property
    def total(self):
        return sum(seconds for _label, seconds in self.stages)

    def met(self):
        if self.target is None:
            return None
        return self.total <= self.target

    def as_dict(self):
        return {
            "name": self.name,
            "cells": self.cells,
            "stages": self.stages,
            "total": self.total,
            "target": self.target,
            "met": self.met(),
            "notes": self.notes,
        }


def digest(array) -> str:
    return hashlib.sha1(np.ascontiguousarray(array).tobytes()).hexdigest()[:16]


def export(case, og, stem):
    case.time("VTKHDF export", og.export_vtkhdf, str(SCRATCH / f"{stem}.vtkhdf"))


# ---------------------------------------------------------------------------- 2D cases
def structured_2d(scale):
    n = round(1000 * scale) + 1
    case = Case(f"2D structured TFI {n} x {n}", 10.0, "< 10 s")
    og = opengrid.Session()
    og.create_arc((0, 0, 0), 1.0, 0.0, math.pi / 2, name="inner")
    og.create_arc((0, 0, 0), 2.0, 0.0, math.pi / 2, name="outer")
    og.create_line((1, 0), (2, 0), name="bottom")
    og.create_line((0, 2), (0, 1), name="left")
    for curve in ("inner", "outer", "bottom", "left"):
        og.create_connector(curve, n, name=f"c-{curve}")
    og.create_domain([["c-bottom"], ["c-outer"], ["c-left"], [("c-inner", True)]], name="d")
    og.create_patch("walls", "wall", ["c-inner", "c-outer"])
    og.create_patch("sides", "symmetry", ["c-bottom", "c-left"])
    grid = case.time("grid (TFI)", og.domain_grid, "d")
    quality = case.time("quality", og.domain_quality, "d")
    export(case, og, "structured_2d")
    case.cells = (n - 1) ** 2
    case.notes.append(f"min scaled Jacobian {quality['scaled_jacobian'].min():.3f}")
    assert grid.shape == (n, n, 3)
    return case


def unstructured_2d(scale, algorithm):
    h = 1.52e-3 / math.sqrt(scale)
    if algorithm == "advancing_front":
        h = 1.36e-3 / math.sqrt(scale)
    n = round(1.0 / h) + 1
    case = Case(f"2D unstructured ({algorithm.replace('_', ' ')})", 10.0, "< 10 s")
    og = opengrid.Session()
    corners = [(0, 0), (1, 0), (1, 1), (0, 1)]
    for k in range(4):
        og.create_line(corners[k], corners[(k + 1) % 4], name=f"e{k}")
        og.create_connector(f"e{k}", n, name=f"c{k}")
    og.create_unstructured_domain(
        connectors=[f"c{k}" for k in range(4)],
        generator={"algorithm": algorithm, "max_size": h, "max_points": 5_000_000},
        name="u",
    )
    og.create_patch("walls", "wall", [f"c{k}" for k in range(4)])
    mesh = case.time("mesh", og.unstructured_mesh, "u")
    quality = case.time("quality", og.domain_quality, "u")
    export(case, og, f"unstructured_2d_{algorithm}")
    case.cells = len(mesh["triangles"]) + len(mesh["quads"])
    case.notes.append(f"min angle {quality['min_angle'].min():.1f} deg")
    return case


def hybrid_2d(scale):
    h = 5.4e-3 / math.sqrt(scale)
    case = Case("2D hybrid (boundary layers + triangles)", 10.0, "< 10 s")
    og = opengrid.Session()
    og.create_circle((0, 0), 0.5, name="body")
    corners = [(-2, -2), (2, -2), (2, 2), (-2, 2)]
    for k in range(4):
        og.create_line(corners[k], corners[(k + 1) % 4], name=f"e{k}")
        og.create_connector(f"e{k}", round(4.0 / h) + 1, name=f"c{k}")
    og.create_connector("body", round(2000 * math.sqrt(scale)), name="wall")
    og.create_unstructured_domain(
        connectors=["wall", *[f"c{k}" for k in range(4)]],
        generator={
            "max_size": h,
            "max_points": 5_000_000,
            "inflation": {
                "connectors": ["wall"],
                "first_height": 1e-4,
                "growth_rate": 1.15,
                "layers": 30,
            },
        },
        name="flow",
    )
    og.create_patch("p-wall", "wall", ["wall"])
    og.create_patch("p-far", "farfield", [f"c{k}" for k in range(4)])
    mesh = case.time("mesh", og.unstructured_mesh, "flow")
    case.time("quality", og.domain_quality, "flow")
    export(case, og, "hybrid_2d")
    case.cells = len(mesh["triangles"]) + len(mesh["quads"])
    case.notes.append(f"{len(mesh['quads'])} layer quads, {len(mesh['triangles'])} triangles")
    return case


def cylinder_ogrid(ni, nj, acceleration, max_iterations=20000):
    og = opengrid.Session()
    og.create_circle((0, 0), 0.5, name="cyl")
    og.create_circle((0, 0), 20.0, name="far")
    og.create_line((0.5, 0), (20, 0), name="cut")
    og.create_connector("cyl", ni, name="c-cyl")
    og.create_connector("far", ni, name="c-far")
    og.create_connector("cut", nj, {"type": "tanh", "begin_spacing": 1e-3}, name="c-cut")
    og.create_domain([[("c-cyl", True)], ["c-cut"], ["c-far"], [("c-cut", True)]], name="d")
    og.set_domain_generator(
        "d",
        elliptic={
            "walls": [{"edge": 0, "spacing": 1e-3}],
            "acceleration": acceleration,
            "max_iterations": max_iterations,
        },
    )
    return og


def elliptic_2d(scale):
    ni = round(500 * math.sqrt(scale))
    nj = round(200 * math.sqrt(scale))
    case = Case(f"2D elliptic smoothing {ni} x {nj} (wall control)", 2.0, "converged < 2 s")
    og = cylinder_ogrid(ni, nj, 5)
    result = case.time("smoothing", og.smooth_domain, "d")
    case.cells = (ni - 1) * (nj - 1)
    spacing = np.linalg.norm(result.grid[1] - result.grid[0], axis=1)
    case.notes.append(
        f"converged {result.converged} in {result.iterations} iterations; wall spacing "
        f"{spacing.min():.6g} .. {spacing.max():.6g} (requested 0.001)"
    )
    if not result.converged:
        case.target = 0.0  # not converged: the target is missed whatever the time
    reference = digest(result.grid)
    # The same on one thread: identical grid.
    threads = opengrid.threads()
    opengrid.set_threads(1)
    serial = cylinder_ogrid(ni, nj, 5)
    start = time.perf_counter()
    one = serial.smooth_domain("d")
    seconds = time.perf_counter() - start
    opengrid.set_threads(threads)
    same = "DIFFERENT GRID"
    if digest(one.grid) == reference:
        same = "identical grid"
    case.notes.append(f"1 thread: {seconds:.2f} s, {same}")
    # Plain relaxation (no acceleration) for comparison, capped at 4000 iterations.
    plain = cylinder_ogrid(ni, nj, 0, 4000)
    start = time.perf_counter()
    p = plain.smooth_domain("d")
    seconds = time.perf_counter() - start
    ratio = p.residuals[-1] / p.residuals[0]
    text = f"without acceleration: {p.iterations} iterations in {seconds:.2f} s"
    if p.converged:
        case.notes.append(text + ", converged")
    else:
        case.notes.append(text + f", not converged (residual {ratio:.1e} of the first)")
    return case


def elliptic_3d(scale):
    n = round(41 * scale ** (1 / 3))
    case = Case(f"3D elliptic block smoothing {2 * n - 1} x {n} x {n} (wall control)")
    og = opengrid.Session()
    # A box extruded from its base with clustering at the base; the smoother enforces the wall
    # spacing and orthogonality there.
    base = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
    names = []
    dims = [2 * n - 1, n, 2 * n - 1, n]
    for k in range(4):
        og.create_line(base[k], base[(k + 1) % 4], name=f"b{k}")
        names.append(og.create_connector(f"b{k}", dims[k], name=f"cb{k}"))
    og.create_domain([[names[0]], [names[1]], [names[2]], [names[3]]], name="base")
    og.create_block(
        extrusion={
            "domain": "base",
            "mode": "translate",
            "vector": (0.0, 0.0, 1.0),
            "dimension": n,
            "distribution": {"type": "tanh", "begin_spacing": 0.005, "end_spacing": 0.1},
        },
        smoothing={"control": "laplace", "walls": [{"face": "kmin", "spacing": 0.005}]},
        name="block",
    )
    case.time("smoothing", og.block_grid, "block")
    case.cells = (2 * n - 2) * (n - 1) ** 2
    case.notes.append(og.block_status("block"))
    return case


# ---------------------------------------------------------------------------- 3D case
def box_faces(og, n):
    """Six unstructured faces of the unit cube, n nodes per edge."""
    corners = {}
    for c in ((i, j, k) for k in (0, 1) for j in (0, 1) for i in (0, 1)):
        corners[c] = tuple(float(v) for v in c)
    edges = {}

    def connector(a, b):
        key = tuple(sorted((a, b)))
        if key not in edges:
            name = f"e{len(edges)}"
            og.create_line(corners[key[0]], corners[key[1]], name=name)
            edges[key] = og.create_connector(name, n, name=f"c-{name}")
        return edges[key]

    faces = []
    for axis in range(3):
        for side in (0, 1):
            a, b = [x for x in range(3) if x != axis]

            def vertex(i, j, axis=axis, side=side, a=a, b=b):
                c = [0, 0, 0]
                c[axis] = side
                c[a] = i
                c[b] = j
                return tuple(c)

            loop = [vertex(0, 0), vertex(1, 0), vertex(1, 1), vertex(0, 1)]
            names = [connector(loop[k], loop[(k + 1) % 4]) for k in range(4)]
            faces.append(og.create_unstructured_domain(connectors=names, name=f"f{axis}{side}"))
    return faces


def tetrahedra_3d(scale):
    n = round(111 * scale ** (1 / 3))
    case = Case("3D tetrahedra via Netgen", 120.0, "volume < 120 s")
    og = opengrid.Session()
    faces = box_faces(og, n)
    og.create_volume(faces, generator={"max_size": 1.0 / (n - 1)}, name="v")
    og.create_patch("walls", "wall", faces)
    case.time("surface meshes", lambda: [og.unstructured_mesh(f) for f in faces])
    mesh = case.time("volume (Netgen + optimisation)", og.volume_mesh, "v")
    # The target is on the volume generation; quality and export are reported alongside.
    volume_seconds = case.total
    quality = case.time("quality", og.volume_quality, "v")
    export(case, og, "tetrahedra_3d")
    case.cells = len(mesh["tets"])
    case.target_total = volume_seconds
    case.notes.append(
        f"volume generation {volume_seconds:.1f} s; min scaled Jacobian "
        f"{quality['scaled_jacobian'].min():.3f}; {og.volume_status('v')}"
    )
    return case


# ---------------------------------------------------------------------------- report
def machine() -> str:
    cpu = platform.processor()
    with contextlib.suppress(OSError):
        cpu = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True
        ).stdout.strip()
    return (
        f"{cpu}, {platform.system()} {platform.release()}, Python {platform.python_version()}, "
        f"OpenGrid {opengrid.__version__}, OpenMP {_core.openmp_enabled()}, "
        f"{opengrid.threads()} threads"
    )


def report(cases, path_txt, path_json):
    lines = ["OpenGrid performance (Phase 9c)", machine(), ""]
    for case in cases:
        total = case.total
        verdict = ""
        measured = getattr(case, "target_total", total)
        if case.target is not None:
            if measured <= case.target:
                verdict = f"  [target {case.target_label}: met]"
            else:
                verdict = f"  [target {case.target_label}: MISSED]"
        lines.append(f"{case.name}: {case.cells:,} cells, {total:.2f} s{verdict}")
        for label, seconds in case.stages:
            lines.append(f"    {label:32s} {seconds:8.2f} s")
        for note in case.notes:
            lines.append(f"    - {note}")
        lines.append("")
    text = "\n".join(lines)
    print(text)
    path_txt.write_text(text + "\n")
    path_json.write_text(
        json.dumps({"machine": machine(), "cases": [c.as_dict() for c in cases]}, indent=1)
    )


def chart(cases, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 0.7 * len(cases) + 1.5))
    names = [c.name for c in cases]
    y = np.arange(len(cases))[::-1]
    for row, case in zip(y, cases, strict=True):
        colour = "#7f8c8d"
        measured = getattr(case, "target_total", case.total)
        if case.target is not None and measured <= case.target:
            colour = "#2a9d8f"
        elif case.target is not None:
            colour = "#e63946"
        ax.barh(row, case.total, color=colour)
        stages = ", ".join(f"{label} {seconds:.2f}" for label, seconds in case.stages)
        right = case.total
        if case.target is not None:
            right = max(right, case.target)
        ax.text(right * 1.15, row, f"{case.total:.2f} s ({stages})", va="center", fontsize=7)
        if case.target is not None and case.target > 0:
            ax.plot([case.target, case.target], [row - 0.4, row + 0.4], color="k", lw=2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{n}\n{c.cells:,} cells" for n, c in zip(names, cases, strict=True)])
    ax.set_xscale("log")
    ax.set_xlim(0.5 * min(c.total for c in cases), 300 * max(c.total for c in cases))
    ax.set_xlabel("seconds (log scale); black marks: targets; green met, red missed")
    ax.set_title(machine(), fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def run_gui(screenshot):
    """The smoother dock on the 500 x 200 O-grid: live residuals, converged."""
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication

    from opengrid.app.main_window import MainWindow

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = MainWindow()
    og = cylinder_ogrid(500, 200, 5)
    path = OUTPUT / "09c_ogrid.ogp"
    og.save(str(path))
    window.open_project(str(path))
    window.show()
    window.set_selection([window.session.entity("d").id])
    window.run_smoother("d")

    def finish():
        if window.smoother_dock.running:
            QTimer.singleShot(200, finish)
            return
        window.quality_action.setChecked(True)
        window.viewport.camera_position = [(0.0, 0.0, 4.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
        window.viewport.render()
        if screenshot is not None:
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

    QTimer.singleShot(1000, finish)
    return app.exec()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--quick", action="store_true", help="1/20 of the sizes (smoke test)")
    parser.add_argument("--skip-3d", action="store_true", help="skip the 5M-tetrahedra case")
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    if args.gui:
        status = run_gui(args.screenshot)
        shutil.rmtree(SCRATCH, ignore_errors=True)
        return status
    scale = 1.0
    stem = "09c_benchmark"
    if args.quick:
        scale = 0.05
        stem = "09c_benchmark_quick"
    cases = [
        structured_2d(scale),
        unstructured_2d(scale, "delaunay"),
        unstructured_2d(scale, "advancing_front"),
        hybrid_2d(scale),
        elliptic_2d(scale),
        elliptic_3d(scale),
    ]
    if not args.skip_3d:
        cases.append(tetrahedra_3d(scale))
    report(cases, OUTPUT / f"{stem}.txt", OUTPUT / f"{stem}.json")
    chart(cases, OUTPUT / f"{stem}.png")
    shutil.rmtree(SCRATCH, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
