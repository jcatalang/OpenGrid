"""Phase 3 demo: NACA 0012 O-grid and C-grid with prescribed wall spacing.

Both grids use elliptic smoothing (Thomas-Middlecoff control) with wall control on the
airfoil (j = 0): first-cell height WALL = 1e-4 chord and orthogonality. Exports VTKHDF and
Plot3D, and prints quality statistics, achieved wall spacing and orthogonality.

    python examples/03_airfoil_grids.py [--gui [--screenshot out.png]]
"""

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid.project import to_json

OUTPUT = Path(__file__).resolve().parent / "output"
WALL = 1e-4
N_SURFACE = 101


def naca0012(x):
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def airfoil(og, prefix):
    xs = [0.5 * (1.0 - math.cos(math.pi * i / 100)) for i in range(101)]
    og.create_spline(
        [(x, naca0012(x)) for x in xs], "clamped", start_tangent=(0, 1), name=f"{prefix}upper"
    )
    og.create_spline(
        [(x, -naca0012(x)) for x in xs], "clamped", start_tangent=(0, -1), name=f"{prefix}lower"
    )
    surface = {"type": "tanh", "begin_spacing": 2e-3, "end_spacing": 1e-3}
    og.create_connector(f"{prefix}upper", N_SURFACE, surface, name=f"{prefix}c-upper")
    og.create_connector(f"{prefix}lower", N_SURFACE, surface, name=f"{prefix}c-lower")


def o_grid(og):
    """Inner loop: airfoil (clockwise from the TE); outer loop: circle R = 15; cut along y = 0."""
    airfoil(og, "o-")
    og.create_circle((0.5, 0.0), 15.0, name="o-farfield")  # starts at (15.5, 0)
    og.create_line((1.0, 0.0), (15.5, 0.0), name="o-cut")
    og.create_connector("o-farfield", 2 * N_SURFACE - 1, name="o-c-farfield")
    og.create_connector(
        "o-cut", 81, {"type": "tanh", "begin_spacing": WALL, "end_spacing": 1.0}, name="o-c-cut"
    )
    # Explicit edges: j = 0 runs clockwise around the airfoil from the TE (lower surface TE -> LE,
    # then upper surface LE -> TE); edge 2 is listed in loop order, i.e. counter-clockwise
    # along the farfield, which is the circle's own direction.
    edges = [
        [("o-c-lower", True), "o-c-upper"],
        ["o-c-cut"],
        ["o-c-farfield"],
        [("o-c-cut", True)],
    ]
    return og.create_domain(edges, generator=elliptic_generator("hermite_tfi"), name="o-grid")


def elliptic_generator(method):
    # The O-grid starts from Hermite TFI (orthogonal at the airfoil); for the C-grid the
    # wall-normal Hermite derivatives fold the grid near the trailing-edge cusp, so it starts
    # from plain TFI and relies on the elliptic wall control.
    return {
        "method": method,
        "elliptic": {
            "control": "thomas_middlecoff",
            "walls": [{"edge": 0, "spacing": WALL}],
            "max_iterations": 20000,
            "tolerance": 1e-6,
        },
    }


def c_grid(og):
    """j = 0: wake (outflow -> TE), lower and upper surface, wake (TE -> outflow)."""
    airfoil(og, "c-")
    x_out = 16.0
    radius = 15.0
    og.create_line((1.0, 0.0), (x_out, 0.0), name="c-wake")
    og.create_line((x_out, 0.0), (x_out, radius), name="c-outflow-top")
    og.create_line((x_out, -radius), (x_out, 0.0), name="c-outflow-bottom")
    og.create_line((1.0, radius), (x_out, radius), name="c-top")
    og.create_line((1.0, -radius), (x_out, -radius), name="c-bottom")
    og.create_arc((1.0, 0.0), radius, math.pi / 2, 3 * math.pi / 2, name="c-front")
    wake = {"type": "tanh", "begin_spacing": 1e-3, "end_spacing": 1.0}
    normal = {"type": "tanh", "begin_spacing": WALL, "end_spacing": 1.0}
    og.create_connector("c-wake", 61, wake, name="c-c-wake")
    og.create_connector("c-outflow-top", 81, normal, name="c-c-outflow-top")
    og.create_connector(
        "c-outflow-bottom",
        81,
        {"type": "tanh", "begin_spacing": 1.0, "end_spacing": WALL},
        name="c-c-outflow-bottom",
    )
    og.create_connector("c-top", 61, wake, name="c-c-top")
    og.create_connector("c-bottom", 61, wake, name="c-c-bottom")
    og.create_connector("c-front", 2 * N_SURFACE - 1, name="c-c-front")
    edges = [
        [("c-c-wake", True), ("c-c-lower", True), "c-c-upper", "c-c-wake"],
        ["c-c-outflow-top"],
        [("c-c-top", True), "c-c-front", "c-c-bottom"],
        ["c-c-outflow-bottom"],
    ]
    return og.create_domain(edges, generator=elliptic_generator("tfi"), name="c-grid")


def wall_report(og, name):
    """First-cell heights and deviation from orthogonality along j = 0 (interior wall nodes),
    plus the maximum deviation away from the trailing-edge corner (x < 0.95)."""
    grid = og.domain_grid(name)
    wall_points = grid[0]
    first = grid[1]
    tangent = np.zeros_like(wall_points)
    tangent[1:-1] = wall_points[2:] - wall_points[:-2]
    normal_vec = first - wall_points
    spacing = np.linalg.norm(normal_vec, axis=1)
    interior = slice(1, -1)
    cosine = np.abs(np.sum(tangent[interior] * normal_vec[interior], axis=1)) / (
        np.linalg.norm(tangent[interior], axis=1) * spacing[interior]
    )
    deviation = np.degrees(np.arcsin(np.clip(cosine, 0, 1)))
    away = wall_points[interior, 0] < 0.95
    return spacing[interior], deviation, float(np.max(deviation[away]))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)

    og = opengrid.Session()
    names = [o_grid(og), c_grid(og)]
    for name in names:
        ni, nj = og.domain_dimensions(name)
        start = time.perf_counter()
        result = og.smooth_domain(name)
        elapsed = time.perf_counter() - start
        q = og.domain_quality(name)
        spacing, deviation, deviation_away = wall_report(og, name)
        print(
            f"{name}: {ni} x {nj}, smoothing {result.iterations} iterations in {elapsed:.1f} s, "
            f"converged={result.converged}"
        )
        for warning in result.warnings:
            print(f"  WARNING {warning}")
        print(
            f"  min scaled Jacobian {q['scaled_jacobian'].min():.3f}, "
            f"max skewness {q['skewness'].max():.3f}, "
            f"max non-orthogonality {q['non_orthogonality'].max():.1f} deg"
        )
        print(
            f"  wall spacing min/max {spacing.min():.4e} / {spacing.max():.4e} "
            f"(requested {WALL:g}); "
            f"max deviation from orthogonal {deviation.max():.2f} deg "
            f"({deviation_away:.2f} deg away from the trailing edge)"
        )
        sj = q["scaled_jacobian"]
        j, i = np.unravel_index(np.argmin(sj), sj.shape)
        print(f"  worst cell (i={i}, j={j}) at x = {og.domain_grid(name)[j, i, 0]:.4f}")
    og.save(OUTPUT / "03_airfoil_grids.ogp")
    og.save_journal(OUTPUT / "03_airfoil_grids_journal.py")
    summary = og.export_vtkhdf(OUTPUT / "03_airfoil_grids.vtkhdf")
    og.export_plot3d(OUTPUT / "03_airfoil_grids.xyz")
    print(f"VTKHDF: {summary}")
    replayed = opengrid.replay((OUTPUT / "03_airfoil_grids_journal.py").read_text())
    assert to_json(replayed.store) == to_json(og.store)
    render(og, names, OUTPUT / "03_airfoil_grids.png")
    if args.gui:
        return run_gui(OUTPUT / "03_airfoil_grids.ogp", args.screenshot)
    return 0


def render(og, names, path):
    import pyvista as pv

    from opengrid.app.scene import quad_mesh

    views = [
        ("o-grid", (0.5, 0.0), 0.7),
        ("o-grid", (1.0, 0.0), 0.02),
        ("c-grid", (0.5, 0.0), 0.7),
        ("c-grid", (-0.005, 0.0), 0.02),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(1800, 1200), border=False)
    for k, (name, center, scale) in enumerate(views):
        plotter.subplot(k // 2, k % 2)
        plotter.set_background("#1f2329")
        mesh = quad_mesh(og.domain_grid(name))
        plotter.add_mesh(mesh, style="wireframe", color="#9fb3c8", line_width=1)
        plotter.add_text(f"{name}  ({center[0]:g}, {center[1]:g})", font_size=10, color="white")
        plotter.enable_parallel_projection()
        plotter.camera.focal_point = (center[0], center[1], 0.0)
        plotter.camera.position = (center[0], center[1], 10.0)
        plotter.camera.up = (0, 1, 0)
        plotter.camera.parallel_scale = scale
    plotter.screenshot(str(path))
    plotter.close()


def run_gui(project, screenshot):
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication

    from opengrid.app.main_window import MainWindow

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = MainWindow()
    window.open_project(str(project))
    window.show()
    for name in window.session.names("domain"):
        window.session.smooth_domain(name)
    window.scene.rebuild()
    window.set_selection([window.session.entity("c-grid").id])
    window.quality_dock.metric.setCurrentIndex(5)  # non-orthogonality
    window.quality_action.setChecked(True)
    camera = window.viewport.camera
    camera.focal_point = (0.5, 0.0, 0.0)
    camera.position = (0.5, 0.0, 10.0)
    camera.parallel_scale = 0.8
    if screenshot is not None:

        def save_and_quit():
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(2000, save_and_quit)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
