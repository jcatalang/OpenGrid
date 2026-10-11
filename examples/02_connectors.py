"""Phase 2 demo: airfoil connectors with leading- and trailing-edge clustering.

Geometry: NACA 0012 (upper/lower splines, LE -> TE), a wake line from the trailing edge to a
circular farfield (R = 20 chords). Connectors:
    c-upper, c-lower   tanh, 121 points, LE spacing 1e-3 c, TE spacing 5e-4 c
    c-wake             tanh, 81 points, end spacing 1 c; the TE spacing is *not* given: it is
                       propagated from the airfoil connectors at the shared node
    c-farfield         161 points, breakpoint at the upstream point (index 80) with tanh
                       spacings 1 c (wake junction) -> 0.5 c (upstream) -> 1 c

Writes to examples/output/: 02_connectors.ogp, 02_connectors_journal.py, 02_connectors.png
(nodes coloured by stretching ratio) and prints each connector's achieved spacings.

    python examples/02_connectors.py [--gui [--screenshot out.png]]
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np

import opengrid
from opengrid.project import to_json

OUTPUT = Path(__file__).resolve().parent / "output"
CHORD = 1.0
R_FARFIELD = 20.0


def naca0012(x: float) -> float:
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def build(og: opengrid.Session) -> dict[str, dict]:
    xs = [0.5 * (1.0 - math.cos(math.pi * i / 100)) for i in range(101)]
    og.create_spline(
        [(x, naca0012(x)) for x in xs], "clamped", start_tangent=(0, 1), name="naca-upper"
    )
    og.create_spline(
        [(x, -naca0012(x)) for x in xs], "clamped", start_tangent=(0, -1), name="naca-lower"
    )
    og.create_circle((0.5, 0.0), R_FARFIELD, name="farfield")  # starts at (20.5, 0)
    og.create_line((CHORD, 0.0), (0.5 + R_FARFIELD, 0.0), name="wake")

    airfoil = {"type": "tanh", "begin_spacing": 1e-3, "end_spacing": 5e-4}
    og.create_connector("naca-upper", 121, airfoil, name="c-upper")
    og.create_connector("naca-lower", 121, airfoil, name="c-lower")
    og.create_connector("wake", 81, {"type": "tanh", "end_spacing": 1.0}, name="c-wake")
    og.create_connector(
        "farfield",
        161,
        {"type": "tanh", "begin_spacing": 1.0, "end_spacing": 0.5},
        breakpoints=[
            {
                "position": 0.5,
                "index": 80,
                "distribution": {"type": "tanh", "begin_spacing": 0.5, "end_spacing": 1.0},
            }
        ],
        name="c-farfield",
    )
    # Requested spacings for the check below.
    return {
        "c-upper": {"begin": 1e-3, "end": 5e-4},
        "c-lower": {"begin": 1e-3, "end": 5e-4},
        "c-wake": {"begin": 5e-4, "end": 1.0},
        "c-farfield": {"begin": 1.0, "end": 1.0},
    }


def check(og: opengrid.Session, requested: dict[str, dict]) -> None:
    for name, spacing in requested.items():
        print(og.connector_report(name))
        nodes = og.connector_nodes(name)
        for end, value in (("begin", nodes.begin_spacing), ("end", nodes.end_spacing)):
            error = abs(value - spacing[end]) / spacing[end]
            print(f"  {end} spacing relative error {error:.1e}")
            assert error < 1e-8, f"{name} {end} spacing off by {error}"
    far = og.connector_nodes("c-farfield").arc_positions
    upstream = far[80] / far[-1]
    print(f"c-farfield breakpoint node 80 at {upstream:.12f} of the circumference (requested 0.5)")


def render_png(og: opengrid.Session, path: Path, threshold: float) -> None:
    import pyvista as pv

    views = [
        ((0.5, 0.0), 22.0, "Overview"),
        ((0.03, 0.0), 0.05, "Leading edge"),
        ((1.0, 0.0), 0.02, "Trailing edge"),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2100, 700), border=False)
    for column, (center, scale, title) in enumerate(views):
        plotter.subplot(0, column)
        plotter.set_background("#1f2329")
        for name in og.names("connector"):
            nodes = og.connector_nodes(name)
            n = len(nodes.points)
            cells = np.column_stack([np.full(n - 1, 2), np.arange(n - 1), np.arange(1, n)]).ravel()
            line = pv.PolyData(nodes.points, lines=cells)
            line.cell_data["stretching"] = nodes.interval_ratios()
            plotter.add_mesh(
                line,
                scalars="stretching",
                cmap="coolwarm",
                clim=(1.0, threshold),
                line_width=4,
                show_scalar_bar=column == 2,
                scalar_bar_args={"title": "stretching ratio", "color": "white"},
            )
            plotter.add_points(
                nodes.points, color="#e8e8e8", point_size=6, render_points_as_spheres=True
            )
        plotter.add_text(title, font_size=12, color="white")
        plotter.enable_parallel_projection()
        plotter.camera.focal_point = (center[0], center[1], 0.0)
        plotter.camera.position = (center[0], center[1], 10.0 * scale)
        plotter.camera.up = (0.0, 1.0, 0.0)
        plotter.camera.parallel_scale = scale
    plotter.screenshot(str(path))
    plotter.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()

    OUTPUT.mkdir(exist_ok=True)
    og = opengrid.Session()
    requested = build(og)
    check(og, requested)

    project = OUTPUT / "02_connectors.ogp"
    journal = OUTPUT / "02_connectors_journal.py"
    og.save(project)
    og.save_journal(journal)
    replayed = opengrid.replay(journal.read_text())
    assert to_json(replayed.store) == to_json(og.store), "journal replay mismatch"
    for name in og.names("connector"):
        assert np.array_equal(
            replayed.connector_nodes(name).points, og.connector_nodes(name).points
        )
    render_png(og, OUTPUT / "02_connectors.png", og.stretching_threshold)
    print(f"wrote {project.name}, {journal.name}, 02_connectors.png")
    if args.gui:
        return run_gui(project, args.screenshot)
    return 0


def run_gui(project: Path, screenshot: Path | None) -> int:
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication

    from opengrid.app.main_window import MainWindow

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = MainWindow()
    window.open_project(str(project))
    window.show()
    window.set_selection([window.session.entity("c-upper").id])
    camera = window.viewport.camera
    camera.focal_point = (0.12, 0.0, 0.0)
    camera.position = (0.12, 0.0, 10.0)
    camera.parallel_scale = 0.16
    if screenshot is not None:

        def save_and_quit():
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(1500, save_and_quit)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
