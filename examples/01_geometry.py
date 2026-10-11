"""Phase 1 demo: NACA 0012 airfoil + circular farfield built entirely from a script.

Writes to examples/output/:
    01_geometry.ogp          native project
    01_geometry_journal.py   replayable journal of the session
    01_geometry.png          overview + leading-edge close-up (headless render)

    python examples/01_geometry.py                 # headless
    python examples/01_geometry.py --gui           # also open the result in the GUI
    python examples/01_geometry.py --gui --screenshot out.png
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np

import opengrid
from opengrid.project import to_json

OUTPUT = Path(__file__).resolve().parent / "output"
N_SURFACE = 101  # points per surface, cosine-clustered at LE and TE
FARFIELD_RADIUS = 20.0  # chords


def naca0012(x: float) -> float:
    """Half thickness of the NACA 0012 with closed trailing edge (-0.1036 x^4 term)."""
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def build(og: opengrid.Session) -> None:
    xs = [0.5 * (1.0 - math.cos(math.pi * i / (N_SURFACE - 1))) for i in range(N_SURFACE)]
    upper = [(x, naca0012(x)) for x in xs]
    lower = [(x, -naca0012(x)) for x in xs]
    # Both surfaces run LE -> TE; a vertical LE tangent keeps the two splines G1 at the nose.
    og.create_spline(upper, method="clamped", start_tangent=(0, 1), name="naca-upper")
    og.create_spline(lower, method="clamped", start_tangent=(0, -1), name="naca-lower")
    og.create_point((0.0, 0.0), name="leading-edge")
    og.create_point((1.0, 0.0), name="trailing-edge")
    og.create_circle((0.5, 0.0), FARFIELD_RADIUS, name="farfield")


def report(og: opengrid.Session) -> None:
    for name in og.names():
        entity = og.entity(name)
        if entity.kind == "point":
            print(f"  {name:14s} point   {entity.params['xyz']}")
        else:
            print(f"  {name:14s} {entity.kind:7s} length = {og.length(name):.12f}")
    perimeter = og.length("naca-upper") + og.length("naca-lower")
    print(f"  airfoil perimeter = {perimeter:.10f} chords")
    exact = 2.0 * math.pi * FARFIELD_RADIUS
    error = abs(og.length("farfield") - exact) / exact
    print(f"  farfield relative length error vs 2*pi*R = {error:.2e}")


def render_png(og: opengrid.Session, path: Path) -> None:
    import pyvista as pv

    plotter = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1600, 700), border=False)
    for column, (center, scale, title) in enumerate(
        [((0.5, 0.0), 22.0, "Overview"), ((0.08, 0.0), 0.12, "Leading edge")]
    ):
        plotter.subplot(0, column)
        plotter.set_background("#1f2329")
        for name in og.names():
            entity = og.entity(name)
            if entity.kind == "point":
                plotter.add_points(
                    np.array([entity.params["xyz"]]),
                    color="#f2f2f2",
                    point_size=12,
                    render_points_as_spheres=True,
                )
                continue
            points = og.geometry(name).tessellate_points(1e-5)
            lines = np.concatenate([[len(points)], np.arange(len(points))])
            plotter.add_mesh(pv.PolyData(points, lines=lines), color="#5ad1ff", line_width=3)
            if name.startswith("naca"):
                knots = np.asarray(og.entity(name).params["points"])
                plotter.add_points(knots, color="#ff9f43", point_size=5)
        plotter.add_text(title, font_size=12, color="white")
        plotter.view_xy()
        plotter.enable_parallel_projection()
        plotter.camera.focal_point = (center[0], center[1], 0.0)
        plotter.camera.position = (center[0], center[1], 10.0 * scale)
        plotter.camera.parallel_scale = scale
    plotter.screenshot(str(path))
    plotter.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true", help="open the result in the GUI")
    parser.add_argument(
        "--screenshot", type=Path, help="with --gui: save a window screenshot and quit"
    )
    args = parser.parse_args()

    OUTPUT.mkdir(exist_ok=True)
    og = opengrid.Session()
    build(og)
    report(og)

    project = OUTPUT / "01_geometry.ogp"
    journal = OUTPUT / "01_geometry_journal.py"
    og.save(project)
    og.save_journal(journal)
    # The journal must rebuild exactly the same project.
    replayed = opengrid.replay(journal.read_text())
    assert to_json(replayed.store) == to_json(og.store), "journal replay mismatch"
    render_png(og, OUTPUT / "01_geometry.png")
    print(f"wrote {project.name}, {journal.name}, 01_geometry.png to {OUTPUT}")

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
    upper = window.session.entity("naca-upper").id
    window.set_selection([upper])
    camera = window.viewport.camera
    camera.focal_point = (0.5, 0.0, 0.0)
    camera.position = (0.5, 0.0, 10.0)
    camera.parallel_scale = 0.45
    window.console.run_line("og.length('naca-upper') + og.length('naca-lower')")
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
