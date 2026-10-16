"""Phase 10e demo: topology templates and automatic 2D topology (cross fields).

  * Templates (opengrid.model.templates): each on its geometry, with no vertex placed by hand:
    2D O-grid (ellipse), butterfly (disk), C-grid (NACA 0012), H-grid (channel with a bump),
    body in channel (cylinder); 3D pipe, cylinder in channel, bump channel, wing in box (C-H,
    extruded; and tapered + swept, smoothed in 3D) and sphere in box. Reported: blocks, time,
    min scaled Jacobian (target >= 0.4), and the C-grid's wall control down to y+ ~ 1.
  * Automatic topology (opengrid.model.crossfield): the cross field of an airfoil in a box, a
    two-element airfoil (main element and flap), a plate with three holes and the 9a
    bracket's section; singularities and the index theorem, the separatrices, and the sketch of
    four-sided blocks they cut.
Writes 10e_templates_auto.png, 10e_templates_auto_report.txt and a project of the plate.

    python examples/10e_templates_auto.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid.validation import ValidationError

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
sys.path.insert(0, str(HERE.parent / "tests" / "python"))
cases = importlib.import_module("test_templates")


def worst(og) -> float:
    return cases.worst(og)


def template_runs(report: list[str]) -> dict:
    """Each template on its geometry; returns the 2D sessions for the picture."""
    runs = {
        "ogrid (ellipse)": ("ogrid", lambda og: {"body": cases.ellipse(og), "wall": 0.002}),
        "butterfly (disk)": (
            "butterfly",
            lambda og: {"boundary": [og.create_circle((0, 0), 1.0, name="c")]},
        ),
        "cgrid (NACA 0012)": ("cgrid", lambda og: {"body": cases.airfoil(og), "wall": 0.001}),
        "hgrid (bump)": (
            "hgrid",
            lambda og: {**cases.channel(og, bump=True), "blocks": (3, 1), "wall": 0.005},
        ),
        "body_in_channel (cylinder)": (
            "body_in_channel",
            lambda og: {
                "body": [og.create_circle((1.0, 0.5), 0.2, name="cyl")],
                **cases.channel(og),
                "wall": 0.005,
            },
        ),
        "pipe": (
            "pipe",
            lambda og: {
                "boundary": [og.create_circle((0, 0), 0.5, name="c")],
                "length": 2.0,
                "span_points": 11,
            },
        ),
        "cylinder_in_channel": (
            "cylinder_in_channel",
            lambda og: {
                "body": [og.create_circle((1.0, 0.5), 0.2, name="cyl")],
                **cases.channel(og),
                "wall": 0.005,
                "length": 1.0,
                "span_points": 7,
            },
        ),
        "bump_channel": (
            "bump_channel",
            lambda og: {
                **cases.channel(og, bump=True),
                "blocks": (3, 1),
                "wall": 0.005,
                "length": 1.0,
            },
        ),
        "wing_in_box (C-H)": (
            "wing_in_box",
            lambda og: {
                "body": cases.airfoil(og),
                "wall": 0.001,
                "length": 1.0,
                "span_points": 9,
            },
        ),
        "wing_in_box tapered + swept": (
            "wing_in_box",
            lambda og: {
                "body": cases.airfoil(og),
                "wall": 0.001,
                "length": 2.0,
                "span_points": 13,
                "taper": 0.6,
                "sweep": 20.0,
            },
        ),
        "sphere_in_box": ("sphere_in_box", lambda og: {"radius": 0.5, "wall": 0.01}),
    }
    report.append("templates (min scaled Jacobian after the group smoothing; target >= 0.4):")
    sessions = {}
    for label, (kind, setup) in runs.items():
        og = opengrid.Session()
        options = setup(og)
        t = time.perf_counter()
        sketch = og.create_template(kind, **options)
        og.smooth_group(og.names("smoothing_group")[0])
        elapsed = time.perf_counter() - t
        quality = worst(og)
        blocks = len(og.entity(sketch).params["blocks"])
        verdict = "met"
        if quality < 0.4:
            verdict = "MISSED"
        report.append(
            f"  {label:28s} {blocks:3d} blocks  {elapsed:5.1f} s  min SJ {quality:.3f} ({verdict})"
        )
        if kind in ("ogrid", "butterfly", "cgrid", "hgrid", "body_in_channel"):
            sessions[label] = og
    report.append(
        "  (wing_in_box: the smoothed 2D C-grid extruded; tapered + swept: a 3D sketch smoothed"
        " in 3D, the wing on ruled surfaces, the end planes symmetry planes)"
    )
    wall_control(report)
    return sessions


def wall_control(report: list[str]) -> None:
    """The C-grid at a y+ ~ 1 wall spacing: first-cell height and orthogonality on the body
    and on the wake cut (an internal wall with twice the spacing)."""
    report.append(
        "  C-grid wall control (first cell; deviation of the first grid line from normal):"
    )
    for wall in (1e-3, 1e-5, 1e-6):
        og = opengrid.Session()
        og.create_template("cgrid", body=cases.airfoil(og), wall=wall, name="c")
        result = og.smooth_group("c-group")
        parts = []
        for label, domain in (("body", "c-d1"), ("wake cut", "c-d0")):
            height, deviation = cases.first_layer(og.domain_grid(domain))
            parts.append(
                f"{label} {height.min():.3g}..{height.max():.3g}, max {deviation.max():.2f} deg"
            )
        state = f"converged in {result.iterations} iterations"
        if not result.converged:
            state = f"NOT converged after {result.iterations} iterations"
        report.append(f"    wall {wall:g}: {state}; " + "; ".join(parts))


def auto_runs(report: list[str]) -> dict:
    report += ["", "automatic topology (cross field):"]
    out = {}
    for label, build in cases.CASES.items():
        og = opengrid.Session()
        curves = build(og)
        t = time.perf_counter()
        layout = og.cross_field(curves)
        field = layout.field
        indices = sorted(s.index for s in field.singularities)
        line = (
            f"  {label:18s} {len(field.singularities)} singularities (indices "
            f"{', '.join(f'{i:+.2f}' for i in indices) or 'none'}), corners "
            f"{sum((2 - c.quarters) / 4 for c in field.corners):+.2f}: sum "
            f"{field.index_sum():+.2f} = Euler characteristic {field.euler():+d}; "
            f"{len(layout.separatrices)} separatrices"
        )
        try:
            sketch = og.auto_topology(curves)
            og.smooth_group(og.names("smoothing_group")[0])
            blocks = og.entity(sketch).params["blocks"]
            line += (
                f"; {len(blocks)} four-sided blocks ({time.perf_counter() - t:.1f} s), "
                f"min scaled Jacobian {worst(og):.3f}"
            )
            out[label] = (og, layout)
        except ValidationError as error:
            line += f"; FAILED: {error.issues[0].message}"
            out[label] = (og, layout)
        report.append(line)
    return out


def render(templates: dict, autos: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 3, figsize=(20, 15), constrained_layout=True)
    panels = list(templates.items())
    views = {
        "ogrid (ellipse)": (-1.5, 1.5, -1.2, 1.2),
        "butterfly (disk)": (-1.1, 1.1, -1.1, 1.1),
        "cgrid (NACA 0012)": (-0.4, 1.8, -0.8, 0.8),
        "hgrid (bump)": (-0.1, 3.1, -0.1, 1.1),
        "body_in_channel (cylinder)": (-0.1, 3.1, -0.1, 1.1),
    }
    for ax, (label, og) in zip(axes.ravel()[:5], panels, strict=True):
        for d in og.names("domain"):
            g = og.domain_grid(d)
            ax.plot(g[..., 0], g[..., 1], color="#305080", lw=0.3)
            ax.plot(g[..., 0].T, g[..., 1].T, color="#305080", lw=0.3)
        lo_x, hi_x, lo_y, hi_y = views[label]
        ax.set_xlim(lo_x, hi_x)
        ax.set_ylim(lo_y, hi_y)
        ax.set_title(f"template {label}")
    for ax, (label, (og, layout)) in zip(axes.ravel()[5:], autos.items(), strict=False):
        for d in og.names("domain"):
            g = og.domain_grid(d)
            ax.plot(g[..., 0], g[..., 1], color="#b0b8c8", lw=0.25)
            ax.plot(g[..., 0].T, g[..., 1].T, color="#b0b8c8", lw=0.25)
        for name in og.names():
            entity = og.entity(name)
            if entity.kind in ("line", "arc", "spline") and not name.startswith("auto"):
                p = np.asarray(og.geometry(name).tessellate_points(1e-4))
                ax.plot(p[:, 0], p[:, 1], color="#202020", lw=1.2)
        for line in layout.separatrices:
            ax.plot(line[:, 0], line[:, 1], color="#d01c8b", lw=1.2)
        for s in layout.field.singularities:
            color = "#1f78b4"
            if s.index > 0:
                color = "#33a02c"
            ax.plot(s.point[0], s.point[1], "o", color=color, ms=6)
        title = f"automatic: {label}"
        if not og.names("domain"):
            title += " (FAILED: field and separatrices only)"
        ax.set_title(title)
    for ax in axes.ravel():
        ax.set_aspect("equal")
    fig.savefig(path, dpi=90)
    plt.close(fig)


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
    window.quality_action.setChecked(True)
    window.run_smoother(window.session.names("smoothing_group")[0])
    if screenshot is not None:

        def save_and_quit():
            window.smoother_dock.wait()
            app.processEvents()
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(4000, save_and_quit)
    return app.exec()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    project = OUTPUT / "10e_plate_auto.ogp"
    if args.gui:
        og = opengrid.Session()
        og.auto_topology(cases.plate(og), name="auto")
        og.save(project)
        return run_gui(project, args.screenshot)
    report: list[str] = []
    templates = template_runs(report)
    autos = auto_runs(report)
    plate, _layout = autos["plate_three_holes"]
    plate.save(project)
    plate.save_journal(OUTPUT / "10e_plate_auto_journal.py")
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "10e_templates_auto_report.txt").write_text(text)
    render(templates, autos, OUTPUT / "10e_templates_auto.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
