"""Phase 11a demo: wall layers - grid lines normal to the wall over N layers.

  * 2D: the C-grid template on a NACA 0012 at a y+ ~ 1 first-cell height (1e-5 chords), with
    wall control only and with 30 wall layers (growth 1.2): deviation of the grid lines from
    the wall normal layer by layer, first-cell heights, quality.
  * 3D: the tapered (0.6) and swept (20 degrees) wing in a box, the same comparison over the
    wing's faces.
Writes 11a_wall_layers.png, 11a_wall_layers_report.txt and a project of the 2D case.

    python examples/11a_wall_layers.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib
import sys
import time
from pathlib import Path

import numpy as np

import opengrid

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
sys.path.insert(0, str(HERE.parent / "tests" / "python"))
cases = importlib.import_module("test_templates")
LAYERS = (0, 1, 4, 8, 16, 24, 29, 32)
WALL = 1e-5


def key(p):
    return tuple(np.round(np.asarray(p, dtype=float), 9))


def angle(a, b):
    c = np.sum(a * b, axis=-1) / (np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1))
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def deviations_2d(og, group) -> dict:
    """Per layer k: deviation (degrees) of the grid line from row k to k + 1 off every wall
    side from the normal of row k (its central-difference tangent, from the two neighbours)."""
    wall = set()
    for w in og.entity(group).params["walls"]:
        wall |= {key(p) for p in og.connector_nodes(og.store.get(w["connector"]).name).points}
    out = {}
    for d in og.names("domain"):
        g = og.domain_grid(d)
        for s in (g, g[::-1], np.swapaxes(g, 0, 1), np.swapaxes(g, 0, 1)[::-1]):
            if not all(key(p) in wall for p in s[0]):
                continue
            for k in LAYERS:
                if k + 1 < s.shape[0]:
                    t = s[k, 2:] - s[k, :-2]
                    v = s[k + 1, 1:-1] - s[k, 1:-1]
                    out.setdefault(k, []).append(np.abs(90.0 - angle(t, v)))
    return {k: np.concatenate(v) for k, v in out.items()}


def deviations_3d(og, group) -> dict:
    wall = set()
    for w in og.entity(group).params["walls"]:
        wall |= {key(p) for p in og.domain_grid(og.store.get(w["domain"]).name).reshape(-1, 3)}
    out = {}
    for b in og.names("block"):
        g = og.block_grid(b)
        for axis in range(3):
            m = np.moveaxis(g, axis, 0)
            for s in (m, m[::-1]):
                if not all(key(p) in wall for p in s[0].reshape(-1, 3)):
                    continue
                for k in LAYERS:
                    if k + 1 < s.shape[0]:
                        row = s[k]
                        n = np.cross(row[1:-1, 2:] - row[1:-1, :-2], row[2:, 1:-1] - row[:-2, 1:-1])
                        a = angle(n, s[k + 1, 1:-1, 1:-1] - row[1:-1, 1:-1])
                        out.setdefault(k, []).append(np.minimum(a, 180.0 - a).ravel())
    return {k: np.concatenate(v) for k, v in out.items()}


def table(report, label, plain, layered):
    report.append(f"  {label}: deviation from the wall normal (mean / max degrees)")
    report.append("    layer   wall control only   30 wall layers")
    for k in LAYERS:
        if k in plain and k in layered:
            a, b = plain[k], layered[k]
            plain_text = f"{a.mean():6.2f} / {a.max():6.2f}"
            layers_text = f"{b.mean():6.2f} / {b.max():6.2f}"
            report.append(f"    {k:5d}   {plain_text}     {layers_text}")


def cgrid(layers):
    og = opengrid.Session()
    options = {"radial_points": (61, 31)}
    if layers:
        options.update(layers=layers, growth=1.2)
    og.create_template("cgrid", body=cases.airfoil(og), wall=WALL, name="c", **options)
    t = time.perf_counter()
    result = og.smooth_group("c-group")
    return og, result, time.perf_counter() - t


def wing(layers):
    og = opengrid.Session()
    options = {"radial_points": (61, 31)}
    if layers:
        options.update(layers=layers, growth=1.2)
    og.create_template(
        "wing_in_box",
        body=cases.airfoil(og),
        wall=WALL,
        length=2.0,
        span_points=13,
        taper=0.6,
        sweep=20.0,
        name="w",
        **options,
    )
    t = time.perf_counter()
    result = og.smooth_group("w-group")
    return og, result, time.perf_counter() - t


def line(og, group, result, elapsed):
    state = f"converged in {result.iterations} iterations"
    if not result.converged:
        state = f"NOT converged after {result.iterations} iterations"
    return (
        f"{og.group_status(group).split(', interface')[0]}; {state}, {elapsed:.1f} s; "
        f"min scaled Jacobian {cases.worst(og):.3f}"
    )


def render(plain, layered, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(20, 12), constrained_layout=True)
    views = {
        "leading edge": (-0.01, 0.03, -0.02, 0.02),
        "trailing edge": (0.97, 1.03, -0.03, 0.03),
        "upper surface": (0.45, 0.55, 0.04, 0.08),
    }
    for row, (label, og) in enumerate((("wall control only", plain), ("30 wall layers", layered))):
        for ax, (where, (x0, x1, y0, y1)) in zip(axes[row], views.items(), strict=True):
            for d in og.names("domain"):
                g = og.domain_grid(d)
                ax.plot(g[..., 0], g[..., 1], color="#305080", lw=0.3)
                ax.plot(g[..., 0].T, g[..., 1].T, color="#305080", lw=0.3)
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_aspect("equal")
            ax.set_title(f"{label}: {where}")
    fig.savefig(path, dpi=80)
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
    window.run_smoother("c-group")
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
    project = OUTPUT / "11a_cgrid_layers.ogp"
    if args.gui:
        return run_gui(project, args.screenshot)
    report = [f"wall layers at a first-cell height of {WALL:g} chords (y+ ~ 1), growth 1.2", ""]
    plain, plain_result, plain_time = cgrid(0)
    layered, layered_result, layered_time = cgrid(30)
    report.append("2D C-grid (NACA 0012), 61 + 31 radial points:")
    report.append(f"  wall control only: {line(plain, 'c-group', plain_result, plain_time)}")
    report.append(f"  30 wall layers:    {line(layered, 'c-group', layered_result, layered_time)}")
    table(
        report,
        "body and wake cut",
        deviations_2d(plain, "c-group"),
        deviations_2d(layered, "c-group"),
    )
    height, _ = cases.first_layer(layered.domain_grid("c-d1"))
    report.append(
        f"  first cell with layers {height.min():.4g} .. {height.max():.4g} (requested {WALL:g})"
    )
    layered.save(project)
    layered.save_journal(OUTPUT / "11a_cgrid_layers_journal.py")
    report.append("")
    w_plain, w_plain_result, w_plain_time = wing(0)
    w_layered, w_layered_result, w_layered_time = wing(30)
    report.append("3D wing in a box (taper 0.6, sweep 20 degrees), 13 span points:")
    report.append(f"  wall control only: {line(w_plain, 'w-group', w_plain_result, w_plain_time)}")
    report.append(
        f"  30 wall layers:    {line(w_layered, 'w-group', w_layered_result, w_layered_time)}"
    )
    table(
        report, "wing faces", deviations_3d(w_plain, "w-group"), deviations_3d(w_layered, "w-group")
    )
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "11a_wall_layers_report.txt").write_text(text)
    render(plain, layered, OUTPUT / "11a_wall_layers.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
