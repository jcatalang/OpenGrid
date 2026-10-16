"""Phase 10c demo: freehand curves, domains from a curve network, live dragging.

  * Strokes: a channel with a bump, drawn freehand (simulated pointer samples with jitter):
    three horizontal strokes (the bottom one over the bump) and four vertical ones crossing
    them. Each stroke is fitted by NURBS curves (Session.create_freehand).
  * Network: the strokes are cut where they cross and the six four-cornered regions become
    structured domains (Session.domains_from_curves) - conforming, one undo step.
  * Smoothing: the six domains as a smoothing group (10a), the bottom as a wall.
  * Live editing: a control point of the bottom curve is dragged upwards with the GUI's
    handle logic (offscreen window): the domains follow during the drag (TFI preview, timed),
    one Modify command on release, then the group is smoothed again.
Writes 10c_freehand_live.png and 10c_freehand_live_report.txt, a project and its journal.

    python examples/10c_freehand_live.py [--gui [--screenshot out.png]]
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

import opengrid

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"


def strokes(rng):
    """Pointer samples of the seven strokes (with hand jitter)."""
    out = {}
    x = np.linspace(-0.15, 3.15, 120)
    bump = 0.25 * np.exp(-(((x - 1.5) / 0.35) ** 2))
    for name, y in (("bottom", bump), ("middle", 0.6 + 0.5 * bump), ("top", 1.2 + 0 * x)):
        jitter = 0.003 * rng.standard_normal(len(x))
        out[name] = np.stack([x, y + jitter, 0 * x], axis=1)
    for k, x0 in enumerate((0.0, 1.0, 2.0, 3.0)):
        y = np.linspace(-0.15, 1.35, 70)
        jitter = 0.003 * rng.standard_normal(len(y))
        out[f"v{k}"] = np.stack([x0 + 0.04 * np.sin(2.5 * y + k) + jitter, y, 0 * y], axis=1)
    return out


def build(report: list[str]):
    og = opengrid.Session()
    rng = np.random.default_rng(10)
    raw = strokes(rng)
    curves = []
    for name, points in raw.items():
        made = og.create_freehand(points, tolerance=0.004, name=name)
        controls = sum(len(og.entity(n).params["control_points"]) for n in made)
        report.append(
            f"stroke {name}: {len(points)} samples -> {len(made)} NURBS, {controls} control points"
        )
        curves.extend(made)
    t = time.perf_counter()
    domains = og.domains_from_curves(curves, spacing=0.04)
    report.append(
        f"domains from curves: {len(domains)} structured domains in "
        f"{time.perf_counter() - t:.2f} s, "
        f"{len(og.topology().interfaces())} interfaces, warnings: {og.last_warnings or 'none'}"
    )
    tfi = {d: og.domain_grid(d).copy() for d in domains}
    walls = [
        c
        for c in og.names("connector")
        if og.store.get(og.entity(c).params["segments"][0]["curve"]).name == "bottom"
    ]
    group = og.create_smoothing_group(domains, walls=walls, name="channel")
    result = og.smooth_group(group)
    report.append(f"smoothing group: {og.group_status(group)}")
    worst = min(float(og.domain_quality(d)["scaled_jacobian"].min()) for d in domains)
    report.append(f"  min scaled Jacobian {worst:.3f}")
    smoothed = {d: g.copy() for d, g in ((d, og.domain_grid(d)) for d in domains)}
    del result
    return og, raw, domains, tfi, smoothed


def drag(og, domains, report: list[str]):
    from PySide6.QtWidgets import QApplication

    from opengrid.app.main_window import MainWindow

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = MainWindow(off_screen=True, session=og)
    window.set_selection([og.entity("bottom").id])
    handles = window.visible_handles()
    # The control point closest to the top of the bump.
    target = min(handles, key=lambda h: np.hypot(h.xyz[0] - 1.5, h.xyz[1] - 0.25))
    assert window.handle_press(target.xyz)
    times = []
    for k in range(1, 21):
        xyz = (target.xyz[0], target.xyz[1] + 0.01 * k, 0.0)
        t = time.perf_counter()
        window.handle_drag(xyz)
        times.append(time.perf_counter() - t)
    final = (target.xyz[0], target.xyz[1] + 0.2, 0.0)
    window.handle_release(final)
    window.smoother_dock.wait()
    app.processEvents()
    og.smooth_group("channel")
    report.append(
        f"drag of bottom control point {target.index}: {len(times)} previews, "
        f"{1000 * np.mean(times):.1f} ms each ({1 / np.mean(times):.0f} frames/s), "
        f"one command: {og.journal_text().splitlines()[-1][:60]}..."
    )
    report.append(f"after the drag: {og.group_status('channel')}")
    window.close()
    return {d: og.domain_grid(d).copy() for d in domains}


def render(raw, og, domains, tfi, smoothed, dragged, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(16, 7.5), constrained_layout=True)
    ax = axes[0, 0]
    for points in raw.values():
        ax.plot(points[:, 0], points[:, 1], ".", ms=2, color="#999999")
    for name in og.names("nurbs"):
        p = np.asarray(og.geometry(name).tessellate_points(1e-4))
        ax.plot(p[:, 0], p[:, 1], lw=1.5)
    ax.set_title("freehand strokes (grey samples) and fitted NURBS curves (bottom: after the drag)")

    def grids(ax, grid_map, title):
        for d in domains:
            g = grid_map[d]
            ax.plot(g[..., 0], g[..., 1], color="#305080", lw=0.4)
            ax.plot(g[..., 0].T, g[..., 1].T, color="#305080", lw=0.4)
        ax.set_title(title)

    grids(axes[0, 1], tfi, "domains from the curve network (TFI)")
    grids(axes[1, 0], smoothed, "smoothing group (bottom = wall)")
    grids(axes[1, 1], dragged, "after dragging a control point of the bump, smoothed again")
    for ax in axes.ravel():
        ax.set_aspect("equal")
        ax.set_xlim(-0.2, 3.2)
        ax.set_ylim(-0.2, 1.4)
    fig.savefig(path, dpi=110)
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
    window.set_selection([window.session.entity("bottom").id])
    window.quality_action.setChecked(True)
    window.run_smoother("channel")
    if screenshot is not None:

        def save_and_quit():
            window.smoother_dock.wait()
            window.set_selection([window.session.entity("bottom").id])
            app.processEvents()
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(3000, save_and_quit)
    return app.exec()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    report: list[str] = []
    og, raw, domains, tfi, smoothed = build(report)
    dragged = drag(og, domains, report)
    text = "\n".join(report) + "\n"
    print(text, end="")
    if not args.gui:
        (OUTPUT / "10c_freehand_live_report.txt").write_text(text)
        render(raw, og, domains, tfi, smoothed, dragged, OUTPUT / "10c_freehand_live.png")
    og.save(OUTPUT / "10c_freehand.ogp")
    og.save_journal(OUTPUT / "10c_freehand_journal.py")
    if args.gui:
        return run_gui(OUTPUT / "10c_freehand.ogp", args.screenshot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
