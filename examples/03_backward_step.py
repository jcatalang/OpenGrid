"""Phase 3 demo: backward-facing step, three structured blocks sharing connectors.

Channel height H = 1 upstream of the step at x = 0, step height h = 1 (expansion ratio 2):
    inlet  block: x in [-5, 0],  y in [1, 2]
    lower  block: x in [0, 20],  y in [0, 1]   (behind the step)
    upper  block: x in [0, 20],  y in [1, 2]
Walls are clustered (tanh, first cell 1e-3 H); the step corner and the shear layer behind it
are refined in x. Blocks are plain TFI (rectangles). Exports VTKHDF and Plot3D multi-block and
checks that the shared connectors produce identical (merged) points.

    python examples/03_backward_step.py [--gui [--screenshot out.png]]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

import opengrid
from opengrid.io import read_plot3d
from opengrid.project import to_json

OUTPUT = Path(__file__).resolve().parent / "output"
WALL = 1e-3


def build(og: opengrid.Session) -> list[str]:
    def line(a, b, name):
        og.create_line(a, b, name=name)
        return name

    # Points: inlet (-5, 1)-(0, 1)-(0, 2)-(-5, 2); step (0, 0)-(0, 1); outlet x = 20.
    names = {
        "in-bottom": line((-5, 1), (0, 1), "in-bottom"),
        "in-top": line((-5, 2), (0, 2), "in-top"),
        "inlet": line((-5, 1), (-5, 2), "inlet"),
        "interface": line((0, 1), (0, 2), "interface"),
        "step": line((0, 0), (0, 1), "step"),
        "lo-bottom": line((0, 0), (20, 0), "lo-bottom"),
        "mid": line((0, 1), (20, 1), "mid"),
        "up-top": line((0, 2), (20, 2), "up-top"),
        "out-lo": line((20, 0), (20, 1), "out-lo"),
        "out-up": line((20, 1), (20, 2), "out-up"),
    }
    del names
    wall_both = {"type": "tanh", "begin_spacing": WALL, "end_spacing": WALL}
    toward_step = {"type": "tanh", "begin_spacing": 0.2, "end_spacing": 0.01}
    behind_step = {"type": "tanh", "begin_spacing": 0.01, "end_spacing": 0.5}
    og.create_connector("in-bottom", 41, toward_step, name="c-in-bottom")
    og.create_connector("in-top", 41, toward_step, name="c-in-top")
    og.create_connector("inlet", 51, wall_both, name="c-inlet")
    og.create_connector("interface", 51, wall_both, name="c-interface")
    og.create_connector("step", 51, wall_both, name="c-step")
    og.create_connector("lo-bottom", 161, behind_step, name="c-lo-bottom")
    # The shear-layer line y = 1 behind the step: same x distribution as the walls.
    og.create_connector("mid", 161, behind_step, name="c-mid")
    og.create_connector("up-top", 161, behind_step, name="c-up-top")
    og.create_connector("out-lo", 51, wall_both, name="c-out-lo")
    og.create_connector("out-up", 51, wall_both, name="c-out-up")

    inlet = og.create_domain(
        [["c-in-bottom"], ["c-interface"], [("c-in-top", True)], [("c-inlet", True)]],
        name="inlet-block",
    )
    lower = og.create_domain(
        [["c-lo-bottom"], ["c-out-lo"], [("c-mid", True)], [("c-step", True)]], name="lower-block"
    )
    upper = og.create_domain(
        [["c-mid"], ["c-out-up"], [("c-up-top", True)], [("c-interface", True)]], name="upper-block"
    )
    return [inlet, lower, upper]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)

    og = opengrid.Session()
    blocks = build(og)
    total_area = 0.0
    for name in blocks:
        ni, nj = og.domain_dimensions(name)
        q = og.domain_quality(name)
        total_area += float(q["area"].sum())
        print(
            f"{name}: {ni} x {nj}, min scaled Jacobian {q['scaled_jacobian'].min():.6f}, "
            f"max aspect ratio {q['aspect_ratio'].max():.1f}, "
            f"max area ratio {q['area_ratio'].max():.3f}"
        )
    print(f"total area {total_area:.12f} (exact 45)")
    assert abs(total_area - 45.0) < 1e-10

    # Shared connectors: the inlet block's right column equals the upper block's left column.
    inlet_grid = og.domain_grid("inlet-block")
    upper_grid = og.domain_grid("upper-block")
    assert np.array_equal(inlet_grid[:, -1], upper_grid[:, 0])
    lower_grid = og.domain_grid("lower-block")
    assert np.array_equal(lower_grid[-1, :], upper_grid[0, :])

    og.save(OUTPUT / "03_backward_step.ogp")
    og.save_journal(OUTPUT / "03_backward_step_journal.py")
    summary = og.export_vtkhdf(OUTPUT / "03_backward_step.vtkhdf")
    og.export_plot3d(OUTPUT / "03_backward_step.xyz")
    grids = read_plot3d(OUTPUT / "03_backward_step.xyz")
    assert all(np.array_equal(g, og.domain_grid(n)) for g, n in zip(grids, blocks, strict=True))
    print(f"VTKHDF: {summary}")
    merged = summary["points_before_merge"] - summary["points"]
    print(f"merged {merged} shared points (expected {51 + 161 - 0})")
    assert to_json(
        opengrid.replay((OUTPUT / "03_backward_step_journal.py").read_text()).store
    ) == to_json(og.store)
    render(og, blocks, OUTPUT / "03_backward_step.png")
    if args.gui:
        return run_gui(OUTPUT / "03_backward_step.ogp", args.screenshot)
    return 0


def render(og, blocks, path):
    import pyvista as pv

    from opengrid.app.scene import quad_mesh

    plotter = pv.Plotter(off_screen=True, shape=(2, 1), window_size=(1800, 1000), border=False)
    for row, (center, scale) in enumerate([((7.5, 1.0), 4.0), ((0.0, 1.0), 0.6)]):
        plotter.subplot(row, 0)
        plotter.set_background("#1f2329")
        colors = ["#9fb3c8", "#ffb86b", "#7ee787"]
        for name, color in zip(blocks, colors, strict=True):
            plotter.add_mesh(
                quad_mesh(og.domain_grid(name)), style="wireframe", color=color, line_width=1
            )
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
