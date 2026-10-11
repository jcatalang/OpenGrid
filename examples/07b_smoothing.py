"""Phase 7b demo: hyperbolic extrusion and 3D elliptic smoothing.

Case 1 - "cylinder_ogrid": O-grid around a cylinder (r = 0.5, span 1) grown by hyperbolic
         extrusion from the wall to r = 10.5: 129 x 21 x 61 points, first cell 1e-3 (2e-3 D),
         periodic around the cylinder (the seam is a block-internal interface), constant-z ends.
Case 2 - "bump_channel": channel over a Gaussian bump, extruded spanwise, then smoothed in 3D
         (Thomas-Middlecoff control, wall control on the bump: first spacing 0.005 and
         orthogonality); compared with the plain extrusion.

Exports CGNS 3D (cgnscheck when available), VTKHDF and Gmsh; writes 07b_smoothing.png and the
projects/journals.

    python examples/07b_smoothing.py [--gui [--screenshot out.png]]
"""

import argparse
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid import _core
from opengrid.project import to_json

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"


def cylinder_ogrid() -> opengrid.Session:
    og = opengrid.Session()
    og.create_circle((0, 0, 0), 0.5, name="circle0")
    og.create_circle((0, 0, 1), 0.5, name="circle1")
    og.create_line((0.5, 0, 0), (0.5, 0, 1), name="seam")
    og.create_connector("circle0", 129, name="c-circle0")
    og.create_connector("circle1", 129, name="c-circle1")
    og.create_connector("seam", 21, name="c-seam")
    og.create_domain(
        [["c-circle0"], ["c-seam"], [("c-circle1", True)], [("c-seam", True)]], name="wall"
    )
    og.create_extruded_block(
        "wall",
        "hyperbolic",
        distance=10.0,
        dimension=61,
        distribution={"type": "tanh", "begin_spacing": 1e-4},
        edges=("periodic", "periodic", "constant_z", "constant_z"),
        name="ogrid",
    )
    og.create_patch("cylinder", "wall", ["wall"])
    og.create_patch("farfield", "farfield", [{"block": "ogrid", "face": "kmax"}])
    og.create_patch(
        "ends",
        "symmetry",
        [{"block": "ogrid", "face": "jmin"}, {"block": "ogrid", "face": "jmax"}],
    )
    return og


def bump_channel(smoothed: bool) -> opengrid.Session:
    og = opengrid.Session()
    xs = np.linspace(0.0, 3.0, 61)
    bump = [(x, 0.0, 0.15 * math.exp(-8.0 * (x - 1.5) ** 2)) for x in xs]
    og.create_spline(bump, name="bump")
    og.create_line((3, 0, bump[-1][2]), (3, 0, 1), name="outlet")
    og.create_line((3, 0, 1), (0, 0, 1), name="top")
    og.create_line((0, 0, 1), (0, 0, bump[0][2]), name="inlet")
    og.create_connector("bump", 61, name="c-bump")
    og.create_connector("outlet", 31, {"type": "tanh", "begin_spacing": 0.005}, name="c-outlet")
    og.create_connector("top", 61, name="c-top")
    og.create_connector("inlet", 31, {"type": "tanh", "end_spacing": 0.005}, name="c-inlet")
    og.create_domain([["c-bump"], ["c-outlet"], ["c-top"], ["c-inlet"]], name="side")
    smoothing = None
    if smoothed:
        smoothing = {
            "control": "thomas_middlecoff",
            "walls": [{"face": "jmin", "spacing": 0.005, "decay": 0.3}],
            "max_iterations": 5000,
            "tolerance": 1e-9,
        }
    og.create_extruded_block(
        "side", "translate", vector=(0, 1, 0), dimension=11, name="channel", smoothing=smoothing
    )
    og.create_patch("bump-wall", "wall", [{"block": "channel", "face": "jmin"}])
    og.create_patch("upper", "symmetry", [{"block": "channel", "face": "jmax"}])
    og.create_patch("in", "inlet", [{"block": "channel", "face": "imin"}])
    og.create_patch("out", "outlet", [{"block": "channel", "face": "imax"}])
    og.create_patch("sides", "symmetry", ["side", {"block": "channel", "face": "kmax"}])
    return og


def wall_metrics(grid, spacing):
    worst_spacing = 0.0
    worst_angle = 0.0
    nk, _nj, ni, _ = grid.shape
    for k in range(1, nk - 1):
        for i in range(1, ni - 1):
            line = grid[k, 1, i] - grid[k, 0, i]
            n = np.cross(
                grid[k, 0, i + 1] - grid[k, 0, i - 1], grid[k + 1, 0, i] - grid[k - 1, 0, i]
            )
            worst_spacing = max(worst_spacing, abs(np.linalg.norm(line) / spacing - 1.0))
            c = abs(line @ n) / (np.linalg.norm(line) * np.linalg.norm(n))
            worst_angle = max(worst_angle, math.degrees(math.acos(min(1.0, c))))
    return worst_spacing, worst_angle


def export(og: opengrid.Session, label: str) -> None:
    result = og.boundary_report()
    assert result.ok, result.issues
    print(f"  interfaces: {', '.join(result.interfaces) or 'none'}")
    if _core.cgns_available():
        cgns = OUTPUT / f"07b_{label}.cgns"
        print(f"  CGNS {og.export_cgns(cgns)}")
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run([checker, str(cgns)], capture_output=True, text=True).stdout
            status = "clean"
            if "error" in out.lower() or "warning" in out.lower():
                status = "PROBLEMS:\n" + out
            print(f"  cgnscheck: {status}")
    print(f"  VTKHDF {og.export_vtkhdf(OUTPUT / f'07b_{label}.vtkhdf')}")
    print(f"  Gmsh {og.export_gmsh(OUTPUT / f'07b_{label}.msh')}")
    og.save(OUTPUT / f"07b_{label}.ogp")
    og.save_journal(OUTPUT / f"07b_{label}_journal.py")
    replayed = opengrid.replay((OUTPUT / f"07b_{label}_journal.py").read_text())
    assert to_json(replayed.store) == to_json(og.store)


def structured(grid):
    import pyvista as pv

    return pv.StructuredGrid(
        np.ascontiguousarray(grid[..., 0].transpose(2, 1, 0)),
        np.ascontiguousarray(grid[..., 1].transpose(2, 1, 0)),
        np.ascontiguousarray(grid[..., 2].transpose(2, 1, 0)),
    )


def render(ogrid, plain, smooth, path):
    import pyvista as pv

    plotter = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2400, 900), border=False)
    plotter.subplot(0, 0)
    mesh = structured(ogrid.block_grid("ogrid"))
    cut = mesh.clip(normal=(0, 0, 1), origin=(0, 0, 0.5), crinkle=True)
    plotter.add_mesh(cut, color="#7fd1b9", show_edges=True, edge_color="#30363d")
    plotter.add_text("cylinder O-grid (hyperbolic), cut z = 0.5", font_size=11, color="white")
    plotter.camera_position = "xy"
    plotter.camera.parallel_projection = True
    plotter.reset_camera()
    plotter.camera.zoom(4.0)
    for k, (og, title) in enumerate(
        ((plain, "bump channel: extrusion (TFI)"), (smooth, "3D elliptic + wall control"))
    ):
        plotter.subplot(0, k + 1)
        grid = og.block_grid("channel")
        mesh = structured(grid)
        q = og.block_quality("channel")["non_orthogonality"].reshape(-1)
        mesh.cell_data["non-orthogonality"] = q
        cut = mesh.clip(normal=(0, 1, 0), origin=(0, 0.5, 0), crinkle=True)
        plotter.add_mesh(
            cut,
            scalars="non-orthogonality",
            clim=(0, 30),
            cmap="viridis",
            show_edges=True,
            edge_color="#30363d",
        )
        plotter.add_text(title, font_size=11, color="white")
        plotter.camera.parallel_projection = True
        plotter.camera.focal_point = (1.5, 0.5, 0.5)
        plotter.camera.position = (1.5, -10.0, 0.5)
        plotter.camera.up = (0.0, 0.0, 1.0)
        plotter.camera.parallel_scale = 1.8
    for k in range(3):
        plotter.subplot(0, k)
        plotter.set_background("#1f2329")
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
    print("cylinder_ogrid:")
    start = time.perf_counter()
    ogrid = cylinder_ogrid()
    grid = ogrid.block_grid("ogrid")
    q = ogrid.block_quality("ogrid")
    r = np.hypot(grid[..., 0], grid[..., 1])
    print(
        f"  {ogrid.block_status('ogrid')} ({time.perf_counter() - start:.2f} s); first cell "
        f"{(r[1] - r[0]).mean():.3e}, outer radius {r[-1].min():.6f}..{r[-1].max():.6f}, "
        f"scaled Jacobian min {q['scaled_jacobian'].min():.4f}, "
        f"max non-orthogonality {q['non_orthogonality'].max():.2f} deg"
    )
    export(ogrid, "cylinder_ogrid")
    print("bump_channel:")
    plain = bump_channel(False)
    start = time.perf_counter()
    smooth = bump_channel(True)
    smooth.block_grid("channel")
    print(f"  smoothing: {smooth.block_status('channel')} ({time.perf_counter() - start:.1f} s)")
    for label, og in (("extrusion", plain), ("smoothed", smooth)):
        q = og.block_quality("channel")
        spacing, angle = wall_metrics(og.block_grid("channel"), 0.005)
        print(
            f"  {label:9s}: first-spacing error {spacing:.3f}, wall angle {angle:.2f} deg, "
            f"max non-orthogonality {q['non_orthogonality'].max():.2f} deg, scaled Jacobian min "
            f"{q['scaled_jacobian'].min():.3f}"
        )
    export(smooth, "bump_channel")
    render(ogrid, plain, smooth, OUTPUT / "07b_smoothing.png")
    if args.gui:
        return run_gui(OUTPUT / "07b_bump_channel.ogp", args.screenshot)
    return 0


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
    window.set_selection([window.session.entity("channel").id])
    window.quality_dock.metric.setCurrentIndex(5)  # non-orthogonality
    window.quality_action.setChecked(True)
    window.cut_action.setChecked(True)
    window.viewport.camera_position = "iso"
    window.viewport.reset_camera()
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
