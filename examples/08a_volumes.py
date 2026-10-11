"""Phase 8a demo: unstructured 3D volumes (tetrahedra via the Netgen backend).

A finite cylinder (radius 0.5, length 1) across a channel [-3, 8] x [-3, 3] x [-2, 2]:
  * the cylinder side is two unstructured domains on extruded (half-cylinder) surfaces, the end
    caps planar unstructured disks;
  * a refinement box [-1.5, 3] x [-1.5, 1.5] x [-1.5, 1.5] around it is the volume "near"
    (fine tetrahedra, max size 0.12), the rest of the channel the volume "far" (max size 0.6);
    the two volumes share the refinement box's faces, so the mesh is conforming across them;
  * the channel walls are structured faces (their quads are split into triangles: no block
    shares them).
Exports CGNS (unstructured zones, BCs, interface connectivity; cgnscheck when available),
VTKHDF and Gmsh 4.1; writes 08a_volumes.png and the project / journal.

    python examples/08a_volumes.py [--gui [--screenshot out.png]]
"""

import argparse
import itertools
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid import _core
from opengrid.model.volumes import total_volume
from opengrid.project import to_json

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"


def box(og, lo, hi, n, tag, kind):
    corners = {
        c: tuple(lo[k] + c[k] * (hi[k] - lo[k]) for k in range(3))
        for c in itertools.product((0, 1), repeat=3)
    }
    edges = {}

    def connector(c, d):
        if (c, d) not in edges:
            name = f"{tag}-e{len(edges)}"
            og.create_line(corners[c], corners[d], name=name)
            length = math.dist(corners[c], corners[d])
            edges[(c, d)] = og.create_connector(
                name, max(3, round(length * n) + 1), name=f"c-{name}"
            )
        return edges[(c, d)]

    faces = []
    for axis in range(3):
        for side in (0, 1):
            a, b = [x for x in range(3) if x != axis]

            def v(i, j, axis=axis, side=side, a=a, b=b):
                c = [0, 0, 0]
                c[axis] = side
                c[a] = i
                c[b] = j
                return tuple(c)

            names = [
                connector(v(0, 0), v(1, 0)),
                connector(v(1, 0), v(1, 1)),
                connector(v(0, 1), v(1, 1)),
                connector(v(0, 0), v(0, 1)),
            ]
            face = f"{tag}-{'xyz'[axis]}{('min', 'max')[side]}"
            if kind == "structured":
                og.create_domain(
                    [[names[0]], [names[1]], [(names[2], True)], [(names[3], True)]], name=face
                )
            else:
                og.create_unstructured_domain(connectors=names, name=face)
            faces.append(face)
    return faces


def cylinder(og, n_half=33, n_axial=11):
    """Cylinder x^2 + y^2 = 0.25, z in [-0.5, 0.5]: two half sides on surfaces, two caps."""
    for z, s in ((-0.5, "0"), (0.5, "1")):
        og.create_arc((0, 0, z), 0.5, 0.0, math.pi, name=f"top{s}")
        og.create_arc((0, 0, z), 0.5, math.pi, 2 * math.pi, name=f"bot{s}")
        og.create_connector(f"top{s}", n_half, name=f"c-top{s}")
        og.create_connector(f"bot{s}", n_half, name=f"c-bot{s}")
    for x, s in ((0.5, "e"), (-0.5, "w")):
        og.create_line((x, 0, -0.5), (x, 0, 0.5), name=f"seam-{s}")
        og.create_connector(f"seam-{s}", n_axial, name=f"c-seam-{s}")
    og.create_extruded_surface("top0", (0, 0, 1), name="cyl-top")
    og.create_extruded_surface("bot0", (0, 0, 1), name="cyl-bot")
    og.create_unstructured_domain(
        connectors=["c-top0", "c-top1", "c-seam-e", "c-seam-w"], surface="cyl-top", name="side-top"
    )
    og.create_unstructured_domain(
        connectors=["c-bot0", "c-bot1", "c-seam-e", "c-seam-w"], surface="cyl-bot", name="side-bot"
    )
    og.create_unstructured_domain(connectors=["c-top0", "c-bot0"], name="cap0")
    og.create_unstructured_domain(connectors=["c-top1", "c-bot1"], name="cap1")
    return ["side-top", "side-bot", "cap0", "cap1"]


def build() -> opengrid.Session:
    og = opengrid.Session()
    body = cylinder(og)
    near = box(og, (-1.5, -1.5, -1.5), (3.0, 1.5, 1.5), 8, "near", "unstructured")
    far = box(og, (-3.0, -3.0, -2.0), (8.0, 3.0, 2.0), 2, "chan", "structured")
    og.create_volume(near + body, generator={"max_size": 0.12, "grading": 0.3}, name="near")
    og.create_volume(far + near, generator={"max_size": 0.6, "grading": 0.3}, name="far")
    og.create_patch("cylinder", "wall", body)
    og.create_patch("inlet", "inlet", ["chan-xmin"])
    og.create_patch("outlet", "outlet", ["chan-xmax"])
    og.create_patch("sides", "symmetry", ["chan-ymin", "chan-ymax", "chan-zmin", "chan-zmax"])
    og.create_zone("fluid", ["near", "far"])
    return og


def unstructured_grid(mesh):
    import pyvista as pv

    cells = [np.hstack([np.full((len(mesh["tets"]), 1), 4), mesh["tets"]]).ravel()]
    types = [np.full(len(mesh["tets"]), pv.CellType.TETRA, dtype=np.uint8)]
    if len(mesh["pyramids"]):
        cells.append(np.hstack([np.full((len(mesh["pyramids"]), 1), 5), mesh["pyramids"]]).ravel())
        types.append(np.full(len(mesh["pyramids"]), pv.CellType.PYRAMID, dtype=np.uint8))
    return pv.UnstructuredGrid(np.concatenate(cells), np.concatenate(types), mesh["points"])


def render(og, path):
    import pyvista as pv

    plotter = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(2400, 1000), border=False)
    colours = {"near": "#f5a623", "far": "#7fd1b9"}
    for k, (title, zoom) in enumerate(
        (("tetrahedra cut at z = 0", 0.55), ("near the cylinder", 3.2))
    ):
        plotter.subplot(0, k)
        plotter.set_background("#1f2329")
        for name in ("near", "far"):
            grid = unstructured_grid(og.volume_mesh(name))
            cut = grid.clip(normal=(0, 0, 1), origin=(0, 0, 0.0), crinkle=True)
            plotter.add_mesh(cut, color=colours[name], show_edges=True, edge_color="#30363d")
        for name in ("side-top", "side-bot"):
            m = og.unstructured_mesh(name)
            from opengrid.app.scene import triangle_mesh

            plotter.add_mesh(triangle_mesh(m["points"], m["triangles"]), color="#c0c0c0")
        plotter.add_text(title, font_size=11, color="white")
        plotter.camera_position = [(2.5, -9.0, 9.0), (2.5, 0.0, 0.0), (0, 0, 1)]
        plotter.camera.zoom(zoom)
        if k == 1:
            plotter.camera.focal_point = (0.3, 0.0, 0.0)
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
    og = build()
    print("cylinder in a channel:")
    total = 0.0
    for name in ("near", "far"):
        start = time.perf_counter()
        mesh = og.volume_mesh(name)
        q = og.volume_quality(name)
        n = len(mesh["tets"])
        total += total_volume(mesh)
        print(
            f"  {name}: {og.volume_status(name)} ({time.perf_counter() - start:.1f} s); "
            f"scaled Jacobian min {q['scaled_jacobian'].min():.3f}, min dihedral "
            f"{q['min_angle'][0, :n].min():.1f} deg, max dihedral {q['max_angle'][0, :n].max():.1f}"
            f" deg"
        )
    analytic = 11 * 6 * 4 - math.pi * 0.25 * 1.0
    print(
        f"  volume {total:.6f}, analytic {analytic:.6f} (cylinder as a polygon: "
        f"{(total - analytic) / analytic:+.2e})"
    )
    report = og.boundary_report()
    assert report.ok, report.issues
    print(f"  interfaces: {len(report.interfaces)} (refinement box faces)")
    if _core.cgns_available():
        cgns = OUTPUT / "08a_volumes.cgns"
        print(f"  CGNS {og.export_cgns(cgns)}")
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run([checker, str(cgns)], capture_output=True, text=True).stdout
            status = "clean"
            if "error" in out.lower() or "warning" in out.lower():
                status = "PROBLEMS:\n" + out
            print(f"  cgnscheck: {status}")
    print(f"  VTKHDF {og.export_vtkhdf(OUTPUT / '08a_volumes.vtkhdf')}")
    print(f"  Gmsh {og.export_gmsh(OUTPUT / '08a_volumes.msh')}")
    og.save(OUTPUT / "08a_volumes.ogp")
    og.save_journal(OUTPUT / "08a_volumes_journal.py")
    replayed = opengrid.replay((OUTPUT / "08a_volumes_journal.py").read_text())
    assert to_json(replayed.store) == to_json(og.store)
    render(og, OUTPUT / "08a_volumes.png")
    if args.gui:
        return run_gui(OUTPUT / "08a_volumes.ogp", args.screenshot)
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
    window.set_selection([window.session.entity("near").id])
    window.quality_dock.metric.setCurrentIndex(0)
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
