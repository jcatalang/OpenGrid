"""Phase 8b demo: 3D boundary-layer inflation (prisms / hexahedra) in unstructured volumes.

Case "sphere": a sphere (radius 0.5) in a box [-2.5, 5] x [-2.5, 2.5]^2. The sphere surface is a
cubed sphere: six unstructured domains, each on a revolved surface whose poles and seam lie
outside the patch. Prism layers (first height 2e-3, growth 1.2, at most 15 layers, stop at the
core size) and Netgen tetrahedra in the core.

Case "wing": a NACA 0012 wing stub (root chord 1, tip chord 0.6, LE sweep 0.3 over a span of 1.5)
on a symmetry plane, in a box [-2, 4] x [-2.5, 2.5] x [0, 3]:
  * the upper and lower surfaces are structured domains on ruled surfaces: hexahedral layers;
  * the tip cap is a planar unstructured domain: prism layers;
  * the symmetry plane is a planar unstructured domain with the root section as a hole: the
    layers slide along it and it is re-meshed around them;
  * at the sharp trailing edge a single marching normal is not valid: those columns stop and
    the layers taper towards them (pyramids and tetrahedra at the steps).
Both cases: exports CGNS (mixed element sections and BCs; cgnscheck when available), VTKHDF and
Gmsh 4.1; writes 08b_inflation.png, projects and journals.

    python examples/08b_inflation.py [--gui [--screenshot out.png]]
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
from opengrid.model.volumes import cell_volumes, total_volume, volume_cells
from opengrid.project import to_json

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"


def box_faces(og, lo, hi, per_unit, tag, plane_hole=None):
    """Six faces of a box (structured); with plane_hole = (name, connectors) the zmin face is an
    unstructured domain with those connectors as a hole."""
    corners = {
        c: tuple(lo[k] + c[k] * (hi[k] - lo[k]) for k in range(3))
        for c in itertools.product((0, 1), repeat=3)
    }
    edges = {}

    def edge(a, b):
        key = tuple(sorted((a, b)))
        if key not in edges:
            name = f"{tag}{len(edges)}"
            og.create_line(corners[key[0]], corners[key[1]], name=name)
            n = max(3, round(math.dist(corners[key[0]], corners[key[1]]) * per_unit) + 1)
            edges[key] = og.create_connector(name, n, name=f"c{name}")
        return edges[key], key != (a, b)

    loops = {
        "xmin": [(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)],
        "xmax": [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)],
        "ymin": [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)],
        "ymax": [(0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1)],
        "zmin": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)],
        "zmax": [(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)],
    }
    faces = {}
    for side, cs in loops.items():
        ed = [edge(cs[k], cs[(k + 1) % 4]) for k in range(4)]
        if side == "zmin" and plane_hole is not None:
            name, hole = plane_hole
            og.create_unstructured_domain(connectors=[e[0] for e in ed] + hole, name=name)
            faces[side] = name
            continue
        og.create_domain([[e] for e in ed], name=f"{tag}_{side}")
        faces[side] = f"{tag}_{side}"
    return faces


def sphere_case() -> opengrid.Session:
    og = opengrid.Session()
    r = 0.5
    # Meridians (south to north pole) about z and about x, revolved into three sphere surfaces:
    # +-x patches use the z-axis sphere with its seam at +y, +-y patches the one with its seam
    # at +x, +-z patches the x-axis sphere (poles at +-x, seam at +y).
    og.create_arc(
        (0, 0, 0), r, -math.pi / 2, math.pi / 2, normal=(0, -1, 0), x_axis=(1, 0, 0), name="mer_z"
    )
    og.create_arc(
        (0, 0, 0), r, -math.pi / 2, math.pi / 2, normal=(0, 0, -1), x_axis=(0, 1, 0), name="mer_x"
    )
    og.create_revolved_surface(
        "mer_z", 2.5 * math.pi, axis=(0, 0, 1), start_angle=0.5 * math.pi, name="sphere_x"
    )
    og.create_revolved_surface("mer_z", 2 * math.pi, axis=(0, 0, 1), name="sphere_y")
    og.create_revolved_surface("mer_x", 2 * math.pi, axis=(1, 0, 0), name="sphere_z")
    corner = {
        c: tuple(r * v / math.sqrt(3) for v in c) for c in itertools.product((-1, 1), repeat=3)
    }
    arcs = {}

    def arc(a, b):
        key = tuple(sorted((a, b)))
        if key not in arcs:
            pa = np.array(corner[key[0]])
            pb = np.array(corner[key[1]])
            mid = r * (pa + pb) / np.linalg.norm(pa + pb)
            name = f"arc{len(arcs)}"
            og.create_arc_3pt(tuple(pa), tuple(mid), tuple(pb), name=name)
            arcs[key] = og.create_connector(name, 17, name=f"c{name}")
        return arcs[key]

    patches = []
    for axis, surface in ((0, "sphere_x"), (1, "sphere_y"), (2, "sphere_z")):
        for side in (-1, 1):
            a, b = [k for k in range(3) if k != axis]
            cs = []
            for i, j in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                c = [0, 0, 0]
                c[axis] = side
                c[a] = i
                c[b] = j
                cs.append(tuple(c))
            name = f"sphere_{'xyz'[axis]}{('min', 'max')[side > 0]}"
            og.create_unstructured_domain(
                connectors=[arc(cs[k], cs[(k + 1) % 4]) for k in range(4)],
                surface=surface,
                name=name,
            )
            patches.append(name)
    outer = box_faces(og, (-2.5, -2.5, -2.5), (5.0, 2.5, 2.5), 1.6, "far")
    og.create_volume(
        list(outer.values()) + patches,
        generator={
            "max_size": 0.6,
            "inflation": {
                "faces": patches,
                "first_height": 2e-3,
                "growth_rate": 1.2,
                "layers": 15,
            },
        },
        name="fluid",
    )
    og.create_patch("sphere", "wall", patches)
    og.create_patch("inlet", "inlet", [outer["xmin"]])
    og.create_patch("outlet", "outlet", [outer["xmax"]])
    og.create_patch("sides", "farfield", [outer[k] for k in ("ymin", "ymax", "zmin", "zmax")])
    return og


def naca0012(x):
    return 0.6 * (
        0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4
    )


def wing_case() -> opengrid.Session:
    og = opengrid.Session()
    span, tip, sweep = 1.5, 0.6, 0.3
    xs = [0.5 * (1.0 - math.cos(math.pi * i / 100)) for i in range(101)]
    for z, chord, x0, tag in ((0.0, 1.0, 0.0, "root"), (span, tip, sweep, "tip")):
        for sign, side in ((1, "upper"), (-1, "lower")):
            og.create_spline(
                [(x0 + chord * x, sign * chord * naca0012(x), z) for x in xs],
                "clamped",
                start_tangent=(0, sign, 0),
                name=f"{tag}_{side}",
            )
            distribution = {
                "type": "tanh",
                "begin_spacing": 0.004 * chord,
                "end_spacing": 0.006 * chord,
            }
            og.create_connector(f"{tag}_{side}", 41, distribution, name=f"c_{tag}_{side}")
    og.create_line((0, 0, 0), (sweep, 0, span), name="le")
    og.create_connector("le", 25, name="c_le")
    og.create_line((1, 0, 0), (sweep + tip, 0, span), name="te")
    og.create_connector("te", 25, name="c_te")
    for side in ("upper", "lower"):
        og.create_ruled_surface(f"root_{side}", f"tip_{side}", name=f"surface_{side}")
        og.create_domain(
            [[f"c_root_{side}"], ["c_te"], [(f"c_tip_{side}", True)], [("c_le", True)]],
            surface=f"surface_{side}",
            name=f"wing_{side}",
        )
    og.create_unstructured_domain(connectors=["c_tip_upper", "c_tip_lower"], name="wing_tip")
    outer = box_faces(
        og,
        (-2.0, -2.5, 0.0),
        (4.0, 2.5, 3.0),
        3,
        "far",
        plane_hole=("symmetry_plane", ["c_root_upper", "c_root_lower"]),
    )
    walls = ["wing_upper", "wing_lower", "wing_tip"]
    og.create_volume(
        list(outer.values()) + walls,
        generator={
            "max_size": 0.5,
            "inflation": {"faces": walls, "first_height": 1e-3, "growth_rate": 1.2, "layers": 20},
        },
        name="fluid",
    )
    og.create_patch("wing", "wall", walls)
    og.create_patch("symmetry", "symmetry", [outer["zmin"]])
    og.create_patch("inlet", "inlet", [outer["xmin"]])
    og.create_patch("outlet", "outlet", [outer["xmax"]])
    og.create_patch("sides", "farfield", [outer[k] for k in ("ymin", "ymax", "zmax")])
    return og


def summary(og, label, analytic):
    start = time.perf_counter()
    mesh = og.volume_mesh("fluid")
    elapsed = time.perf_counter() - start
    q = og.volume_quality("fluid")
    info = mesh["inflation"]
    print(f"{label}: {og.volume_status('fluid')} ({elapsed:.1f} s)")
    print(f"  layers per wall node: max {info['max_layers']}; stops {info['stops']}")
    if info["sliding"]:
        names = ", ".join(og.store.get(d).name for d in info["sliding"])
        print(f"  layers slide on: {names} (re-meshed around them)")
    start_index = 0
    for key, cells, _ in volume_cells(mesh):
        n = len(cells)
        if n:
            part = slice(start_index, start_index + n)
            assert np.all(cell_volumes(mesh["points"], key, cells) > 0)
            sj = q["scaled_jacobian"][0, part].min()
            lo = q["min_angle"][0, part].min()
            hi = q["max_angle"][0, part].max()
            print(
                f"  {key:9s} {n:7d}: scaled Jacobian min {sj:.3f}, dihedral {lo:.1f}..{hi:.1f} deg"
            )
        start_index += n
    filled = total_volume(mesh)
    text = f"  volume {filled:.9f} = enclosed by the boundary {mesh['enclosed_volume']:.9f}"
    if analytic is not None:
        text += f"; analytic {analytic:.9f} (faceted wall: {(filled - analytic) / analytic:+.2e})"
    print(text)
    report = og.boundary_report()
    assert report.ok, report.issues
    return mesh


def export(og, stem):
    if _core.cgns_available():
        path = OUTPUT / f"{stem}.cgns"
        print(f"  CGNS {og.export_cgns(path)}")
        checker = shutil.which("cgnscheck")
        if checker:
            out = subprocess.run([checker, str(path)], capture_output=True, text=True).stdout
            status = "clean"
            if "error" in out.lower() or "warning" in out.lower():
                status = "PROBLEMS:\n" + out
            print(f"  cgnscheck: {status}")
    print(f"  VTKHDF {og.export_vtkhdf(OUTPUT / f'{stem}.vtkhdf')}")
    print(f"  Gmsh {og.export_gmsh(OUTPUT / f'{stem}.msh')}")
    og.save(OUTPUT / f"{stem}.ogp")
    og.save_journal(OUTPUT / f"{stem}_journal.py")
    replayed = opengrid.replay((OUTPUT / f"{stem}_journal.py").read_text())
    assert to_json(replayed.store) == to_json(og.store)


def render(cases, path):
    import pyvista as pv

    from opengrid.app.scene import volume_grid

    colours = ["#7fd1b9", "#e8c547", "#f5a623", "#d95f5f"]
    views = [
        ("sphere", (0, 0, 1), (0.3, 0.0, 0.0), (0, 1, 0), 2.6, "sphere: cut z = 0"),
        ("wing", (0, 0, 1), (0.45, 0.0, 0.5), (0, 1, 0), 3.2, "wing: cut z = 0.5"),
        ("wing", (0, 0, 1), (0.95, 0.0, 0.5), (0, 1, 0), 22.0, "trailing edge"),
        (
            "wing",
            (1, 0, 0),
            (0.4, 0.0, 1.0),
            (0, 1, 0),
            2.0,
            "wing: cut x = 0.4 (root on the right, tip)",
        ),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 1600), border=False)
    for k, (case, normal, origin, up, zoom, title) in enumerate(views):
        plotter.subplot(k // 2, k % 2)
        plotter.set_background("#1f2329")
        mesh = cases[case].volume_mesh("fluid")
        grid = volume_grid(mesh)
        grid.cell_data["type"] = np.concatenate(
            [np.full(len(c), t) for t, (_, c, _) in enumerate(volume_cells(mesh))]
        )
        cut = grid.slice(normal=normal, origin=origin)
        plotter.add_mesh(
            cut,
            scalars="type",
            cmap=colours,
            clim=[0, 3],
            show_edges=True,
            edge_color="#30363d",
            show_scalar_bar=False,
        )
        plotter.add_text(title, font_size=11, color="white")
        camera = np.asarray(origin) + 10.0 * np.asarray(normal)
        plotter.camera_position = [tuple(camera), origin, up]
        plotter.camera.zoom(zoom)
    plotter.add_text(
        "tetrahedra / pyramids / prisms / hexahedra", position="lower_right", font_size=9
    )
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
    cases = {"sphere": sphere_case(), "wing": wing_case()}
    summary(cases["sphere"], "sphere in a box", 7.5 * 5 * 5 - 4.0 / 3.0 * math.pi * 0.5**3)
    export(cases["sphere"], "08b_sphere")
    wing = cases["wing"]
    mesh = summary(wing, "wing stub", None)
    walls = [wing.entity(n).id for n in ("wing_upper", "wing_lower")]
    nodes = np.unique(np.concatenate([np.concatenate(mesh["faces"][w]) for w in walls]))
    layers = mesh["inflation"]["layers"][nodes]
    print(
        f"  wing nodes: {np.mean(layers == mesh['inflation']['max_layers']) * 100:.0f}% with all "
        f"{mesh['inflation']['max_layers']} layers, {np.sum(layers == 0)} without (the sharp "
        "trailing edge; the layers taper towards it)"
    )
    export(wing, "08b_wing")
    render(cases, OUTPUT / "08b_inflation.png")
    if args.gui:
        return run_gui(OUTPUT / "08b_wing.ogp", args.screenshot)
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
    window.set_selection([window.session.entity("fluid").id])
    window.quality_action.setChecked(True)
    window.cut_action.setChecked(True)
    window.viewport.camera_position = [(3.2, -2.6, 2.4), (0.6, 0.0, 0.5), (0.0, 0.0, 1.0)]
    if screenshot is not None:

        def save_and_quit():
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(3000, save_and_quit)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
