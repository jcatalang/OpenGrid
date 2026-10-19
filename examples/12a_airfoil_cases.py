"""Two complete 2D aerodynamic cases built from the standard shapes, as a user would script them.

Case "far field": a NACA 2412 at 4 degrees in a 30 x 20 chord rectangle.
  * create_airfoil gives the upper and lower splines, create_rectangle the four far-field lines;
  * the airfoil connectors cluster at the leading and trailing edges (tanh), the far field is
    coarse;
  * one unstructured domain between the two loops (the airfoil is the hole), with inflation
    layers from the airfoil (first height 2e-5 chord, growth 1.15) in a triangle core;
  * boundary types in one step (wall, farfield), the boundary report checked, CGNS and Gmsh
    written.

Case "cascade": a NACA 6412 at 30 degrees of stagger in a periodic passage of pitch 1.
  * inlet on the left, outlet on the right, the bottom and top paired as periodic with
    make_periodic_pair (the translation is found from the two sides);
  * the 2D model is then extruded one cell thick with extrude_model, the front and back a
    second periodic pair: the 3D mesh that solvers without a 2D mode want.

Writes 12a_airfoil_cases.png, the projects, journals and CGNS / Gmsh files.

    python examples/12a_airfoil_cases.py
"""

import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid import _core

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"

FIRST_HEIGHT = 2e-5  # chords: y+ of about 1 at Re = 3e6
LAYERS = {"first_height": FIRST_HEIGHT, "growth_rate": 1.15, "layers": 40}


def airfoil_connectors(og, name: str, points: int = 121) -> list[str]:
    """Connectors on the airfoil's upper and lower splines, clustered at both edges."""
    clustering = {"type": "tanh", "begin_spacing": 2e-3, "end_spacing": 1e-3}
    return [
        og.create_connector(f"{name}-{side}", points, clustering, name=f"c-{side}")
        for side in ("upper", "lower")
    ]


def side_connectors(og, rectangle: str, counts: dict) -> dict:
    """One uniform connector on each side of a rectangle made by create_rectangle."""
    return {
        side: og.create_connector(f"{rectangle}-{side}", n, name=f"c-{side}")
        for side, n in counts.items()
    }


def far_field_case() -> opengrid.Session:
    og = opengrid.Session()
    og.create_airfoil("2412", chord=1.0, angle=4.0, name="naca")
    og.create_rectangle((0.5, 0.0), 30.0, 20.0, name="far")
    wall = airfoil_connectors(og, "naca")
    far = side_connectors(og, "far", {"bottom": 31, "right": 21, "top": 31, "left": 21})
    og.create_unstructured_domain(
        connectors=[*wall, *far.values()],
        generator={
            "max_size": 2.0,
            "growth_rate": 1.15,
            "inflation": {"connectors": wall, **LAYERS},
        },
        name="flow",
    )
    og.set_boundary_type(wall, "wall")
    og.set_boundary_type(list(far.values()), "farfield")
    return og


def cascade_case() -> opengrid.Session:
    og = opengrid.Session()
    og.create_airfoil("6412", chord=1.0, angle=30.0, name="blade")
    # The passage: one pitch high, centred on the blade's mid-chord.
    og.create_rectangle((0.45, -0.25), 4.0, 1.0, name="passage")
    wall = airfoil_connectors(og, "blade", 101)
    sides = side_connectors(og, "passage", {"bottom": 61, "right": 16, "top": 61, "left": 16})
    og.create_unstructured_domain(
        connectors=[*wall, *sides.values()],
        generator={"max_size": 0.06, "inflation": {"connectors": wall, **LAYERS}},
        name="passage-flow",
    )
    og.set_boundary_type(wall, "wall")
    og.set_boundary_type([sides["left"]], "inlet")
    og.set_boundary_type([sides["right"]], "outlet")
    og.make_periodic_pair([sides["bottom"]], [sides["top"]], name="pitch")
    return og


def summary(og, title: str) -> dict:
    report = og.boundary_report()
    if not report.ok:
        for issue in report.issues:
            print("  ", issue)
        raise SystemExit(f"{title}: the boundary is not ready")
    domain = (og.names("unstructured_domain") or og.names("domain"))[0]
    mesh = og.unstructured_mesh(domain)
    quality = og.domain_quality(domain)
    jacobian = np.asarray(quality["scaled_jacobian"]).ravel()
    print(
        f"{title}: {len(mesh['triangles'])} triangles, {len(mesh['quads'])} quads, "
        f"{mesh['layers']} inflation layers; min scaled Jacobian {jacobian.min():.3f}; "
        f"{len(report.boundary)} boundaries typed, {len(report.periodic)} periodic pair(s)"
    )
    return mesh


def export(og, stem: str) -> None:
    og.save(OUTPUT / f"{stem}.ogp")
    (OUTPUT / f"{stem}_journal.py").write_text(og.journal_text())
    og.export_gmsh(OUTPUT / f"{stem}.msh")
    if _core.cgns_available():
        og.export_cgns(OUTPUT / f"{stem}.cgns")


def polydata(mesh):
    import pyvista as pv

    cells = []
    for kind in ("triangles", "quads"):
        block = np.asarray(mesh[kind], dtype=np.int64)
        if len(block):
            cells.append(np.hstack([np.full((len(block), 1), block.shape[1]), block]).ravel())
    return pv.PolyData(np.asarray(mesh["points"], dtype=float), np.concatenate(cells))


def render(far: dict, cascade: dict, path: Path) -> None:
    import pyvista as pv

    views = [
        (far, (0.5, 0.0), 1.6, (0,), "NACA 2412 at 4 deg (far field at 15 chords)"),
        (far, (0.02, 0.0), 0.07, (0,), "leading edge: 40 layers from 2e-5 chord"),
        (cascade, (0.45, -0.25), 1.3, (-1, 0, 1), "NACA 6412 cascade: one pitch, periodic copies"),
        (cascade, (0.87, -0.5), 0.05, (0,), "cascade trailing edge: layer fan"),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 1600), border=False)
    for k, (mesh, centre, half, pitches, title) in enumerate(views):
        plotter.subplot(k // 2, k % 2)
        plotter.set_background("#1f2329")
        for shift in pitches:
            colour = "#e8c547"
            if shift != 0:
                colour = "#7a6a30"  # the periodic neighbours (pitch 1 above and below)
            piece = polydata(mesh).translate((0.0, float(shift), 0.0), inplace=False)
            plotter.add_mesh(piece, color=colour, show_edges=True, edge_color="#30363d")
        plotter.add_text(title, font_size=11, color="white")
        plotter.camera_position = [
            (centre[0], centre[1], 10.0),
            (centre[0], centre[1], 0.0),
            (0, 1, 0),
        ]
        plotter.enable_parallel_projection()
        plotter.camera.parallel_scale = half
    plotter.screenshot(str(path))
    plotter.close()
    save_compact(path)


def save_compact(path: Path) -> None:
    """Re-save a screenshot as a 128-colour PNG: a third of the size, no visible change."""
    from PIL import Image

    image = Image.open(path).convert("RGB").quantize(colors=128, method=Image.Quantize.MEDIANCUT)
    image.save(path, optimize=True)


def main() -> int:
    OUTPUT.mkdir(exist_ok=True)
    started = time.perf_counter()
    far = far_field_case()
    far_mesh = summary(far, "far field")
    export(far, "12a_far_field")
    cascade = cascade_case()
    cascade_mesh = summary(cascade, "cascade")
    export(cascade, "12a_cascade_2d")
    names = cascade.extrude_model(0.05)  # one cell thick, front / back periodic
    report = cascade.boundary_report()
    print(
        f"cascade extruded: {', '.join(names)}; {len(report.periodic)} periodic pairs, "
        f"boundary ready: {report.ok}"
    )
    export(cascade, "12a_cascade_3d")
    render(far_mesh, cascade_mesh, OUTPUT / "12a_airfoil_cases.png")
    print(f"done in {time.perf_counter() - started:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
