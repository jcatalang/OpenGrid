"""Phase 9a demo: STEP / IGES import through OpenCASCADE, with connectors and domains created
automatically from the CAD edges and faces.

The part: a block 1.2 x 0.8 x 0.4 with its four vertical edges filleted (radius 0.12) and a
through hole (radius 0.15), modelled with OCCT (netgen.occ) and written to
output/09a_bracket.step. Then:
  * "solid": imported with domains="auto" (four-sided faces become structured TFI domains, on
    the face's NURBS surface when curved; the others unstructured) and filled with tetrahedra;
  * "fluid": imported with unstructured domains and placed in a box [-1.2, 3] x [-1, 1.8] x
    [-0.8, 1.2]; boundary layers (first height 2e-3, growth 1.2, 12 layers) on every part face
    and tetrahedra in the core; patches, CGNS / VTKHDF / Gmsh exports, project and journal
    (replayed without the STEP file);
  * an IGES file of unsewn trimmed faces (tests/data/block_hole.igs) is sewn at import.
Writes 09a_cad_import.png.

    python examples/09a_cad_import.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid.model.volumes import total_volume, volume_cells

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
DATA = HERE.parent / "tests" / "data"
sys.path.insert(0, str(HERE))
phase8b = importlib.import_module("08b_inflation")

STEP = OUTPUT / "09a_bracket.step"


def make_part(path: Path) -> float:
    """Write the bracket to STEP with OCCT (netgen.occ); returns its exact volume."""
    import netgen.occ as occ

    box = occ.Box((0, 0, 0), (1.2, 0.8, 0.4))
    vertical = [e for e in box.edges if e.bounding_box[1][2] - e.bounding_box[0][2] > 0.3]
    body = box.MakeFillet(vertical, 0.12)
    hole = occ.Cylinder((0.6, 0.4, -0.1), (0, 0, 1), r=0.15, h=0.6)
    part = body - hole
    part.WriteStep(str(path))
    return part.mass


def solid_case(exact_volume: float) -> opengrid.Session:
    og = opengrid.Session()
    result = og.import_cad(STEP, name="bracket", domains="auto", spacing=0.04)
    print(result.report())
    og.create_volume(result.domains, generator={"max_size": 0.06}, name="solid")
    og.create_patch("bracket", "wall", result.domains)
    start = time.perf_counter()
    mesh = og.volume_mesh("solid")
    elapsed = time.perf_counter() - start
    print(f"solid: {og.volume_status('solid')} ({elapsed:.1f} s)")
    filled = total_volume(mesh)
    print(
        f"  volume {filled:.6f}; CAD volume {exact_volume:.6f} "
        f"(faceted fillets and hole: {(filled - exact_volume) / exact_volume:+.2e})"
    )
    return og


def fluid_case() -> tuple[opengrid.Session, list[str]]:
    og = opengrid.Session()
    result = og.import_cad(STEP, name="bracket", spacing=0.04)
    walls = result.domains
    outer = phase8b.box_faces(og, (-1.2, -1.0, -0.8), (3.0, 1.8, 1.2), 4, "far")
    og.create_volume(
        list(outer.values()) + walls,
        generator={
            "max_size": 0.3,
            "inflation": {
                "faces": walls,
                "first_height": 2e-3,
                "growth_rate": 1.2,
                "layers": 12,
            },
        },
        name="fluid",
    )
    og.create_patch("bracket", "wall", walls)
    og.create_patch("inlet", "inlet", [outer["xmin"]])
    og.create_patch("outlet", "outlet", [outer["xmax"]])
    og.create_patch("sides", "farfield", [outer[k] for k in ("ymin", "ymax", "zmin", "zmax")])
    return og, walls


def fluid_report(og, walls) -> None:
    start = time.perf_counter()
    mesh = og.volume_mesh("fluid")
    elapsed = time.perf_counter() - start
    print(f"fluid: {og.volume_status('fluid')} ({elapsed:.1f} s)")
    info = mesh["inflation"]
    ids = [og.entity(n).id for n in walls]
    nodes = np.unique(np.concatenate([np.concatenate(mesh["faces"][w]) for w in ids]))
    layers = info["layers"][nodes]
    print(
        f"  wall nodes: mean {layers.mean():.2f} of {info['max_layers']} layers, "
        f"{np.sum(layers == 0)} without; stops {info['stops']}"
    )
    q = mesh["quality"]
    print(
        f"  quality floor: {q['below_floor']} of {q['cells']} cells below it; worst scaled "
        f"Jacobian {q['worst_jacobian']:.3f}, skewness {q['worst_skewness']:.3f}"
    )
    report = og.boundary_report()
    assert report.ok, report.issues


def iges_case() -> None:
    og = opengrid.Session()
    result = og.import_cad(DATA / "block_hole.igs", spacing=0.1)
    unsewn = og.import_cad(DATA / "block_hole.igs", name="unsewn", domains="none", sew=False)
    sewing = [line for line in result.log if "sewing" in line]
    print(f"IGES block_hole.igs: {sewing[0]}")
    print(
        f"  {len(result.domains)} domains, {len(result.free_edges)} free edges after sewing "
        f"(without sewing: {len(unsewn.free_edges)})"
    )


def render(solid, fluid, path) -> None:
    import pyvista as pv

    from opengrid.app.scene import volume_grid

    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 1600), border=False)
    camera = [(2.4, -1.6, 1.6), (0.6, 0.4, 0.2), (0, 0, 1)]

    # Imported faces: structured (blue) and unstructured (orange) domains with their grids.
    plotter.subplot(0, 0)
    plotter.set_background("#1f2329")
    for name in solid.names("domain"):
        grid = solid.domain_grid(name)
        surface = pv.StructuredGrid(grid[..., 0], grid[..., 1], grid[..., 2])
        plotter.add_mesh(surface, color="#6fa8dc", show_edges=True, edge_color="#20303f")
    for name in solid.names("unstructured_domain"):
        mesh = solid.unstructured_mesh(name)
        faces = np.hstack([np.full((len(mesh["triangles"]), 1), 3), mesh["triangles"]])
        surface = pv.PolyData(np.asarray(mesh["points"]), faces.ravel())
        plotter.add_mesh(surface, color="#f5a623", show_edges=True, edge_color="#4a3510")
    plotter.add_text("imported faces: structured / unstructured", font_size=11, color="white")
    plotter.camera_position = camera

    colours = ["#7fd1b9", "#e8c547", "#f5a623", "#d95f5f"]
    views = [
        (solid, "solid", (0, 0, 1), (0.6, 0.4, 0.2), 3.0, "solid: tetrahedra, cut z = 0.2"),
        (fluid, "fluid", (0, 0, 1), (0.6, 0.4, 0.2), 2.2, "fluid: cut z = 0.2"),
        (fluid, "fluid", (0, 0, 1), (1.12, 0.72, 0.2), 9.0, "fillet: layers, cut z = 0.2"),
    ]
    for k, (og, volume, normal, origin, zoom, title) in enumerate(views, start=1):
        plotter.subplot(k // 2, k % 2)
        plotter.set_background("#1f2329")
        mesh = og.volume_mesh(volume)
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
        eye = np.asarray(origin) + 10.0 * np.asarray(normal)
        plotter.camera_position = [tuple(eye), origin, (0, 1, 0)]
        plotter.camera.zoom(zoom)
    plotter.add_text(
        "tetrahedra / pyramids / prisms / hexahedra", position="lower_right", font_size=9
    )
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
    window.set_selection([window.session.entity("fluid").id])
    window.cut_action.setChecked(True)
    window.viewport.camera_position = [(2.6, -1.8, 1.8), (0.6, 0.4, 0.2), (0.0, 0.0, 1.0)]
    if screenshot is not None:

        def save_and_quit():
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
    exact = make_part(STEP)
    solid = solid_case(exact)
    phase8b.export(solid, "09a_solid")
    fluid, walls = fluid_case()
    fluid_report(fluid, walls)
    phase8b.export(fluid, "09a_fluid")
    iges_case()
    render(solid, fluid, OUTPUT / "09a_cad_import.png")
    if args.gui:
        return run_gui(OUTPUT / "09a_fluid.ogp", args.screenshot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
