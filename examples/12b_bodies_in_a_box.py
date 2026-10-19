"""Bodies in a box: 3D volume meshes around the standard shapes, as a user would script them.

Case "sphere": a sphere of radius 0.5 (create_sphere: two halves) in a 5 x 4 x 4 box.
Case "wing": a closed NACA 0012 wing (create_wing: root chord 1, tip chord 0.6, span 1.5,
15 degrees of sweep, root and tip caps) in an 8 x 6 x 6 box.

For both:
  * mesh_surfaces puts surface meshes on the body (fine) and on the box faces (coarse): the edge
    curves are found or made, the connectors shared, so each set is a closed shell;
  * create_volume fills the space between the shells with Netgen tetrahedra and grows prism
    layers off the body (pyramids where they meet quads, fans at the wing's trailing edge);
  * set_boundary_type gives the body a wall and the box a far field, the boundary report is
    checked, CGNS and Gmsh are written.

Writes 12b_bodies_in_a_box.png (cuts coloured by cell type), the projects, journals and files.

    python examples/12b_bodies_in_a_box.py
"""

import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid import _core

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
BOX_FACES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")


def body_in_a_box(og, body_surfaces, box: str, body_spacing, box_spacing, layers) -> str:
    """Surface meshes on the body and the box, the volume between them with prism layers off
    the body, the boundary types; returns the volume's name."""
    body = og.mesh_surfaces(body_surfaces, spacing=body_spacing)
    outer = og.mesh_surfaces([f"{box}-{s}" for s in BOX_FACES], spacing=box_spacing)
    volume = og.create_volume(
        outer + body,
        generator={"max_size": box_spacing, "inflation": {"faces": body, **layers}},
        name="fluid",
    )
    og.set_boundary_type(body, "wall")
    og.set_boundary_type(outer, "farfield")
    return volume


def sphere_case() -> opengrid.Session:
    og = opengrid.Session()
    parts = og.create_sphere(radius=0.5, name="sphere")
    og.create_box(center=(0.5, 0.0, 0.0), size=(5.0, 4.0, 4.0), name="box")
    surfaces = [n for n in parts if og.entity(n).kind == "revolved_surface"]
    layers = {"first_height": 2e-3, "growth_rate": 1.2, "layers": 12}
    body_in_a_box(og, surfaces, "box", 0.04, 0.5, layers)
    return og


def wing_case() -> opengrid.Session:
    og = opengrid.Session()
    parts = og.create_wing("0012", root_chord=1.0, tip_chord=0.6, span=1.5, sweep=15.0, name="wing")
    og.create_box(center=(1.0, 0.75, 0.0), size=(8.0, 6.0, 6.0), name="box")
    surfaces = [n for n in parts if og.entity(n).kind == "ruled_surface"]
    layers = {"first_height": 2e-3, "growth_rate": 1.25, "layers": 8}
    body_in_a_box(og, surfaces, "box", 0.04, 0.6, layers)
    return og


def summary(og, title: str) -> None:
    report = og.boundary_report()
    if not report.ok:
        for issue in report.issues:
            print("  ", issue)
        raise SystemExit(f"{title}: the boundary is not ready")
    quality = og.volume_quality("fluid")
    jacobian = np.concatenate([np.asarray(v).ravel() for v in [quality["scaled_jacobian"]]])
    print(f"{title}: {og.volume_status('fluid')}")
    print(f"  min scaled Jacobian {jacobian.min():.3f}; {len(report.boundary)} boundaries typed")


def export(og, stem: str) -> None:
    og.save(OUTPUT / f"{stem}.ogp")
    (OUTPUT / f"{stem}_journal.py").write_text(og.journal_text())
    og.export_gmsh(OUTPUT / f"{stem}.msh")
    if _core.cgns_available():
        og.export_cgns(OUTPUT / f"{stem}.cgns")


def render(cases: dict, path: Path) -> None:
    import pyvista as pv

    from opengrid.app.scene import volume_grid
    from opengrid.model.volumes import volume_cells

    colours = ["#7fd1b9", "#e8c547", "#f5a623", "#d95f5f"]
    views = [
        ("sphere", (0, 0, 1), (0.5, 0.0, 0.0), (0, 1, 0), 1.6, "sphere: cut z = 0"),
        ("sphere", (0, 0, 1), (0.0, 0.5, 0.0), (0, 1, 0), 0.12, "sphere: 12 prism layers"),
        # Mid-span (y = 0.75): the leading edge swept back to x = 0.2, the chord 0.8.
        # Looking along +y (from the root side): x to the right, the flow left to right.
        ("wing", (0, -1, 0), (0.6, 0.75, 0.0), (0, 0, 1), 0.9, "wing: cut at mid-span"),
        ("wing", (0, -1, 0), (1.0, 0.75, 0.0), (0, 0, 1), 0.08, "wing: trailing edge fan"),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 1600), border=False)
    for k, (case, normal, origin, up, half, title) in enumerate(views):
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
        plotter.enable_parallel_projection()
        plotter.camera.parallel_scale = half
    plotter.add_text(
        "tetrahedra / pyramids / prisms / hexahedra", position="lower_right", font_size=9
    )
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
    cases = {"sphere": sphere_case(), "wing": wing_case()}
    for name, og in cases.items():
        summary(og, name)
        export(og, f"12b_{name}")
    render(cases, OUTPUT / "12b_bodies_in_a_box.png")
    print(f"done in {time.perf_counter() - started:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
