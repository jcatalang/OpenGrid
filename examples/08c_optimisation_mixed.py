"""Phase 8c demo: ridge fans, volume-mesh optimisation and mixed structured / unstructured models.

Case "wing" (the Phase 8b wing stub, examples/08b_inflation.py): at the sharp trailing edge
(turning 164 degrees) and the 90-degree rim of the tip cap the layers now grow in fans: several
columns from each ridge node, joined by faces of zero area at the wall (prisms lying on their
side in the first layer, hexahedra above); where the trailing edge meets the tip rim a cap around
a centre column closes the corner. The core is then optimised (node moves of interior nodes and
pyramid apexes, 2-3 / 3-2 swaps of tetrahedra); the comparison without optimisation is printed.

Case "duct": a structured block (hexahedra, x in [0, 1]) feeds an unstructured volume (x in
[1, 3]) through a shared structured face (pyramids on the volume side). The volume holds a cube
obstacle with structured faces: hexahedral layers, fans along its edges and caps at its corners.
Exported to CGNS as one base with a structured and an unstructured zone joined by vertex
connectivity (cgnscheck when available), VTKHDF and Gmsh.

Writes 08c_optimisation_mixed.png, projects and journals.

    python examples/08c_optimisation_mixed.py [--gui [--screenshot out.png]]
"""

import argparse
import importlib
import itertools
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
sys.path.insert(0, str(HERE))
phase8b = importlib.import_module("08b_inflation")


def box_domains(og, lo, hi, n, tag, shared=None):
    """Six structured domains of an axis-aligned box, n nodes per edge; returns them by face
    ("xmin", ..., "zmax"). ``shared``: existing domains to use for some faces."""
    corners = {
        c: tuple(lo[k] + c[k] * (hi[k] - lo[k]) for k in range(3))
        for c in itertools.product((0, 1), repeat=3)
    }
    edges = {}

    def connector(a, b):
        key = tuple(sorted((corners[a], corners[b])))
        if key not in edges:
            existing = [
                name
                for name in og.names("line")
                if {og.entity(name).params["p0"], og.entity(name).params["p1"]} == set(key)
            ]
            if existing:
                edges[key] = "c_" + existing[0]
            else:
                name = f"{tag}_e{len(edges)}"
                og.create_line(key[0], key[1], name=name)
                edges[key] = og.create_connector(name, n, name=f"c_{name}")
        line = og.entity(edges[key][2:])
        return edges[key], line.params["p0"] != corners[a]

    faces = {}
    for axis, label in enumerate("xyz"):
        for side, end in ((0, "min"), (1, "max")):
            name = f"{label}{end}"
            if shared and name in shared:
                faces[name] = shared[name]
                continue
            a, b = [k for k in range(3) if k != axis]

            def v(i, j, axis=axis, side=side, a=a, b=b):
                c = [0, 0, 0]
                c[axis] = side
                c[a] = i
                c[b] = j
                return tuple(c)

            loop = [v(0, 0), v(1, 0), v(1, 1), v(0, 1)]
            sides = [connector(loop[k], loop[(k + 1) % 4]) for k in range(4)]
            faces[name] = og.create_domain([[s] for s in sides], name=f"{tag}_{name}")
    return faces


def duct_case() -> opengrid.Session:
    og = opengrid.Session()
    block = box_domains(og, (0, 0, 0), (1, 1, 1), 9, "inlet")
    og.create_block(faces=list(block.values()), name="inlet_block")
    duct = box_domains(og, (1, 0, 0), (3, 1, 1), 9, "duct", shared={"xmin": block["xmax"]})
    cube = box_domains(og, (1.8, 0.35, 0.35), (2.1, 0.65, 0.65), 11, "cube")
    og.create_volume(
        list(duct.values()) + list(cube.values()),
        generator={
            "max_size": 0.08,
            "inflation": {
                "faces": list(cube.values()),
                "first_height": 2e-3,
                "growth_rate": 1.2,
                "layers": 10,
            },
        },
        name="fluid",
    )
    og.create_patch("inlet", "inlet", [block["xmin"]])
    og.create_patch("outlet", "outlet", [duct["xmax"]])
    walls = [block[k] for k in ("ymin", "ymax", "zmin", "zmax")]
    walls += [duct[k] for k in ("ymin", "ymax", "zmin", "zmax")]
    og.create_patch("duct", "wall", walls)
    og.create_patch("obstacle", "wall", list(cube.values()))
    return og


def quality_table(mesh):
    """Per cell type, layer and core cells: count, min scaled Jacobian, cells below 0.1 / 0.3."""
    points = np.ascontiguousarray(mesh["points"])
    layer = (mesh["inflation"] or {"cells": {}})["cells"]
    rows = []
    for key, cells, _ in volume_cells(mesh):
        q = np.array([_core.cell_quality_3d(points, c.tolist()) for c in cells])
        n_layer = layer.get(key, 0)
        for part, values in (("layer", q[:n_layer]), ("core", q[n_layer:])):
            if len(values):
                rows.append(
                    (
                        key,
                        part,
                        len(values),
                        values.min(),
                        np.sum(values < 0.1),
                        np.sum(values < 0.3),
                    )
                )
    return rows


def print_quality(rows, indent="  "):
    for key, part, n, worst, below01, below03 in rows:
        print(
            f"{indent}{part:5s} {key:9s} {n:7d}: scaled Jacobian min {worst:.3f}, "
            f"{below01} below 0.1, {below03} below 0.3"
        )


def wing_report(og):
    start = time.perf_counter()
    mesh = og.volume_mesh("fluid")
    elapsed = time.perf_counter() - start
    info = mesh["inflation"]
    print(f"wing stub: {og.volume_status('fluid')} ({elapsed:.1f} s)")
    print(f"  fan columns {info['fan_columns']}; column stops {info['stops']}")
    n_surface = mesh["n_surface"]
    reason = np.array(info["reason"][:n_surface], dtype=object)
    wall = reason != ""
    layers = info["layers"][:n_surface][wall]
    points = mesh["points"][:n_surface][wall]
    print(
        f"  wall nodes: {np.mean(layers == info['max_layers']) * 100:.0f}% with all "
        f"{info['max_layers']} layers (Phase 8b: 43%), mean {layers.mean():.1f}"
    )
    # Trailing edge: the line from (1, 0, 0) to (0.9, 0, 1.5).
    a = np.array([1.0, 0.0, 0.0])
    d = np.array([-0.1, 0.0, 1.5])
    t = np.clip((points - a) @ d / (d @ d), 0.0, 1.0)
    on_edge = np.linalg.norm(points - (a + t[:, None] * d), axis=1) < 1e-9
    print(
        f"  trailing-edge nodes: {np.sum(layers[on_edge] == info['max_layers'])} of "
        f"{on_edge.sum()} with all layers (Phase 8b: none)"
    )
    print(
        f"  volume {total_volume(mesh):.9f} = enclosed by the boundary "
        f"{mesh['enclosed_volume']:.9f}"
    )
    stats = mesh["optimization"]
    print(
        f"  optimisation: {stats['moved']} node moves, {stats['swaps_23']} 2-3 and "
        f"{stats['swaps_32']} 3-2 swaps; cells below 0.3: {stats['below_before']} -> "
        f"{stats['below_after']}"
    )
    rows = quality_table(mesh)
    print("  with optimisation:")
    print_quality(rows, "    ")
    og.set_volume_generator("fluid", optimize=False)
    print("  without optimisation (core cells):")
    print_quality([r for r in quality_table(og.volume_mesh("fluid")) if r[1] == "core"], "    ")
    og.undo()
    for key, cells, _ in volume_cells(mesh):
        assert np.all(cell_volumes(mesh["points"], key, cells) > 0), key
    return mesh


def duct_report(og):
    start = time.perf_counter()
    mesh = og.volume_mesh("fluid")
    elapsed = time.perf_counter() - start
    info = mesh["inflation"]
    print(f"duct: block inlet_block {og.block_status('inlet_block')}")
    print(f"  volume fluid: {og.volume_status('fluid')} ({elapsed:.1f} s)")
    print(f"  fan columns {info['fan_columns']}; column stops {info['stops'] or 'none'}")
    print(f"  optimisation {mesh['optimization']}")
    print_quality(quality_table(mesh))
    expected = 2.0 - 0.3**3
    print(f"  volume {total_volume(mesh):.12f} (analytic {expected:.12f})")
    for key, cells, _ in volume_cells(mesh):
        assert np.all(cell_volumes(mesh["points"], key, cells) > 0), key
    report = og.boundary_report()
    assert report.ok, report.issues
    print(f"  interfaces: {report.interfaces}")
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


def render(wing, duct, path):
    import pyvista as pv

    from opengrid.app.scene import volume_grid

    colours = ["#7fd1b9", "#e8c547", "#f5a623", "#d95f5f"]

    def typed(mesh):
        grid = volume_grid(mesh)
        grid.cell_data["type"] = np.concatenate(
            [np.full(len(c), t) for t, (_, c, _) in enumerate(volume_cells(mesh))]
        )
        return grid

    wing_grid = typed(wing.volume_mesh("fluid"))
    duct_grid = typed(duct.volume_mesh("fluid"))
    block = pv.StructuredGrid(*np.moveaxis(duct.block_grid("inlet_block"), -1, 0)[::-1])
    block = block.cast_to_unstructured_grid()
    block.cell_data["type"] = np.full(block.n_cells, 3)
    views = [
        (
            wing_grid,
            (0, 0, 1),
            (0.95, 0.0, 0.5),
            (0, 1, 0),
            9.0,
            "wing: trailing-edge fan, z = 0.5",
        ),
        (wing_grid, (0, 0, 1), (0.92, 0.0, 1.2), (0, 1, 0), 6.0, "wing: z = 1.2 (towards the tip)"),
        (wing_grid, (1, 0, 0), (0.83, 0.0, 1.5), (0, 1, 0), 9.0, "wing: x = 0.83, tip rim fans"),
        (None, (0, 0, 1), (1.5, 0.5, 0.5), (0, 1, 0), 2.6, "duct: z = 0.5 (block | volume)"),
    ]
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 1600), border=False)
    for k, (grid, normal, origin, up, zoom, title) in enumerate(views):
        plotter.subplot(k // 2, k % 2)
        plotter.set_background("#1f2329")
        grids = [grid]
        if grid is None:
            grids = [duct_grid, block]
        for g in grids:
            plotter.add_mesh(
                g.slice(normal=normal, origin=origin),
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
    wing = phase8b.wing_case()
    wing_report(wing)
    export(wing, "08c_wing")
    duct = duct_case()
    duct_report(duct)
    export(duct, "08c_duct")
    render(wing, duct, OUTPUT / "08c_optimisation_mixed.png")
    if args.gui:
        return phase8b.run_gui(OUTPUT / "08c_duct.ogp", args.screenshot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
