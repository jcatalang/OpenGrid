"""Phase 11f demo: straight-sided high-order elements of every cell type, Gmsh / CGNS / VTKHDF
(``curved=False``; the curved boundary elements of 11g: examples/11g_curved_elements.py).

  * Node orders: OpenGrid's generated node orders against each format's own reference, for
    every element kind (line, triangle, quadrilateral, tetrahedron, pyramid, prism,
    hexahedron): Gmsh (orders 1-6, reference data from Gmsh), VTK Lagrange cells (orders 1-6,
    live) and CGNS (orders 1-4: Gmsh's CGNS writer and the SIDS node descriptions).
  * 2D hybrid (structured O-ring of quadrilaterals around a cylinder, unstructured triangles
    outside): orders 1-4, node counts, Gmsh / CGNS / VTKHDF files; the distance of the wall's
    high-order nodes from the cylinder (straight-sided: the chord sagitta at every order) next
    to the opt-in curved elements of the structured ring (Phase 11e).
  * 3D body in a box with hexahedral layers on its walls (tetrahedra, pyramids, prisms,
    hexahedra): orders 2-4, counts, build and write times.
Writes 11f_high_order_elements.png, 11f_high_order_elements_report.txt and the meshes
11f_cylinder_p4.{msh,cgns,vtkhdf}, 11f_layers_p3.{msh,cgns,vtkhdf}.

    python examples/11f_high_order_elements.py
"""

import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from opengrid.io import lagrange

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
TESTS = HERE.parent / "tests" / "python"
sys.path.insert(0, str(TESTS))
checks = importlib.import_module("test_high_order_elements")
unstructured_cases = importlib.import_module("test_unstructured")
layer_cases = importlib.import_module("test_inflation3d")
curved_cases = importlib.import_module("test_high_order")


def node_order_report() -> list[str]:
    out = ["node orders against the formats' own references"]
    ref = json.loads((TESTS / "data" / "gmsh_reference_nodes.json").read_text())["nodes"]
    gmsh_ok = 0
    for kind in lagrange.KINDS:
        for p in range(1, 7):
            entry = ref[f"{checks.GMSH_FAMILY[kind]}:{p}"]
            d = lagrange.DIM[kind]
            coords = np.array(entry["coords"]).reshape(entry["nodes"], -1)[:, :d]
            want = checks._exact(kind, checks._from_gmsh(kind, coords), p)
            gmsh_ok += int(list(lagrange.nodes("gmsh", kind, p)) == want)
    out.append(f"  Gmsh: {gmsh_ok} of 42 elements (7 kinds x orders 1-6) match")
    vtk_ok = 0
    for kind in ("line", "triangle", "quad", "tetra", "hexa", "prism"):
        try:
            checks.test_vtk_node_orders(kind)
            vtk_ok += 6
        except AssertionError:
            pass
    out.append(
        f"  VTK: {vtk_ok} of 36 Lagrange cells (orders 1-6) match; VTK has no Lagrange pyramid "
        "(pyramids are written as linear pyramids, exact for straight-sided ones)"
    )
    checks.test_cgns_node_orders()
    checks.test_cgns_orders_match_the_sids_descriptions()
    out.append(
        "  CGNS: 26 of 26 element types written by Gmsh's CGNS writer match; PENTA_40 / PENTA_75 "
        "(not written by Gmsh) and TETRA_35 / PYRA_55 agree with all 46 SIDS node descriptions"
    )
    return out


def wall_distance(mesh, radius: float) -> float:
    """Largest distance of the wall's boundary-element nodes from the circle."""
    lines = mesh.faces["line"][mesh.face_patch["line"] == mesh.patch_names.index("p-wall")]
    r = np.linalg.norm(mesh.points[lines.reshape(-1), :2], axis=1)
    return float(np.abs(r - radius).max())


def curved_wall_distance(og, p: int) -> float:
    """Largest distance of the curved elements' wall edges from the cylinder."""
    return curved_cases.arc_error(og.curved_high_order("ring", order=p))


def main() -> int:
    OUTPUT.mkdir(exist_ok=True)
    report = ["straight-sided high-order elements (OpenGrid 11f)", ""]
    report += node_order_report()
    report.append("")
    og = unstructured_cases.hybrid()
    report.append(
        "2D hybrid: structured O-ring (r 1 .. 2, quadrilaterals) + unstructured triangles "
        "(r 2 .. 6); patches p-wall, p-farfield"
    )
    for p in (1, 2, 3, 4):
        t = time.perf_counter()
        mesh = og.high_order(p, curved=False)
        built = time.perf_counter() - t
        counts = mesh.counts()
        report.append(
            f"  order {p}: {counts['quad']} quadrilaterals, {counts['triangle']} triangles, "
            f"{counts['boundary line']} boundary lines, {len(mesh.points)} nodes "
            f"({built * 1e3:.0f} ms); wall nodes off the cylinder by up to "
            f"{wall_distance(mesh, 1.0):.2e}"
        )
    mesh4 = og.high_order(4, curved=False)
    info = {}
    for suffix in ("msh", "cgns", "vtkhdf"):
        info[suffix] = og.export_high_order(
            OUTPUT / f"11f_cylinder_p4.{suffix}", order=4, curved=False
        )
    report.append(
        f"  order 4 written: Gmsh {info['msh']['elements']} elements / {info['msh']['nodes']} "
        f"nodes, CGNS {info['cgns']['sections']} sections in {info['cgns']['zones']} zone, "
        f"VTKHDF {info['vtkhdf']['cells']} Lagrange cells"
    )
    curved = curved_wall_distance(og, 4)
    report.append(
        "  the straight-sided wall is the polygon of the linear mesh at every order: its edges "
        f"are chords, off the cylinder by up to the sagitta 1 - cos(pi / 80) = "
        f"{1 - np.cos(np.pi / 80):.2e}; the curved order-4 elements of the structured ring "
        f"(opt-in, Phase 11e: 4 x 4 grid cells per element) follow it to {curved:.1e}"
    )
    report.append("")
    layers, *_ = layer_cases.body_in_box("structured")
    report.append("3D: body in a box, hexahedral layers on the structured walls")
    first = ", with the linear volume mesh"
    for p in (2, 3, 4):
        t = time.perf_counter()
        mesh = layers.high_order(p, curved=False)
        built = time.perf_counter() - t
        c = mesh.counts()
        report.append(
            f"  order {p}: {c['tetra']} tetrahedra, {c['pyramid']} pyramids, {c['prism']} "
            f"prisms, {c['hexa']} hexahedra, {c['boundary triangle'] + c['boundary quad']} "
            f"boundary faces; {len(mesh.points)} nodes ({built:.2f} s{first})"
        )
        first = ""
    times = {}
    for suffix in ("msh", "cgns", "vtkhdf"):
        t = time.perf_counter()
        layers.export_high_order(OUTPUT / f"11f_layers_p3.{suffix}", order=3, curved=False)
        times[suffix] = time.perf_counter() - t
    report.append(
        "  order 3 written: "
        + ", ".join(f"{k} {v:.1f} s" for k, v in times.items())
        + " (CGNS: TETRA_20, PYRA_30, PENTA_40, HEXA_64 sections; TRI_10 / QUAD_16 boundaries)"
    )
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "11f_high_order_elements_report.txt").write_text(text)
    render(mesh4, OUTPUT / "11f_high_order_elements.png")
    return 0


def render(mesh, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(20, 7), constrained_layout=True)
    ax = fig.add_subplot(1, 3, 1)
    x0, x1, y0, y1 = 1.55, 2.6, -0.45, 0.45
    p = mesh.order
    for kind, color in (("quad", "#305080"), ("triangle", "#c0392b")):
        conn = mesh.cells[kind]
        nv = len(lagrange.VERTICES[kind])
        corners = mesh.points[conn[:, :nv]]
        centre = corners.mean(axis=1)
        inside = (
            (centre[:, 0] > x0) & (centre[:, 0] < x1) & (centre[:, 1] > y0) & (centre[:, 1] < y1)
        )
        for element, nodes in zip(corners[inside], conn[inside], strict=True):
            loop = np.vstack([element, element[:1]])
            ax.plot(loop[:, 0], loop[:, 1], color=color, lw=0.8)
            pts = mesh.points[nodes]
            ax.plot(pts[:, 0], pts[:, 1], ".", color=color, ms=2.5)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_title(f"order {p}: quadrilaterals (ring) and triangles (outside), straight sides")
    colors = {1: "#000000", 2: "#2e86c1", 3: "#27ae60", 4: "#e67e22"}
    for slot, kind in zip((2, 3), ("pyramid", "prism"), strict=True):
        ax3 = fig.add_subplot(1, 3, slot, projection="3d")
        nodes = lagrange.nodes("cgns", kind, 4)
        x = lagrange.as_array(nodes)
        n_vertices = len(lagrange.VERTICES[kind])
        for k, (node, point) in enumerate(zip(nodes, x, strict=True)):
            support = sum(1 for w in lagrange.weights(kind, node) if w)
            level = 4  # interior
            if support <= 2:
                level = support  # vertex, edge
            elif support < n_vertices:
                level = 3  # face
            ax3.scatter(*point, color=colors[level], s=18)
            ax3.text(*point, str(k + 1), fontsize=6)
        spec = lagrange.CGNS[kind]
        corners = x[:n_vertices]
        for a, b in spec["edges"]:
            ax3.plot(*np.stack([corners[a], corners[b]]).T, color="#888888", lw=0.6)
        ax3.set_title(
            f"{lagrange.cgns_type(kind, 4)}: CGNS node numbers (vertex, edge, face, interior)"
        )
        ax3.set_axis_off()
    fig.savefig(path, dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
