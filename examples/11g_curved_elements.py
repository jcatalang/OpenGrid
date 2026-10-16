"""Phase 11g demo: curved boundary elements - high-order nodes on the geometry, the elements
next to it following, every element checked.

  * Cylinder ring: the element edges' distance from the circle, straight-sided against
    curved, for p = 1 .. 5 with 16 and 32 wall edges, and the rates.
  * A coarse cylinder (16 wall edges) with a 1e-5 first cell: the chord sagitta is 2000 cells
    thick, so the boundary layer has to bend with the wall.
  * NACA 0012 C-grid template (1e-5 wall, 24 wall layers): order 4, wall nodes on the airfoil;
    Gmsh and CGNS files.
  * Sphere in a box, unstructured volume with wall layers: hexahedra on structured sphere
    faces, prisms on triangulated ones (with tetrahedra and pyramids); order 3.
  * CAD bracket (STEP through OpenCASCADE, Phase 9a): order 2 in the fluid volume, its nodes on
    the NURBS faces.
Writes 11g_curved_elements.png, 11g_curved_elements_report.txt, 11g_naca_p4.{msh,cgns},
11g_sphere_p3.msh.

    python examples/11g_curved_elements.py
"""

import importlib
import math
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid.model import curving, shape

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
sys.path.insert(0, str(HERE.parent / "tests" / "python"))
tests = importlib.import_module("test_curved_elements")
cases = importlib.import_module("test_templates")


def quality(report) -> str:
    return (
        f"scaled Jacobian of the curved elements {report['min_scaled_jacobian']:.3f} "
        f"(straight-sided {report['min_scaled_jacobian_straight']:.3f})"
    )


def ring_report() -> list[str]:
    out = ["cylinder ring: element edges off the circle, straight-sided -> curved"]
    for p in (1, 2, 3, 4, 5):
        rows = []
        for n in (16, 32):
            og = tests.ring(n)
            rows.append(
                (
                    tests.edge_error(og.high_order(p, curved=False), "p-wall", 1.0),
                    tests.edge_error(og.high_order(p), "p-wall", 1.0),
                )
            )
        rate = math.log2(rows[0][1] / rows[1][1])
        out.append(
            f"  p = {p}: 16 edges {rows[0][0]:.1e} -> {rows[0][1]:.1e}, 32 edges {rows[1][0]:.1e} "
            f"-> {rows[1][1]:.1e} (rate {rate:.1f}; straight-sided 2.0)"
        )
    return out


def main() -> int:
    OUTPUT.mkdir(exist_ok=True)
    report = ["curved boundary elements (OpenGrid 11g)", ""]
    report += ring_report()
    report.append("")
    thin_case = tests.ring(16, first=1e-5, radial=60)
    thin = thin_case.high_order(4)
    flat = thin_case.high_order(4, curved=False)
    r = thin.curving
    wall = curving.boundary_nodes(thin)
    wall = wall[np.linalg.norm(flat.points[wall, :2], axis=1) < 1.5]
    moved = float(np.linalg.norm(thin.points[wall] - flat.points[wall], axis=1).max())
    report.append(
        f"coarse cylinder (16 wall edges), 1e-5 first cell, order 4: wall nodes moved up to "
        f"{moved:.1e} ({moved / 1e-5:.0f} first cells); {r['solved']} nodes followed; "
        f"{quality(r)}"
    )
    report.append("")
    og = opengrid.Session()
    body = cases.airfoil(og)
    og.create_template("cgrid", body=body, wall=1e-5, layers=24, radial_points=(61, 33), name="c")
    og.smooth_group("c-group")
    t = time.perf_counter()
    naca = og.high_order(4)
    built = time.perf_counter() - t
    straight = og.high_order(4, curved=False)
    curves = [og.store.geometry(og.store.get(b).id) for b in body]
    wall = curving.boundary_nodes(naca)

    def off(points, ids):
        return max(min(c.project(tuple(points[i]))[-1] for c in curves) for i in ids)

    on = [i for i in wall if min(c.project(tuple(straight.points[i]))[-1] for c in curves) < 1e-3]
    on = [i for i in on if abs(straight.points[i, 0] - 0.5) < 0.6]
    r = naca.curving
    report.append(
        f"NACA 0012 C-grid, 1e-5 wall, 24 layers, order 4 ({built:.2f} s): {naca.elements} "
        f"quadrilaterals; wall nodes off the airfoil {off(straight.points, on):.1e} "
        f"(straight-sided) -> {off(naca.points, on):.1e}; {r['curved_elements']} curved "
        f"elements; {quality(r)}"
    )
    og.create_patch("airfoil", "wall", body_connectors(og, body))
    og.create_patch(
        "far",
        "farfield",
        [c for c in boundary_connectors(og) if c not in body_connectors(og, body)],
    )
    info_msh = og.export_high_order(OUTPUT / "11g_naca_p4.msh", order=4)
    info_cgns = og.export_high_order(OUTPUT / "11g_naca_p4.cgns", order=4)
    report.append(
        f"  written: Gmsh {info_msh['elements']} elements / {info_msh['nodes']} nodes, CGNS "
        f"{info_cgns['sections']} sections"
    )
    report.append("")
    for kind, label in (("structured", "hexahedral layers"), ("unstructured", "prism layers")):
        sphere, radius = tests.sphere_in_box(kind)
        t = time.perf_counter()
        mesh = sphere.high_order(3)
        built = time.perf_counter() - t
        flat = sphere.high_order(3, curved=False)
        r = mesh.curving
        counts = ", ".join(f"{v} {k}" for k, v in mesh.counts().items() if "boundary" not in k)
        report.append(
            f"sphere in a box, {label}, order 3 ({built:.1f} s): {counts}; sphere faces off the "
            f"sphere {tests.sphere_error(flat, radius):.1e} -> "
            f"{tests.sphere_error(mesh, radius):.1e}; {quality(r)}"
        )
        if kind == "structured":
            sphere.export_high_order(OUTPUT / "11g_sphere_p3.msh", order=3)
    report.append("")
    report += cad_report()
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "11g_curved_elements_report.txt").write_text(text)
    render(straight, naca, thin, OUTPUT / "11g_curved_elements.png")
    return 0


def boundary_connectors(og) -> list[str]:
    return sorted(og.store.get(c).name for c in og.topology().boundary())


def body_connectors(og, body) -> list[str]:
    ids = {og.store.get(b).id for b in body}
    out = []
    for c in boundary_connectors(og):
        segments = og.store.get(c).params["segments"]
        if all(s["curve"] in ids for s in segments):
            out.append(c)
    return out


def cad_report() -> list[str]:
    spec = importlib.util.spec_from_file_location("ex09a", HERE / "09a_cad_import.py")
    ex = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ex)
    if not ex.STEP.exists():
        ex.make_part(ex.STEP)
    og, _walls = ex.fluid_case()
    t = time.perf_counter()
    mesh = og.high_order(2)
    built = time.perf_counter() - t
    flat = og.high_order(2, curved=False)
    assembly = og._linear_assembly()
    worst_flat = worst = 0.0
    faces = 0
    for source in assembly.boundary_sources:
        if source["domain"] is None:
            continue
        entity = og.store.get(source["domain"])
        if entity.params.get("surface") is None:
            continue
        surface = og.store.geometry(entity.params["surface"])
        tris = np.array([c for c in source["cells"] if len(c) == 3], dtype=np.int64)
        if len(tris) == 0:
            continue
        faces += 1
        ids = curving._face_node_ids(mesh, "triangle", tris)[:, 3:].reshape(-1)
        ids = ids[ids >= 0]
        worst_flat = max(worst_flat, max(surface.project(tuple(flat.points[i]))[-1] for i in ids))
        worst = max(worst, max(surface.project(tuple(mesh.points[i]))[-1] for i in ids))
    counts = ", ".join(f"{v} {k}" for k, v in mesh.counts().items() if "boundary" not in k)
    r = mesh.curving
    return [
        f"CAD bracket (STEP), fluid volume, order 2 ({built:.1f} s): {counts}; nodes on its "
        f"{faces} NURBS faces off them {worst_flat:.1e} -> {worst:.1e}; {quality(r)}"
    ]


def render(straight, naca, thin, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(21, 7), constrained_layout=True)
    views = [
        (axes[0], straight, "straight-sided", (-0.012, 0.03, -0.021, 0.021)),
        (axes[1], naca, "curved", (-0.012, 0.03, -0.021, 0.021)),
        (axes[2], thin, "curved, 16 wall edges, 1e-5 first cell", (0.90, 1.04, 0.13, 0.27)),
    ]
    t = np.linspace(0, 1, 25)[:, None]
    values, _ = shape.basis("line", naca.order, t)
    for ax, mesh, title, (x0, x1, y0, y1) in views:
        conn = mesh.cells["quad"]
        p = mesh.order
        nodes = mesh.points[conn]
        inside = (nodes[..., 0] > x0) & (nodes[..., 0] < x1) & (nodes[..., 1] > y0)
        inside = (inside & (nodes[..., 1] < y1)).any(axis=1)
        for element in conn[inside]:
            # Gmsh quad order: vertices, then the nodes of edges 01, 12, 23, 30.
            for k, (a, b) in enumerate(((0, 1), (1, 2), (2, 3), (3, 0))):
                inner = element[4 + k * (p - 1) : 4 + (k + 1) * (p - 1)]
                ids = np.concatenate([[element[a], element[b]], inner])  # a Gmsh line
                curve = values @ mesh.points[ids]
                ax.plot(curve[:, 0], curve[:, 1], color="#305080", lw=0.5)
        if mesh is thin:
            angle = np.linspace(0.0, np.pi / 8, 200)
            ax.plot(np.cos(angle), np.sin(angle), "--", color="#d01c8b", lw=1.0, label="circle")
            ax.plot(
                [1.0, np.cos(np.pi / 8)],
                [0.0, np.sin(np.pi / 8)],
                ":",
                color="#555555",
                label="chord (straight-sided wall)",
            )
            ax.legend(loc="lower left")
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_title(f"order {p}, {title}")
    fig.savefig(path, dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
