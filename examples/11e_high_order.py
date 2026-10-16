"""Phase 11e demo: curved high-order elements for continuous and discontinuous Galerkin.

  * A quarter annulus (one domain, arcs of radius 1 and 2): the distance of the curved element
    edges from the inner circle for p = 1 .. 6 with 4 and 8 elements per side, and the rate.
  * NACA 0012 C-grid (template, 1e-5 wall, 24 wall layers, intervals multiples of 4): order-4
    quadrilaterals - count, validity (scaled Jacobian of the curved elements), Gauss-Lobatto
    nodes snapped onto the airfoil; written as VTKHDF Lagrange cells (equispaced nodes) and
    HDF5 (GLL and Gauss node sets).
  * Pipe (3D template, extruded): order-2 hexahedra, VTKHDF.
Writes 11e_high_order.png, 11e_high_order_report.txt, 11e_naca_p4.vtkhdf, 11e_naca_p4_gll.h5,
11e_pipe_p2.vtkhdf.

    python examples/11e_high_order.py
"""

import importlib
import math
import sys
import time
from pathlib import Path

import numpy as np

import opengrid
from opengrid import _core

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
sys.path.insert(0, str(HERE.parent / "tests" / "python"))
cases = importlib.import_module("test_templates")
annulus_case = importlib.import_module("test_high_order")


def lagrange(nodes, t) -> np.ndarray:
    out = np.ones((len(t), len(nodes)))
    for j in range(len(nodes)):
        for k in range(len(nodes)):
            if k != j:
                out[:, j] *= (t - nodes[k]) / (nodes[j] - nodes[k])
    return out


def edges(mesh, samples=41):
    """Sampled element edges (curves) of a 2D high-order mesh: (edges, samples, 3)."""
    p = mesh.order
    basis = lagrange(_core.nodes_1d(mesh.nodes, p + 1), np.linspace(-1, 1, samples))
    out = []
    for e in mesh.xyz:
        e = e.reshape(p + 1, p + 1, 3)
        for side in (e[0, :], e[-1, :], e[:, 0], e[:, -1]):
            out.append(basis @ side)
    return np.array(out)


def main() -> int:
    OUTPUT.mkdir(exist_ok=True)
    report = ["curved high-order elements from structured grids (OpenGrid 11e)", ""]
    report.append("quarter annulus: distance of the element edges from the inner circle (r = 1)")
    for p in (1, 2, 3, 4, 6):
        ec = annulus_case.arc_error(annulus_case.annulus(4, p).curved_high_order("d", order=p))
        ef = annulus_case.arc_error(annulus_case.annulus(8, p).curved_high_order("d", order=p))
        report.append(
            f"  p = {p}: 4 x 4 elements {ec:.2e}, 8 x 8 elements {ef:.2e} "
            f"(rate {math.log2(ec / ef):.1f}; h^(p+1) would be {p + 1})"
        )
    report.append("")
    og = opengrid.Session()
    t = time.perf_counter()
    og.create_template(
        "cgrid",
        body=cases.airfoil(og),
        wall=1e-5,
        layers=24,
        radial_points=(61, 33),
        interval_multiple=4,
        name="c",
    )
    og.smooth_group("c-group")
    built = time.perf_counter() - t
    t = time.perf_counter()
    mesh = og.curved_high_order(order=4, nodes="gll")
    report.append(
        f"NACA 0012 C-grid, 1e-5 wall, 24 wall layers: grid {built:.1f} s; order 4: "
        f"{mesh.elements} quadrilaterals ({time.perf_counter() - t:.2f} s), scaled Jacobian of "
        f"the curved elements {mesh.scaled_jacobian.min():.3f} .. {mesh.scaled_jacobian.max():.3f}"
    )
    worst = int(np.argmin(mesh.scaled_jacobian))
    centre = mesh.xyz[worst].mean(axis=0)
    unsnapped = og.curved_high_order(order=4, nodes="gll", snap=False)
    move = np.linalg.norm(mesh.xyz - unsnapped.xyz, axis=-1)
    body = (np.abs(mesh.xyz[..., 1]) < 0.2) & (mesh.xyz[..., 0] > -0.05) & (mesh.xyz[..., 0] < 1.05)
    report.append(
        f"  worst element at ({centre[0]:.3f}, {centre[1]:.3f}) (the trailing edge); GLL nodes "
        f"snapped onto the curves: {mesh.snapped}, moved at most {move[body].max():.1e} chords "
        f"on the airfoil and {move.max():.1e} on the far-field arc (the degree-4 mappings' own "
        "distance from the curves)"
    )
    info = og.export_curved_high_order(OUTPUT / "11e_naca_p4.vtkhdf", order=4)
    og.export_curved_high_order(OUTPUT / "11e_naca_p4_gll.h5", order=4, nodes="gll")
    og.export_curved_high_order(OUTPUT / "11e_naca_p4_gauss.h5", order=4, nodes="gauss")
    report.append(
        f"  VTKHDF: {info['cells']} Lagrange quadrilaterals, {info['points']} points; HDF5: "
        "GLL and Gauss node sets"
    )
    report.append("")
    pipe = opengrid.Session()
    pipe.create_template(
        "pipe",
        boundary=[pipe.create_circle((0, 0), 0.5, name="c")],
        length=2.0,
        span_points=9,
        interval_multiple=2,
    )
    pipe.smooth_group(pipe.names("smoothing_group")[0])
    hexes = pipe.curved_high_order(order=2)
    info = pipe.export_curved_high_order(OUTPUT / "11e_pipe_p2.vtkhdf", order=2)
    report.append(
        f"pipe (3D, extruded): order 2: {hexes.elements} hexahedra, scaled Jacobian "
        f"{hexes.scaled_jacobian.min():.3f} .. {hexes.scaled_jacobian.max():.3f}; VTKHDF "
        f"{info['points']} points"
    )
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "11e_high_order_report.txt").write_text(text)
    render(og, mesh, OUTPUT / "11e_high_order.png")
    return 0


def render(og, mesh, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(20, 7), constrained_layout=True)
    curves = edges(mesh)
    views = {
        "leading edge": (-0.02, 0.05, -0.035, 0.035),
        "trailing edge": (0.96, 1.04, -0.04, 0.04),
        "C-grid": (-1.0, 3.0, -1.6, 1.6),
    }
    for ax, (label, (x0, x1, y0, y1)) in zip(axes, views.items(), strict=True):
        for c in curves:
            ax.plot(c[:, 0], c[:, 1], color="#305080", lw=0.5)
        if label != "C-grid":
            pts = mesh.xyz.reshape(-1, 3)
            inside = (pts[:, 0] > x0) & (pts[:, 0] < x1) & (pts[:, 1] > y0) & (pts[:, 1] < y1)
            ax.plot(pts[inside, 0], pts[inside, 1], ".", color="#d01c8b", ms=2)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_title(f"order 4, GLL nodes: {label}")
    fig.savefig(path, dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
