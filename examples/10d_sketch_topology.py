"""Phase 10d demo: sketch topology (a floating wireframe of blocks) around CAD geometry.

  * "airfoil": a NACA 0012 (chord 1) in a C-shaped far field (radius 5, outflow at x = 6),
    modelled with OCCT (netgen.occ), written to output/10d_airfoil.step and imported (curves
    only). A C-topology of 12 blocks (19 vertices: an inner layer around the airfoil and its
    wake, an outer layer to the far field) is drawn by hand, roughly; the curves near its
    edges are assigned (Session.assign_sketch_geometry), the airfoil is a wall and the wake
    cut an internal wall. The vertices on the airfoil, the far field and the outflow slide
    along their curves until the grid stops moving them; the inner ring floats in the group.
    Checks: min scaled Jacobian; the floating vertices moved at random by 20 % of the local
    block size (three trials) and the grid compared, in local spacings; wall nodes on the CAD
    curves; undo, journal replay and .ogp round trip; a ring vertex dragged with the GUI's
    handle logic (offscreen window).
  * "sphere": a sphere (radius 0.5) in a box (side 3), 3D: a cubed-sphere O-grid (6 blocks)
    inside a box layer of 26 blocks (72 vertices, 32 blocks); the sphere is a wall, the box
    faces slide on their planes. The inner cube's vertices float; they are moved at random
    by 20 % of the block size and the grids compared.
  * "wing": a wing stub between two side walls, 3D: the fluid region (the 2D fluid face
    extruded by one chord, STEP), and the C-sketch extruded likewise (38 vertices, 12
    blocks); the CAD edges and faces near the sketch's edges and faces are assigned; upper
    and lower wing faces are walls, the wake cut an internal wall; the vertices on the side
    walls slide on them (those on edges, e.g. where the outflow meets a side wall, stay where
    drawn, projected onto their CAD edge).
Writes 10d_sketch_topology.png, 10d_sketch_topology_report.txt, projects and journals.

    python examples/10d_sketch_topology.py [--quick] [--gui [--screenshot out.png]]
"""

import argparse
import itertools
import logging
import math
import sys
import time
from pathlib import Path

import numpy as np

import opengrid

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
AIRFOIL_STEP = OUTPUT / "10d_airfoil.step"
WING_STEP = OUTPUT / "10d_wing.step"
RADIUS = 5.0  # far field
OUTFLOW = 6.0
SPAN = 1.0
WALL = 1e-3
WAKE = 2e-3


# --------------------------------------------------------------------------- geometry
def naca(x):
    return 0.6 * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)


def airfoil_face():
    import netgen.occ as occ

    beta = np.linspace(0.0, math.pi, 61)
    x = 0.5 * (1.0 - np.cos(beta))
    upper = occ.SplineInterpolation([occ.Pnt(a, naca(a), 0) for a in x])
    lower = occ.SplineInterpolation([occ.Pnt(a, -naca(a), 0) for a in x])
    return occ.Face(occ.Wire([upper, lower]))


def write_cad() -> None:
    """The 2D fluid face (far field minus the airfoil) and the 3D fluid region around the
    wing stub (that face extruded by the span), as STEP."""
    import netgen.occ as occ

    foil = airfoil_face()
    arc = occ.ArcOfCircle(occ.Pnt(1, RADIUS, 0), occ.Pnt(1 - RADIUS, 0, 0), occ.Pnt(1, -RADIUS, 0))
    bottom = occ.Segment(occ.Pnt(1, -RADIUS, 0), occ.Pnt(OUTFLOW, -RADIUS, 0))
    outflow = occ.Segment(occ.Pnt(OUTFLOW, -RADIUS, 0), occ.Pnt(OUTFLOW, RADIUS, 0))
    top = occ.Segment(occ.Pnt(OUTFLOW, RADIUS, 0), occ.Pnt(1, RADIUS, 0))
    fluid = occ.Face(occ.Wire([arc, bottom, outflow, top])) - foil
    fluid.WriteStep(str(AIRFOIL_STEP))
    fluid.Extrude(occ.Vec(0, 0, SPAN)).WriteStep(str(WING_STEP))


# ------------------------------------------------------------------------ the C-sketch
# Along the C: outflow (on the cut), trailing edge, lower middle, leading edge, upper middle,
# trailing edge, outflow; a ring around the airfoil and its wake; the far field.
BASE = (0, 1, 2, 3, 4, 1, 0)
RING = (5, 6, 7, 8, 9, 10, 11)
FAR = (12, 13, 14, 15, 16, 17, 18)
C_BLOCKS = [(BASE[k], BASE[k + 1], RING[k + 1], RING[k]) for k in range(6)] + [
    (RING[k], RING[k + 1], FAR[k + 1], FAR[k]) for k in range(6)
]
# Drawn by hand: roughly around the geometry, nowhere exactly on it.
C_DRAWN = [
    (6.05, 0.04),  # 0 outflow end of the wake cut
    (1.01, 0.0),  # 1 trailing edge
    (0.47, -0.11),  # 2 lower middle
    (-0.03, 0.02),  # 3 leading edge
    (0.55, 0.12),  # 4 upper middle
    (5.9, -0.7),  # 5..11 ring
    (1.1, -0.5),
    (0.4, -0.75),
    (-0.7, 0.1),
    (0.6, 0.8),
    (1.0, 0.65),
    (6.1, 0.55),
    (5.9, -5.1),  # 12..18 far field
    (1.05, -4.9),
    (-2.4, -3.8),
    (-4.1, 0.3),
    (-2.6, 3.3),
    (0.95, 5.05),
    (6.05, 4.9),
]
C_FLOATING = (2, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16)  # slide on curves or float
C_DIMENSIONS = {(0, 1): 41, (1, 2): 31, (2, 3): 31, (3, 4): 31, (4, 1): 31, (1, 6): 31, (6, 13): 31}
C_DISTRIBUTIONS = {
    (3, 2): {"type": "tanh", "begin_spacing": 0.004, "end_spacing": 0.02},
    (3, 4): {"type": "tanh", "begin_spacing": 0.004, "end_spacing": 0.02},
    (1, 2): {"type": "tanh", "begin_spacing": 0.005, "end_spacing": 0.02},
    (1, 4): {"type": "tanh", "begin_spacing": 0.005, "end_spacing": 0.02},
    (1, 0): {"type": "tanh", "begin_spacing": 0.005, "end_spacing": 0.4},
}
C_WALLS = {(1, 2): WALL, (2, 3): WALL, (3, 4): WALL, (4, 1): WALL, (0, 1): WAKE}


def c_edges() -> list[tuple[int, int]]:
    edges = set()
    for block in C_BLOCKS:
        for k in range(4):
            a, b = block[k], block[(k + 1) % 4]
            edges.add((min(a, b), max(a, b)))
    return sorted(edges)


def block_sizes(points: np.ndarray, edges) -> np.ndarray:
    """Local block size per vertex: its shortest sketch edge."""
    size = np.full(len(points), np.inf)
    for a, b in edges:
        length = float(np.linalg.norm(points[a] - points[b]))
        size[a] = min(size[a], length)
        size[b] = min(size[b], length)
    return size


def local_spacing(grid: np.ndarray) -> np.ndarray:
    """Shortest grid edge at every node (any dimension)."""
    shape = grid.shape[:-1]
    out = np.full(shape, np.inf)
    for axis in range(len(shape)):
        h = np.linalg.norm(np.diff(grid, axis=axis), axis=-1)
        lo = [slice(None)] * len(shape)
        hi = [slice(None)] * len(shape)
        lo[axis] = slice(0, -1)
        hi[axis] = slice(1, None)
        out[tuple(lo)] = np.minimum(out[tuple(lo)], h)
        out[tuple(hi)] = np.minimum(out[tuple(hi)], h)
    return out


def change(a: dict, b: dict) -> float:
    """Largest node displacement between two grid sets, in local spacings."""
    worst = 0.0
    for key, grid in a.items():
        moved = np.linalg.norm(b[key] - grid, axis=-1) / local_spacing(grid)
        worst = max(worst, float(moved.max()))
    return worst


class PassCounter(logging.Handler):
    def __init__(self):
        super().__init__()
        self.lines = []

    def emit(self, record):
        message = record.getMessage()
        if "vertex passes" in message:
            self.lines.append(message)


def airfoil_session(drawn, curves=None) -> opengrid.Session:
    """The airfoil sketch; ``curves``: the assignment to use (default: the CAD curves near
    the edges, Session.assign_sketch_geometry, as the GUI's "assign selected geometry")."""
    og = opengrid.Session()
    og.import_cad(AIRFOIL_STEP, name="cad", domains="none")
    og.delete(og.names("connector"))  # the sketch makes its own
    # Set up with the vertex passes off; the last edit turns them on (one relaxation).
    og.create_sketch(
        drawn,
        C_BLOCKS,
        dimensions=C_DIMENSIONS,
        distributions=C_DISTRIBUTIONS,
        relax=False,
        name="c",
    )
    if curves is None:
        og.assign_sketch_geometry("c", og.names("nurbs"))
    else:
        og.edit_sketch("c", curves=curves)
    og.edit_sketch("c", walls=C_WALLS, relax=True)
    return og


def named_curves(og, sketch: str) -> dict:
    return {
        item["edge"]: og.store.get(item["curve"]).name
        for item in og.entity(sketch).params["curves"]
    }


def domain_grids(og) -> dict:
    return {d: og.domain_grid(d).copy() for d in og.names("domain")}


def min_jacobian(og, names, kind="domain") -> float:
    if kind == "domain":
        return min(float(og.domain_quality(n)["scaled_jacobian"].min()) for n in names)
    return min(float(og.block_quality(n)["scaled_jacobian"].min()) for n in names)


def airfoil_case(report: list[str], counter: PassCounter, quick: bool):
    t = time.perf_counter()
    og = airfoil_session(C_DRAWN)
    elapsed = time.perf_counter() - t
    t = time.perf_counter()
    og.smooth_group("c-group")
    smooth = time.perf_counter() - t
    grids = domain_grids(og)
    domains = og.names("domain")
    assigned = og.entity("c").params["curves"]
    report += [
        f"airfoil: C-topology sketch, {len(C_DRAWN)} vertices, {len(C_BLOCKS)} blocks, "
        f"{len(assigned)} edges assigned to CAD curves by assign_sketch_geometry",
        f"  generated in {elapsed:.1f} s ({counter.lines[-1]}), group smoothed in {smooth:.2f} s",
        f"  {og.group_status('c-group')}",
        f"  min scaled Jacobian {min_jacobian(og, domains):.3f} (target >= 0.5)",
    ]
    report.append(
        "  (the kink is at the trailing edge, where the C-lines turn by the airfoil's 14 deg"
        " wedge; it decays to 8 deg six cells away)"
    )
    # Wall nodes on the CAD curves.
    worst = 0.0
    for item in assigned:
        curve = og.store.geometry(item["curve"])
        a, b = item["edge"]
        nodes = og.display_connector_nodes(f"c-c{a}-{b}").points
        worst = max(worst, max(curve.project(tuple(p))[2] for p in nodes))
    report.append(f"  nodes on assigned CAD curves within {worst:.1e} (chord 1)")
    # The floating vertices moved by 20 % of the local block size.
    edges = c_edges()
    drawn = np.array([(*p, 0.0) for p in C_DRAWN])
    size = block_sizes(drawn, edges)
    rng = np.random.default_rng(10)
    trials = 3
    if quick:
        trials = 1
    changes = []
    for _trial in range(trials):
        moved = drawn.copy()
        for v in C_FLOATING:
            angle = rng.uniform(0.0, 2.0 * math.pi)
            moved[v, :2] += 0.2 * size[v] * np.array([math.cos(angle), math.sin(angle)])
        # The same assignment (moved vertices could bring other curves within reach).
        other = airfoil_session([tuple(p) for p in moved], named_curves(og, "c"))
        other.smooth_group("c-group")
        changes.append(change(grids, domain_grids(other)))
    report.append(
        f"  floating vertices moved by 20 % of the block size ({trials} trials): the grid "
        f"changes by at most {max(changes):.1e} local spacings (target < 1e-2)"
    )
    # Undo, journal, project.
    og.save(OUTPUT / "10d_airfoil.ogp")
    og.save_journal(OUTPUT / "10d_airfoil_journal.py")
    replayed = opengrid.replay(og.journal_text())
    replayed.smooth_group("c-group")
    loaded = opengrid.Session()
    loaded.load(OUTPUT / "10d_airfoil.ogp")
    loaded.smooth_group("c-group")
    report.append(
        f"  journal replay: grid identical {change(grids, domain_grids(replayed)) == 0.0}; "
        f".ogp round trip: identical {change(grids, domain_grids(loaded)) == 0.0}"
    )
    before = len(og.names())
    og.undo()  # the walls
    walls = len(og.entity("c").params["walls"])
    og.undo()  # the assignment
    curves = len(og.entity("c").params["curves"])
    og.redo()
    og.redo()
    report.append(
        f"  undo: one step per edit (after one undo {walls} walls, after two {curves} assigned "
        f"curves; redo restores both; {before} entities)"
    )
    try:
        og.rename("c-c3-8", "mine")
    except opengrid.validation.ValidationError as error:
        report.append(f"  renaming an output directly is refused: {error.issues[0].message}")
    return og, grids


def drag_ring_vertex(og, report: list[str]) -> dict:
    """Drag ring vertex 9 with the GUI's handle logic: previews, one command, then compare."""
    from PySide6.QtWidgets import QApplication

    from opengrid.app.main_window import MainWindow

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = MainWindow(off_screen=True, session=og)
    window.set_selection([og.entity("c").id])
    target = next(h for h in window.visible_handles() if h.entity == "c" and h.index == 9)
    assert window.handle_press(target.xyz)
    times = []
    for k in range(1, 11):
        xyz = (target.xyz[0] + 0.02 * k, target.xyz[1] - 0.015 * k, 0.0)
        t = time.perf_counter()
        window.handle_drag(xyz)
        times.append(time.perf_counter() - t)
    t = time.perf_counter()
    window.handle_release((target.xyz[0] + 0.2, target.xyz[1] - 0.15, 0.0))
    release = time.perf_counter() - t
    window.smoother_dock.wait()
    app.processEvents()
    window.close()
    og.smooth_group("c-group")
    report.append(
        f"  ring vertex 9 dragged by (0.2, -0.15): {len(times)} previews at "
        f"{1000 * np.mean(times):.0f} ms ({1 / np.mean(times):.0f} frames/s), release "
        f"(regeneration with the vertex passes) {release:.1f} s; journal: "
        f"{og.journal_text().splitlines()[-1][:58]}..."
    )
    return domain_grids(og)


# --------------------------------------------------------------------- sphere in a box
SPHERE = 0.5
CUBE = 0.9
BOX = 1.5


def sphere_geometry(og) -> dict:
    og.create_arc(
        (0, 0, 0),
        SPHERE,
        -math.pi / 2,
        math.pi / 2,
        normal=(0, -1, 0),
        x_axis=(1, 0, 0),
        name="mer_z",
    )
    og.create_arc(
        (0, 0, 0),
        SPHERE,
        -math.pi / 2,
        math.pi / 2,
        normal=(0, 0, -1),
        x_axis=(0, 1, 0),
        name="mer_x",
    )
    og.create_revolved_surface(
        "mer_z", 2.5 * math.pi, axis=(0, 0, 1), start_angle=0.5 * math.pi, name="sphere_x"
    )
    og.create_revolved_surface("mer_z", 2 * math.pi, axis=(0, 0, 1), name="sphere_y")
    og.create_revolved_surface("mer_x", 2 * math.pi, axis=(1, 0, 0), name="sphere_z")
    return {0: "sphere_x", 1: "sphere_y", 2: "sphere_z"}


def sphere_session(jitter=None, n=13, n_r=17, n_o=9):
    """Box layer: a 4 x 4 x 4 lattice of vertices without its centre cell; inside it six
    shells from the cubed sphere's patches to the inner cube. ``jitter``: (rng, size) moves
    the inner cube's vertices."""
    og = opengrid.Session()
    surfaces = sphere_geometry(og)
    coords = (-BOX, -CUBE, CUBE, BOX)
    vertices = []
    index = {}
    for k, j, i in itertools.product(range(4), repeat=3):
        point = np.array([coords[i], coords[j], coords[k]])
        if jitter is not None and 0 < i < 3 and 0 < j < 3 and 0 < k < 3:
            rng, size = jitter
            direction = rng.normal(size=3)
            point = point + 0.2 * size * direction / np.linalg.norm(direction)
        index[(i, j, k)] = len(vertices)
        vertices.append(tuple(point.tolist()))
    on_sphere = {}
    for c in itertools.product((-1, 1), repeat=3):
        on_sphere[c] = len(vertices)
        vertices.append(tuple((SPHERE * np.array(c) / math.sqrt(3)).tolist()))
    corners = (
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    )
    blocks = []
    for k, j, i in itertools.product(range(3), repeat=3):
        if (i, j, k) != (1, 1, 1):
            blocks.append(tuple(index[(i + a, j + b, k + c)] for a, b, c in corners))
    faces = {}
    arcs = {}
    for axis in range(3):
        for side in (-1, 1):
            a, b = [q for q in range(3) if q != axis]
            cs = []
            for u, v in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                c = [0, 0, 0]
                c[axis], c[a], c[b] = side, u, v
                cs.append(tuple(c))
            inner = [index[tuple(1 + (s > 0) for s in c)] for c in cs]
            blocks.append(tuple(on_sphere[c] for c in cs) + tuple(inner))
            faces[tuple(sorted(on_sphere[c] for c in cs))] = surfaces[axis]
            for m in range(4):
                edge = tuple(sorted((on_sphere[cs[m]], on_sphere[cs[(m + 1) % 4]])))
                if edge not in arcs:
                    pa = np.array(vertices[edge[0]])
                    pb = np.array(vertices[edge[1]])
                    mid = SPHERE * (pa + pb) / np.linalg.norm(pa + pb)
                    arcs[edge] = og.create_arc_3pt(
                        tuple(pa), tuple(mid), tuple(pb), name=f"arc{len(arcs)}"
                    )
    dimensions = {next(iter(arcs)): n, (blocks[-1][0], blocks[-1][4]): n_r}
    for axis in range(3):
        for segment, count in ((0, n_o), (1, n), (2, n_o)):
            lo = [0, 0, 0]
            hi = [0, 0, 0]
            lo[axis] = segment
            hi[axis] = segment + 1
            dimensions[(index[tuple(lo)], index[tuple(hi)])] = count
    og.create_sketch(
        vertices,
        blocks,
        dimensions=dimensions,
        curves=arcs,
        surfaces=faces,
        walls=dict.fromkeys(faces, 0.01),
        name="box",
    )
    inner = [index[(i, j, k)] for k, j, i in itertools.product((1, 2), repeat=3)]
    return og, inner


def block_grids(og) -> dict:
    return {b: og.block_grid(b).copy() for b in og.names("block")}


def sphere_case(report: list[str], quick: bool):
    t = time.perf_counter()
    og, inner = sphere_session()
    elapsed = time.perf_counter() - t
    t = time.perf_counter()
    og.smooth_group("box-group")
    smooth = time.perf_counter() - t
    blocks = og.names("block")
    grids = block_grids(og)
    points = np.concatenate([g.reshape(-1, 3) for g in grids.values()])
    radius = np.linalg.norm(points, axis=1)
    near = radius < SPHERE * (1 + 1e-6)
    report += [
        "",
        f"sphere in a box: {len(og.entity('box').params['vertices'])} vertices, {len(blocks)} "
        f"blocks (cubed-sphere O-grid in a 26-block box layer), generated in {elapsed:.1f} s",
        f"  {og.group_status('box-group')}, {smooth:.1f} s",
        f"  min scaled Jacobian {min_jacobian(og, blocks, 'block'):.3f} (target >= 0.5)",
        f"  {int(near.sum())} nodes on the sphere within {np.abs(radius[near] - SPHERE).max():.1e}"
        " (wall faces projected onto their surface)",
    ]
    rng = np.random.default_rng(3)
    trials = 2
    if quick:
        trials = 1
    changes = []
    for _trial in range(trials):
        other, _inner = sphere_session(jitter=(rng, CUBE - SPHERE / math.sqrt(3)))
        other.smooth_group("box-group")
        changes.append(change(grids, block_grids(other)))
    report.append(
        f"  inner cube vertices moved by 20 % of the block size ({trials} trials): the grid "
        f"changes by at most {max(changes):.1e} local spacings"
    )
    og.save(OUTPUT / "10d_sphere.ogp")
    og.save_journal(OUTPUT / "10d_sphere_journal.py")
    return og, inner


# ---------------------------------------------------------------------------- wing stub
def wing_session():
    og = opengrid.Session()
    result = og.import_cad(WING_STEP, name="wing", domains="auto")
    # Only the CAD curves and surfaces are kept: the sketch makes the grid.
    og.delete(result.domains + og.names("connector"))
    n = len(C_DRAWN)
    vertices = [(*p, 0.0) for p in C_DRAWN] + [(*p, SPAN) for p in C_DRAWN]
    blocks = [(*b, *(v + n for v in b)) for b in C_BLOCKS]
    dimensions = dict(C_DIMENSIONS)
    dimensions[(0, n)] = 13
    distributions = dict(C_DISTRIBUTIONS)
    distributions.update({(a + n, b + n): d for (a, b), d in C_DISTRIBUTIONS.items()})
    og.create_sketch(
        vertices,
        blocks,
        dimensions=dimensions,
        distributions=distributions,
        relax=False,
        name="w",
    )
    geometry = og.names("nurbs") + og.names("nurbs_surface") + og.names("line")
    og.assign_sketch_geometry("w", geometry)
    walls = {}
    for a, b in ((1, 2), (2, 3), (3, 4), (4, 1)):
        walls[(a, b, b + n, a + n)] = WALL
    walls[(0, 1, 1 + n, n)] = WAKE  # the wake cut: an internal wall
    og.edit_sketch("w", walls=walls, relax=True)
    return og


WAKE_BLOCKS = ("w-b0", "w-b5", "w-b6", "w-b11")


def wing_case(report: list[str]):
    t = time.perf_counter()
    og = wing_session()
    elapsed = time.perf_counter() - t
    t = time.perf_counter()
    og.smooth_group("w-group")
    smooth = time.perf_counter() - t
    blocks = og.names("block")
    params = og.entity("w").params
    report += [
        "",
        f"wing stub: the C-sketch extruded by one chord ({len(params['vertices'])} vertices, "
        f"{len(blocks)} blocks), {len(params['curves'])} edges and {len(params['surfaces'])} "
        f"faces assigned to the CAD model, generated in {elapsed:.1f} s",
        f"  {og.group_status('w-group')}, {smooth:.1f} s",
        f"  min scaled Jacobian {min_jacobian(og, blocks, 'block'):.3f}; outside the four wake"
        f" blocks {min_jacobian(og, [b for b in blocks if b not in WAKE_BLOCKS], 'block'):.3f}",
        "  (the wake blocks' vertices on the outflow edges stay where they were drawn - 3D"
        " vertices on curves do not slide - and the cells between the drawn ring and the"
        " internal wall are sheared there)",
    ]
    og.save(OUTPUT / "10d_wing.ogp")
    og.save_journal(OUTPUT / "10d_wing_journal.py")
    return og


# ------------------------------------------------------------------------------ picture
def render(airfoil, sphere, wing, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    fig, axes = plt.subplots(2, 3, figsize=(19, 11), constrained_layout=True)
    og, grids = airfoil
    sketch = og.entity("c").params
    drawn = np.asarray(sketch["vertices"])

    def wireframe(ax, points, edges):
        for a, b in edges:
            ax.plot(points[[a, b], 0], points[[a, b], 1], "-", color="#d01c8b", lw=1.6)
        ax.plot(points[:, 0], points[:, 1], "o", color="#d01c8b", ms=4)

    def geometry(ax):
        for name in og.names("nurbs"):
            p = np.asarray(og.geometry(name).tessellate_points(1e-4))
            ax.plot(p[:, 0], p[:, 1], color="#202020", lw=1.4)

    def grid_lines(ax, grid_map, lw=0.35):
        for grid in grid_map.values():
            ax.plot(grid[..., 0], grid[..., 1], color="#305080", lw=lw)
            ax.plot(grid[..., 0].T, grid[..., 1].T, color="#305080", lw=lw)

    ax = axes[0, 0]
    geometry(ax)
    wireframe(ax, drawn, c_edges())
    ax.set_title("airfoil: CAD curves and the sketch as drawn")
    ax.set_xlim(-5.5, 6.8)
    ax.set_ylim(-5.6, 5.6)
    ax = axes[0, 1]
    grid_lines(ax, grids)
    ax.set_title("generated and smoothed (C-topology, 12 blocks)")
    ax.set_xlim(-5.5, 6.8)
    ax.set_ylim(-5.6, 5.6)
    ax = axes[0, 2]
    grid_lines(ax, grids, 0.3)
    geometry(ax)
    wireframe(ax, drawn, c_edges())
    ax.set_title("near the airfoil (sketch in magenta)")
    ax.set_xlim(-0.35, 1.45)
    ax.set_ylim(-0.6, 0.6)

    def cut(ax, session, title, level=0.0, axis=2):
        """Cells of all blocks crossing the plane x[axis] = level, coloured by quality."""
        polys = []
        values = []
        others = [c for c in range(3) if c != axis]
        for name in session.names("block"):
            grid = session.block_grid(name)
            quality = session.block_quality(name)["scaled_jacobian"]
            nk, nj, ni = grid.shape[:3]
            coord = grid[..., axis]
            for k, j, i in itertools.product(range(nk - 1), range(nj - 1), range(ni - 1)):
                values_here = coord[k : k + 2, j : j + 2, i : i + 2]
                if not values_here.min() <= level < values_here.max():
                    continue
                polygon = _section(grid[k : k + 2, j : j + 2, i : i + 2], axis, level, others)
                if polygon is not None:
                    polys.append(polygon)
                    values.append(quality[k, j, i])
        collection = PolyCollection(
            polys, array=np.array(values), cmap="viridis", edgecolors="k", linewidths=0.15
        )
        collection.set_clim(0.0, 1.0)
        ax.add_collection(collection)
        ax.autoscale_view()
        ax.set_title(title)
        return collection

    # Just off z = 0, a plane of nodes (cells meet it in faces there).
    collection = cut(axes[1, 0], sphere, "sphere in a box: cut z = 0 (scaled Jacobian)", level=1e-6)
    axes[1, 0].set_xlim(-1.6, 1.6)
    axes[1, 0].set_ylim(-1.6, 1.6)
    cut(axes[1, 1], wing, "wing stub: cut at mid-span (scaled Jacobian)", level=0.5 * SPAN)
    axes[1, 1].set_xlim(-0.6, 1.8)
    axes[1, 1].set_ylim(-0.9, 0.9)
    cut(axes[1, 2], wing, "wing stub: the whole cut at mid-span", level=0.5 * SPAN)
    fig.colorbar(collection, ax=axes[1, :], shrink=0.6, label="scaled Jacobian")
    for ax in axes.ravel():
        ax.set_aspect("equal")
    fig.savefig(path, dpi=100)
    plt.close(fig)


def _section(cell: np.ndarray, axis: int, level: float, others) -> np.ndarray | None:
    """Polygon where a hexahedral cell (2, 2, 2, 3) meets the plane x[axis] = level (its
    four edges along the most aligned direction, interpolated)."""
    best = None
    for direction in range(3):
        a = np.take(cell, 0, axis=direction).reshape(-1, 3)
        b = np.take(cell, 1, axis=direction).reshape(-1, 3)
        span = np.abs(b[:, axis] - a[:, axis])
        if best is None or span.min() > best[0]:
            best = (span.min(), a, b)
    _span, a, b = best
    if _span <= 0.0:
        return None
    t = (level - a[:, axis]) / (b[:, axis] - a[:, axis])
    if np.any(t < -1e-9) or np.any(t > 1 + 1e-9):
        return None
    p = a + t[:, None] * (b - a)
    return p[[0, 1, 3, 2]][:, others]


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
    window.set_selection([window.session.entity("c").id])
    window.quality_action.setChecked(True)
    window.run_smoother("c-group")
    if screenshot is not None:

        def save_and_quit():
            window.smoother_dock.wait()
            window.set_selection([window.session.entity("c").id])
            app.processEvents()
            window.viewport.render()
            window.grab().save(str(screenshot))
            print(f"saved {screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(4000, save_and_quit)
    return app.exec()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    counter = PassCounter()
    logging.getLogger("opengrid").addHandler(counter)
    logging.getLogger("opengrid").setLevel(logging.INFO)
    write_cad()
    if args.gui:
        airfoil_session(C_DRAWN).save(OUTPUT / "10d_airfoil.ogp")
        return run_gui(OUTPUT / "10d_airfoil.ogp", args.screenshot)
    report: list[str] = []
    airfoil = airfoil_case(report, counter, args.quick)
    dragged = drag_ring_vertex(airfoil[0], report)
    report.append(
        f"    after the drag the grid differs from before by {change(airfoil[1], dragged):.1e} "
        "local spacings (the ring floats)"
    )
    sphere, _inner = sphere_case(report, args.quick)
    wing = wing_case(report)
    text = "\n".join(report) + "\n"
    print(text, end="")
    (OUTPUT / "10d_sketch_topology_report.txt").write_text(text)
    render(airfoil, sphere, wing, OUTPUT / "10d_sketch_topology.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
