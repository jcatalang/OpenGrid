"""Automatic 2D block topology from a cross field (frame field): a planar region bounded by
closed loops of curves is cut into four-sided blocks, which become a 10d sketch.

Method:
1. Triangulate the region (the 5a Delaunay mesher, in a scratch session).
2. Cross field: the complex value z = exp(4 i theta) is fixed on the boundary nodes (theta the
   boundary tangent; at corners the mean of the two legs' values) and interpolated
   harmonically inside (graph Laplacian, conjugate gradients). A cross is the four directions
   arg(z) / 4 + k pi / 2.
3. Singularities: triangles around which arg(z) winds by 2 pi k, k != 0 (index k / 4); a
   boundary corner of interior angle n pi / 2 counts (2 - n) / 4. Poincare-Hopf for cross
   fields: the indices add up to the Euler characteristic 1 - holes (checked and reported).
4. Separatrices: from each singularity, the directions where a cross direction points
   radially (4 - 4 k of them: 3 for index +1/4, 5 for -1/4), and from each boundary corner or
   curve junction the n - 1 cross directions into the region; each is traced as a streamline
   of the cross field (midpoint rule, the cross direction closest to the previous one) until
   it leaves the region (ended exactly on the boundary curve there) or reaches another
   singularity (ended at it; that singularity's separatrix in that direction is the same).
5. Blocks: the boundary curves and the separatrices (as polylines) form a planar network
   (opengrid.model.network: crossings, faces); every face has four corners, because each
   separatrix runs on to the boundary. The faces become the blocks of a sketch whose edges
   are assigned to the pieces of the curves and separatrices between the vertices.
References: N. Kowalski, F. Ledoux, P. Frey, "A PDE based approach to multi-domain
partitioning and quadrilateral meshing", IMR 2012; D. Bommes et al., "Quad-mesh generation
and processing: a survey", Computer Graphics Forum 32 (2013) 51-76; M. Reberol et al.,
"Quasi-structured quadrilateral meshing in Gmsh" (2021) for the boundary-corner counts.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from opengrid.model.connectors import CURVE_KINDS
from opengrid.validation import Issue, ValidationError

CORNER_ANGLE = 30.0  # degrees of boundary turning that make a corner


@dataclass
class Singularity:
    point: np.ndarray
    index: float  # +-1/4 (or multiples)
    triangle: int


@dataclass
class Corner:
    point: np.ndarray
    node: int
    quarters: int  # n: interior angle n pi / 2
    curve_junction: bool


@dataclass
class Field:
    points: np.ndarray  # (n, 3)
    triangles: np.ndarray  # (m, 3)
    z: np.ndarray  # complex, per node
    boundary: list  # loops of node ids (region on the left)
    corners: list[Corner]
    singularities: list[Singularity]
    holes: int
    size: float  # mesh size

    def index_sum(self) -> float:
        return sum(s.index for s in self.singularities) + sum(
            (2 - c.quarters) / 4.0 for c in self.corners
        )

    def euler(self) -> int:
        return 1 - self.holes


@dataclass
class Layout:
    field: Field
    separatrices: list  # polylines (k, 3)
    warnings: list = field(default_factory=list)
    marks: list = field(default_factory=list)  # boundary points that are corners of faces


def _fail(message: str):
    raise ValidationError([Issue("automatic topology", message)])


# ------------------------------------------------------------------------- triangulation
def _triangulate(store, curve_ids, size: float):
    """Points, triangles and boundary loops of the region (a scratch session's 5a mesh)."""
    import opengrid

    scratch = opengrid.Session()
    names = []
    for c in curve_ids:
        entity = store.get(c)
        scratch._create(entity.kind, dict(entity.params), entity.name)
        positions = _boundary_positions(store.geometry(c), size)
        names.append(
            scratch.create_connector(
                entity.name,
                len(positions),
                {"type": "explicit", "positions": positions},
                name=f"c-{entity.name}",
            )
        )
    domain = scratch.create_unstructured_domain(connectors=names, name="region")
    mesh = scratch.unstructured_mesh(domain)
    points = np.asarray(mesh["points"], dtype=float)
    triangles = np.asarray(mesh["triangles"], dtype=np.int64)
    return points, triangles


def _boundary_positions(curve, size: float, max_turn: float = 3.0) -> tuple:
    """Node positions (fractions of the arc length) along a boundary curve: a new node where
    the tangent has turned by ``max_turn`` degrees or the length reached ``size`` - so the
    field's 4 theta turns by well under half a turn between neighbours, even at a leading
    edge."""
    t = np.linspace(0.0, 1.0, 20001)
    p = np.asarray(curve.eval_many(t))[:, :2]
    d = np.diff(p, axis=0)
    seg = np.linalg.norm(d, axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    a = np.arctan2(d[:, 1], d[:, 0])
    turn = np.concatenate([[0.0], np.abs(_wrap(np.diff(a)))])
    out = [0.0]
    last_s = 0.0
    acc = 0.0
    for k in range(1, len(seg)):
        acc += math.degrees(float(turn[k]))
        if acc >= max_turn or s[k] - last_s >= size:
            out.append(float(s[k] / s[-1]))
            last_s = s[k]
            acc = 0.0
    if 1.0 - out[-1] < 0.25 * size / s[-1] and len(out) > 1:
        out.pop()
    out.append(1.0)
    if len(out) < 4:
        out = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
    return tuple(out)


def _boundary_loops(triangles: np.ndarray) -> list[list[int]]:
    """Boundary loops (each with the region on its left)."""
    edges: dict = {}
    for t in triangles:
        for k in range(3):
            a, b = int(t[k]), int(t[(k + 1) % 3])
            edges[(a, b)] = edges.get((a, b), 0) + 1
    directed = {}
    for a, b in edges:
        if (b, a) not in edges:
            directed[a] = b  # oriented as in a counter-clockwise triangle: region on the left
    loops = []
    seen = set()
    for start in sorted(directed):
        if start in seen:
            continue
        loop = [start]
        seen.add(start)
        current = directed[start]
        while current != start:
            loop.append(current)
            seen.add(current)
            current = directed[current]
        loops.append(loop)
    return loops


def _orient(points, triangles) -> np.ndarray:
    p = points[triangles]
    area = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (p[:, 2, 0] - p[:, 0, 0]) * (
        p[:, 1, 1] - p[:, 0, 1]
    )
    out = triangles.copy()
    flip = area < 0
    out[flip] = out[flip][:, [0, 2, 1]]
    return out


# --------------------------------------------------------------------------- the field
def _laplace(n, triangles, fixed_mask, values) -> np.ndarray:
    """Harmonic interpolation (graph Laplacian) of complex ``values`` from the fixed nodes."""
    edges = set()
    for t in triangles:
        for k in range(3):
            a, b = int(t[k]), int(t[(k + 1) % 3])
            edges.add((min(a, b), max(a, b)))
    edges = np.array(sorted(edges), dtype=np.int64)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    degree = np.bincount(rows, minlength=n).astype(float)
    free = ~fixed_mask
    x = values.astype(complex).copy()
    x[free] = 0.0

    def apply(v):  # (L v) on all nodes, L = D - A
        out = degree * v
        np.subtract.at(out, rows, v[cols])
        return out

    # Solve L_ff x_f = -L_fc x_c by conjugate gradients (real and imaginary parts together).
    rhs = -apply(np.where(free, 0.0, x))
    b = np.where(free, rhs, 0.0)
    xf = np.zeros(n, dtype=complex)
    r = b - np.where(free, apply(xf), 0.0)
    p = r.copy()
    rr = np.vdot(r, r).real
    norm_b = math.sqrt(np.vdot(b, b).real) or 1.0
    for _iteration in range(20 * n):
        if math.sqrt(rr) <= 1e-12 * norm_b:
            break
        ap = np.where(free, apply(p), 0.0)
        alpha = rr / np.vdot(p, ap).real
        xf += alpha * p
        r -= alpha * ap
        new = np.vdot(r, r).real
        p = r + (new / rr) * p
        rr = new
    x[free] = xf[free]
    return x


def _smooth_closed(curve) -> bool:
    """A closed curve whose tangent is continuous at its seam."""
    a = np.asarray(curve.eval(0.0))
    b = np.asarray(curve.eval(1.0))
    if np.linalg.norm(a - b) > 1e-9 * max(1.0, float(np.max(np.abs(a)))):
        return False
    step = 1e-6
    u = np.asarray(curve.eval(step)) - a
    v = b - np.asarray(curve.eval(1.0 - step))
    cos = float(u @ v) / (float(np.linalg.norm(u) * np.linalg.norm(v)) or 1.0)
    return cos > math.cos(math.radians(CORNER_ANGLE))


def _mesh_directions(prev, here, nxt):
    u = here - prev
    v = nxt - here
    return math.atan2(u[1], u[0]), math.atan2(v[1], v[0])


def _curve_directions(curves, prev, here, nxt, a_in, a_out):
    """Directions of the boundary arriving at and leaving a curve junction, from the curves'
    tangents: the curve holding the previous (next) boundary node, near its end here."""

    def tangent(neighbour, leaving):
        # The curve end at this point whose inside lies towards the neighbouring boundary
        # node (a closed curve has both of its ends here).
        scale = max(1.0, float(np.max(np.abs(here))))
        best = None
        step = 1e-6
        for curve in curves.values():
            for t_end, inside in ((0.0, step), (1.0, 1.0 - step)):
                if np.linalg.norm(np.asarray(curve.eval(t_end)) - here) > 1e-9 * scale:
                    continue
                q = np.asarray(curve.eval(inside))
                _t, _q, d = curve.project(tuple(float(c) for c in neighbour))
                score = (d, float(np.linalg.norm(q - neighbour)))
                if best is None or score < best[0]:
                    best = (score, q)
        if best is None:
            return None
        q = best[1]
        if leaving:
            v = q - here
        else:
            v = here - q
        return math.atan2(v[1], v[0])

    incoming = tangent(prev, False)
    outgoing = tangent(nxt, True)
    if incoming is None:
        incoming = a_in
    if outgoing is None:
        outgoing = a_out
    return incoming, outgoing


def _curve_tangent(curves, point, along) -> float:
    """Direction of the nearest curve's tangent at ``point``, oriented like ``along``."""
    best = None
    for curve in curves.values():
        t, _q, d = curve.project(tuple(float(c) for c in point))
        if best is None or d < best[0]:
            best = (d, curve, t)
    _d, curve, t = best
    step = 1e-6
    a = np.asarray(curve.eval(max(t - step, 0.0)))
    b = np.asarray(curve.eval(min(t + step, 1.0)))
    v = b - a
    if float(v[:2] @ along[:2]) < 0.0:
        v = -v
    return math.atan2(v[1], v[0])


NORMALIZED_SWEEPS = 400


def _normalize_smooth(n, triangles, fixed, z, sweeps: int):
    """Unit-length smoothing of the harmonic field: each free value becomes the normalised
    mean of its neighbours, ``sweeps`` times (Jacobi). Minimises a Ginzburg-Landau-type
    energy, which keeps the singularities off the boundary (where the harmonic solution
    leaves them next to strongly curved parts)."""
    edges = set()
    for t in triangles:
        for k in range(3):
            a, b = int(t[k]), int(t[(k + 1) % 3])
            edges.add((min(a, b), max(a, b)))
    edges = np.array(sorted(edges), dtype=np.int64)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    free = ~fixed
    out = z.copy()
    for _sweep in range(sweeps):
        total = np.zeros(n, dtype=complex)
        np.add.at(total, rows, out[cols])
        magnitude = np.abs(total)
        ok = free & (magnitude > 0)
        out[ok] = total[ok] / magnitude[ok]
    return out


def _wrap(a):
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def compute_field(store, curve_ids, size: float | None = None) -> Field:
    """The cross field of the region bounded by the curves, with its singularities."""
    lengths = [store.geometry(c).length() for c in curve_ids]
    if size is None:
        size = sum(lengths) / 200.0
    points, triangles = _triangulate(store, curve_ids, size)
    triangles = _orient(points, triangles)
    loops = _boundary_loops(triangles)
    n = len(points)
    values = np.zeros(n, dtype=complex)
    fixed = np.zeros(n, dtype=bool)
    corners = []
    junction_points = []
    for c in curve_ids:
        curve = store.geometry(c)
        if _smooth_closed(curve):
            continue  # its seam is no corner (sketch edges run across it)
        for t in (0.0, 1.0):
            junction_points.append(np.asarray(curve.eval(t)))
    junction_points = np.asarray(junction_points).reshape(-1, 3)
    curves = {c: store.geometry(c) for c in curve_ids}
    scale = max(1.0, float(np.max(np.abs(points))))
    for loop in loops:
        m = len(loop)
        for k, node in enumerate(loop):
            prev = points[loop[k - 1]]
            here = points[node]
            nxt = points[loop[(k + 1) % m]]
            a_in, a_out = _mesh_directions(prev, here, nxt)
            junction = bool(np.min(np.linalg.norm(junction_points - here, axis=1)) < 1e-9 * scale)
            if junction:
                # Exact directions from the curves meeting there (the mesh polygon cuts
                # strongly curved ends short).
                a_in, a_out = _curve_directions(curves, prev, here, nxt, a_in, a_out)
            else:
                # The curve's own tangent at the node (boundary nodes the mesher inserted lie
                # on chords, where the polygon kinks).
                a_in = a_out = _curve_tangent(curves, here, nxt - prev)
            turn = math.degrees(_wrap(a_out - a_in))
            mean = np.exp(4j * a_in) + np.exp(4j * a_out)
            if abs(mean) < 1e-6:
                mean = np.exp(4j * a_out)
            values[node] = mean / abs(mean)
            fixed[node] = True
            if junction:
                interior = 180.0 - turn  # degrees, inside the region (on the left)
                quarters = max(1, round(interior / 90.0))
                corners.append(Corner(here.copy(), node, quarters, True))
    z = _laplace(n, triangles, fixed, values)
    z = _normalize_smooth(n, triangles, fixed, z, NORMALIZED_SWEEPS)
    singularities = []
    angle = np.angle(z)
    for index, t in enumerate(triangles):
        winding = 0.0
        for k in range(3):
            winding += _wrap(angle[t[(k + 1) % 3]] - angle[t[k]])
        k = round(winding / (2.0 * math.pi))
        if k != 0:
            singularities.append(Singularity(points[t].mean(axis=0), k / 4.0, index))
    return Field(points, triangles, z, loops, corners, singularities, len(loops) - 1, size)


# ---------------------------------------------------------------------------- tracing
class _Locator:
    """Triangle containing a point (uniform bins)."""

    def __init__(self, points, triangles, size):
        self.points = points
        self.triangles = triangles
        self.lo = points[:, :2].min(axis=0)
        self.cell = 2.0 * size
        corners = points[triangles][:, :, :2]
        tmin = np.floor((corners.min(axis=1) - self.lo) / self.cell).astype(int)
        tmax = np.floor((corners.max(axis=1) - self.lo) / self.cell).astype(int)
        edges = corners - np.roll(corners, 1, axis=1)
        self.sizes = np.linalg.norm(edges, axis=2).max(axis=1)  # longest edge per triangle
        self.bins: dict = {}
        for t in range(len(triangles)):
            for i in range(tmin[t, 0], tmax[t, 0] + 1):
                for j in range(tmin[t, 1], tmax[t, 1] + 1):
                    self.bins.setdefault((i, j), []).append(t)

    def locate(self, p):
        key = tuple(np.floor((np.asarray(p[:2]) - self.lo) / self.cell).astype(int))
        for t in self.bins.get(key, ()):
            a, b, c = self.points[self.triangles[t]][:, :2]
            m = np.array([[b[0] - a[0], c[0] - a[0]], [b[1] - a[1], c[1] - a[1]]])
            det = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]
            if det == 0.0:
                continue
            d = np.asarray(p[:2]) - a
            s = (m[1, 1] * d[0] - m[0, 1] * d[1]) / det
            r = (-m[1, 0] * d[0] + m[0, 0] * d[1]) / det
            if s >= -1e-12 and r >= -1e-12 and s + r <= 1.0 + 1e-12:
                return t, np.array([1.0 - s - r, s, r])
        return None, None


def _local_size(fld: Field, locator: _Locator, p) -> float:
    t, _w = locator.locate(p)
    if t is None:
        return fld.size
    return float(locator.sizes[t])


def _field_at(fld: Field, locator: _Locator, p):
    t, w = locator.locate(p)
    if t is None:
        return None
    return complex(np.dot(w, fld.z[fld.triangles[t]]))


def _closest_direction(z: complex, direction: np.ndarray) -> np.ndarray:
    base = np.angle(z) / 4.0
    best = None
    best_dot = -2.0
    for k in range(4):
        a = base + k * math.pi / 2.0
        d = np.array([math.cos(a), math.sin(a), 0.0])
        dot = float(d @ direction)
        if dot > best_dot:
            best = d
            best_dot = dot
    return best


def _directions_at(fld: Field, locator, centre, radius) -> list[np.ndarray]:
    """Directions from ``centre`` where a cross direction points radially."""
    samples = 720
    values = []
    for k in range(samples):
        a = 2.0 * math.pi * k / samples
        q = centre + radius * np.array([math.cos(a), math.sin(a), 0.0])
        z = _field_at(fld, locator, q)
        if z is None or z == 0:
            values.append(None)
            continue
        values.append(_wrap(4.0 * a - np.angle(z)) / 4.0)  # in (-pi/4, pi/4]
    out = []
    for k in range(samples):
        f0 = values[k]
        f1 = values[(k + 1) % samples]
        if f0 is None or f1 is None:
            continue
        if f0 <= 0.0 < f1 and f1 - f0 < 0.5:
            a = 2.0 * math.pi * (k + f0 / (f0 - f1)) / samples
            out.append(np.array([math.cos(a), math.sin(a), 0.0]))
    return out


def _segment_hit(p, q, a, b):
    """Parameter s in [0, 1] along p -> q where it crosses segment a -> b, or None."""
    r = q[:2] - p[:2]
    s = b[:2] - a[:2]
    den = r[0] * s[1] - r[1] * s[0]
    if den == 0.0:
        return None
    d = a[:2] - p[:2]
    t = (d[0] * s[1] - d[1] * s[0]) / den
    u = (d[0] * r[1] - d[1] * r[0]) / den
    if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
        return t
    return None


def _trace(fld: Field, locator, start, direction, singularities, curves, start_index, straight=0.0):
    """A separatrix from ``start`` along ``direction`` (straight for the first ``straight``
    length): (points, end kind, end index)."""
    h = 0.5 * _local_size(fld, locator, start)
    p = np.asarray(start, dtype=float) + h * direction
    d = direction
    path = [np.asarray(start, dtype=float), p.copy()]
    travelled = h
    while travelled < straight:
        q = p + h * d
        if locator.locate(q)[0] is None:
            break
        p = q
        path.append(p.copy())
        travelled += h
    diameter = float(np.max(np.ptp(fld.points[:, :2], axis=0)))
    boundary = [
        (fld.points[loop[k]], fld.points[loop[(k + 1) % len(loop)]])
        for loop in fld.boundary
        for k in range(len(loop))
    ]
    length = h
    while length < 6.0 * diameter:
        for s_index, s in enumerate(singularities):
            if s_index != start_index and np.linalg.norm(s.point - p) < s.capture:
                path.append(s.point.copy())
                return np.asarray(path), "source", s_index
        z = _field_at(fld, locator, p)
        if z is None:
            break
        # Steps of half the local mesh size (fine near strongly curved boundaries).
        h = 0.5 * _local_size(fld, locator, p)
        d1 = _closest_direction(z, d)
        mid = p + 0.5 * h * d1
        zm = _field_at(fld, locator, mid)
        if zm is None:
            q = p + h * d1
        else:
            d1 = _closest_direction(zm, d1)
            q = p + h * d1
        if locator.locate(q)[0] is None:
            # Leaving the region: end exactly on the boundary curve there.
            best = None
            for a, b in boundary:
                t = _segment_hit(p, q, a, b)
                if t is not None and (best is None or t < best):
                    best = t
            hit = q
            if best is not None:
                hit = p + best * (q - p)
            end, _curve = _on_curves(curves, hit)
            # Landing next to a boundary source (a junction or corner): end there; its own
            # separatrix in this direction is this one.
            for s_index, src in enumerate(singularities):
                if (
                    s_index != start_index
                    and src.junction
                    and _boundary_directions(src)
                    and np.linalg.norm(src.point - end) < 4.0 * src.capture
                ):
                    path.append(src.point.copy())
                    return np.asarray(path), "source", s_index
            path.append(end)
            return np.asarray(path), "boundary", None
        d = d1
        p = q
        path.append(p.copy())
        length += h
    end, _curve = _on_curves(curves, p)
    path.append(end)
    return np.asarray(path), "lost", None


def _on_curves(curves, point):
    best = None
    for c, curve in curves.items():
        _t, q, dist = curve.project(tuple(float(v) for v in point))
        if best is None or dist < best[2]:
            best = (np.asarray(q), c, dist)
    return best[0], best[1]


def _simplify(path: np.ndarray, spacing: float) -> np.ndarray:
    """Every ~spacing a point, the ends kept exactly."""
    keep = [path[0]]
    for q in path[1:-1]:
        if np.linalg.norm(q - keep[-1]) >= spacing:
            keep.append(q)
    if len(keep) > 1 and np.linalg.norm(path[-1] - keep[-1]) < 0.5 * spacing:
        keep.pop()
    keep.append(path[-1])
    return np.asarray(keep)


@dataclass
class Source:
    """Where separatrices start: an interior singularity (or several close ones), or a point
    on the boundary - a corner or curve junction, and singularities snapped onto it."""

    point: np.ndarray
    radius: float  # direction search circle and capture radius
    legs: list = field(default_factory=list)  # boundary directions at a boundary source
    quarters: int = 0  # boundary source: interior angle in quarter turns
    index: float = 0.0  # singularities merged into it
    capture: float = 0.0  # separatrices coming this close end here
    junction: bool = False  # a corner or curve junction (not a singularity snapped on)


SNAP = 4.0  # local mesh sizes within which singularities join the boundary or each other
STRAIGHT = 2.0  # boundary separatrices start straight for this many source radii
LEG_ANGLE = 15.0  # degrees: separatrix directions closer to a boundary leg are dropped


def _sources(fld: Field, locator, curves) -> list[Source]:
    sources = []
    for corner in fld.corners:
        loop = next(loop for loop in fld.boundary if corner.node in loop)
        k = loop.index(corner.node)
        prev = fld.points[loop[k - 1]] - corner.point
        nxt = fld.points[loop[(k + 1) % len(loop)]] - corner.point
        size = float(locator.sizes[_node_triangle(fld, corner.node)])
        sources.append(
            Source(
                corner.point.copy(),
                2.0 * size,
                [prev / np.linalg.norm(prev), nxt / np.linalg.norm(nxt)],
            )
        )
        sources[-1].quarters = corner.quarters
        sources[-1].junction = True
    boundary_points = np.concatenate([fld.points[loop] for loop in fld.boundary])
    for s in fld.singularities:
        size = float(locator.sizes[s.triangle])
        gaps = np.linalg.norm(boundary_points - s.point, axis=1)
        if float(gaps.min()) < SNAP * size:
            # Next to the boundary: onto it, joining a corner there if there is one.
            near = [
                b
                for b in sources
                if b.legs
                and np.linalg.norm(b.point - s.point) < 2.0 * SNAP * size
                # the corner is (about) the nearest boundary point, not one further along
                and np.linalg.norm(b.point - s.point) < 2.5 * float(gaps.min())
                # a +1/4 is a convex corner of its own (a mark), never joined to a junction
                and s.index < 0
            ]
            if near:
                target = min(near, key=lambda b: float(np.linalg.norm(b.point - s.point)))
                target.index += s.index
                target.radius = max(
                    target.radius, 2.0 * float(np.linalg.norm(target.point - s.point)) + 2.0 * size
                )
                continue
            point, _curve = _on_curves(curves, s.point)
            k = int(np.argmin(np.linalg.norm(boundary_points - point, axis=1)))
            legs = _legs_at(fld, point)
            sources.append(Source(point, 2.0 * size + float(np.linalg.norm(point - s.point)), legs))
            sources[-1].quarters = 2
            sources[-1].index = s.index
            del k
            continue
        near = [
            b for b in sources if not b.legs and np.linalg.norm(b.point - s.point) < SNAP * size
        ]
        if near:
            target = near[0]
            target.radius = max(
                target.radius, 2.0 * float(np.linalg.norm(target.point - s.point)) + 2.0 * size
            )
            continue
        sources.append(Source(s.point.copy(), 2.0 * size))
    for source in sources:
        source.capture = 1.5 * _local_size(fld, locator, source.point)
    return sources


def _boundary_directions(source: Source) -> list:
    """A boundary source's separatrices: quarters - 1 - 4 index of them, evenly spaced over
    the interior angle between its legs (counter-clockwise from the outgoing leg)."""
    count = round(source.quarters - 1 - 4.0 * source.index)
    if count <= 0 or len(source.legs) < 2:
        return []
    incoming, outgoing = source.legs  # towards the previous and the next boundary node
    a0 = math.atan2(outgoing[1], outgoing[0])
    a1 = math.atan2(incoming[1], incoming[0])
    sweep = (a1 - a0) % (2.0 * math.pi)  # the region lies on the left of the boundary
    out = []
    for j in range(1, count + 1):
        a = a0 + sweep * j / (count + 1)
        out.append(np.array([math.cos(a), math.sin(a), 0.0]))
    return out


def _node_triangle(fld: Field, node: int) -> int:
    return int(np.nonzero(np.any(fld.triangles == node, axis=1))[0][0])


def _legs_at(fld: Field, point) -> list:
    """Boundary directions at a boundary point: towards the boundary nodes on either side."""
    best = None
    for loop in fld.boundary:
        pts = fld.points[loop]
        for k in range(len(loop)):
            a = pts[k]
            b = pts[(k + 1) % len(loop)]
            ab = b - a
            t = float(np.clip((point - a)[:2] @ ab[:2] / max(float(ab[:2] @ ab[:2]), 1e-300), 0, 1))
            d = float(np.linalg.norm(a + t * ab - point))
            if best is None or d < best[0]:
                best = (d, a, b)
    _d, a, b = best
    out = []
    for q in (a, b):
        v = q - point
        if np.linalg.norm(v) > 0:
            out.append(v / np.linalg.norm(v))
    return out


def trace_layout(store, curve_ids, size=None) -> Layout:
    """The cross field and its separatrices (polylines)."""
    fld = compute_field(store, curve_ids, size)
    locator = _Locator(fld.points, fld.triangles, fld.size)
    curves = {c: store.geometry(c) for c in curve_ids}
    warnings = []
    expected = fld.euler()
    if abs(fld.index_sum() - expected) > 1e-9:
        warnings.append(
            f"cross field: singularity indices add up to {fld.index_sum():.2f}, the Euler "
            f"characteristic is {expected} (the field is not resolved: refine)"
        )
    sources = _sources(fld, locator, curves)
    separatrices = []
    leaving: dict[int, list[np.ndarray]] = {}
    cosine = math.cos(math.radians(20.0))
    marks = []
    for index, source in enumerate(sources):
        if source.legs and abs(source.index) > 0 and not _boundary_directions(source):
            # A singularity on a smooth boundary that emits nothing (+1/4 there: a convex
            # corner): still a corner of the face around it.
            marks.append(source.point.copy())
        if source.legs:
            directions = _boundary_directions(source)
        else:
            directions = _directions_at(fld, locator, source.point, source.radius)
        for d in directions:
            if any(float(d @ leg) > math.cos(math.radians(LEG_ANGLE)) for leg in source.legs):
                continue  # along the boundary
            if any(float(d @ a) > cosine for a in leaving.get(index, ())):
                continue  # already traced (from the other end, or a duplicate root)
            leaving.setdefault(index, []).append(d)
            straight = 0.0
            if source.legs:
                straight = STRAIGHT * source.radius
            path, kind, other = _trace(
                fld, locator, source.point, d, sources, curves, index, straight
            )
            if kind == "lost":
                warnings.append("cross field: a separatrix did not reach the boundary (stopped)")
            if kind == "source":
                # The arrival direction from a few steps back (the last step can be the
                # jump onto the source).
                back = path[max(0, len(path) - 8)] - path[-1]
                back = back / np.linalg.norm(back)
                if any(
                    float(back @ a) > math.cos(math.radians(30.0)) for a in leaving.get(other, ())
                ):
                    continue  # the other source traced this one already (from its end)
                leaving.setdefault(other, []).append(back)
            separatrices.append(path)
    return Layout(fld, [_simplify(p, fld.size) for p in separatrices], warnings, marks)


# ------------------------------------------------------------------------- the sketch
def auto_topology(store, curves, size=None, spacing=None, name=None, relax=False):
    """Separatrix polylines and a sketch of the four-sided faces (inside the open
    transaction). Returns ([sketch name], warnings)."""
    from opengrid.model.entities import normalize_params
    from opengrid.model.network import build_network
    from opengrid.model.sketch import create_sketch
    from opengrid.model.topology_edit import fresh_name

    if isinstance(curves, (str, int)):
        curves = [curves]
    ids = []
    for c in curves:
        if c not in store or store.get(c).kind not in CURVE_KINDS:
            _fail(f"{c!r} is not a curve")
        ids.append(store.get(c).id)
    layout = trace_layout(store, ids, size)
    polylines = []
    for k, path in enumerate(layout.separatrices):
        base = fresh_name(store, f"separatrix-{k + 1}")
        params = normalize_params("polyline", {"points": tuple(map(tuple, path))}, base, store)
        polylines.append(store.add("polyline", params, base).id)
    net = build_network(store, ids + polylines, corner_angle=60.0)
    mark_vertices = set()
    for mark in layout.marks:
        split = _split_at(store, net, mark, set(ids))
        if split is not None:
            mark_vertices.add(split)
    locator = _Locator(layout.field.points, layout.field.triangles, layout.field.size)
    vertices = [v.copy() for v in net.vertices]
    blocks = []
    assigned = {}
    issues = []
    # A vertex where only the two pieces of one closed curve meet (its seam) is no corner.
    degree: dict = {}
    for edge in net.edges:
        for v in (edge.v0, edge.v1):
            degree.setdefault(v, set()).add(edge.curve)
    count_at: dict = {}
    for edge in net.edges:
        for v in (edge.v0, edge.v1):
            count_at[v] = count_at.get(v, 0) + 1
    seams = {
        v for v, n in count_at.items() if n == 2 and len(degree[v]) == 1 and v not in mark_vertices
    }
    used: dict = {}
    for face in net.faces:
        if locator.locate(_interior_point(net, face))[0] is None:
            continue  # a hole
        corners = []
        curves_of = []
        for e, forward in face:
            edge = net.edges[e]
            start_vertex = edge.v0
            if not forward:
                start_vertex = edge.v1
            if start_vertex not in seams:
                corners.append(start_vertex)
                curves_of.append(edge.curve)
        if len(corners) != 4:
            where = np.round(np.mean([vertices[v] for v in corners], axis=0)[:2], 3)
            issues.append(
                Issue(
                    "automatic topology",
                    f"a face near ({where[0]}, {where[1]}) has {len(corners)} corners, not 4",
                )
            )
            continue
        block = []
        for k, v in enumerate(corners):
            if v not in used:
                used[v] = len(used)
            block.append(used[v])
            nxt = corners[(k + 1) % 4]
            assigned[(v, nxt)] = curves_of[k]
        blocks.append(tuple(block))
    if issues:
        raise ValidationError(issues)
    if spacing is None:
        spacing = 2.0 * layout.field.size
    if name is None:
        name = fresh_name(store, "auto")
    order = sorted(used, key=lambda v: used[v])
    curves_by_edge = {(used[a], used[b]): c for (a, b), c in assigned.items()}
    names, warnings = create_sketch(
        store,
        name=name,
        vertices=tuple(tuple(float(c) for c in vertices[v]) for v in order),
        blocks=tuple(blocks),
        curves=curves_by_edge,
        spacing=float(spacing),
        relax=relax,
    )
    return names, [*layout.warnings, *warnings]


def _split_at(store, net, point, boundary) -> int | None:
    """Make a boundary point a vertex of the network: split the boundary edge through it.
    Returns the new vertex."""
    from opengrid.model.network import Edge

    best = None
    for e, edge in enumerate(net.edges):
        if edge.curve not in boundary:
            continue
        curve = store.geometry(edge.curve)
        t, q, d = curve.project(tuple(float(c) for c in point))
        lo, hi = sorted((edge.t0, edge.t1))
        if lo < t < hi and (best is None or d < best[0]):
            best = (d, e, t, np.asarray(q))
    if best is None:
        return None
    _d, e, t, q = best
    edge = net.edges[e]
    net.vertices.append(q)
    v = len(net.vertices) - 1
    curve = store.geometry(edge.curve)

    def piece(t0, t1, v0, v1):
        sample = np.asarray(curve.eval_many(np.linspace(t0, t1, 33)))
        length = float(np.sum(np.linalg.norm(np.diff(sample, axis=0), axis=1)))
        return Edge(edge.curve, t0, t1, v0, v1, length, sample)

    net.edges[e] = piece(edge.t0, t, edge.v0, v)
    net.edges.append(piece(t, edge.t1, v, edge.v1))
    new = len(net.edges) - 1
    for f, face in enumerate(net.faces):
        out = []
        for half in face:
            if half[0] != e:
                out.append(half)
            elif half[1]:
                out.extend([(e, True), (new, True)])
            else:
                out.extend([(new, False), (e, False)])
        net.faces[f] = out
    return v


def _interior_point(net, face) -> np.ndarray:
    """A point just inside a face: off the middle of its longest edge, to the left."""
    best = None
    for e, forward in face:
        pts = net.edges[e].points
        if not forward:
            pts = pts[::-1]
        length = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        if best is None or length > best[0]:
            best = (length, pts)
    length, pts = best
    k = len(pts) // 2
    t = pts[min(k + 1, len(pts) - 1)] - pts[max(k - 1, 0)]
    left = np.array([-t[1], t[0], 0.0])
    left /= np.linalg.norm(left)
    return pts[k] + 1e-3 * length * left


__all__ = ["auto_topology", "compute_field", "trace_layout"]
