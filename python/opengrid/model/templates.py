"""Topology templates: parametrised generators that build a 10d sketch (and the far-field
geometry they need) from the selected geometry, with no user-placed vertices.

2D (curves in a plane z = const):
  * ``ogrid``: an O-grid around a closed body: ``blocks`` (>= 4) blocks between the body and
    a circle around it; the body is cut at its curve junctions and, to reach the block count,
    where rays from its centroid meet it.
  * ``butterfly``: the inside of a closed curve, a core block and four around it (O-H).
  * ``cgrid``: a C-grid around a body with a sharp trailing edge (the sharpest corner): two
    layers of blocks around the body and its wake to a C-shaped far field (an arc of
    ``radius`` around the trailing edge and an outflow ``downstream`` of it, along the chord);
    the wake cut is an internal wall.
  * ``hgrid``: a channel from its four sides (``bottom``, ``right``, ``top``, ``left``), an
    ``nx`` x ``ny`` lattice of blocks (TFI of the sides).
  * ``body_in_channel``: a body in a channel: a 4 x 4 lattice of the channel without its
    centre cell, which an O-grid of four blocks around the body replaces (each channel corner
    in one block).
3D: ``pipe`` (butterfly), ``cylinder_in_channel`` (body in channel), ``bump_channel``
(H-grid) and ``wing_in_box`` (C-grid) are the 2D template with ``extrude``: its smoothed grid
extruded along z by ``length`` (``span_points``, ``span_distribution`` with spacings as
lengths, e.g. clustered towards the side walls), the same in every section; with
``three_d_blocks=True`` a 3D sketch of six-face blocks instead (smoothed in 3D with the end
planes as symmetry planes, for editing the layout in 3D); ``wing_in_box`` takes ``taper``
(tip / root chord, about the leading edge) and ``sweep`` (degrees of the leading edge),
which make it such a 3D sketch with the wing on ruled surfaces between root and tip;
``sphere_in_box`` builds a cubed-sphere O-grid inside a 26-block box layer around a sphere it
creates.
Boundary edges are assigned to the curves (surfaces) that contain them; the template places
the boundary vertices exactly on the geometry (fixed) and the interior ones float in the
smoothing group (10d), so nothing is placed by hand. ``layers`` / ``growth`` (every template):
wall layers on all its walls (11a: the first layers marched normal to the wall and kept).

Method: the standard multi-block topologies of Thompson, Warsi and Mastin;
references: J. F. Thompson, Z. U. A. Warsi, C. W. Mastin, Numerical Grid Generation (1985),
ch. VI-VII; P. R. Eiseman, NASA CP-3291 (1995).
"""

import itertools
import math

import numpy as np

from opengrid.model.connectors import CURVE_KINDS
from opengrid.model.entities import normalize_params
from opengrid.validation import Issue, ValidationError

KINDS_2D = ("ogrid", "butterfly", "cgrid", "hgrid", "body_in_channel")
KINDS_3D = ("pipe", "cylinder_in_channel", "bump_channel", "wing_in_box", "sphere_in_box")
EXTRUDED = {
    "pipe": "butterfly",
    "cylinder_in_channel": "body_in_channel",
    "bump_channel": "hgrid",
    "wing_in_box": "cgrid",
}
SAMPLES = 400  # loop samples per curve


def _fail(where: str, message: str):
    raise ValidationError([Issue(where, message)])


# ------------------------------------------------------------------------------- geometry
def _curve_ids(store, names, where) -> list[int]:
    if isinstance(names, (str, int)):
        names = [names]
    out = []
    for name in names:
        if name not in store or store.get(name).kind not in CURVE_KINDS:
            _fail(where, f"{name!r} is not a curve")
        out.append(store.get(name).id)
    if not out:
        _fail(where, "no curves given")
    return out


class Loop:
    """A closed loop of curves, counter-clockwise: pieces (curve id, reversed), samples."""

    def __init__(self, store, curve_ids, where: str):
        curves = {c: store.geometry(c) for c in curve_ids}
        ends = {
            c: (np.asarray(curves[c].eval(0.0)), np.asarray(curves[c].eval(1.0))) for c in curve_ids
        }
        scale = max(1.0, max(float(np.max(np.abs(e))) for pair in ends.values() for e in pair))
        tol = 1e-6 * scale
        order = [(curve_ids[0], False)]
        left = list(curve_ids[1:])
        while left:
            c, rev = order[-1]
            tail = ends[c][1]
            if rev:
                tail = ends[c][0]
            for other in left:
                if np.linalg.norm(ends[other][0] - tail) <= tol:
                    order.append((other, False))
                    break
                if np.linalg.norm(ends[other][1] - tail) <= tol:
                    order.append((other, True))
                    break
            else:
                _fail(where, "the curves do not form one closed loop")
            left.remove(order[-1][0])
        first, rev = order[0]
        head = ends[first][0]
        if rev:
            head = ends[first][1]
        last, rev = order[-1]
        tail = ends[last][1]
        if rev:
            tail = ends[last][0]
        if np.linalg.norm(head - tail) > tol:
            _fail(where, "the curves do not form a closed loop")
        self.curves = curves
        self.pieces = order
        self.tol = tol
        self._sample()
        if self.area < 0.0:
            self.pieces = [(c, not rev) for c, rev in reversed(self.pieces)]
            self._sample()

    def _sample(self):
        points = []
        owner = []
        starts = []
        for c, rev in self.pieces:
            t = np.linspace(0.0, 1.0, SAMPLES + 1)
            if rev:
                t = t[::-1]
            p = np.asarray(self.curves[c].eval_many(t))
            starts.append(len(points))
            points.extend(p[:-1])
            owner.extend([c] * SAMPLES)
        self.points = np.asarray(points)
        self.owner = owner
        self.starts = starts
        closed = np.concatenate([self.points, self.points[:1]])
        self.s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))])
        self.length = float(self.s[-1])
        x = self.points[:, 0]
        y = self.points[:, 1]
        self.area = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
        cx = float(np.sum((x + np.roll(x, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)))
        cy = float(np.sum((y + np.roll(y, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)))
        self.centroid = np.array(
            [cx / (6 * self.area), cy / (6 * self.area), float(self.points[0, 2])]
        )

    def junctions(self) -> list[int]:
        """Sample indices where the curves meet (none for a single closed curve)."""
        if len(self.pieces) == 1:
            return []
        return list(self.starts)

    def turn(self, k: int) -> float:
        """Turning angle (degrees, signed) of the loop at sample k."""
        n = len(self.points)
        a = self.points[(k - 1) % n]
        b = self.points[k]
        c = self.points[(k + 1) % n]
        u = b - a
        v = c - b
        return math.degrees(math.atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1]))

    def ray(self, angle: float) -> int:
        """Sample whose direction from the centroid is closest to ``angle``."""
        d = self.points - self.centroid
        angles = np.arctan2(d[:, 1], d[:, 0])
        diff = np.abs((angles - angle + math.pi) % (2 * math.pi) - math.pi)
        return int(np.argmin(diff))

    def normal(self, k: int) -> np.ndarray:
        """Outward unit normal (the loop is counter-clockwise)."""
        n = len(self.points)
        t = self.points[(k + 1) % n] - self.points[(k - 1) % n]
        out = np.array([t[1], -t[0], 0.0])
        return out / np.linalg.norm(out)


def _choose(loop: Loop, required: list[int], candidates: list[int], count: int) -> list[int]:
    """Loop sample indices: all ``required`` ones, then candidates far (in arc length) from
    those already taken, until ``count``; in loop order."""
    taken = sorted(set(required))
    for c in candidates:
        if len(taken) >= count:
            break
        gaps = [abs(loop.s[c] - loop.s[t]) for t in taken]
        gaps = [min(g, loop.length - g) for g in gaps]
        if not taken or min(gaps) > 0.08 * loop.length:
            taken.append(c)
    while len(taken) < count:
        # Split the longest arc between taken samples.
        taken.sort()
        arcs = []
        for a, b in zip(taken, [*taken[1:], taken[0]], strict=True):
            length = (loop.s[b] - loop.s[a]) % loop.length
            arcs.append((length, a, b))
        _length, a, b = max(arcs)
        target = (loop.s[a] + 0.5 * ((loop.s[b] - loop.s[a]) % loop.length)) % loop.length
        taken.append(int(np.argmin(np.abs(loop.s[:-1] - target))))
    return sorted(set(taken))


def _arc_between(loop: Loop, pa, pb) -> float:
    """Length along the loop from pa to pb (counter-clockwise)."""
    ka = int(np.argmin(np.linalg.norm(loop.points - pa, axis=1)))
    kb = int(np.argmin(np.linalg.norm(loop.points - pb, axis=1)))
    return float((loop.s[kb] - loop.s[ka]) % loop.length)


def _add(store, kind: str, raw: dict, base: str) -> int:
    from opengrid.model.topology_edit import fresh_name

    name = fresh_name(store, base)
    params = normalize_params(kind, raw, name, store)
    return store.add(kind, params, name).id


def _assign(store, vertices: np.ndarray, edges, curve_ids) -> dict:
    """Boundary edge -> the curve holding both its ends (and nearest at its middle)."""
    out = {}
    scale = max(1.0, float(np.max(np.abs(vertices))))
    tol = 1e-6 * scale
    for a, b in edges:
        best = None
        for c in curve_ids:
            curve = store.geometry(c)
            ends = max(curve.project(tuple(vertices[v]))[2] for v in (a, b))
            if ends > tol:
                continue
            middle = curve.project(tuple(0.5 * (vertices[a] + vertices[b])))[2]
            if middle > 0.5 * float(np.linalg.norm(vertices[b] - vertices[a])):
                continue
            if best is None or middle < best[1]:
                best = (c, middle)
        if best is not None:
            out[(a, b)] = best[0]
    return out


def _boundary_edges(blocks) -> list[tuple[int, int]]:
    uses: dict = {}
    for block in blocks:
        for k in range(4):
            key = tuple(sorted((block[k], block[(k + 1) % 4])))
            uses[key] = uses.get(key, 0) + 1
    return sorted(k for k, n in uses.items() if n == 1)


def _option(options: dict, key: str, default):
    value = options.pop(key, None)
    if value is None:
        return default
    return value


# ------------------------------------------------------------------------- 2D templates
def _ogrid(store, options) -> dict:
    loop = Loop(store, _curve_ids(store, options.pop("body"), "ogrid body"), "ogrid body")
    size = math.sqrt(abs(loop.area))
    radius = float(_option(options, "radius", 10.0 * size))
    count = int(_option(options, "blocks", 4))
    if count < 4:
        _fail("ogrid", "at least 4 blocks")
    sharpest = max(abs(loop.turn(k)) for k in range(len(loop.points)))
    if sharpest > 120.0:
        _fail(
            "ogrid",
            f"the body has a cusp (it turns by {sharpest:.0f} degrees): an O-grid folds there; "
            "use the cgrid template",
        )
    first = 0.25 * math.pi
    if loop.junctions():
        d = loop.points[loop.junctions()[0]] - loop.centroid
        first = math.atan2(d[1], d[0])
    rays = [loop.ray(first + 2 * math.pi * k / count) for k in range(count)]
    body = _choose(loop, loop.junctions(), rays, count)
    c = loop.centroid
    circle = _add(
        store,
        "arc",
        {"center": tuple(c), "radius": radius, "start_angle": 0.0, "end_angle": 2 * math.pi},
        "farfield",
    )
    vertices = [loop.points[k] for k in body]
    for k in body:
        d = loop.points[k] - c
        angle = math.atan2(d[1], d[0])
        vertices.append(c + radius * np.array([math.cos(angle), math.sin(angle), 0.0]))
    n = len(body)
    blocks = [(k, (k + 1) % n, n + (k + 1) % n, n + k) for k in range(n)]
    return _finish_2d(
        store,
        options,
        vertices,
        blocks,
        [p[0] for p in loop.pieces] + [circle],
        {
            "radial": ((0, n), _option(options, "radial_points", 31)),
            "walls": [(k, (k + 1) % n) for k in range(n)],
            "fixed": range(2 * n),
        },
    )


def _butterfly(store, options) -> dict:
    loop = Loop(store, _curve_ids(store, options.pop("boundary"), "butterfly"), "butterfly")
    core = float(_option(options, "core", 0.5))
    first = 0.25 * math.pi
    if loop.junctions():
        d = loop.points[loop.junctions()[0]] - loop.centroid
        first = math.atan2(d[1], d[0])
    rays = [loop.ray(first + 0.5 * math.pi * k) for k in range(4)]
    body = _choose(loop, loop.junctions(), rays, 4)
    if len(body) != 4:
        _fail("butterfly", f"the boundary has {len(body)} curve junctions (at most 4)")
    c = loop.centroid
    vertices = [loop.points[k] for k in body]
    vertices += [c + core * (loop.points[k] - c) for k in body]
    blocks = [(4, 5, 6, 7)] + [(k, (k + 1) % 4, 4 + (k + 1) % 4, 4 + k) for k in range(4)]
    return _finish_2d(
        store,
        options,
        vertices,
        blocks,
        [p[0] for p in loop.pieces],
        {
            "radial": ((0, 4), _option(options, "radial_points", 11)),
            "walls": [(k, (k + 1) % 4) for k in range(4)],
            "fixed": range(4),
        },
    )


def _cgrid(store, options) -> dict:
    loop = Loop(store, _curve_ids(store, options.pop("body"), "cgrid body"), "cgrid body")
    turns = [abs(loop.turn(k)) for k in range(len(loop.points))]
    te = int(np.argmax(turns))
    if turns[te] < 30.0:
        _fail("cgrid", "the body has no sharp trailing edge (a corner of 30 degrees or more)")
    le = int(np.argmax(np.linalg.norm(loop.points - loop.points[te], axis=1)))
    p_te = loop.points[te]
    chord = float(np.linalg.norm(loop.points[le] - p_te))
    e = (p_te - loop.points[le]) / chord
    nrm = np.array([-e[1], e[0], 0.0])
    radius = float(_option(options, "radius", 5.0 * chord))
    downstream = float(_option(options, "downstream", 5.0 * chord))
    ring = float(_option(options, "ring", 0.6 * chord))
    n_points = len(loop.points)
    # Body vertices from the trailing edge around (counter-clockwise): junctions, the leading
    # edge, and the middle of each side when it has none.
    required = sorted(set(loop.junctions()) | {te, le})
    order = sorted(required, key=lambda k: (k - te) % n_points)
    rel = [(k - te) % n_points for k in order]
    le_pos = rel.index((le - te) % n_points)
    sides = [order[1:le_pos], order[le_pos + 1 :]]
    middles = []
    for side, (lo, hi) in zip(sides, ((te, le), (le, te)), strict=True):
        if not side:
            target = (loop.s[lo] + 0.5 * ((loop.s[hi] - loop.s[lo]) % loop.length)) % loop.length
            middles.append(int(np.argmin(np.abs(loop.s[:-1] - target))))
    body = sorted(set(required) | set(middles), key=lambda k: (k - te) % n_points)
    outflow = p_te + downstream * e
    # C along the body: outflow, trailing edge, the body around, trailing edge, outflow.
    base_points = [outflow, p_te] + [loop.points[k] for k in body[1:]] + [p_te, outflow]
    m = len(base_points)
    ring_points = [outflow + ring * nrm, p_te + ring * nrm]
    far_points = [outflow + radius * nrm, p_te + radius * nrm]
    for k in body[1:]:
        r = loop.points[k] + ring * loop.normal(k)
        ring_points.append(r)
        d = r - p_te
        angle = math.atan2(float(d @ nrm), float(d @ e))
        if angle >= 0:
            angle = min(max(angle, 0.5 * math.pi), math.pi)
        else:
            angle = max(min(angle, -0.5 * math.pi), -math.pi)
        far_points.append(p_te + radius * (math.cos(angle) * e + math.sin(angle) * nrm))
    ring_points += [p_te - ring * nrm, outflow - ring * nrm]
    far_points += [p_te - radius * nrm, outflow - radius * nrm]
    # Vertex numbering: outflow end of the cut (one vertex), trailing edge (one), the body,
    # then the ring and the far field.
    base_ids = [0, 1, *range(2, m - 2), 1, 0]
    vertices = [outflow, p_te, *base_points[2 : m - 2]]
    ring_ids = list(range(len(vertices), len(vertices) + m))
    vertices += ring_points
    far_ids = list(range(len(vertices), len(vertices) + m))
    vertices += far_points
    blocks = [(base_ids[k], base_ids[k + 1], ring_ids[k + 1], ring_ids[k]) for k in range(m - 1)]
    blocks += [(ring_ids[k], ring_ids[k + 1], far_ids[k + 1], far_ids[k]) for k in range(m - 1)]
    arc = _add(
        store,
        "arc",
        {
            "center": tuple(p_te),
            "radius": radius,
            "normal": (0.0, 0.0, 1.0),
            "x_axis": tuple(nrm),
            "start_angle": 0.0,
            "end_angle": math.pi,
        },
        "farfield-arc",
    )
    top = _add(
        store,
        "line",
        {"p0": tuple(p_te + radius * nrm), "p1": tuple(far_points[0])},
        "farfield-top",
    )
    bottom = _add(
        store,
        "line",
        {"p0": tuple(p_te - radius * nrm), "p1": tuple(far_points[-1])},
        "farfield-bottom",
    )
    out = _add(store, "line", {"p0": tuple(far_points[-1]), "p1": tuple(far_points[0])}, "outflow")
    body_edges = [(base_ids[k], base_ids[k + 1]) for k in range(1, m - 2)]
    inner, outer = _option(options, "radial_points", (31, 31))
    wake_points = int(_option(options, "wake_points", 41))
    body_spacing = float(_option(options, "spacing", chord / 30.0))
    options["spacing"] = body_spacing
    # Clustering along the body at the trailing and leading edges and along the wake cut at
    # the trailing edge (by default 1/200 and 1/250 of the chord).
    te_spacing = float(_option(options, "te_spacing", chord / 200.0))
    le_spacing = float(_option(options, "le_spacing", chord / 250.0))
    le_id = base_ids[1 + body.index(le)]
    distributions = {(1, 0): {"type": "tanh", "begin_spacing": te_spacing, "end_spacing": None}}
    for k in range(1, m - 2):
        a, b = base_ids[k], base_ids[k + 1]
        ends = {1: te_spacing, le_id: le_spacing}
        begin = ends.get(a)
        end = ends.get(b)
        if begin is None and end is None:
            continue
        body_length = _arc_between(loop, vertices[a], vertices[b])
        # The ring and far-field edges across the blocks get the same relative clustering,
        # so the blocks' TFI starts nearly orthogonal (the smoothing then moves little).
        for p0, p1 in ((a, b), (ring_ids[k], ring_ids[k + 1]), (far_ids[k], far_ids[k + 1])):
            ratio = 1.0
            if (p0, p1) != (a, b):
                ratio = float(np.linalg.norm(vertices[p1] - vertices[p0])) / body_length
            scaled = {"type": "tanh", "begin_spacing": None, "end_spacing": None}
            if begin is not None:
                scaled["begin_spacing"] = begin * ratio
            if end is not None:
                scaled["end_spacing"] = end * ratio
            distributions[(p0, p1)] = scaled
    plan = {
        "radial": ((1, ring_ids[1]), inner),
        "radial2": ((ring_ids[1], far_ids[1]), outer),
        "wake": ((0, 1), wake_points),
        # Points along the body from the body's own lengths (each family also holds the much
        # longer ring and far-field edges).
        "extra": {
            (a, b): max(
                5, math.ceil(_arc_between(loop, vertices[a], vertices[b]) / body_spacing) + 1
            )
            for a, b in body_edges
        },
        "walls": body_edges,
        "internal": [(0, 1)],
        # Body and far-field vertices stay at their places (sliding, they drift towards the
        # trailing edge); the ring and the outflow vertices relax.
        "fixed": list(range(2, m - 2)) + far_ids[2 : m - 2],
        "distributions": distributions,
    }
    curves = [p[0] for p in loop.pieces] + [arc, top, bottom, out]
    params = _finish_2d(store, options, vertices, blocks, curves, plan)
    params["_body"] = tuple(p[0] for p in loop.pieces)
    params["_leading_edge"] = tuple(loop.points[le].tolist())
    params["_chord"] = tuple(e.tolist())
    return params


def _corners(store, sides: dict, where: str):
    """Corner points and the four sides oriented bottom: c00 -> c10, right: c10 -> c11, top:
    c01 -> c11, left: c00 -> c01."""
    curves = {k: store.geometry(_curve_ids(store, sides[k], where)[0]) for k in sides}
    ends = {k: (np.asarray(c.eval(0.0)), np.asarray(c.eval(1.0))) for k, c in curves.items()}

    def meet(a, b):
        for p in ends[a]:
            for q in ends[b]:
                if np.linalg.norm(p - q) <= 1e-6 * max(1.0, float(np.max(np.abs(p)))):
                    return p
        _fail(where, f"{a} and {b} do not meet")

    c00 = meet("bottom", "left")
    c10 = meet("bottom", "right")
    c11 = meet("top", "right")
    c01 = meet("top", "left")
    oriented = {}
    for key, start in (("bottom", c00), ("right", c10), ("top", c01), ("left", c00)):
        curve = curves[key]
        flip = np.linalg.norm(ends[key][0] - start) > np.linalg.norm(ends[key][1] - start)
        oriented[key] = (curve, flip)
    return (c00, c10, c11, c01), oriented


def _side_point(side, fraction: float) -> np.ndarray:
    curve, flip = side
    if flip:
        fraction = 1.0 - fraction
    return np.asarray(curve.eval(curve.t_at_length(fraction * curve.length())))


def _tfi(oriented, corners, u: float, v: float) -> np.ndarray:
    c00, c10, c11, c01 = corners
    bottom = _side_point(oriented["bottom"], u)
    top = _side_point(oriented["top"], u)
    left = _side_point(oriented["left"], v)
    right = _side_point(oriented["right"], v)
    return (
        (1 - v) * bottom
        + v * top
        + (1 - u) * left
        + u * right
        - ((1 - u) * (1 - v) * c00 + u * (1 - v) * c10 + u * v * c11 + (1 - u) * v * c01)
    )


def _channel(store, options, where):
    sides = {k: options.pop(k) for k in ("bottom", "right", "top", "left")}
    corners, oriented = _corners(store, sides, where)
    ids = [_curve_ids(store, sides[k], where)[0] for k in ("bottom", "right", "top", "left")]
    return corners, oriented, ids


def _lattice(oriented, corners, us, vs):
    index = {}
    vertices = []
    for j, v in enumerate(vs):
        for i, u in enumerate(us):
            index[(i, j)] = len(vertices)
            vertices.append(_tfi(oriented, corners, u, v))
    return index, vertices


def _hgrid(store, options) -> dict:
    corners, oriented, ids = _channel(store, options, "hgrid")
    nx, ny = _option(options, "blocks", (3, 1))
    us = [i / nx for i in range(nx + 1)]
    vs = [j / ny for j in range(ny + 1)]
    index, vertices = _lattice(oriented, corners, us, vs)
    blocks = [
        (index[(i, j)], index[(i + 1, j)], index[(i + 1, j + 1)], index[(i, j + 1)])
        for j in range(ny)
        for i in range(nx)
    ]
    walls = []
    for side in _option(options, "wall_sides", ("bottom",)):
        if side == "bottom":
            walls += [(index[(i, 0)], index[(i + 1, 0)]) for i in range(nx)]
        elif side == "top":
            walls += [(index[(i, ny)], index[(i + 1, ny)]) for i in range(nx)]
        else:
            _fail("hgrid", f"wall_sides: bottom or top, got {side!r}")
    return _finish_2d(
        store,
        options,
        vertices,
        blocks,
        ids,
        {
            "radial": ((index[(0, 0)], index[(0, 1)]), _option(options, "normal_points", None)),
            "walls": walls,
        },
    )


def _body_in_channel(store, options) -> dict:
    loop = Loop(store, _curve_ids(store, options.pop("body"), "body_in_channel"), "body")
    corners, oriented, ids = _channel(store, options, "body_in_channel")
    lo = loop.points.min(axis=0)
    hi = loop.points.max(axis=0)
    size = float(np.linalg.norm(hi - lo))
    margin = float(_option(options, "margin", 0.5 * size))
    c00, c10, _c11, c01 = corners
    width = float(np.linalg.norm(c10 - c00))
    height = float(np.linalg.norm(c01 - c00))
    us = [0.0, (lo[0] - margin - c00[0]) / width, (hi[0] + margin - c00[0]) / width, 1.0]
    vs = [0.0, (lo[1] - margin - c00[1]) / height, (hi[1] + margin - c00[1]) / height, 1.0]
    if not (0.0 < us[1] < us[2] < 1.0 and 0.0 < vs[1] < vs[2] < 1.0):
        _fail("body_in_channel", "the body (with its margin) does not fit in the channel")
    index, vertices = _lattice(oriented, corners, us, vs)
    blocks = []
    for j, i in itertools.product(range(3), range(3)):
        if (i, j) != (1, 1):
            blocks.append(
                (index[(i, j)], index[(i + 1, j)], index[(i + 1, j + 1)], index[(i, j + 1)])
            )
    inner = [index[(1, 1)], index[(2, 1)], index[(2, 2)], index[(1, 2)]]
    angles = []
    for q in inner:
        d = vertices[q] - loop.centroid
        angles.append(math.atan2(d[1], d[0]))
    rays = [loop.ray(a) for a in angles]
    junctions = loop.junctions()
    if len(junctions) > 4:
        _fail("body_in_channel", f"the body has {len(junctions)} curve junctions (at most 4)")
    chosen = list(rays)
    for j in junctions:
        # Each junction replaces the ray point nearest to it.
        k = min(range(4), key=lambda r: abs(loop.s[chosen[r]] - loop.s[j]))
        chosen[k] = j
    if len(set(chosen)) != 4:
        _fail("body_in_channel", "could not place four vertices on the body")
    order = sorted(range(4), key=lambda r: chosen[r])
    body = [chosen[r] for r in order]
    inner = [inner[r] for r in order]
    start = len(vertices)
    vertices += [loop.points[k] for k in body]
    for k in range(4):
        blocks.append((start + k, start + (k + 1) % 4, inner[(k + 1) % 4], inner[k]))
    return _finish_2d(
        store,
        options,
        vertices,
        blocks,
        ids + [p[0] for p in loop.pieces],
        {
            "radial": ((start, inner[0]), _option(options, "radial_points", 21)),
            "walls": [(start + k, start + (k + 1) % 4) for k in range(4)],
        },
    )


def _finish_2d(store, options, vertices, blocks, curve_ids, plan) -> dict:
    """Sketch parameters: assignment, dimensions, walls."""
    points = np.asarray(vertices, dtype=float)
    boundary = _boundary_edges(blocks)
    curves = _assign(store, points, boundary, curve_ids)
    lengths = [float(np.linalg.norm(points[a] - points[b])) for a, b in boundary]
    spacing = options.pop("spacing", None)
    if spacing is None:
        spacing = float(np.median(lengths)) / 20.0
    dimensions = {}
    for key in ("radial", "radial2", "wake"):
        if key in plan and plan[key][1] is not None:
            dimensions[plan[key][0]] = int(plan[key][1])
    dimensions.update(plan.get("extra", {}))
    wall = options.pop("wall", None)
    walls = {}
    if wall is not None:
        for edge in plan.get("walls", ()):
            walls[edge] = float(wall)
        wake = options.pop("wake", None)
        if wake is None:
            wake = 2.0 * float(wall)
        for edge in plan.get("internal", ()):
            walls[edge] = float(wake)
    else:
        options.pop("wake", None)
    # Vertices the template fixes (where sliding has a neutral mode, e.g. rotating along a
    # circle, or drifts); the others relax (10d).
    fixed = sorted(plan.get("fixed", ()))
    return {
        "vertices": tuple(tuple(float(c) for c in p) for p in points),
        "blocks": tuple(tuple(b) for b in blocks),
        "dimensions": dimensions,
        "spacing": float(spacing),
        "curves": curves,
        "walls": walls,
        "fixed": tuple(fixed),
        "distributions": plan.get("distributions", {}),
    }


# ------------------------------------------------------------------------- 3D templates
def _translated(store, curve_id: int, offset) -> int:
    """A copy of a curve moved by ``offset``."""
    d = np.asarray(offset, dtype=float)
    return _mapped(store, curve_id, lambda p: p + d, 1.0)


def _mapped(store, curve_id: int, mapping, scale: float) -> int:
    """A copy of a curve under a similarity ``mapping`` (point -> point; ``scale`` its length
    factor): exact for lines, arcs, polylines, interpolating splines and NURBS."""
    entity = store.get(curve_id)

    def moved(point):
        return tuple(np.asarray(mapping(np.asarray(point, dtype=float)), dtype=float).tolist())

    raw = dict(entity.params)
    if entity.kind == "line":
        raw = {"p0": moved(raw["p0"]), "p1": moved(raw["p1"])}
    elif entity.kind == "arc":
        raw["center"] = moved(raw["center"])
        raw["radius"] = raw["radius"] * scale
    elif entity.kind in ("polyline", "spline"):
        raw["points"] = tuple(moved(q) for q in raw["points"])
    elif entity.kind == "nurbs":
        raw["control_points"] = tuple(moved(q) for q in raw["control_points"])
    return _add(store, entity.kind, raw, f"{entity.name}-end")


def _relaxed_vertices(store, params2d: dict) -> np.ndarray:
    """The 2D template's vertices after its 10d relaxation and group smoothing (interior
    vertices float in the group): a temporary sketch generated, smoothed and removed again
    inside the transaction."""
    from opengrid.model.sketch import Lattice, _evaluators, create_sketch
    from opengrid.model.topology_edit import fresh_name

    name = fresh_name(store, "template-2d")
    create_sketch(store, name=name, **params2d)
    sketch = store.get(name)
    outputs = dict(sketch.params["outputs"])
    lattice = Lattice.build(sketch.params["vertices"], sketch.params["blocks"], name)
    result = _evaluators(store).smooth(outputs["group"])
    points = np.asarray(params2d["vertices"], dtype=float)
    for index, block in enumerate(lattice.blocks):
        grid = result.grids[outputs[f"d{index}"]]
        for corner, (j, i) in zip(block, ((0, 0), (0, -1), (-1, -1), (-1, 0)), strict=True):
            points[corner] = grid[j, i]
    order = {"smoothing_group": 0, "block": 1, "domain": 2, "connector": 3}
    for _role, entity_id in sorted(
        outputs.items(), key=lambda item: order.get(store.get(item[1]).kind, 4)
    ):
        store.remove(entity_id)
    store.remove(sketch.id)
    return points


def _blocks_3d(
    store, params2d: dict, length: float, span_points: int, section=None, body=()
) -> dict:
    """A 3D sketch of six-face blocks: the 2D sketch (relaxed in 2D first) swept along z by
    ``length``, edges at both ends on the curves (copies at the far end), the faces between on
    extruded surfaces of the curves, walls from wall edges, the end planes symmetry planes.
    ``section``: (mapping, scale) - a similarity applied to the ``body`` curves and to the
    vertices not on other curves at the far end (a tapered or swept body; its faces are then
    on ruled surfaces between the two sections). The plain extrusion is ``extrude``."""
    vertices = list(_relaxed_vertices(store, params2d))
    n = len(vertices)
    offset = (0.0, 0.0, length)
    body = set(body)
    other = {v for (a, b), c in params2d["curves"].items() if c not in body for v in (a, b)}
    top = []
    for k, v in enumerate(vertices):
        if section is not None and k not in other:
            top.append(np.asarray(section[0](v), dtype=float) + np.asarray(offset))
        else:
            top.append(v + np.asarray(offset))
    blocks = []
    for b in params2d["blocks"]:
        blocks.append((*b, *(v + n for v in b)))
    surfaces = {}
    curves = {}
    made = {}
    ends = {}
    for (a, b), curve in params2d["curves"].items():
        if curve not in made and section is not None and curve in body:
            mapping, scale = section
            d = np.asarray(offset)
            ends[curve] = _mapped(store, curve, lambda p, m=mapping, d=d: m(p) + d, scale)
            made[curve] = _add(
                store,
                "ruled_surface",
                {"curve0": curve, "curve1": ends[curve]},
                f"{store.get(curve).name}-ruled",
            )
        if curve not in made:
            made[curve] = _add(
                store,
                "extruded_surface",
                {"curve": curve, "direction": offset},
                f"{store.get(curve).name}-swept",
            )
            ends[curve] = _translated(store, curve, offset)
        surfaces[tuple(sorted((a, b, b + n, a + n)))] = made[curve]
        curves[(a, b)] = curve
        curves[(a + n, b + n)] = ends[curve]
    if section is not None:
        # Where two body curves meet (the leading edge), the spanwise edge lies between two
        # ruled surfaces: a straight line from the root to the tip section (exact for a
        # similarity of the sections).
        meets: dict[int, set] = {}
        for (a, b), curve in params2d["curves"].items():
            if curve in body:
                meets.setdefault(a, set()).add(curve)
                meets.setdefault(b, set()).add(curve)
        for v in sorted(meets):
            if len(meets[v]) >= 2 and v not in other:
                curves[(v, v + n)] = _add(
                    store,
                    "line",
                    {
                        "p0": tuple(float(c) for c in vertices[v]),
                        "p1": tuple(float(c) for c in top[v]),
                    },
                    "span-edge",
                )
    walls = {}
    for (a, b), spacing in params2d["walls"].items():
        walls[(a, b, b + n, a + n)] = spacing
    dimensions = dict(params2d["dimensions"])
    dimensions[(0, n)] = int(span_points)
    distributions = {}
    for (a, b), d in params2d.get("distributions", {}).items():
        distributions[(a, b)] = d
        distributions[(a + n, b + n)] = d
    fixed = [v for v in params2d["fixed"]] + [v + n for v in params2d["fixed"]]
    # The end planes are symmetry planes: every section is the same.
    symmetry = [tuple(sorted(b)) for b in params2d["blocks"]]
    symmetry += [tuple(sorted(v + n for v in b)) for b in params2d["blocks"]]
    return {
        "fixed": tuple(fixed),
        "symmetry": tuple(symmetry),
        "vertices": tuple(tuple(float(c) for c in p) for p in vertices + top),
        "blocks": tuple(blocks),
        "dimensions": dimensions,
        "spacing": params2d["spacing"],
        "curves": curves,
        "surfaces": surfaces,
        "walls": walls,
        "distributions": distributions,
    }


def _sphere_in_box(store, options) -> dict:
    center = np.asarray(_option(options, "center", (0.0, 0.0, 0.0)), dtype=float)
    radius = float(_option(options, "radius", 0.5))
    box = float(_option(options, "box", 3.0 * radius))
    cube = float(_option(options, "cube", 0.6 * (radius / math.sqrt(3) + box)))
    n = int(_option(options, "points", 13))
    n_r = int(_option(options, "radial_points", 17))
    n_o = int(_option(options, "outer_points", 9))
    wall = options.pop("wall", None)
    mer_z = _add(
        store,
        "arc",
        {
            "center": tuple(center),
            "radius": radius,
            "normal": (0.0, -1.0, 0.0),
            "x_axis": (1.0, 0.0, 0.0),
            "start_angle": -0.5 * math.pi,
            "end_angle": 0.5 * math.pi,
        },
        "sphere-meridian",
    )
    mer_x = _add(
        store,
        "arc",
        {
            "center": tuple(center),
            "radius": radius,
            "normal": (0.0, 0.0, -1.0),
            "x_axis": (0.0, 1.0, 0.0),
            "start_angle": -0.5 * math.pi,
            "end_angle": 0.5 * math.pi,
        },
        "sphere-meridian",
    )
    axis_point = tuple(center)
    surf = {
        0: _add(
            store,
            "revolved_surface",
            {
                "curve": mer_z,
                "axis_point": axis_point,
                "axis": (0.0, 0.0, 1.0),
                "start_angle": 0.5 * math.pi,
                "end_angle": 2.5 * math.pi,
            },
            "sphere",
        ),
        1: _add(
            store,
            "revolved_surface",
            {
                "curve": mer_z,
                "axis_point": axis_point,
                "axis": (0.0, 0.0, 1.0),
                "end_angle": 2 * math.pi,
            },
            "sphere",
        ),
        2: _add(
            store,
            "revolved_surface",
            {
                "curve": mer_x,
                "axis_point": axis_point,
                "axis": (1.0, 0.0, 0.0),
                "end_angle": 2 * math.pi,
            },
            "sphere",
        ),
    }
    coords = (-box, -cube, cube, box)
    vertices = []
    index = {}
    for k, j, i in itertools.product(range(4), repeat=3):
        index[(i, j, k)] = len(vertices)
        vertices.append(center + np.array([coords[i], coords[j], coords[k]]))
    on_sphere = {}
    for c in itertools.product((-1, 1), repeat=3):
        on_sphere[c] = len(vertices)
        vertices.append(center + radius * np.array(c) / math.sqrt(3))
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
            faces[tuple(sorted(on_sphere[c] for c in cs))] = surf[axis]
            for m in range(4):
                edge = tuple(sorted((on_sphere[cs[m]], on_sphere[cs[(m + 1) % 4]])))
                if edge in arcs:
                    continue
                pa = vertices[edge[0]] - center
                pb = vertices[edge[1]] - center
                normal = np.cross(pa, pb)
                normal /= np.linalg.norm(normal)
                x_axis = pa / np.linalg.norm(pa)
                sweep = math.acos(float(np.clip(pa @ pb / (radius * radius), -1.0, 1.0)))
                arcs[edge] = _add(
                    store,
                    "arc",
                    {
                        "center": tuple(center),
                        "radius": radius,
                        "normal": tuple(normal),
                        "x_axis": tuple(x_axis),
                        "start_angle": 0.0,
                        "end_angle": sweep,
                    },
                    "sphere-edge",
                )
    dimensions = {next(iter(arcs)): n, (blocks[-1][0], blocks[-1][4]): n_r}
    for axis in range(3):
        for segment, count in ((0, n_o), (1, n), (2, n_o)):
            lo = [0, 0, 0]
            hi = [0, 0, 0]
            lo[axis] = segment
            hi[axis] = segment + 1
            dimensions[(index[tuple(lo)], index[tuple(hi)])] = count
    walls = {}
    if wall is not None:
        walls = dict.fromkeys(faces, float(wall))
    return {
        "vertices": tuple(tuple(float(c) for c in p) for p in vertices),
        "blocks": tuple(blocks),
        "dimensions": dimensions,
        "curves": arcs,
        "surfaces": faces,
        "walls": walls,
    }


BUILDERS_2D = {
    "ogrid": _ogrid,
    "butterfly": _butterfly,
    "cgrid": _cgrid,
    "hgrid": _hgrid,
    "body_in_channel": _body_in_channel,
}


def create_template(store, kind: str, name=None, **options):
    """Build a template's geometry and sketch (inside the open transaction). Returns
    ([sketch name], warnings)."""
    from opengrid.model.sketch import _canonical_input, create_sketch
    from opengrid.model.topology_edit import fresh_name

    options = {k: _canonical_input(v) for k, v in options.items()}
    if kind not in KINDS_2D and kind not in KINDS_3D:
        _fail("create_template", f"unknown template {kind!r} ({', '.join(KINDS_2D + KINDS_3D)})")
    relax = options.pop("relax", True)
    layers = options.pop("layers", 0)
    growth = options.pop("growth", 1.2)
    if layers != 0 and "wall" in options and "wake" not in options:
        # A wake cut with layers takes the wall's spacing: twice it (the default without
        # layers) steps the layer front by 2x where the cut meets the body.
        options["wake"] = options["wall"]
    if kind == "sphere_in_box":
        params = _sphere_in_box(store, options)
    elif kind in EXTRUDED:
        length = options.pop("length", None)
        span_points = int(options.pop("span_points", 11))
        span_distribution = options.pop("span_distribution", "uniform")
        three_d_blocks = bool(options.pop("three_d_blocks", False))
        taper = float(options.pop("taper", 1.0))
        sweep = float(options.pop("sweep", 0.0))
        if kind != "wing_in_box" and (taper != 1.0 or sweep != 0.0):
            _fail(kind, "taper and sweep are options of wing_in_box")
        params = BUILDERS_2D[EXTRUDED[kind]](store, options)
        body = params.pop("_body", ())
        leading = np.asarray(params.pop("_leading_edge", (0.0, 0.0, 0.0)))
        chord = np.asarray(params.pop("_chord", (1.0, 0.0, 0.0)))
        if length is None:
            points = np.asarray(params["vertices"])
            length = float(np.max(np.ptp(points[:, :2], axis=0))) / 4.0
        section = None
        if taper != 1.0 or sweep != 0.0:
            # The tip section: the root scaled by taper about the leading edge, moved back by
            # the sweep (degrees of the leading edge).
            shift = float(length) * math.tan(math.radians(sweep)) * chord

            def mapping(p, leading=leading, taper=taper, shift=shift):
                return leading + taper * (np.asarray(p) - leading) + shift

            section = (mapping, taper)
            three_d_blocks = True
        if three_d_blocks:
            params = _blocks_3d(
                store, params, float(length), span_points, section=section, body=body
            )
            relax = False  # relaxed in 2D already
        else:
            # The smoothed 2D grid extruded along z: the same in every section.
            params["extrude"] = {
                "vector": (0.0, 0.0, float(length)),
                "points": span_points,
                "distribution": span_distribution,
            }
    else:
        params = BUILDERS_2D[kind](store, options)
    # Builder notes for the 3D variants (body curves, leading edge, chord) are not sketch
    # parameters.
    params = {k: v for k, v in params.items() if not k.startswith("_")}
    if layers != 0:
        # Wall layers on every wall (the wake cut of a C-grid too).
        params["walls"] = {
            key: {"spacing": spacing, "layers": layers, "growth": growth}
            for key, spacing in params["walls"].items()
        }
    if options:
        _fail("create_template", f"unexpected options for {kind}: {', '.join(sorted(options))}")
    if name is None:
        name = fresh_name(store, kind)
    return create_sketch(store, name=name, relax=relax, **params)


__all__ = ["KINDS_2D", "KINDS_3D", "create_template"]
