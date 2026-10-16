"""Structured domains from a network of planar curves (strokes drawn across each other).

Method:
1. Vertices: the curves' end points, their sharp corners (turning angle of the tessellation
   above ``corner_angle``), points where an end lies on another curve, and crossings - found
   on the tessellations (segment-segment tests), then refined by Newton's method on
   c_a(t_a) = c_b(t_b) in the plane. Points closer than the tolerance are one vertex.
2. Edges: the pieces of each curve between consecutive vertices.
3. Faces: the bounded faces of the planar arrangement, traced with the half-edge rule "leave
   each vertex along the next edge clockwise from the one you arrived by" (the face lies on the
   left); faces with negative area (the outer boundary) are dropped.
4. Corners of a face: vertices where its boundary turns by more than ``corner_angle``. A
   face with four corners becomes a structured domain whose sides are the edge chains between
   them. Other faces (e.g. a smooth closed curve, no corners) are reported, or become
   unstructured domains on request.
5. Point counts: each edge gets length / spacing; edges on opposite single-edge sides of a
   face are made equal (union-find, the largest count of each class), so the domains are
   conforming; a face whose opposite sides are chains of several edges keeps its counts and
   is reported if they do not match.
The connectors are pieces (t0, t1) of the original curves, so the curves are not modified.
"""

import itertools
import math
from dataclasses import dataclass, field

import numpy as np

from opengrid.validation import Issue, ValidationError

TESSELLATION = 2e-4


@dataclass
class Edge:
    curve: int  # entity id
    t0: float
    t1: float
    v0: int
    v1: int
    length: float
    points: np.ndarray  # tessellation from t0 to t1


@dataclass
class Network:
    vertices: list[np.ndarray] = field(default_factory=list)
    edges: list[Edge] = field(default_factory=list)
    faces: list[list[tuple[int, bool]]] = field(default_factory=list)  # (edge, forward)
    issues: list[Issue] = field(default_factory=list)


def _xy(p) -> np.ndarray:
    return np.asarray(p, dtype=float)[..., :2]


def _turn(a, b, c) -> float:
    u = _xy(b) - _xy(a)
    v = _xy(c) - _xy(b)
    return math.degrees(math.atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1]))


def _segment_hits(pa: np.ndarray, pb: np.ndarray):
    """(i, j, s, r): segment i of polyline pa crosses segment j of pb at fractions s, r."""
    a0 = _xy(pa[:-1])[:, None, :]
    a1 = _xy(pa[1:])[:, None, :]
    b0 = _xy(pb[:-1])[None, :, :]
    b1 = _xy(pb[1:])[None, :, :]
    da = a1 - a0
    db = b1 - b0
    denom = da[..., 0] * db[..., 1] - da[..., 1] * db[..., 0]
    w = b0 - a0
    with np.errstate(divide="ignore", invalid="ignore"):
        s = (w[..., 0] * db[..., 1] - w[..., 1] * db[..., 0]) / denom
        r = (w[..., 0] * da[..., 1] - w[..., 1] * da[..., 0]) / denom
    hit = (denom != 0) & (s >= 0) & (s <= 1) & (r >= 0) & (r <= 1)
    ii, jj = np.nonzero(hit)
    return [(int(i), int(j), float(s[i, j]), float(r[i, j])) for i, j in zip(ii, jj, strict=True)]


def _refine(ca, cb, ta: float, tb: float) -> tuple[float, float]:
    """Newton on ca(ta) - cb(tb) = 0 (xy), finite-difference derivatives."""
    h = 1e-7
    for _ in range(30):
        f = _xy(ca.eval(ta)) - _xy(cb.eval(tb))
        if float(np.linalg.norm(f)) < 1e-15:
            break
        da = (_xy(ca.eval(min(ta + h, 1.0))) - _xy(ca.eval(max(ta - h, 0.0)))) / (
            min(ta + h, 1.0) - max(ta - h, 0.0)
        )
        db = (_xy(cb.eval(min(tb + h, 1.0))) - _xy(cb.eval(max(tb - h, 0.0)))) / (
            min(tb + h, 1.0) - max(tb - h, 0.0)
        )
        jac = np.array([[da[0], -db[0]], [da[1], -db[1]]])
        if abs(np.linalg.det(jac)) < 1e-300:
            break
        step = np.linalg.solve(jac, -f)
        ta = float(np.clip(ta + step[0], 0.0, 1.0))
        tb = float(np.clip(tb + step[1], 0.0, 1.0))
    return ta, tb


def build_network(store, curve_ids, corner_angle: float = 30.0) -> Network:
    curves = {c: store.geometry(c) for c in curve_ids}
    params = {c: np.asarray(curves[c].tessellate(TESSELLATION)) for c in curve_ids}
    points = {c: np.asarray(curves[c].eval_many(params[c])) for c in curve_ids}
    all_pts = np.concatenate(list(points.values()))
    scale = float(np.max(np.ptp(all_pts, axis=0))) or 1.0
    tol = 1e-7 * scale
    net = Network()
    z = float(all_pts[0, 2])
    if float(np.ptp(all_pts[:, 2])) > 1e-9 * scale:
        raise ValidationError([Issue("curve network", "curves must lie in a plane z = const")])
    # Cut parameters per curve: ends, sharp corners, crossings and ends on other curves.
    cuts: dict[int, list[float]] = {c: [0.0, 1.0] for c in curve_ids}
    for c in curve_ids:
        p = points[c]
        for k in range(1, len(p) - 1):
            if abs(_turn(p[k - 1], p[k], p[k + 1])) > corner_angle:
                cuts[c].append(float(params[c][k]))
    ids = list(curve_ids)
    for a_index, a in enumerate(ids):
        for b in ids[a_index + 1 :]:
            for i, j, s, r in _segment_hits(points[a], points[b]):
                ta = params[a][i] + s * (params[a][i + 1] - params[a][i])
                tb = params[b][j] + r * (params[b][j + 1] - params[b][j])
                ta, tb = _refine(curves[a], curves[b], ta, tb)
                if np.linalg.norm(_xy(curves[a].eval(ta)) - _xy(curves[b].eval(tb))) > tol:
                    continue
                cuts[a].append(ta)
                cuts[b].append(tb)
    for a in ids:
        for end in (0.0, 1.0):
            p = curves[a].eval(end)
            for b in ids:
                if b == a:
                    continue
                t, _q, d = curves[b].project(tuple(p))
                if d <= tol and 0.0 < t < 1.0:
                    cuts[b].append(float(t))

    # Vertices (merged within the tolerance) and edges.
    def vertex(point) -> int:
        q = np.asarray(point, dtype=float)
        for k, v in enumerate(net.vertices):
            if np.linalg.norm(v - q) <= 1e3 * tol:
                return k
        net.vertices.append(q)
        return len(net.vertices) - 1

    for c in ids:
        ts = sorted(cuts[c])
        unique = [ts[0]]
        for t in ts[1:]:
            if t - unique[-1] > 1e-9:
                unique.append(t)
        for t0, t1 in itertools.pairwise(unique):
            v0 = vertex(curves[c].eval(t0))
            v1 = vertex(curves[c].eval(t1))
            sample = np.asarray(curves[c].eval_many(np.linspace(t0, t1, 33)))
            if v0 == v1 and float(np.linalg.norm(sample[16] - sample[0])) <= 1e3 * tol:
                continue  # degenerate
            length = float(np.sum(np.linalg.norm(np.diff(sample, axis=0), axis=1)))
            net.edges.append(Edge(c, t0, t1, v0, v1, length, sample))
    # Faces: half-edge traversal.
    out: dict[int, list[tuple[float, int, bool]]] = {}
    for e, edge in enumerate(net.edges):
        for forward in (True, False):
            pts = edge.points
            if not forward:
                pts = pts[::-1]
            start = edge.v0
            if not forward:
                start = edge.v1
            d = _xy(pts[1]) - _xy(pts[0])
            out.setdefault(start, []).append((math.atan2(d[1], d[0]), e, forward))
    for v in out:
        out[v].sort()
    used = set()
    for e in range(len(net.edges)):
        for forward in (True, False):
            if (e, forward) in used:
                continue
            cycle = []
            current = (e, forward)
            while current not in used:
                used.add(current)
                cycle.append(current)
                edge = net.edges[current[0]]
                end = edge.v1
                if not current[1]:
                    end = edge.v0
                # Arrived at `end` along the reverse of the outgoing half-edge (current[0],
                # not current[1]); take the next one clockwise from it.
                around = out[end]
                index = around.index(
                    next(item for item in around if item[1] == current[0] and item[2] != current[1])
                )
                nxt = around[(index - 1) % len(around)]
                current = (nxt[1], nxt[2])
            polygon = np.concatenate([_half(net, h)[:-1] for h in cycle])
            x = polygon[:, 0]
            y = polygon[:, 1]
            area = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
            if area > 0:
                net.faces.append(cycle)
    del z
    return net


def _half(net: Network, half) -> np.ndarray:
    edge = net.edges[half[0]]
    if half[1]:
        return edge.points
    return edge.points[::-1]


def face_sides(net: Network, face, corner_angle: float = 30.0):
    """The face's sides (lists of half-edges) between its corners, or None."""
    n = len(face)
    corners = []
    for k in range(n):
        incoming = _half(net, face[k - 1])
        outgoing = _half(net, face[k])
        if abs(_turn(incoming[-2], incoming[-1], outgoing[1])) > corner_angle:
            corners.append(k)
    if len(corners) != 4:
        return None, len(corners)
    sides = []
    for m in range(4):
        a = corners[m]
        b = corners[(m + 1) % 4]
        side = []
        k = a
        while True:
            side.append(face[k])
            k = (k + 1) % n
            if k == b:
                break
        sides.append(side)
    return sides, 4


class _UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, a):
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def plan_domains(store, curve_ids, spacing: float, corner_angle: float = 30.0, unstructured=False):
    """Entities to create for the network: (items for CreateEntities, issues).

    Items: one connector per edge used by a face, one structured domain per four-sided face
    (unstructured domains for the others when ``unstructured``)."""
    net = build_network(store, curve_ids, corner_angle)
    issues = list(net.issues)
    if not net.faces:
        issues.append(Issue("curve network", "the curves enclose no region"))
        return [], issues, net
    counts = [max(2, round(e.length / spacing) + 1) for e in net.edges]
    uf = _UnionFind(len(net.edges))
    structured = []
    others = []
    for f, face in enumerate(net.faces):
        sides, corners = face_sides(net, face, corner_angle)
        if sides is None:
            others.append((f, corners))
            continue
        structured.append((f, sides))
        for m in (0, 1):
            a, b = sides[m], sides[m + 2]
            if len(a) == 1 and len(b) == 1:
                uf.union(a[0][0], b[0][0])
    best: dict[int, int] = {}
    for e in range(len(net.edges)):
        r = uf.find(e)
        best[r] = max(best.get(r, 2), counts[e])
    counts = [best[uf.find(e)] for e in range(len(net.edges))]
    names = {}
    items = []
    prefix = 1
    while any(entity.name.startswith(f"net{prefix}-") for entity in store):
        prefix += 1

    def connector(e: int) -> str:
        if e not in names:
            edge = net.edges[e]
            names[e] = f"net{prefix}-c{len(names) + 1}"
            items.append(
                (
                    "connector",
                    {
                        "segments": ({"curve": edge.curve, "t0": edge.t0, "t1": edge.t1},),
                        "dimension": counts[e],
                    },
                    names[e],
                )
            )
        return names[e]

    for f, sides in structured:
        lengths = [sum(counts[h[0]] - 1 for h in side) for side in sides]
        if lengths[0] != lengths[2] or lengths[1] != lengths[3]:
            issues.append(
                Issue(
                    f"face {f}",
                    f"opposite sides have {lengths[0]}/{lengths[2]} and {lengths[1]}/{lengths[3]} "
                    "intervals (sides made of several edges); adjust the connectors",
                )
            )
            continue
        edges = []
        for side in sides:
            edges.append(
                tuple({"connector": connector(e), "reversed": not forward} for e, forward in side)
            )
        items.append(("domain", {"edges": tuple(edges)}, f"net{prefix}-d{f + 1}"))
    for f, corners in others:
        if unstructured:
            loop = tuple(
                {"connector": connector(e), "reversed": not forward} for e, forward in net.faces[f]
            )
            items.append(("unstructured_domain", {"loops": (loop,)}, f"net{prefix}-u{f + 1}"))
            continue
        issues.append(Issue(f"face {f}", f"has {corners} corner(s), not 4 (no structured domain)"))
    return items, issues, net
