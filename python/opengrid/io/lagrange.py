"""Node sets of straight-sided Lagrange elements of any order and their orderings in Gmsh,
CGNS and VTK (Phase 11f).

A straight-sided element of order p is the linear element (straight edges, planar triangles,
bilinear quadrilaterals) with the equispaced Lagrange nodes of order p placed by its own
geometric map: barycentric for triangles and tetrahedra, bilinear / trilinear for
quadrilaterals and hexahedra, triangle x segment for prisms, and the collapsed map for
pyramids (the base square shrinks linearly to the apex). Nodes are generated as exact
rational coordinates in one canonical reference element per kind (VTK's parametric
elements: unit simplex, unit square / cube, prism = unit triangle x [0, 1], pyramid = unit
square base with apex (1/2, 1/2, 1)), so the orderings of the three formats can be compared
and every node maps to weights over the element's vertices (the geometry, and the key that
merges the nodes shared by neighbouring elements).

The orderings follow one pattern - vertices, then the interior nodes of each edge (from its
first vertex to its second), of each face (in the face's own vertex order) and of the volume -
with each format's own edge directions, face vertex orders and interior modes written as data
below (FORMATS). Interior modes: "recursive" (the interior lattice as a lower-order element of
the same kind), "grid" (rows, the first direction fastest), "ring" (quadrilaterals: ring by
ring from the outside in), "layers" (hexahedra, prisms: the interior of each layer from the
bottom) and "columns" (prisms: each triangle node's column, its two ends first). The tables
are checked in the tests against Gmsh's own element definitions, VTK's vtkLagrange* cells and
CGNS through Gmsh's CGNS writer and the SIDS node descriptions (complete elements up to
order 4).

References: Gmsh reference manual, "Node ordering"; CGNS SIDS, section 3.3 "Element
Connectivity"; W. J. Schroeder et al., "Definition of higher-order Lagrange cells in VTK"
(Kitware, 2018).
"""

from fractions import Fraction
from functools import cache

import numpy as np

F = Fraction  # exact reference coordinates

KINDS = ("line", "triangle", "quad", "tetra", "hexa", "prism", "pyramid")
DIM = {"line": 1, "triangle": 2, "quad": 2, "tetra": 3, "hexa": 3, "prism": 3, "pyramid": 3}

# Canonical reference vertices, in the vertex order shared by VTK, Gmsh and CGNS.
VERTICES = {
    "line": ((0,), (1,)),
    "triangle": ((0, 0), (1, 0), (0, 1)),
    "quad": ((0, 0), (1, 0), (1, 1), (0, 1)),
    "tetra": ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
    "hexa": (
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    ),
    "prism": ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)),
    "pyramid": ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (F(1, 2), F(1, 2), 1)),
}


def _vec(point) -> tuple:
    return tuple(F(c) for c in point)


def _lerp(a, b, t) -> tuple:
    return tuple(x + (y - x) * t for x, y in zip(a, b, strict=True))


# ------------------------------------------------------------------ patterns
def _segment_interior(a, b, p) -> list:
    return [_lerp(a, b, F(k, p)) for k in range(1, p)]


def _centroid(corners) -> tuple:
    return tuple(sum(x) / len(corners) for x in zip(*corners, strict=True))


def _bilinear(corners, u, v) -> tuple:
    a, b, c, d = corners
    return tuple(
        (1 - u) * (1 - v) * x + u * (1 - v) * y + u * v * z + (1 - u) * v * w
        for x, y, z, w in zip(a, b, c, d, strict=True)
    )


def _inner_triangle(corners, p) -> list:
    """Corners of the interior lattice (an order p - 3 triangle) of an order-p triangle."""
    s = F(1, p)
    out = []
    for k in range(3):
        v, u, w = corners[k], corners[(k + 1) % 3], corners[(k + 2) % 3]
        out.append(tuple(x + (y - x) * s + (z - x) * s for x, y, z in zip(v, u, w, strict=True)))
    return out


def _inner_quad(corners, p) -> list:
    """Corners of the interior lattice (an order p - 2 quadrilateral) of an order-p one."""
    s = F(1, p)
    return [_bilinear(corners, u, v) for u, v in ((s, s), (1 - s, s), (1 - s, 1 - s), (s, 1 - s))]


def _polygon(kind, corners, p, fmt) -> list:
    """All nodes of an order-p triangle / quadrilateral on the given corners."""
    corners = [_vec(c) for c in corners]
    if p == 0:
        return [_centroid(corners)]
    spec = fmt[kind]
    out = list(corners)
    for i, j in spec["edges"]:
        out += _segment_interior(corners[i], corners[j], p)
    return out + _polygon_interior(kind, corners, p, fmt, spec["interior"])


def _polygon_interior(kind, corners, p, fmt, mode) -> list:
    if kind == "triangle":
        if p < 3:
            return []
        if mode == "grid":
            a, b, c = corners
            return [
                tuple(
                    x + (y - x) * F(i, p) + (z - x) * F(j, p)
                    for x, y, z in zip(a, b, c, strict=True)
                )
                for j in range(1, p - 1)
                for i in range(1, p - j)
            ]
        return _polygon(kind, _inner_triangle(corners, p), p - 3, fmt)
    if p < 2:
        return []
    if mode == "grid":
        return [_bilinear(corners, F(i, p), F(j, p)) for j in range(1, p) for i in range(1, p)]
    inner = _inner_quad(corners, p)
    if mode == "ring":
        return _rings(inner, p - 2)
    return _polygon(kind, inner, p - 2, fmt)


def _rings(corners, q) -> list:
    """An order-q quadrilateral lattice ring by ring from the outside in, each ring
    counter-clockwise from its first corner."""
    if q == 0:
        return [_centroid(corners)]
    out = []
    for k in range(4):
        a, b = corners[k], corners[(k + 1) % 4]
        out += [a, *_segment_interior(a, b, q)]
    if q >= 2:
        out += _rings(_inner_quad(corners, q), q - 2)
    return out


def _inner_tetra(corners, p) -> list:
    s = F(1, p)
    out = []
    for k in range(4):
        others = [corners[m] for m in range(4) if m != k]
        out.append(tuple(x + sum((o[i] - x) * s for o in others) for i, x in enumerate(corners[k])))
    return out


def _inner_hexa(corners, p) -> list:
    s = F(1, p)
    out = []
    for w in (s, 1 - s):
        bottom = [_lerp(corners[k], corners[k + 4], w) for k in range(4)]
        out += _inner_quad(bottom, p)
    return out


def _inner_pyramid(corners, p) -> list:
    """Corners of the inner pyramid (order p - 3) of an order-p pyramid: its base is the
    interior of the first layer above the base (lattice indices 1 .. p - 2 of p - 1), its
    apex the centre of layer p - 2."""
    base = corners[:4]
    level = [_lerp(b, corners[4], F(1, p)) for b in base]
    n = p - 1
    out = [
        _bilinear(level, F(u, n), F(v, n))
        for u, v in ((1, 1), (n - 1, 1), (n - 1, n - 1), (1, n - 1))
    ]
    out.append(_lerp(_centroid(base), corners[4], F(p - 2, p)))
    return out


def _element(kind, corners, p, fmt) -> list:
    """All nodes of an order-p element of the given kind on the given corners, in the order
    `fmt` (one of the format tables below) defines."""
    if DIM[kind] < 3:
        if kind == "line":
            a, b = (_vec(c) for c in corners)
            return [a, b, *_segment_interior(a, b, p)]
        return _polygon(kind, corners, p, fmt)
    corners = [_vec(c) for c in corners]
    if p == 0:
        return [_centroid(corners)]
    spec = fmt[kind]
    out = list(corners)
    for i, j in spec["edges"]:
        out += _segment_interior(corners[i], corners[j], p)
    for face in spec["faces"]:
        face_kind = "quad"
        if len(face) == 3:
            face_kind = "triangle"
        mode = spec_face_mode(fmt, kind, face_kind)
        out += _polygon_interior(face_kind, [corners[v] for v in face], p, fmt, mode)
    return out + _volume_interior(kind, corners, p, fmt, spec["interior"])


def spec_face_mode(fmt, kind, face_kind) -> str:
    return fmt[kind].get("face_interior", {}).get(face_kind, fmt[face_kind]["interior"])


def _volume_interior(kind, corners, p, fmt, mode) -> list:
    if kind == "tetra":
        if p < 4:
            return []
        return _element(kind, _inner_tetra(corners, p), p - 4, fmt)
    if kind == "pyramid":
        if p < 3:
            return []
        return _element(kind, _inner_pyramid(corners, p), p - 3, fmt)
    if kind == "hexa":
        if p < 2:
            return []
        if mode == "grid":
            return [
                _lerp(
                    _bilinear(corners[:4], F(i, p), F(j, p)),
                    _bilinear(corners[4:], F(i, p), F(j, p)),
                    F(k, p),
                )
                for k in range(1, p)
                for j in range(1, p)
                for i in range(1, p)
            ]
        if mode == "layers":
            out = []
            quad_mode = spec_face_mode(fmt, kind, "quad")
            for k in range(1, p):
                layer = [_lerp(corners[m], corners[m + 4], F(k, p)) for m in range(4)]
                out += _polygon_interior("quad", layer, p, fmt, quad_mode)
            return out
        return _element(kind, _inner_hexa(corners, p), p - 2, fmt)
    # prism: the interior triangle lattice of each layer, by layer ("layers") or the
    # triangle node outer and the height inner with the two end layers first ("columns").
    if p < 3:
        return []
    heights = [F(k, p) for k in range(1, p)]
    triangle = _polygon_interior(
        "triangle", corners[:3], p, fmt, spec_face_mode(fmt, kind, "triangle")
    )
    lift = tuple(c - b for c, b in zip(corners[3], corners[0], strict=True))

    def at(point, w):
        return tuple(x + d * w for x, d in zip(point, lift, strict=True))

    if mode == "layers":
        return [at(t, w) for w in heights for t in triangle]
    ends_first = [heights[0], heights[-1], *heights[1:-1]]
    return [at(t, w) for t in triangle for w in ends_first]


# ------------------------------------------------------------------ format definitions
# Edges (from, to), faces (vertex lists, in each face's own order) and interior modes per
# format and kind.
GMSH = {
    "triangle": {"edges": ((0, 1), (1, 2), (2, 0)), "interior": "recursive"},
    "quad": {"edges": ((0, 1), (1, 2), (2, 3), (3, 0)), "interior": "recursive"},
    "tetra": {
        "edges": ((0, 1), (1, 2), (2, 0), (3, 0), (3, 2), (3, 1)),
        "faces": ((0, 2, 1), (0, 1, 3), (0, 3, 2), (3, 1, 2)),
        "interior": "recursive",
    },
    "hexa": {
        "edges": (
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 5),
            (2, 3),
            (2, 6),
            (3, 7),
            (4, 5),
            (4, 7),
            (5, 6),
            (6, 7),
        ),
        "faces": (
            (0, 3, 2, 1),
            (0, 1, 5, 4),
            (0, 4, 7, 3),
            (1, 2, 6, 5),
            (2, 3, 7, 6),
            (4, 5, 6, 7),
        ),
        "interior": "recursive",
    },
    "prism": {
        "edges": ((0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (3, 5), (4, 5)),
        "faces": ((0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (0, 3, 5, 2), (1, 2, 5, 4)),
        "interior": "columns",
    },
    "pyramid": {
        "edges": ((0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4)),
        "faces": ((0, 1, 4), (3, 0, 4), (1, 2, 4), (2, 3, 4), (0, 3, 2, 1)),
        "interior": "recursive",
    },
}

VTK = {
    "triangle": {"edges": ((0, 1), (1, 2), (2, 0)), "interior": "recursive"},
    "quad": {"edges": ((0, 1), (1, 2), (3, 2), (0, 3)), "interior": "grid"},
    "tetra": {
        "edges": ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),
        "faces": ((0, 1, 3), (2, 3, 1), (0, 3, 2), (0, 2, 1)),
        "interior": "recursive",
    },
    "hexa": {
        "edges": (
            (0, 1),
            (1, 2),
            (3, 2),
            (0, 3),
            (4, 5),
            (5, 6),
            (7, 6),
            (4, 7),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ),
        "faces": (
            (0, 3, 7, 4),
            (1, 2, 6, 5),
            (0, 1, 5, 4),
            (3, 2, 6, 7),
            (0, 1, 2, 3),
            (4, 5, 6, 7),
        ),
        "interior": "grid",
    },
    "prism": {
        "edges": ((0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)),
        "faces": ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
        "face_interior": {"triangle": "grid", "quad": "grid"},
        "interior": "layers",
    },
}

CGNS = {
    "triangle": {"edges": ((0, 1), (1, 2), (2, 0)), "interior": "recursive"},
    "quad": {"edges": ((0, 1), (1, 2), (2, 3), (3, 0)), "interior": "ring"},
    "tetra": {
        "edges": ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),
        "faces": ((0, 1, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3)),
        "interior": "recursive",
    },
    "hexa": {
        "edges": (
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
        ),
        "faces": (
            (0, 1, 2, 3),
            (0, 1, 5, 4),
            (1, 2, 6, 5),
            (2, 3, 7, 6),
            (3, 0, 4, 7),
            (4, 5, 6, 7),
        ),
        "interior": "layers",
    },
    # PENTA_40 / PENTA_75 from the SIDS node descriptions (Gmsh does not write them).
    "prism": {
        "edges": ((0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5), (5, 3)),
        "faces": ((0, 1, 2), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5), (3, 4, 5)),
        "interior": "layers",
    },
    "pyramid": {
        "edges": ((0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)),
        "faces": ((0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
        "interior": "recursive",
    },
}

FORMATS = {"gmsh": GMSH, "vtk": VTK, "cgns": CGNS}


@cache
def nodes(fmt: str, kind: str, p: int) -> tuple:
    """Canonical reference coordinates (exact fractions) of the order-p nodes in the order of
    format `fmt` ("gmsh", "vtk", "cgns")."""
    return tuple(_element(kind, VERTICES[kind], p, FORMATS[fmt]))


@cache
def _permutation(fmt: str, kind: str, p: int) -> tuple:
    canonical = {x: k for k, x in enumerate(nodes("gmsh", kind, p))}
    return tuple(canonical[x] for x in nodes(fmt, kind, p))


def permutation(fmt: str, kind: str, p: int) -> np.ndarray:
    """perm such that conn[:, perm] lists an element's nodes (held in Gmsh's order) in the
    order of format ``fmt``."""
    if fmt not in FORMATS or (kind not in FORMATS[fmt] and kind != "line"):
        raise ValueError(f"{fmt} has no order-{p} {kind} elements")
    return np.array(_permutation(fmt, kind, p), dtype=np.int64)


# Element type identifiers per kind, orders 1 .. 6 (Gmsh), the SIDS names up to order 4 (CGNS)
# and the VTK Lagrange cell types (VTK 9.7 has no Lagrange pyramid class).
GMSH_TYPES = {
    "line": (1, 8, 26, 27, 28, 62),
    "triangle": (2, 9, 21, 23, 25, 42),
    "quad": (3, 10, 36, 37, 38, 47),
    "tetra": (4, 11, 29, 30, 31, 71),
    "hexa": (5, 12, 92, 93, 94, 95),
    "prism": (6, 13, 90, 91, 106, 107),
    "pyramid": (7, 14, 118, 119, 120, 121),
}
CGNS_PREFIX = {
    "line": "BAR",
    "triangle": "TRI",
    "quad": "QUAD",
    "tetra": "TETRA",
    "pyramid": "PYRA",
    "prism": "PENTA",
    "hexa": "HEXA",
}
CGNS_MAX_ORDER = 4
VTK_LAGRANGE = {"line": 68, "triangle": 69, "quad": 70, "tetra": 71, "hexa": 72, "prism": 73}


def gmsh_type(kind: str, p: int) -> int:
    return GMSH_TYPES[kind][p - 1]


def cgns_type(kind: str, p: int) -> str:
    if p > CGNS_MAX_ORDER:
        raise ValueError(f"CGNS defines elements up to order {CGNS_MAX_ORDER}, not {p}")
    return f"{CGNS_PREFIX[kind]}_{count(kind, p)}"


def weights(kind: str, point) -> tuple:
    """Exact values of the linear element's shape functions (one per vertex) at a canonical
    reference point: the node's position as a combination of the element's vertices."""
    point = _vec(point)
    if kind == "line":
        (u,) = point
        return (1 - u, u)
    if kind == "triangle":
        u, v = point
        return (1 - u - v, u, v)
    if kind == "quad":
        u, v = point
        return ((1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v)
    if kind == "tetra":
        u, v, w = point
        return (1 - u - v - w, u, v, w)
    if kind == "hexa":
        u, v, w = point
        base = ((1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v)
        return tuple(b * (1 - w) for b in base) + tuple(b * w for b in base)
    if kind == "prism":
        u, v, w = point
        base = (1 - u - v, u, v)
        return tuple(b * (1 - w) for b in base) + tuple(b * w for b in base)
    u, v, w = point  # pyramid: the bilinear map of the layer square at height w
    if w == 1:
        return (0, 0, 0, 0, F(1))
    s = (u - w / 2) / (1 - w)
    t = (v - w / 2) / (1 - w)
    base = ((1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t)
    return (*(b * (1 - w) for b in base), w)


def count(kind: str, p: int) -> int:
    if kind == "line":
        return p + 1
    if kind == "triangle":
        return (p + 1) * (p + 2) // 2
    if kind == "quad":
        return (p + 1) ** 2
    if kind == "tetra":
        return (p + 1) * (p + 2) * (p + 3) // 6
    if kind == "hexa":
        return (p + 1) ** 3
    if kind == "prism":
        return (p + 1) * (p + 1) * (p + 2) // 2
    return (p + 1) * (p + 2) * (2 * p + 3) // 6  # pyramid


def as_array(nodes) -> np.ndarray:
    return np.array([[float(c) for c in n] for n in nodes])
