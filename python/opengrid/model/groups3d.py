"""Smoothing groups of blocks (3D): six-face blocks smoothed together as one hexahedral grid.

The member block grids (their own TFI, without their own smoothing) are merged into one node
set: block faces are copies of their domain grids, so shared nodes have bitwise identical
coordinates and are merged by exact coordinates (as the 3D topology matches interfaces). The
solver is ``opengrid._core.smooth_multiblock_3d``.

Roles of the member blocks' faces (domains): "interface" (two member blocks), "outside" (also
used by a block or volume outside the group: fixed), "wall" (fixed, wall control), "sliding"
(nodes slide on the face's surface, or on its plane when it is planar; a curved face without a
surface stays fixed: "fixed (curved)"), "fixed" (also faces listed in ``pinned``). Nodes of wall
faces on a surface are projected onto it. A wall between two member blocks is an internal wall:
fixed, with wall control on both sides. Nodes on boundary edges: fixed when a fixed,
wall or outside face holds them, otherwise on the surface shared by their faces, or along their
connector between two sliding faces on different surfaces; connector end nodes are fixed.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from opengrid import _core
from opengrid.model.blocks import FACE_AXIS, FACE_NAMES
from opengrid.validation import Issue, ValidationError

PLANE_TOLERANCE = 1e-9


def _clean(points: np.ndarray) -> np.ndarray:
    """Coordinates with -0.0 turned into 0.0 (so equal values have equal bytes)."""
    return np.ascontiguousarray(points, dtype=float) + 0.0


def face_ids(ids: np.ndarray, face: str) -> np.ndarray:
    """Node ids (nk, nj, ni) of a block face, in the layout of blocks.face_array."""
    if face == "imin":
        return ids[:, :, 0]
    if face == "imax":
        return ids[:, :, -1]
    if face == "jmin":
        return ids[:, 0, :]
    if face == "jmax":
        return ids[:, -1, :]
    if face == "kmin":
        return ids[0]
    return ids[-1]


def _plane(grid: np.ndarray):
    """(PlaneSurface, key) for a planar face grid (or points), else (None, None)."""
    points = grid.reshape(-1, 3)
    centre = points.mean(axis=0)
    _u, _singular, vt = np.linalg.svd(points - centre, full_matrices=False)
    normal = vt[2]
    scale = max(1.0, float(np.max(np.abs(points))))
    if float(np.max(np.abs((points - centre) @ normal))) > PLANE_TOLERANCE * scale:
        return None, None
    a = vt[0]
    b = vt[1]
    extent = float(np.max(np.linalg.norm(points - centre, axis=1)))
    # A parallelogram patch three times the face's size around it (parameters are clamped).
    origin = centre - 3.0 * extent * (a + b)
    surface = _core.PlaneSurface(
        tuple(origin.tolist()),
        tuple((6.0 * extent * a).tolist()),
        tuple((6.0 * extent * b).tolist()),
    )
    if normal[np.argmax(np.abs(normal))] < 0:
        normal = -normal
    key = ("plane", *np.round(normal, 9).tolist(), round(float(normal @ centre), 9))
    return surface, key


@dataclass
class Layout3D:
    ids: dict[int, np.ndarray]  # block id -> (nk, nj, ni) global ids
    domain_ids: dict[int, np.ndarray]  # face domain id -> (nj, ni) ids in the domain's layout
    points: np.ndarray
    hexes: np.ndarray
    fixed: np.ndarray
    surfaces: list
    chains: list
    slides: list[tuple[int, str, int]]
    walls: list[tuple[int, float, float, float, float]]
    roles: dict[int, str]  # face domain id -> role
    interfaces: list[tuple[int, str, int, str]]  # (block, face, block, face) per interface
    symmetric: np.ndarray = None  # (n,) uint8: nodes of symmetry-plane faces


@dataclass
class GroupResult3D:
    grids: dict[int, np.ndarray]  # block id -> (nk, nj, ni, 3)
    faces: dict[int, np.ndarray]  # face domain id -> (nj, ni, 3) in the domain's layout
    residuals: np.ndarray
    iterations: int
    converged: bool
    cancelled: bool
    node_counts: dict[str, int]
    lines: int
    colours: int
    roles: dict[int, str]
    max_kink: float
    warnings: list[str] = field(default_factory=list)


def face_roles(store, params: dict, faces: dict[int, dict]) -> dict[int, str]:
    """Role of every face domain of the member blocks (``faces``: block -> {face: domain})."""
    members = set(params["blocks"])
    uses: dict[int, int] = {}
    for block in params["blocks"]:
        for d in faces[block].values():
            uses[d] = uses.get(d, 0) + 1
    walls = {w["domain"] for w in params["walls"]}
    slide = params["slide"]
    pinned = set(params.get("pinned", ()))
    roles = {}
    for d in sorted(uses):
        outside = [
            e for e in store.dependents(d) if e.kind in ("block", "volume") and e.id not in members
        ]
        if outside:
            roles[d] = "outside"
        elif d in walls:
            roles[d] = "wall"
        elif d in pinned:
            roles[d] = "fixed"
        elif uses[d] >= 2:
            roles[d] = "interface"
        elif slide == "auto" or d in slide:
            roles[d] = "sliding"
        else:
            roles[d] = "fixed"
    return roles


def build_layout_3d(store, connectors, domains, blocks, params: dict) -> Layout3D:
    """Merge the member blocks (their unsmoothed grids) into one hexahedral grid."""
    members = params["blocks"]
    grids = {}
    faces = {}
    for b in members:
        result = blocks.evaluate(b, raw=True)
        grids[b] = np.asarray(result["grid"])
        faces[b] = result["faces"]
    roles = face_roles(store, params, faces)
    all_points = _clean(np.concatenate([grids[b].reshape(-1, 3) for b in members]))
    unique, first, inverse = np.unique(all_points, axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    order = np.argsort(first, kind="stable")
    rank = np.empty(len(unique), dtype=np.int64)
    rank[order] = np.arange(len(unique))
    flat = rank[inverse]
    points = all_points[first[order]]
    ids = {}
    offset = 0
    hexes = []
    for b in members:
        nk, nj, ni = grids[b].shape[:3]
        block_ids = flat[offset : offset + nk * nj * ni].reshape(nk, nj, ni)
        offset += nk * nj * ni
        ids[b] = block_ids
        corners = []
        for c in range(8):
            di = c & 1
            dj = (c >> 1) & 1
            dk = (c >> 2) & 1
            corners.append(block_ids[dk : nk - 1 + dk, dj : nj - 1 + dj, di : ni - 1 + di])
        cells = np.stack(corners, axis=-1).reshape(-1, 8)
        g = grids[b]
        o = g[0, 0, 0]
        volume = np.dot(np.cross(g[0, 0, 1] - o, g[0, 1, 0] - o), g[1, 0, 0] - o)
        if volume < 0:
            cells = cells[:, [1, 0, 3, 2, 5, 4, 7, 6]]  # mirrored in i
        hexes.append(cells)
    # Coordinates -> id for boundary nodes (face domains and their connectors).
    lookup: dict[bytes, int] = {}
    for b in members:
        for face in FACE_NAMES:
            fid = face_ids(ids[b], face).reshape(-1)
            fpts = points[fid]
            for k in range(len(fid)):
                lookup[fpts[k].tobytes()] = int(fid[k])

    def id_of(point) -> int:
        key = _clean(np.asarray(point, dtype=float)).tobytes()
        if key not in lookup:
            where = f"node {tuple(point)}"
            raise ValidationError([Issue("smoothing group", f"{where} is not on a block face")])
        return lookup[key]

    domain_ids = {}
    for d in roles:
        grid = _clean(domains.grid(d))
        flat_ids = np.array([id_of(p) for p in grid.reshape(-1, 3)], dtype=np.int64)
        domain_ids[d] = flat_ids.reshape(grid.shape[:2])
    n = len(points)
    fixed = np.zeros(n, dtype=np.uint8)
    surfaces = []
    surface_of_face = {}
    for d in sorted(roles):
        if roles[d] != "sliding":
            continue
        surface_id = store.get(d).params.get("surface")
        if surface_id is not None:
            surface_of_face[d] = (store.geometry(surface_id), ("surface", surface_id))
            continue
        plane, key = _plane(domains.grid(d))
        if plane is None:
            roles[d] = "fixed (curved)"
            continue
        surface_of_face[d] = (plane, key)
    # Nodes of wall faces on a surface go onto it (the curves between them may only approximate
    # the surface, e.g. the edges of a sketch), unless they also lie on another surface.
    on_surface: dict[int, set] = {}
    for d in sorted(roles):
        surface_id = store.get(d).params.get("surface")
        if roles[d] == "wall" and surface_id is not None:
            for node in domain_ids[d].reshape(-1):
                on_surface.setdefault(int(node), set()).add(surface_id)
    for node, ids_on in sorted(on_surface.items()):
        if len(ids_on) == 1:
            surface = store.geometry(next(iter(ids_on)))
            points[node] = surface.project(tuple(points[node].tolist()))[2]
    # Coplanar sliding faces share one plane patch, fitted to all of them (a patch around one
    # face only would clamp the nodes of the others).
    coplanar: dict[tuple, list[int]] = {}
    for d in sorted(surface_of_face):
        key = surface_of_face[d][1]
        if key[0] == "plane":
            coplanar.setdefault(key, []).append(d)
    for key, sharing in coplanar.items():
        if len(sharing) > 1:
            grid = np.concatenate([domains.grid(d).reshape(-1, 3) for d in sharing])
            plane, _key = _plane(grid)
            if plane is not None:
                for d in sharing:
                    surface_of_face[d] = (plane, key)
    surface_index = {}
    for d in sorted(surface_of_face):
        surface, key = surface_of_face[d]
        if key not in surface_index:
            surface_index[key] = len(surfaces)
            surfaces.append(surface)
    boundary = [d for d in sorted(roles) if roles[d] != "interface"]
    # Fixed faces: every node.
    for d in boundary:
        if roles[d] != "sliding":
            fixed[domain_ids[d].reshape(-1)] = 1
    # Sliding faces: interior nodes on the face's surface.
    slide_of: dict[int, tuple[str, int]] = {}
    for d in boundary:
        if roles[d] != "sliding":
            continue
        index = surface_index[surface_of_face[d][1]]
        for node in domain_ids[d][1:-1, 1:-1].reshape(-1):
            slide_of[int(node)] = ("surface", index)
    # Boundary edges (connectors of boundary faces).
    edge_faces: dict[int, list[int]] = {}
    for d in boundary:
        for edge in store.get(d).params["edges"]:
            for item in edge:
                edge_faces.setdefault(item["connector"], []).append(d)
    chains = []
    inside_surface: dict[int, bool] = {}  # connector end node -> every connector there is
    for c in sorted(edge_faces):
        nodes = [id_of(p) for p in connectors.nodes(c).points]
        holders = edge_faces[c]
        fixed[nodes[0]] = 1
        fixed[nodes[-1]] = 1
        inside = False
        if all(roles[d] == "sliding" for d in holders):
            inside = len({surface_of_face[d][1] for d in holders}) == 1
        for end in (nodes[0], nodes[-1]):
            inside_surface[end] = inside_surface.get(end, True) and inside
        if any(roles[d] != "sliding" for d in holders):
            fixed[nodes] = 1
            continue
        keys = {surface_of_face[d][1] for d in holders}
        if len(keys) == 1:
            index = surface_index[keys.pop()]
            for node in nodes[1:-1]:
                slide_of[node] = ("surface", index)
            continue
        if len(holders) != 2:
            fixed[nodes] = 1
            continue
        chain_index = len(chains)
        chains.append(connectors.chain(c))
        for node in nodes[1:-1]:
            slide_of[node] = ("chain", chain_index)
    if params.get("joints") == "slide":
        # Face corners inside one surface (all faces there slide on it, all connectors ending
        # there lie inside it) slide on it too.
        corner_faces: dict[int, set] = {}
        for d in boundary:
            g = domain_ids[d]
            for node in (g[0, 0], g[0, -1], g[-1, 0], g[-1, -1]):
                corner_faces.setdefault(int(node), set()).add(d)
        for node, holding in sorted(corner_faces.items()):
            if not inside_surface.get(node, False):
                continue
            if any(roles[d] != "sliding" for d in holding):
                continue
            keys = {surface_of_face[d][1] for d in holding}
            if len(keys) == 1:
                fixed[node] = 0
                slide_of[node] = ("surface", surface_index[keys.pop()])
    slides = [
        (node, kind, index) for node, (kind, index) in sorted(slide_of.items()) if not fixed[node]
    ]
    walls = []
    seen = set()
    for wall in params["walls"]:
        d = wall["domain"]
        spacing = wall["spacing"]
        if spacing is None:
            spacing = 0.0
        for node in domain_ids[d].reshape(-1):
            if int(node) in seen:
                continue
            seen.add(int(node))
            walls.append(
                (int(node), spacing, wall["decay"], wall["relaxation"], wall["max_source"])
            )
    interfaces = []
    holders: dict[int, list[tuple[int, str]]] = {}
    for b in members:
        for face, d in faces[b].items():
            holders.setdefault(d, []).append((b, face))
    for d in sorted(holders):
        if roles[d] == "interface" and len(holders[d]) == 2:
            (ba, fa), (bb, fb) = holders[d]
            interfaces.append((ba, fa, bb, fb))
    symmetric = np.zeros(n, dtype=np.uint8)
    for d in params.get("symmetry", ()):
        if d in domain_ids and roles.get(d) == "sliding":
            symmetric[domain_ids[d].reshape(-1)] = 1
    return Layout3D(
        symmetric=symmetric,
        ids=ids,
        domain_ids=domain_ids,
        points=np.ascontiguousarray(points),
        hexes=np.ascontiguousarray(np.concatenate(hexes).astype(np.int64)),
        fixed=fixed,
        surfaces=surfaces,
        chains=chains,
        slides=slides,
        walls=walls,
        roles=roles,
        interfaces=interfaces,
    )


def run_group_3d(layout: Layout3D, params: dict, progress=None):
    """The solver call alone (no model state)."""
    return _core.smooth_multiblock_3d(
        layout.points,
        layout.hexes,
        layout.fixed,
        layout.surfaces,
        layout.chains,
        layout.slides,
        layout.walls,
        omega=params["omega"],
        max_iterations=params["max_iterations"],
        tolerance=params["tolerance"],
        acceleration=params["acceleration"],
        progress=progress,
        slide_limit=params["slide_limit"],
        absolute_tolerance=params.get("absolute_tolerance", 0.0),
        symmetric=layout.symmetric,
    )


def _turn(a, b, c) -> float:
    u = b - a
    v = c - b
    norm = float(np.linalg.norm(u) * np.linalg.norm(v))
    if norm == 0.0:
        return 0.0
    return math.degrees(math.acos(float(np.clip(np.dot(u, v) / norm, -1.0, 1.0))))


def interface_kinks_3d(layout: Layout3D, points: np.ndarray, irregular: np.ndarray, margin=3):
    """Largest kink (degrees) of the grid lines crossing the interfaces: turn at the interface
    node minus the mean turn at the neighbouring nodes, over the face nodes farther than
    ``margin`` lines from a singular node (``irregular`` mask) and from the face's edges."""
    worst = 0.0
    for ba, fa, bb, fb in layout.interfaces:
        sides = []
        for b, face in ((ba, fa), (bb, fb)):
            ids = layout.ids[b]
            axis, side = FACE_AXIS[face]
            # Move the face axis to the front, inward first: (depth, a, b) arrays of ids.
            moved = np.moveaxis(ids, 2 - axis, 0)
            if side > 0:
                moved = moved[::-1]
            sides.append(moved)
        a_ids, b_ids = sides
        face_a = a_ids[0]
        # Matching node positions in b's face.
        where_b = {int(v): idx for idx, v in np.ndenumerate(b_ids[0])}
        rows, cols = face_a.shape
        near = np.zeros(face_a.shape, dtype=bool)
        # Lines next to the face's edges follow the grid of the faces there (a fixed wall
        # grid, e.g. the patches of a cubed sphere, has its own turns): left out too.
        near[: margin + 1] = True
        near[rows - margin - 1 :] = True
        near[:, : margin + 1] = True
        near[:, cols - margin - 1 :] = True
        for r, c in np.argwhere(irregular[face_a]):
            near[max(0, r - margin) : r + margin + 1, max(0, c - margin) : c + margin + 1] = True
        for r in range(1, rows - 1):
            for c in range(1, cols - 1):
                if near[r, c] or a_ids.shape[0] < 3 or b_ids.shape[0] < 3:
                    continue
                node = int(face_a[r, c])
                rb, cb = where_b[node]
                g = points[node]
                a1 = points[a_ids[1, r, c]]
                a2 = points[a_ids[2, r, c]]
                b1 = points[b_ids[1, rb, cb]]
                b2 = points[b_ids[2, rb, cb]]
                here = _turn(a1, g, b1)
                around = 0.5 * (_turn(a2, a1, g) + _turn(g, b1, b2))
                worst = max(worst, abs(here - around))
    return worst


def finish_3d(layout: Layout3D, params: dict, raw) -> GroupResult3D:
    from opengrid.model.groups import NODE_CLASSES

    points, classes, residuals, iterations, converged, cancelled, lines, colours, warnings = raw
    classes = np.asarray(classes)
    counts = np.bincount(classes.astype(np.int64), minlength=4)
    # Singular nodes: any number of hexahedra around them other than a regular interior (8),
    # boundary face (4), boundary edge (2) or boundary corner (1) - irregular edges, and
    # boundary points where three blocks meet (e.g. a box corner); grid lines must turn there.
    around = np.bincount(layout.hexes.reshape(-1), minlength=len(points))
    irregular = ~np.isin(around, (1, 2, 4, 8))
    return GroupResult3D(
        grids={b: points[layout.ids[b]] for b in params["blocks"]},
        faces={d: points[ids] for d, ids in layout.domain_ids.items()},
        residuals=np.asarray(residuals),
        iterations=iterations,
        converged=converged,
        cancelled=cancelled,
        node_counts={name: int(counts[k]) for k, name in enumerate(NODE_CLASSES)},
        lines=lines,
        colours=colours,
        roles=dict(layout.roles),
        max_kink=interface_kinks_3d(layout, points, irregular),
        warnings=list(warnings),
    )
