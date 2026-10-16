"""Structured 3D blocks and the 3D topology.

A block is either
  * ``faces``: six structured domains enclosing a hexahedral region. Orientation and point
    matching are automatic: the eight corners are identified (exact coordinates: faces share
    connector nodes), a right-handed (i, j, k) frame is chosen at the first face's first corner,
    and every face grid is transposed / flipped onto its block face; the interior is 3D TFI
    (opengrid._core.tfi_block). Or
  * ``extrusion``: a structured domain extruded (translate / rotate / normal / hyperbolic) with a
    connector style distribution of ``dimension`` layers (opengrid._core.extrude,
    opengrid._core.extrude_hyperbolic). The domain is the block's kmin face (mirrored in i when
    needed for a right-handed block). Hyperbolic extrusion marches ``distance`` along smoothed
    layer normals with edge conditions per base edge (imin, imax, jmin, jmax: free, constant_x /
    _y / _z, periodic) and tangential ``dissipation``.
Either kind can be smoothed afterwards (``smoothing``: 3D elliptic, Laplace or Thomas-Middlecoff,
with wall control on chosen faces; opengrid._core.smooth_elliptic_3d); faces stay fixed.
Grids are arrays (nk, nj, ni, 3). Block faces are named imin, imax, jmin, jmax, kmin, kmax; their
arrays are (nk, nj, 3), (nk, ni, 3) and (nj, ni, 3) respectively.

3D topology: block faces with identical nodes are interfaces (1-to-1, with the CGNS transform),
the others are boundary faces. Patches name boundary faces through domains ({"domain": id}, the
faces of six-face blocks and the base of extruded blocks) or directly ({"block": id, "face"}).
"""

import hashlib
import math
from collections import OrderedDict
from dataclasses import dataclass, field

import numpy as np

from opengrid import _core
from opengrid.model.domains import MAX_ACCELERATION
from opengrid.validation import Issue, ValidationError

FACE_NAMES = ("imin", "imax", "jmin", "jmax", "kmin", "kmax")
EXTRUSION_MODES = ("translate", "rotate", "normal", "hyperbolic")
EDGE_CONDITIONS = ("free", "constant_x", "constant_y", "constant_z", "periodic")
SMOOTHING_CONTROLS = ("laplace", "thomas_middlecoff")
DEFAULT_WALL = {
    "face": "kmin",
    "spacing": None,
    "decay": 0.1,
    "relaxation": 0.1,
    "max_source": 100.0,
}
DEFAULT_SMOOTHING = {
    "control": "laplace",
    "walls": (),
    "omega": 1.0,
    "max_iterations": 500,
    "tolerance": 1e-6,
    "acceleration": 5,
}
DEFAULT_DISTRIBUTION = {
    "type": "uniform",
    "begin_spacing": None,
    "end_spacing": None,
    "growth_ratio": None,
}
DEFAULT_EXTRUSION = {
    "domain": None,
    "mode": "translate",
    "vector": (0.0, 0.0, 1.0),
    "axis_point": (0.0, 0.0, 0.0),
    "angle": 0.0,
    "distance": 0.0,
    "dimension": 11,
    "distribution": DEFAULT_DISTRIBUTION,
    "edges": ("free", "free", "free", "free"),
    "dissipation": 0.2,
}
# Block faces: (corner at the face origin, along the first face index, along the second), in
# terms of block corners (a, b, c) = (i, j, k) in {0, 1}.
FACE_CORNERS = {
    "imin": ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
    "imax": ((1, 0, 0), (1, 1, 0), (1, 0, 1)),
    "jmin": ((0, 0, 0), (1, 0, 0), (0, 0, 1)),
    "jmax": ((0, 1, 0), (1, 1, 0), (0, 1, 1)),
    "kmin": ((0, 0, 0), (1, 0, 0), (0, 1, 0)),
    "kmax": ((0, 0, 1), (1, 0, 1), (0, 1, 1)),
}
# Block index axis held constant on each face and its side (-1 min, +1 max).
FACE_AXIS = {
    "imin": (0, -1),
    "imax": (0, 1),
    "jmin": (1, -1),
    "jmax": (1, 1),
    "kmin": (2, -1),
    "kmax": (2, 1),
}


# ----------------------------------------------------------------------------- parameters
def _vec3(value, label, issues):
    if isinstance(value, (str, bytes)) or not hasattr(value, "__len__") or len(value) not in (2, 3):
        issues.append(Issue(label, f"expected 2 or 3 coordinates, got {value!r}"))
        return None
    out = []
    for v in value:
        if hasattr(v, "item") and not isinstance(v, (int, float)):
            v = v.item()
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
            issues.append(Issue(label, f"expected numbers, got {value!r}"))
            return None
        out.append(float(v))
    if len(out) == 2:
        out.append(0.0)
    return tuple(out)


def _number(value, label, issues, optional=False):
    if value is None and optional:
        return None
    if hasattr(value, "item") and not isinstance(value, (int, float)):
        value = value.item()
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        issues.append(Issue(label, f"expected a finite number, got {value!r}"))
        return None
    return float(value)


def _domain_ref(ref, label, store, issues):
    if isinstance(ref, bool) or not isinstance(ref, (str, int)) or ref not in store:
        issues.append(Issue(label, f"no domain {ref!r}"))
        return None
    entity = store.get(ref)
    if entity.kind != "domain":
        issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}, not a structured domain"))
        return None
    return entity.id


def normalize_block(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    from opengrid.model.quality_floor import normalize_quality_floor

    for key in sorted(set(raw) - {"faces", "extrusion", "smoothing", "quality_floor"}):
        issues.append(Issue(f"{where}: {key}", "unexpected parameter for a block"))
    smoothing = _normalize_smoothing(raw.get("smoothing"), f"{where}: smoothing", issues)
    floor = normalize_quality_floor(raw.get("quality_floor"), f"{where}: quality_floor", issues)
    faces_raw = raw.get("faces")
    extrusion_raw = raw.get("extrusion")
    if (faces_raw is None) == (extrusion_raw is None):
        raise ValidationError(
            [*issues, Issue(where, "give either six faces or an extrusion (not both)")]
        )
    faces = None
    extrusion = None
    if faces_raw is not None:
        if not isinstance(faces_raw, (list, tuple)) or len(faces_raw) != 6:
            issues.append(Issue(f"{where}: faces", "exactly six structured domains are required"))
        else:
            faces = []
            for k, ref in enumerate(faces_raw):
                faces.append(_domain_ref(ref, f"{where}: faces[{k}]", store, issues))
            if len({f for f in faces if f is not None}) != len([f for f in faces if f is not None]):
                issues.append(Issue(f"{where}: faces", "a domain is listed twice"))
            faces = tuple(faces)
    else:
        extrusion = _normalize_extrusion(extrusion_raw, f"{where}: extrusion", store, issues)
    if issues:
        raise ValidationError(issues)
    return {"faces": faces, "extrusion": extrusion, "smoothing": smoothing, "quality_floor": floor}


def _normalize_smoothing(raw, label, issues) -> dict | None:
    if raw is None:
        return None
    if not isinstance(raw, dict):
        issues.append(Issue(label, f"expected a dict or None, got {raw!r}"))
        return None
    for key in sorted(set(raw) - set(DEFAULT_SMOOTHING)):
        issues.append(Issue(f"{label}.{key}", "unexpected key"))
    merged = dict(DEFAULT_SMOOTHING)
    merged.update(raw)
    if merged["control"] not in SMOOTHING_CONTROLS:
        issues.append(Issue(f"{label}.control", f"must be one of {', '.join(SMOOTHING_CONTROLS)}"))
    walls = []
    raw_walls = merged["walls"]
    if isinstance(raw_walls, dict):
        raw_walls = [raw_walls]
    for k, w in enumerate(raw_walls):
        wl = f"{label}.walls[{k}]"
        if isinstance(w, str):
            w = {"face": w}
        if not isinstance(w, dict):
            issues.append(Issue(wl, "expected a dict or a face name"))
            continue
        for key in sorted(set(w) - set(DEFAULT_WALL)):
            issues.append(Issue(f"{wl}.{key}", "unexpected key"))
        wall = dict(DEFAULT_WALL)
        wall.update(w)
        if wall["face"] not in FACE_NAMES:
            issues.append(Issue(f"{wl}.face", f"must be one of {', '.join(FACE_NAMES)}"))
        out = {"face": wall["face"]}
        for key in ("spacing", "decay", "relaxation", "max_source"):
            out[key] = _number(wall[key], f"{wl}.{key}", issues, optional=key == "spacing")
        walls.append(out)
    iterations = merged["max_iterations"]
    if isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 1:
        issues.append(Issue(f"{label}.max_iterations", "must be an integer >= 1"))
    acceleration = merged["acceleration"]
    if (
        isinstance(acceleration, bool)
        or not isinstance(acceleration, int)
        or not 0 <= acceleration <= MAX_ACCELERATION
    ):
        issues.append(
            Issue(
                f"{label}.acceleration",
                f"must be an integer from 0 (off) to {MAX_ACCELERATION}, got {acceleration!r}",
            )
        )
    return {
        "control": merged["control"],
        "walls": tuple(walls),
        "omega": _number(merged["omega"], f"{label}.omega", issues),
        "max_iterations": iterations,
        "tolerance": _number(merged["tolerance"], f"{label}.tolerance", issues),
        "acceleration": acceleration,
    }


def _normalize_extrusion(raw, label, store, issues) -> dict | None:
    if not isinstance(raw, dict):
        issues.append(Issue(label, f"expected a dict, got {raw!r}"))
        return None
    for key in sorted(set(raw) - set(DEFAULT_EXTRUSION)):
        issues.append(Issue(f"{label}.{key}", "unexpected key"))
    merged = dict(DEFAULT_EXTRUSION)
    merged.update(raw)
    domain = _domain_ref(merged["domain"], f"{label}.domain", store, issues)
    mode = merged["mode"]
    if mode not in EXTRUSION_MODES:
        issues.append(
            Issue(f"{label}.mode", f"must be one of {', '.join(EXTRUSION_MODES)}, got {mode!r}")
        )
    dimension = merged["dimension"]
    if isinstance(dimension, bool) or not isinstance(dimension, int) or dimension < 2:
        issues.append(Issue(f"{label}.dimension", f"must be an integer >= 2, got {dimension!r}"))
    distribution = dict(DEFAULT_DISTRIBUTION)
    if not isinstance(merged["distribution"], dict):
        issues.append(Issue(f"{label}.distribution", "expected a dict"))
    else:
        for key in sorted(set(merged["distribution"]) - set(DEFAULT_DISTRIBUTION)):
            issues.append(Issue(f"{label}.distribution.{key}", "unexpected key"))
        distribution.update(merged["distribution"])
        for key in ("begin_spacing", "end_spacing", "growth_ratio"):
            distribution[key] = _number(
                distribution[key], f"{label}.distribution.{key}", issues, optional=True
            )
    out = {
        "domain": domain,
        "mode": mode,
        "vector": _vec3(merged["vector"], f"{label}.vector", issues),
        "axis_point": _vec3(merged["axis_point"], f"{label}.axis_point", issues),
        "angle": _number(merged["angle"], f"{label}.angle", issues),
        "distance": _number(merged["distance"], f"{label}.distance", issues),
        "dimension": dimension,
        "distribution": distribution,
        "edges": _edges(merged["edges"], label, issues),
        "dissipation": _number(merged["dissipation"], f"{label}.dissipation", issues),
    }
    if mode == "translate" and out["vector"] is not None and not any(out["vector"]):
        issues.append(Issue(f"{label}.vector", "must not be zero for a translation"))
    if mode == "rotate" and out["angle"] == 0.0:
        issues.append(Issue(f"{label}.angle", "must not be zero for a rotation (radians)"))
    if mode in ("normal", "hyperbolic") and out["distance"] == 0.0:
        issues.append(Issue(f"{label}.distance", "must not be zero for a normal extrusion"))
    return out


def _edges(raw, label, issues) -> tuple:
    if not isinstance(raw, (list, tuple)) or len(raw) != 4:
        issues.append(Issue(f"{label}.edges", "four conditions (imin, imax, jmin, jmax)"))
        return ("free", "free", "free", "free")
    for value in raw:
        if value not in EDGE_CONDITIONS:
            issues.append(
                Issue(
                    f"{label}.edges", f"conditions are {', '.join(EDGE_CONDITIONS)}, got {value!r}"
                )
            )
    return tuple(raw)


def block_references(params: dict) -> set[int]:
    if params["faces"] is not None:
        return set(params["faces"])
    return {params["extrusion"]["domain"]}


def block_script_params(params: dict, store) -> dict:
    result = dict(params)
    if params["faces"] is not None:
        result["faces"] = tuple(store.get(f).name for f in params["faces"])
    else:
        extrusion = dict(params["extrusion"])
        extrusion["domain"] = store.get(extrusion["domain"]).name
        result["extrusion"] = extrusion
    return result


# ------------------------------------------------------------------- six-face orientation
def _key(point) -> tuple:
    return tuple(float(c) for c in point)


def _transforms(grid: np.ndarray):
    """The eight index transforms of a face grid (n2, n1, 3)."""
    for transpose in (False, True):
        g = grid
        if transpose:
            g = np.swapaxes(grid, 0, 1)
        for flip_rows in (False, True):
            for flip_cols in (False, True):
                h = g
                if flip_rows:
                    h = h[::-1]
                if flip_cols:
                    h = h[:, ::-1]
                yield np.ascontiguousarray(h)


def orient_faces(grids: list[np.ndarray], names: list[str], where: str) -> dict:
    """Map six face grids onto the faces of a right-handed block.

    Returns {block face: (index into grids, oriented array)}."""
    corners = []
    for g in grids:
        corners.append((_key(g[0, 0]), _key(g[0, -1]), _key(g[-1, 0]), _key(g[-1, -1])))
    vertices = sorted({c for face in corners for c in face})
    if len(vertices) != 8:
        raise ValidationError(
            [Issue(where, f"the six faces have {len(vertices)} distinct corners, not 8 (no box)")]
        )
    neighbours: dict[tuple, set] = {v: set() for v in vertices}
    for c00, c01, c10, c11 in corners:
        for a, b in ((c00, c01), (c01, c11), (c11, c10), (c10, c00)):
            neighbours[a].add(b)
            neighbours[b].add(a)
    for v, near in neighbours.items():
        if len(near) != 3:
            raise ValidationError(
                [Issue(where, f"corner {v} joins {len(near)} edges, not 3 (faces do not close)")]
            )
    origin, along_i, along_j = corners[0][0], corners[0][1], corners[0][2]
    (along_k,) = neighbours[origin] - {along_i, along_j}
    o = np.array(origin)
    if np.linalg.det(np.array([along_i, along_j, along_k]) - o) < 0:
        along_i, along_j = along_j, along_i
    vertex = {(0, 0, 0): origin, (1, 0, 0): along_i, (0, 1, 0): along_j, (0, 0, 1): along_k}

    def fourth(a, b, c):
        for face in corners:
            if {a, b, c} <= set(face):
                (d,) = set(face) - {a, b, c}
                return d
        raise ValidationError([Issue(where, "faces do not form a box")])

    vertex[(1, 1, 0)] = fourth(origin, along_i, along_j)
    vertex[(1, 0, 1)] = fourth(origin, along_i, along_k)
    vertex[(0, 1, 1)] = fourth(origin, along_j, along_k)
    (last,) = set(vertices) - set(vertex.values())
    vertex[(1, 1, 1)] = last
    out = {}
    used = set()
    for face, (c0, c1, c2) in FACE_CORNERS.items():
        o, a, b = vertex[c0], vertex[c1], vertex[c2]
        wanted = {o, a, b}
        match = None
        for index, g in enumerate(grids):
            if index in used or not wanted <= set(corners[index]):
                continue
            for t in _transforms(g):
                if _key(t[0, 0]) == o and _key(t[0, -1]) == a and _key(t[-1, 0]) == b:
                    match = (index, t)
                    break
            if match is not None:
                break
        if match is None:
            raise ValidationError([Issue(where, f"no face fits the block's {face} side")])
        used.add(match[0])
        out[face] = match
    return out


# ------------------------------------------------------------------------------ evaluator
class BlockEvaluator:
    """Grids of blocks, cached by a key over parameters and face grids."""

    CACHE_SIZE = 16

    def __init__(self, store, domains):
        self.store = store
        self.domains = domains
        self.groups = None  # GroupEvaluator: members of 3D smoothing groups take their grids there
        self._cache: OrderedDict[str, dict] = OrderedDict()

    def _inputs(self, entity_id: int, raw: bool = False) -> tuple[str, list]:
        params = self.store.get(entity_id).params
        digest = hashlib.sha256(repr(params).encode())
        digest.update(repr(raw).encode())
        grids = []
        refs = params["faces"]
        if refs is None:
            refs = (params["extrusion"]["domain"],)
        for ref in refs:
            g = np.ascontiguousarray(self.domains.grid(ref))
            grids.append(g)
            digest.update(g.tobytes())
        return digest.hexdigest(), grids

    def evaluate(self, entity_id: int, raw: bool = False, compute: bool = True) -> dict:
        """{"grid", "faces": {face: domain id or None}, "mirrored", "warnings", "smoothing",
        "group"}. A member of a 3D smoothing group gets its part of the group's result (run now
        unless ``compute`` is False: then the unsmoothed grid); ``raw``: the block's own grid
        without its smoothing (the group's input)."""
        if not raw and self.groups is not None:
            group = self.groups.group_of_block(entity_id)
            if group is not None:
                base = self.evaluate(entity_id, raw=True)
                cached = self.groups.result(group)
                if cached is None and compute:
                    cached = self.groups.smooth(group)
                if cached is None:
                    return dict(base, group=group)
                smoothing = {
                    "iterations": cached.iterations,
                    "converged": cached.converged,
                    "residuals": cached.residuals,
                }
                return dict(
                    base,
                    grid=cached.grids[entity_id],
                    smoothing=smoothing,
                    warnings=list(base["warnings"]) + list(cached.warnings),
                    group=group,
                )
        key, grids = self._inputs(entity_id, raw)
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        entity = self.store.get(entity_id)
        params = entity.params
        warnings: list[str] = []
        if params["faces"] is not None:
            names = [self.store.get(f).name for f in params["faces"]]
            oriented = orient_faces(grids, names, entity.name)
            try:
                grid = _core.tfi_block(*(oriented[f][1] for f in FACE_NAMES))
            except ValueError as error:
                raise ValidationError([Issue(entity.name, str(error))]) from None
            faces = {f: params["faces"][oriented[f][0]] for f in FACE_NAMES}
            mirrored = False
        else:
            e = params["extrusion"]
            d = e["distribution"]
            fractions, _warnings = _core.distribute(
                e["dimension"], d["type"], d["begin_spacing"], d["end_spacing"], d["growth_ratio"]
            )
            vector = e["vector"]
            try:
                if e["mode"] == "hyperbolic":
                    grid, mirrored, _inverted, messages = _core.extrude_hyperbolic(
                        grids[0],
                        e["distance"],
                        list(fractions),
                        edges=list(e["edges"]),
                        dissipation=e["dissipation"],
                    )
                    warnings.extend(messages)
                else:
                    grid, mirrored = _core.extrude(
                        grids[0],
                        e["mode"],
                        vector=vector,
                        axis_point=e["axis_point"],
                        angle=e["angle"],
                        distance=e["distance"],
                        fractions=list(fractions),
                    )
            except ValueError as error:
                raise ValidationError([Issue(entity.name, str(error))]) from None
            faces = {f: None for f in FACE_NAMES}
            faces["kmin"] = e["domain"]
        smoothing = None
        if params["smoothing"] is not None and not raw:
            sm = params["smoothing"]
            walls = [
                (
                    FACE_NAMES.index(w["face"]),
                    w["spacing"],
                    w["decay"],
                    w["relaxation"],
                    w["max_source"],
                )
                for w in sm["walls"]
            ]
            try:
                grid, residuals, iterations, converged, _cancelled, messages = (
                    _core.smooth_elliptic_3d(
                        np.ascontiguousarray(grid),
                        sm["control"],
                        walls,
                        sm["omega"],
                        sm["max_iterations"],
                        sm["tolerance"],
                        acceleration=sm["acceleration"],
                    )
                )
            except ValueError as error:
                raise ValidationError([Issue(entity.name, str(error))]) from None
            warnings.extend(messages)
            smoothing = {"iterations": iterations, "converged": converged, "residuals": residuals}
        result = {
            "grid": grid,
            "faces": faces,
            "mirrored": bool(mirrored),
            "warnings": warnings,
            "smoothing": smoothing,
            "group": None,
        }
        self._cache[key] = result
        while len(self._cache) > self.CACHE_SIZE:
            self._cache.popitem(last=False)
        return result

    def grid(self, entity_id: int) -> np.ndarray:
        return self.evaluate(entity_id)["grid"]


def face_array(grid: np.ndarray, face: str) -> np.ndarray:
    if face == "imin":
        return grid[:, :, 0]
    if face == "imax":
        return grid[:, :, -1]
    if face == "jmin":
        return grid[:, 0, :]
    if face == "jmax":
        return grid[:, -1, :]
    if face == "kmin":
        return grid[0]
    return grid[-1]


def face_range(grid: np.ndarray, face: str) -> tuple:
    """0-based (imin, jmin, kmin, imax, jmax, kmax) of a block face."""
    nk, nj, ni, _ = grid.shape
    lo = [0, 0, 0]
    hi = [ni - 1, nj - 1, nk - 1]
    axis, side = FACE_AXIS[face]
    if side < 0:
        hi[axis] = 0
    else:
        lo[axis] = hi[axis]
    return (*lo, *hi)


# ------------------------------------------------------------------------------ topology
@dataclass
class BlockFace:
    block: int  # entity id
    face: str
    domain: int | None  # domain the face comes from (six-face blocks, extrusion base)


@dataclass
class Topology3D:
    boundary: list[BlockFace] = field(default_factory=list)
    interfaces: list[tuple[BlockFace, BlockFace]] = field(default_factory=list)
    problems: list[Issue] = field(default_factory=list)


def _face_index(grid, face) -> dict:
    """Coordinates -> (i, j, k) of the nodes of a block face."""
    lo = face_range(grid, face)
    out = {}
    for k in range(lo[2], lo[5] + 1):
        for j in range(lo[1], lo[4] + 1):
            for i in range(lo[0], lo[3] + 1):
                out[_key(grid[k, j, i])] = (i, j, k)
    return out


def analyze_blocks(store, blocks: BlockEvaluator) -> Topology3D:
    topo = Topology3D()
    groups: dict[tuple, list[tuple[BlockFace, np.ndarray]]] = {}
    for entity in store.entities("block"):
        result = blocks.evaluate(entity.id)
        grid = result["grid"]
        for face in FACE_NAMES:
            a = face_array(grid, face)
            key = (
                frozenset({_key(a[0, 0]), _key(a[0, -1]), _key(a[-1, 0]), _key(a[-1, -1])}),
                tuple(sorted(a.shape[:2])),
            )
            groups.setdefault(key, []).append(
                (BlockFace(entity.id, face, result["faces"][face]), a)
            )
    for items in groups.values():
        if len(items) == 1:
            topo.boundary.append(items[0][0])
            continue
        if len(items) > 2:
            names = ", ".join(f"{store.get(f.block).name}.{f.face}" for f, _ in items)
            topo.problems.append(Issue(names, "more than two block faces coincide"))
            continue
        (fa, a), (fb, b) = items
        if any(np.array_equal(a, t) for t in _transforms(b)):
            topo.interfaces.append((fa, fb))
        else:
            topo.problems.append(
                Issue(
                    f"{store.get(fa.block).name}.{fa.face}",
                    f"shares its corners with {store.get(fb.block).name}.{fb.face} but the nodes "
                    "differ (point mismatch)",
                )
            )
            topo.boundary.extend([fa, fb])
    topo.boundary.sort(key=lambda f: (f.block, FACE_NAMES.index(f.face)))
    topo.interfaces.sort(key=lambda p: (p[0].block, FACE_NAMES.index(p[0].face)))
    return topo


def one_to_one_3d(blocks: BlockEvaluator, this: BlockFace, donor: BlockFace) -> dict:
    """CGNS 1-to-1 data (0-based ranges) from ``this`` face to ``donor``."""
    g = blocks.grid(this.block)
    h = blocks.grid(donor.block)
    lo = face_range(g, this.face)
    index = _face_index(h, donor.face)
    begin = index[_key(g[lo[2], lo[1], lo[0]])]
    end = index[_key(g[lo[5], lo[4], lo[3]])]
    axis, side = FACE_AXIS[this.face]
    d_axis, d_side = FACE_AXIS[donor.face]
    transform = [0, 0, 0]
    transform[axis] = -side * d_side * (d_axis + 1)
    for a in range(3):
        if a == axis:
            continue
        step = list(lo[:3])
        step[a] += 1
        p = index[_key(g[step[2], step[1], step[0]])]
        delta = [p[m] - begin[m] for m in range(3)]
        (m,) = [m for m in range(3) if delta[m] != 0]
        transform[a] = int(np.sign(delta[m])) * (m + 1)
    return {"range": lo, "donor_range": (*begin, *end), "transform": tuple(transform)}


def patch_faces(store, topo: Topology3D) -> tuple[dict, list[Issue]]:
    """{(block id, face): patch name} for boundary faces, and membership problems."""
    issues = []
    by_domain: dict[int, list[BlockFace]] = {}
    by_face: dict[tuple, BlockFace] = {}
    for f in topo.boundary:
        if f.domain is not None:
            by_domain.setdefault(f.domain, []).append(f)
        by_face[(f.block, f.face)] = f
    assigned: dict[tuple, list[str]] = {}
    for patch in store.entities("patch"):
        for item in patch.params["members"]:
            if "domain" in item:
                faces = by_domain.get(item["domain"], [])
                if not faces:
                    issues.append(
                        Issue(
                            store.get(item["domain"]).name,
                            f"is in patch {patch.name} but is not a boundary face of any block",
                        )
                    )
                for f in faces:
                    assigned.setdefault((f.block, f.face), []).append(patch.name)
            elif "block" in item:
                key = (item["block"], item["face"])
                if key not in by_face:
                    issues.append(
                        Issue(
                            f"{store.get(item['block']).name}.{item['face']}",
                            f"is in patch {patch.name} but is not a boundary face",
                        )
                    )
                    continue
                assigned.setdefault(key, []).append(patch.name)
    result = {}
    for key, names in assigned.items():
        if len(names) > 1:
            issues.append(
                Issue(
                    f"{store.get(key[0]).name}.{key[1]}",
                    f"belongs to several patches: {', '.join(names)}",
                )
            )
        result[key] = names[0]
    for f in topo.boundary:
        if (f.block, f.face) not in result:
            issues.append(
                Issue(f"{store.get(f.block).name}.{f.face}", "is a boundary face in no patch")
            )
    return result, issues
