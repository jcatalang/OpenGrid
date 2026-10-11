"""Unstructured 3D volumes: tetrahedra (and pyramids on quad faces) filling the region enclosed
by a set of domains.

A volume lists its boundary ``faces`` (structured or unstructured domains, planar or on
surfaces; any number of closed shells: the outer boundary and inner bodies). The face meshes are
merged into one surface by exact node coordinates (faces share connector nodes, so they merge
bitwise) and checked:
  * watertight: every surface edge is used by exactly two face cells (else the open or
    non-manifold edges are reported with the domains that use them);
  * orientable: each connected shell can be oriented consistently (cells flipped as needed);
  * the outer shell (largest enclosed volume) is oriented outward and the inner shells inward,
    so every surface normal points out of the meshed region.
Quad faces: with ``pyramids`` "auto" (default) quads of faces shared with structured blocks
become pyramids, so the interface stays conforming, and other quads are split into two triangles;
"always" / "never" force either. A pyramid's apex is placed inside
the region at 0.3 times the quad's mean edge length from its centroid along the inward normal
(side faces at ~31 degrees: pyramids on faces meeting at a right angle stay clear of each other;
much sharper concave edges are not supported; along edges where two pyramid faces meet, the
tetrahedra between the pyramids are thin - minimum dihedral angles of a few degrees - which is why
"auto" splits boundary quads), and the
quad's four side triangles replace it in the surface given to the tetrahedral backend; the
pyramid shares the quad with the face, so neighbouring structured blocks stay conforming.
Backends implement ``TetBackend.mesh(points, triangles, max_size, grading)``: ``points`` are the
surface nodes, ``triangles`` the surface oriented out of the region; they return all points
(the surface nodes first, unchanged) and the tetrahedra. Netgen (LGPL-2.1, conda-forge package
``netgen``) is the first backend, used through its Python module at run time; mmg (optimisation)
and an own 3D Delaunay mesher can be added behind the same interface. The backend is a Python
class rather than the C++ interface sketched in ARCHITECTURE: Netgen's Python module is the
supported, dynamically linked way to drive it, and keeps it an optional run-time dependency.
"""

import hashlib
import math
from collections import OrderedDict

import numpy as np

from opengrid.validation import Issue, ValidationError

BACKENDS = ("netgen",)
DEFAULT_VOLUME_GENERATOR = {
    "backend": "netgen",
    "max_size": None,
    "grading": 0.3,
    "pyramids": "auto",
}
PYRAMID_MODES = ("auto", "always", "never")


def normalize_volume(raw: dict, store, where: str) -> dict:
    issues: list[Issue] = []
    for key in sorted(set(raw) - {"faces", "generator"}):
        issues.append(Issue(f"{where}: {key}", "unexpected parameter for a volume"))
    faces = []
    refs = raw.get("faces")
    if isinstance(refs, (str, int)):
        refs = [refs]
    if not isinstance(refs, (list, tuple)) or len(refs) < 2:
        issues.append(Issue(f"{where}: faces", "at least two domains are needed"))
        refs = []
    for k, ref in enumerate(refs):
        label = f"{where}: faces[{k}]"
        if isinstance(ref, bool) or not isinstance(ref, (str, int)) or ref not in store:
            issues.append(Issue(label, f"no domain {ref!r}"))
            continue
        entity = store.get(ref)
        if entity.kind not in ("domain", "unstructured_domain"):
            issues.append(Issue(label, f"{entity.name!r} is a {entity.kind}, not a domain"))
            continue
        if entity.id in faces:
            issues.append(Issue(label, f"{entity.name!r} listed twice"))
            continue
        faces.append(entity.id)
    generator = dict(DEFAULT_VOLUME_GENERATOR)
    graw = raw.get("generator") or {}
    if not isinstance(graw, dict):
        issues.append(Issue(f"{where}: generator", "expected a dict"))
        graw = {}
    for key in sorted(set(graw) - set(DEFAULT_VOLUME_GENERATOR)):
        issues.append(Issue(f"{where}: generator.{key}", "unexpected key"))
    generator.update(graw)
    if generator["backend"] not in BACKENDS:
        issues.append(Issue(f"{where}: generator.backend", f"must be one of {', '.join(BACKENDS)}"))
    for key, optional in (("max_size", True), ("grading", False)):
        value = generator[key]
        if value is None and optional:
            continue
        if hasattr(value, "item") and not isinstance(value, (int, float)):
            value = value.item()
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
            issues.append(Issue(f"{where}: generator.{key}", f"must be positive, got {value!r}"))
            continue
        generator[key] = float(value)
    if generator["pyramids"] not in PYRAMID_MODES:
        issues.append(
            Issue(f"{where}: generator.pyramids", f"must be one of {', '.join(PYRAMID_MODES)}")
        )
    if issues:
        raise ValidationError(issues)
    return {"faces": tuple(faces), "generator": generator}


def volume_script_params(params: dict, store) -> dict:
    result = dict(params)
    result["faces"] = tuple(store.get(f).name for f in params["faces"])
    return result


# ---------------------------------------------------------------------------- the surface
def face_cells(session, domain_id: int) -> tuple[np.ndarray, list[np.ndarray]]:
    """Points and cells (triangles and quads, as index arrays) of a domain's mesh."""
    entity = session.store.get(domain_id)
    if entity.kind == "unstructured_domain":
        mesh = session.unstructured_mesh(entity.name)
        cells = [np.asarray(t) for t in mesh["triangles"]] + [np.asarray(q) for q in mesh["quads"]]
        return mesh["points"], cells
    grid = session.domain_grid(entity.name)
    nj, ni, _ = grid.shape
    index = np.arange(ni * nj).reshape(nj, ni)
    quads = np.stack(
        [index[:-1, :-1], index[:-1, 1:], index[1:, 1:], index[1:, :-1]], axis=-1
    ).reshape(-1, 4)
    return grid.reshape(-1, 3), list(quads)


def _signed_volume(points, cells) -> float:
    total = 0.0
    for cell in cells:
        p = points[cell]
        for k in range(1, len(cell) - 1):
            total += np.dot(p[0], np.cross(p[k], p[k + 1]))
    return total / 6.0


def assemble_shell(session, face_ids, where: str) -> dict:
    """Merged, oriented boundary surface: {"points", "cells", "cell_face", "shells"}."""
    store = session.store
    pool = []
    raw_cells = []
    cell_face = []
    offset = 0
    for f in face_ids:
        points, cells = face_cells(session, f)
        pool.append(points)
        for cell in cells:
            raw_cells.append(cell + offset)
            cell_face.append(f)
        offset += len(points)
    allpoints = np.concatenate(pool)
    points, inverse = np.unique(allpoints, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    cells = [inverse[c] for c in raw_cells]
    # Edges -> cells.
    edges: dict[tuple, list[int]] = {}
    for c, cell in enumerate(cells):
        n = len(cell)
        for k in range(n):
            a, b = int(cell[k]), int(cell[(k + 1) % n])
            edges.setdefault((min(a, b), max(a, b)), []).append(c)
    bad = {e: u for e, u in edges.items() if len(u) != 2}
    if bad:
        names = sorted({store.get(cell_face[c]).name for u in bad.values() for c in u})
        open_edges = sum(1 for u in bad.values() if len(u) == 1)
        non_manifold = len(bad) - open_edges
        raise ValidationError(
            [
                Issue(
                    where,
                    f"the faces do not form closed shells: {open_edges} open and {non_manifold} "
                    f"non-manifold edges (domains {', '.join(names)})",
                )
            ]
        )
    # Shells (connected components) and consistent orientation by breadth-first propagation:
    # neighbours across an edge must traverse it in opposite directions.
    shell_of = [-1] * len(cells)
    shells = []
    for start in range(len(cells)):
        if shell_of[start] >= 0:
            continue
        shell_id = len(shells)
        members = [start]
        shell_of[start] = shell_id
        queue = [start]
        while queue:
            c = queue.pop()
            cell = cells[c]
            n = len(cell)
            for k in range(n):
                a, b = int(cell[k]), int(cell[(k + 1) % n])
                for o in edges[(min(a, b), max(a, b))]:
                    if o == c:
                        continue
                    other = cells[o]
                    m = len(other)
                    same = any(
                        int(other[j]) == a and int(other[(j + 1) % m]) == b for j in range(m)
                    )
                    if shell_of[o] < 0:
                        if same:
                            cells[o] = other[::-1].copy()
                        shell_of[o] = shell_id
                        members.append(o)
                        queue.append(o)
                    elif same:
                        raise ValidationError(
                            [Issue(where, "a shell is not orientable (Moebius-like faces)")]
                        )
        shells.append(members)
    volumes = [_signed_volume(points, [cells[c] for c in m]) for m in shells]
    outer = int(np.argmax([abs(v) for v in volumes]))
    for s, members in enumerate(shells):
        want_positive = s == outer
        if (volumes[s] > 0) != want_positive:
            for c in members:
                cells[c] = cells[c][::-1].copy()
    if len(shells) > 1:
        # Inner shells must lie inside the outer one (by bounding box).
        lo = points[np.unique(np.concatenate([cells[c] for c in shells[outer]]))].min(axis=0)
        hi = points[np.unique(np.concatenate([cells[c] for c in shells[outer]]))].max(axis=0)
        for s, members in enumerate(shells):
            if s == outer:
                continue
            ids = np.unique(np.concatenate([cells[c] for c in members]))
            if np.any(points[ids].min(axis=0) < lo) or np.any(points[ids].max(axis=0) > hi):
                raise ValidationError([Issue(where, "an inner shell lies outside the outer one")])
    return {
        "points": points,
        "cells": cells,
        "cell_face": cell_face,
        "shells": len(shells),
        "volume": abs(volumes[outer]) - sum(abs(v) for s, v in enumerate(volumes) if s != outer),
    }


# ------------------------------------------------------------------------------ backends
class TetBackend:
    name = "abstract"

    def mesh(self, points, triangles, max_size, grading):  # pragma: no cover - interface
        raise NotImplementedError


class NetgenBackend(TetBackend):
    """Netgen's volume mesher from a fixed surface mesh (nglib through the Python module)."""

    name = "netgen"

    def mesh(self, points, triangles, max_size, grading):
        try:
            import netgen.meshing as ng
        except ImportError as error:
            raise ValidationError(
                [Issue("netgen", f"the Netgen backend is not installed ({error})")]
            ) from None
        m = ng.Mesh(dim=3)
        ids = [m.Add(ng.MeshPoint(ng.Pnt(*map(float, p)))) for p in points]
        fd = m.Add(ng.FaceDescriptor(bc=1, domin=1, domout=0, surfnr=1))
        for t in triangles:
            m.Add(ng.Element2D(fd, [ids[int(t[0])], ids[int(t[1])], ids[int(t[2])]]))
        options = {"grading": grading}
        if max_size is not None:
            options["maxh"] = max_size
        m.GenerateVolumeMesh(ng.MeshingParameters(**options))
        coords = np.array([tuple(p) for p in m.Points()], dtype=float)
        tets = np.array([[v.nr - 1 for v in e.vertices] for e in m.Elements3D()], dtype=np.int64)
        if len(tets) == 0:
            raise ValidationError([Issue("netgen", "the volume mesher produced no tetrahedra")])
        return coords, tets


def backend(name: str) -> TetBackend:
    if name == "netgen":
        return NetgenBackend()
    raise ValidationError([Issue("backend", f"unknown backend {name!r}")])


# ------------------------------------------------------------------------------ evaluator
def _tet_volume(p):
    return (
        np.einsum("ij,ij->i", p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))
        / 6.0
    )


class VolumeEvaluator:
    """Volume meshes, cached by a key over the generator settings and the face meshes."""

    CACHE_SIZE = 8

    def __init__(self, session):
        self.session = session
        self._cache: OrderedDict[str, dict] = OrderedDict()

    def mesh(self, entity_id: int) -> dict:
        """{"points", "tets", "pyramids", "surface_points", "faces": {domain: cells},
        "warnings", "volume"} (points: the merged surface nodes first)."""
        session = self.session
        entity = session.store.get(entity_id)
        generator = entity.params["generator"]
        shell = assemble_shell(session, entity.params["faces"], entity.name)
        digest = hashlib.sha256(repr(generator).encode())
        digest.update(shell["points"].tobytes())
        for cell in shell["cells"]:
            digest.update(np.asarray(cell, dtype=np.int64).tobytes())
        key = digest.hexdigest()
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        points = shell["points"]
        triangles = []
        pyramids = []
        apexes = []
        n_surface = len(points)
        # Faces that are also faces of structured blocks keep their quads (pyramids); with "auto"
        # other quads are split into triangles.
        block_faces = set()
        for block in session.store.entities("block"):
            block_faces.update(session.blocks.evaluate(block.id)["faces"].values())
        # Boundary cells as meshed (split quads become their two triangles), per face domain.
        faces = {}
        for cell, face in zip(shell["cells"], shell["cell_face"], strict=True):
            if len(cell) == 3:
                triangles.append(cell)
                faces.setdefault(face, []).append(cell)
                continue
            mode = generator["pyramids"]
            if mode == "never" or (mode == "auto" and face not in block_faces):
                for half in (cell[[0, 1, 2]], cell[[0, 2, 3]]):
                    triangles.append(half)
                    faces.setdefault(face, []).append(half)
                continue
            faces.setdefault(face, []).append(cell)
            p = points[cell]
            normal = np.cross(p[2] - p[0], p[3] - p[1])  # outward of the region
            normal = normal / np.linalg.norm(normal)
            size = np.mean(np.linalg.norm(p - np.roll(p, -1, axis=0), axis=1))
            # Apex height 0.3 h: the side faces rise at ~31 degrees, so pyramids on two faces
            # meeting at a right angle do not touch (at 0.5 h their sides would coincide).
            apex = p.mean(axis=0) - 0.3 * size * normal
            index = n_surface + len(apexes)
            apexes.append(apex)
            pyramids.append([cell[0], cell[3], cell[2], cell[1], index])
            for k in range(4):
                # Side triangle of the pyramid, oriented out of the remaining region (that is,
                # into the pyramid): the face edge in its own direction, then the apex.
                triangles.append([cell[k], cell[(k + 1) % 4], index])
        if apexes:
            points = np.concatenate([points, np.asarray(apexes)])
        triangles = np.asarray(triangles, dtype=np.int64)
        coords, tets = backend(generator["backend"]).mesh(
            points, triangles, generator["max_size"], generator["grading"]
        )
        if len(coords) < len(points) or not np.array_equal(coords[: len(points)], points):
            raise ValidationError(
                [Issue(entity.name, "the backend moved or reordered boundary points")]
            )
        coords[: len(points)] = points
        volumes = _tet_volume(coords[tets])
        flip = volumes < 0
        tets[flip] = tets[flip][:, [0, 2, 1, 3]]
        pyramids = np.asarray(pyramids, dtype=np.int64).reshape(-1, 5)
        warnings = []
        if np.any(np.abs(volumes) <= 1e-14 * shell["volume"] / max(len(tets), 1)):
            warnings.append("some tetrahedra are degenerate")
        result = {
            "points": coords,
            "tets": tets,
            "pyramids": pyramids,
            "n_surface": n_surface,
            "faces": faces,
            "shells": shell["shells"],
            "enclosed_volume": shell["volume"],
            "warnings": warnings,
        }
        self._cache[key] = result
        while len(self._cache) > self.CACHE_SIZE:
            self._cache.popitem(last=False)
        return result


def pyramid_volume(p) -> np.ndarray:
    """Volumes of pyramids (base 0-3, apex 4) as two tetrahedra on the base diagonal."""
    t1 = _tet_volume(p[:, [0, 1, 2, 4]])
    t2 = _tet_volume(p[:, [0, 2, 3, 4]])
    return t1 + t2


def total_volume(mesh) -> float:
    v = float(_tet_volume(mesh["points"][mesh["tets"]]).sum())
    if len(mesh["pyramids"]):
        v += float(pyramid_volume(mesh["points"][mesh["pyramids"]]).sum())
    return v


# ------------------------------------------------------------------------------- quality
def volume_quality(mesh) -> dict:
    """Per-cell metrics of a volume mesh (tetrahedra first, then pyramids), shape (1, m):
    scaled_jacobian (min corner det(e1, e2, e3) / (|e1| |e2| |e3|), scaled to 1 for the regular
    tetrahedron / square pyramid), aspect_ratio (tetrahedra: circumradius / (3 inradius);
    pyramids: longest / shortest edge), min_angle / max_angle (dihedral angles of tetrahedra,
    face corner angles of pyramids), skewness (from those angles: ideal 70.53 degrees for
    tetrahedra), non_orthogonality (face area vector vs centroid link), volume_ratio, volume."""
    pts = mesh["points"]
    tets = mesh["tets"]
    pyrs = mesh["pyramids"]
    out = {}
    p = pts[tets]
    vol = _tet_volume(p)
    corner = []
    for o, (a, b, c) in enumerate(((1, 2, 3), (2, 0, 3), (0, 1, 3), (1, 0, 2))):
        e1 = p[:, a] - p[:, o]
        e2 = p[:, b] - p[:, o]
        e3 = p[:, c] - p[:, o]
        det = np.einsum("ij,ij->i", e1, np.cross(e2, e3))
        lengths = (
            np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) * np.linalg.norm(e3, axis=1)
        )
        corner.append(det / np.maximum(lengths, 1e-300) * math.sqrt(2.0))
    tet_sj = np.min(np.abs(np.stack(corner, axis=1)), axis=1) * np.sign(vol)
    # Faces (opposite vertex k), areas and inradius / circumradius.
    faces = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    normals = [np.cross(p[:, b] - p[:, a], p[:, c] - p[:, a]) for a, b, c in faces]
    area_sum = sum(0.5 * np.linalg.norm(n, axis=1) for n in normals)
    inradius = 3.0 * np.abs(vol) / np.maximum(area_sum, 1e-300)
    a = p[:, 1] - p[:, 0]
    b = p[:, 2] - p[:, 0]
    c = p[:, 3] - p[:, 0]
    num = (
        np.einsum("ij,ij->i", c, c)[:, None] * np.cross(a, b)
        + np.einsum("ij,ij->i", a, a)[:, None] * np.cross(b, c)
        + np.einsum("ij,ij->i", b, b)[:, None] * np.cross(c, a)
    )
    circum = np.linalg.norm(num, axis=1) / np.maximum(12.0 * np.abs(vol), 1e-300)
    tet_aspect = circum / np.maximum(3.0 * inradius, 1e-300)
    # Dihedral angle along edge (i, j): angle between the opposite vertices k and m seen
    # perpendicular to the edge.
    dihedral = []
    for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        k, m = [v for v in range(4) if v not in (i, j)]
        e = p[:, j] - p[:, i]
        e = e / np.maximum(np.linalg.norm(e, axis=1), 1e-300)[:, None]
        u = p[:, k] - p[:, i]
        v = p[:, m] - p[:, i]
        u = u - np.einsum("ij,ij->i", u, e)[:, None] * e
        v = v - np.einsum("ij,ij->i", v, e)[:, None] * e
        cos = np.einsum("ij,ij->i", u, v) / np.maximum(
            np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1), 1e-300
        )
        dihedral.append(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
    dihedral = np.stack(dihedral, axis=1)
    ideal = math.degrees(math.acos(1.0 / 3.0))
    tet_min = dihedral.min(axis=1)
    tet_max = dihedral.max(axis=1)
    tet_skew = np.maximum((tet_max - ideal) / (180.0 - ideal), (ideal - tet_min) / ideal)
    blocks = {
        "scaled_jacobian": [tet_sj],
        "aspect_ratio": [tet_aspect],
        "min_angle": [tet_min],
        "max_angle": [tet_max],
        "skewness": [tet_skew],
        "volume": [vol],
    }
    if len(pyrs):
        q = pts[pyrs]
        pvol = pyramid_volume(q)
        corners = []
        for k in range(4):
            e1 = q[:, (k + 1) % 4] - q[:, k]
            e2 = q[:, (k + 3) % 4] - q[:, k]
            e3 = q[:, 4] - q[:, k]
            det = np.einsum("ij,ij->i", e1, np.cross(e2, e3))
            lengths = (
                np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) * np.linalg.norm(e3, axis=1)
            )
            corners.append(det / np.maximum(lengths, 1e-300) * math.sqrt(2.0))
        edges = [q[:, (k + 1) % 4] - q[:, k] for k in range(4)] + [
            q[:, 4] - q[:, k] for k in range(4)
        ]
        lengths = np.stack([np.linalg.norm(e, axis=1) for e in edges], axis=1)
        angles = []
        for k in range(4):
            u = q[:, (k + 1) % 4] - q[:, k]
            v = q[:, (k + 3) % 4] - q[:, k]
            angles.append(
                np.degrees(
                    np.arccos(
                        np.clip(
                            np.einsum("ij,ij->i", u, v)
                            / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)),
                            -1.0,
                            1.0,
                        )
                    )
                )
            )
        angles = np.stack(angles, axis=1)
        blocks["scaled_jacobian"].append(np.stack(corners, axis=1).min(axis=1))
        blocks["aspect_ratio"].append(lengths.max(axis=1) / np.maximum(lengths.min(axis=1), 1e-300))
        blocks["min_angle"].append(angles.min(axis=1))
        blocks["max_angle"].append(angles.max(axis=1))
        blocks["skewness"].append(
            np.maximum((angles.max(axis=1) - 90.0) / 90.0, (90.0 - angles.min(axis=1)) / 90.0)
        )
        blocks["volume"].append(pvol)
    for key, parts in blocks.items():
        out[key] = np.concatenate(parts)
    # Face neighbours (vectorised): triangle faces of tetrahedra and pyramid sides are paired by
    # their sorted vertex triples; pyramid bases (quads) lie on the boundary.
    n_tets = len(tets)
    n_cells = n_tets + len(pyrs)
    centroid = np.zeros((n_cells, 3))
    centroid[:n_tets] = pts[tets].mean(axis=1)
    tri_faces = [tets[:, [1, 2, 3]], tets[:, [0, 3, 2]], tets[:, [0, 1, 3]], tets[:, [0, 2, 1]]]
    tri_owner = [np.arange(n_tets)] * 4
    quad_faces = np.zeros((0, 4), dtype=np.int64)
    quad_owner = np.zeros(0, dtype=np.int64)
    if len(pyrs):
        centroid[n_tets:] = pts[pyrs].mean(axis=1)
        owners = np.arange(n_tets, n_cells)
        for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
            tri_faces.append(np.stack([pyrs[:, a], pyrs[:, b], pyrs[:, 4]], axis=1))
            tri_owner.append(owners)
        quad_faces = pyrs[:, :4]
        quad_owner = owners
    faces = np.concatenate(tri_faces)
    owner = np.concatenate(tri_owner)
    fp = pts[faces]
    area = 0.5 * np.cross(fp[:, 1] - fp[:, 0], fp[:, 2] - fp[:, 0])
    face_centre = fp.mean(axis=1)
    keys = np.sort(faces, axis=1)
    order = np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0]))
    sk = keys[order]
    same = np.all(sk[1:] == sk[:-1], axis=1)
    other = np.full(len(faces), -1, dtype=np.int64)
    first = order[:-1][same]
    second = order[1:][same]
    other[first] = owner[second]
    other[second] = owner[first]
    d = face_centre - centroid[owner]
    paired = other >= 0
    d[paired] = centroid[other[paired]] - centroid[owner[paired]]
    na = np.linalg.norm(area, axis=1)
    nd = np.linalg.norm(d, axis=1)
    cosine = np.abs(np.einsum("ij,ij->i", area, d)) / np.maximum(na * nd, 1e-300)
    angle = np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0)))
    non_orth = np.zeros(n_cells)
    np.maximum.at(non_orth, owner, angle)
    if len(quad_faces):
        qp = pts[quad_faces]
        qa = 0.5 * np.cross(qp[:, 2] - qp[:, 0], qp[:, 3] - qp[:, 1])
        qd = qp.mean(axis=1) - centroid[quad_owner]
        qc = np.abs(np.einsum("ij,ij->i", qa, qd)) / np.maximum(
            np.linalg.norm(qa, axis=1) * np.linalg.norm(qd, axis=1), 1e-300
        )
        np.maximum.at(non_orth, quad_owner, np.degrees(np.arccos(np.clip(qc, 0.0, 1.0))))
    volume = np.abs(out["volume"])
    ratio = np.ones(n_cells)
    va = volume[owner[paired]]
    vb = volume[other[paired]]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.maximum(va / vb, vb / va)
    r = np.nan_to_num(r, nan=1.0, posinf=1.0)
    np.maximum.at(ratio, owner[paired], r)
    out["volume_ratio"] = ratio
    out["non_orthogonality"] = non_orth
    return {key: value.reshape(1, -1) for key, value in out.items()}


# ------------------------------------------------------------------------- 3D topology
def analyze_3d(session) -> dict:
    """Boundary / interface faces of a 3D model (blocks and volumes) and their patches.

    Returns {"entries": [{"label", "block": (block id, face) or None, "volume": (volume id,
    domain id) or None, "keys": set, "patch": name or None}], "interfaces": [labels],
    "issues": [Issue]}. Keys are ("domain", id) and ("block", id, face): what patch members name.
    """
    from opengrid.model.blocks import Topology3D, analyze_blocks

    store = session.store
    issues: list[Issue] = []
    if store.entities("block"):
        topo = analyze_blocks(store, session.blocks)
    else:
        topo = Topology3D()
    issues.extend(topo.problems)
    volume_uses: dict[int, list[int]] = {}
    for volume in store.entities("volume"):
        for d in volume.params["faces"]:
            volume_uses.setdefault(d, []).append(volume.id)
    entries = []
    interfaces = []

    def name(entity_id):
        return store.get(entity_id).name

    for a, b in topo.interfaces:
        interfaces.append(f"{name(a.block)}.{a.face} | {name(b.block)}.{b.face}")
    block_domains = set()
    for f in topo.boundary:
        label = f"{name(f.block)}.{f.face}"
        if f.domain is not None:
            block_domains.add(f.domain)
        if f.domain is not None and f.domain in volume_uses:
            vols = ", ".join(name(v) for v in volume_uses[f.domain])
            interfaces.append(f"{label} | {vols} ({name(f.domain)})")
            continue
        keys = {("block", f.block, f.face)}
        if f.domain is not None:
            keys.add(("domain", f.domain))
        entries.append({"label": label, "block": (f.block, f.face), "volume": None, "keys": keys})
    for d, vols in sorted(volume_uses.items()):
        count = len(vols) + int(d in block_domains)
        if count == 1:
            entries.append(
                {
                    "label": f"{name(vols[0])}.{name(d)}",
                    "block": None,
                    "volume": (vols[0], d),
                    "keys": {("domain", d)},
                }
            )
        elif count == 2 and d not in block_domains:
            interfaces.append(f"{name(vols[0])} | {name(vols[1])} ({name(d)})")
        elif count > 2:
            issues.append(Issue(name(d), f"bounds {count} blocks / volumes (at most 2)"))
    # Patches.
    by_key: dict[tuple, list[int]] = {}
    for index, entry in enumerate(entries):
        entry["patch"] = None
        for key in entry["keys"]:
            by_key.setdefault(key, []).append(index)
    assigned: dict[int, list[str]] = {}
    for patch in store.entities("patch"):
        for item in patch.params["members"]:
            if "connector" in item:
                issues.append(
                    Issue(
                        name(item["connector"]),
                        f"is a connector in patch {patch.name}; 3D patches contain domains or "
                        "block faces",
                    )
                )
                continue
            if "block" in item:
                key = ("block", item["block"], item["face"])
                label = f"{name(item['block'])}.{item['face']}"
            else:
                key = ("domain", item["domain"])
                label = name(item["domain"])
            hits = by_key.get(key, [])
            if not hits:
                issues.append(Issue(label, f"is in patch {patch.name} but is not a boundary face"))
            for index in hits:
                assigned.setdefault(index, []).append(patch.name)
        if patch.params["periodic"] is not None:
            issues.append(Issue(patch.name, "periodic pairing of 3D faces is not supported yet"))
    for index, entry in enumerate(entries):
        names = sorted(set(assigned.get(index, [])))
        if not names:
            issues.append(Issue(entry["label"], "is a boundary face in no patch"))
        elif len(names) > 1:
            issues.append(Issue(entry["label"], f"belongs to several patches: {', '.join(names)}"))
        else:
            entry["patch"] = names[0]
    return {"entries": entries, "interfaces": interfaces, "issues": issues, "blocks": topo}
