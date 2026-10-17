"""Periodic pairing of 3D patches: block faces, volume faces, the ends and side walls of
extruded volumes.

Method. The nodes of the source patch's boundary faces are transformed (opengrid.model.periodic)
and matched to the nodes of the target patch's faces within MATCH_TOLERANCE times the model
size by a uniform grid hash (cells of twice the tolerance: a match lies in the point's own
cell or one of its 26 neighbours). The nodes must match one to one, and every source face must
map onto a target face with the same nodes, so the two face meshes correspond cell by cell,
not only node by node. Exporters get the correspondence in their own numbering (the assembled
mesh; a CGNS zone's vertex numbers or (i, j, k) triples) by the same matching on that
numbering; two block faces that map corner to corner give a CGNS 1-to-1 periodic connection
(``one_to_one_periodic``, the transform as in blocks.one_to_one_3d).
"""

from dataclasses import dataclass

import numpy as np

from opengrid.model.blocks import FACE_AXIS, face_array, face_range
from opengrid.model.periodic import MATCH_TOLERANCE, apply
from opengrid.validation import Issue, ValidationError

NEIGHBOURS = np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)])
PRIMES = np.array([73856093, 19349663, 83492791], dtype=np.int64)


@dataclass(frozen=True)
class FaceMatch:
    """A matched 3D periodic pair: how many nodes and faces map."""

    nodes: int
    faces: int


def tolerance_of(*point_sets) -> float:
    """The matching tolerance for these points: relative to the largest coordinate."""
    scale = 1.0
    for points in point_sets:
        if len(points):
            scale = max(scale, float(np.max(np.abs(points))))
    return MATCH_TOLERANCE * scale


def _cell_keys(cells: np.ndarray) -> np.ndarray:
    # Integer overflow wraps (numpy array arithmetic); a collision only adds a candidate that
    # the distance check rejects.
    return (cells[:, 0] * PRIMES[0]) ^ (cells[:, 1] * PRIMES[1]) ^ (cells[:, 2] * PRIMES[2])


def match_points(image, targets, tolerance: float) -> tuple[np.ndarray, np.ndarray]:
    """For every image point: the index of the nearest target point within ``tolerance`` (-1
    when there is none) and the distance to the nearest target point found in its cell and
    the neighbouring ones (inf when there is none)."""
    image = np.asarray(image, dtype=float).reshape(-1, 3)
    targets = np.asarray(targets, dtype=float).reshape(-1, 3)
    found = np.full(len(image), -1, dtype=np.int64)
    distance = np.full(len(image), np.inf)
    if len(image) == 0 or len(targets) == 0:
        return found, distance
    h = 2.0 * tolerance
    keys = _cell_keys(np.floor(targets / h).astype(np.int64))
    order = np.argsort(keys, kind="stable")
    sorted_keys = keys[order]
    cells = np.floor(image / h).astype(np.int64)
    for offset in NEIGHBOURS:
        probe = _cell_keys(cells + offset)
        left = np.searchsorted(sorted_keys, probe, side="left")
        right = np.searchsorted(sorted_keys, probe, side="right")
        count = right - left
        single = np.flatnonzero(count == 1)
        if len(single):
            candidate = order[left[single]]
            d = np.linalg.norm(targets[candidate] - image[single], axis=1)
            better = d < distance[single]
            distance[single[better]] = d[better]
            found[single[better]] = candidate[better]
        for k in np.flatnonzero(count > 1):
            candidates = order[left[k] : right[k]]
            d = np.linalg.norm(targets[candidates] - image[k], axis=1)
            j = int(np.argmin(d))
            if d[j] < distance[k]:
                distance[k] = d[j]
                found[k] = candidates[j]
    found[distance > tolerance] = -1
    return found, distance


def _face_nodes(grid: np.ndarray, face: str) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """A block face's nodes (n, 3), their (i, j, k) (n, 3) and its quads as rows of node
    indices, in the order of blocks.face_array."""
    lo = face_range(grid, face)
    kk, jj, ii = np.meshgrid(
        np.arange(lo[2], lo[5] + 1),
        np.arange(lo[1], lo[4] + 1),
        np.arange(lo[0], lo[3] + 1),
        indexing="ij",
    )
    ijk = np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1)
    a = face_array(grid, face)
    n2, n1 = a.shape[:2]
    index = np.arange(n2 * n1).reshape(n2, n1)
    quads = np.stack([index[:-1, :-1], index[:-1, 1:], index[1:, 1:], index[1:, :-1]], axis=-1)
    return a.reshape(-1, 3), ijk, list(quads.reshape(-1, 4))


def entry_faces(session, entry) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
    """The nodes (n, 3) and faces (rows of node indices) of a boundary entry of
    volumes.analyze_3d, and where its nodes are in their zone: (i, j, k) triples (n, 3) of a
    block face, vertex numbers (n,) of a volume mesh; all 0-based."""
    store = session.store
    if entry["block"] is not None:
        block_id, face = entry["block"]
        points, ijk, quads = _face_nodes(session.block_grid(store.get(block_id).name), face)
        return points, quads, ijk
    volume_id, key = entry["volume"]
    mesh = session.volume_mesh(store.get(volume_id).name)
    faces = [np.asarray(f, dtype=np.int64) for f in mesh["faces"][key]]
    if not faces:
        return np.zeros((0, 3)), [], np.zeros(0, dtype=np.int64)
    used = np.unique(np.concatenate(faces))
    return mesh["points"][used], [np.searchsorted(used, f) for f in faces], used


def patch_faces(session, entries, patch_name: str) -> tuple[np.ndarray, list[np.ndarray]]:
    """All nodes and faces of the entries in patch ``patch_name`` (face rows index the pooled
    nodes; a node shared by two entries appears twice and is merged by the matching)."""
    pool = []
    faces = []
    offset = 0
    for entry in entries:
        if entry["patch"] != patch_name:
            continue
        points, rows, _local = entry_faces(session, entry)
        pool.append(points)
        faces.extend(row + offset for row in rows)
        offset += len(points)
    if not pool:
        return np.zeros((0, 3)), []
    return np.concatenate(pool), faces


def patch_centroid(session, patch) -> np.ndarray:
    """Centroid of the distinct nodes of a 3D patch's boundary faces."""
    from opengrid.model.volumes import analyze_3d

    points, faces = patch_faces(session, analyze_3d(session)["entries"], patch.name)
    if not faces:
        raise ValueError(f"patch {patch.name} has no boundary faces: give the translation")
    return np.unique(points, axis=0).mean(axis=0)


def match_faces(session, source, target, entries) -> tuple[FaceMatch | None, list[Issue]]:
    """Node-by-node and face-by-face correspondence of a 3D periodic pair (``entries`` from
    analyze_3d with their patches assigned); the problems found instead of a match."""
    periodic = source.params["periodic"]
    s_points, s_faces = patch_faces(session, entries, source.name)
    t_points, t_faces = patch_faces(session, entries, target.name)
    issues = []
    for patch, faces in ((source, s_faces), (target, t_faces)):
        if not faces:
            issues.append(Issue(patch.name, "has no boundary faces to pair"))
    if issues:
        return None, issues
    s_unique, s_inverse = np.unique(s_points, axis=0, return_inverse=True)
    t_unique, t_inverse = np.unique(t_points, axis=0, return_inverse=True)
    s_inverse = s_inverse.reshape(-1)
    t_inverse = t_inverse.reshape(-1)
    if len(s_unique) != len(t_unique):
        issues.append(
            Issue(source.name, f"has {len(s_unique)} nodes, {target.name} has {len(t_unique)}")
        )
    if len(s_faces) != len(t_faces):
        issues.append(
            Issue(source.name, f"has {len(s_faces)} faces, {target.name} has {len(t_faces)}")
        )
    image = apply(s_unique, periodic)
    found, _distance = match_points(image, t_unique, tolerance_of(s_unique, t_unique))
    missing = np.flatnonzero(found < 0)
    if len(missing):
        first = s_unique[missing[0]]
        nearest = float(np.min(np.linalg.norm(t_unique - image[missing[0]], axis=1)))
        issues.append(
            Issue(
                source.name,
                f"{len(missing)} of its {len(s_unique)} nodes have no periodic image in "
                f"{target.name} (its node at ({first[0]:.6g}, {first[1]:.6g}, {first[2]:.6g}) "
                f"maps {nearest:.3g} from the nearest node of {target.name})",
            )
        )
    elif len(np.unique(found)) < len(found):
        issues.append(
            Issue(source.name, f"several of its nodes map onto the same node of {target.name}")
        )
    if issues:
        return None, issues
    target_sets = {frozenset(t_inverse[row].tolist()) for row in t_faces}
    unmatched = 0
    for row in s_faces:
        if frozenset(found[s_inverse[row]].tolist()) not in target_sets:
            unmatched += 1
    if unmatched:
        issues.append(
            Issue(
                source.name,
                f"{unmatched} of its {len(s_faces)} faces map onto no face of {target.name} "
                "with the same nodes (the face meshes differ)",
            )
        )
        return None, issues
    return FaceMatch(len(s_unique), len(s_faces)), issues


def _patch_points(boundary_cells, boundary_patch: np.ndarray, index: int) -> np.ndarray:
    rows = [np.asarray(boundary_cells[k]) for k in np.flatnonzero(boundary_patch == index)]
    if not rows:
        return np.zeros(0, dtype=np.int64)
    return np.unique(np.concatenate(rows))


def assembly_links(points, boundary_cells, boundary_patch, patch_names, periodic) -> list:
    """Periodic links of an assembled 3D mesh (io/assembly3d.py) in its point numbering:
    (source patch index, target patch index, transform, target point ids, matching source
    point ids), as the 2D assembly lists them."""
    links = []
    boundary_patch = np.asarray(boundary_patch)
    for source, target, _match in periodic:
        s = patch_names.index(source.name)
        t = patch_names.index(target.name)
        s_ids = _patch_points(boundary_cells, boundary_patch, s)
        t_ids = _patch_points(boundary_cells, boundary_patch, t)
        transform = dict(source.params["periodic"])
        tolerance = tolerance_of(points[s_ids], points[t_ids])
        found, _distance = match_points(apply(points[s_ids], transform), points[t_ids], tolerance)
        if len(s_ids) != len(t_ids) or np.any(found < 0):
            raise ValidationError(
                [
                    Issue(
                        source.name,
                        f"its nodes do not map onto those of {target.name} in the assembled mesh",
                    )
                ]
            )
        links.append((s, t, transform, t_ids[found], s_ids))
    return links


def one_to_one_periodic(grid, face, donor_grid, donor_face, periodic, tolerance) -> dict | None:
    """CGNS 1-to-1 data (0-based ranges, as blocks.one_to_one_3d) of block face ``face`` of
    ``grid`` whose periodic image is the whole block face ``donor_face`` of ``donor_grid``, or
    None when the two faces do not map corner to corner (their nodes then go into point
    lists)."""
    lo = face_range(grid, face)
    donor_lo = face_range(donor_grid, donor_face)
    donor_points, donor_ijk, _quads = _face_nodes(donor_grid, donor_face)
    probes = [grid[lo[2], lo[1], lo[0]], grid[lo[5], lo[4], lo[3]]]
    axis, side = FACE_AXIS[face]
    running = [a for a in range(3) if a != axis]
    for a in running:
        step = list(lo[:3])
        step[a] += 1
        probes.append(grid[step[2], step[1], step[0]])
    found, _distance = match_points(apply(np.array(probes), periodic), donor_points, tolerance)
    if np.any(found < 0):
        return None
    begin = donor_ijk[found[0]]
    end = donor_ijk[found[1]]
    for m in range(3):
        if {int(begin[m]), int(end[m])} != {donor_lo[m], donor_lo[m + 3]}:
            return None  # the corners are not opposite corners of the donor face
    d_axis, d_side = FACE_AXIS[donor_face]
    transform = [0, 0, 0]
    transform[axis] = -side * d_side * (d_axis + 1)
    for n, a in enumerate(running):
        delta = donor_ijk[found[2 + n]] - begin
        moved = [m for m in range(3) if delta[m] != 0]
        if len(moved) != 1:
            return None
        transform[a] = int(np.sign(delta[moved[0]])) * (moved[0] + 1)
    return {
        "range": lo,
        "donor_range": (*begin.tolist(), *end.tolist()),
        "transform": tuple(transform),
    }
