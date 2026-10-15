"""NumPy reference implementation of the volume-mesh quality metrics (the implementation used
before Phase 9c moved them to C++, opengrid._core.volume_quality); tests compare the two."""

import math

import numpy as np

from opengrid.model.volumes import CELL_FACES, _area_vector, cell_volumes, volume_cells

# Corners (node, three neighbours) whose triple product is positive for a valid cell, and the
# factor scaling the ideal cell's corner value to 1.
CELL_CORNERS = {
    "tets": (((0, 1, 2, 3), (1, 2, 0, 3), (2, 0, 1, 3), (3, 0, 2, 1)), math.sqrt(2.0)),
    "pyramids": (((0, 1, 3, 4), (1, 2, 0, 4), (2, 3, 1, 4), (3, 0, 2, 4)), math.sqrt(2.0)),
    "prisms": (
        ((0, 1, 2, 3), (1, 2, 0, 4), (2, 0, 1, 5), (3, 5, 4, 0), (4, 3, 5, 1), (5, 4, 3, 2)),
        2.0 / math.sqrt(3.0),
    ),
    "hexes": (
        (
            (0, 1, 3, 4),
            (1, 2, 0, 5),
            (2, 3, 1, 6),
            (3, 0, 2, 7),
            (4, 7, 5, 0),
            (5, 4, 6, 1),
            (6, 5, 7, 2),
            (7, 6, 4, 3),
        ),
        1.0,
    ),
}
# Ideal dihedral angles (smallest, largest) of each type, for the equiangle skewness.
CELL_IDEAL_ANGLES = {
    "tets": (math.degrees(math.acos(1.0 / 3.0)),) * 2,
    "pyramids": (math.degrees(math.atan(math.sqrt(2.0))), 180.0 - math.degrees(math.acos(1 / 3))),
    "prisms": (60.0, 90.0),
    "hexes": (90.0, 90.0),
}


def _edge_face_pairs(faces):
    """(face a, face b) for every edge: the two faces traversing it in opposite directions."""
    owner = {}
    for f, face in enumerate(faces):
        for k in range(len(face)):
            owner[(face[k], face[(k + 1) % len(face)])] = f
    return tuple(sorted({tuple(sorted((owner[(a, b)], owner[(b, a)]))) for (a, b) in owner}))


CELL_EDGE_FACES = {key: _edge_face_pairs(faces) for key, faces in CELL_FACES.items()}


def volume_quality_reference(mesh) -> dict:
    """Per-cell metrics of a volume mesh (tetrahedra, pyramids, prisms, hexahedra in that order),
    arrays of shape (1, m):
      * scaled_jacobian: the smallest corner det(e1, e2, e3) / (|e1| |e2| |e3|), scaled to 1 for
        the ideal cell (regular tetrahedron, square pyramid with equilateral sides, right prism
        on an equilateral triangle, cube) and capped at 1; negative for inverted cells;
      * aspect_ratio: tetrahedra circumradius / (3 inradius), other cells longest / shortest edge;
      * min_angle / max_angle: dihedral angles between adjacent faces (degrees);
      * skewness: equiangle skewness of those angles against the type's ideal angles;
      * non_orthogonality: largest angle between a face's area vector and the link of the two
        cell centroids (boundary faces: face centre);
      * volume_ratio: largest volume ratio with a face neighbour; volume."""
    pts = mesh["points"]
    blocks: dict[str, list] = {
        k: []
        for k in ("scaled_jacobian", "aspect_ratio", "min_angle", "max_angle", "skewness", "volume")
    }
    all_faces = {3: [], 4: []}
    all_owner = {3: [], 4: []}
    centroids = []
    offset = 0
    for key, cells, _vtk in volume_cells(mesh):
        m = len(cells)
        if m == 0:
            continue
        p = pts[cells]
        vol = cell_volumes(pts, key, cells)
        corners_table, factor = CELL_CORNERS[key]
        corner = []
        for o, a, b, c in corners_table:
            e1 = p[:, a] - p[:, o]
            e2 = p[:, b] - p[:, o]
            e3 = p[:, c] - p[:, o]
            det = np.einsum("ij,ij->i", e1, np.cross(e2, e3))
            lengths = (
                np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) * np.linalg.norm(e3, axis=1)
            )
            # Scaled so the ideal cell's corners are 1; capped at 1 (taller pyramid / prism corners
            # than the ideal one are not better).
            corner.append(np.minimum(det / np.maximum(lengths, 1e-300) * factor, 1.0))
        corner = np.stack(corner, axis=1)
        if key == "tets":
            sj = np.min(np.abs(corner), axis=1) * np.sign(vol)
            aspect = _tet_aspect(p, vol)
        else:
            sj = corner.min(axis=1)
            edges = sorted({tuple(sorted(e)) for e in _cell_edges(key)})
            lengths = np.stack(
                [np.linalg.norm(p[:, b] - p[:, a], axis=1) for a, b in edges], axis=1
            )
            aspect = lengths.max(axis=1) / np.maximum(lengths.min(axis=1), 1e-300)
        normals = []
        for face in CELL_FACES[key]:
            normals.append(_area_vector(p[:, list(face)]))
        dihedral = []
        for fa, fb in CELL_EDGE_FACES[key]:
            na = normals[fa]
            nb = normals[fb]
            cos = np.einsum("ij,ij->i", na, nb) / np.maximum(
                np.linalg.norm(na, axis=1) * np.linalg.norm(nb, axis=1), 1e-300
            )
            dihedral.append(180.0 - np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
        dihedral = np.stack(dihedral, axis=1)
        lo = dihedral.min(axis=1)
        hi = dihedral.max(axis=1)
        ideal_lo, ideal_hi = CELL_IDEAL_ANGLES[key]
        skew = np.maximum((hi - ideal_hi) / (180.0 - ideal_hi), (ideal_lo - lo) / ideal_lo)
        blocks["scaled_jacobian"].append(sj)
        blocks["aspect_ratio"].append(aspect)
        blocks["min_angle"].append(lo)
        blocks["max_angle"].append(hi)
        blocks["skewness"].append(np.maximum(skew, 0.0))
        blocks["volume"].append(vol)
        centroids.append(p.mean(axis=1))
        for face in CELL_FACES[key]:
            all_faces[len(face)].append(cells[:, list(face)])
            all_owner[len(face)].append(offset + np.arange(m))
        offset += m
    out = {key: np.concatenate([*parts, np.zeros(0)]) for key, parts in blocks.items()}
    centroid = np.concatenate([*centroids, np.zeros((0, 3))])
    n_cells = offset
    non_orth = np.zeros(n_cells)
    ratio = np.ones(n_cells)
    volume = np.abs(out["volume"])
    for size in (3, 4):
        if not all_faces[size]:
            continue
        faces = np.concatenate(all_faces[size])
        owner = np.concatenate(all_owner[size]).astype(np.int64)
        fp = pts[faces]
        area = _area_vector(fp)
        keys = np.sort(faces, axis=1)
        order = np.lexsort(keys.T[::-1])
        sk = keys[order]
        same = np.all(sk[1:] == sk[:-1], axis=1)
        other = np.full(len(faces), -1, dtype=np.int64)
        first = order[:-1][same]
        second = order[1:][same]
        other[first] = owner[second]
        other[second] = owner[first]
        d = fp.mean(axis=1) - centroid[owner]
        paired = other >= 0
        d[paired] = centroid[other[paired]] - centroid[owner[paired]]
        cosine = np.abs(np.einsum("ij,ij->i", area, d)) / np.maximum(
            np.linalg.norm(area, axis=1) * np.linalg.norm(d, axis=1), 1e-300
        )
        np.maximum.at(non_orth, owner, np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0))))
        va = volume[owner[paired]]
        vb = volume[other[paired]]
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.maximum(va / vb, vb / va)
        np.maximum.at(ratio, owner[paired], np.nan_to_num(r, nan=1.0, posinf=1.0))
    out["volume_ratio"] = ratio
    out["non_orthogonality"] = non_orth
    return {key: value.reshape(1, -1) for key, value in out.items()}


def _cell_edges(key):
    return [
        (face[k], face[(k + 1) % len(face)]) for face in CELL_FACES[key] for k in range(len(face))
    ]


def _tet_aspect(p, vol):
    """Circumradius / (3 inradius) of tetrahedra (1 for the regular one)."""
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
    return circum / np.maximum(3.0 * inradius, 1e-300)
