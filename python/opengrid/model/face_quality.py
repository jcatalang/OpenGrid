"""Quality metrics of triangle / quad cells in 3D (faces of 3D shells, domains on surfaces), with
the definitions of the planar metrics (opengrid._core.quality_cells / quality_2d): corner angles
are measured in 3D and signed against the cell's Newell normal (the normal of the oriented
polygon), so a folded cell has a negative scaled Jacobian.

Keys: scaled_jacobian, skewness, aspect_ratio, min_angle, max_angle, non_orthogonality,
area_ratio, area; arrays of shape (1, m), triangles first, then quads.
"""

import numpy as np


def _cells(triangles, quads):
    cells = [np.asarray(t, dtype=np.int64) for t in triangles]
    cells += [np.asarray(q, dtype=np.int64) for q in quads]
    return cells


def face_quality(points, triangles, quads) -> dict:
    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    quads = np.asarray(quads, dtype=np.int64).reshape(-1, 4)
    out = {}
    blocks = []
    for cells, ideal in ((triangles, 60.0), (quads, 90.0)):
        if len(cells) == 0:
            continue
        p = points[cells]  # (m, n, 3)
        n = cells.shape[1]
        nxt = np.roll(p, -1, axis=1)
        prv = np.roll(p, 1, axis=1)
        e1 = nxt - p
        e2 = prv - p
        newell = np.cross(p, nxt).sum(axis=1)  # 2 x area vector
        area = 0.5 * np.linalg.norm(newell, axis=1)
        normal = newell / np.maximum(np.linalg.norm(newell, axis=1), 1e-300)[:, None]
        cross = np.cross(e1, e2)
        sin_signed = np.einsum("mkc,mc->mk", cross, normal)
        l1 = np.linalg.norm(e1, axis=2)
        l2 = np.linalg.norm(e2, axis=2)
        with np.errstate(divide="ignore", invalid="ignore"):
            jac = sin_signed / (l1 * l2)
        jac = np.nan_to_num(jac)
        if n == 3:
            jac = jac / (np.sqrt(3.0) / 2.0)
        angle = np.degrees(np.arctan2(sin_signed, np.einsum("mkc,mkc->mk", e1, e2)))
        angle = np.where(angle < 0.0, angle + 360.0, angle)
        min_angle = angle.min(axis=1)
        max_angle = angle.max(axis=1)
        skew = np.maximum((max_angle - ideal) / (180.0 - ideal), (ideal - min_angle) / ideal)
        if n == 3:
            la, lb, lc = l1[:, 1], l1[:, 2], l1[:, 0]
            with np.errstate(divide="ignore", invalid="ignore"):
                circ = la * lb * lc / (4.0 * area)
                inr = 2.0 * area / (la + lb + lc)
                aspect = circ / (2.0 * inr)
        else:
            a0 = 0.5 * (l1[:, 0] + l1[:, 2])
            a1 = 0.5 * (l1[:, 1] + l1[:, 3])
            with np.errstate(divide="ignore", invalid="ignore"):
                aspect = np.maximum(a0, a1) / np.minimum(a0, a1)
        blocks.append(
            {
                "scaled_jacobian": jac.min(axis=1),
                "skewness": skew,
                "aspect_ratio": np.nan_to_num(aspect, nan=np.inf),
                "min_angle": min_angle,
                "max_angle": max_angle,
                "area": area,
                "centroid": p.mean(axis=1),
                "normal": normal,
            }
        )
    keys = ("scaled_jacobian", "skewness", "aspect_ratio", "min_angle", "max_angle", "area")
    centroid = np.zeros((0, 3))
    normal = np.zeros((0, 3))
    for key in keys:
        out[key] = np.zeros(0)
    if blocks:
        for key in keys:
            out[key] = np.concatenate([b[key] for b in blocks])
        centroid = np.concatenate([b["centroid"] for b in blocks])
        normal = np.concatenate([b["normal"] for b in blocks])
    # Edge neighbours: area ratio and non-orthogonality (in-surface edge normal vs centroid link).
    cells = _cells(triangles, quads)
    m = len(cells)
    owner: dict[tuple, list] = {}
    for c, cell in enumerate(cells):
        k = len(cell)
        for e in range(k):
            a, b = int(cell[e]), int(cell[(e + 1) % k])
            owner.setdefault((min(a, b), max(a, b)), []).append((c, a, b))
    area_ratio = np.ones(m)
    non_orth = np.zeros(m)
    area = out["area"]
    for (_a, _b), uses in owner.items():
        for c, a, b in uses:
            edge = points[b] - points[a]
            outward = np.cross(edge, normal[c])
            others = [u for u in uses if u[0] != c]
            if others:
                d = centroid[others[0][0]] - centroid[c]
                o = others[0][0]
                if area[c] > 0 and area[o] > 0:
                    area_ratio[c] = max(area_ratio[c], area[c] / area[o], area[o] / area[c])
            else:
                d = 0.5 * (points[a] + points[b]) - centroid[c]
            na = np.linalg.norm(outward)
            nd = np.linalg.norm(d)
            if na > 0 and nd > 0:
                cosine = min(1.0, abs(outward @ d) / (na * nd))
                non_orth[c] = max(non_orth[c], np.degrees(np.arccos(cosine)))
    out["area_ratio"] = area_ratio
    out["non_orthogonality"] = non_orth
    return {key: value.reshape(1, -1) for key, value in out.items()}
