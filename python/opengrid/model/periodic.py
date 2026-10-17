"""Periodic pairing of patches.

A periodic pair is stored on its *source* patch: ``periodic = {"partner": target patch id,
"rotation_center", "rotation_axis", "rotation_angle", "translation"}``. The transform maps
source points onto target points as a rotation about ``rotation_axis`` (a unit vector, z by
default) through ``rotation_center`` by ``rotation_angle`` (radians, right-hand rule) followed
by ``translation`` (the CGNS GridConnectivityProperty convention):
    x_target = R(axis, angle) (x_source - center) + center + translation.
2D (``match``): every connector of the source must map node by node (in either direction) onto
a distinct connector of the target. 3D: the faces of the two patches must map node by node and
face by face (opengrid.model.periodic_faces).
"""

import math
from dataclasses import dataclass

import numpy as np

from opengrid.validation import Issue

MATCH_TOLERANCE = 1e-8  # relative to the model size
DEFAULT_AXIS = (0.0, 0.0, 1.0)


@dataclass(frozen=True)
class ConnectorPair:
    source: int  # connector id
    target: int
    same_direction: bool  # source first node -> target first node


def rotation_matrix(periodic: dict) -> np.ndarray:
    """3 x 3 rotation about ``rotation_axis`` by ``rotation_angle`` (Rodrigues' formula)."""
    axis = np.asarray(periodic.get("rotation_axis", DEFAULT_AXIS), dtype=float)
    length = float(np.linalg.norm(axis))
    if length == 0.0:
        raise ValueError("the rotation axis must not be zero")
    k = axis / length
    a = float(periodic["rotation_angle"])
    cross = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    return np.eye(3) + math.sin(a) * cross + (1.0 - math.cos(a)) * (cross @ cross)


def angle_vector(periodic: dict) -> tuple[float, float, float]:
    """The rotation as a vector: axis times angle (what CGNS stores as RotationAngle)."""
    axis = np.asarray(periodic.get("rotation_axis", DEFAULT_AXIS), dtype=float)
    length = float(np.linalg.norm(axis))
    if length == 0.0:
        raise ValueError("the rotation axis must not be zero")
    vector = axis / length * float(periodic["rotation_angle"])
    return (float(vector[0]), float(vector[1]), float(vector[2]))


def cgns_property(periodic: dict) -> tuple:
    """(rotation centre, rotation angle vector, translation) as the CGNS writers take it."""
    return (
        tuple(float(v) for v in periodic["rotation_center"]),
        angle_vector(periodic),
        tuple(float(v) for v in periodic["translation"]),
    )


def apply(points: np.ndarray, periodic: dict) -> np.ndarray:
    c = np.asarray(periodic["rotation_center"], dtype=float)
    t = np.asarray(periodic["translation"], dtype=float)
    rotation = rotation_matrix(periodic)
    return (np.asarray(points, dtype=float) - c) @ rotation.T + c + t


def inverse(periodic: dict) -> dict:
    """Transform from target back to source in the same convention: the same centre and axis,
    the opposite angle, and the translation turned back (x_s = R^-1 (x_t - c) + c - R^-1 t)."""
    t = np.asarray(periodic["translation"], dtype=float)
    back = -(rotation_matrix(periodic).T @ t)
    return {
        "rotation_center": tuple(float(v) for v in periodic["rotation_center"]),
        "rotation_axis": tuple(float(v) for v in periodic.get("rotation_axis", DEFAULT_AXIS)),
        "rotation_angle": -periodic["rotation_angle"],
        "translation": tuple(float(v) for v in back),
    }


def affine_matrix(periodic: dict) -> np.ndarray:
    """4 x 4 homogeneous matrix of the transform."""
    c = np.asarray(periodic["rotation_center"], dtype=float)
    t = np.asarray(periodic["translation"], dtype=float)
    rotation = rotation_matrix(periodic)
    matrix = np.eye(4)
    matrix[:3, :3] = rotation
    matrix[:3, 3] = c + t - rotation @ c
    return matrix


def match(session, source, target) -> tuple[list[ConnectorPair], list[Issue]]:
    """Connector-by-connector correspondence of a periodic pair (entities ``source``/``target``)."""
    periodic = source.params["periodic"]
    nodes = {}
    scale = 1.0
    for patch in (source, target):
        for item in patch.params["members"]:
            if "connector" not in item:
                continue
            points = session.connectors.nodes(item["connector"]).points
            nodes[item["connector"]] = points
            scale = max(scale, float(np.max(np.abs(points))))
    tolerance = MATCH_TOLERANCE * scale
    pairs = []
    issues = []
    used = set()
    for item in source.params["members"]:
        if "connector" not in item:
            continue
        a = item["connector"]
        image = apply(nodes[a], periodic)
        best = (math.inf, None, True)  # (max node offset, target connector, same direction)
        counts = []
        found = None
        for candidate in target.params["members"]:
            if "connector" not in candidate:
                continue
            b = candidate["connector"]
            if b in used:
                continue
            counts.append(len(nodes[b]))
            if len(nodes[b]) != len(image):
                continue
            forward = float(np.max(np.linalg.norm(nodes[b] - image, axis=1)))
            backward = float(np.max(np.linalg.norm(nodes[b][::-1] - image, axis=1)))
            if forward < best[0]:
                best = (forward, b, True)
            if backward < best[0]:
                best = (backward, b, False)
            if forward <= tolerance:
                found = ConnectorPair(a, b, True)
                break
            if backward <= tolerance:
                found = ConnectorPair(a, b, False)
                break
        name = session.store.get(a).name
        if found is None:
            issues.append(
                Issue(
                    source.name,
                    f"connector {name} has no periodic image in {target.name}: "
                    + _mismatch_reason(session, image, nodes, best, counts, periodic),
                )
            )
            continue
        used.add(found.target)
        pairs.append(found)
    for item in target.params["members"]:
        if "connector" not in item:
            continue
        if item["connector"] not in used and not issues:
            name = session.store.get(item["connector"]).name
            issues.append(
                Issue(
                    target.name,
                    f"connector {name} is not the image of any connector of {source.name}",
                )
            )
    return pairs, issues


def _mismatch_reason(session, image, nodes, best, counts, periodic) -> str:
    """Why a source connector's image matches no target connector: the point counts differ,
    the sides are not images of each other (the ends are off), or the ends coincide and the
    node distributions differ."""
    offset, b, same_direction = best
    if b is None:
        have = ", ".join(str(c) for c in sorted(set(counts)))
        if not counts:
            return "the other side has no connector left to match"
        return (
            f"it has {len(image)} points and the other side's connectors have {have}: give "
            "both sides the same number of points"
        )
    other = nodes[b]
    if not same_direction:
        other = other[::-1]
    ends = max(
        float(np.linalg.norm(other[0] - image[0])), float(np.linalg.norm(other[-1] - image[-1]))
    )
    target_name = session.store.get(b).name
    transform = "translation"
    if periodic["rotation_angle"] != 0.0:
        transform = "rotation"
    if ends > 0.01 * offset + 1e-300 and ends >= offset * 0.5:
        return (
            f"its ends land {ends:.3g} away from the ends of {target_name}: the two sides are "
            f"not images of each other under the {transform} (draw them equal: snap to the "
            "reference grid, or correct the translation / rotation in the patch's panel)"
        )
    return (
        f"its ends coincide with those of {target_name} but the nodes differ by up to "
        f"{offset:.3g}: give both connectors the same distribution"
    )


def automatic_translation(session, source, target) -> tuple[float, float, float]:
    """Translation that maps the centroid of the source nodes onto that of the target (2D: the
    connector nodes; 3D: the nodes of the patches' boundary faces)."""
    if session.is_3d():
        from opengrid.model.periodic_faces import patch_centroid

        delta = patch_centroid(session, target) - patch_centroid(session, source)
        return (float(delta[0]), float(delta[1]), float(delta[2]))

    def centroid(patch):
        points = [
            session.connectors.nodes(i["connector"]).points
            for i in patch.params["members"]
            if "connector" in i
        ]
        if not points:
            raise ValueError(f"patch {patch.name} has no connectors: give the translation")
        return np.concatenate(points).mean(axis=0)

    delta = centroid(target) - centroid(source)
    return (float(delta[0]), float(delta[1]), float(delta[2]))
