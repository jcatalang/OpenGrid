"""Periodic pairing of patches.

A periodic pair is stored on its *source* patch: ``periodic = {"partner": target patch id,
"rotation_center", "rotation_angle", "translation"}``. The transform maps source points onto
target points as a rotation about the z axis through ``rotation_center`` by ``rotation_angle``
(radians) followed by ``translation`` (the CGNS GridConnectivityProperty convention):
    x_target = R(angle) (x_source - center) + center + translation.
Every connector of the source must map node by node (in either direction) onto a distinct
connector of the target.
"""

import math
from dataclasses import dataclass

import numpy as np

from opengrid.validation import Issue

MATCH_TOLERANCE = 1e-8  # relative to the model size


@dataclass(frozen=True)
class ConnectorPair:
    source: int  # connector id
    target: int
    same_direction: bool  # source first node -> target first node


def apply(points: np.ndarray, periodic: dict) -> np.ndarray:
    c = np.asarray(periodic["rotation_center"], dtype=float)
    t = np.asarray(periodic["translation"], dtype=float)
    a = periodic["rotation_angle"]
    rotation = np.array(
        [[math.cos(a), -math.sin(a), 0.0], [math.sin(a), math.cos(a), 0.0], [0.0, 0.0, 1.0]]
    )
    return (points - c) @ rotation.T + c + t


def inverse(periodic: dict) -> dict:
    """Transform from target back to source in the same convention."""
    c = np.asarray(periodic["rotation_center"], dtype=float)
    t = np.asarray(periodic["translation"], dtype=float)
    return {
        "rotation_center": tuple(float(v) for v in c + t),
        "rotation_angle": -periodic["rotation_angle"],
        "translation": tuple(float(v) for v in -t),
    }


def affine_matrix(periodic: dict) -> np.ndarray:
    """4 x 4 homogeneous matrix of the transform."""
    c = np.asarray(periodic["rotation_center"], dtype=float)
    t = np.asarray(periodic["translation"], dtype=float)
    a = periodic["rotation_angle"]
    rotation = np.array(
        [[math.cos(a), -math.sin(a), 0.0], [math.sin(a), math.cos(a), 0.0], [0.0, 0.0, 1.0]]
    )
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
            points = session.connectors.nodes(item["connector"]).points
            nodes[item["connector"]] = points
            scale = max(scale, float(np.max(np.abs(points))))
    tolerance = MATCH_TOLERANCE * scale
    pairs = []
    issues = []
    used = set()
    for item in source.params["members"]:
        a = item["connector"]
        image = apply(nodes[a], periodic)
        best = math.inf
        found = None
        for candidate in target.params["members"]:
            b = candidate["connector"]
            if b in used or len(nodes[b]) != len(image):
                continue
            forward = float(np.max(np.linalg.norm(nodes[b] - image, axis=1)))
            backward = float(np.max(np.linalg.norm(nodes[b][::-1] - image, axis=1)))
            best = min(best, forward, backward)
            if forward <= tolerance:
                found = ConnectorPair(a, b, True)
                break
            if backward <= tolerance:
                found = ConnectorPair(a, b, False)
                break
        name = session.store.get(a).name
        if found is None:
            detail = "no target connector with the same number of points"
            if math.isfinite(best):
                detail = f"closest target connector is {best:.3g} away"
            issues.append(
                Issue(
                    source.name,
                    f"connector {name} has no periodic image in {target.name} ({detail})",
                )
            )
            continue
        used.add(found.target)
        pairs.append(found)
    for item in target.params["members"]:
        if item["connector"] not in used and not issues:
            name = session.store.get(item["connector"]).name
            issues.append(
                Issue(
                    target.name,
                    f"connector {name} is not the image of any connector of {source.name}",
                )
            )
    return pairs, issues


def automatic_translation(session, source, target) -> tuple[float, float, float]:
    """Translation that maps the centroid of the source nodes onto that of the target."""

    def centroid(patch):
        points = np.concatenate(
            [session.connectors.nodes(i["connector"]).points for i in patch.params["members"]]
        )
        return points.mean(axis=0)

    delta = centroid(target) - centroid(source)
    return (float(delta[0]), float(delta[1]), float(delta[2]))
