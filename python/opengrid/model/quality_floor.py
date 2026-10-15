"""The quality floor of 3D cells: every tetrahedron, pyramid, prism and hexahedron of a volume or
block must have a scaled Jacobian of at least ``min_jacobian`` and an equiangle skewness of at
most ``max_skewness`` (definitions in core quality/quality_cells3d.hpp and
opengrid.model.volumes.volume_quality).

Volumes enforce it while they are built (layer cells that would miss it are not made; the
optimiser rescues core cells below it); whatever still misses it, in volumes or blocks, is
reported (status, warnings, boundary report), so exports refuse such a mesh. The default (0.05,
0.95) follows common solver guidance (skewness above 0.95 is "bad").
"""

import numpy as np

from opengrid import _core
from opengrid.validation import Issue

DEFAULT_QUALITY_FLOOR = {"min_jacobian": 0.05, "max_skewness": 0.95}


def normalize_quality_floor(raw, label: str, issues: list) -> dict:
    """Canonical floor from None (the default) or a dict with some of its keys."""
    floor = dict(DEFAULT_QUALITY_FLOOR)
    if raw is None:
        return floor
    if not isinstance(raw, dict):
        issues.append(Issue(label, "expected a dict with min_jacobian and max_skewness"))
        return floor
    for key in sorted(set(raw) - set(DEFAULT_QUALITY_FLOOR)):
        issues.append(Issue(f"{label}.{key}", "unexpected key"))
    checks = (
        ("min_jacobian", lambda v: 0.0 <= v < 1.0, "in [0, 1)"),
        ("max_skewness", lambda v: 0.0 < v <= 1.0, "in (0, 1]"),
    )
    for key, valid, text in checks:
        if key not in raw:
            continue
        value = raw[key]
        if hasattr(value, "item") and not isinstance(value, (int, float)):
            value = value.item()
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not valid(value):
            issues.append(Issue(f"{label}.{key}", f"must be {text}, got {value!r}"))
            continue
        floor[key] = float(value)
    return floor


def floor_summary(scaled_jacobian, skewness, floor: dict) -> dict:
    """{"floor", "cells", "below_floor", "worst_jacobian", "worst_skewness", "failing" (cell
    indices)} of per-cell metric arrays."""
    sj = np.asarray(scaled_jacobian, dtype=float).reshape(-1)
    sk = np.asarray(skewness, dtype=float).reshape(-1)
    failing = np.flatnonzero((sj < floor["min_jacobian"]) | (sk > floor["max_skewness"]))
    worst_jacobian = 1.0
    worst_skewness = 0.0
    if len(sj):
        worst_jacobian = float(sj.min())
        worst_skewness = float(sk.max())
    return {
        "floor": dict(floor),
        "cells": len(sj),
        "below_floor": len(failing),
        "worst_jacobian": worst_jacobian,
        "worst_skewness": worst_skewness,
        "failing": failing,
    }


def block_floor(grid, floor: dict) -> dict:
    """Floor check of a block grid (nk, nj, ni, 3) with its hexahedra's metrics."""
    q = _core.quality_3d(np.ascontiguousarray(grid, dtype=float))
    return floor_summary(q["scaled_jacobian"], q["skewness"], floor)


def floor_issue(name: str, summary: dict) -> Issue | None:
    """An issue when cells miss the floor, else None."""
    if summary["below_floor"] == 0:
        return None
    floor = summary["floor"]
    return Issue(
        name,
        f"{summary['below_floor']} of {summary['cells']} cells miss the quality floor (scaled "
        f"Jacobian >= {floor['min_jacobian']:g}, skewness <= {floor['max_skewness']:g}); worst "
        f"scaled Jacobian {summary['worst_jacobian']:.3g}, skewness "
        f"{summary['worst_skewness']:.3g}",
    )


CELL_KEYS = ("tets", "pyramids", "prisms", "hexes")
