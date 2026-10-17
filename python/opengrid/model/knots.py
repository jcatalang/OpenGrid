"""Knot insertion and removal of NURBS curves (a control point added or taken away).

Method: Piegl and Tiller, The NURBS Book (2nd ed.), A5.1 (Boehm's knot insertion: one knot
inserted at u gives one more control point and leaves the curve unchanged) and A5.8 (Tiller's
knot removal: the control points of the affected span are recomputed by running the insertion
recurrence backwards from both ends of the span; the removal is exact when the two sweeps
meet, otherwise the curve changes and the caller measures by how much). Both work in
homogeneous coordinates (w x, w y, w z, w), so rational curves are handled.
"""

import numpy as np


def _homogeneous(control_points, weights) -> np.ndarray:
    points = np.asarray(control_points, dtype=float).reshape(-1, 3)
    w = np.asarray(weights, dtype=float).reshape(-1)
    return np.column_stack([points * w[:, None], w])


def _cartesian(pw: np.ndarray) -> tuple[tuple, tuple]:
    w = pw[:, 3]
    points = pw[:, :3] / w[:, None]
    return tuple(tuple(float(c) for c in p) for p in points), tuple(float(v) for v in w)


def find_span(degree: int, knots, u: float) -> int:
    """Index k with knots[k] <= u < knots[k + 1] (the last non-empty span at the end)."""
    n = len(knots) - degree - 1  # control points
    if u >= knots[n]:
        return n - 1
    k = degree
    while knots[k + 1] <= u:
        k += 1
    return k


def multiplicity(knots, u: float, tolerance: float = 0.0) -> int:
    return int(sum(1 for value in knots if abs(value - u) <= tolerance))


def insert_knot(degree: int, knots, control_points, weights, u: float) -> tuple:
    """(knots, control_points, weights) with ``u`` inserted once; the curve is unchanged."""
    knots = [float(value) for value in knots]
    pw = _homogeneous(control_points, weights)
    n = len(pw)
    k = find_span(degree, knots, u)
    out = np.empty((n + 1, 4))
    for i in range(n + 1):
        if i <= k - degree:
            out[i] = pw[i]
        elif i > k:
            out[i] = pw[i - 1]
        else:
            alpha = (u - knots[i]) / (knots[i + degree] - knots[i])
            out[i] = alpha * pw[i] + (1.0 - alpha) * pw[i - 1]
    new_knots = (*knots[: k + 1], float(u), *knots[k + 1 :])
    points, w = _cartesian(out)
    return new_knots, points, w


def removable_control_points(degree: int, knots) -> dict[int, int]:
    """{control point index: knot index r} of the control points a single knot removal takes
    away: every distinct interior knot value (r = its last index, s = its multiplicity)
    removes control point (2 r - s - degree) // 2. The end control points are never removable
    (they fix the curve's ends)."""
    knots = [float(value) for value in knots]
    n = len(knots) - degree - 1
    out = {}
    r = degree + 1
    while r < n:
        value = knots[r]
        last = r
        while last + 1 < n and knots[last + 1] == value:
            last += 1
        s = multiplicity(knots, value)
        out[(2 * last - s - degree) // 2] = last
        r = last + 1
    return out


def remove_knot(degree: int, knots, control_points, weights, r: int) -> tuple:
    """(knots, control_points, weights) with knot ``r`` (the last index of its value) removed
    once; exact only when the curve allows it (the caller compares the curves)."""
    knots = [float(value) for value in knots]
    pw = _homogeneous(control_points, weights)
    u = knots[r]
    s = multiplicity(knots, u)
    order = degree + 1
    fout = (2 * r - s - degree) // 2
    last = r - s
    first = r - degree
    temp = np.zeros((2 * degree + 1, 4))
    off = first - 1
    temp[0] = pw[off]
    temp[last + 1 - off] = pw[last + 1]
    i = first
    j = last
    ii = 1
    jj = last - off
    while j - i > 0:
        alpha_i = (u - knots[i]) / (knots[i + order] - knots[i])
        alpha_j = (u - knots[j]) / (knots[j + order] - knots[j])
        temp[ii] = (pw[i] - (1.0 - alpha_i) * temp[ii - 1]) / alpha_i
        temp[jj] = (pw[j] - alpha_j * temp[jj + 1]) / (1.0 - alpha_j)
        i += 1
        ii += 1
        j -= 1
        jj -= 1
    i = first
    j = last
    while j - i > 0:
        pw[i] = temp[i - off]
        pw[j] = temp[j - off]
        i += 1
        j -= 1
    del knots[r]
    pw = np.delete(pw, fout, axis=0)
    points, w = _cartesian(pw)
    return tuple(knots), points, w
