"""Freehand strokes to NURBS curves.

Method: the stroke (pointer samples) is split at corners - samples where the direction over a
window of 3 % of the stroke length on either side turns by at least ``corner_angle`` degrees,
the sharpest sample of each run - and every piece is fitted by a clamped cubic B-spline (unit
weights): least squares over chord-length parameters with the end points interpolated, a
small second-difference penalty against wiggles from pointer noise, and knot insertion at the
middle of the span with the largest error until every sample lies within ``tolerance`` (a
fraction of the stroke length) or the control-point budget is used. The parameters are then
corrected twice by projecting the samples onto the curve (Hoschek's parameter correction).
References: L. Piegl, W. Tiller, The NURBS Book, 2nd ed. (1997), ch. 9.4 (least-squares
approximation, eq. 9.63-9.67); J. Hoschek, "Intrinsic parametrization for approximation",
Computer Aided Geometric Design 5 (1988) 27-31.
"""

import itertools
import math

import numpy as np

DEGREE = 3


def _basis(knots: np.ndarray, degree: int, u: np.ndarray) -> np.ndarray:
    """B-spline basis matrix (len(u), number of control points), Cox-de Boor recursion."""
    n = len(knots) - degree - 1
    u = np.clip(u, knots[degree], knots[-degree - 1])
    # Degree 0: the span containing u (the last span includes its right end).
    basis = np.zeros((len(u), len(knots) - 1))
    for k in range(len(knots) - 1):
        if knots[k] < knots[k + 1]:
            inside = (u >= knots[k]) & (u < knots[k + 1])
            if knots[k + 1] == knots[-degree - 1]:
                inside |= u == knots[k + 1]
            basis[inside, k] = 1.0
    for p in range(1, degree + 1):
        nxt = np.zeros((len(u), len(knots) - p - 1))
        for k in range(len(knots) - p - 1):
            left = knots[k + p] - knots[k]
            right = knots[k + p + 1] - knots[k + 1]
            term = np.zeros(len(u))
            if left > 0:
                term += (u - knots[k]) / left * basis[:, k]
            if right > 0:
                term += (knots[k + p + 1] - u) / right * basis[:, k + 1]
            nxt[:, k] = term
        basis = nxt
    return basis[:, :n]


def _clean(points) -> np.ndarray:
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    keep = [0]
    for k in range(1, len(pts)):
        if np.linalg.norm(pts[k] - pts[keep[-1]]) > 0.0:
            keep.append(k)
    return pts[keep]


def _arc(points: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))])


def find_corners(points, corner_angle: float = 60.0, window: float = 0.03) -> list[int]:
    """Indices of corner samples (end points excluded)."""
    pts = _clean(points)
    s = _arc(pts)
    if len(pts) < 5 or s[-1] == 0.0:
        return []
    w = window * s[-1]
    angles = np.zeros(len(pts))
    for k in range(1, len(pts) - 1):
        back = int(np.searchsorted(s, s[k] - w, side="right")) - 1
        ahead = int(np.searchsorted(s, s[k] + w))
        back = max(back, 0)
        ahead = min(ahead, len(pts) - 1)
        if back == k or ahead == k:
            continue
        u = pts[k] - pts[back]
        v = pts[ahead] - pts[k]
        norm = float(np.linalg.norm(u) * np.linalg.norm(v))
        if norm == 0.0:
            continue
        angles[k] = math.degrees(math.acos(float(np.clip(np.dot(u, v) / norm, -1.0, 1.0))))
    corners = []
    k = 1
    while k < len(pts) - 1:
        if angles[k] < corner_angle:
            k += 1
            continue
        run_end = k
        while run_end + 1 < len(pts) - 1 and angles[run_end + 1] >= corner_angle:
            run_end += 1
        best = k + int(np.argmax(angles[k : run_end + 1]))
        # Corners closer than a window to an end point or to each other are dropped.
        if s[best] > w and s[-1] - s[best] > w and (not corners or s[best] - s[corners[-1]] > w):
            corners.append(best)
        k = run_end + 1
    return corners


def _uniform_knots(spans: int, degree: int = DEGREE) -> np.ndarray:
    inner = np.linspace(0.0, 1.0, spans + 1)
    return np.concatenate([[0.0] * degree, inner, [1.0] * degree])


def _solve(points: np.ndarray, u: np.ndarray, knots: np.ndarray, degree: int, smooth: float):
    """Control points: least squares with the end points fixed, plus a second-difference
    penalty ``smooth`` (relative to the data term)."""
    a = _basis(knots, degree, u)
    m = a.shape[1]
    control = np.zeros((m, 3))
    control[0] = points[0]
    control[-1] = points[-1]
    if m <= 2:
        return control
    rhs = points - np.outer(a[:, 0], points[0]) - np.outer(a[:, -1], points[-1])
    inner = a[:, 1:-1]
    lhs = inner.T @ inner
    second = np.zeros((max(m - 2, 0), m))
    for k in range(m - 2):
        second[k, k : k + 3] = (1.0, -2.0, 1.0)
    if len(second):
        scale = smooth * float(np.trace(lhs)) / max(1, m)
        d_inner = second[:, 1:-1]
        lhs = lhs + scale * d_inner.T @ d_inner
        rhs_pen = (
            -scale
            * d_inner.T
            @ (np.outer(second[:, 0], points[0]) + np.outer(second[:, -1], points[-1]))
        )
    else:
        rhs_pen = 0.0
    control[1:-1] = np.linalg.solve(lhs, inner.T @ rhs + rhs_pen)
    return control


def _evaluate(knots, control, degree, u):
    return _basis(knots, degree, u) @ control


def _project(points, knots, control, degree, u):
    """Parameter correction: one Newton step per sample on (C(u) - p) . C'(u) = 0."""
    h = 1e-6
    c = _evaluate(knots, control, degree, u)
    cp = (
        _evaluate(knots, control, degree, np.clip(u + h, 0, 1))
        - _evaluate(knots, control, degree, np.clip(u - h, 0, 1))
    ) / (2 * h)
    num = np.sum((c - points) * cp, axis=1)
    den = np.sum(cp * cp, axis=1)
    step = np.zeros_like(u)
    ok = den > 0
    step[ok] = num[ok] / den[ok]
    out = np.clip(u - step, 0.0, 1.0)
    out[0] = 0.0
    out[-1] = 1.0
    return np.maximum.accumulate(out)


def fit_piece(points, tolerance: float, max_control: int = 64, smooth: float = 1e-6) -> dict:
    """One clamped cubic B-spline through the end points of ``points``, within the absolute
    ``tolerance`` where the budget allows. Returns create_nurbs parameters plus "error"."""
    pts = _clean(points)
    if len(pts) < 2:
        raise ValueError("a stroke needs two distinct points")
    s = _arc(pts)
    u = s / s[-1]
    degree = min(DEGREE, len(pts) - 1)
    spans = 1
    knots = _uniform_knots(spans, degree)
    while True:
        control = _solve(pts, u, knots, degree, smooth)
        for _ in range(2):
            u = _project(pts, knots, control, degree, u)
            control = _solve(pts, u, knots, degree, smooth)
        error = np.linalg.norm(_evaluate(knots, control, degree, u) - pts, axis=1)
        worst = float(error.max())
        if worst <= tolerance or len(control) >= max_control or len(control) >= len(pts):
            break
        # Insert a knot in the middle of the span holding the largest error.
        k = int(np.argmax(error))
        span = int(np.searchsorted(knots, u[k], side="right")) - 1
        span = min(max(span, degree), len(knots) - degree - 2)
        knots = np.sort(np.append(knots, 0.5 * (knots[span] + knots[span + 1])))
    return {
        "degree": degree,
        "knots": tuple(float(v) for v in knots),
        "control_points": tuple(tuple(float(c) for c in p) for p in control),
        "weights": tuple(1.0 for _ in control),
        "error": worst,
    }


def fit_stroke(
    points,
    tolerance: float = 0.005,
    corner_angle: float | None = 60.0,
    max_control: int = 64,
    close_tolerance: float = 0.0,
) -> list[dict]:
    """NURBS pieces for a pointer stroke: split at corners (``corner_angle`` degrees; None:
    no splitting) and fitted within ``tolerance`` times the stroke length. A stroke that ends
    within ``close_tolerance`` (absolute; 0: never) of its start is closed first: its start
    point is appended, so the fitted pieces form a loop whose ends meet exactly (strokes
    shorter than four tolerances are left alone)."""
    pts = _clean(points)
    if len(pts) < 2:
        raise ValueError("a stroke needs two distinct points")
    length = float(_arc(pts)[-1])
    if close_tolerance > 0.0 and length > 4.0 * close_tolerance:
        gap = float(np.linalg.norm(pts[-1] - pts[0]))
        if 0.0 < gap <= close_tolerance:
            pts = np.vstack([pts, pts[:1]])
            length = float(_arc(pts)[-1])
    cuts = []
    if corner_angle is not None:
        cuts = find_corners(pts, corner_angle)
    bounds = [0, *cuts, len(pts) - 1]
    pieces = []
    for a, b in itertools.pairwise(bounds):
        pieces.append(fit_piece(pts[a : b + 1], tolerance * length, max_control))
    return pieces
