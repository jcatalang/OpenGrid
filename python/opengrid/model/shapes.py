"""Standard shapes as geometry specifications (no store): NACA 4-digit airfoils, ellipses as
rational NURBS, rectangles, boxes, cylinders, spheres and tapered swept wings.

Methods and references
- NACA 4-digit sections: thickness y_t = 5 t (0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2
  + 0.2843 x^3 - 0.1036 x^4) (the last coefficient 0.1015 for the open trailing edge) and the
  two-parabola camber line y_c, the surface offset normal to the camber line: Abbott &
  von Doenhoff, "Theory of Wing Sections" (Dover, 1959), §6.4. Points are cosine spaced in
  x, clustering at both edges.
- Ellipses: the rational quadratic NURBS circle of nine control points with weights
  (1, sqrt(2)/2, ...) on a triple-knot vector (Piegl & Tiller, "The NURBS Book", §7.5),
  scaled by the semi-axes; exact.
- Cylinders and spheres: surfaces of revolution of a line and of a meridian arc.
- Wings: ruled surfaces between a root and a tip section (linear taper, sweep, dihedral).
"""

import math


def naca4(code: str, points: int = 101, closed_trailing_edge: bool = True):
    """(upper, lower) point lists of a NACA 4-digit section of unit chord, both from the
    leading edge (0, 0) to the trailing edge (1, 0), ``points`` each, cosine spaced."""
    digits = str(code).strip()
    if len(digits) != 4 or not digits.isdigit():
        raise ValueError(f"a NACA 4-digit code has four digits, got {code!r}")
    if points < 3:
        raise ValueError(f"at least 3 points per side, got {points}")
    m = int(digits[0]) / 100.0
    p = int(digits[1]) / 10.0
    t = int(digits[2:]) / 100.0
    a4 = 0.1036
    if not closed_trailing_edge:
        a4 = 0.1015
    upper = []
    lower = []
    for i in range(points):
        x = 0.5 * (1.0 - math.cos(math.pi * i / (points - 1)))
        yt = (
            5.0
            * t
            * (0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - a4 * x**4)
        )
        if p == 0.0 or m == 0.0:
            yc = 0.0
            slope = 0.0
        elif x < p:
            yc = m / p**2 * (2.0 * p * x - x**2)
            slope = 2.0 * m / p**2 * (p - x)
        else:
            yc = m / (1.0 - p) ** 2 * ((1.0 - 2.0 * p) + 2.0 * p * x - x**2)
            slope = 2.0 * m / (1.0 - p) ** 2 * (p - x)
        theta = math.atan(slope)
        upper.append((x - yt * math.sin(theta), yc + yt * math.cos(theta)))
        lower.append((x + yt * math.sin(theta), yc - yt * math.cos(theta)))
    if closed_trailing_edge:
        # Exact closure: y_t(1) = 0 with the 0.1036 term and y_c(1) = 0, up to rounding.
        upper[-1] = (1.0, 0.0)
        lower[-1] = (1.0, 0.0)
    upper[0] = (0.0, 0.0)
    lower[0] = (0.0, 0.0)
    return upper, lower


def placed(points, chord: float, leading_edge, angle_deg: float):
    """Unit-chord section points scaled by ``chord``, rotated nose-up by the angle of attack
    (the section turns clockwise) and moved to the leading edge (2D, z = 0)."""
    c = math.cos(math.radians(-angle_deg))
    s = math.sin(math.radians(-angle_deg))
    x0, y0 = float(leading_edge[0]), float(leading_edge[1])
    out = []
    for x, y in points:
        xs, ys = chord * x, chord * y
        out.append((x0 + c * xs - s * ys, y0 + s * xs + c * ys, 0.0))
    return out


def ellipse_nurbs(center, a: float, b: float):
    """(degree, knots, control points, weights) of the exact ellipse with semi-axes a (x)
    and b (y) about ``center`` in the z = 0 plane (a == b: a circle)."""
    w = math.sqrt(0.5)
    square = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0)]
    cx, cy = float(center[0]), float(center[1])
    control = [(cx + a * x, cy + b * y, 0.0) for x, y in square]
    weights = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0]
    knots = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]
    return 2, knots, control, weights


def rectangle_corners(center, width: float, height: float):
    """The four corners counter-clockwise from the lower left (z = 0)."""
    cx, cy = float(center[0]), float(center[1])
    hw, hh = 0.5 * float(width), 0.5 * float(height)
    return [
        (cx - hw, cy - hh, 0.0),
        (cx + hw, cy - hh, 0.0),
        (cx + hw, cy + hh, 0.0),
        (cx - hw, cy + hh, 0.0),
    ]


def frame(axis):
    """An orthonormal frame (u, v, w) with w along ``axis``."""
    ax, ay, az = (float(c) for c in axis)
    norm = math.sqrt(ax * ax + ay * ay + az * az)
    if norm == 0.0:
        raise ValueError("the axis must not be zero")
    w = (ax / norm, ay / norm, az / norm)
    # u: the x axis (the y axis for an axis near x) made orthogonal to w, so a z axis gives
    # the usual (x, y, z) frame; v = w x u closes the right-handed triad.
    helper = (1.0, 0.0, 0.0)
    if abs(w[0]) > 0.9:
        helper = (0.0, 1.0, 0.0)
    dot = sum(h * c for h, c in zip(helper, w, strict=True))
    u = tuple(h - dot * c for h, c in zip(helper, w, strict=True))
    norm_u = math.sqrt(sum(c * c for c in u))
    u = (u[0] / norm_u, u[1] / norm_u, u[2] / norm_u)
    v = (w[1] * u[2] - w[2] * u[1], w[2] * u[0] - w[0] * u[2], w[0] * u[1] - w[1] * u[0])
    return u, v, w


def along(point, direction, distance: float):
    return tuple(float(point[k]) + distance * direction[k] for k in range(3))


def box_planes(center, size):
    """The six faces of an axis-aligned box as (suffix, origin, u vector, v vector), outward
    normals u x v."""
    cx, cy, cz = (float(c) for c in center)
    sx, sy, sz = (float(c) for c in size)
    x0, y0, z0 = cx - 0.5 * sx, cy - 0.5 * sy, cz - 0.5 * sz
    x1, y1, z1 = cx + 0.5 * sx, cy + 0.5 * sy, cz + 0.5 * sz
    return [
        ("xmin", (x0, y0, z0), (0.0, 0.0, sz), (0.0, sy, 0.0)),
        ("xmax", (x1, y0, z0), (0.0, sy, 0.0), (0.0, 0.0, sz)),
        ("ymin", (x0, y0, z0), (sx, 0.0, 0.0), (0.0, 0.0, sz)),
        ("ymax", (x0, y1, z0), (0.0, 0.0, sz), (sx, 0.0, 0.0)),
        ("zmin", (x0, y0, z0), (0.0, sy, 0.0), (sx, 0.0, 0.0)),
        ("zmax", (x0, y0, z1), (sx, 0.0, 0.0), (0.0, sy, 0.0)),
    ]


def wing_sections(
    code: str,
    root_chord: float,
    tip_chord: float,
    span: float,
    sweep_deg: float,
    dihedral_deg: float,
    points: int,
):
    """The root and tip sections of a straight tapered wing: chord along +x, span along +y,
    thickness along z; the tip's leading edge swept back by span tan(sweep) and raised by
    span tan(dihedral). Returns {"root": (upper, lower), "tip": (upper, lower)} of 3D
    points, each side from the leading edge to the trailing edge."""
    upper, lower = naca4(code, points)
    tip_x = span * math.tan(math.radians(sweep_deg))
    tip_z = span * math.tan(math.radians(dihedral_deg))

    def section(chord, y, x_offset, z_offset):
        up = [(x_offset + chord * x, y, z_offset + chord * z) for x, z in upper]
        low = [(x_offset + chord * x, y, z_offset + chord * z) for x, z in lower]
        return up, low

    return {
        "root": section(root_chord, 0.0, 0.0, 0.0),
        "tip": section(tip_chord, span, tip_x, tip_z),
    }
