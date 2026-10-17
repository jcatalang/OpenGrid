"""2D frames of unstructured domains: how a domain's 3D boundary is mapped to the plane in which
the 2D meshers run, and back.

  * z = const (2D models): the identity on (x, y).
  * Planar in another plane (faces of 3D shells): an orthonormal frame of the best-fit plane
    (SVD of the boundary nodes); the boundary must be planar to 1e-9 of its size.
  * On a surface: the parameters (u, v) of the projected boundary nodes, scaled by the mean
    lengths of the tangent vectors |x_u|, |x_v| over the domain (a constant metric correction,
    so sizes in the 2D frame approximate physical sizes; strongly varying metrics are only
    approximated). Interior points are evaluated on the surface.
  * On a surface of revolution with a pole (the generating curve meets the axis, as a
    sphere's halves do): a conformal map, so the pole is an ordinary boundary point of the
    plane and triangles keep their shape. With s the arc length along the generating curve
    and r(s) its distance to the axis, sigma = integral ds / r is the isothermal coordinate
    (Gauss; for a sphere the Mercator latitude), w = exp(sigma + i theta) the conformal
    plane (for a sphere, the stereographic projection from the far pole: the near pole at 0,
    the far one at infinity), and z = (w - 1) / (w + 1) a Mobius map that brings the far
    pole to +1, the near pole to -1 and the meridian opposite the patch's middle to
    infinity, so a patch of less than a full turn is a bounded lens; z is scaled to 3D
    length units on average. The map is exact in both directions (sigma is tabulated along
    the curve and inverted by interpolation).
Planar and surface frames are oriented so that the outer loop (largest enclosed area) runs
counter-clockwise; z = const frames keep the (x, y) orientation (2D convention).
"""

import numpy as np

from opengrid.validation import Issue, ValidationError

PLANARITY = 1e-9


def _signed_area(points2: np.ndarray) -> float:
    x = points2[:, 0]
    y = points2[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


class Frame:
    kind = "z"

    def __init__(self, z0: float):
        self.z0 = z0

    def to2d(self, points: np.ndarray) -> np.ndarray:
        return np.ascontiguousarray(np.asarray(points, dtype=float)[:, :2])

    def to3d(self, uv: np.ndarray) -> np.ndarray:
        out = np.zeros((len(uv), 3))
        out[:, :2] = uv
        out[:, 2] = self.z0
        return out


class PlaneFrame(Frame):
    kind = "plane"

    def __init__(self, origin, e1, e2):
        self.origin = np.asarray(origin, dtype=float)
        self.e1 = np.asarray(e1, dtype=float)
        self.e2 = np.asarray(e2, dtype=float)

    def to2d(self, points):
        d = np.asarray(points, dtype=float) - self.origin
        return np.ascontiguousarray(np.stack([d @ self.e1, d @ self.e2], axis=1))

    def to3d(self, uv):
        uv = np.asarray(uv, dtype=float)
        return self.origin + np.outer(uv[:, 0], self.e1) + np.outer(uv[:, 1], self.e2)

    def flip(self) -> None:
        self.e2 = -self.e2


class SurfaceFrame(Frame):
    kind = "surface"

    def __init__(self, surface, points):
        self.surface = surface
        self.flip_v = 1.0
        uv = []
        for p in np.asarray(points, dtype=float):
            u, v, _q, _d = surface.project(tuple(float(c) for c in p))
            uv.append((u, v))
        uv = np.asarray(uv)
        self._cache = {
            tuple(np.asarray(p, dtype=float).tolist()): tuple(q)
            for p, q in zip(points, uv, strict=True)
        }
        lo = uv.min(axis=0)
        hi = uv.max(axis=0)
        su = []
        sv = []
        for a in np.linspace(0.0, 1.0, 5):
            for b in np.linspace(0.0, 1.0, 5):
                u = lo[0] + a * (hi[0] - lo[0])
                v = lo[1] + b * (hi[1] - lo[1])
                h = 1e-6
                xu = np.array(surface.eval(min(u + h, 1.0), v)) - np.array(
                    surface.eval(max(u - h, 0.0), v)
                )
                xv = np.array(surface.eval(u, min(v + h, 1.0))) - np.array(
                    surface.eval(u, max(v - h, 0.0))
                )
                su.append(np.linalg.norm(xu) / (min(u + h, 1.0) - max(u - h, 0.0)))
                sv.append(np.linalg.norm(xv) / (min(v + h, 1.0) - max(v - h, 0.0)))
        self.su = float(np.mean(su)) or 1.0
        self.sv = float(np.mean(sv)) or 1.0

    def to2d(self, points):
        out = []
        for p in np.asarray(points, dtype=float):
            key = tuple(p.tolist())
            if key in self._cache:
                u, v = self._cache[key]
            else:
                u, v, _q, _d = self.surface.project(key)
            out.append((u * self.su, self.flip_v * v * self.sv))
        return np.ascontiguousarray(np.asarray(out, dtype=float))

    def to3d(self, uv):
        uv = np.asarray(uv, dtype=float)
        out = np.zeros((len(uv), 3))
        for k, (a, b) in enumerate(uv):
            u = min(max(a / self.su, 0.0), 1.0)
            v = min(max(self.flip_v * b / self.sv, 0.0), 1.0)
            out[k] = self.surface.eval(u, v)
        return out

    def flip(self) -> None:
        self.flip_v = -self.flip_v


class RevolutionFrame(Frame):
    """The conformal plane of a surface of revolution with a pole (module docstring)."""

    kind = "revolution"
    TABLE = 4001

    def __init__(self, surface, axis_point, axis, start_angle, end_angle, curve):
        self.surface = surface
        self.axis_point = np.asarray(axis_point, dtype=float)
        w = np.asarray(axis, dtype=float)
        self.w = w / np.linalg.norm(w)
        self.a0 = float(start_angle)
        self.span = float(end_angle) - self.a0
        self.sign = 1.0
        ts = np.linspace(0.0, 1.0, self.TABLE)
        points = np.asarray([curve.eval(float(t)) for t in ts], dtype=float)
        rel = points - self.axis_point
        height = rel @ self.w
        r = np.linalg.norm(rel - np.outer(height, self.w), axis=1)
        scale = max(1.0, float(np.max(np.abs(points))))
        self.pole0 = r[0] <= 1e-9 * scale
        self.pole1 = r[-1] <= 1e-9 * scale
        ds = np.linalg.norm(np.diff(points, axis=0), axis=1)
        r_mid = 0.5 * (r[:-1] + r[1:])
        sigma = np.concatenate([[0.0], np.cumsum(ds / np.maximum(r_mid, 1e-300))])
        sigma -= sigma[self.TABLE // 2]
        rho = np.exp(sigma)
        if self.pole0:
            rho[0] = 0.0
        self.ts = ts
        self.rho = rho
        # 3D length per unit of the conformal plane, averaged (geometric mean) over the patch.
        self.length = 1.0
        ratios = []
        for u in np.linspace(0.1, 0.9, 9):
            for v in np.linspace(0.05, 0.95, 7):
                p = np.asarray(surface.eval(float(u), float(v)))
                q = np.asarray(surface.eval(float(u) + 1e-4, float(v)))
                d3 = float(np.linalg.norm(q - p))
                d2 = float(np.linalg.norm(self._plane(u + 1e-4, v) - self._plane(u, v)))
                if d3 > 0 and d2 > 0:
                    ratios.append(d3 / d2)
        if ratios:
            self.length = float(np.exp(np.mean(np.log(ratios))))

    def _rho(self, u: float) -> float:
        if self.pole1 and u >= self.ts[-2]:
            # rho ~ 1 / (1 - t) at a pole at the curve's end.
            remaining = max(1.0 - float(u), 1e-300)
            return float(self.rho[-2]) * (1.0 - float(self.ts[-2])) / remaining
        return float(np.interp(u, self.ts, self.rho))

    def _u(self, rho: float) -> float:
        if self.pole1 and rho >= self.rho[-2]:
            return 1.0 - (1.0 - float(self.ts[-2])) * float(self.rho[-2]) / max(rho, 1e-300)
        return float(np.interp(rho, self.rho[:-1], self.ts[:-1]))

    def _plane(self, u: float, v: float) -> np.ndarray:
        """The unscaled conformal plane point of surface parameters (u, v)."""
        theta = (float(v) - 0.5) * self.span
        w = self._rho(u) * complex(np.cos(theta), np.sin(theta))
        if not np.isfinite(w.real) or abs(w) > 1e150:
            return np.asarray([1.0, 0.0])
        z = (w - 1.0) / (w + 1.0)
        return np.asarray([z.real, z.imag])

    def to2d(self, points):
        out = []
        for p in np.asarray(points, dtype=float):
            u, v, _q, _d = self.surface.project(tuple(float(c) for c in p))
            x, y = self._plane(u, v) * self.length
            out.append((x, self.sign * y))
        return np.ascontiguousarray(np.asarray(out, dtype=float))

    def to3d(self, uv):
        uv = np.asarray(uv, dtype=float)
        out = np.zeros((len(uv), 3))
        for k, (x, y) in enumerate(uv):
            z = complex(x, self.sign * y) / self.length
            if abs(z - 1.0) < 1e-12:
                u, v = 1.0, 0.5
            else:
                w = (1.0 + z) / (1.0 - z)
                u = min(max(self._u(abs(w)), 0.0), 1.0)
                theta = float(np.angle(w))
                v = min(max(theta / self.span + 0.5, 0.0), 1.0)
            out[k] = self.surface.eval(u, v)
        return out

    def flip(self) -> None:
        self.sign = -self.sign


def revolution_frame(store, entity, surface):
    """A RevolutionFrame for a revolved surface whose generating curve meets the axis at an
    end (a pole) over less than a full turn, else None."""
    if entity is None or entity.kind != "revolved_surface":
        return None
    p = entity.params
    span = float(p["end_angle"]) - float(p["start_angle"])
    if abs(abs(span) - 2.0 * np.pi) < 1e-9 or span == 0.0:
        return None
    curve = store.geometry(store.get(p["curve"]).name)
    axis_point = np.asarray(p["axis_point"], dtype=float)
    w = np.asarray(p["axis"], dtype=float)
    w = w / np.linalg.norm(w)
    ends = [np.asarray(curve.eval(0.0)), np.asarray(curve.eval(1.0))]
    scale = max(1.0, float(np.max(np.abs(ends))))
    radii = []
    for end in ends:
        rel = end - axis_point
        radii.append(float(np.linalg.norm(rel - float(rel @ w) * w)))
    if min(radii) > 1e-9 * scale:
        return None
    return RevolutionFrame(surface, axis_point, w, p["start_angle"], p["end_angle"], curve)


def make_frame(
    loops: list[np.ndarray], surface=None, where: str = "", surface_entity=None, store=None
) -> Frame:
    """Frame of a domain from its boundary loops (3D node arrays, closing node not repeated);
    with the surface's entity and the store, a surface of revolution with a pole gets its
    conformal frame."""
    everything = np.concatenate(loops)
    size = float(np.max(np.ptp(everything, axis=0))) or 1.0
    if surface is None and np.ptp(everything[:, 2]) <= PLANARITY * max(size, 1.0):
        return Frame(float(everything[0, 2]))
    if surface is not None:
        frame = None
        if store is not None:
            frame = revolution_frame(store, surface_entity, surface)
        if frame is None:
            frame = SurfaceFrame(surface, everything)
    else:
        centre = everything.mean(axis=0)
        _u, s, vt = np.linalg.svd(everything - centre)
        if s[-1] > PLANARITY * max(size, 1.0) * np.sqrt(len(everything)):
            raise ValidationError(
                [
                    Issue(
                        where,
                        "the boundary is not planar: unstructured domains in 3D must be planar "
                        "or lie on a surface (set 'surface')",
                    )
                ]
            )
        frame = PlaneFrame(centre, vt[0], vt[1])
    areas = [_signed_area(frame.to2d(loop)) for loop in loops]
    outer = int(np.argmax([abs(a) for a in areas]))
    if areas[outer] < 0:
        frame.flip()
    return frame
