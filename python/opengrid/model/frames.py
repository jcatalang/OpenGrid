"""2D frames of unstructured domains: how a domain's 3D boundary is mapped to the plane in which
the 2D meshers run, and back.

  * z = const (2D models): the identity on (x, y).
  * Planar in another plane (faces of 3D shells): an orthonormal frame of the best-fit plane
    (SVD of the boundary nodes); the boundary must be planar to 1e-9 of its size.
  * On a surface: the parameters (u, v) of the projected boundary nodes, scaled by the mean
    lengths of the tangent vectors |x_u|, |x_v| over the domain (a constant metric correction,
    so sizes in the 2D frame approximate physical sizes; strongly varying metrics are only
    approximated). Interior points are evaluated on the surface.
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


def make_frame(loops: list[np.ndarray], surface=None, where: str = "") -> Frame:
    """Frame of a domain from its boundary loops (3D node arrays, closing node not repeated)."""
    everything = np.concatenate(loops)
    size = float(np.max(np.ptp(everything, axis=0))) or 1.0
    if surface is None and np.ptp(everything[:, 2]) <= PLANARITY * max(size, 1.0):
        return Frame(float(everything[0, 2]))
    if surface is not None:
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
