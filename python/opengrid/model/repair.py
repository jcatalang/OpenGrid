"""Repair of curved high-order elements and rigorous validity (Phase 12c).

Validity. The Jacobian determinant J of an order-p element is a polynomial on the reference
element: degree 2(p - 1) on triangles, 3(p - 1) on tetrahedra, 2p - 1 per direction on
quadrilaterals, 3p - 1 per direction on hexahedra, triangle degree 3p - 2 times line degree
3p - 1 on prisms. Interpolated at the lattice of that degree it is represented exactly, and
its Bernstein-Bezier coefficients bound it on the element (convex hull): min B <= J <= max B.
All lattice values positive but a negative coefficient is undetermined; the lattice of a
higher degree (degree elevation) tightens the bound, three rounds here. Pyramids have a
rational Jacobian and are checked by sampling only.

Repair. Elements that fold or lose most of their straight-sided quality after curving are
repaired by moving their free high-order nodes (not on the geometry, not on the boundary,
never the linear mesh's vertices) to minimise, over them and their neighbours (whose other
nodes stay), the distortion sum_q (r - 1)^2 + alpha / h(r) with r = J / J_straight at the
sample points and h(r) = (r + sqrt(r^2 + delta^2)) / 2 the regularised Jacobian of Escobar
et al., which untangles folded elements (delta is reduced in rounds). L-BFGS with an Armijo
line search, analytic gradients (d det A / dA = cof A). The result is kept when the worst
element improved.

References: A. Johnen, J.-F. Remacle, C. Geuzaine, "Geometrical validity of curvilinear
finite elements", J. Comput. Phys. 233 (2013) 359-372; T. Toulorge, C. Geuzaine, J.-F.
Remacle, J. Lambrechts, "Robust untangling of curvilinear meshes", J. Comput. Phys. 254
(2013) 8-26; J. M. Escobar et al., "Simultaneous untangling and smoothing of tetrahedral
meshes", Comput. Methods Appl. Mech. Engrg. 192 (2003) 2775-2787.
"""

from functools import cache
from math import comb, factorial

import numpy as np

from opengrid.io import lagrange
from opengrid.model import shape

ELEVATIONS = (0, 2, 4)  # extra degrees tried when the Bezier bound is inconclusive
REPAIR_RATIO = 0.25  # repair elements whose scaled Jacobian is below this ...
REPAIR_LOSS = 0.7  # ... or below this fraction of their straight-sided one
BEZIER_RATIO = 0.3  # rigorous check for elements whose sampled scaled Jacobian is below this
ALPHA = 0.05  # weight of the barrier in the distortion
DELTAS = (0.1, 0.01, 0.001)  # regularisation rounds
MAX_ITERATIONS = 60  # per round
MAX_EVALUATIONS = 400  # energy evaluations in all (a hopeless patch stops early)


# ------------------------------------------------------------------ Bernstein lattices
@cache
def line_lattice(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Lobatto points of degree n on [0, 1] (the equispaced Bernstein Vandermonde is too
    ill-conditioned beyond degree 10) and the matrix taking function values on them to
    Bernstein coefficients."""
    if n == 0:
        return np.zeros(1), np.ones((1, 1))
    from opengrid import _core

    t = (np.asarray(_core.nodes_1d("gll", n + 1)) + 1.0) / 2.0
    m = np.array([[comb(n, i) * x**i * (1 - x) ** (n - i) for i in range(n + 1)] for x in t])
    return t, np.linalg.inv(m)


@cache
def simplex_lattice(d: int, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Lattice of degree n on the unit simplex of dimension d (points (m, d)) and the matrix
    taking function values on it to Bernstein coefficients. The lattice is equispaced,
    warped towards the Gauss-Lobatto points along every barycentric direction (the 1D warp
    of Warburton), which keeps the Vandermonde well conditioned at the degrees used."""
    if n == 0:
        return np.full((1, d), 1.0 / (d + 1)), np.ones((1, 1))
    indices = []

    def recurse(prefix, remaining, depth):
        if depth == d:
            indices.append((*prefix, remaining))
            return
        for k in range(remaining + 1):
            recurse((*prefix, k), remaining - k, depth + 1)

    recurse((), n, 0)  # (i_1 .. i_d, i_0): i_0 for the barycentric coordinate 1 - sum
    gll, _ = line_lattice(n)
    # Every lattice point moves to the GLL position of each of its barycentric lattice
    # coordinates, renormalised to sum to one.
    pts = []
    for index in indices:
        lam = np.array([index[-1], *index[:d]], dtype=float) / n
        w = np.array([gll[round(value * n)] for value in lam])
        w = w / w.sum()
        pts.append(w[1:])
    points = np.array(pts)
    rows = []
    for point in points:
        lam = [1.0 - point.sum(), *point]
        row = []
        for index in indices:
            multi = (index[-1], *index[:d])
            coefficient = factorial(n)
            for k in multi:
                coefficient //= factorial(k)
            value = float(coefficient)
            for k, coordinate in zip(multi, lam, strict=True):
                value *= coordinate**k
            row.append(value)
        rows.append(row)
    return points, np.linalg.inv(np.array(rows))


def jacobian_degree(kind: str, p: int) -> tuple:
    """Degrees of the Jacobian determinant polynomial: (line degree,) per direction for tensor
    kinds, (total degree,) for simplices, (triangle degree, line degree) for prisms."""
    if kind == "triangle":
        return (2 * (p - 1),)
    if kind == "tetra":
        return (3 * (p - 1),)
    if kind == "quad":
        return (2 * p - 1,)
    if kind == "hexa":
        return (3 * p - 1,)
    if kind == "prism":
        return (3 * p - 2, 3 * p - 1)
    raise ValueError(f"no polynomial Jacobian for {kind}")


@cache
def _lattice(kind: str, degrees: tuple) -> tuple[np.ndarray, tuple]:
    """Lattice points (m, d) of the given degrees and the conversion operators."""
    if kind in ("triangle", "tetra"):
        points, matrix = simplex_lattice(lagrange.DIM[kind], degrees[0])
        return points, (matrix,)
    if kind in ("quad", "hexa"):
        t, matrix = line_lattice(degrees[0])
        d = lagrange.DIM[kind]
        grid = np.meshgrid(*([t] * d), indexing="ij")
        points = np.column_stack([g.reshape(-1) for g in grid])  # last axis fastest
        return points, (matrix,)
    tri, tri_matrix = simplex_lattice(2, degrees[0])
    t, line_matrix = line_lattice(degrees[1])
    points = np.column_stack(
        [np.repeat(tri, len(t), axis=0), np.tile(t, len(tri))]
    )  # triangle outer, line inner
    return points, (tri_matrix, line_matrix)


def bezier_coefficients(kind: str, degrees: tuple, values: np.ndarray) -> np.ndarray:
    """Bernstein coefficients (e, m) of polynomials given by their lattice values (e, m)."""
    _points, operators = _lattice(kind, degrees)
    if kind in ("triangle", "tetra"):
        return values @ operators[0].T
    if kind in ("quad", "hexa"):
        m = operators[0]
        n1 = m.shape[0]
        d = lagrange.DIM[kind]
        out = values.reshape(len(values), *([n1] * d))
        for axis in range(1, d + 1):
            out = np.moveaxis(np.tensordot(out, m, axes=([axis], [1])), -1, axis)
        return out.reshape(len(values), -1)
    tri, line = operators
    out = values.reshape(len(values), tri.shape[0], line.shape[0])
    out = np.einsum("ij,ejk->eik", tri, out)
    out = np.einsum("kl,eil->eik", line, out)
    return out.reshape(len(values), -1)


def rigorous(kind: str, p: int) -> bool:
    """Whether the Bezier bound is accurate for this kind and order: the lattice-to-Bernstein
    conversion of a tensor element amplifies rounding by cond^dim (6.8e4 per direction at
    degree 17: a 0.1 error for hexahedra of order 6)."""
    if kind == "pyramid":
        return False
    degrees = jacobian_degree(kind, p)
    if kind in ("quad", "hexa"):
        return degrees[0] <= 14
    if kind == "prism":
        return degrees[1] <= 14 and degrees[0] <= 13
    return degrees[0] <= 15


def validity(kind: str, p: int, xyz: np.ndarray, dim: int, sign=None) -> tuple:
    """Rigorous validity of elements (e, n, 3): (status (e,): 1 valid, -1 invalid, 0
    undetermined after degree elevation; lower bound (e,) of J relative to its lattice
    maximum). ``sign``: the orientation of each element (+1 / -1), default from its
    straight-sided centre. Pyramids: sampled at the (p + 2) lattice (status from the samples,
    never 0)."""
    e = len(xyz)
    if e == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    if sign is None:
        sign = np.sign(shape.jacobians_at_centre(kind, p, xyz, dim))
    if kind == "pyramid":
        det = shape.jacobians(kind, p, xyz, dim) * sign[:, None]
        status = np.where(det.min(axis=1) > 0, 1, -1)
        return status, det.min(axis=1) / np.maximum(det.max(axis=1), 1e-300)
    base = jacobian_degree(kind, p)
    if not rigorous(kind, p):  # the conversion is too ill-conditioned: sampled
        det = shape.jacobians(kind, p, xyz, dim) * sign[:, None]
        status = np.where(det.min(axis=1) > 0, 1, -1)
        return status, det.min(axis=1) / np.maximum(det.max(axis=1), 1e-300)
    status = np.zeros(e, dtype=np.int64)
    bound = np.zeros(e)
    todo = np.arange(e)
    for extra in ELEVATIONS:
        degrees = tuple(d + extra for d in base)
        points, _operators = _lattice(kind, degrees)
        det = shape.jacobians(kind, p, xyz[todo], dim, points) * sign[todo, None]
        scale = np.maximum(det.max(axis=1), 1e-300)
        coefficients = bezier_coefficients(kind, degrees, det)
        low = coefficients.min(axis=1) / scale
        bound[todo] = low
        invalid = det.min(axis=1) <= 0.0
        valid = (low > 0.0) & ~invalid
        status[todo[invalid]] = -1
        status[todo[valid]] = 1
        todo = todo[~(invalid | valid)]
        if len(todo) == 0:
            break
    return status, bound


# ------------------------------------------------------------------ repair
def _cofactors(a: np.ndarray) -> np.ndarray:
    """Cofactor matrices of (..., d, d) with d = 2 or 3: d det / d a."""
    if a.shape[-1] == 2:
        out = np.empty_like(a)
        out[..., 0, 0] = a[..., 1, 1]
        out[..., 0, 1] = -a[..., 1, 0]
        out[..., 1, 0] = -a[..., 0, 1]
        out[..., 1, 1] = a[..., 0, 0]
        return out
    out = np.empty_like(a)
    for i in range(3):
        for j in range(3):
            rows = [k for k in range(3) if k != i]
            cols = [k for k in range(3) if k != j]
            minor = (
                a[..., rows[0], cols[0]] * a[..., rows[1], cols[1]]
                - a[..., rows[0], cols[1]] * a[..., rows[1], cols[0]]
            )
            out[..., i, j] = (-1) ** (i + j) * minor
    return out


class Patch:
    """The elements of a repair (flagged ones and their neighbours), their sample-point
    gradients, straight-sided Jacobians and free nodes, with the distortion and its gradient
    in the free nodes' coordinates."""

    def __init__(self, mesh, straight, elements: dict, free: np.ndarray):
        self.dim = mesh.dim
        self.free = np.asarray(free, dtype=np.int64)
        nodes = np.unique(
            np.concatenate([mesh.cells[k][r].reshape(-1) for k, r in elements.items()])
        )
        self.nodes = nodes
        local = {int(n): k for k, n in enumerate(nodes)}
        self.local_free = np.array([local[int(n)] for n in self.free], dtype=np.int64)
        self.positions = mesh.points[nodes].copy()
        self.kinds = []
        for kind, rows in elements.items():
            conn = mesh.cells[kind][rows]
            conn_local = np.vectorize(local.__getitem__)(conn)
            _v, grads = shape.basis(kind, mesh.order, shape.sample_points(kind, mesh.order))
            j0 = shape.jacobians(kind, mesh.order, straight[conn], mesh.dim)
            sign = np.sign(shape.jacobians_at_centre(kind, mesh.order, straight[conn], mesh.dim))
            j0 = j0 * sign[:, None]
            j0 = np.maximum(j0, 1e-12 * np.abs(j0).max(axis=1, keepdims=True))
            self.kinds.append((kind, conn_local, grads, j0, sign))
        self.delta = DELTAS[0]
        self.evaluations = 0

    def ratios(self, positions=None) -> np.ndarray:
        """Scaled Jacobians (min / max over the samples) of the patch's elements."""
        if positions is None:
            positions = self.positions
        out = []
        for _kind, conn, grads, _j0, sign in self.kinds:
            det = self._determinants(positions[conn], grads) * sign[:, None]
            out.append(det.min(axis=1) / np.maximum(det.max(axis=1), 1e-300))
        return np.concatenate(out)

    def _matrices(self, x: np.ndarray, grads: np.ndarray) -> np.ndarray:
        """Jacobian matrices (e, q, d, d) of elements with nodes x (e, n, 3)."""
        q, n, d = grads.shape
        flat = grads.transpose(1, 0, 2).reshape(n, q * d)
        a = (np.ascontiguousarray(x.transpose(0, 2, 1)) @ flat).reshape(len(x), 3, q, d)
        return a.transpose(0, 2, 1, 3)[:, :, : self.dim, : self.dim]

    def _determinants(self, x: np.ndarray, grads: np.ndarray) -> np.ndarray:
        return np.linalg.det(self._matrices(x, grads))

    def energy(self, x_free: np.ndarray) -> tuple[float, np.ndarray]:
        positions = self.positions.copy()
        positions[self.local_free] = x_free
        total = 0.0
        gradient = np.zeros_like(positions)
        d = self.dim
        delta = self.delta
        self.evaluations += 1
        for _kind, conn, grads, j0, sign in self.kinds:
            x = positions[conn]  # (e, n, 3)
            a = self._matrices(x, grads)  # (e, q, d, d)
            det = np.linalg.det(a) * sign[:, None]
            r = det / j0
            root = np.sqrt(r * r + delta * delta)
            h = 0.5 * (r + root)
            total += float(np.sum((r - 1.0) ** 2 + ALPHA / h))
            dh = 0.5 * (1.0 + r / root)
            df = 2.0 * (r - 1.0) - ALPHA * dh / (h * h)
            coefficient = df * sign[:, None] / j0  # dE / d det
            cof = _cofactors(a) * coefficient[:, :, None, None]  # (e, q, d, d)
            q, n, _ = grads.shape
            flat = grads[:, :, :d].transpose(1, 0, 2).reshape(n, q * d)  # (n, q d)
            dx = (cof.transpose(0, 2, 1, 3).reshape(len(x), d, q * d) @ flat.T).transpose(0, 2, 1)
            if d == 2:
                dx = np.concatenate([dx, np.zeros((*dx.shape[:2], 1))], axis=2)
            np.add.at(gradient, conn, dx)
        return total, gradient[self.local_free]

    def minimise(self) -> int:
        """L-BFGS; the regularisation rounds only while something is folded. Returns the
        iteration count."""
        x = self.positions[self.local_free].copy()
        iterations = 0
        self.evaluations = 0
        for delta in DELTAS:
            self.delta = delta
            x, used = _lbfgs(self.energy, x, MAX_ITERATIONS)
            iterations += used
            self.positions[self.local_free] = x
            if self.ratios().min() > 0.0 or self.evaluations > MAX_EVALUATIONS:
                break
        return iterations


def _lbfgs(function, x0, max_iterations: int, memory: int = 10):
    x = x0.copy()
    f, g = function(x)
    s_list, y_list = [], []
    scale = max(float(np.abs(x).max()), 1.0)
    for iteration in range(max_iterations):
        norm = float(np.linalg.norm(g))
        if norm < 1e-12:
            return x, iteration
        # Two-loop recursion.
        q = g.copy()
        alphas = []
        for s, y in reversed(list(zip(s_list, y_list, strict=True))):
            rho = 1.0 / float(np.sum(y * s))
            alpha = rho * float(np.sum(s * q))
            alphas.append((alpha, rho, s, y))
            q = q - alpha * y
        if s_list:
            s, y = s_list[-1], y_list[-1]
            q = q * (float(np.sum(s * y)) / float(np.sum(y * y)))
        else:
            q = q * (0.01 * scale / norm)
        for alpha, rho, s, y in reversed(alphas):
            beta = rho * float(np.sum(y * q))
            q = q + s * (alpha - beta)
        direction = -q
        slope = float(np.sum(g * direction))
        if slope >= 0.0:  # not a descent direction: restart
            s_list, y_list = [], []
            direction = -g * (0.01 * scale / norm)
            slope = float(np.sum(g * direction))
        step = 1.0
        accepted = False
        for _ in range(20):
            x_new = x + step * direction
            f_new, g_new = function(x_new)
            if np.isfinite(f_new) and f_new <= f + 1e-4 * step * slope:
                accepted = True
                break
            step *= 0.5
        if not accepted:
            return x, iteration
        s = x_new - x
        y = g_new - g
        if float(np.sum(s * y)) > 1e-14 * float(np.sum(s * s) + 1e-300):
            s_list.append(s)
            y_list.append(y)
            if len(s_list) > memory:
                s_list.pop(0)
                y_list.pop(0)
        decrease = f - f_new
        x, f, g = x_new, f_new, g_new
        if decrease < 1e-10 * max(abs(f), 1.0):
            return x, iteration + 1
    return x, max_iterations


def flagged_elements(ratio: np.ndarray, ratio_straight: np.ndarray) -> np.ndarray:
    """Mask of the elements to repair: folded, poor, or having lost most of their quality."""
    return (ratio <= 0.0) | (ratio < REPAIR_RATIO) | (ratio < REPAIR_LOSS * ratio_straight)


def repair(mesh, straight, flagged: dict, fixed: np.ndarray) -> dict:
    """Repair the flagged elements ({kind: element rows}) of ``mesh`` in place; ``fixed``
    marks the nodes that must not move (boundary, geometry). Returns the report."""
    flagged_nodes = np.unique(
        np.concatenate([mesh.cells[k][r].reshape(-1) for k, r in flagged.items()])
    )
    free = flagged_nodes[(flagged_nodes >= mesh.vertices) & (fixed[flagged_nodes] == 0)]
    report = {"elements": int(sum(len(r) for r in flagged.values())), "free_nodes": len(free)}
    if len(free) == 0:
        report["iterations"] = 0
        return report
    # The patch: flagged elements and every element sharing a node with a free node.
    touched = np.zeros(len(mesh.points), dtype=bool)
    touched[free] = True
    elements = {}
    for kind, conn in mesh.cells.items():
        rows = np.flatnonzero(touched[conn].any(axis=1))
        if len(rows):
            elements[kind] = rows
    patch = Patch(mesh, straight, elements, free)
    before = patch.ratios()
    iterations = patch.minimise()
    after = patch.ratios()
    report["iterations"] = int(iterations)
    report["min_scaled_jacobian_before"] = float(before.min())
    report["min_scaled_jacobian_after"] = float(after.min())
    report["patch_elements"] = len(before)
    if after.min() > before.min():
        mesh.points[patch.nodes] = patch.positions
        report["kept"] = True
    else:
        report["kept"] = False
    return report
