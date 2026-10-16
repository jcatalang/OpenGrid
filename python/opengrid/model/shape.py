"""Lagrange shape functions of order p on the reference elements and the Jacobian validity of
high-order elements (Phase 11g).

Bases (equispaced nodes, canonical reference elements of ``io/lagrange.py``):

* line, quadrilateral, hexahedron: tensor products of 1D Lagrange polynomials;
* triangle, tetrahedron: Silvester's formula, N = prod_m R_(i_m)(lambda_m) with
  R_i(l) = prod_(s < i) (p l - s) / (s + 1) over the node's barycentric lattice indices;
* prism: the triangle's basis times a 1D Lagrange polynomial in the height;
* pyramid: the rational space of the collapsed coordinates s = (x - z/2) / (1 - z),
  t = (y - z/2) / (1 - z): s^a t^b (1 - z)^max(a, b) z^c with 0 <= c <= p, 0 <= a, b <= p - c
  (it contains the linear pyramid's map), through the inverse of its Vandermonde matrix at the
  nodes. Solvers may use another pyramid space; the nodes, and so the geometry at them, are the
  same.

Validity: the Jacobian determinant of the element map at the nodes of the order-(p + 2)
lattice (for pyramids without the apex, where the collapsed map is singular), relative to the
straight-sided element's orientation; the scaled Jacobian is min / max over those points.

References: P. Silvester, "High-order polynomial triangular finite elements for potential
problems", Int. J. Eng. Sci. 7 (1969) 849-861; M. Bergot, G. Cohen, M. Durufle, "Higher-order
finite elements for hybrid meshes using new nodal pyramidal elements", J. Sci. Comput. 42 (2010)
345-381.
"""

from functools import cache

import numpy as np

from opengrid.io import lagrange


def _lattice(kind, p) -> np.ndarray:
    """(n, d) integer lattice indices of the Gmsh-ordered nodes (coordinates x p; pyramids:
    x 2p)."""
    scale = p
    if kind == "pyramid":
        scale = 2 * p
    out = [[int(c * scale) for c in x] for x in lagrange.nodes("gmsh", kind, p)]
    return np.array(out, dtype=np.int64)


def _lagrange_1d(p, t) -> tuple[np.ndarray, np.ndarray]:
    """Values and derivatives of the p + 1 equispaced Lagrange polynomials on [0, 1] at t."""
    t = np.asarray(t, dtype=float)
    nodes = np.arange(p + 1) / p
    values = np.ones((len(t), p + 1))
    derivs = np.zeros((len(t), p + 1))
    for j in range(p + 1):
        for k in range(p + 1):
            if k == j:
                continue
            factor = (t - nodes[k]) / (nodes[j] - nodes[k])
            derivs[:, j] = derivs[:, j] * factor + values[:, j] / (nodes[j] - nodes[k])
            values[:, j] *= factor
    return values, derivs


def _silvester(p, lam) -> tuple[np.ndarray, np.ndarray]:
    """R_i(lambda) and dR_i / dlambda for i = 0 .. p at the given barycentric values."""
    lam = np.asarray(lam, dtype=float)
    values = np.ones((len(lam), p + 1))
    derivs = np.zeros((len(lam), p + 1))
    for i in range(1, p + 1):
        factor = (p * lam - (i - 1)) / i
        derivs[:, i] = derivs[:, i - 1] * factor + values[:, i - 1] * p / i
        values[:, i] = values[:, i - 1] * factor
    return values, derivs


def _simplex(p, points, lattice) -> tuple[np.ndarray, np.ndarray]:
    """Silvester basis: values (q, n) and gradients (q, n, d) for a triangle or tetrahedron."""
    q, d = points.shape
    lam = np.column_stack([1.0 - points.sum(axis=1), points])  # (q, d + 1)
    bary = np.column_stack([p - lattice.sum(axis=1), lattice])  # (n, d + 1)
    r = []
    dr = []
    for m in range(d + 1):
        v, dv = _silvester(p, lam[:, m])
        r.append(v[:, bary[:, m]])  # (q, n)
        dr.append(dv[:, bary[:, m]])
    values = np.prod(np.stack(r), axis=0)
    grads = np.zeros((q, len(lattice), d))
    for m in range(d + 1):
        others = np.prod(np.stack([r[k] for k in range(d + 1) if k != m]), axis=0)
        term = dr[m] * others  # d N / d lambda_m
        if m == 0:
            grads -= term[:, :, None]  # lambda_0 = 1 - sum x
        else:
            grads[:, :, m - 1] += term
    return values, grads


def _tensor(p, points, lattice) -> tuple[np.ndarray, np.ndarray]:
    q, d = points.shape
    one = [_lagrange_1d(p, points[:, a]) for a in range(d)]
    factors = [one[a][0][:, lattice[:, a]] for a in range(d)]
    dfactors = [one[a][1][:, lattice[:, a]] for a in range(d)]
    values = np.prod(np.stack(factors), axis=0)
    grads = np.zeros((q, len(lattice), d))
    for a in range(d):
        rest = np.prod(
            np.stack([factors[b] for b in range(d) if b != a] or [np.ones_like(values)]), axis=0
        )
        grads[:, :, a] = dfactors[a] * rest
    return values, grads


def _prism(p, points, lattice) -> tuple[np.ndarray, np.ndarray]:
    tv, tg = _simplex(p, points[:, :2], lattice[:, :2])
    zv, zd = _lagrange_1d(p, points[:, 2])
    zv = zv[:, lattice[:, 2]]
    zd = zd[:, lattice[:, 2]]
    values = tv * zv
    grads = np.concatenate([tg * zv[:, :, None], (tv * zd)[:, :, None]], axis=2)
    return values, grads


def _pyramid_monomials(p) -> list:
    return [(a, b, c) for c in range(p + 1) for b in range(p - c + 1) for a in range(p - c + 1)]


def _pyramid_basis(p, points) -> tuple[np.ndarray, np.ndarray]:
    """phi = s^a t^b w^max(a, b) z^c (s = (x - z/2) / w, t = (y - z/2) / w, w = 1 - z) and its
    gradients in (x, y, z); at the apex (w = 0) phi is its limit (0, or z^c for a = b = 0)."""
    x, y, z = points.T
    w = 1.0 - z
    apex = w < 1e-14
    safe = np.where(apex, 1.0, w)
    s = np.where(apex, 0.0, (x - z / 2) / safe)
    t = np.where(apex, 0.0, (y - z / 2) / safe)
    ds_dz = (s - 0.5) / safe
    dt_dz = (t - 0.5) / safe
    terms = _pyramid_monomials(p)
    values = np.zeros((len(points), len(terms)))
    grads = np.zeros((len(points), len(terms), 3))
    for k, (a, b, c) in enumerate(terms):
        m = max(a, b)
        sa = s**a
        tb = t**b
        wm = w**m
        zc = z**c
        values[:, k] = sa * tb * wm * zc
        d_s = a * s ** max(a - 1, 0) * tb * wm * zc  # d phi / d s
        d_t = b * sa * t ** max(b - 1, 0) * wm * zc
        d_w = m * sa * tb * w ** max(m - 1, 0) * zc
        d_z = c * sa * tb * wm * z ** max(c - 1, 0)
        grads[:, k, 0] = d_s / safe
        grads[:, k, 1] = d_t / safe
        grads[:, k, 2] = d_s * ds_dz + d_t * dt_dz - d_w + d_z
    return values, grads


@cache
def _pyramid_coefficients(p) -> np.ndarray:
    nodes = lagrange.as_array(lagrange.nodes("gmsh", "pyramid", p))
    vandermonde, _ = _pyramid_basis(p, nodes)
    return np.linalg.inv(vandermonde)


def basis(kind: str, p: int, points) -> tuple[np.ndarray, np.ndarray]:
    """Shape functions of the order-p ``kind`` element (nodes in Gmsh's order) at reference
    points (q, d): values (q, n) and gradients (q, n, d)."""
    points = np.atleast_2d(np.asarray(points, dtype=float))
    if kind in ("line", "quad", "hexa"):
        return _tensor(p, points, _lattice(kind, p))
    if kind in ("triangle", "tetra"):
        return _simplex(p, points, _lattice(kind, p))
    if kind == "prism":
        return _prism(p, points, _lattice(kind, p))
    values, grads = _pyramid_basis(p, points)
    c = _pyramid_coefficients(p)
    return values @ c, np.einsum("qkd,kn->qnd", grads, c)


@cache
def sample_points(kind: str, p: int) -> np.ndarray:
    """Points where validity is checked: the order-(p + 2) lattice (pyramids: without the
    apex)."""
    x = lagrange.as_array(lagrange.nodes("gmsh", kind, p + 2))
    if kind == "pyramid":
        x = x[x[:, 2] < 1.0 - 1e-12]
    return x


def jacobians(kind: str, p: int, xyz: np.ndarray, dim: int) -> np.ndarray:
    """Jacobian determinants (elements, samples) of elements with node coordinates ``xyz``
    (elements, n, 3); 2D: the in-plane (x, y) determinant."""
    _values, grads = basis(kind, p, sample_points(kind, p))
    q, n, d = grads.shape
    flat = grads.transpose(1, 0, 2).reshape(n, q * d)  # (n, q d)
    out = []
    for start in range(0, len(xyz), 2048):
        chunk = np.ascontiguousarray(xyz[start : start + 2048].transpose(0, 2, 1))  # (e, 3, n)
        jac = (chunk @ flat).reshape(len(chunk), 3, q, d).transpose(0, 2, 1, 3)  # (e, q, 3, d)
        if dim == 2:
            out.append(np.linalg.det(jac[:, :, :2, :2]))
        else:
            out.append(np.linalg.det(jac))
    if not out:
        return np.zeros((0, len(sample_points(kind, p))))
    return np.concatenate(out)


def jacobians_at_centre(kind: str, p: int, xyz: np.ndarray, dim: int) -> np.ndarray:
    """Jacobian determinant (elements,) at the reference element's vertex centroid."""
    centre = lagrange.as_array(lagrange.VERTICES[kind]).mean(axis=0)
    _values, grads = basis(kind, p, centre[None, :])
    jac = np.einsum("qnd,enc->eqcd", grads, xyz)[:, 0]
    if dim == 2:
        return np.linalg.det(jac[:, :2, :2])
    return np.linalg.det(jac)
