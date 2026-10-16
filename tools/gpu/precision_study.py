"""Precision study for a GPU elliptic smoother (11d).

O-grid around a unit circle, first cell 1e-6, outer radius 20, 81 radial x 128 angular nodes.
Inverse control functions make the exact polar grid the discrete solution; the start shears
the interior angularly. Radial line Jacobi in delta form (tridiagonal along each radial line,
angular neighbours explicit), as the CPU smoother does. Variants:
  fp64   - everything in float64
  fp32   - everything in float32 (positions too)
  delta  - positions stored in float64 (stands for two float32 per coordinate on the GPU);
           neighbour differences, metrics, coefficients, residual, line solve in float32;
           only x += omega * delta in float64
Reports the error against the exact grid in units of the local radial cell height.
"""

import time

import numpy as np

NI, NJ = 128, 81
R0, R1, H0 = 1.0, 20.0, 1e-6


def radial():
    # geometric distribution: H0 * (q^n - 1) / (q - 1) = R1 - R0
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        if H0 * (q ** (NJ - 1) - 1) / (q - 1) > R1 - R0:
            hi = q
        else:
            lo = q
    r = R0 + H0 * (q ** np.arange(NJ) - 1) / (q - 1)
    r[-1] = R1
    return r


def exact():
    r = radial()
    t = 2 * np.pi * np.arange(NI) / NI
    x = np.stack([np.outer(np.cos(t), r), np.outer(np.sin(t), r)], -1)  # (NI, NJ, 2)
    return x


def diffs(x, dt):
    """Neighbour differences relative to each interior node (periodic in i), in dtype dt."""
    c = x[:, 1:-1]
    e = np.roll(x, -1, 0)[:, 1:-1]
    w = np.roll(x, 1, 0)[:, 1:-1]
    n = x[:, 2:]
    s = x[:, :-2]
    ne = np.roll(x, -1, 0)[:, 2:]
    nw = np.roll(x, 1, 0)[:, 2:]
    se = np.roll(x, -1, 0)[:, :-2]
    sw = np.roll(x, 1, 0)[:, :-2]
    return [(v - c).astype(dt) for v in (e, w, n, s, ne, nw, se, sw)]


def operator(d, pc, qc, dt):
    """Residual of alpha (x_ss + pc x_s) - 2 beta x_st + gamma (x_tt + qc x_t) from differences."""
    e, w, n, s, ne, nw, se, sw = d
    half = dt(0.5)
    xs = half * (e - w)
    xt = half * (n - s)
    xss = e + w
    xtt = n + s
    xst = dt(0.25) * (ne - nw - se + sw)
    alpha = np.sum(xt * xt, -1, keepdims=True)
    beta = np.sum(xs * xt, -1, keepdims=True)
    gamma = np.sum(xs * xs, -1, keepdims=True)
    res = alpha * (xss + pc * xs) - dt(2) * beta * xst + gamma * (xtt + qc * xt)
    return res, alpha, gamma


def controls(x):
    """pc, qc (per node) with the exact grid as the solution: alpha pc x_s + gamma qc x_t = -L."""
    d = diffs(x, np.float64)
    zero = np.zeros(x[:, 1:-1, :1].shape)
    lhs, alpha, gamma = operator(d, zero, zero, np.float64)
    e, w, n, s = d[:4]
    xs = 0.5 * (e - w)
    xt = 0.5 * (n - s)
    a = np.stack([alpha * xs, gamma * xt], -1)  # (NI, NJ-2, 2, 2): columns pc, qc
    pq = np.linalg.solve(a, -lhs[..., None])[..., 0]
    return pq[..., :1], pq[..., 1:]


def thomas(lo, di, up, rhs):
    """Batched tridiagonal solve along axis 1 (interior unknowns, zero Dirichlet deltas)."""
    n = di.shape[1]
    c = np.zeros_like(up)
    d = np.zeros_like(rhs)
    c[:, 0] = up[:, 0] / di[:, 0]
    d[:, 0] = rhs[:, 0] / di[:, 0][..., None]
    for k in range(1, n):
        m = di[:, k] - lo[:, k] * c[:, k - 1]
        c[:, k] = up[:, k] / m
        d[:, k] = (rhs[:, k] - lo[:, k][..., None] * d[:, k - 1]) / m[..., None]
    out = np.zeros_like(rhs)
    out[:, -1] = d[:, -1]
    for k in range(n - 2, -1, -1):
        out[:, k] = d[:, k] - c[:, k][..., None] * out[:, k + 1]
    return out


def solve(variant, x_exact, iterations=6000, omega=1.0):
    pc, qc = controls(x_exact)
    dt = np.float64
    if variant != "fp64":
        dt = np.float32
    pc = pc.astype(dt)
    qc = qc.astype(dt)
    # Start: interior sheared in angle by up to 0.3 rad (zero at both walls).
    r = radial()
    s = (r - R0) / (R1 - R0)
    t = 2 * np.pi * np.arange(NI) / NI
    shear = 0.3 * np.sin(np.pi * np.sqrt(s))[None, :] * np.cos(2 * t)[:, None]
    tt = t[:, None] + shear
    x = np.stack([np.cos(tt) * r[None, :], np.sin(tt) * r[None, :]], -1)
    store = np.float64
    if variant == "fp32":
        store = np.float32
    x = x.astype(store)
    first = None
    history = []
    for _ in range(iterations):
        d = diffs(x, dt)
        res, alpha, gamma = operator(d, pc, qc, dt)
        # Delta form, along j: -gamma (1 - qc / 2) d_{j-1} + 2 (alpha + gamma) d_j
        #   - gamma (1 + qc / 2) d_{j+1} = res
        g = gamma[..., 0]
        q = qc[..., 0]
        lo = -g * (dt(1) - dt(0.5) * q)
        up = -g * (dt(1) + dt(0.5) * q)
        di = dt(2) * (alpha[..., 0] + g)
        delta = thomas(lo, di, up, res)
        x[:, 1:-1] += (dt(omega) * delta).astype(store)
        rms = float(np.sqrt(np.mean(np.sum(delta.astype(np.float64) ** 2, -1))))
        if first is None:
            first = rms
        history.append(rms)
        if rms < 1e-13 * first:
            break
    return x.astype(np.float64), history


def report(variant, x, x_exact, history, elapsed):
    r = radial()
    h = np.diff(r)  # radial cell heights
    hloc = np.minimum(np.concatenate([[h[0]], h]), np.concatenate([h, [h[-1]]]))
    err = np.linalg.norm(x - x_exact, axis=-1) / hloc[None, :]
    first_cell = np.linalg.norm(x[:, 1] - x[:, 0], axis=-1)
    print(
        f"{variant:5s}: {len(history)} iterations ({elapsed:.1f} s), residual "
        f"{history[-1] / history[0]:.1e} of the first; error / local cell height: max "
        f"{err.max():.1e}, wall row {err[:, 1].max():.1e}, outer rows {err[:, -10:].max():.1e}; "
        f"first cell {first_cell.min():.6e} .. {first_cell.max():.6e} (exact {H0:g})"
    )


if __name__ == "__main__":
    xe = exact()
    for variant in ("fp64", "delta", "fp32"):
        t0 = time.perf_counter()
        x, hist = solve(variant, xe)
        report(variant, x, xe, hist, time.perf_counter() - t0)
