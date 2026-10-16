"""GPU benchmark for 11d: 3D Winslow sweeps on an ni x nj x nk structured block (OpenCL).

Kernels (one work-item per node, or per k-line):
  point_fp32   - point Jacobi, positions and arithmetic in float32
  point_df     - point Jacobi, positions as float-float (hi + lo), differences formed from
                 both parts, metrics and update in float32, x += omega delta with two-sum
  line_df      - k-line Jacobi in delta form (Thomas along k per work-item), float-float
                 positions, float32 arithmetic (the CPU smoother relaxes lines)
Prints ms per sweep and the change per sweep for a sanity check.
Needs pyopencl (not in the opengrid environment):
    python tools/gpu/opencl_bench.py ni nj nk [sweeps [dump.npy]]
"""

import sys
import time

import numpy as np
import pyopencl as cl

SRC = r"""
#define IDX(i, j, k) ((k) * (NI * NJ) + (j) * NI + (i))

// float-float difference a - b with a = ah + al, b = bh + bl (accurate to float32 rounding
// of the result; ah - bh is exact for neighbouring nodes by Sterbenz).
inline float3 dd(float3 ah, float3 al, float3 bh, float3 bl) { return (ah - bh) + (al - bl); }

inline void cof(float3 xi, float3 xj, float3 xk, float* a) {
  float g11 = dot(xi, xi), g22 = dot(xj, xj), g33 = dot(xk, xk);
  float g12 = dot(xi, xj), g13 = dot(xi, xk), g23 = dot(xj, xk);
  a[0] = g22 * g33 - g23 * g23;  // A11
  a[1] = g11 * g33 - g13 * g13;  // A22
  a[2] = g11 * g22 - g12 * g12;  // A33
  a[3] = g13 * g23 - g12 * g33;  // A12
  a[4] = g12 * g23 - g13 * g22;  // A13
  a[5] = g12 * g13 - g11 * g23;  // A23
}

#define LOAD(dst, i, j, k) dst = vload3(IDX(i, j, k), xh)
#define LOADL(dst, i, j, k) dst = vload3(IDX(i, j, k), xl)

// Differences of the 18 needed neighbours relative to the centre (float-float or float).
#define D(i, j, k) dd(vload3(IDX(i, j, k), xh), vload3(IDX(i, j, k), xl), ch, cl_)
#define DF(i, j, k) (vload3(IDX(i, j, k), xh) - ch)

inline float3 winslow_delta(float3 dip, float3 dim, float3 djp, float3 djm, float3 dkp,
                            float3 dkm, float3 xij, float3 xik, float3 xjk, float* aout) {
  float3 xi = 0.5f * (dip - dim), xj = 0.5f * (djp - djm), xk = 0.5f * (dkp - dkm);
  float a[6];
  cof(xi, xj, xk, a);
  float3 num = a[0] * (dip + dim) + a[1] * (djp + djm) + a[2] * (dkp + dkm) +
               2.0f * (a[3] * xij + a[4] * xik + a[5] * xjk);
  float diag = 2.0f * (a[0] + a[1] + a[2]);
  if (aout) { for (int m = 0; m < 6; ++m) aout[m] = a[m]; }
  return num / diag;
}

__kernel void point_fp32(__global const float* xh, __global float* out, float omega,
                         __global float* change) {
  int i = get_global_id(0) + 1, j = get_global_id(1) + 1, k = get_global_id(2) + 1;
  float3 ch = vload3(IDX(i, j, k), xh);
  float3 xij = 0.25f * (DF(i+1,j+1,k) - DF(i+1,j-1,k) - DF(i-1,j+1,k) + DF(i-1,j-1,k));
  float3 xik = 0.25f * (DF(i+1,j,k+1) - DF(i+1,j,k-1) - DF(i-1,j,k+1) + DF(i-1,j,k-1));
  float3 xjk = 0.25f * (DF(i,j+1,k+1) - DF(i,j+1,k-1) - DF(i,j-1,k+1) + DF(i,j-1,k-1));
  float3 d = winslow_delta(DF(i+1,j,k), DF(i-1,j,k), DF(i,j+1,k), DF(i,j-1,k), DF(i,j,k+1),
                           DF(i,j,k-1), xij, xik, xjk, 0);
  vstore3(ch + omega * d, IDX(i, j, k), out);
  change[IDX(i, j, k)] = dot(d, d);
}

// two-sum: (h, l) + d -> renormalised (h, l)
inline void add_ff(float3* h, float3* l, float3 d) {
  float3 s = *h + d;
  float3 bb = s - *h;
  float3 err = (*h - (s - bb)) + (d - bb);
  float3 lo = *l + err;
  float3 hi = s + lo;
  *l = lo - (hi - s);
  *h = hi;
}

__kernel void point_df(__global const float* xh, __global const float* xl, __global float* oh,
                       __global float* ol, float omega, __global float* change) {
  int i = get_global_id(0) + 1, j = get_global_id(1) + 1, k = get_global_id(2) + 1;
  float3 ch = vload3(IDX(i, j, k), xh), cl_ = vload3(IDX(i, j, k), xl);
  float3 xij = 0.25f * (D(i+1,j+1,k) - D(i+1,j-1,k) - D(i-1,j+1,k) + D(i-1,j-1,k));
  float3 xik = 0.25f * (D(i+1,j,k+1) - D(i+1,j,k-1) - D(i-1,j,k+1) + D(i-1,j,k-1));
  float3 xjk = 0.25f * (D(i,j+1,k+1) - D(i,j+1,k-1) - D(i,j-1,k+1) + D(i,j-1,k-1));
  float3 d = winslow_delta(D(i+1,j,k), D(i-1,j,k), D(i,j+1,k), D(i,j-1,k), D(i,j,k+1),
                           D(i,j,k-1), xij, xik, xjk, 0);
  float3 h = ch, l = cl_;
  add_ff(&h, &l, omega * d);
  vstore3(h, IDX(i, j, k), oh);
  vstore3(l, IDX(i, j, k), ol);
  change[IDX(i, j, k)] = dot(d, d);
}

// k-line Jacobi in delta form: A33 implicit along k, the rest explicit. One work-item per
// (i, j) column; Thomas in private memory (NK <= 256).
__kernel void line_df(__global const float* xh, __global const float* xl, __global float* oh,
                      __global float* ol, float omega, __global float* change) {
  int i = get_global_id(0) + 1, j = get_global_id(1) + 1;
  float cp[NK];
  float3 dp[NK];
  float prev_c = 0.0f;
  float3 prev_d = (float3)(0.0f);
  for (int k = 1; k < NK - 1; ++k) {
    float3 ch = vload3(IDX(i, j, k), xh), cl_ = vload3(IDX(i, j, k), xl);
    float3 xij = 0.25f * (D(i+1,j+1,k) - D(i+1,j-1,k) - D(i-1,j+1,k) + D(i-1,j-1,k));
    float3 xik = 0.25f * (D(i+1,j,k+1) - D(i+1,j,k-1) - D(i-1,j,k+1) + D(i-1,j,k-1));
    float3 xjk = 0.25f * (D(i,j+1,k+1) - D(i,j+1,k-1) - D(i,j-1,k+1) + D(i,j-1,k-1));
    float a[6];
    float3 r = winslow_delta(D(i+1,j,k), D(i-1,j,k), D(i,j+1,k), D(i,j-1,k), D(i,j,k+1),
                             D(i,j,k-1), xij, xik, xjk, a);
    // residual * diag = num; tridiagonal: -A33 d_{k-1} + diag d_k - A33 d_{k+1} = num
    float diag = 2.0f * (a[0] + a[1] + a[2]);
    float3 num = r * diag;
    float lo = -a[2], up = -a[2];
    float m = diag - lo * prev_c;
    prev_c = up / m;
    prev_d = (num - lo * prev_d) / m;
    cp[k] = prev_c;
    dp[k] = prev_d;
  }
  float3 next = (float3)(0.0f);
  float sum = 0.0f;
  for (int k = NK - 2; k >= 1; --k) {
    float3 d = dp[k] - cp[k] * next;
    next = d;
    float3 h = vload3(IDX(i, j, k), xh), l = vload3(IDX(i, j, k), xl);
    add_ff(&h, &l, omega * d);
    vstore3(h, IDX(i, j, k), oh);
    vstore3(l, IDX(i, j, k), ol);
    sum += dot(d, d);
  }
  change[j * NI + i] = sum;
}
"""


def grid(ni, nj, nk):
    # A stretched, sheared box: k clustered towards k = 0 (first cell 1e-5), mild shear.
    i = np.linspace(0, 1, ni)
    j = np.linspace(0, 1, nj)
    q = 1.12
    kk = np.cumsum(np.concatenate([[0.0], 1e-5 * q ** np.arange(nk - 1)]))
    kk /= kk[-1]
    gk, gj, gi = np.meshgrid(kk, j, i, indexing="ij")
    x = np.stack([gi + 0.1 * np.sin(3 * gj) * gk, gj, gk + 0.05 * gi * (1 - gk)], -1)
    return x.astype(np.float64)


def main():
    ni, nj, nk = (int(v) for v in sys.argv[1:4])
    sweeps_arg = 50
    if len(sys.argv) > 4:
        sweeps_arg = int(sys.argv[4])
    dump = None
    if len(sys.argv) > 5:
        dump = sys.argv[5]
    ctx = cl.create_some_context(interactive=False)
    queue = cl.CommandQueue(ctx)
    dev = ctx.devices[0]
    print(f"device {dev.name}, {dev.max_compute_units} compute units, OpenCL {dev.version}")
    prog = cl.Program(ctx, SRC).build(options=f"-DNI={ni} -DNJ={nj} -DNK={nk}")
    x = grid(ni, nj, nk)
    hi = x.astype(np.float32)
    lo = (x - hi.astype(np.float64)).astype(np.float32)
    mf = cl.mem_flags
    n = ni * nj * nk
    bh = [
        cl.Buffer(ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=hi.reshape(-1).copy())
        for _ in range(2)
    ]
    bl = [
        cl.Buffer(ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=lo.reshape(-1).copy())
        for _ in range(2)
    ]
    change = cl.Buffer(ctx, mf.READ_WRITE, size=4 * n)
    interior = (ni - 2, nj - 2, nk - 2)
    sweeps = sweeps_arg
    results = {}
    for name in ("point_fp32", "point_df", "line_df"):
        kernel = getattr(prog, name)
        if name == "line_df":
            size = (ni - 2, nj - 2)
        else:
            size = interior
        for _warm in range(2):
            t0 = time.perf_counter()
            for s in range(sweeps):
                a, b = s % 2, 1 - s % 2
                if name == "point_fp32":
                    kernel(queue, size, None, bh[a], bh[b], np.float32(1.0), change)
                else:
                    kernel(queue, size, None, bh[a], bl[a], bh[b], bl[b], np.float32(1.0), change)
            queue.finish()
            elapsed = time.perf_counter() - t0
        results[name] = 1e3 * elapsed / sweeps
        if dump and name == "line_df":
            # final positions after the second (timed) run: start from the initial grid again
            for b in (0, 1):
                cl.enqueue_copy(queue, bh[b], hi.reshape(-1).copy())
                cl.enqueue_copy(queue, bl[b], lo.reshape(-1).copy())
            for s in range(sweeps):
                a, b = s % 2, 1 - s % 2
                kernel(queue, size, None, bh[a], bl[a], bh[b], bl[b], np.float32(1.0), change)
            h = np.empty(n * 3, np.float32)
            low = np.empty(n * 3, np.float32)
            last = sweeps % 2
            cl.enqueue_copy(queue, h, bh[last])
            cl.enqueue_copy(queue, low, bl[last])
            queue.finish()
            np.save(dump, h.astype(np.float64) + low.astype(np.float64))
        print(
            f"{name:10s}: {results[name]:7.2f} ms per sweep ({n / 1e6:.2f} M nodes, "
            f"{1e9 * elapsed / sweeps / n:.2f} ns per node)"
        )


if __name__ == "__main__":
    main()
