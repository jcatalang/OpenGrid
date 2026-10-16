// OpenCL C sources of the GPU kernels (compiled once by the runtime).
#pragma once

namespace og::gpu {

inline constexpr const char* kKernels = R"CLC(
// No contraction into fused multiply-adds: the two-sum below needs plain float32 rounding, and
// the results stay the same across compilers.
#pragma OPENCL FP_CONTRACT OFF

// a - b for float-float positions (x = hi + lo): exact hi difference for neighbouring nodes
// (Sterbenz), rounded once to float32.
inline float3 ffdiff(__global const float* xh, __global const float* xl, int a, int b) {
  return (vload3(a, xh) - vload3(b, xh)) + (vload3(a, xl) - vload3(b, xl));
}

// x[g] += d in float-float (two-sum, then renormalisation).
inline void ffadd(__global float* xh, __global float* xl, int g, float3 d) {
  float3 h = vload3(g, xh);
  float3 l = vload3(g, xl);
  float3 s = h + d;
  float3 bb = s - h;
  float3 err = (h - (s - bb)) + (d - bb);
  float3 lo = l + err;
  float3 hi = s + lo;
  vstore3(hi, g, xh);
  vstore3(lo - (hi - s), g, xl);
}

#define SLOT(i, j, k) (13 + (i) + 3 * (j) + 9 * (k))

// Self-test: out[q] = (hi, lo) of a[q] + d[q] and the float-float difference a[q] - a[q + 1].
__kernel void fftest(__global float* xh, __global float* xl, __global const float* d,
                     __global float* diff) {
  int q = get_global_id(0);
  vstore3(ffdiff(xh, xl, 2 * q, 2 * q + 1), q, diff);
  ffadd(xh, xl, 2 * q, vload3(q, d));
}

// Frozen coefficients of the regular nodes (as Coefficients in multiblock_elliptic3d.cpp):
// coef[16 r] = a[6] by direction (2 m: +m, 2 m + 1: -m), diag, cross term (3).
__kernel void coefficients(int count, __global const int* regular, __global const int* stencil,
                           __global const float* src, __global const float* xh,
                           __global const float* xl, __global float* coef) {
  int r = get_global_id(0);
  if (r >= count) {
    return;
  }
  int g = regular[r];
  __global const int* s = stencil + 27 * r;
#define D(i, j, k) ffdiff(xh, xl, s[SLOT(i, j, k)], g)
  float3 d0 = 0.5f * (D(1, 0, 0) - D(-1, 0, 0));
  float3 d1 = 0.5f * (D(0, 1, 0) - D(0, -1, 0));
  float3 d2 = 0.5f * (D(0, 0, 1) - D(0, 0, -1));
  float g11 = dot(d0, d0);
  float g22 = dot(d1, d1);
  float g33 = dot(d2, d2);
  float g12 = dot(d0, d1);
  float g13 = dot(d0, d2);
  float g23 = dot(d1, d2);
  float a0 = g22 * g33 - g23 * g23;
  float a1 = g11 * g33 - g13 * g13;
  float a2 = g11 * g22 - g12 * g12;
  float s0 = src[3 * r];
  float s1 = src[3 * r + 1];
  float s2 = src[3 * r + 2];
  __global float* c = coef + 16 * r;
  c[0] = a0 * (1.0f + 0.5f * s0);
  c[1] = a0 * (1.0f - 0.5f * s0);
  c[2] = a1 * (1.0f + 0.5f * s1);
  c[3] = a1 * (1.0f - 0.5f * s1);
  c[4] = a2 * (1.0f + 0.5f * s2);
  c[5] = a2 * (1.0f - 0.5f * s2);
  c[6] = 2.0f * a0 + 2.0f * a1 + 2.0f * a2;
  float a12 = g13 * g23 - g12 * g33;
  float a13 = g12 * g23 - g13 * g22;
  float a23 = g12 * g13 - g11 * g23;
  float3 c01 = 0.25f * (D(1, 1, 0) - D(1, -1, 0) - D(-1, 1, 0) + D(-1, -1, 0));
  float3 c02 = 0.25f * (D(1, 0, 1) - D(1, 0, -1) - D(-1, 0, 1) + D(-1, 0, -1));
  float3 c12 = 0.25f * (D(0, 1, 1) - D(0, 1, -1) - D(0, -1, 1) + D(0, -1, -1));
  float3 cross = 2.0f * (a12 * c01 + a13 * c02 + a23 * c12);
  c[7] = cross.x;
  c[8] = cross.y;
  c[9] = cross.z;
#undef D
}

// Line relaxation of one colour in delta form: per line, the tridiagonal system along the
// line (the line's two directions implicit, the others explicit, end nodes fixed) for the
// update of its nodes; residual = cross + sum over the six directions a[d] (x_d - x_g).
// Lines of one colour share no node; they read the positions at the start of the colour
// (xh, xl) and write their nodes into (wh, wl), a copy of them.
__kernel void relax_lines(int count, int first, __global const int* lines,
                          __global const int* line_offset, __global const int* line_size,
                          __global const int* line_nodes, __global const char* line_prev,
                          __global const char* line_next, __global const int* coef_of,
                          __global const int* stencil, __global const float* coef,
                          __global const float* xh, __global const float* xl,
                          __global float* wh, __global float* wl, __global float* scratch_c,
                          __global float* scratch_d, float omega, __global float* line_sum) {
  int q = get_global_id(0);
  if (q >= count) {
    return;
  }
  int l = lines[first + q];
  int off = line_offset[l];
  int m = line_size[l];
  float pc = 0.0f;
  float3 pd = (float3)(0.0f);
  for (int k = 0; k < m; ++k) {
    int g = line_nodes[off + k];
    int r = coef_of[g];
    __global const float* c = coef + 16 * r;
    __global const int* s = stencil + 27 * r;
    float3 res = (float3)(c[7], c[8], c[9]);
    res += c[0] * ffdiff(xh, xl, s[SLOT(1, 0, 0)], g);
    res += c[1] * ffdiff(xh, xl, s[SLOT(-1, 0, 0)], g);
    res += c[2] * ffdiff(xh, xl, s[SLOT(0, 1, 0)], g);
    res += c[3] * ffdiff(xh, xl, s[SLOT(0, -1, 0)], g);
    res += c[4] * ffdiff(xh, xl, s[SLOT(0, 0, 1)], g);
    res += c[5] * ffdiff(xh, xl, s[SLOT(0, 0, -1)], g);
    float lo = -c[(int)line_prev[off + k]];
    float up = -c[(int)line_next[off + k]];
    float mm = c[6] - lo * pc;
    pc = up / mm;
    pd = (res - lo * pd) / mm;
    scratch_c[off + k] = pc;
    vstore3(pd, off + k, scratch_d);
  }
  float3 next = (float3)(0.0f);
  float sum = 0.0f;
  for (int k = m - 1; k >= 0; --k) {
    float3 dlt = vload3(off + k, scratch_d) - scratch_c[off + k] * next;
    next = dlt;
    ffadd(wh, wl, line_nodes[off + k], omega * dlt);
    sum += dot(dlt, dlt);
  }
  line_sum[l] = sum;
}
)CLC";

}  // namespace og::gpu
