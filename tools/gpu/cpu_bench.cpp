// CPU baseline for 11d: the same 3D Winslow sweeps as opencl_bench.py, float64, OpenMP.
// Build: clang++ -std=c++20 -O3 -ffp-contract=off -Xpreprocessor -fopenmp -I$CONDA_PREFIX/include
//   -L$CONDA_PREFIX/lib -lomp -Wl,-rpath,$CONDA_PREFIX/lib cpu_bench.cpp -o cpu_bench
// Run: ./cpu_bench ni nj nk [sweeps [dump.bin]]
#include <omp.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct V {
  double x, y, z;
};
static inline V operator+(V a, V b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
static inline V operator-(V a, V b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
static inline V operator*(double s, V a) {
  return {s * a.x, s * a.y, s * a.z};
}
static inline double dot(V a, V b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

int NI, NJ, NK;
static inline long idx(int i, int j, int k) {
  return (long(k) * NJ + j) * NI + i;
}

static inline V delta(const V* x, int i, int j, int k, double* a6) {
  const V c = x[idx(i, j, k)];
  auto d = [&](int a, int b, int e) { return x[idx(a, b, e)] - c; };
  V dip = d(i + 1, j, k), dim = d(i - 1, j, k), djp = d(i, j + 1, k), djm = d(i, j - 1, k);
  V dkp = d(i, j, k + 1), dkm = d(i, j, k - 1);
  V xij =
      0.25 * (d(i + 1, j + 1, k) - d(i + 1, j - 1, k) - d(i - 1, j + 1, k) + d(i - 1, j - 1, k));
  V xik =
      0.25 * (d(i + 1, j, k + 1) - d(i + 1, j, k - 1) - d(i - 1, j, k + 1) + d(i - 1, j, k - 1));
  V xjk =
      0.25 * (d(i, j + 1, k + 1) - d(i, j + 1, k - 1) - d(i, j - 1, k + 1) + d(i, j - 1, k - 1));
  V xi = 0.5 * (dip - dim), xj = 0.5 * (djp - djm), xk = 0.5 * (dkp - dkm);
  double g11 = dot(xi, xi), g22 = dot(xj, xj), g33 = dot(xk, xk);
  double g12 = dot(xi, xj), g13 = dot(xi, xk), g23 = dot(xj, xk);
  double a[6] = {g22 * g33 - g23 * g23, g11 * g33 - g13 * g13, g11 * g22 - g12 * g12,
                 g13 * g23 - g12 * g33, g12 * g23 - g13 * g22, g12 * g13 - g11 * g23};
  V num = a[0] * (dip + dim) + a[1] * (djp + djm) + a[2] * (dkp + dkm) +
          2.0 * (a[3] * xij + a[4] * xik + a[5] * xjk);
  double diag = 2.0 * (a[0] + a[1] + a[2]);
  if (a6) {
    for (int m = 0; m < 6; ++m)
      a6[m] = a[m];
  }
  return (1.0 / diag) * num;
}

int main(int argc, char** argv) {
  NI = std::atoi(argv[1]);
  NJ = std::atoi(argv[2]);
  NK = std::atoi(argv[3]);
  const long n = long(NI) * NJ * NK;
  std::vector<V> x(n), y(n);
  for (int k = 0; k < NK; ++k) {
    for (int j = 0; j < NJ; ++j) {
      for (int i = 0; i < NI; ++i) {
        double s = double(k) / (NK - 1);
        double kk = (std::pow(1.12, s * (NK - 1)) - 1) / (std::pow(1.12, NK - 1) - 1);
        double I = double(i) / (NI - 1), J = double(j) / (NJ - 1);
        x[idx(i, j, k)] = {I + 0.1 * std::sin(3 * J) * kk, J, kk + 0.05 * I * (1 - kk)};
      }
    }
  }
  y = x;
  int sweeps = 20;
  if (argc > 4)
    sweeps = std::atoi(argv[4]);
  std::printf("threads %d, %.2f M nodes\n", omp_get_max_threads(), n / 1e6);
  const bool dumping = argc > 5;
  const std::vector<V> initial = x;
  for (int warm = 0; warm < 2 && !dumping; ++warm) {
    auto t0 = std::chrono::steady_clock::now();
    for (int s = 0; s < sweeps; ++s) {
#pragma omp parallel for collapse(2) schedule(static)
      for (int k = 1; k < NK - 1; ++k) {
        for (int j = 1; j < NJ - 1; ++j) {
          for (int i = 1; i < NI - 1; ++i) {
            y[idx(i, j, k)] = x[idx(i, j, k)] + delta(x.data(), i, j, k, nullptr);
          }
        }
      }
      std::swap(x, y);
    }
    double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (warm == 1)
      std::printf("point_fp64: %7.2f ms per sweep (%.2f ns per node)\n", 1e3 * dt / sweeps,
                  1e9 * dt / sweeps / n);
  }
  x = initial;
  y = initial;
  for (int warm = 0; warm < 2; ++warm) {
    if (dumping && warm == 0)
      continue;
    auto t0 = std::chrono::steady_clock::now();
    for (int s = 0; s < sweeps; ++s) {
#pragma omp parallel
      {
        std::vector<double> cp(long(NK) * NI);
        std::vector<V> dp(long(NK) * NI);
#pragma omp for schedule(static)
        for (int j = 1; j < NJ - 1; ++j) {
          // Thomas along k for all i of row j at once (i innermost: contiguous).
          for (int k = 1; k < NK - 1; ++k) {
            for (int i = 1; i < NI - 1; ++i) {
              double a[6];
              V r = delta(x.data(), i, j, k, a);
              double diag = 2.0 * (a[0] + a[1] + a[2]);
              V num = diag * r;
              double lo = -a[2], up = -a[2];
              double pc = 0.0;
              V pd{0, 0, 0};
              if (k > 1) {
                pc = cp[long(k - 1) * NI + i];
                pd = dp[long(k - 1) * NI + i];
              }
              double m = diag - lo * pc;
              cp[long(k) * NI + i] = up / m;
              dp[long(k) * NI + i] = (1.0 / m) * (num - lo * pd);
            }
          }
          for (int k = NK - 2; k >= 1; --k) {
            for (int i = 1; i < NI - 1; ++i) {
              V next{0, 0, 0};
              if (k < NK - 2)
                next = dp[long(k + 1) * NI + i];
              V d = dp[long(k) * NI + i] - cp[long(k) * NI + i] * next;
              dp[long(k) * NI + i] = d;
              y[idx(i, j, k)] = x[idx(i, j, k)] + d;
            }
          }
        }
      }
      std::swap(x, y);
    }
    double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (warm == 1)
      std::printf("line_fp64 : %7.2f ms per sweep (%.2f ns per node)\n", 1e3 * dt / sweeps,
                  1e9 * dt / sweeps / n);
  }
  if (argc > 5) {
    FILE* f = std::fopen(argv[5], "wb");
    std::fwrite(x.data(), sizeof(V), n, f);
    std::fclose(f);
  }
  return 0;
}
