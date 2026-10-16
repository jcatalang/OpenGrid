#include "opengrid/core/highorder/elements.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

#include "opengrid/core/util/parallel.hpp"

namespace og {

namespace {

// out[a, b, c] = sum A[a, i] B[b, j] C[c, k] in[i, j, k] (sum factorisation; i fastest in both,
// in: n^3 values, A: na x n, ...). For 2D pass nc = 1 and C = {1} with n3 = 1.
void contract(const std::vector<double>& a, std::size_t na, const std::vector<double>& b,
              std::size_t nb, const std::vector<double>& c, std::size_t nc, const Vec3* in,
              std::size_t n, std::size_t n3, Vec3* out, std::vector<Vec3>& t1,
              std::vector<Vec3>& t2) {
  // t1[a, j, k] = sum_i A[a, i] in[i, j, k]
  t1.assign(na * n * n3, Vec3{});
  for (std::size_t k = 0; k < n3; ++k) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t q = 0; q < na; ++q) {
        Vec3 s;
        for (std::size_t i = 0; i < n; ++i) {
          s += a[q * n + i] * in[(k * n + j) * n + i];
        }
        t1[(k * n + j) * na + q] = s;
      }
    }
  }
  // t2[a, b, k] = sum_j B[b, j] t1[a, j, k]
  t2.assign(na * nb * n3, Vec3{});
  for (std::size_t k = 0; k < n3; ++k) {
    for (std::size_t q = 0; q < nb; ++q) {
      for (std::size_t p = 0; p < na; ++p) {
        Vec3 s;
        for (std::size_t j = 0; j < n; ++j) {
          s += b[q * n + j] * t1[(k * n + j) * na + p];
        }
        t2[(k * nb + q) * na + p] = s;
      }
    }
  }
  // out[a, b, c] = sum_k C[c, k] t2[a, b, k]
  for (std::size_t r = 0; r < nc; ++r) {
    for (std::size_t q = 0; q < nb; ++q) {
      for (std::size_t p = 0; p < na; ++p) {
        Vec3 s;
        for (std::size_t k = 0; k < n3; ++k) {
          s += c[r * n3 + k] * t2[(k * nb + q) * na + p];
        }
        out[(r * nb + q) * na + p] = s;
      }
    }
  }
}

}  // namespace

HighOrderElements high_order_from_grid(const std::vector<Vec3>& grid, int ni, int nj, int nk,
                                       int order, NodeSet set) {
  if (order < 1) {
    throw std::invalid_argument("high-order elements: the order must be at least 1");
  }
  if (ni < 2 || nj < 2 || nk < 1 ||
      grid.size() != static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj) *
                         static_cast<std::size_t>(nk)) {
    throw std::invalid_argument("high-order elements: the grid size does not match");
  }
  HighOrderElements out;
  out.dim = 2;
  if (nk > 1) {
    out.dim = 3;
  }
  out.order = order;
  const std::array<int, 3> counts = {ni, nj, nk};
  for (int d = 0; d < out.dim; ++d) {
    const int intervals = counts[static_cast<std::size_t>(d)] - 1;
    if (intervals % order != 0) {
      throw std::invalid_argument("high-order elements: " + std::to_string(intervals) +
                                  " intervals along " + std::string(1, "ijk"[d]) +
                                  " are not a multiple of the order " + std::to_string(order));
    }
    out.ne[static_cast<std::size_t>(d)] = intervals / order;
  }
  if (out.dim == 2) {
    out.ne[2] = 1;
  }
  out.count = out.ne[0] * out.ne[1] * out.ne[2];
  const auto n = static_cast<std::size_t>(order) + 1;
  std::size_t n3 = 1;
  if (out.dim == 3) {
    n3 = n;
  }
  const std::size_t per = n * n * n3;
  const std::vector<double> eq = nodes_1d(NodeSet::Equispaced, static_cast<int>(n));
  const std::vector<double> target = nodes_1d(set, static_cast<int>(n));
  const std::vector<double> to_set = interpolation_matrix(eq, target);
  const std::vector<double> check = nodes_1d(NodeSet::GaussLobatto, static_cast<int>(n) + 1);
  const std::size_t nq = check.size();
  std::size_t nq3 = 1;
  if (out.dim == 3) {
    nq3 = nq;
  }
  const std::vector<double> m = interpolation_matrix(eq, check);
  const std::vector<double> dm = differentiation_matrix(eq, check);
  const std::vector<double> one = {1.0};
  out.xyz.assign(static_cast<std::size_t>(out.count) * per * 3, 0.0);
  out.min_jacobian.assign(static_cast<std::size_t>(out.count), 0.0);
  out.max_jacobian.assign(static_cast<std::size_t>(out.count), 0.0);
  const auto sni = static_cast<std::size_t>(ni);
  const auto snj = static_cast<std::size_t>(nj);
  const auto p = static_cast<std::size_t>(order);
  parallel_for(
      0, out.count,
      [&](std::int64_t e) {
        const auto ue = static_cast<std::size_t>(e);
        const std::size_t ei = ue % static_cast<std::size_t>(out.ne[0]);
        const std::size_t ej =
            (ue / static_cast<std::size_t>(out.ne[0])) % static_cast<std::size_t>(out.ne[1]);
        const std::size_t ek = ue / static_cast<std::size_t>(out.ne[0] * out.ne[1]);
        std::vector<Vec3> local(per);
        for (std::size_t c = 0; c < n3; ++c) {
          for (std::size_t b = 0; b < n; ++b) {
            for (std::size_t a = 0; a < n; ++a) {
              const std::size_t g = ((ek * p + c) * snj + (ej * p + b)) * sni + (ei * p + a);
              local[(c * n + b) * n + a] = grid[g];
            }
          }
        }
        std::vector<Vec3> t1;
        std::vector<Vec3> t2;
        std::vector<Vec3> nodes(per);
        if (set == NodeSet::Equispaced) {
          nodes = local;
        } else {
          const std::vector<double>* third = &one;
          if (out.dim == 3) {
            third = &to_set;
          }
          contract(to_set, n, to_set, n, *third, n3, local.data(), n, n3, nodes.data(), t1, t2);
        }
        double* dst = &out.xyz[ue * per * 3];
        for (std::size_t q = 0; q < per; ++q) {
          dst[3 * q] = nodes[q].x;
          dst[3 * q + 1] = nodes[q].y;
          dst[3 * q + 2] = nodes[q].z;
        }
        // Jacobian at the check points.
        const std::size_t count = nq * nq * nq3;
        std::vector<Vec3> du(count);
        std::vector<Vec3> dv(count);
        std::vector<Vec3> dw(count);
        if (out.dim == 3) {
          contract(dm, nq, m, nq, m, nq, local.data(), n, n, du.data(), t1, t2);
          contract(m, nq, dm, nq, m, nq, local.data(), n, n, dv.data(), t1, t2);
          contract(m, nq, m, nq, dm, nq, local.data(), n, n, dw.data(), t1, t2);
        } else {
          contract(dm, nq, m, nq, one, 1, local.data(), n, 1, du.data(), t1, t2);
          contract(m, nq, dm, nq, one, 1, local.data(), n, 1, dv.data(), t1, t2);
        }
        std::vector<double> det(count);
        Vec3 normal;
        for (std::size_t q = 0; q < count; ++q) {
          normal += cross(du[q], dv[q]);
        }
        const double nn = norm(normal);
        double sum = 0.0;
        for (std::size_t q = 0; q < count; ++q) {
          if (out.dim == 3) {
            det[q] = dot(cross(du[q], dv[q]), dw[q]);
          } else if (nn > 0.0) {
            det[q] = dot(cross(du[q], dv[q]), (1.0 / nn) * normal);
          } else {
            det[q] = 0.0;
          }
          sum += det[q];
        }
        // Orientation of the grid (a left-handed block): the sign of the mean.
        double sign = 1.0;
        if (sum < 0.0) {
          sign = -1.0;
        }
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (std::size_t q = 0; q < count; ++q) {
          lo = std::min(lo, sign * det[q]);
          hi = std::max(hi, sign * det[q]);
        }
        out.min_jacobian[ue] = lo;
        out.max_jacobian[ue] = hi;
      },
      16);
  return out;
}

}  // namespace og
