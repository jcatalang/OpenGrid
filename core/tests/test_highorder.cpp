// High-order elements (highorder/nodes.hpp, highorder/elements.hpp).
#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "opengrid/core/highorder/elements.hpp"
#include "opengrid/core/highorder/nodes.hpp"

namespace {

double integrate(const std::vector<double>& x, const std::vector<double>& w, int degree) {
  double sum = 0.0;
  for (std::size_t k = 0; k < x.size(); ++k) {
    sum += w[k] * std::pow(x[k], degree);
  }
  return sum;
}

double exact_moment(int degree) {
  if (degree % 2 == 1) {
    return 0.0;
  }
  return 2.0 / (degree + 1);
}

}  // namespace

TEST(HighOrderNodes, QuadratureIsExactToItsDegree) {
  for (int n = 2; n <= 12; ++n) {
    const std::vector<double> gx = og::nodes_1d(og::NodeSet::GaussLegendre, n);
    const std::vector<double> gw = og::weights_1d(og::NodeSet::GaussLegendre, n);
    const std::vector<double> lx = og::nodes_1d(og::NodeSet::GaussLobatto, n);
    const std::vector<double> lw = og::weights_1d(og::NodeSet::GaussLobatto, n);
    EXPECT_DOUBLE_EQ(lx.front(), -1.0);
    EXPECT_DOUBLE_EQ(lx.back(), 1.0);
    for (int d = 0; d <= 2 * n - 1; ++d) {
      EXPECT_NEAR(integrate(gx, gw, d), exact_moment(d), 1e-13) << "GL n " << n << " d " << d;
    }
    for (int d = 0; d <= 2 * n - 3; ++d) {
      EXPECT_NEAR(integrate(lx, lw, d), exact_moment(d), 1e-13) << "GLL n " << n << " d " << d;
    }
  }
  EXPECT_THROW(og::weights_1d(og::NodeSet::Equispaced, 4), std::invalid_argument);
  EXPECT_THROW(og::node_set_from("chebyshev"), std::invalid_argument);
}

TEST(HighOrderNodes, InterpolationAndDerivativesAreExactForPolynomials) {
  const std::vector<double> from = og::nodes_1d(og::NodeSet::Equispaced, 6);  // degree 5
  const std::vector<double> to = {-0.93, -0.4, 0.0, 0.2, 0.71, 1.0};
  const std::vector<double> m = og::interpolation_matrix(from, to);
  const std::vector<double> d = og::differentiation_matrix(from, to);
  const std::vector<double> dn = og::differentiation_matrix(from, from);
  auto f = [](double x) { return 2.0 - x + 3.0 * x * x * x - 0.5 * std::pow(x, 5); };
  auto df = [](double x) { return -1.0 + 9.0 * x * x - 2.5 * std::pow(x, 4); };
  for (std::size_t i = 0; i < to.size(); ++i) {
    double v = 0.0;
    double dv = 0.0;
    for (std::size_t j = 0; j < from.size(); ++j) {
      v += m[i * from.size() + j] * f(from[j]);
      dv += d[i * from.size() + j] * f(from[j]);
    }
    EXPECT_NEAR(v, f(to[i]), 1e-13);
    EXPECT_NEAR(dv, df(to[i]), 1e-12);
  }
  for (std::size_t i = 0; i < from.size(); ++i) {
    double dv = 0.0;
    for (std::size_t j = 0; j < from.size(); ++j) {
      dv += dn[i * from.size() + j] * f(from[j]);
    }
    EXPECT_NEAR(dv, df(from[i]), 1e-12);  // at the nodes themselves
  }
}

TEST(HighOrderElements, ReproducePolynomialMappingsAndTheirJacobian) {
  // A grid sampled from a quadratic map: order-2 elements represent it exactly, so the GLL
  // nodes are the map at the GLL points and the Jacobian is the map's.
  auto map = [](double s, double t) { return og::Vec3{s + 0.1 * t * t, t + 0.2 * s * t, 0.0}; };
  const int ni = 9;
  const int nj = 5;
  std::vector<og::Vec3> grid;
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      grid.push_back(map(static_cast<double>(i) / (ni - 1), static_cast<double>(j) / (nj - 1)));
    }
  }
  const og::HighOrderElements e =
      og::high_order_from_grid(grid, ni, nj, 1, 2, og::NodeSet::GaussLobatto);
  ASSERT_EQ(e.count, 8);
  const std::vector<double> gll = og::nodes_1d(og::NodeSet::GaussLobatto, 3);
  for (std::int64_t k = 0; k < e.count; ++k) {
    const std::int64_t ei = k % 4;
    const std::int64_t ej = k / 4;
    for (int b = 0; b < 3; ++b) {
      for (int a = 0; a < 3; ++a) {
        const double s =
            (static_cast<double>(ei) + 0.5 * (gll[static_cast<std::size_t>(a)] + 1.0)) / 4.0;
        const double t =
            (static_cast<double>(ej) + 0.5 * (gll[static_cast<std::size_t>(b)] + 1.0)) / 2.0;
        const og::Vec3 want = map(s, t);
        const std::size_t q = static_cast<std::size_t>(k) * 9 + static_cast<std::size_t>(b) * 3 +
                              static_cast<std::size_t>(a);
        EXPECT_NEAR(e.xyz[3 * q], want.x, 1e-14);
        EXPECT_NEAR(e.xyz[3 * q + 1], want.y, 1e-14);
      }
    }
    // Jacobian (in xi in [-1, 1]: ds/dxi = 1/8, dt/deta = 1/4): det = (1 + 0.2 s)/32 - ...
    EXPECT_GT(e.min_jacobian[static_cast<std::size_t>(k)], 0.0);
    EXPECT_LE(e.min_jacobian[static_cast<std::size_t>(k)],
              e.max_jacobian[static_cast<std::size_t>(k)]);
  }
  EXPECT_THROW(og::high_order_from_grid(grid, ni, nj, 1, 3, og::NodeSet::Equispaced),
               std::invalid_argument);  // 8 intervals are not a multiple of 3
}

TEST(HighOrderElements, FoldedElementIsReported) {
  // A 3 x 3 grid (one order-2 element) whose middle node is pushed past a corner.
  std::vector<og::Vec3> grid;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      grid.emplace_back(0.5 * i, 0.5 * j, 0.0);
    }
  }
  grid[4] = {1.4, 1.4, 0.0};
  const og::HighOrderElements e =
      og::high_order_from_grid(grid, 3, 3, 1, 2, og::NodeSet::Equispaced);
  EXPECT_LT(e.min_jacobian[0], 0.0);
}
