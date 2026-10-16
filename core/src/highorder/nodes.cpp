#include "opengrid/core/highorder/nodes.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>

namespace og {

namespace {

// P_n(x) and P_{n-1}(x) by the three-term recurrence.
void legendre(int n, double x, double& pn, double& pm) {
  double p0 = 1.0;
  double p1 = x;
  if (n == 0) {
    pn = 1.0;
    pm = 0.0;
    return;
  }
  for (int k = 1; k < n; ++k) {
    const double p2 = ((2.0 * k + 1.0) * x * p1 - k * p0) / (k + 1.0);
    p0 = p1;
    p1 = p2;
  }
  pn = p1;
  pm = p0;
}

// P'_n(x) for |x| < 1.
double legendre_derivative(int n, double x) {
  double pn = 0.0;
  double pm = 0.0;
  legendre(n, x, pn, pm);
  return n * (x * pn - pm) / (x * x - 1.0);
}

std::vector<double> barycentric_weights(const std::vector<double>& x) {
  std::vector<double> w(x.size(), 1.0);
  for (std::size_t j = 0; j < x.size(); ++j) {
    for (std::size_t k = 0; k < x.size(); ++k) {
      if (k != j) {
        w[j] /= x[j] - x[k];
      }
    }
  }
  return w;
}

constexpr double kSame = 1e-14;

}  // namespace

NodeSet node_set_from(const std::string& name) {
  if (name == "equispaced") {
    return NodeSet::Equispaced;
  }
  if (name == "gll" || name == "gauss_lobatto") {
    return NodeSet::GaussLobatto;
  }
  if (name == "gauss" || name == "gauss_legendre") {
    return NodeSet::GaussLegendre;
  }
  throw std::invalid_argument("node set must be 'equispaced', 'gll' or 'gauss', got '" + name +
                              "'");
}

std::vector<double> nodes_1d(NodeSet set, int count) {
  if (count < 1) {
    throw std::invalid_argument("node set: need at least one point");
  }
  std::vector<double> x(static_cast<std::size_t>(count), 0.0);
  if (set == NodeSet::Equispaced) {
    if (count == 1) {
      return x;
    }
    for (int k = 0; k < count; ++k) {
      x[static_cast<std::size_t>(k)] = -1.0 + 2.0 * k / (count - 1);
    }
    return x;
  }
  if (set == NodeSet::GaussLegendre) {
    for (int k = 0; k < count; ++k) {
      double t = std::cos(std::numbers::pi * (k + 0.75) / (count + 0.5));
      for (int it = 0; it < 100; ++it) {
        double pn = 0.0;
        double pm = 0.0;
        legendre(count, t, pn, pm);
        const double step = pn / legendre_derivative(count, t);
        t -= step;
        if (std::fabs(step) < 1e-16) {
          break;
        }
      }
      x[static_cast<std::size_t>(k)] = t;
    }
    std::sort(x.begin(), x.end());
    return x;
  }
  if (count < 2) {
    throw std::invalid_argument("Gauss-Lobatto points: need at least two");
  }
  const int n = count - 1;  // the interior points are the roots of P'_n
  x.front() = -1.0;
  x.back() = 1.0;
  for (int k = 1; k < n; ++k) {
    double t = -std::cos(std::numbers::pi * k / n);
    for (int it = 0; it < 100; ++it) {
      double pn = 0.0;
      double pm = 0.0;
      legendre(n, t, pn, pm);
      const double d1 = n * (t * pn - pm) / (t * t - 1.0);
      const double d2 = (2.0 * t * d1 - n * (n + 1.0) * pn) / (1.0 - t * t);
      const double step = d1 / d2;
      t -= step;
      if (std::fabs(step) < 1e-16) {
        break;
      }
    }
    x[static_cast<std::size_t>(k)] = t;
  }
  std::sort(x.begin(), x.end());
  return x;
}

std::vector<double> weights_1d(NodeSet set, int count) {
  const std::vector<double> x = nodes_1d(set, count);
  std::vector<double> w(x.size(), 0.0);
  if (set == NodeSet::GaussLegendre) {
    for (std::size_t k = 0; k < x.size(); ++k) {
      const double d = legendre_derivative(count, x[k]);
      w[k] = 2.0 / ((1.0 - x[k] * x[k]) * d * d);
    }
    return w;
  }
  if (set == NodeSet::GaussLobatto) {
    const int n = count - 1;
    for (std::size_t k = 0; k < x.size(); ++k) {
      double pn = 0.0;
      double pm = 0.0;
      legendre(n, x[k], pn, pm);
      w[k] = 2.0 / (n * (n + 1.0) * pn * pn);
    }
    return w;
  }
  throw std::invalid_argument("quadrature weights: only for Gauss-Lobatto or Gauss points");
}

std::vector<double> interpolation_matrix(const std::vector<double>& from,
                                         const std::vector<double>& to) {
  const std::vector<double> w = barycentric_weights(from);
  const std::size_t n = from.size();
  std::vector<double> m(to.size() * n, 0.0);
  for (std::size_t i = 0; i < to.size(); ++i) {
    double* row = &m[i * n];
    bool exact = false;
    for (std::size_t j = 0; j < n; ++j) {
      if (std::fabs(to[i] - from[j]) < kSame) {
        row[j] = 1.0;
        exact = true;
      }
    }
    if (exact) {
      continue;
    }
    double sum = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      row[j] = w[j] / (to[i] - from[j]);
      sum += row[j];
    }
    for (std::size_t j = 0; j < n; ++j) {
      row[j] /= sum;
    }
  }
  return m;
}

std::vector<double> differentiation_matrix(const std::vector<double>& from,
                                           const std::vector<double>& to) {
  const std::vector<double> w = barycentric_weights(from);
  const std::size_t n = from.size();
  std::vector<double> m(to.size() * n, 0.0);
  for (std::size_t i = 0; i < to.size(); ++i) {
    double* row = &m[i * n];
    std::size_t at = n;
    for (std::size_t j = 0; j < n; ++j) {
      if (std::fabs(to[i] - from[j]) < kSame) {
        at = j;
      }
    }
    if (at < n) {
      // At node `at`: D_aj = (w_j / w_a) / (x_a - x_j), D_aa = -sum of the others.
      double diag = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        if (j != at) {
          row[j] = (w[j] / w[at]) / (from[at] - from[j]);
          diag -= row[j];
        }
      }
      row[at] = diag;
      continue;
    }
    double s = 0.0;
    double ds = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      const double inv = 1.0 / (to[i] - from[j]);
      s += w[j] * inv;
      ds -= w[j] * inv * inv;
    }
    for (std::size_t j = 0; j < n; ++j) {
      const double inv = 1.0 / (to[i] - from[j]);
      const double lj = w[j] * inv / s;
      row[j] = (-w[j] * inv * inv - lj * ds) / s;
    }
  }
  return m;
}

}  // namespace og
