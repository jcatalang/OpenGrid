#include "opengrid/core/structured/tfi.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace og {
namespace {

std::vector<double> fractions(const std::vector<Vec3>& points) {
  std::vector<double> u(points.size(), 0.0);
  for (std::size_t k = 1; k < points.size(); ++k) {
    u[k] = u[k - 1] + distance(points[k], points[k - 1]);
  }
  const double total = u.back();
  if (total > 0.0) {
    for (double& value : u) {
      value /= total;
    }
  }
  u.back() = 1.0;
  return u;
}

double extent(const Boundary2D& b) {
  double size = 0.0;
  for (const auto* edge : {&b.bottom, &b.right, &b.top, &b.left}) {
    for (std::size_t k = 1; k < edge->size(); ++k) {
      size = std::max(size, distance((*edge)[k], edge->front()));
    }
  }
  return size;
}

void validate(const Boundary2D& b) {
  if (b.bottom.size() < 2 || b.left.size() < 2) {
    throw std::invalid_argument("tfi: each direction needs at least 2 points");
  }
  if (b.top.size() != b.bottom.size() || b.right.size() != b.left.size()) {
    throw std::invalid_argument(
        "tfi: opposite boundaries must have the same number of points (" +
        std::to_string(b.bottom.size()) + " vs " + std::to_string(b.top.size()) + ", " +
        std::to_string(b.left.size()) + " vs " + std::to_string(b.right.size()) + ")");
  }
  const double tol = 1e-9 * std::max(1.0, extent(b));
  const bool ok = distance(b.bottom.front(), b.left.front()) <= tol &&
                  distance(b.bottom.back(), b.right.front()) <= tol &&
                  distance(b.top.front(), b.left.back()) <= tol &&
                  distance(b.top.back(), b.right.back()) <= tol;
  if (!ok) {
    throw std::invalid_argument("tfi: boundary corners do not coincide");
  }
}

// Plane normal of the boundary loop (Newell's method); +z for an xy domain traversed CCW.
Vec3 plane_normal(const Boundary2D& b) {
  std::vector<Vec3> loop(b.bottom);
  loop.insert(loop.end(), b.right.begin() + 1, b.right.end());
  loop.insert(loop.end(), b.top.rbegin() + 1, b.top.rend());
  loop.insert(loop.end(), b.left.rbegin() + 1, b.left.rend() - 1);
  Vec3 n;
  for (std::size_t k = 0; k < loop.size(); ++k) {
    n += cross(loop[k], loop[(k + 1) % loop.size()]);
  }
  const double len = norm(n);
  if (len == 0.0) {
    return {0.0, 0.0, 1.0};
  }
  return n / len;
}

struct Params {
  std::vector<double> xi;   // per node
  std::vector<double> eta;  // per node
};

Params blending_params(const Boundary2D& b) {
  const auto ub = fractions(b.bottom);
  const auto ut = fractions(b.top);
  const auto vl = fractions(b.left);
  const auto vr = fractions(b.right);
  const std::size_t ni = b.bottom.size();
  const std::size_t nj = b.left.size();
  Params p;
  p.xi.resize(ni * nj);
  p.eta.resize(ni * nj);
  for (std::size_t j = 0; j < nj; ++j) {
    for (std::size_t i = 0; i < ni; ++i) {
      const double du = ut[i] - ub[i];
      const double dv = vr[j] - vl[j];
      const double xi = (ub[i] + vl[j] * du) / (1.0 - dv * du);
      const double eta = vl[j] + xi * dv;
      p.xi[j * ni + i] = xi;
      p.eta[j * ni + i] = eta;
    }
  }
  return p;
}

// Copy the boundary nodes exactly (TFI reproduces them only up to round-off).
// Sides first, so that the corners carry the bottom/top values.
void stamp_boundary(const Boundary2D& b, Grid2D& g) {
  for (int j = 0; j < g.nj; ++j) {
    g.at(0, j) = b.left[static_cast<std::size_t>(j)];
    g.at(g.ni - 1, j) = b.right[static_cast<std::size_t>(j)];
  }
  for (int i = 0; i < g.ni; ++i) {
    g.at(i, 0) = b.bottom[static_cast<std::size_t>(i)];
    g.at(i, g.nj - 1) = b.top[static_cast<std::size_t>(i)];
  }
}

}  // namespace

Grid2D tfi_linear(const Boundary2D& b) {
  validate(b);
  const int ni = static_cast<int>(b.bottom.size());
  const int nj = static_cast<int>(b.left.size());
  const Params p = blending_params(b);
  Grid2D g(ni, nj);
  const Vec3& c00 = b.bottom.front();
  const Vec3& c10 = b.bottom.back();
  const Vec3& c01 = b.top.front();
  const Vec3& c11 = b.top.back();
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      const std::size_t k = g.index(i, j);
      const double s = p.xi[k];
      const double t = p.eta[k];
      const auto ui = static_cast<std::size_t>(i);
      const auto uj = static_cast<std::size_t>(j);
      g.x[k] =
          (1.0 - t) * b.bottom[ui] + t * b.top[ui] + (1.0 - s) * b.left[uj] + s * b.right[uj] -
          ((1.0 - s) * (1.0 - t) * c00 + s * (1.0 - t) * c10 + (1.0 - s) * t * c01 + s * t * c11);
    }
  }
  stamp_boundary(b, g);
  return g;
}

Grid2D tfi_hermite(const Boundary2D& b, const HermiteOptions& options) {
  validate(b);
  const int ni = static_cast<int>(b.bottom.size());
  const int nj = static_cast<int>(b.left.size());
  const Params p = blending_params(b);
  const Vec3 n = plane_normal(b);

  // dx/deta on the bottom and top rows. At the end columns the side boundaries' own
  // eta-derivatives are used, so that the Boolean-sum correction P_xi - P_xi P_eta vanishes to
  // first order at the wall and does not spread the corner mismatch into the interior.
  const auto vl = fractions(b.left);
  const auto vr = fractions(b.right);
  const std::size_t n_side = b.left.size();
  auto side_derivative = [&](const std::vector<Vec3>& side, const std::vector<double>& v,
                             bool at_bottom) {
    if (at_bottom) {
      return (side[1] - side[0]) / v[1];
    }
    return (side[n_side - 1] - side[n_side - 2]) / (1.0 - v[n_side - 2]);
  };
  // Near the end columns the derivative direction turns smoothly from the side boundary's
  // direction to the wall normal (smoothstep over 10% of the columns), avoiding crossing
  // columns where the wall meets the sides at an oblique angle.
  const auto last_column = static_cast<double>(b.bottom.size() - 1);
  const double zone = std::max(1.0, std::round(0.1 * last_column));
  auto derivative = [&](const std::vector<Vec3>& edge, bool orthogonal, bool at_bottom) {
    std::vector<Vec3> d(edge.size());
    const Vec3 left_d = side_derivative(b.left, vl, at_bottom);
    const Vec3 right_d = side_derivative(b.right, vr, at_bottom);
    for (std::size_t i = 0; i < edge.size(); ++i) {
      const Vec3 chord = b.top[i] - b.bottom[i];
      if (!orthogonal) {
        d[i] = chord;
        continue;
      }
      if (i == 0) {
        d[i] = left_d;
        continue;
      }
      if (i + 1 == edge.size()) {
        d[i] = right_d;
        continue;
      }
      const Vec3 tangent = edge[i + 1] - edge[i - 1];
      Vec3 normal = cross(n, tangent);
      const double len = norm(normal);
      if (len == 0.0) {
        d[i] = chord;
        continue;
      }
      normal = normal / len;
      const auto fi = static_cast<double>(i);
      double w = std::min(1.0, std::min(fi, last_column - fi) / zone);
      w = w * w * (3.0 - 2.0 * w);
      Vec3 side = left_d;
      if (fi > 0.5 * last_column) {
        side = right_d;
      }
      const double side_len = norm(side);
      Vec3 direction = normal;
      if (side_len > 0.0) {
        direction = (1.0 - w) * (side / side_len) + w * normal;
      }
      const double dlen = norm(direction);
      if (dlen == 0.0) {
        d[i] = chord;
        continue;
      }
      d[i] = direction * (options.strength * norm(chord) / dlen);
    }
    return d;
  };
  const std::vector<Vec3> db = derivative(b.bottom, options.bottom, true);
  const std::vector<Vec3> dt = derivative(b.top, options.top, false);

  // P_eta at column i for parameter t.
  auto p_eta = [&](std::size_t i, double t) {
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double a0 = 2.0 * t3 - 3.0 * t2 + 1.0;
    const double a1 = -2.0 * t3 + 3.0 * t2;
    const double b0 = t3 - 2.0 * t2 + t;
    const double b1 = t3 - t2;
    return a0 * b.bottom[i] + a1 * b.top[i] + b0 * db[i] + b1 * dt[i];
  };
  Grid2D g(ni, nj);
  const auto last = static_cast<std::size_t>(ni - 1);
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      const std::size_t k = g.index(i, j);
      const double s = p.xi[k];
      const double t = p.eta[k];
      const auto ui = static_cast<std::size_t>(i);
      const auto uj = static_cast<std::size_t>(j);
      const Vec3 pe = p_eta(ui, t);
      const Vec3 px = (1.0 - s) * b.left[uj] + s * b.right[uj];
      const Vec3 pxe = (1.0 - s) * p_eta(0, t) + s * p_eta(last, t);
      g.x[k] = pe + px - pxe;
    }
  }
  stamp_boundary(b, g);
  return g;
}

}  // namespace og
