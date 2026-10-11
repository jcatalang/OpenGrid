#include "opengrid/core/structured/elliptic3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <optional>
#include <stdexcept>

namespace og {
namespace {

using Idx = std::array<int, 3>;

// Scalar field on the block nodes.
struct Field3 {
  int ni = 0;
  int nj = 0;
  int nk = 0;
  std::vector<double> v;
  Field3(int a, int b, int c)
      : ni(a),
        nj(b),
        nk(c),
        v(static_cast<std::size_t>(a) * static_cast<std::size_t>(b) * static_cast<std::size_t>(c),
          0.0) {}
  double& operator()(int i, int j, int k) {
    return v[(static_cast<std::size_t>(k) * static_cast<std::size_t>(nj) +
              static_cast<std::size_t>(j)) *
                 static_cast<std::size_t>(ni) +
             static_cast<std::size_t>(i)];
  }
};

// Transfinite blend (index space) of a quantity known on the four faces around an axis.
double blend(double f0b, double f1b, double fa0, double fa1, double c00, double c10, double c01,
             double c11, double s, double t) {
  return (1.0 - s) * f0b + s * f1b + (1.0 - t) * fa0 + t * fa1 -
         ((1.0 - s) * (1.0 - t) * c00 + s * (1.0 - t) * c10 + (1.0 - s) * t * c01 + s * t * c11);
}

struct Wall {
  int axis = 2;
  int side = -1;  // -1 min face, +1 max face
  WallControl control;
  int na = 0;  // face size along its first / second tangential axis
  int nb = 0;
  int ta = 0;  // tangential axes
  int tb = 1;
  std::vector<double> spacing;           // per face node
  std::array<std::vector<double>, 3> s;  // sources per block axis, per face node
  bool clamped = false;
};

}  // namespace

Elliptic3DResult smooth_elliptic_3d(const Grid3D& initial, const Elliptic3DOptions& options,
                                    SmoothingProgress* progress) {
  const int ni = initial.ni;
  const int nj = initial.nj;
  const int nk = initial.nk;
  if (ni < 3 || nj < 3 || nk < 3) {
    throw std::invalid_argument("elliptic 3D: the block needs at least 3 x 3 x 3 points");
  }
  if (options.control == ControlType::Ttm) {
    throw std::invalid_argument("elliptic 3D: control must be Laplace or Thomas-Middlecoff");
  }
  const std::array<int, 3> n = {ni, nj, nk};
  Grid3D g = initial;
  auto at = [&](const Idx& p) -> Vec3& { return g.at(p[0], p[1], p[2]); };
  auto shift = [](Idx p, int axis, int d) {
    p[static_cast<std::size_t>(axis)] += d;
    return p;
  };
  double size = 0.0;
  for (const Vec3& p : g.x) {
    size = std::max(size, distance(p, g.x.front()));
  }

  // Background control functions (Thomas-Middlecoff): P[m] on every node.
  std::array<Field3, 3> background = {Field3(ni, nj, nk), Field3(ni, nj, nk), Field3(ni, nj, nk)};
  if (options.control == ControlType::ThomasMiddlecoff) {
    for (int m = 0; m < 3; ++m) {
      // Value of -(r_m . r_mm)/|r_m|^2 at a boundary node (interior along m; ends copied).
      auto tm = [&](Idx p) {
        const int c = p[static_cast<std::size_t>(m)];
        const int nm = n[static_cast<std::size_t>(m)];
        p[static_cast<std::size_t>(m)] = std::clamp(c, 1, nm - 2);
        const Vec3 rp = at(shift(p, m, 1));
        const Vec3 r0 = at(p);
        const Vec3 rm = at(shift(p, m, -1));
        const Vec3 d1 = 0.5 * (rp - rm);
        const Vec3 d2 = rp - 2.0 * r0 + rm;
        const double l2 = dot(d1, d1);
        if (l2 == 0.0) {
          return 0.0;
        }
        return -dot(d1, d2) / l2;
      };
      const int a = (m + 1) % 3;
      const int b = (m + 2) % 3;
      const int na = n[static_cast<std::size_t>(a)];
      const int nb = n[static_cast<std::size_t>(b)];
      for (int k = 0; k < nk; ++k) {
        for (int j = 0; j < nj; ++j) {
          for (int i = 0; i < ni; ++i) {
            const Idx p = {i, j, k};
            const int ia = p[static_cast<std::size_t>(a)];
            const int ib = p[static_cast<std::size_t>(b)];
            auto on = [&](int va, int vb) {
              Idx q = p;
              q[static_cast<std::size_t>(a)] = va;
              q[static_cast<std::size_t>(b)] = vb;
              return tm(q);
            };
            const double s = static_cast<double>(ia) / (na - 1);
            const double t = static_cast<double>(ib) / (nb - 1);
            background[static_cast<std::size_t>(m)](i, j, k) =
                blend(on(0, ib), on(na - 1, ib), on(ia, 0), on(ia, nb - 1), on(0, 0), on(na - 1, 0),
                      on(0, nb - 1), on(na - 1, nb - 1), s, t);
          }
        }
      }
    }
  }

  // Walls.
  std::vector<Wall> walls;
  for (const WallControl3D& w : options.walls) {
    if (w.face < 0 || w.face > 5) {
      throw std::invalid_argument("elliptic 3D: wall face must be 0..5");
    }
    Wall wall;
    wall.axis = w.face / 2;
    wall.side = -1;
    if (w.face % 2 == 1) {
      wall.side = 1;
    }
    wall.control = w.control;
    wall.ta = (wall.axis + 1) % 3;
    wall.tb = (wall.axis + 2) % 3;
    wall.na = n[static_cast<std::size_t>(wall.ta)];
    wall.nb = n[static_cast<std::size_t>(wall.tb)];
    if (n[static_cast<std::size_t>(wall.axis)] < 4) {
      throw std::invalid_argument("elliptic 3D: wall control needs at least 4 points across");
    }
    const auto count = static_cast<std::size_t>(wall.na) * static_cast<std::size_t>(wall.nb);
    for (auto& s : wall.s) {
      s.assign(count, 0.0);
    }
    wall.spacing.assign(count, 0.0);
    walls.push_back(wall);
  }
  auto wall_node = [&](const Wall& w, int ia, int ib, int depth) {
    Idx p{};
    p[static_cast<std::size_t>(w.ta)] = ia;
    p[static_cast<std::size_t>(w.tb)] = ib;
    int c = depth;
    if (w.side > 0) {
      c = n[static_cast<std::size_t>(w.axis)] - 1 - depth;
    }
    p[static_cast<std::size_t>(w.axis)] = c;
    return p;
  };
  for (Wall& w : walls) {
    // Requested first spacing: explicit, or bilinear (index space) from the first intervals of
    // the four block edges leaving the wall face's corners.
    std::array<double, 4> corner{};
    int m = 0;
    for (const int ib : {0, w.nb - 1}) {
      for (const int ia : {0, w.na - 1}) {
        corner[static_cast<std::size_t>(m++)] =
            distance(at(wall_node(w, ia, ib, 1)), at(wall_node(w, ia, ib, 0)));
      }
    }
    for (int ib = 0; ib < w.nb; ++ib) {
      for (int ia = 0; ia < w.na; ++ia) {
        const double s = static_cast<double>(ia) / (w.na - 1);
        const double t = static_cast<double>(ib) / (w.nb - 1);
        double value = (1 - s) * (1 - t) * corner[0] + s * (1 - t) * corner[1] +
                       (1 - s) * t * corner[2] + s * t * corner[3];
        const std::optional<double> requested = w.control.spacing;
        if (requested.has_value()) {
          value = *requested;
        }
        w.spacing[static_cast<std::size_t>(ib) * static_cast<std::size_t>(w.na) +
                  static_cast<std::size_t>(ia)] = value;
      }
    }
  }
  auto update_walls = [&]() {
    for (Wall& w : walls) {
      for (int ib = 1; ib < w.nb - 1; ++ib) {
        for (int ia = 1; ia < w.na - 1; ++ia) {
          const auto id = static_cast<std::size_t>(ib) * static_cast<std::size_t>(w.na) +
                          static_cast<std::size_t>(ia);
          const Idx p0 = wall_node(w, ia, ib, 0);
          const Vec3 line = at(wall_node(w, ia, ib, 1)) - at(p0);
          const double nl = norm(line);
          if (nl == 0.0) {
            continue;
          }
          const double error = (nl - w.spacing[id]) / w.spacing[id];
          const auto toward = static_cast<double>(w.side);  // -1 at a min face
          auto& sa = w.s[static_cast<std::size_t>(w.axis)][id];
          sa += toward * w.control.relaxation * std::tanh(error);
          for (const int t : {w.ta, w.tb}) {
            const Vec3 tangent = at(shift(p0, t, 1)) - at(shift(p0, t, -1));
            const double tl = norm(tangent);
            if (tl == 0.0) {
              continue;
            }
            const double lean = std::acos(std::clamp(dot(tangent, line) / (tl * nl), -1.0, 1.0)) -
                                0.5 * std::numbers::pi;
            auto& st = w.s[static_cast<std::size_t>(t)][id];
            st += w.control.relaxation * std::tanh(lean);
          }
          for (auto& s : w.s) {
            if (std::fabs(s[id]) > w.control.max_source) {
              w.clamped = true;
              s[id] = std::clamp(s[id], -w.control.max_source, w.control.max_source);
            }
          }
        }
      }
    }
  };
  auto sources = [&](const Idx& p, std::array<double, 3>& out) {
    for (int m = 0; m < 3; ++m) {
      out[static_cast<std::size_t>(m)] = background[static_cast<std::size_t>(m)](p[0], p[1], p[2]);
    }
    for (const Wall& w : walls) {
      int depth = p[static_cast<std::size_t>(w.axis)];
      if (w.side > 0) {
        depth = n[static_cast<std::size_t>(w.axis)] - 1 - depth;
      }
      const double weight = std::exp(-w.control.decay * depth);
      const auto id = static_cast<std::size_t>(p[static_cast<std::size_t>(w.tb)]) *
                          static_cast<std::size_t>(w.na) +
                      static_cast<std::size_t>(p[static_cast<std::size_t>(w.ta)]);
      for (int m = 0; m < 3; ++m) {
        out[static_cast<std::size_t>(m)] += weight * w.s[static_cast<std::size_t>(m)][id];
      }
    }
  };

  Elliptic3DResult result;
  double first = 0.0;
  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    if (progress != nullptr && progress->cancelled()) {
      result.cancelled = true;
      break;
    }
    update_walls();
    double sum = 0.0;
    std::size_t count = 0;
    for (int k = 1; k < nk - 1; ++k) {
      for (int j = 1; j < nj - 1; ++j) {
        for (int i = 1; i < ni - 1; ++i) {
          const Idx p = {i, j, k};
          std::array<Vec3, 3> d;
          for (int m = 0; m < 3; ++m) {
            d[static_cast<std::size_t>(m)] = 0.5 * (at(shift(p, m, 1)) - at(shift(p, m, -1)));
          }
          const double g11 = dot(d[0], d[0]);
          const double g22 = dot(d[1], d[1]);
          const double g33 = dot(d[2], d[2]);
          const double g12 = dot(d[0], d[1]);
          const double g13 = dot(d[0], d[2]);
          const double g23 = dot(d[1], d[2]);
          const std::array<double, 3> a = {g22 * g33 - g23 * g23, g11 * g33 - g13 * g13,
                                           g11 * g22 - g12 * g12};
          const double a12 = g13 * g23 - g12 * g33;
          const double a13 = g12 * g23 - g13 * g22;
          const double a23 = g12 * g13 - g11 * g23;
          auto cross_derivative = [&](int m1, int m2) {
            return 0.25 *
                   (at(shift(shift(p, m1, 1), m2, 1)) - at(shift(shift(p, m1, 1), m2, -1)) -
                    at(shift(shift(p, m1, -1), m2, 1)) + at(shift(shift(p, m1, -1), m2, -1)));
          };
          std::array<double, 3> src{};
          sources(p, src);
          Vec3 numerator = 2.0 * (a12 * cross_derivative(0, 1) + a13 * cross_derivative(0, 2) +
                                  a23 * cross_derivative(1, 2));
          double diag = 0.0;
          for (int m = 0; m < 3; ++m) {
            const auto um = static_cast<std::size_t>(m);
            numerator += a[um] * ((1.0 + 0.5 * src[um]) * at(shift(p, m, 1)) +
                                  (1.0 - 0.5 * src[um]) * at(shift(p, m, -1)));
            diag += 2.0 * a[um];
          }
          if (diag <= 0.0) {
            continue;
          }
          const Vec3 target = numerator / diag;
          const Vec3 delta = options.omega * (target - at(p));
          at(p) += delta;
          sum += dot(delta, delta);
          ++count;
        }
      }
    }
    const double rms = std::sqrt(sum / static_cast<double>(std::max<std::size_t>(count, 1)));
    result.residuals.push_back(rms);
    if (progress != nullptr) {
      progress->record(rms);
    }
    result.iterations = iteration + 1;
    if (iteration == 0) {
      first = rms;
    }
    if (rms <= options.tolerance * first || rms <= 1e-13 * std::max(size, 1e-300)) {
      result.converged = true;
      break;
    }
  }
  for (const Wall& w : walls) {
    if (w.clamped) {
      result.warnings.emplace_back("elliptic 3D: wall sources reached max_source (face " +
                                   std::to_string(w.axis * 2 + (w.side + 1) / 2) + ")");
    }
  }
  if (!result.converged && !result.cancelled) {
    result.warnings.emplace_back("elliptic 3D: not converged in " +
                                 std::to_string(options.max_iterations) + " iterations");
  }
  result.grid = std::move(g);
  return result;
}

}  // namespace og
