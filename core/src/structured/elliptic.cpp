#include "opengrid/core/structured/elliptic.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <tuple>

#include "opengrid/core/util/anderson.hpp"
#include "opengrid/core/util/parallel.hpp"

namespace og {
namespace {

struct Field {
  int ni;
  int nj;
  std::vector<double> v;
  Field(int ni_, int nj_)
      : ni(ni_), nj(nj_), v(static_cast<std::size_t>(ni_) * static_cast<std::size_t>(nj_), 0.0) {}
  double& operator()(int i, int j) {
    return v[static_cast<std::size_t>(j) * static_cast<std::size_t>(ni) +
             static_cast<std::size_t>(i)];
  }
  double operator()(int i, int j) const {
    return v[static_cast<std::size_t>(j) * static_cast<std::size_t>(ni) +
             static_cast<std::size_t>(i)];
  }
};

double sign(double value) {
  if (value > 0.0) {
    return 1.0;
  }
  if (value < 0.0) {
    return -1.0;
  }
  return 0.0;
}

// Background control functions P0, Q0 (fixed during the iteration).
void background_control(const Field& x, const Field& y, const EllipticOptions& options, Field& p,
                        Field& q) {
  const int ni = x.ni;
  const int nj = x.nj;
  if (options.control == ControlType::ThomasMiddlecoff) {
    auto phi = [&](int i, int j) {
      const double xs = 0.5 * (x(i + 1, j) - x(i - 1, j));
      const double ys = 0.5 * (y(i + 1, j) - y(i - 1, j));
      const double xss = x(i + 1, j) - 2.0 * x(i, j) + x(i - 1, j);
      const double yss = y(i + 1, j) - 2.0 * y(i, j) + y(i - 1, j);
      const double d = xs * xs + ys * ys;
      if (d == 0.0) {
        return 0.0;
      }
      return -(xs * xss + ys * yss) / d;
    };
    auto psi = [&](int i, int j) {
      const double xt = 0.5 * (x(i, j + 1) - x(i, j - 1));
      const double yt = 0.5 * (y(i, j + 1) - y(i, j - 1));
      const double xtt = x(i, j + 1) - 2.0 * x(i, j) + x(i, j - 1);
      const double ytt = y(i, j + 1) - 2.0 * y(i, j) + y(i, j - 1);
      const double d = xt * xt + yt * yt;
      if (d == 0.0) {
        return 0.0;
      }
      return -(xt * xtt + yt * ytt) / d;
    };
    for (int i = 1; i < ni - 1; ++i) {
      const double bottom = phi(i, 0);
      const double top = phi(i, nj - 1);
      for (int j = 0; j < nj; ++j) {
        const double t = static_cast<double>(j) / (nj - 1);
        p(i, j) = (1.0 - t) * bottom + t * top;
      }
    }
    for (int j = 1; j < nj - 1; ++j) {
      const double left = psi(0, j);
      const double right = psi(ni - 1, j);
      for (int i = 0; i < ni; ++i) {
        const double s = static_cast<double>(i) / (ni - 1);
        q(i, j) = (1.0 - s) * left + s * right;
      }
    }
  } else if (options.control == ControlType::Ttm) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        double pv = 0.0;
        double qv = 0.0;
        for (const LineAttraction& line : options.lines) {
          if (line.xi_line) {
            const double d = i - line.position;
            pv -= line.amplitude * sign(d) * std::exp(-line.decay * std::fabs(d));
          } else {
            const double d = j - line.position;
            qv -= line.amplitude * sign(d) * std::exp(-line.decay * std::fabs(d));
          }
        }
        for (const PointAttraction& point : options.points) {
          const double di = i - point.xi;
          const double dj = j - point.eta;
          const double r = std::sqrt(di * di + dj * dj);
          pv -= point.amplitude * sign(di) * std::exp(-point.decay * r);
          qv -= point.amplitude * sign(dj) * std::exp(-point.decay * r);
        }
        p(i, j) = pv;
        q(i, j) = qv;
      }
    }
  }
}

struct WallState {
  std::vector<double> p;
  std::vector<double> q;
  std::vector<double> spacing;
};

// Requested first spacing along a wall row: explicit, or interpolated from the side
// boundaries' first intervals.
std::vector<double> wall_spacing(const Field& x, const Field& y, const WallControl& control,
                                 bool bottom) {
  const int ni = x.ni;
  const int nj = x.nj;
  std::vector<double> result(static_cast<std::size_t>(ni), 0.0);
  if (control.spacing) {
    std::fill(result.begin(), result.end(), *control.spacing);
    return result;
  }
  int j0 = 0;
  int j1 = 1;
  if (!bottom) {
    j0 = nj - 1;
    j1 = nj - 2;
  }
  const double left = std::hypot(x(0, j1) - x(0, j0), y(0, j1) - y(0, j0));
  const double right = std::hypot(x(ni - 1, j1) - x(ni - 1, j0), y(ni - 1, j1) - y(ni - 1, j0));
  // Arc-length fraction along the wall row.
  std::vector<double> s(static_cast<std::size_t>(ni), 0.0);
  for (int i = 1; i < ni; ++i) {
    s[static_cast<std::size_t>(i)] = s[static_cast<std::size_t>(i - 1)] +
                                     std::hypot(x(i, j0) - x(i - 1, j0), y(i, j0) - y(i - 1, j0));
  }
  for (int i = 0; i < ni; ++i) {
    double f = 0.0;
    if (s.back() > 0.0) {
      f = s[static_cast<std::size_t>(i)] / s.back();
    }
    result[static_cast<std::size_t>(i)] = (1.0 - f) * left + f * right;
  }
  return result;
}

// Wall-source update by direct feedback (Hilgenstock): q is driven by the relative error of the
// first-cell height, p by the deviation of the first grid line from the wall normal. At a
// converged state both errors vanish, so spacing and orthogonality are met exactly (up to
// the clamp on the sources).
void update_wall(const Field& x, const Field& y, const WallControl& control, bool bottom,
                 WallState& state) {
  const int ni = x.ni;
  const int nj = x.nj;
  int j0 = 0;
  int j1 = 1;
  double toward_wall = -1.0;  // sign of q that pulls nodes toward this wall
  if (!bottom) {
    j0 = nj - 1;
    j1 = nj - 2;
    toward_wall = 1.0;
  }
  for (int i = 1; i < ni - 1; ++i) {
    const auto ui = static_cast<std::size_t>(i);
    const double tx = x(i + 1, j0) - x(i - 1, j0);
    const double ty = y(i + 1, j0) - y(i - 1, j0);
    const double nx = x(i, j1) - x(i, j0);
    const double ny = y(i, j1) - y(i, j0);
    const double tl = std::hypot(tx, ty);
    const double nl = std::hypot(nx, ny);
    if (tl == 0.0 || nl == 0.0) {
      continue;
    }
    // Angle between the wall tangent (+xi) and the first grid line, minus 90 degrees (radians):
    // positive when the line leans towards -xi.
    const double lean =
        std::acos(std::clamp((tx * nx + ty * ny) / (tl * nl), -1.0, 1.0)) - 0.5 * std::numbers::pi;
    const double spacing_error = (nl - state.spacing[ui]) / state.spacing[ui];
    state.p[ui] = std::clamp(state.p[ui] + control.relaxation * std::tanh(lean),
                             -control.max_source, control.max_source);
    state.q[ui] =
        std::clamp(state.q[ui] + toward_wall * control.relaxation * std::tanh(spacing_error),
                   -control.max_source, control.max_source);
  }
}

// Thomas algorithm on k = 1..n-2 of a line with Dirichlet ends (k = 0 and k = n-1 hold the
// boundary values in rx/ry and are not modified). Returns false on a zero pivot.
bool solve_line(const double* la, double* lb, const double* lc, double* rx, double* ry,
                std::size_t n) {
  const std::size_t last = n - 2;
  rx[1] -= la[1] * rx[0];
  ry[1] -= la[1] * ry[0];
  rx[last] -= lc[last] * rx[n - 1];
  ry[last] -= lc[last] * ry[n - 1];
  for (std::size_t k = 2; k <= last; ++k) {
    if (!(lb[k - 1] != 0.0)) {
      return false;
    }
    const double w = la[k] / lb[k - 1];
    lb[k] -= w * lc[k - 1];
    rx[k] -= w * rx[k - 1];
    ry[k] -= w * ry[k - 1];
  }
  if (!(lb[last] != 0.0)) {
    return false;
  }
  rx[last] /= lb[last];
  ry[last] /= lb[last];
  for (std::size_t k = last - 1; k >= 1; --k) {
    rx[k] = (rx[k] - lc[k] * rx[k + 1]) / lb[k];
    ry[k] = (ry[k] - lc[k] * ry[k + 1]) / lb[k];
  }
  return true;
}

}  // namespace

EllipticResult smooth_elliptic(const Grid2D& initial, const EllipticOptions& options,
                               SmoothingProgress* progress) {
  const int ni = initial.ni;
  const int nj = initial.nj;
  if (ni < 3 || nj < 3) {
    throw std::invalid_argument("elliptic: the grid needs at least 3 x 3 points");
  }
  const double z0 = initial.x.front().z;
  double size = 0.0;
  for (const Vec3& point : initial.x) {
    size = std::max(size, distance(point, initial.x.front()));
  }
  for (const Vec3& point : initial.x) {
    if (std::fabs(point.z - z0) > 1e-9 * std::max(1.0, size)) {
      throw std::invalid_argument("elliptic: only planar grids in a z = const plane are supported");
    }
  }
  if ((options.bottom.enabled || options.top.enabled) && nj < 4) {
    throw std::invalid_argument("elliptic: wall control needs at least 4 points across");
  }

  Field x(ni, nj);
  Field y(ni, nj);
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      x(i, j) = initial.at(i, j).x;
      y(i, j) = initial.at(i, j).y;
    }
  }
  Field p0(ni, nj);
  Field q0(ni, nj);
  background_control(x, y, options, p0, q0);

  WallState bottom;
  WallState top;
  for (WallState* state : {&bottom, &top}) {
    state->p.assign(static_cast<std::size_t>(ni), 0.0);
    state->q.assign(static_cast<std::size_t>(ni), 0.0);
  }
  if (options.bottom.enabled) {
    bottom.spacing = wall_spacing(x, y, options.bottom, true);
  }
  if (options.top.enabled) {
    top.spacing = wall_spacing(x, y, options.top, false);
  }

  // Decay of the wall sources with the line index (computed once).
  std::vector<double> bottom_weight(static_cast<std::size_t>(nj), 0.0);
  std::vector<double> top_weight(static_cast<std::size_t>(nj), 0.0);
  for (int j = 0; j < nj; ++j) {
    bottom_weight[static_cast<std::size_t>(j)] = std::exp(-options.bottom.decay * j);
    top_weight[static_cast<std::size_t>(j)] = std::exp(-options.top.decay * (nj - 1 - j));
  }
  // Total control at (i, j) = background + decaying wall sources.
  auto control = [&](int i, int j, double& p, double& q) {
    p = p0(i, j);
    q = q0(i, j);
    if (options.bottom.enabled) {
      const double w = bottom_weight[static_cast<std::size_t>(j)];
      p += bottom.p[static_cast<std::size_t>(i)] * w;
      q += bottom.q[static_cast<std::size_t>(i)] * w;
    }
    if (options.top.enabled) {
      const double w = top_weight[static_cast<std::size_t>(j)];
      p += top.p[static_cast<std::size_t>(i)] * w;
      q += top.q[static_cast<std::size_t>(i)] * w;
    }
  };

  // Discrete equation at (i, j), central differences with unit computational spacing:
  //   aw x_{i-1} + ae x_{i+1} + as x_{j-1} + an x_{j+1} - 2 (alpha + gamma) x + cross = 0,
  //   ae/aw = alpha (1 +/- P/2), an/as = gamma (1 +/- Q/2), cross = -2 beta x_xieta.
  // Source terms enter the neighbour coefficients, so line relaxation treats them implicitly.
  struct Stencil {
    double aw;
    double ae;
    double as;
    double an;
    double diag;
    double cross_x;
    double cross_y;
  };
  auto stencil = [&](int i, int j) {
    const double xs = 0.5 * (x(i + 1, j) - x(i - 1, j));
    const double ys = 0.5 * (y(i + 1, j) - y(i - 1, j));
    const double xt = 0.5 * (x(i, j + 1) - x(i, j - 1));
    const double yt = 0.5 * (y(i, j + 1) - y(i, j - 1));
    const double alpha = xt * xt + yt * yt;
    const double beta = xs * xt + ys * yt;
    const double gamma = xs * xs + ys * ys;
    const double xst =
        0.25 * (x(i + 1, j + 1) - x(i + 1, j - 1) - x(i - 1, j + 1) + x(i - 1, j - 1));
    const double yst =
        0.25 * (y(i + 1, j + 1) - y(i + 1, j - 1) - y(i - 1, j + 1) + y(i - 1, j - 1));
    double p = 0.0;
    double q = 0.0;
    control(i, j, p, q);
    Stencil s{};
    s.aw = alpha * (1.0 - 0.5 * p);
    s.ae = alpha * (1.0 + 0.5 * p);
    s.as = gamma * (1.0 - 0.5 * q);
    s.an = gamma * (1.0 + 0.5 * q);
    s.diag = 2.0 * (alpha + gamma);
    s.cross_x = -2.0 * beta * xst;
    s.cross_y = -2.0 * beta * yst;
    return s;
  };

  // Line-solver scratch, one row per line so lines of one colour can be solved concurrently,
  // and the squared corrections per line (or point row), summed in order afterwards.
  const auto line_storage = static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj);
  std::vector<double> la(line_storage);
  std::vector<double> lb(line_storage);
  std::vector<double> lc(line_storage);
  std::vector<double> rx(line_storage);
  std::vector<double> ry(line_storage);
  std::vector<double> line_sum(static_cast<std::size_t>(std::max(ni, nj)), 0.0);

  // Iteration state for the acceleration: x, y and the wall sources.
  const std::size_t field_size = x.v.size();
  const std::size_t state_size = 2 * field_size + 4 * static_cast<std::size_t>(ni);
  auto pack = [&](std::vector<double>& u) {
    u.resize(state_size);
    std::copy(x.v.begin(), x.v.end(), u.begin());
    std::copy(y.v.begin(), y.v.end(), u.begin() + static_cast<std::ptrdiff_t>(field_size));
    auto at = u.begin() + static_cast<std::ptrdiff_t>(2 * field_size);
    for (const std::vector<double>* source : {&bottom.p, &bottom.q, &top.p, &top.q}) {
      at = std::copy(source->begin(), source->end(), at);
    }
  };
  auto unpack = [&](const std::vector<double>& u) {
    std::copy(u.begin(), u.begin() + static_cast<std::ptrdiff_t>(field_size), x.v.begin());
    std::copy(u.begin() + static_cast<std::ptrdiff_t>(field_size),
              u.begin() + static_cast<std::ptrdiff_t>(2 * field_size), y.v.begin());
    auto at = u.begin() + static_cast<std::ptrdiff_t>(2 * field_size);
    for (std::vector<double>* source : {&bottom.p, &bottom.q, &top.p, &top.q}) {
      std::copy(at, at + static_cast<std::ptrdiff_t>(ni), source->begin());
      at += static_cast<std::ptrdiff_t>(ni);
    }
  };
  AndersonAccelerator anderson(state_size, options.acceleration);
  std::vector<double> before;
  std::vector<double> after;

  EllipticResult result;
  double first_residual = 0.0;
  auto interior = static_cast<double>((ni - 2) * (nj - 2));
  if (options.solver == EllipticSolver::LineSor) {
    interior *= 2.0;  // two sweeps per iteration
  }
  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    if (progress != nullptr && progress->cancelled()) {
      result.cancelled = true;
      break;
    }
    if (options.acceleration > 0) {
      pack(before);
    }
    if (options.bottom.enabled) {
      update_wall(x, y, options.bottom, true, bottom);
    }
    if (options.top.enabled) {
      update_wall(x, y, options.top, false, top);
    }
    // Correction of node (i, j); returns its squared length.
    auto apply = [&](int i, int j, double gx, double gy) {
      const double cx = gx - x(i, j);
      const double cy = gy - y(i, j);
      x(i, j) += options.omega * cx;
      y(i, j) += options.omega * cy;
      return cx * cx + cy * cy;
    };
    std::fill(line_sum.begin(), line_sum.end(), 0.0);
    if (options.solver == EllipticSolver::PointSor) {
      // Zebra rows: odd j, then even j (the cross term reaches the rows j +/- 1, so points of
      // one row cannot be split by colour); the points of a row are relaxed in order.
      for (int colour = 1; colour >= 0; --colour) {
        const std::int64_t count = (nj - 2 - (1 - colour) + 1) / 2;
        parallel_for(
            0, count,
            [&](std::int64_t m) {
              const auto j = static_cast<int>(2 * m + 2 - colour);
              double sum = 0.0;
              for (int i = 1; i < ni - 1; ++i) {
                const Stencil s = stencil(i, j);
                if (!(s.diag > 0.0)) {
                  continue;
                }
                const double gx = (s.aw * x(i - 1, j) + s.ae * x(i + 1, j) + s.as * x(i, j - 1) +
                                   s.an * x(i, j + 1) + s.cross_x) /
                                  s.diag;
                const double gy = (s.aw * y(i - 1, j) + s.ae * y(i + 1, j) + s.as * y(i, j - 1) +
                                   s.an * y(i, j + 1) + s.cross_y) /
                                  s.diag;
                sum += apply(i, j, gx, gy);
              }
              line_sum[static_cast<std::size_t>(j)] = sum;
            },
            8);
      }
    } else {
      // Alternating zebra line relaxation: eta-lines (implicit in j) at odd, then even i; then
      // xi-lines (implicit in i) at odd, then even j.
      const auto unj = static_cast<std::size_t>(nj);
      const auto uni = static_cast<std::size_t>(ni);
      auto eta_line = [&](int i) {
        const std::size_t o = static_cast<std::size_t>(i) * unj;
        for (int j = 1; j < nj - 1; ++j) {
          const Stencil s = stencil(i, j);
          const std::size_t k = o + static_cast<std::size_t>(j);
          la[k] = -s.as;
          lb[k] = s.diag;
          lc[k] = -s.an;
          rx[k] = s.aw * x(i - 1, j) + s.ae * x(i + 1, j) + s.cross_x;
          ry[k] = s.aw * y(i - 1, j) + s.ae * y(i + 1, j) + s.cross_y;
        }
        rx[o] = x(i, 0);
        ry[o] = y(i, 0);
        rx[o + unj - 1] = x(i, nj - 1);
        ry[o + unj - 1] = y(i, nj - 1);
        if (!solve_line(&la[o], &lb[o], &lc[o], &rx[o], &ry[o], unj)) {
          return 0.0;
        }
        double sum = 0.0;
        for (int j = 1; j < nj - 1; ++j) {
          sum +=
              apply(i, j, rx[o + static_cast<std::size_t>(j)], ry[o + static_cast<std::size_t>(j)]);
        }
        return sum;
      };
      auto xi_line = [&](int j) {
        const std::size_t o = static_cast<std::size_t>(j) * uni;
        for (int i = 1; i < ni - 1; ++i) {
          const Stencil s = stencil(i, j);
          const std::size_t k = o + static_cast<std::size_t>(i);
          la[k] = -s.aw;
          lb[k] = s.diag;
          lc[k] = -s.ae;
          rx[k] = s.as * x(i, j - 1) + s.an * x(i, j + 1) + s.cross_x;
          ry[k] = s.as * y(i, j - 1) + s.an * y(i, j + 1) + s.cross_y;
        }
        rx[o] = x(0, j);
        ry[o] = y(0, j);
        rx[o + uni - 1] = x(ni - 1, j);
        ry[o + uni - 1] = y(ni - 1, j);
        if (!solve_line(&la[o], &lb[o], &lc[o], &rx[o], &ry[o], uni)) {
          return 0.0;
        }
        double sum = 0.0;
        for (int i = 1; i < ni - 1; ++i) {
          sum +=
              apply(i, j, rx[o + static_cast<std::size_t>(i)], ry[o + static_cast<std::size_t>(i)]);
        }
        return sum;
      };
      std::vector<double> eta_sum(uni, 0.0);
      for (int colour = 1; colour >= 0; --colour) {
        const std::int64_t count = (ni - 2 - (1 - colour) + 1) / 2;
        parallel_for(
            0, count,
            [&](std::int64_t m) {
              const auto i = static_cast<int>(2 * m + 2 - colour);
              eta_sum[static_cast<std::size_t>(i)] = eta_line(i);
            },
            8);
      }
      for (int colour = 1; colour >= 0; --colour) {
        const std::int64_t count = (nj - 2 - (1 - colour) + 1) / 2;
        parallel_for(
            0, count,
            [&](std::int64_t m) {
              const auto j = static_cast<int>(2 * m + 2 - colour);
              line_sum[static_cast<std::size_t>(j)] = xi_line(j);
            },
            8);
      }
      for (const double value : eta_sum) {
        line_sum[0] += value;
      }
    }
    double sum = 0.0;
    for (const double value : line_sum) {
      sum += value;
    }
    const double residual = std::sqrt(sum / interior);
    if (!std::isfinite(residual) && options.acceleration > 0) {
      // An extrapolation blew up within one step: back to the best iterate.
      pack(after);
      if (anderson.guarded_step(before, after, residual)) {
        unpack(after);
        result.iterations = iteration + 1;
        continue;
      }
    }
    if (!std::isfinite(residual)) {
      result.warnings.emplace_back("elliptic: iteration diverged; returning the initial grid");
      result.grid = initial;
      result.iterations = iteration + 1;
      return result;
    }
    result.residuals.push_back(residual);
    if (progress != nullptr) {
      progress->record(residual);
    }
    if (iteration == 0) {
      first_residual = residual;
    }
    if (options.acceleration > 0) {
      pack(after);
      if (anderson.guarded_step(before, after, residual)) {
        unpack(after);
      }
    }
    result.iterations = iteration + 1;
    // Relative convergence, with an absolute floor at round-off level of the grid size (a grid
    // that already solves the equations starts at a round-off residual).
    if (residual <= options.tolerance * first_residual || residual <= 1e-13 * size) {
      result.converged = true;
      break;
    }
  }
  // A source stuck at its clamp means the wall targets could not be met: the iteration can
  // still settle ("converge"), so report it explicitly.
  for (const auto& [wall_state, wall_control, label] :
       {std::tuple{&bottom, &options.bottom, "j = 0"},
        std::tuple{&top, &options.top, "j = jmax"}}) {
    if (!wall_control->enabled) {
      continue;
    }
    int saturated = 0;
    for (std::size_t i = 0; i < wall_state->p.size(); ++i) {
      if (std::fabs(wall_state->p[i]) >= wall_control->max_source ||
          std::fabs(wall_state->q[i]) >= wall_control->max_source) {
        saturated += 1;
      }
    }
    if (saturated > 0) {
      result.warnings.push_back("elliptic: wall control at " + std::string(label) +
                                " saturated at " + std::to_string(saturated) +
                                " node(s); spacing/orthogonality may not be met (raise max_source "
                                "or lower decay)");
    }
  }
  if (anderson.restarts() > 0) {
    result.warnings.push_back("elliptic: " + anderson.restart_note());
  }
  if (!result.converged && !result.cancelled) {
    result.warnings.push_back("elliptic: not converged after " + std::to_string(result.iterations) +
                              " iterations");
  }
  result.grid = initial;
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      result.grid.at(i, j).x = x(i, j);
      result.grid.at(i, j).y = y(i, j);
    }
  }
  return result;
}

}  // namespace og
