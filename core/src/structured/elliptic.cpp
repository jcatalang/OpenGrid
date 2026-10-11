#include "opengrid/core/structured/elliptic.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <tuple>

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

  // Total control at (i, j) = background + decaying wall sources.
  auto control = [&](int i, int j, double& p, double& q) {
    p = p0(i, j);
    q = q0(i, j);
    if (options.bottom.enabled) {
      const double w = std::exp(-options.bottom.decay * j);
      p += bottom.p[static_cast<std::size_t>(i)] * w;
      q += bottom.q[static_cast<std::size_t>(i)] * w;
    }
    if (options.top.enabled) {
      const double w = std::exp(-options.top.decay * (nj - 1 - j));
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

  // Thomas algorithm on k = 1..n-2 of a line with Dirichlet ends (k = 0 and k = n-1 hold the
  // boundary values in rx/ry and are not modified). Returns false on a zero pivot.
  std::vector<double> la;
  std::vector<double> lb;
  std::vector<double> lc;
  std::vector<double> rx;
  std::vector<double> ry;
  auto solve_line = [&](std::size_t n) {
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
  };
  auto resize = [&](int n) {
    for (auto* v : {&la, &lb, &lc, &rx, &ry}) {
      v->assign(static_cast<std::size_t>(n), 0.0);
    }
  };

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
    if (options.bottom.enabled) {
      update_wall(x, y, options.bottom, true, bottom);
    }
    if (options.top.enabled) {
      update_wall(x, y, options.top, false, top);
    }
    double sum = 0.0;
    auto apply = [&](int i, int j, double gx, double gy) {
      const double cx = gx - x(i, j);
      const double cy = gy - y(i, j);
      sum += cx * cx + cy * cy;
      x(i, j) += options.omega * cx;
      y(i, j) += options.omega * cy;
    };
    if (options.solver == EllipticSolver::PointSor) {
      for (int j = 1; j < nj - 1; ++j) {
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
          apply(i, j, gx, gy);
        }
      }
    } else {
      // Alternating line relaxation: eta-lines (implicit in j) for every interior i, then
      // xi-lines (implicit in i) for every interior j.
      resize(nj);
      for (int i = 1; i < ni - 1; ++i) {
        for (int j = 1; j < nj - 1; ++j) {
          const Stencil s = stencil(i, j);
          const auto k = static_cast<std::size_t>(j);
          la[k] = -s.as;
          lb[k] = s.diag;
          lc[k] = -s.an;
          rx[k] = s.aw * x(i - 1, j) + s.ae * x(i + 1, j) + s.cross_x;
          ry[k] = s.aw * y(i - 1, j) + s.ae * y(i + 1, j) + s.cross_y;
        }
        rx[0] = x(i, 0);
        ry[0] = y(i, 0);
        rx[static_cast<std::size_t>(nj - 1)] = x(i, nj - 1);
        ry[static_cast<std::size_t>(nj - 1)] = y(i, nj - 1);
        if (!solve_line(static_cast<std::size_t>(nj))) {
          continue;
        }
        for (int j = 1; j < nj - 1; ++j) {
          apply(i, j, rx[static_cast<std::size_t>(j)], ry[static_cast<std::size_t>(j)]);
        }
      }
      resize(ni);
      for (int j = 1; j < nj - 1; ++j) {
        for (int i = 1; i < ni - 1; ++i) {
          const Stencil s = stencil(i, j);
          const auto k = static_cast<std::size_t>(i);
          la[k] = -s.aw;
          lb[k] = s.diag;
          lc[k] = -s.ae;
          rx[k] = s.as * x(i, j - 1) + s.an * x(i, j + 1) + s.cross_x;
          ry[k] = s.as * y(i, j - 1) + s.an * y(i, j + 1) + s.cross_y;
        }
        rx[0] = x(0, j);
        ry[0] = y(0, j);
        rx[static_cast<std::size_t>(ni - 1)] = x(ni - 1, j);
        ry[static_cast<std::size_t>(ni - 1)] = y(ni - 1, j);
        if (!solve_line(static_cast<std::size_t>(ni))) {
          continue;
        }
        for (int i = 1; i < ni - 1; ++i) {
          apply(i, j, rx[static_cast<std::size_t>(i)], ry[static_cast<std::size_t>(i)]);
        }
      }
    }
    const double residual = std::sqrt(sum / interior);
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
