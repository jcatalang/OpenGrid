#include "opengrid/core/structured/surface_elliptic.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace og {
EllipticResult smooth_on_surface(const Surface& surface, const Grid2D& initial,
                                 const EllipticOptions& options, SmoothingProgress* progress) {
  const int ni = initial.ni;
  const int nj = initial.nj;
  if (ni < 3 || nj < 3) {
    throw std::invalid_argument("surface elliptic: the grid needs at least 3 x 3 points");
  }
  if (options.control == ControlType::Ttm || options.bottom.enabled || options.top.enabled ||
      !options.lines.empty() || !options.points.empty()) {
    throw std::invalid_argument(
        "surface elliptic: only Laplace and Thomas-Middlecoff control are supported on surfaces");
  }
  Grid2D g = initial;
  // Parameters of every node (projection once; then Newton from the previous parameters).
  std::vector<double> u(g.x.size());
  std::vector<double> v(g.x.size());
  double size = 0.0;
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    const SurfaceProjection q = surface.project(g.x[k]);
    u[k] = q.u;
    v[k] = q.v;
    size = std::max(size, distance(g.x[k], g.x.front()));
  }
  // Thomas-Middlecoff control functions from the boundary distributions (3D dot products).
  std::vector<double> p(g.x.size(), 0.0);
  std::vector<double> q(g.x.size(), 0.0);
  if (options.control == ControlType::ThomasMiddlecoff) {
    auto along_i = [&](int i, int j) {
      i = std::clamp(i, 1, ni - 2);
      const Vec3 d1 = 0.5 * (g.at(i + 1, j) - g.at(i - 1, j));
      const Vec3 d2 = g.at(i + 1, j) - 2.0 * g.at(i, j) + g.at(i - 1, j);
      const double l2 = dot(d1, d1);
      if (l2 == 0.0) {
        return 0.0;
      }
      return -dot(d1, d2) / l2;
    };
    auto along_j = [&](int i, int j) {
      j = std::clamp(j, 1, nj - 2);
      const Vec3 d1 = 0.5 * (g.at(i, j + 1) - g.at(i, j - 1));
      const Vec3 d2 = g.at(i, j + 1) - 2.0 * g.at(i, j) + g.at(i, j - 1);
      const double l2 = dot(d1, d1);
      if (l2 == 0.0) {
        return 0.0;
      }
      return -dot(d1, d2) / l2;
    };
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const double s = static_cast<double>(i) / (ni - 1);
        const double t = static_cast<double>(j) / (nj - 1);
        // P from the bottom / top rows, interpolated in eta; Q from left / right in xi.
        p[g.index(i, j)] = (1.0 - t) * along_i(i, 0) + t * along_i(i, nj - 1);
        q[g.index(i, j)] = (1.0 - s) * along_j(0, j) + s * along_j(ni - 1, j);
      }
    }
  }
  EllipticResult result;
  double first = 0.0;
  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    if (progress != nullptr && progress->cancelled()) {
      result.cancelled = true;
      break;
    }
    double sum = 0.0;
    std::size_t count = 0;
    for (int j = 1; j < nj - 1; ++j) {
      for (int i = 1; i < ni - 1; ++i) {
        const Vec3 xs = 0.5 * (g.at(i + 1, j) - g.at(i - 1, j));
        const Vec3 xt = 0.5 * (g.at(i, j + 1) - g.at(i, j - 1));
        const double alpha = dot(xt, xt);
        const double beta = dot(xs, xt);
        const double gamma = dot(xs, xs);
        const Vec3 xst = 0.25 * (g.at(i + 1, j + 1) - g.at(i + 1, j - 1) - g.at(i - 1, j + 1) +
                                 g.at(i - 1, j - 1));
        const std::size_t id = g.index(i, j);
        const double pv = p[id];
        const double qv = q[id];
        const double diag = 2.0 * (alpha + gamma);
        if (diag <= 0.0) {
          continue;
        }
        const Vec3 target =
            (alpha * ((1.0 + 0.5 * pv) * g.at(i + 1, j) + (1.0 - 0.5 * pv) * g.at(i - 1, j)) +
             gamma * ((1.0 + 0.5 * qv) * g.at(i, j + 1) + (1.0 - 0.5 * qv) * g.at(i, j - 1)) -
             2.0 * beta * xst) /
            diag;
        const Vec3 moved = g.at(i, j) + options.omega * (target - g.at(i, j));
        const SurfaceProjection pr = surface.project_from(moved, u[id], v[id]);
        u[id] = pr.u;
        v[id] = pr.v;
        const Vec3 delta = pr.point - g.at(i, j);
        g.at(i, j) = pr.point;
        sum += dot(delta, delta);
        ++count;
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
  if (!result.converged && !result.cancelled) {
    result.warnings.emplace_back("surface elliptic: not converged in " +
                                 std::to_string(options.max_iterations) + " iterations");
  }
  result.grid = std::move(g);
  return result;
}

}  // namespace og
