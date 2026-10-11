#include "opengrid/core/structured/tfi3d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace og {
namespace {

// Arc-length fractions along a polyline given by a node accessor.
template <class F>
std::vector<double> fractions(int n, F&& node) {
  std::vector<double> s(static_cast<std::size_t>(n), 0.0);
  for (int a = 1; a < n; ++a) {
    s[static_cast<std::size_t>(a)] =
        s[static_cast<std::size_t>(a - 1)] + distance(node(a), node(a - 1));
  }
  const double total = s.back();
  for (int a = 0; a < n; ++a) {
    if (total > 0.0) {
      s[static_cast<std::size_t>(a)] /= total;
    } else {
      s[static_cast<std::size_t>(a)] = static_cast<double>(a) / (n - 1);
    }
  }
  return s;
}

void require(bool ok, const std::string& what) {
  if (!ok) {
    throw std::invalid_argument("block TFI: " + what);
  }
}

// Linear TFI in index space of a scalar f(a, b) on [0, na-1] x [0, nb-1] from its boundary
// values: f0b = f(0, b), f1b = f(na-1, b), fa0 = f(a, 0), fa1 = f(a, nb-1).
double blend(double f0b, double f1b, double fa0, double fa1, double c00, double c10, double c01,
             double c11, double s, double t) {
  return (1.0 - s) * f0b + s * f1b + (1.0 - t) * fa0 + t * fa1 -
         ((1.0 - s) * (1.0 - t) * c00 + s * (1.0 - t) * c10 + (1.0 - s) * t * c01 + s * t * c11);
}

}  // namespace

Grid3D tfi_block(const BlockFaces& f) {
  const int ni = f.kmin.ni;
  const int nj = f.kmin.nj;
  const int nk = f.imin.nj;
  require(ni >= 2 && nj >= 2 && nk >= 2, "every direction needs at least 2 points");
  require(f.kmax.ni == ni && f.kmax.nj == nj, "kmin / kmax sizes differ");
  require(f.jmin.ni == ni && f.jmax.ni == ni && f.jmin.nj == nk && f.jmax.nj == nk,
          "jmin / jmax sizes do not match (ni, nk)");
  require(f.imin.ni == nj && f.imax.ni == nj && f.imax.nj == nk,
          "imin / imax sizes do not match (nj, nk)");
  // Shared edges must coincide.
  double size = 0.0;
  for (const auto* face : {&f.imin, &f.imax, &f.jmin, &f.jmax, &f.kmin, &f.kmax}) {
    for (const Vec3& p : face->x) {
      size = std::max(size, norm(p - f.kmin.x.front()));
    }
  }
  const double tol = 1e-9 * std::max(size, 1e-300);
  auto same = [&](const Vec3& a, const Vec3& b, const char* edge) {
    require(distance(a, b) <= tol, std::string("faces disagree on the ") + edge + " edge");
  };
  for (int i = 0; i < ni; ++i) {
    same(f.kmin.at(i, 0), f.jmin.at(i, 0), "(j=0, k=0)");
    same(f.kmin.at(i, nj - 1), f.jmax.at(i, 0), "(j=max, k=0)");
    same(f.kmax.at(i, 0), f.jmin.at(i, nk - 1), "(j=0, k=max)");
    same(f.kmax.at(i, nj - 1), f.jmax.at(i, nk - 1), "(j=max, k=max)");
  }
  for (int j = 0; j < nj; ++j) {
    same(f.kmin.at(0, j), f.imin.at(j, 0), "(i=0, k=0)");
    same(f.kmin.at(ni - 1, j), f.imax.at(j, 0), "(i=max, k=0)");
    same(f.kmax.at(0, j), f.imin.at(j, nk - 1), "(i=0, k=max)");
    same(f.kmax.at(ni - 1, j), f.imax.at(j, nk - 1), "(i=max, k=max)");
  }
  for (int k = 0; k < nk; ++k) {
    same(f.jmin.at(0, k), f.imin.at(0, k), "(i=0, j=0)");
    same(f.jmin.at(ni - 1, k), f.imax.at(0, k), "(i=max, j=0)");
    same(f.jmax.at(0, k), f.imin.at(nj - 1, k), "(i=0, j=max)");
    same(f.jmax.at(ni - 1, k), f.imax.at(nj - 1, k), "(i=max, j=max)");
  }
  // Arc-length fractions of the grid lines on the faces.
  // xi on the four faces containing i-lines: jmin(i, k), jmax(i, k), kmin(i, j), kmax(i, j).
  auto lines = [](const Grid2D& g, bool along_first) {
    std::vector<std::vector<double>> out;
    if (along_first) {
      for (int b = 0; b < g.nj; ++b) {
        out.push_back(fractions(g.ni, [&](int a) { return g.at(a, b); }));
      }
    } else {
      for (int a = 0; a < g.ni; ++a) {
        out.push_back(fractions(g.nj, [&](int b) { return g.at(a, b); }));
      }
    }
    return out;
  };
  const auto xi_jmin = lines(f.jmin, true);     // [k][i]
  const auto xi_jmax = lines(f.jmax, true);     // [k][i]
  const auto xi_kmin = lines(f.kmin, true);     // [j][i]
  const auto xi_kmax = lines(f.kmax, true);     // [j][i]
  const auto eta_imin = lines(f.imin, true);    // [k][j]
  const auto eta_imax = lines(f.imax, true);    // [k][j]
  const auto eta_kmin = lines(f.kmin, false);   // [i][j]
  const auto eta_kmax = lines(f.kmax, false);   // [i][j]
  const auto zeta_imin = lines(f.imin, false);  // [j][k]
  const auto zeta_imax = lines(f.imax, false);  // [j][k]
  const auto zeta_jmin = lines(f.jmin, false);  // [i][k]
  const auto zeta_jmax = lines(f.jmax, false);  // [i][k]
  auto at = [](const std::vector<std::vector<double>>& v, int a, int b) {
    return v[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
  };

  Grid3D g(ni, nj, nk);
  for (int k = 0; k < nk; ++k) {
    const double zb = static_cast<double>(k) / (nk - 1);
    for (int j = 0; j < nj; ++j) {
      const double yb = static_cast<double>(j) / (nj - 1);
      for (int i = 0; i < ni; ++i) {
        const double xb = static_cast<double>(i) / (ni - 1);
        // Blending parameters: TFI in index space of the face fractions.
        const double xi = blend(at(xi_jmin, k, i), at(xi_jmax, k, i), at(xi_kmin, j, i),
                                at(xi_kmax, j, i), at(xi_jmin, 0, i), at(xi_jmax, 0, i),
                                at(xi_jmin, nk - 1, i), at(xi_jmax, nk - 1, i), yb, zb);
        const double eta = blend(at(eta_imin, k, j), at(eta_imax, k, j), at(eta_kmin, i, j),
                                 at(eta_kmax, i, j), at(eta_imin, 0, j), at(eta_imax, 0, j),
                                 at(eta_imin, nk - 1, j), at(eta_imax, nk - 1, j), xb, zb);
        const double zeta = blend(at(zeta_imin, j, k), at(zeta_imax, j, k), at(zeta_jmin, i, k),
                                  at(zeta_jmax, i, k), at(zeta_imin, 0, k), at(zeta_imax, 0, k),
                                  at(zeta_imin, nj - 1, k), at(zeta_imax, nj - 1, k), xb, yb);
        const double p0 = 1.0 - xi;
        const double p1 = xi;
        const double q0 = 1.0 - eta;
        const double q1 = eta;
        const double r0 = 1.0 - zeta;
        const double r1 = zeta;
        const Vec3 u = p0 * f.imin.at(j, k) + p1 * f.imax.at(j, k);
        const Vec3 v = q0 * f.jmin.at(i, k) + q1 * f.jmax.at(i, k);
        const Vec3 w = r0 * f.kmin.at(i, j) + r1 * f.kmax.at(i, j);
        // Edges: k-direction at (i, j) corners (from imin / imax), i-direction at (j, k) corners
        // (from jmin / jmax), j-direction at (i, k) corners (from imin / imax).
        const Vec3 uv = p0 * (q0 * f.imin.at(0, k) + q1 * f.imin.at(nj - 1, k)) +
                        p1 * (q0 * f.imax.at(0, k) + q1 * f.imax.at(nj - 1, k));
        const Vec3 vw = q0 * (r0 * f.jmin.at(i, 0) + r1 * f.jmin.at(i, nk - 1)) +
                        q1 * (r0 * f.jmax.at(i, 0) + r1 * f.jmax.at(i, nk - 1));
        const Vec3 wu = r0 * (p0 * f.imin.at(j, 0) + p1 * f.imax.at(j, 0)) +
                        r1 * (p0 * f.imin.at(j, nk - 1) + p1 * f.imax.at(j, nk - 1));
        const Vec3 uvw =
            r0 * (q0 * (p0 * f.imin.at(0, 0) + p1 * f.imax.at(0, 0)) +
                  q1 * (p0 * f.imin.at(nj - 1, 0) + p1 * f.imax.at(nj - 1, 0))) +
            r1 * (q0 * (p0 * f.imin.at(0, nk - 1) + p1 * f.imax.at(0, nk - 1)) +
                  q1 * (p0 * f.imin.at(nj - 1, nk - 1) + p1 * f.imax.at(nj - 1, nk - 1)));
        g.at(i, j, k) = u + v + w - uv - vw - wu + uvw;
      }
    }
  }
  // Faces exactly (the formula reproduces them up to rounding).
  const BlockFaces& b = f;
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      g.at(0, j, k) = b.imin.at(j, k);
      g.at(ni - 1, j, k) = b.imax.at(j, k);
    }
    for (int i = 0; i < ni; ++i) {
      g.at(i, 0, k) = b.jmin.at(i, k);
      g.at(i, nj - 1, k) = b.jmax.at(i, k);
    }
  }
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      g.at(i, j, 0) = b.kmin.at(i, j);
      g.at(i, j, nk - 1) = b.kmax.at(i, j);
    }
  }
  return g;
}

BlockFaces block_faces(const Grid3D& g) {
  BlockFaces f;
  f.imin = Grid2D(g.nj, g.nk);
  f.imax = Grid2D(g.nj, g.nk);
  f.jmin = Grid2D(g.ni, g.nk);
  f.jmax = Grid2D(g.ni, g.nk);
  f.kmin = Grid2D(g.ni, g.nj);
  f.kmax = Grid2D(g.ni, g.nj);
  for (int k = 0; k < g.nk; ++k) {
    for (int j = 0; j < g.nj; ++j) {
      f.imin.at(j, k) = g.at(0, j, k);
      f.imax.at(j, k) = g.at(g.ni - 1, j, k);
    }
    for (int i = 0; i < g.ni; ++i) {
      f.jmin.at(i, k) = g.at(i, 0, k);
      f.jmax.at(i, k) = g.at(i, g.nj - 1, k);
    }
  }
  for (int j = 0; j < g.nj; ++j) {
    for (int i = 0; i < g.ni; ++i) {
      f.kmin.at(i, j) = g.at(i, j, 0);
      f.kmax.at(i, j) = g.at(i, j, g.nk - 1);
    }
  }
  return f;
}

}  // namespace og
