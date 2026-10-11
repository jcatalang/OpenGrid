#include "opengrid/core/quality/quality3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace og {
namespace {

constexpr double kDeg = 180.0 / std::numbers::pi;

double det3(const Vec3& a, const Vec3& b, const Vec3& c) {
  return dot(a, cross(b, c));
}

double angle_deg(const Vec3& a, const Vec3& b) {
  const double na = norm(a);
  const double nb = norm(b);
  if (na == 0.0 || nb == 0.0) {
    return 0.0;
  }
  return std::atan2(norm(cross(a, b)), dot(a, b)) * kDeg;
}

double acute_angle_deg(const Vec3& a, const Vec3& b) {
  const double na = norm(a);
  const double nb = norm(b);
  if (na == 0.0 || nb == 0.0) {
    return 90.0;
  }
  return std::acos(std::clamp(std::fabs(dot(a, b)) / (na * nb), 0.0, 1.0)) * kDeg;
}

// Linear weight of corner index `bit` at parameter x in [0, 1].
double weight(std::size_t bit, double x) {
  if (bit == 0) {
    return 1.0 - x;
  }
  return x;
}

// Cell corners c[a][b][c] for a, b, c in {0, 1}.
using Corners = std::array<std::array<std::array<Vec3, 2>, 2>, 2>;

// Volume of the trilinear map from the unit cube: 2-point Gauss in each direction is exact.
double trilinear_volume(const Corners& c) {
  const double g = 0.5 * std::numbers::inv_sqrt3;
  const std::array<double, 2> pts = {0.5 - g, 0.5 + g};
  double volume = 0.0;
  for (const double r : pts) {
    for (const double s : pts) {
      for (const double t : pts) {
        // Partial derivatives of x(r, s, t).
        Vec3 xr;
        Vec3 xs;
        Vec3 xt;
        for (std::size_t b = 0; b < 2; ++b) {
          for (std::size_t cc = 0; cc < 2; ++cc) {
            const double ws = weight(b, s);
            const double wt = weight(cc, t);
            xr += ws * wt * (c[1][b][cc] - c[0][b][cc]);
          }
        }
        for (std::size_t a = 0; a < 2; ++a) {
          for (std::size_t cc = 0; cc < 2; ++cc) {
            const double wr = weight(a, r);
            const double wt = weight(cc, t);
            xs += wr * wt * (c[a][1][cc] - c[a][0][cc]);
          }
        }
        for (std::size_t a = 0; a < 2; ++a) {
          for (std::size_t b = 0; b < 2; ++b) {
            const double wr = weight(a, r);
            const double ws = weight(b, s);
            xt += wr * ws * (c[a][b][1] - c[a][b][0]);
          }
        }
        volume += 0.125 * det3(xr, xs, xt);
      }
    }
  }
  return volume;
}

}  // namespace

Quality3D quality_3d(const Grid3D& g) {
  const int ci = g.ni - 1;
  const int cj = g.nj - 1;
  const int ck = g.nk - 1;
  const auto n = static_cast<std::size_t>(std::max(ci, 0)) *
                 static_cast<std::size_t>(std::max(cj, 0)) *
                 static_cast<std::size_t>(std::max(ck, 0));
  Quality3D q;
  for (auto* field : {&q.scaled_jacobian, &q.skewness, &q.aspect_ratio, &q.min_angle, &q.max_angle,
                      &q.non_orthogonality, &q.volume_ratio, &q.volume}) {
    field->assign(n, 0.0);
  }
  std::vector<Vec3> centroid(n);
  auto cell = [&](int i, int j, int k) {
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(cj) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(ci) +
           static_cast<std::size_t>(i);
  };
  auto corners = [&](int i, int j, int k) {
    Corners c{};
    for (std::size_t a = 0; a < 2; ++a) {
      for (std::size_t b = 0; b < 2; ++b) {
        for (std::size_t d = 0; d < 2; ++d) {
          c[a][b][d] =
              g.at(i + static_cast<int>(a), j + static_cast<int>(b), k + static_cast<int>(d));
        }
      }
    }
    return c;
  };
  for (int k = 0; k < ck; ++k) {
    for (int j = 0; j < cj; ++j) {
      for (int i = 0; i < ci; ++i) {
        const std::size_t id = cell(i, j, k);
        const Corners c = corners(i, j, k);
        double min_jac = std::numeric_limits<double>::infinity();
        std::array<double, 3> edge_sum = {0.0, 0.0, 0.0};
        Vec3 centre;
        for (std::size_t a = 0; a < 2; ++a) {
          for (std::size_t b = 0; b < 2; ++b) {
            for (std::size_t d = 0; d < 2; ++d) {
              const Vec3& o = c[a][b][d];
              centre += 0.125 * o;
              const Vec3 e1 = c[1 - a][b][d] - o;
              const Vec3 e2 = c[a][1 - b][d] - o;
              const Vec3 e3 = c[a][b][1 - d] - o;
              double sign = 1.0;
              if (a == 1) {
                sign = -sign;
              }
              if (b == 1) {
                sign = -sign;
              }
              if (d == 1) {
                sign = -sign;
              }
              const double l = norm(e1) * norm(e2) * norm(e3);
              double jac = 0.0;
              if (l > 0.0) {
                jac = sign * det3(e1, e2, e3) / l;
              }
              min_jac = std::min(min_jac, jac);
              // Edges counted from their lower corner only.
              if (a == 0) {
                edge_sum[0] += norm(e1);
              }
              if (b == 0) {
                edge_sum[1] += norm(e2);
              }
              if (d == 0) {
                edge_sum[2] += norm(e3);
              }
            }
          }
        }
        // Faces (corner loops) for angles and skewness.
        const std::array<std::array<Vec3, 4>, 6> faces = {{
            {c[0][0][0], c[0][1][0], c[0][1][1], c[0][0][1]},
            {c[1][0][0], c[1][1][0], c[1][1][1], c[1][0][1]},
            {c[0][0][0], c[1][0][0], c[1][0][1], c[0][0][1]},
            {c[0][1][0], c[1][1][0], c[1][1][1], c[0][1][1]},
            {c[0][0][0], c[1][0][0], c[1][1][0], c[0][1][0]},
            {c[0][0][1], c[1][0][1], c[1][1][1], c[0][1][1]},
        }};
        double min_angle = 360.0;
        double max_angle = 0.0;
        for (const auto& face : faces) {
          for (std::size_t m = 0; m < 4; ++m) {
            const double angle =
                angle_deg(face[(m + 1) % 4] - face[m], face[(m + 3) % 4] - face[m]);
            min_angle = std::min(min_angle, angle);
            max_angle = std::max(max_angle, angle);
          }
        }
        const double small = std::min({edge_sum[0], edge_sum[1], edge_sum[2]});
        const double large = std::max({edge_sum[0], edge_sum[1], edge_sum[2]});
        q.scaled_jacobian[id] = min_jac;
        q.min_angle[id] = min_angle;
        q.max_angle[id] = max_angle;
        q.skewness[id] = std::max((max_angle - 90.0) / 90.0, (90.0 - min_angle) / 90.0);
        if (small > 0.0) {
          q.aspect_ratio[id] = large / small;
        } else {
          q.aspect_ratio[id] = std::numeric_limits<double>::infinity();
        }
        q.volume[id] = trilinear_volume(c);
        centroid[id] = centre;
      }
    }
  }
  // Face neighbours: non-orthogonality and volume ratio.
  for (int k = 0; k < ck; ++k) {
    for (int j = 0; j < cj; ++j) {
      for (int i = 0; i < ci; ++i) {
        const std::size_t id = cell(i, j, k);
        const Corners c = corners(i, j, k);
        double worst_angle = 0.0;
        double worst_ratio = 1.0;
        // (face corner loop, neighbour offset)
        struct Face {
          std::array<Vec3, 4> p;
          int di;
          int dj;
          int dk;
        };
        const std::array<Face, 6> faces = {{
            {{c[0][0][0], c[0][1][0], c[0][1][1], c[0][0][1]}, -1, 0, 0},
            {{c[1][0][0], c[1][1][0], c[1][1][1], c[1][0][1]}, 1, 0, 0},
            {{c[0][0][0], c[1][0][0], c[1][0][1], c[0][0][1]}, 0, -1, 0},
            {{c[0][1][0], c[1][1][0], c[1][1][1], c[0][1][1]}, 0, 1, 0},
            {{c[0][0][0], c[1][0][0], c[1][1][0], c[0][1][0]}, 0, 0, -1},
            {{c[0][0][1], c[1][0][1], c[1][1][1], c[0][1][1]}, 0, 0, 1},
        }};
        for (const Face& face : faces) {
          const Vec3 area = 0.5 * cross(face.p[2] - face.p[0], face.p[3] - face.p[1]);
          const int ni2 = i + face.di;
          const int nj2 = j + face.dj;
          const int nk2 = k + face.dk;
          Vec3 d;
          if (ni2 >= 0 && ni2 < ci && nj2 >= 0 && nj2 < cj && nk2 >= 0 && nk2 < ck) {
            const std::size_t other = cell(ni2, nj2, nk2);
            d = centroid[other] - centroid[id];
            const double vo = std::fabs(q.volume[other]);
            const double vc = std::fabs(q.volume[id]);
            if (vo > 0.0 && vc > 0.0) {
              worst_ratio = std::max({worst_ratio, vo / vc, vc / vo});
            }
          } else {
            d = 0.25 * (face.p[0] + face.p[1] + face.p[2] + face.p[3]) - centroid[id];
          }
          worst_angle = std::max(worst_angle, acute_angle_deg(area, d));
        }
        q.non_orthogonality[id] = worst_angle;
        q.volume_ratio[id] = worst_ratio;
      }
    }
  }
  return q;
}

}  // namespace og
