#include "opengrid/core/structured/hyperbolic.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace og {
namespace {

// Neighbour index along a direction with the edge conditions: periodic wraps (the last line
// duplicates the first, so it is skipped), free / constant edges clamp (one-sided).
int step_index(int c, int d, int n, bool periodic) {
  int r = c + d;
  if (periodic) {
    // Lines 0 and n-1 coincide: the ring has n-1 distinct lines.
    const int m = n - 1;
    r = ((r % m) + m) % m;
    return r;
  }
  return std::clamp(r, 0, n - 1);
}

}  // namespace

HyperbolicResult extrude_hyperbolic(const Grid2D& input, double distance_total,
                                    const std::vector<double>& s,
                                    const HyperbolicOptions& options) {
  if (input.ni < 2 || input.nj < 2) {
    throw std::invalid_argument("hyperbolic extrusion: the base grid needs at least 2 x 2 points");
  }
  if (s.size() < 2 || s.front() != 0.0 || s.back() != 1.0) {
    throw std::invalid_argument("hyperbolic extrusion: layer fractions must run from 0 to 1");
  }
  for (std::size_t k = 1; k < s.size(); ++k) {
    if (!(s[k] > s[k - 1])) {
      throw std::invalid_argument("hyperbolic extrusion: layer fractions must increase");
    }
  }
  if (distance_total == 0.0) {
    throw std::invalid_argument("hyperbolic extrusion: zero distance");
  }
  HyperbolicResult result;
  Grid2D base = input;
  std::array<HyperbolicEdge, 4> edges = options.edges;
  if (distance_total < 0.0) {
    // March to the other side: mirror i (swapping the i edges).
    result.mirrored = true;
    for (int j = 0; j < input.nj; ++j) {
      for (int i = 0; i < input.ni; ++i) {
        base.at(i, j) = input.at(input.ni - 1 - i, j);
      }
    }
    std::swap(edges[0], edges[1]);
  }
  const double total = std::fabs(distance_total);
  const int ni = base.ni;
  const int nj = base.nj;
  const bool periodic_i =
      edges[0] == HyperbolicEdge::Periodic || edges[1] == HyperbolicEdge::Periodic;
  const bool periodic_j =
      edges[2] == HyperbolicEdge::Periodic || edges[3] == HyperbolicEdge::Periodic;
  double size = 0.0;
  for (const Vec3& p : base.x) {
    size = std::max(size, distance(p, base.x.front()));
  }
  const double tol = 1e-9 * std::max(size, 1e-300);
  if (periodic_i) {
    for (int j = 0; j < nj; ++j) {
      if (distance(base.at(0, j), base.at(ni - 1, j)) > tol) {
        throw std::invalid_argument(
            "hyperbolic extrusion: i is periodic but lines 0 and imax differ");
      }
    }
  }
  if (periodic_j) {
    for (int i = 0; i < ni; ++i) {
      if (distance(base.at(i, 0), base.at(i, nj - 1)) > tol) {
        throw std::invalid_argument(
            "hyperbolic extrusion: j is periodic but lines 0 and jmax differ");
      }
    }
  }
  const auto nk = static_cast<int>(s.size());
  result.grid = Grid3D(ni, nj, nk);
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      result.grid.at(i, j, 0) = base.at(i, j);
    }
  }
  Grid2D layer = base;
  std::vector<Vec3> normal(layer.x.size());
  std::vector<Vec3> scratch(layer.x.size());
  auto edge_of = [&](int i, int j) {
    // Edge condition applying to node (i, j) (corners: the i edge wins, then the j edge).
    if (!periodic_i && i == 0 && edges[0] != HyperbolicEdge::Free) {
      return edges[0];
    }
    if (!periodic_i && i == ni - 1 && edges[1] != HyperbolicEdge::Free) {
      return edges[1];
    }
    if (!periodic_j && j == 0 && edges[2] != HyperbolicEdge::Free) {
      return edges[2];
    }
    if (!periodic_j && j == nj - 1 && edges[3] != HyperbolicEdge::Free) {
      return edges[3];
    }
    return HyperbolicEdge::Free;
  };
  auto constrain = [&](int i, int j, Vec3 p) {
    const Vec3& origin = base.at(i, j);
    switch (edge_of(i, j)) {
      case HyperbolicEdge::ConstantX:
        p.x = origin.x;
        break;
      case HyperbolicEdge::ConstantY:
        p.y = origin.y;
        break;
      case HyperbolicEdge::ConstantZ:
        p.z = origin.z;
        break;
      default:
        break;
    }
    return p;
  };
  for (int k = 1; k < nk; ++k) {
    const double h = total * (s[static_cast<std::size_t>(k)] - s[static_cast<std::size_t>(k - 1)]);
    // Layer normals: central (one-sided at free edges) differences.
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const Vec3 ti = layer.at(step_index(i, 1, ni, periodic_i), j) -
                        layer.at(step_index(i, -1, ni, periodic_i), j);
        const Vec3 tj = layer.at(i, step_index(j, 1, nj, periodic_j)) -
                        layer.at(i, step_index(j, -1, nj, periodic_j));
        Vec3 nrm = cross(ti, tj);
        const double l = norm(nrm);
        if (l > 0.0) {
          nrm = nrm / l;
        }
        normal[layer.index(i, j)] = nrm;
      }
    }
    const std::vector<Vec3> raw = normal;
    for (int sweep = 0; sweep < options.normal_smoothing; ++sweep) {
      for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
          Vec3 sum = 2.0 * normal[layer.index(i, j)];
          sum += normal[layer.index(step_index(i, 1, ni, periodic_i), j)];
          sum += normal[layer.index(step_index(i, -1, ni, periodic_i), j)];
          sum += normal[layer.index(i, step_index(j, 1, nj, periodic_j))];
          sum += normal[layer.index(i, step_index(j, -1, nj, periodic_j))];
          const double l = norm(sum);
          if (l > 0.0) {
            sum = sum / l;
          }
          // Keep the direction visible from the raw normal (no folding back).
          if (dot(sum, raw[layer.index(i, j)]) < 0.5) {
            sum = raw[layer.index(i, j)];
          }
          scratch[layer.index(i, j)] = sum;
        }
      }
      normal = scratch;
    }
    // Step: the layer height, corrected for the angle to the raw normal.
    Grid2D next = layer;
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const std::size_t id = layer.index(i, j);
        const double c = std::max(dot(normal[id], raw[id]), 0.5);
        next.at(i, j) = constrain(i, j, layer.at(i, j) + (h / c) * normal[id]);
      }
    }
    // Dissipation: move towards the mean of the layer neighbours, tangentially.
    if (options.dissipation > 0.0) {
      Grid2D smooth = next;
      for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
          const bool interior_i = periodic_i || (i > 0 && i < ni - 1);
          const bool interior_j = periodic_j || (j > 0 && j < nj - 1);
          Vec3 mean;
          int count = 0;
          if (interior_i) {
            mean += next.at(step_index(i, 1, ni, periodic_i), j) +
                    next.at(step_index(i, -1, ni, periodic_i), j);
            count += 2;
          }
          if (interior_j) {
            mean += next.at(i, step_index(j, 1, nj, periodic_j)) +
                    next.at(i, step_index(j, -1, nj, periodic_j));
            count += 2;
          }
          if (count == 0) {
            continue;
          }
          const std::size_t id = layer.index(i, j);
          Vec3 move = options.dissipation * (mean / count - next.at(i, j));
          move -= dot(move, normal[id]) * normal[id];
          smooth.at(i, j) = constrain(i, j, next.at(i, j) + move);
        }
      }
      next = smooth;
    }
    if (periodic_i) {
      for (int j = 0; j < nj; ++j) {
        next.at(ni - 1, j) = next.at(0, j);
      }
    }
    if (periodic_j) {
      for (int i = 0; i < ni; ++i) {
        next.at(i, nj - 1) = next.at(i, 0);
      }
    }
    // Cell check: (x_i x x_j) . x_k > 0 at the cell centres.
    bool inverted = false;
    for (int j = 0; j + 1 < nj && !inverted; ++j) {
      for (int i = 0; i + 1 < ni && !inverted; ++i) {
        const Vec3 a = layer.at(i + 1, j + 1) - layer.at(i, j);
        const Vec3 b = layer.at(i, j + 1) - layer.at(i + 1, j);
        const Vec3 up = 0.25 * (next.at(i, j) + next.at(i + 1, j) + next.at(i, j + 1) +
                                next.at(i + 1, j + 1) - layer.at(i, j) - layer.at(i + 1, j) -
                                layer.at(i, j + 1) - layer.at(i + 1, j + 1));
        inverted = dot(cross(a, b), up) <= 0.0;
      }
    }
    if (inverted) {
      result.inverted_layers += 1;
    }
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        result.grid.at(i, j, k) = next.at(i, j);
      }
    }
    layer = next;
  }
  if (result.inverted_layers > 0) {
    result.warnings.push_back("hyperbolic extrusion: " + std::to_string(result.inverted_layers) +
                              " layers have inverted cells (reduce the distance or increase the "
                              "dissipation)");
  }
  return result;
}

}  // namespace og
