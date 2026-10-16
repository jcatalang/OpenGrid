#include "opengrid/core/structured/wall_layers.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace og {

namespace {

constexpr std::int64_t kNone = -1;
constexpr double kBlend = 0.1;  // weight of the neighbours' mean direction

struct Corner {
  std::size_t cell;
  unsigned corner;
};

std::string where(const Vec3& p, int dim) {
  std::ostringstream out;
  out << std::setprecision(4) << "(" << p.x << ", " << p.y;
  if (dim == 3) {
    out << ", " << p.z;
  }
  out << ")";
  return out.str();
}

Vec3 unit_or_zero(const Vec3& v) {
  const double length = norm(v);
  if (length > 0.0) {
    return (1.0 / length) * v;
  }
  return {};
}

class Complex {
 public:
  Complex(const WallLayerInput& input, std::size_t n) : input_(input), incident_(n) {
    const unsigned corners = corner_count();
    for (std::size_t c = 0; c < input.cells.size(); ++c) {
      for (unsigned k = 0; k < corners; ++k) {
        const std::int64_t g = input.cells[c][k];
        if (g < 0 || std::cmp_greater_equal(g, n)) {
          throw std::invalid_argument("wall layers: cell " + std::to_string(c) +
                                      " references a node out of range");
        }
        incident_[static_cast<std::size_t>(g)].push_back({c, k});
      }
    }
  }

  [[nodiscard]] unsigned corner_count() const {
    if (input_.dim == 3) {
      return 8U;
    }
    return 4U;
  }

  [[nodiscard]] std::int64_t node(const Corner& at, unsigned corner) const {
    return input_.cells[at.cell][corner];
  }

  // Corners edge-adjacent to `corner` in one cell.
  [[nodiscard]] std::vector<unsigned> adjacent(unsigned corner) const {
    if (input_.dim == 3) {
      return {corner ^ 1U, corner ^ 2U, corner ^ 4U};
    }
    return {(corner + 1U) % 4U, (corner + 3U) % 4U};
  }

  [[nodiscard]] const std::vector<Corner>& incident(std::int64_t g) const {
    return incident_[static_cast<std::size_t>(g)];
  }

  // Edge neighbours of g, sorted.
  [[nodiscard]] std::vector<std::int64_t> neighbours(std::int64_t g) const {
    std::vector<std::int64_t> out;
    for (const Corner& at : incident(g)) {
      for (const unsigned c : adjacent(at.corner)) {
        out.push_back(node(at, c));
      }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  // The corner of `other` in the cell of `at`, if it is edge-adjacent to at's corner.
  [[nodiscard]] int corner_of(const Corner& at, std::int64_t other) const {
    for (const unsigned c : adjacent(at.corner)) {
      if (node(at, c) == other) {
        return static_cast<int>(c);
      }
    }
    return -1;
  }

  // The node after `cur` on the grid line arriving from `prev`: the one edge neighbour of cur
  // that shares no cell with the edge (prev, cur). kNone if there is not exactly one.
  [[nodiscard]] std::int64_t continuation(std::int64_t prev, std::int64_t cur) const {
    std::vector<std::int64_t> beside;
    for (const Corner& at : incident(cur)) {
      if (corner_of(at, prev) < 0) {
        continue;
      }
      for (const unsigned c : adjacent(at.corner)) {
        beside.push_back(node(at, c));
      }
    }
    std::int64_t found = kNone;
    int count = 0;
    for (const std::int64_t q : neighbours(cur)) {
      if (q == prev || std::find(beside.begin(), beside.end(), q) != beside.end()) {
        continue;
      }
      found = q;
      count += 1;
    }
    if (count == 1) {
      return found;
    }
    return kNone;
  }

 private:
  const WallLayerInput& input_;
  std::vector<std::vector<Corner>> incident_;
};

}  // namespace

WallLayerResult march_wall_layers(const WallLayerInput& input, std::vector<Vec3>& x) {
  const std::size_t n = x.size();
  if (input.movable.size() != n || input.on_wall.size() != n ||
      (!input.plane.empty() && input.plane.size() != n)) {
    throw std::invalid_argument("wall layers: per-node arrays do not match the nodes");
  }
  const int dim = input.dim;
  const Complex complex(input, n);
  const std::vector<Vec3> x0 = x;
  auto free_node = [&](std::int64_t g) {
    return g >= 0 && input.movable[static_cast<std::size_t>(g)] != 0 &&
           input.on_wall[static_cast<std::size_t>(g)] == 0;
  };
  auto plane_of = [&](std::int64_t g) -> const WallLayerPlane* {
    if (input.plane.empty()) {
      return nullptr;
    }
    const std::int32_t p = input.plane[static_cast<std::size_t>(g)];
    if (p < 0 || std::cmp_greater_equal(p, input.planes.size())) {
      return nullptr;
    }
    return &input.planes[static_cast<std::size_t>(p)];
  };

  // ------------------------------------------------------------------- the lines
  WallLayerResult result;
  std::vector<int> layer_of(n, -1);  // 0: wall node, k: layer k
  for (std::size_t g = 0; g < n; ++g) {
    if (input.on_wall[g] != 0) {
      layer_of[g] = 0;
    }
  }
  std::vector<std::int64_t> owner(n, kNone);  // line of a layer node
  std::vector<double> spacing;
  std::vector<double> growth;
  for (const WallLayerNode& wall : input.walls) {
    if (wall.node < 0 || std::cmp_greater_equal(wall.node, n) || wall.layers <= 0) {
      continue;
    }
    if (!(wall.growth >= 1.0)) {
      throw std::invalid_argument("wall layers: growth must be at least 1");
    }
    for (const std::int64_t first : complex.neighbours(wall.node)) {
      if (!free_node(first)) {
        continue;
      }
      WallLayerLine line;
      line.base = wall.node;
      line.source = wall.source;
      std::int64_t prev = wall.node;
      std::int64_t cur = first;
      while (true) {
        const auto line_id = static_cast<std::int64_t>(result.lines.size());
        if (owner[static_cast<std::size_t>(cur)] != kNone &&
            owner[static_cast<std::size_t>(cur)] != line_id) {
          throw std::invalid_argument(
              "wall layers: the layers of two walls overlap at " +
              where(x[static_cast<std::size_t>(cur)], dim) +
              " (fewer grid lines between them than their layers): lower `layers`");
        }
        owner[static_cast<std::size_t>(cur)] = line_id;
        line.nodes.push_back(cur);
        if (std::cmp_greater_equal(line.nodes.size(), wall.layers)) {
          break;
        }
        const std::int64_t next = complex.continuation(prev, cur);
        if (!free_node(next)) {
          throw std::invalid_argument(
              "wall layers: the grid line leaving the wall at " +
              where(x[static_cast<std::size_t>(wall.node)], dim) + " ends after " +
              std::to_string(line.nodes.size()) + " of " + std::to_string(wall.layers) +
              " layers (at a block corner, a boundary or a wall): lower `layers`");
        }
        prev = cur;
        cur = next;
      }
      double s = wall.spacing;
      if (!(s > 0.0)) {
        s = distance(x[static_cast<std::size_t>(first)], x[static_cast<std::size_t>(wall.node)]);
      }
      if (!(s > 0.0)) {
        throw std::invalid_argument("wall layers: zero first height at " +
                                    where(x[static_cast<std::size_t>(wall.node)], dim));
      }
      line.next_spacing = s * std::pow(wall.growth, static_cast<double>(wall.layers));
      for (std::size_t k = 0; k < line.nodes.size(); ++k) {
        layer_of[static_cast<std::size_t>(line.nodes[k])] = static_cast<int>(k) + 1;
      }
      spacing.push_back(s);
      growth.push_back(wall.growth);
      result.lines.push_back(std::move(line));
    }
  }
  const std::size_t count = result.lines.size();
  if (count == 0) {
    return result;
  }
  std::size_t deepest = 0;
  for (const WallLayerLine& line : result.lines) {
    deepest = std::max(deepest, line.nodes.size());
  }
  auto node_at = [&](std::size_t l, std::size_t k) {  // k = 0: the base
    if (k == 0) {
      return result.lines[l].base;
    }
    return result.lines[l].nodes[k - 1];
  };

  // ------------------------------------------------------------------ marching
  // Facets at node p across the edge (p, next) in each cell around it, with the cell's corner
  // of p and the edge's axis bit (3D) or orientation (2D).
  auto facet_normal = [&](const Corner& at, int next_corner, Vec3& normal,
                          std::vector<std::int64_t>& facet, std::vector<std::int64_t>& above) {
    facet.clear();
    above.clear();
    const unsigned c = at.corner;
    const auto nc = static_cast<unsigned>(next_corner);
    auto pos = [&](unsigned k) { return x[static_cast<std::size_t>(complex.node(at, k))]; };
    if (dim == 2) {
      unsigned other = (c + 1U) % 4U;
      bool clockwise = false;  // rotate the facet clockwise to point at next
      if (nc == (c + 1U) % 4U) {
        other = (c + 3U) % 4U;
        clockwise = true;
      }
      const Vec3 f = pos(other) - pos(c);
      if (clockwise) {
        normal = {f.y, -f.x, 0.0};
      } else {
        normal = {-f.y, f.x, 0.0};
      }
      facet.push_back(complex.node(at, other));
      above.push_back(complex.node(at, (c + 2U) % 4U));
      return;
    }
    const unsigned edge = c ^ nc;  // the bit of the line's axis
    int a = 0;
    while ((1U << static_cast<unsigned>(a)) != edge) {
      ++a;
    }
    const unsigned b1 = 1U << static_cast<unsigned>((a + 1) % 3);
    const unsigned b2 = 1U << static_cast<unsigned>((a + 2) % 3);
    double sign = 1.0;
    for (const unsigned b : {edge, b1, b2}) {
      if ((c & b) != 0U) {
        sign = -sign;
      }
    }
    normal = sign * cross(pos(c ^ b1 ^ b2) - pos(c), pos(c ^ b2) - pos(c ^ b1));
    facet = {complex.node(at, c ^ b1), complex.node(at, c ^ b2), complex.node(at, c ^ b1 ^ b2)};
    above = {complex.node(at, nc ^ b1), complex.node(at, nc ^ b2)};
  };
  std::vector<Vec3> direction(count);
  std::vector<std::vector<std::size_t>> beside(count);
  std::vector<std::int64_t> facet;
  std::vector<std::int64_t> above;
  for (std::size_t k = 1; k <= deepest; ++k) {
    for (std::size_t l = 0; l < count; ++l) {
      beside[l].clear();
      if (result.lines[l].nodes.size() < k) {
        continue;
      }
      const std::int64_t p = node_at(l, k - 1);
      const std::int64_t next = node_at(l, k);
      Vec3 sum;
      for (const Corner& at : complex.incident(p)) {
        const int nc = complex.corner_of(at, next);
        if (nc < 0) {
          continue;
        }
        Vec3 normal;
        facet_normal(at, nc, normal, facet, above);
        bool front = true;
        for (const std::int64_t q : facet) {
          if (layer_of[static_cast<std::size_t>(q)] != static_cast<int>(k) - 1) {
            front = false;
          }
        }
        if (!front) {
          continue;
        }
        sum += unit_or_zero(normal);
        for (const std::int64_t q : above) {
          const std::int64_t m = owner[static_cast<std::size_t>(q)];
          if (m >= 0 && std::cmp_equal(layer_of[static_cast<std::size_t>(q)], k)) {
            beside[l].push_back(static_cast<std::size_t>(m));
          }
        }
      }
      Vec3 d = unit_or_zero(sum);
      if (norm(d) == 0.0) {
        d = unit_or_zero(x0[static_cast<std::size_t>(next)] - x0[static_cast<std::size_t>(p)]);
      }
      direction[l] = d;
      std::sort(beside[l].begin(), beside[l].end());
      beside[l].erase(std::unique(beside[l].begin(), beside[l].end()), beside[l].end());
    }
    for (std::size_t l = 0; l < count; ++l) {
      if (result.lines[l].nodes.size() < k) {
        continue;
      }
      Vec3 d = direction[l];
      if (!beside[l].empty()) {
        Vec3 mean;
        for (const std::size_t m : beside[l]) {
          mean += direction[m];
        }
        d = unit_or_zero((1.0 - kBlend) * d +
                         (kBlend / static_cast<double>(beside[l].size())) * mean);
      }
      const std::int64_t p = node_at(l, k - 1);
      const std::int64_t next = node_at(l, k);
      const WallLayerPlane* plane = plane_of(next);
      if (plane != nullptr) {
        d = unit_or_zero(d - dot(d, plane->normal) * plane->normal);
      }
      const double height = spacing[l] * std::pow(growth[l], static_cast<double>(k) - 1.0);
      Vec3 q = x[static_cast<std::size_t>(p)] + height * d;
      if (plane != nullptr) {
        q -= dot(q - plane->point, plane->normal) * plane->normal;
      }
      x[static_cast<std::size_t>(next)] = q;
    }
  }

  // ------------------------------------------------------------- folded cells
  const unsigned corners = complex.corner_count();
  for (const auto& cell : input.cells) {
    bool layer_cell = true;
    bool any_layer = false;
    for (unsigned k = 0; k < corners; ++k) {
      const int layer = layer_of[static_cast<std::size_t>(cell[k])];
      if (layer < 0) {
        layer_cell = false;
      }
      if (layer > 0) {
        any_layer = true;
      }
    }
    if (!layer_cell || !any_layer) {
      continue;
    }
    for (unsigned k = 0; k < corners; ++k) {
      const Vec3& p = x[static_cast<std::size_t>(cell[k])];
      double det = 0.0;
      if (dim == 2) {
        const Vec3 a = x[static_cast<std::size_t>(cell[(k + 1U) % 4U])] - p;
        const Vec3 b = x[static_cast<std::size_t>(cell[(k + 3U) % 4U])] - p;
        det = a.x * b.y - a.y * b.x;
      } else {
        const Vec3 a = x[static_cast<std::size_t>(cell[k ^ 1U])] - p;
        const Vec3 b = x[static_cast<std::size_t>(cell[k ^ 2U])] - p;
        const Vec3 e = x[static_cast<std::size_t>(cell[k ^ 4U])] - p;
        det = dot(cross(a, b), e);
        for (const unsigned bit : {1U, 2U, 4U}) {
          if ((k & bit) != 0U) {
            det = -det;
          }
        }
      }
      if (!(det > 0.0)) {
        int layer = 0;
        for (unsigned m = 0; m < corners; ++m) {
          layer = std::max(layer, layer_of[static_cast<std::size_t>(cell[m])]);
        }
        throw std::invalid_argument("wall layers: layer " + std::to_string(layer) + " folds at " +
                                    where(p, dim) +
                                    " (a concave corner or too thick layers): lower `layers` or "
                                    "`growth`");
      }
    }
  }

  // ------------------------------------------- the grid beyond the layers follows
  std::vector<Vec3> shift(n);
  for (std::size_t l = 0; l < count; ++l) {
    const WallLayerLine& line = result.lines[l];
    const std::int64_t last = line.nodes.back();
    std::int64_t prev = node_at(l, line.nodes.size() - 1);
    const Vec3 moved = x[static_cast<std::size_t>(last)] - x0[static_cast<std::size_t>(last)];
    std::vector<std::int64_t> beyond;
    std::vector<double> arc;
    double length = 0.0;
    std::int64_t cur = last;
    double step = 0.0;
    while (true) {
      const std::int64_t next = complex.continuation(prev, cur);
      if (next < 0) {
        length += step;  // no further node known: one more interval
        break;
      }
      step = distance(x0[static_cast<std::size_t>(next)], x0[static_cast<std::size_t>(cur)]);
      length += step;
      if (!free_node(next) || layer_of[static_cast<std::size_t>(next)] >= 0 || beyond.size() >= n) {
        break;
      }
      beyond.push_back(next);
      arc.push_back(length);
      prev = cur;
      cur = next;
    }
    if (!(length > 0.0)) {
      continue;
    }
    for (std::size_t j = 0; j < beyond.size(); ++j) {
      shift[static_cast<std::size_t>(beyond[j])] += (1.0 - arc[j] / length) * moved;
    }
  }
  for (std::size_t g = 0; g < n; ++g) {
    if (layer_of[g] >= 0 || norm(shift[g]) == 0.0) {
      continue;
    }
    x[g] += shift[g];
    const WallLayerPlane* plane = plane_of(static_cast<std::int64_t>(g));
    if (plane != nullptr) {
      x[g] -= dot(x[g] - plane->point, plane->normal) * plane->normal;
    }
  }
  return result;
}

}  // namespace og
