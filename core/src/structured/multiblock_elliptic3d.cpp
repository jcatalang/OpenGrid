#include "opengrid/core/structured/multiblock_elliptic3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>

#include "opengrid/core/util/anderson.hpp"
#include "opengrid/core/util/parallel.hpp"

namespace og {
namespace {

constexpr std::int64_t kNone = -1;
constexpr int kCentre = 13;
constexpr std::array<int, 3> kStride = {1, 3, 9};

// Face directions 0..5: axis d / 2, sign + for even d.
int dir_axis(int d) {
  return d / 2;
}
int dir_sign(int d) {
  if (d % 2 == 0) {
    return 1;
  }
  return -1;
}
int dir_opposite(int d) {
  return static_cast<int>(static_cast<unsigned>(d) ^ 1U);
}
// Bit of corner index c = i + 2 j + 4 k for local axis a.
unsigned bit(int a) {
  return 1U << static_cast<unsigned>(a);
}
int dir_of(int axis, int sign) {
  if (sign > 0) {
    return 2 * axis;
  }
  return 2 * axis + 1;
}
// Stencil slot of a frame offset (each component -1, 0, 1).
int slot_of(const std::array<int, 3>& offset) {
  return kCentre + offset[0] * kStride[0] + offset[1] * kStride[1] + offset[2] * kStride[2];
}
int slot_of_dir(int d) {
  return kCentre + dir_sign(d) * kStride[static_cast<std::size_t>(dir_axis(d))];
}

using Stencil27 = std::array<std::int64_t, 27>;

enum class Shape : std::uint8_t { None, Interior, Face, Edge };

struct NodeTopo {
  Stencil27 nb{};
  Shape shape = Shape::None;
  int inward = -1;                   // Face: direction into the domain; Edge: unused
  int edge_axis = -1;                // Edge: frame axis along the boundary edge
  std::array<int, 2> inward_axes{};  // Edge: the two axes pointing inside, with signs below
  std::array<int, 2> inward_signs{};
  int hexes = 0;
  bool closed = false;     // every face around the node shared by two hexahedra
  bool structural = true;  // one fan of hexahedra, no face shared by three
  bool consistent = true;  // the frames of all hexahedra agree (a stencil exists)
};

struct Incident {
  std::int64_t hex;
  unsigned corner;
};

// Walks the hexahedra around node g across shared faces and assembles the 3 x 3 x 3 stencil.
NodeTopo node_topology(std::int64_t g, const std::vector<Incident>& incident,
                       const std::vector<std::array<std::int64_t, 8>>& hexes) {
  NodeTopo t;
  t.nb.fill(kNone);
  const std::size_t k = incident.size();
  t.hexes = static_cast<int>(k);
  if (k == 0) {
    t.structural = false;
    t.consistent = false;
    return t;
  }
  // Start: the hexahedron in which g is corner 0 (lowest index), else the lowest index.
  std::size_t start = k;
  for (std::size_t e = 0; e < k; ++e) {
    if (incident[e].corner == 0U && (start == k || incident[e].hex < incident[start].hex)) {
      start = e;
    }
  }
  if (start == k) {
    start = 0;
    for (std::size_t e = 1; e < k; ++e) {
      if (incident[e].hex < incident[start].hex) {
        start = e;
      }
    }
  }
  struct Map {
    std::array<int, 3> axis{};
    std::array<int, 3> sign{};
  };
  std::vector<Map> maps(k);
  std::vector<std::uint8_t> seen(k, 0);
  maps[start] = {{0, 1, 2}, {1, 1, 1}};
  seen[start] = 1;
  std::vector<std::size_t> queue{start};
  bool boundary = false;
  auto corner_node = [&](std::size_t e, unsigned bits) {
    return hexes[static_cast<std::size_t>(incident[e].hex)][static_cast<std::size_t>(bits)];
  };
  for (std::size_t q = 0; q < queue.size(); ++q) {
    const std::size_t e = queue[q];
    const unsigned b = incident[e].corner;
    const Map& m = maps[e];
    for (int a = 0; a < 3; ++a) {
      const int a1 = (a + 1) % 3;
      const int a2 = (a + 2) % 3;
      const std::int64_t n1 = corner_node(e, b ^ bit(a1));
      const std::int64_t n2 = corner_node(e, b ^ bit(a2));
      const std::int64_t n12 = corner_node(e, b ^ bit(a1) ^ bit(a2));
      std::size_t other = k;
      for (std::size_t f = 0; f < k; ++f) {
        if (f == e) {
          continue;
        }
        const auto& h = hexes[static_cast<std::size_t>(incident[f].hex)];
        const bool has1 = std::find(h.begin(), h.end(), n1) != h.end();
        const bool has2 = std::find(h.begin(), h.end(), n2) != h.end();
        const bool has12 = std::find(h.begin(), h.end(), n12) != h.end();
        if (has1 && has2 && has12) {
          if (other != k) {
            t.structural = false;  // three hexahedra on one face
          }
          other = f;
        }
      }
      if (other == k) {
        boundary = true;
        continue;
      }
      Map next;
      int fresh = 0;
      const unsigned bo = incident[other].corner;
      for (int ao = 0; ao < 3; ++ao) {
        const std::int64_t node = corner_node(other, bo ^ bit(ao));
        if (node == n1) {
          next.axis[static_cast<std::size_t>(ao)] = m.axis[static_cast<std::size_t>(a1)];
          next.sign[static_cast<std::size_t>(ao)] = m.sign[static_cast<std::size_t>(a1)];
        } else if (node == n2) {
          next.axis[static_cast<std::size_t>(ao)] = m.axis[static_cast<std::size_t>(a2)];
          next.sign[static_cast<std::size_t>(ao)] = m.sign[static_cast<std::size_t>(a2)];
        } else {
          next.axis[static_cast<std::size_t>(ao)] = m.axis[static_cast<std::size_t>(a)];
          next.sign[static_cast<std::size_t>(ao)] = -m.sign[static_cast<std::size_t>(a)];
          fresh += 1;
        }
      }
      if (fresh != 1) {
        t.structural = false;
        continue;
      }
      if (seen[other] != 0) {
        if (maps[other].axis != next.axis || maps[other].sign != next.sign) {
          t.consistent = false;
        }
        continue;
      }
      seen[other] = 1;
      maps[other] = next;
      queue.push_back(other);
    }
  }
  if (queue.size() != k) {
    t.structural = false;  // several fans meet at the node
  }
  t.closed = !boundary;
  if (!t.consistent || !t.structural) {
    // Irregular (frames rotate around an edge of odd valence) or malformed: no stencil.
    t.nb.fill(kNone);
    return t;
  }
  // Octants and the stencil.
  std::array<int, 3> lo = {1, 1, 1};
  std::array<int, 3> hi = {-1, -1, -1};
  std::vector<std::array<int, 3>> octants;
  for (std::size_t e = 0; e < k; ++e) {
    const Map& m = maps[e];
    std::array<int, 3> octant{};
    for (int a = 0; a < 3; ++a) {
      octant[static_cast<std::size_t>(m.axis[static_cast<std::size_t>(a)])] =
          m.sign[static_cast<std::size_t>(a)];
    }
    for (int a = 0; a < 3; ++a) {
      lo[static_cast<std::size_t>(a)] =
          std::min(lo[static_cast<std::size_t>(a)], octant[static_cast<std::size_t>(a)]);
      hi[static_cast<std::size_t>(a)] =
          std::max(hi[static_cast<std::size_t>(a)], octant[static_cast<std::size_t>(a)]);
    }
    if (std::find(octants.begin(), octants.end(), octant) != octants.end()) {
      t.consistent = false;
      t.nb.fill(kNone);
      return t;
    }
    octants.push_back(octant);
    const unsigned b = incident[e].corner;
    for (unsigned c = 0; c < 8; ++c) {
      const unsigned d = c ^ b;
      std::array<int, 3> offset{};
      for (int a = 0; a < 3; ++a) {
        if ((d & bit(a)) != 0U) {
          offset[static_cast<std::size_t>(m.axis[static_cast<std::size_t>(a)])] =
              m.sign[static_cast<std::size_t>(a)];
        }
      }
      const auto s = static_cast<std::size_t>(slot_of(offset));
      const std::int64_t node = corner_node(e, c);
      if (t.nb[s] != kNone && t.nb[s] != node) {
        t.consistent = false;
        t.nb.fill(kNone);
        return t;
      }
      t.nb[s] = node;
    }
  }
  (void)g;
  if (t.closed && k == 8) {
    t.shape = Shape::Interior;
  } else if (!t.closed && k == 4) {
    // A half space: one axis with a single sign.
    int fixed_axes = 0;
    for (int a = 0; a < 3; ++a) {
      if (lo[static_cast<std::size_t>(a)] == hi[static_cast<std::size_t>(a)]) {
        fixed_axes += 1;
        t.inward = dir_of(a, lo[static_cast<std::size_t>(a)]);
      }
    }
    if (fixed_axes == 1) {
      t.shape = Shape::Face;
    }
  } else if (!t.closed && k == 2) {
    int varying = 0;
    int m = 0;
    for (int a = 0; a < 3; ++a) {
      if (lo[static_cast<std::size_t>(a)] != hi[static_cast<std::size_t>(a)]) {
        varying += 1;
        t.edge_axis = a;
      } else if (m < 2) {
        t.inward_axes[static_cast<std::size_t>(m)] = a;
        t.inward_signs[static_cast<std::size_t>(m)] = lo[static_cast<std::size_t>(a)];
        ++m;
      }
    }
    if (varying == 1) {
      t.shape = Shape::Edge;
    }
  }
  return t;
}

int dir_to(const Stencil27& nb, std::int64_t node) {
  for (int d = 0; d < 6; ++d) {
    if (nb[static_cast<std::size_t>(slot_of_dir(d))] == node) {
      return d;
    }
  }
  return -1;
}

bool solve_line3(const double* la, double* lb, const double* lc, std::array<double*, 3> r,
                 std::size_t n) {
  const std::size_t last = n - 2;
  for (double* x : r) {
    x[1] -= la[1] * x[0];
    x[last] -= lc[last] * x[n - 1];
  }
  for (std::size_t k = 2; k <= last; ++k) {
    if (!(lb[k - 1] != 0.0)) {
      return false;
    }
    const double w = la[k] / lb[k - 1];
    lb[k] -= w * lc[k - 1];
    for (double* x : r) {
      x[k] -= w * x[k - 1];
    }
  }
  if (!(lb[last] != 0.0)) {
    return false;
  }
  for (double* x : r) {
    x[last] /= lb[last];
    for (std::size_t k = last - 1; k >= 1; --k) {
      x[k] = (x[k] - lc[k] * x[k + 1]) / lb[k];
    }
  }
  return true;
}

struct Line {
  std::size_t offset;
  std::size_t size;
  std::int64_t end0;
  std::int64_t end1;
};

struct WallReach {
  std::size_t wall;
  double weight;
  int normal;
  std::array<int, 2> tangent;
};

}  // namespace

Multiblock3DResult smooth_multiblock_3d(const Multiblock3DInput& input,
                                        const MultiblockOptions& options,
                                        SmoothingProgress* progress) {
  const std::size_t n = input.xyz.size() / 3;
  if (input.xyz.size() != 3 * n || input.fixed.size() != n) {
    throw std::invalid_argument("multiblock 3D: xyz and fixed sizes do not match");
  }
  std::vector<std::vector<Incident>> incident(n);
  for (std::size_t h = 0; h < input.hexes.size(); ++h) {
    for (unsigned c = 0; c < 8; ++c) {
      const std::int64_t g = input.hexes[h][c];
      if (g < 0 || std::cmp_greater_equal(g, n)) {
        throw std::invalid_argument("multiblock 3D: hexahedron " + std::to_string(h) +
                                    " references a node out of range");
      }
      incident[static_cast<std::size_t>(g)].push_back({static_cast<std::int64_t>(h), c});
    }
  }
  std::vector<NodeTopo> topo(n);
  parallel_for(
      0, static_cast<std::int64_t>(n),
      [&](std::int64_t g) {
        topo[static_cast<std::size_t>(g)] =
            node_topology(g, incident[static_cast<std::size_t>(g)], input.hexes);
      },
      256);

  Multiblock3DResult result;
  auto& cls = result.node_class;
  cls.assign(n, static_cast<std::uint8_t>(MultiblockNode::Fixed));
  std::vector<std::int64_t> slide_of(n, -1);
  for (std::size_t s = 0; s < input.slides.size(); ++s) {
    const Multiblock3DSlide& slide = input.slides[s];
    bool bad = slide.node < 0 || std::cmp_greater_equal(slide.node, n);
    if (slide.target == SlideTarget::Surface) {
      bad = bad || slide.index >= input.surfaces.size() || !input.surfaces[slide.index];
    } else {
      bad = bad || slide.index >= input.chains.size() || !input.chains[slide.index];
    }
    if (bad) {
      throw std::invalid_argument("multiblock 3D: slide " + std::to_string(s) + " is out of range");
    }
    slide_of[static_cast<std::size_t>(slide.node)] = static_cast<std::int64_t>(s);
  }
  int kept = 0;
  int bad_slides = 0;
  for (std::size_t g = 0; g < n; ++g) {
    const NodeTopo& t = topo[g];
    if (input.fixed[g] != 0) {
      continue;
    }
    if (slide_of[g] >= 0) {
      const Multiblock3DSlide& slide = input.slides[static_cast<std::size_t>(slide_of[g])];
      // On a surface: inside a boundary face. Along a curve: on a boundary edge, or on a
      // smooth seam between two faces (four hexahedra, e.g. two patches of one sphere).
      const bool face = slide.target == SlideTarget::Surface && t.shape == Shape::Face;
      const bool edge =
          slide.target == SlideTarget::Chain && (t.shape == Shape::Edge || t.shape == Shape::Face);
      if (t.consistent && (face || edge)) {
        cls[g] = static_cast<std::uint8_t>(MultiblockNode::Sliding);
      } else {
        bad_slides += 1;
      }
      continue;
    }
    if (!t.structural || !t.closed) {
      kept += 1;
      continue;
    }
    if (t.consistent && t.shape == Shape::Interior) {
      cls[g] = static_cast<std::uint8_t>(MultiblockNode::Regular);
    } else if (t.hexes >= 3) {
      cls[g] = static_cast<std::uint8_t>(MultiblockNode::Irregular);
    } else {
      kept += 1;
    }
  }
  if (kept > 0) {
    result.warnings.push_back("multiblock 3D: " + std::to_string(kept) +
                              " free node(s) on the boundary or without a regular fan of "
                              "hexahedra kept fixed");
  }
  if (bad_slides > 0) {
    result.warnings.push_back("multiblock 3D: " + std::to_string(bad_slides) +
                              " sliding node(s) not inside a boundary face or edge kept fixed");
  }
  auto is = [&](std::int64_t g, MultiblockNode c) {
    return g >= 0 && cls[static_cast<std::size_t>(g)] == static_cast<std::uint8_t>(c);
  };
  auto nb_dir = [&](std::size_t g, int d) {
    return topo[g].nb[static_cast<std::size_t>(slot_of_dir(d))];
  };

  std::vector<Vec3> x(n);
  for (std::size_t g = 0; g < n; ++g) {
    x[g] = {input.xyz[3 * g], input.xyz[3 * g + 1], input.xyz[3 * g + 2]};
  }
  double size = 0.0;
  for (const Vec3& p : x) {
    size = std::max(size, distance(p, x.front()));
  }

  // ------------------------------------------------------------------- grid lines
  std::vector<std::int64_t> slot(3 * n, -1);  // line through (node, axis); -2: opening node
  std::vector<std::int64_t> line_nodes;
  std::vector<int> line_prev;
  std::vector<int> line_next;
  std::vector<std::pair<int, std::int64_t>> line_side;  // 4 per line node
  std::vector<Line> lines;
  for (std::size_t s = 0; s < n; ++s) {
    if (cls[s] != static_cast<std::uint8_t>(MultiblockNode::Regular)) {
      continue;
    }
    for (int axis = 0; axis < 3; ++axis) {
      if (slot[3 * s + static_cast<std::size_t>(axis)] != -1) {
        continue;
      }
      const auto line_id = static_cast<std::int64_t>(lines.size());
      auto walk = [&](int dir, std::vector<std::int64_t>& nodes, bool& closed) {
        auto cur = static_cast<std::int64_t>(s);
        int d = dir;
        closed = false;
        while (true) {
          const std::int64_t next = nb_dir(static_cast<std::size_t>(cur), d);
          if (std::cmp_equal(next, s)) {
            closed = true;
            return next;
          }
          if (!is(next, MultiblockNode::Regular)) {
            return next;
          }
          const int back = dir_to(topo[static_cast<std::size_t>(next)].nb, cur);
          if (back < 0) {
            return next;
          }
          const auto next_slot =
              3 * static_cast<std::size_t>(next) + static_cast<std::size_t>(dir_axis(back));
          if (slot[next_slot] != -1) {
            return next;
          }
          nodes.push_back(next);
          slot[next_slot] = line_id;
          cur = next;
          d = dir_opposite(back);
        }
      };
      slot[3 * s + static_cast<std::size_t>(axis)] = line_id;
      std::vector<std::int64_t> forward;
      bool closed = false;
      const std::int64_t end1 = walk(dir_of(axis, 1), forward, closed);
      std::vector<std::int64_t> ordered;
      std::int64_t end0 = kNone;
      if (closed) {
        slot[3 * s + static_cast<std::size_t>(axis)] = -2;
        ordered = forward;
        end0 = static_cast<std::int64_t>(s);
      } else {
        std::vector<std::int64_t> backward;
        bool unused = false;
        end0 = walk(dir_of(axis, -1), backward, unused);
        ordered.assign(backward.rbegin(), backward.rend());
        ordered.push_back(static_cast<std::int64_t>(s));
        ordered.insert(ordered.end(), forward.begin(), forward.end());
      }
      if (ordered.empty()) {
        continue;
      }
      lines.push_back({line_nodes.size(), ordered.size(), end0, end1});
      for (std::size_t m = 0; m < ordered.size(); ++m) {
        std::int64_t before = end0;
        std::int64_t after = end1;
        if (m > 0) {
          before = ordered[m - 1];
        }
        if (m + 1 < ordered.size()) {
          after = ordered[m + 1];
        }
        const auto& nb = topo[static_cast<std::size_t>(ordered[m])].nb;
        line_nodes.push_back(ordered[m]);
        const int prev = dir_to(nb, before);
        const int next = dir_to(nb, after);
        line_prev.push_back(prev);
        line_next.push_back(next);
        // The four neighbours off the line (direction, node), stored in line order so a sweep
        // reads them sequentially.
        for (int d = 0; d < 6; ++d) {
          if (d != prev && d != next) {
            line_side.emplace_back(d, nb[static_cast<std::size_t>(slot_of_dir(d))]);
          }
        }
      }
    }
  }
  std::vector<int> colour(lines.size(), -1);
  int colours = 0;
  for (std::size_t l = 0; l < lines.size(); ++l) {
    std::vector<std::uint8_t> used(static_cast<std::size_t>(colours) + 1, 0);
    const Line& line = lines[l];
    for (std::size_t m = 0; m < line.size; ++m) {
      const std::int64_t node = line_nodes[line.offset + m];
      for (const std::int64_t s : topo[static_cast<std::size_t>(node)].nb) {
        if (s < 0) {
          continue;
        }
        for (int axis = 0; axis < 3; ++axis) {
          const std::int64_t other =
              slot[3 * static_cast<std::size_t>(s) + static_cast<std::size_t>(axis)];
          if (other >= 0 && std::cmp_not_equal(other, l) &&
              colour[static_cast<std::size_t>(other)] >= 0) {
            used[static_cast<std::size_t>(colour[static_cast<std::size_t>(other)])] = 1;
          }
        }
      }
    }
    int c = 0;
    while (c < colours && used[static_cast<std::size_t>(c)] != 0) {
      ++c;
    }
    colour[l] = c;
    colours = std::max(colours, c + 1);
  }
  std::vector<std::vector<std::size_t>> by_colour(static_cast<std::size_t>(colours));
  for (std::size_t l = 0; l < lines.size(); ++l) {
    by_colour[static_cast<std::size_t>(colour[l])].push_back(l);
  }
  result.lines = static_cast<int>(lines.size());
  result.colours = colours;

  // ------------------------------------------------------------------ wall control
  struct WallNode {
    std::size_t node;
    double spacing;
    double decay;
    double relaxation;
    double max_source;
    int normal;                  // inward direction
    std::array<int, 2> tangent;  // positive directions of the two other axes
  };
  std::vector<WallNode> walls;
  // Symmetry-plane nodes (flagged face nodes sliding on a plane): interior nodes for the
  // smoothing, so wall control reaches them too.
  std::vector<std::uint8_t> symmetric(n, 0);
  for (std::size_t g = 0; g < n; ++g) {
    if (!is(static_cast<std::int64_t>(g), MultiblockNode::Sliding) ||
        topo[g].shape != Shape::Face || slide_of[g] < 0) {
      continue;
    }
    const Multiblock3DSlide& slide = input.slides[static_cast<std::size_t>(slide_of[g])];
    if (g < input.symmetric.size() && input.symmetric[g] != 0 &&
        slide.target == SlideTarget::Surface &&
        dynamic_cast<const PlaneSurface*>(input.surfaces[slide.index].get()) != nullptr) {
      symmetric[g] = 1;
    }
  }
  auto is_symmetric = [&](std::int64_t g) {
    return g >= 0 && symmetric[static_cast<std::size_t>(g)] != 0;
  };
  std::vector<std::uint8_t> on_wall(n, 0);
  for (const MultiblockWall& wall : input.walls) {
    if (wall.node < 0 || std::cmp_greater_equal(wall.node, n)) {
      throw std::invalid_argument("multiblock 3D: wall node out of range");
    }
    on_wall[static_cast<std::size_t>(wall.node)] = 1;
  }
  for (const MultiblockWall& wall : input.walls) {
    const auto g = static_cast<std::size_t>(wall.node);
    const NodeTopo& t = topo[g];
    // A wall node on the boundary controls the grid line leaving it; a fixed regular node
    // inside the grid (an internal wall, e.g. a pinned wake cut) controls the lines leaving
    // it on both sides, along the one axis whose neighbours are not wall nodes.
    std::vector<int> inward;
    if (t.consistent && t.shape == Shape::Face) {
      inward.push_back(t.inward);
    } else if (t.consistent && t.shape == Shape::Interior && input.fixed[g] != 0) {
      int normal_axis = -1;
      int count = 0;
      for (int a = 0; a < 3; ++a) {
        const std::int64_t plus = nb_dir(g, dir_of(a, 1));
        const std::int64_t minus = nb_dir(g, dir_of(a, -1));
        const bool wall_plus = plus >= 0 && on_wall[static_cast<std::size_t>(plus)] != 0;
        const bool wall_minus = minus >= 0 && on_wall[static_cast<std::size_t>(minus)] != 0;
        if (!wall_plus && !wall_minus) {
          normal_axis = a;
          count += 1;
        }
      }
      if (count == 1) {
        inward.push_back(dir_of(normal_axis, 1));
        inward.push_back(dir_of(normal_axis, -1));
      }
    } else if (t.consistent && t.shape == Shape::Edge) {
      // Where a wall meets a symmetry plane: the grid line leaving the wall in the plane.
      for (int m = 0; m < 2; ++m) {
        const int d = dir_of(t.inward_axes[static_cast<std::size_t>(m)],
                             t.inward_signs[static_cast<std::size_t>(m)]);
        if (is_symmetric(nb_dir(g, d))) {
          inward.push_back(d);
        }
      }
    }
    for (const int d : inward) {
      const int axis = dir_axis(d);
      const std::array<int, 2> tangent = {dir_of((axis + 1) % 3, 1), dir_of((axis + 2) % 3, 1)};
      double spacing = wall.spacing;
      if (!(spacing > 0.0)) {
        spacing = distance(x[static_cast<std::size_t>(nb_dir(g, d))], x[g]);
      }
      if (!(spacing > 0.0)) {
        continue;
      }
      walls.push_back({g, spacing, wall.decay, wall.relaxation, wall.max_source, d, tangent});
    }
  }
  std::vector<std::vector<WallReach>> reach(n);
  for (std::size_t w = 0; w < walls.size(); ++w) {
    auto cur = static_cast<std::int64_t>(walls[w].node);
    int normal = walls[w].normal;
    std::array<int, 2> tangent = walls[w].tangent;
    for (std::size_t steps = 1; steps <= n; ++steps) {
      const auto& nb = topo[static_cast<std::size_t>(cur)].nb;
      const std::int64_t next = nb[static_cast<std::size_t>(slot_of_dir(normal))];
      if (!is(next, MultiblockNode::Regular) && !is_symmetric(next)) {
        break;
      }
      const auto& next_nb = topo[static_cast<std::size_t>(next)].nb;
      const int back = dir_to(next_nb, cur);
      if (back < 0) {
        break;
      }
      std::array<int, 2> next_tangent{};
      bool ok = true;
      for (int m = 0; m < 2; ++m) {
        std::array<int, 3> offset{};
        offset[static_cast<std::size_t>(dir_axis(normal))] = dir_sign(normal);
        offset[static_cast<std::size_t>(dir_axis(tangent[static_cast<std::size_t>(m)]))] =
            dir_sign(tangent[static_cast<std::size_t>(m)]);
        std::int64_t side = nb[static_cast<std::size_t>(slot_of(offset))];
        bool flipped = false;
        if (side < 0) {
          // Across a symmetry plane: take the mirror side and turn the direction round.
          offset[static_cast<std::size_t>(dir_axis(tangent[static_cast<std::size_t>(m)]))] =
              -dir_sign(tangent[static_cast<std::size_t>(m)]);
          side = nb[static_cast<std::size_t>(slot_of(offset))];
          flipped = true;
        }
        int d = dir_to(next_nb, side);
        if (side < 0 || d < 0) {
          ok = false;
          break;
        }
        if (flipped) {
          d = dir_opposite(d);
        }
        next_tangent[static_cast<std::size_t>(m)] = d;
      }
      if (!ok) {
        break;
      }
      const double weight = std::exp(-walls[w].decay * static_cast<double>(steps));
      if (weight < 1e-14) {
        break;
      }
      normal = dir_opposite(back);
      tangent = next_tangent;
      reach[static_cast<std::size_t>(next)].push_back({w, weight, normal, tangent});
      cur = next;
    }
  }
  // Sources per wall node: normal, first and second tangential.
  std::vector<std::array<double, 3>> wall_s(walls.size(), {0.0, 0.0, 0.0});
  auto update_walls = [&] {
    for (std::size_t w = 0; w < walls.size(); ++w) {
      const WallNode& wall = walls[w];
      const Vec3 line = x[static_cast<std::size_t>(nb_dir(wall.node, wall.normal))] - x[wall.node];
      const double nl = norm(line);
      if (nl == 0.0) {
        continue;
      }
      const double error = (nl - wall.spacing) / wall.spacing;
      auto& s = wall_s[w];
      s[0] -= wall.relaxation * std::tanh(error);
      for (int m = 0; m < 2; ++m) {
        const int d = wall.tangent[static_cast<std::size_t>(m)];
        const std::int64_t plus = nb_dir(wall.node, d);
        const std::int64_t minus = nb_dir(wall.node, dir_opposite(d));
        if (plus < 0 || minus < 0) {
          continue;
        }
        const Vec3 tangent = x[static_cast<std::size_t>(plus)] - x[static_cast<std::size_t>(minus)];
        const double tl = norm(tangent);
        if (tl == 0.0) {
          continue;
        }
        const double lean = std::acos(std::clamp(dot(tangent, line) / (tl * nl), -1.0, 1.0)) -
                            0.5 * std::numbers::pi;
        s[static_cast<std::size_t>(m) + 1] += wall.relaxation * std::tanh(lean);
      }
      for (double& value : s) {
        value = std::clamp(value, -wall.max_source, wall.max_source);
      }
    }
  };

  // Coefficients of a regular node (frozen for one iteration: computed for all regular
  // nodes before the line sweeps; the cross-derivative term is evaluated fresh per line).
  struct Coefficients {
    std::array<double, 6> a;  // per direction
    double diag;
    double a12;
    double a13;
    double a23;
    Vec3 cross;
  };
  auto at_offset = [&](std::size_t g, const std::array<int, 3>& offset) -> const Vec3& {
    return x[static_cast<std::size_t>(topo[g].nb[static_cast<std::size_t>(slot_of(offset))])];
  };
  auto unit = [](int axis, int sign) {
    std::array<int, 3> o{};
    o[static_cast<std::size_t>(axis)] = sign;
    return o;
  };
  // Coefficients from the neighbours' positions given by `pos` (an offset -> point accessor:
  // the stored neighbours, or mirror images across a symmetry plane).
  auto coefficients_with = [&](std::size_t g, const auto& pos) {
    std::array<Vec3, 3> d;
    for (int m = 0; m < 3; ++m) {
      d[static_cast<std::size_t>(m)] = 0.5 * (pos(unit(m, 1)) - pos(unit(m, -1)));
    }
    const double g11 = dot(d[0], d[0]);
    const double g22 = dot(d[1], d[1]);
    const double g33 = dot(d[2], d[2]);
    const double g12 = dot(d[0], d[1]);
    const double g13 = dot(d[0], d[2]);
    const double g23 = dot(d[1], d[2]);
    const std::array<double, 3> a = {g22 * g33 - g23 * g23, g11 * g33 - g13 * g13,
                                     g11 * g22 - g12 * g12};
    std::array<double, 3> src{};
    for (const WallReach& r : reach[g]) {
      const auto& s = wall_s[r.wall];
      src[static_cast<std::size_t>(dir_axis(r.normal))] += dir_sign(r.normal) * s[0] * r.weight;
      for (int m = 0; m < 2; ++m) {
        const int t = r.tangent[static_cast<std::size_t>(m)];
        src[static_cast<std::size_t>(dir_axis(t))] +=
            dir_sign(t) * s[static_cast<std::size_t>(m) + 1] * r.weight;
      }
    }
    Coefficients c{};
    c.diag = 0.0;
    for (int m = 0; m < 3; ++m) {
      const auto um = static_cast<std::size_t>(m);
      c.a[static_cast<std::size_t>(dir_of(m, 1))] = a[um] * (1.0 + 0.5 * src[um]);
      c.a[static_cast<std::size_t>(dir_of(m, -1))] = a[um] * (1.0 - 0.5 * src[um]);
      c.diag += 2.0 * a[um];
    }
    c.a12 = g13 * g23 - g12 * g33;
    c.a13 = g12 * g23 - g13 * g22;
    c.a23 = g12 * g13 - g11 * g23;
    return c;
  };
  auto cross_with = [&](const Coefficients& c, const auto& pos) {
    auto cross_derivative = [&](int m1, int m2) {
      auto o = [&](int s1, int s2) {
        std::array<int, 3> v{};
        v[static_cast<std::size_t>(m1)] = s1;
        v[static_cast<std::size_t>(m2)] = s2;
        return v;
      };
      return 0.25 * (pos(o(1, 1)) - pos(o(1, -1)) - pos(o(-1, 1)) + pos(o(-1, -1)));
    };
    return 2.0 * (c.a12 * cross_derivative(0, 1) + c.a13 * cross_derivative(0, 2) +
                  c.a23 * cross_derivative(1, 2));
  };
  auto coefficients = [&](std::size_t g) {
    return coefficients_with(g,
                             [&](const std::array<int, 3>& o) -> Vec3 { return at_offset(g, o); });
  };
  auto cross_term = [&](std::size_t g, const Coefficients& c) {
    return cross_with(c, [&](const std::array<int, 3>& o) -> Vec3 { return at_offset(g, o); });
  };
  std::vector<std::size_t> regular;
  for (std::size_t g = 0; g < n; ++g) {
    if (is(static_cast<std::int64_t>(g), MultiblockNode::Regular)) {
      regular.push_back(g);
    }
  }
  std::vector<Coefficients> frozen(n);

  // Point updates.
  std::vector<std::size_t> irregular;
  std::vector<std::size_t> sliding;
  for (std::size_t g = 0; g < n; ++g) {
    if (is(static_cast<std::int64_t>(g), MultiblockNode::Irregular)) {
      irregular.push_back(g);
    } else if (is(static_cast<std::int64_t>(g), MultiblockNode::Sliding)) {
      sliding.push_back(g);
    }
  }
  // Sliding state: surface (u, v) or chain parameter; initial position and allowed radius
  // (half the smallest initial distance to the boundary neighbours).
  struct SlideState {
    const Multiblock3DSlide* slide;
    double u = 0.0;
    double v = 0.0;
    Vec3 origin;
    double radius = 0.0;
    // A face node on a plane is a symmetry-plane node: updated like an interior node with the
    // neighbours across the plane mirrored, and not limited (it is the interior of the grid's
    // sections, e.g. the end plane of an extrusion).
    bool planar = false;
    Vec3 normal;
  };
  std::vector<SlideState> slides(sliding.size());
  for (std::size_t s = 0; s < sliding.size(); ++s) {
    const std::size_t g = sliding[s];
    SlideState& st = slides[s];
    st.slide = &input.slides[static_cast<std::size_t>(slide_of[g])];
    st.origin = x[g];
    const NodeTopo& t = topo[g];
    double nearest = std::numeric_limits<double>::infinity();
    if (st.slide->target == SlideTarget::Surface) {
      const SurfaceProjection p = input.surfaces[st.slide->index]->project(x[g]);
      st.u = p.u;
      st.v = p.v;
    } else {
      st.u = input.chains[st.slide->index]->project_param(x[g]);
    }
    // Boundary neighbours: the in-face directions (face nodes) or along the edge (edge nodes).
    std::vector<int> along;
    if (t.shape == Shape::Face) {
      const int axis = dir_axis(t.inward);
      along = {dir_of((axis + 1) % 3, 1), dir_of((axis + 1) % 3, -1), dir_of((axis + 2) % 3, 1),
               dir_of((axis + 2) % 3, -1)};
    } else {
      along = {dir_of(t.edge_axis, 1), dir_of(t.edge_axis, -1)};
    }
    for (const int d : along) {
      const std::int64_t other = nb_dir(g, d);
      if (other >= 0) {
        nearest = std::min(nearest, distance(x[static_cast<std::size_t>(other)], x[g]));
      }
    }
    if (!std::isfinite(nearest)) {
      nearest = 0.0;
    }
    st.radius = options.slide_limit * nearest;
    if (st.slide->target == SlideTarget::Chain && t.shape == Shape::Edge) {
      // Between two symmetry planes: not limited either.
      const std::int64_t on1 = nb_dir(g, dir_of(t.inward_axes[1], t.inward_signs[1]));
      const std::int64_t on2 = nb_dir(g, dir_of(t.inward_axes[0], t.inward_signs[0]));
      if (is_symmetric(on1) && is_symmetric(on2)) {
        st.radius = std::numeric_limits<double>::infinity();
      }
    }
    if (symmetric[g] != 0) {
      const Surface& surface = *input.surfaces[st.slide->index];
      {
        const Vec3 nrm = cross(surface.du(st.u, st.v), surface.dv(st.u, st.v));
        const double length = norm(nrm);
        if (length > 0.0) {
          st.planar = true;
          st.normal = (1.0 / length) * nrm;
          st.radius = std::numeric_limits<double>::infinity();
        }
      }
    }
  }
  // Initial parameters and the metric there, so the limit acts in parameter space: the
  // limited position lies exactly on the geometry and limiting it again changes nothing (an
  // idempotent clamp; a clamp in space followed by a projection is not, and a smooth limiter
  // pulls every node back - both stall the accelerated iteration).
  std::vector<std::array<double, 2>> origin_param(sliding.size());
  std::vector<std::array<Vec3, 2>> origin_metric(sliding.size());
  for (std::size_t s = 0; s < sliding.size(); ++s) {
    const SlideState& st = slides[s];
    origin_param[s] = {st.u, st.v};
    if (st.slide->target == SlideTarget::Surface) {
      const Surface& surface = *input.surfaces[st.slide->index];
      origin_metric[s] = {surface.du(st.u, st.v), surface.dv(st.u, st.v)};
    } else {
      const CurveChain& chain = *input.chains[st.slide->index];
      const double h = 1e-7;
      const double a = std::max(0.0, st.u - h);
      const double b = std::min(1.0, st.u + h);
      origin_metric[s] = {(1.0 / (b - a)) * (chain.point_at_param(b) - chain.point_at_param(a)),
                          Vec3{0.0, 0.0, 0.0}};
    }
  }
  // New position of sliding node s heading for `wanted`: projected onto its geometry, the
  // parameter step from the initial parameters scaled down to the allowed radius.
  // Point of sliding node s at parameters (ut, vt), the step from the initial parameters
  // scaled down to the allowed radius; (u, v) receive the parameters used.
  auto limit = [&](std::size_t s, double ut, double vt, double& u, double& v) {
    const SlideState& st = slides[s];
    const auto& [u0, v0] = origin_param[s];
    const auto& [mu, mv] = origin_metric[s];
    double du = ut - u0;
    double dv = vt - v0;
    if (st.slide->target == SlideTarget::Chain) {
      dv = 0.0;
    }
    const double length = norm(du * mu + dv * mv);
    if (length > st.radius) {
      double scale = 0.0;
      if (length > 0.0) {
        scale = st.radius / length;
      }
      du *= scale;
      dv *= scale;
    }
    if (st.slide->target == SlideTarget::Surface) {
      u = u0 + du;
      v = v0 + dv;
      return input.surfaces[st.slide->index]->eval(u, v);
    }
    u = std::clamp(u0 + du, 0.0, 1.0);
    v = 0.0;
    return input.chains[st.slide->index]->point_at_param(u);
  };
  // Parameters of the point of node s's geometry closest to p.
  auto project = [&](std::size_t s, const Vec3& p, double& u, double& v) {
    const SlideState& st = slides[s];
    if (st.slide->target == SlideTarget::Surface) {
      const SurfaceProjection q = input.surfaces[st.slide->index]->project_from(p, st.u, st.v);
      u = q.u;
      v = q.v;
      return;
    }
    u = input.chains[st.slide->index]->project_param(p);
    v = 0.0;
  };
  auto place = [&](std::size_t s, const Vec3& wanted, double& u, double& v) {
    double ut = 0.0;
    double vt = 0.0;
    project(s, wanted, ut, vt);
    return limit(s, ut, vt, u, v);
  };
  // Unit normal of the plane a symmetry-plane node slides on.
  auto plane_normal = [&](std::size_t g) {
    const Multiblock3DSlide& slide = input.slides[static_cast<std::size_t>(slide_of[g])];
    const Surface& surface = *input.surfaces[slide.index];
    const Vec3 nrm = cross(surface.du(0.0, 0.0), surface.dv(0.0, 0.0));
    return (1.0 / norm(nrm)) * nrm;
  };
  // Point of chain node s at parameter u, not limited (edges between symmetry planes).
  auto limit_free = [&](std::size_t s, double ut, double& u, double& v) {
    u = std::clamp(ut, 0.0, 1.0);
    v = 0.0;
    return input.chains[slides[s].slide->index]->point_at_param(u);
  };
  std::vector<Vec3> slide_target(sliding.size());
  std::vector<std::array<double, 2>> slide_param(sliding.size());
  auto reproject = [&] {
    parallel_for(
        0, static_cast<std::int64_t>(sliding.size()),
        [&](std::int64_t s) {
          const auto us = static_cast<std::size_t>(s);
          double u = 0.0;
          double v = 0.0;
          slide_target[us] = place(us, x[sliding[us]], u, v);
          slide_param[us] = {u, v};
        },
        16);
    for (std::size_t s = 0; s < sliding.size(); ++s) {
      x[sliding[s]] = slide_target[s];
      slides[s].u = slide_param[s][0];
      slides[s].v = slide_param[s][1];
    }
  };

  const std::size_t scratch = line_nodes.size() + 2 * lines.size();
  std::vector<double> la(scratch);
  std::vector<double> lb(scratch);
  std::vector<double> lc(scratch);
  std::array<std::vector<double>, 3> rhs = {
      std::vector<double>(scratch), std::vector<double>(scratch), std::vector<double>(scratch)};
  std::vector<double> line_sum(lines.size(), 0.0);
  std::vector<double> point_sum(irregular.size() + sliding.size(), 0.0);
  std::vector<Vec3> irregular_target(irregular.size());

  auto relax_line = [&](std::size_t l) {
    const Line& line = lines[l];
    const std::size_t o = line.offset + 2 * l;
    const std::size_t m = line.size;
    for (std::size_t k = 0; k < m; ++k) {
      const auto g = static_cast<std::size_t>(line_nodes[line.offset + k]);
      const Coefficients& c = frozen[g];
      const int prev = line_prev[line.offset + k];
      const int next = line_next[line.offset + k];
      Vec3 b = c.cross;
      const std::size_t side = 4 * (line.offset + k);
      for (std::size_t q = side; q < side + 4; ++q) {
        const auto& [d, node] = line_side[q];
        b += c.a[static_cast<std::size_t>(d)] * x[static_cast<std::size_t>(node)];
      }
      la[o + k + 1] = -c.a[static_cast<std::size_t>(prev)];
      lb[o + k + 1] = c.diag;
      lc[o + k + 1] = -c.a[static_cast<std::size_t>(next)];
      rhs[0][o + k + 1] = b.x;
      rhs[1][o + k + 1] = b.y;
      rhs[2][o + k + 1] = b.z;
    }
    const Vec3& e0 = x[static_cast<std::size_t>(line.end0)];
    const Vec3& e1 = x[static_cast<std::size_t>(line.end1)];
    rhs[0][o] = e0.x;
    rhs[1][o] = e0.y;
    rhs[2][o] = e0.z;
    rhs[0][o + m + 1] = e1.x;
    rhs[1][o + m + 1] = e1.y;
    rhs[2][o + m + 1] = e1.z;
    if (!solve_line3(&la[o], &lb[o], &lc[o], {&rhs[0][o], &rhs[1][o], &rhs[2][o]}, m + 2)) {
      return 0.0;
    }
    double sum = 0.0;
    for (std::size_t k = 0; k < m; ++k) {
      const auto g = static_cast<std::size_t>(line_nodes[line.offset + k]);
      const Vec3 delta = Vec3{rhs[0][o + k + 1], rhs[1][o + k + 1], rhs[2][o + k + 1]} - x[g];
      x[g] += options.omega * delta;
      sum += dot(delta, delta);
    }
    return sum;
  };

  const std::size_t state_size = 3 * n + 3 * walls.size();
  auto pack = [&](std::vector<double>& state) {
    state.resize(state_size);
    std::size_t k = 0;
    for (const Vec3& p : x) {
      state[k++] = p.x;
      state[k++] = p.y;
      state[k++] = p.z;
    }
    for (const auto& s : wall_s) {
      for (const double value : s) {
        state[k++] = value;
      }
    }
  };
  auto unpack = [&](const std::vector<double>& state) {
    std::size_t k = 0;
    for (Vec3& p : x) {
      p.x = state[k++];
      p.y = state[k++];
      p.z = state[k++];
    }
    for (auto& s : wall_s) {
      for (double& value : s) {
        value = state[k++];
      }
    }
  };
  AndersonAccelerator anderson(state_size, options.acceleration);
  std::vector<double> before;
  std::vector<double> after;

  auto updates = static_cast<double>(line_nodes.size() + irregular.size() + sliding.size());
  updates = std::max(updates, 1.0);
  double first_residual = 0.0;
  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    if (progress != nullptr && progress->cancelled()) {
      result.cancelled = true;
      break;
    }
    if (options.acceleration > 0) {
      pack(before);
    }
    update_walls();
    parallel_for(
        0, static_cast<std::int64_t>(regular.size()),
        [&](std::int64_t k) {
          const std::size_t g = regular[static_cast<std::size_t>(k)];
          frozen[g] = coefficients(g);
          frozen[g].cross = cross_term(g, frozen[g]);
        },
        1024);
    for (const std::vector<std::size_t>& group : by_colour) {
      parallel_for(
          0, static_cast<std::int64_t>(group.size()),
          [&](std::int64_t k) {
            const std::size_t l = group[static_cast<std::size_t>(k)];
            line_sum[l] = relax_line(l);
          },
          4);
    }
    // Irregular nodes (Jacobi).
    parallel_for(
        0, static_cast<std::int64_t>(irregular.size()),
        [&](std::int64_t k) {
          const std::size_t g = irregular[static_cast<std::size_t>(k)];
          Vec3 edges{0.0, 0.0, 0.0};
          Vec3 completions{0.0, 0.0, 0.0};
          int count_edges = 0;
          for (const Incident& inc : incident[g]) {
            const auto& h = input.hexes[static_cast<std::size_t>(inc.hex)];
            auto at = [&](unsigned bits) {
              return x[static_cast<std::size_t>(h[static_cast<std::size_t>(inc.corner ^ bits)])];
            };
            const Vec3 a = at(1U);
            const Vec3 b = at(2U);
            const Vec3 c = at(4U);
            edges += a + b + c;
            count_edges += 3;
            completions += a + b + c - at(3U) - at(5U) - at(6U) + at(7U);
          }
          const auto hexes = static_cast<double>(incident[g].size());
          irregular_target[static_cast<std::size_t>(k)] =
              0.5 * ((1.0 / count_edges) * edges + (1.0 / hexes) * completions);
        },
        64);
    // Sliding nodes (Jacobi): the projection of the inward neighbour (faces) or of the inward
    // diagonal neighbour (edges), halfway.
    parallel_for(
        0, static_cast<std::int64_t>(sliding.size()),
        [&](std::int64_t s) {
          const auto us = static_cast<std::size_t>(s);
          const std::size_t g = sliding[us];
          const NodeTopo& t = topo[g];
          if (slides[us].planar) {
            // Mirror images across the plane through the node stand in for the missing
            // neighbours on the far side.
            const int axis = dir_axis(t.inward);
            const Vec3 here = x[g];
            const Vec3 nrm = slides[us].normal;
            auto pos = [&](const std::array<int, 3>& o) -> Vec3 {
              const std::int64_t id = t.nb[static_cast<std::size_t>(slot_of(o))];
              if (id >= 0) {
                return x[static_cast<std::size_t>(id)];
              }
              std::array<int, 3> inside = o;
              inside[static_cast<std::size_t>(axis)] = -inside[static_cast<std::size_t>(axis)];
              const Vec3 q =
                  x[static_cast<std::size_t>(t.nb[static_cast<std::size_t>(slot_of(inside))])];
              return q - (2.0 * dot(q - here, nrm)) * nrm;
            };
            const Coefficients c = coefficients_with(g, pos);
            Vec3 sum = cross_with(c, pos);
            for (int d = 0; d < 6; ++d) {
              sum += c.a[static_cast<std::size_t>(d)] * pos(unit(dir_axis(d), dir_sign(d)));
            }
            double u = 0.0;
            double v = 0.0;
            slide_target[us] = place(us, (1.0 / c.diag) * sum, u, v);
            slide_param[us] = {u, v};
            return;
          }
          std::int64_t inner = kNone;
          if (t.shape == Shape::Edge) {
            // An edge on a symmetry plane: between two planes both sides are mirrored (an
            // interior node of the sections' boundary); between a plane and another surface
            // it slides as a 2D boundary node, to the foot of its neighbour in the plane.
            const int a1 = t.inward_axes[0];
            const int a2 = t.inward_axes[1];
            const std::int64_t on1 = nb_dir(g, dir_of(a2, t.inward_signs[1]));  // face normal to a1
            const std::int64_t on2 = nb_dir(g, dir_of(a1, t.inward_signs[0]));  // face normal to a2
            const bool sym1 = is_symmetric(on1);
            const bool sym2 = is_symmetric(on2);
            if (sym1 && sym2) {
              const Vec3 here = x[g];
              const Vec3 n1 = plane_normal(static_cast<std::size_t>(on1));
              const Vec3 n2 = plane_normal(static_cast<std::size_t>(on2));
              auto pos = [&](const std::array<int, 3>& o) -> Vec3 {
                const std::int64_t id = t.nb[static_cast<std::size_t>(slot_of(o))];
                if (id >= 0) {
                  return x[static_cast<std::size_t>(id)];
                }
                std::array<int, 3> inside = o;
                const bool out1 = o[static_cast<std::size_t>(a1)] == -t.inward_signs[0];
                const bool out2 = o[static_cast<std::size_t>(a2)] == -t.inward_signs[1];
                if (out1) {
                  inside[static_cast<std::size_t>(a1)] = t.inward_signs[0];
                }
                if (out2) {
                  inside[static_cast<std::size_t>(a2)] = t.inward_signs[1];
                }
                const std::int64_t mirror_of = t.nb[static_cast<std::size_t>(slot_of(inside))];
                if (mirror_of < 0) {
                  return here;
                }
                Vec3 q = x[static_cast<std::size_t>(mirror_of)];
                if (out1) {
                  q = q - (2.0 * dot(q - here, n1)) * n1;
                }
                if (out2) {
                  q = q - (2.0 * dot(q - here, n2)) * n2;
                }
                return q;
              };
              const Coefficients c = coefficients_with(g, pos);
              Vec3 sum = cross_with(c, pos);
              for (int d = 0; d < 6; ++d) {
                sum += c.a[static_cast<std::size_t>(d)] * pos(unit(dir_axis(d), dir_sign(d)));
              }
              double uf = slides[us].u;
              double vf = slides[us].v;
              project(us, (1.0 / c.diag) * sum, uf, vf);
              double u = 0.0;
              double v = 0.0;
              slide_target[us] = limit_free(us, uf, u, v);
              slide_param[us] = {u, v};
              return;
            }
            if (sym1 || sym2) {
              inner = on1;
              if (sym2) {
                inner = on2;
              }
            }
          }
          if (inner >= 0) {
            // set above (an edge on one symmetry plane)
          } else if (t.shape == Shape::Face) {
            inner = nb_dir(g, t.inward);
          } else {
            std::array<int, 3> offset{};
            for (int m = 0; m < 2; ++m) {
              offset[static_cast<std::size_t>(t.inward_axes[static_cast<std::size_t>(m)])] =
                  t.inward_signs[static_cast<std::size_t>(m)];
            }
            inner = t.nb[static_cast<std::size_t>(slot_of(offset))];
          }
          // Halfway to the foot, in parameter space (one projection per node).
          double uf = slides[us].u;
          double vf = slides[us].v;
          if (inner >= 0) {
            project(us, x[static_cast<std::size_t>(inner)], uf, vf);
          }
          const double ut = slides[us].u + 0.5 * (uf - slides[us].u);
          const double vt = slides[us].v + 0.5 * (vf - slides[us].v);
          double u = 0.0;
          double v = 0.0;
          slide_target[us] = limit(us, ut, vt, u, v);
          slide_param[us] = {u, v};
        },
        16);
    for (std::size_t k = 0; k < irregular.size(); ++k) {
      const std::size_t g = irregular[k];
      const Vec3 delta = irregular_target[k] - x[g];
      x[g] += options.omega * delta;
      point_sum[k] = dot(delta, delta);
    }
    for (std::size_t s = 0; s < sliding.size(); ++s) {
      const std::size_t g = sliding[s];
      const Vec3 delta = slide_target[s] - x[g];
      x[g] = slide_target[s];
      slides[s].u = slide_param[s][0];
      slides[s].v = slide_param[s][1];
      point_sum[irregular.size() + s] = dot(delta, delta);
    }
    double sum = 0.0;
    for (const double value : line_sum) {
      sum += value;
    }
    for (const double value : point_sum) {
      sum += value;
    }
    const double residual = std::sqrt(sum / updates);
    if (!std::isfinite(residual)) {
      result.warnings.emplace_back("multiblock 3D: iteration diverged; returning the initial grid");
      result.xyz = input.xyz;
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
        reproject();
      }
    }
    result.iterations = iteration + 1;
    if (residual <= options.tolerance * first_residual || residual <= 1e-13 * size ||
        residual <= options.absolute_tolerance) {
      result.converged = true;
      break;
    }
  }
  int saturated = 0;
  for (std::size_t w = 0; w < walls.size(); ++w) {
    for (const double value : wall_s[w]) {
      if (std::fabs(value) >= walls[w].max_source) {
        saturated += 1;
        break;
      }
    }
  }
  if (saturated > 0) {
    result.warnings.push_back("multiblock 3D: wall control saturated at " +
                              std::to_string(saturated) +
                              " node(s); spacing/orthogonality may not be met (raise max_source "
                              "or lower decay)");
  }
  if (anderson.restarts() > 0) {
    result.warnings.push_back("multiblock: " + anderson.restart_note());
  }
  if (!result.converged && !result.cancelled) {
    result.warnings.push_back("multiblock 3D: not converged after " +
                              std::to_string(result.iterations) + " iterations");
  }
  result.xyz.resize(3 * n);
  for (std::size_t g = 0; g < n; ++g) {
    result.xyz[3 * g] = x[g].x;
    result.xyz[3 * g + 1] = x[g].y;
    result.xyz[3 * g + 2] = x[g].z;
  }
  return result;
}

}  // namespace og
