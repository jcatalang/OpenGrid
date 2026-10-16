#include "opengrid/core/structured/multiblock_elliptic.hpp"

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

// Stencil labels in a node's local frame: edge neighbours first, then diagonals.
constexpr int kE = 0;
constexpr int kN = 1;
constexpr int kW = 2;
constexpr int kS = 3;
constexpr int kNE = 4;
constexpr int kNW = 5;
constexpr int kSW = 6;
constexpr int kSE = 7;
constexpr std::int64_t kNone = -1;

int opposite(int label) {
  return (label + 2) % 4;
}

// Diagonal between two perpendicular edge labels.
int diagonal(int a, int b) {
  const int lo = std::min(a, b);
  const int hi = std::max(a, b);
  if (lo == kE && hi == kN) {
    return kNE;
  }
  if (lo == kN && hi == kW) {
    return kNW;
  }
  if (lo == kW && hi == kS) {
    return kSW;
  }
  if (lo == kE && hi == kS) {
    return kSE;
  }
  return -1;
}

// One incident quad seen from a node g: corners (g, p, d, m) counter-clockwise.
struct Fan {
  std::int64_t quad;
  std::int64_t p;
  std::int64_t d;
  std::int64_t m;
};

struct Topology {
  std::vector<std::array<std::int64_t, 8>> nb;  // stencil per node (kNone where missing)
  std::vector<std::vector<Fan>> fans;           // chained counter-clockwise
  std::vector<std::uint8_t> closed;
  std::vector<std::uint8_t> manifold;
};

Topology build_topology(const MultiblockInput& input, std::size_t n) {
  std::vector<std::vector<Fan>> incident(n);
  std::vector<std::vector<int>> corner(n);
  for (std::size_t q = 0; q < input.quads.size(); ++q) {
    const auto& quad = input.quads[q];
    for (int c = 0; c < 4; ++c) {
      const std::int64_t g = quad[static_cast<std::size_t>(c)];
      if (g < 0 || std::cmp_greater_equal(g, n)) {
        throw std::invalid_argument("multiblock: quad " + std::to_string(q) +
                                    " references a node out of range");
      }
      incident[static_cast<std::size_t>(g)].push_back(
          {static_cast<std::int64_t>(q), quad[static_cast<std::size_t>((c + 1) % 4)],
           quad[static_cast<std::size_t>((c + 2) % 4)],
           quad[static_cast<std::size_t>((c + 3) % 4)]});
      corner[static_cast<std::size_t>(g)].push_back(c);
    }
  }
  Topology topo;
  topo.nb.assign(n, {kNone, kNone, kNone, kNone, kNone, kNone, kNone, kNone});
  topo.fans.resize(n);
  topo.closed.assign(n, 0);
  topo.manifold.assign(n, 1);
  for (std::size_t g = 0; g < n; ++g) {
    const std::vector<Fan>& fans = incident[g];
    const std::size_t k = fans.size();
    if (k == 0) {
      topo.manifold[g] = 0;
      continue;
    }
    // Next fan after f: the one whose p is f's m.
    auto find_p = [&](std::int64_t node) {
      std::size_t found = k;
      for (std::size_t e = 0; e < k; ++e) {
        if (fans[e].p == node) {
          if (found != k) {
            return k + 1;  // two quads on the same side: not a fan
          }
          found = e;
        }
      }
      return found;
    };
    std::vector<std::size_t> starts;
    for (std::size_t e = 0; e < k; ++e) {
      bool preceded = false;
      for (std::size_t f = 0; f < k; ++f) {
        if (fans[f].m == fans[e].p) {
          preceded = true;
        }
      }
      if (!preceded) {
        starts.push_back(e);
      }
    }
    std::size_t start = 0;
    bool closed = false;
    if (starts.empty()) {
      closed = true;
      // Prefer the quad in which g is the first corner (so a block's interior nodes keep its
      // (i, j) frame), then the lowest quad index.
      std::size_t best = k;
      for (std::size_t e = 0; e < k; ++e) {
        if (corner[g][e] == 0 && (best == k || fans[e].quad < fans[best].quad)) {
          best = e;
        }
      }
      if (best == k) {
        best = 0;
        for (std::size_t e = 1; e < k; ++e) {
          if (fans[e].quad < fans[best].quad) {
            best = e;
          }
        }
      }
      start = best;
    } else if (starts.size() == 1) {
      start = starts.front();
    } else {
      topo.manifold[g] = 0;  // several fans meet at a point
      continue;
    }
    std::vector<Fan> chain;
    std::size_t e = start;
    while (chain.size() <= k) {
      chain.push_back(fans[e]);
      const std::size_t next = find_p(fans[e].m);
      if (next >= k || next == start) {
        break;
      }
      e = next;
    }
    if (chain.size() != k) {
      topo.manifold[g] = 0;
      continue;
    }
    topo.closed[g] = static_cast<std::uint8_t>(closed);
    auto& nb = topo.nb[g];
    if (closed && k == 4) {
      nb[kE] = chain[0].p;
      nb[kNE] = chain[0].d;
      nb[kN] = chain[1].p;
      nb[kNW] = chain[1].d;
      nb[kW] = chain[2].p;
      nb[kSW] = chain[2].d;
      nb[kS] = chain[3].p;
      nb[kSE] = chain[3].d;
    } else if (!closed && k == 2) {
      nb[kE] = chain[0].p;
      nb[kNE] = chain[0].d;
      nb[kN] = chain[0].m;
      nb[kNW] = chain[1].d;
      nb[kW] = chain[1].m;
    }
    topo.fans[g] = std::move(chain);
  }
  return topo;
}

int label_of(const std::array<std::int64_t, 8>& nb, std::int64_t node) {
  for (int l = 0; l < 8; ++l) {
    if (nb[static_cast<std::size_t>(l)] == node) {
      return l;
    }
  }
  return -1;
}

// Thomas algorithm on k = 1..n-2 with Dirichlet values at k = 0 and k = n-1 (as in the
// single-block smoother). Returns false on a zero pivot.
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

struct Line {
  std::size_t offset;  // into line_nodes / prev / next
  std::size_t size;    // interior nodes
  std::int64_t end0;
  std::int64_t end1;
};

// Wall source carried to a node: weight and the node's labels of the wall-normal (away from
// the wall) and wall-tangent directions.
struct WallReach {
  std::size_t wall;
  double weight;
  int normal;
  int tangent;
};

void add_source(int label, double value, double& p, double& q) {
  if (label == kE) {
    p += value;
  } else if (label == kW) {
    p -= value;
  } else if (label == kN) {
    q += value;
  } else if (label == kS) {
    q -= value;
  }
}

}  // namespace

MultiblockResult smooth_multiblock(const MultiblockInput& input, const MultiblockOptions& options,
                                   SmoothingProgress* progress) {
  const std::size_t n = input.xy.size() / 2;
  if (input.xy.size() != 2 * n || input.fixed.size() != n) {
    throw std::invalid_argument("multiblock: xy and fixed sizes do not match");
  }
  const Topology topo = build_topology(input, n);
  MultiblockResult result;
  auto& cls = result.node_class;
  cls.assign(n, static_cast<std::uint8_t>(MultiblockNode::Fixed));

  // Sliding nodes: chain and the boundary neighbours' parameters.
  std::vector<std::int64_t> slide_of(n, -1);
  for (std::size_t s = 0; s < input.slides.size(); ++s) {
    const MultiblockSlide& slide = input.slides[s];
    if (slide.node < 0 || std::cmp_greater_equal(slide.node, n) ||
        slide.chain >= input.chains.size() || !input.chains[slide.chain]) {
      throw std::invalid_argument("multiblock: slide " + std::to_string(s) + " is out of range");
    }
    slide_of[static_cast<std::size_t>(slide.node)] = static_cast<std::int64_t>(s);
  }
  int unchained = 0;
  int open_free = 0;
  int bad_slides = 0;
  for (std::size_t g = 0; g < n; ++g) {
    const std::size_t k = topo.fans[g].size();
    if (input.fixed[g] != 0) {
      continue;
    }
    if (topo.manifold[g] == 0) {
      unchained += 1;
      continue;
    }
    if (slide_of[g] >= 0) {
      if (topo.closed[g] == 0 && k == 2) {
        cls[g] = static_cast<std::uint8_t>(MultiblockNode::Sliding);
      } else {
        bad_slides += 1;
      }
      continue;
    }
    if (topo.closed[g] == 0) {
      open_free += 1;
      continue;
    }
    if (k == 4) {
      cls[g] = static_cast<std::uint8_t>(MultiblockNode::Regular);
    } else if (k >= 3) {
      cls[g] = static_cast<std::uint8_t>(MultiblockNode::Irregular);
    } else {
      unchained += 1;
    }
  }
  if (unchained > 0) {
    result.warnings.push_back("multiblock: " + std::to_string(unchained) +
                              " node(s) without a regular fan of quads kept fixed");
  }
  if (open_free > 0) {
    result.warnings.push_back("multiblock: " + std::to_string(open_free) +
                              " boundary node(s) neither fixed nor sliding kept fixed");
  }
  if (bad_slides > 0) {
    result.warnings.push_back("multiblock: " + std::to_string(bad_slides) +
                              " sliding node(s) not inside a boundary edge kept fixed");
  }
  auto is = [&](std::int64_t g, MultiblockNode c) {
    return g >= 0 && cls[static_cast<std::size_t>(g)] == static_cast<std::uint8_t>(c);
  };

  std::vector<double> x(n);
  std::vector<double> y(n);
  for (std::size_t g = 0; g < n; ++g) {
    x[g] = input.xy[2 * g];
    y[g] = input.xy[2 * g + 1];
  }
  double size = 0.0;
  for (std::size_t g = 0; g < n; ++g) {
    size = std::max(size, std::hypot(x[g] - x[0], y[g] - y[0]));
  }

  // ---------------------------------------------------------------- grid lines
  std::vector<std::int64_t> slot(2 * n, -1);  // line through (node, axis); -2: opening node
  std::vector<std::int64_t> line_nodes;
  std::vector<int> line_prev;
  std::vector<int> line_next;
  std::vector<Line> lines;
  for (std::size_t s = 0; s < n; ++s) {
    if (cls[s] != static_cast<std::uint8_t>(MultiblockNode::Regular)) {
      continue;
    }
    for (int axis = 0; axis < 2; ++axis) {
      if (slot[2 * s + static_cast<std::size_t>(axis)] != -1) {
        continue;
      }
      const auto line_id = static_cast<std::int64_t>(lines.size());
      // Walk from s in direction `label`; returns the nodes passed and the end node.
      auto walk = [&](int label, std::vector<std::int64_t>& nodes, std::vector<int>& axes,
                      bool& closed) {
        auto cur = static_cast<std::int64_t>(s);
        int dir = label;
        closed = false;
        while (true) {
          const std::int64_t next =
              topo.nb[static_cast<std::size_t>(cur)][static_cast<std::size_t>(dir)];
          if (std::cmp_equal(next, s)) {
            closed = true;
            return next;
          }
          if (!is(next, MultiblockNode::Regular)) {
            return next;
          }
          const int back = label_of(topo.nb[static_cast<std::size_t>(next)], cur);
          if (back < 0 || back > kS) {
            return next;
          }
          const int next_axis = back % 2;
          if (slot[2 * static_cast<std::size_t>(next) + static_cast<std::size_t>(next_axis)] !=
              -1) {
            return next;  // already on a line (only through inconsistent frames): an end
          }
          nodes.push_back(next);
          axes.push_back(next_axis);
          slot[2 * static_cast<std::size_t>(next) + static_cast<std::size_t>(next_axis)] = line_id;
          cur = next;
          dir = opposite(back);
        }
      };
      slot[2 * s + static_cast<std::size_t>(axis)] = line_id;
      std::vector<std::int64_t> forward;
      std::vector<int> forward_axes;
      bool closed = false;
      const std::int64_t end1 = walk(axis, forward, forward_axes, closed);
      std::vector<std::int64_t> ordered;
      std::int64_t end0 = kNone;
      if (closed) {
        // Open the loop at s: s is a Dirichlet end on both sides (relaxed by its other line).
        slot[2 * s + static_cast<std::size_t>(axis)] = -2;
        ordered = forward;
        end0 = static_cast<std::int64_t>(s);
      } else {
        std::vector<std::int64_t> backward;
        std::vector<int> backward_axes;
        bool unused = false;
        end0 = walk(opposite(axis), backward, backward_axes, unused);
        ordered.assign(backward.rbegin(), backward.rend());
        ordered.push_back(static_cast<std::int64_t>(s));
        ordered.insert(ordered.end(), forward.begin(), forward.end());
      }
      if (ordered.empty()) {
        continue;
      }
      Line line{line_nodes.size(), ordered.size(), end0, end1};
      for (std::size_t m = 0; m < ordered.size(); ++m) {
        std::int64_t before = end0;
        std::int64_t after = end1;
        if (m > 0) {
          before = ordered[m - 1];
        }
        if (m + 1 < ordered.size()) {
          after = ordered[m + 1];
        }
        const auto& nb = topo.nb[static_cast<std::size_t>(ordered[m])];
        line_nodes.push_back(ordered[m]);
        line_prev.push_back(label_of(nb, before));
        line_next.push_back(label_of(nb, after));
      }
      lines.push_back(line);
    }
  }
  // Greedy colouring: a line conflicts with every line through its nodes or their stencils.
  std::vector<int> colour(lines.size(), -1);
  int colours = 0;
  for (std::size_t l = 0; l < lines.size(); ++l) {
    std::vector<std::uint8_t> used(static_cast<std::size_t>(colours) + 1, 0);
    const Line& line = lines[l];
    for (std::size_t m = 0; m < line.size; ++m) {
      const std::int64_t node = line_nodes[line.offset + m];
      const auto& nb = topo.nb[static_cast<std::size_t>(node)];
      for (int k = -1; k < 8; ++k) {
        std::int64_t s = node;
        if (k >= 0) {
          s = nb[static_cast<std::size_t>(k)];
        }
        if (s < 0) {
          continue;
        }
        for (int axis = 0; axis < 2; ++axis) {
          const std::int64_t other =
              slot[2 * static_cast<std::size_t>(s) + static_cast<std::size_t>(axis)];
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

  // ------------------------------------------------------------- wall control
  // A wall node on the boundary (two quads) controls the grid line leaving it; a fixed node
  // inside the grid (four quads: an internal wall such as a pinned wake cut) controls the
  // lines leaving it on both sides, as two walls.
  struct WallNode {
    std::size_t node;
    double spacing;
    double decay;
    double relaxation;
    double max_source;
    int normal;   // label of the wall-normal direction at the wall node
    int tangent;  // label of the wall-tangent direction (normal turned clockwise)
  };
  std::vector<WallNode> walls;
  int bad_walls = 0;
  for (const MultiblockWall& wall : input.walls) {
    if (wall.node < 0 || std::cmp_greater_equal(wall.node, n)) {
      throw std::invalid_argument("multiblock: wall node out of range");
    }
    const auto g = static_cast<std::size_t>(wall.node);
    std::vector<std::pair<int, int>> sides;
    if (topo.manifold[g] != 0 && topo.closed[g] == 0 && topo.fans[g].size() == 2) {
      sides.emplace_back(kN, kE);
    } else if (topo.manifold[g] != 0 && topo.closed[g] != 0 && topo.fans[g].size() == 4 &&
               input.fixed[g] != 0) {
      sides.emplace_back(kN, kE);
      sides.emplace_back(kS, kW);
    } else {
      bad_walls += 1;  // a block corner on the wall: no grid line leaves it normally
      continue;
    }
    for (const auto& [normal, tangent] : sides) {
      double spacing = wall.spacing;
      if (!(spacing > 0.0)) {
        const auto inward = static_cast<std::size_t>(topo.nb[g][static_cast<std::size_t>(normal)]);
        spacing = std::hypot(x[inward] - x[g], y[inward] - y[g]);
      }
      if (!(spacing > 0.0)) {
        bad_walls += 1;
        continue;
      }
      walls.push_back({g, spacing, wall.decay, wall.relaxation, wall.max_source, normal, tangent});
    }
  }
  (void)bad_walls;  // corners of blocks on walls are expected (as i = 0, imax in one block)
  // Grid lines leaving each wall, through regular free nodes (across interfaces).
  std::vector<std::vector<WallReach>> reach(n);
  for (std::size_t w = 0; w < walls.size(); ++w) {
    auto cur = static_cast<std::int64_t>(walls[w].node);
    int normal = walls[w].normal;
    int tangent = walls[w].tangent;
    for (std::size_t steps = 1; steps <= n; ++steps) {
      const auto& nb = topo.nb[static_cast<std::size_t>(cur)];
      const std::int64_t next = nb[static_cast<std::size_t>(normal)];
      if (!is(next, MultiblockNode::Regular)) {
        break;
      }
      const int diag = diagonal(normal, tangent);
      if (diag < 0) {
        break;
      }
      const std::int64_t side = nb[static_cast<std::size_t>(diag)];
      const auto& next_nb = topo.nb[static_cast<std::size_t>(next)];
      const int back = label_of(next_nb, cur);
      const int next_tangent = label_of(next_nb, side);
      if (back < 0 || back > kS || next_tangent < 0 || next_tangent > kS) {
        break;
      }
      const double weight = std::exp(-walls[w].decay * static_cast<double>(steps));
      if (weight < 1e-14) {
        break;
      }
      normal = opposite(back);
      tangent = next_tangent;
      reach[static_cast<std::size_t>(next)].push_back({w, weight, normal, tangent});
      cur = next;
    }
  }
  std::vector<double> wall_p(walls.size(), 0.0);
  std::vector<double> wall_q(walls.size(), 0.0);
  auto update_walls = [&] {
    for (std::size_t w = 0; w < walls.size(); ++w) {
      const WallNode& wall = walls[w];
      const auto& nb = topo.nb[wall.node];
      const auto e = static_cast<std::size_t>(nb[static_cast<std::size_t>(wall.tangent)]);
      const auto west =
          static_cast<std::size_t>(nb[static_cast<std::size_t>(opposite(wall.tangent))]);
      const auto in = static_cast<std::size_t>(nb[static_cast<std::size_t>(wall.normal)]);
      const double tx = x[e] - x[west];
      const double ty = y[e] - y[west];
      const double nx = x[in] - x[wall.node];
      const double ny = y[in] - y[wall.node];
      const double tl = std::hypot(tx, ty);
      const double nl = std::hypot(nx, ny);
      if (tl == 0.0 || nl == 0.0) {
        continue;
      }
      const double lean = std::acos(std::clamp((tx * nx + ty * ny) / (tl * nl), -1.0, 1.0)) -
                          0.5 * std::numbers::pi;
      const double error = (nl - wall.spacing) / wall.spacing;
      wall_p[w] = std::clamp(wall_p[w] + wall.relaxation * std::tanh(lean), -wall.max_source,
                             wall.max_source);
      wall_q[w] = std::clamp(wall_q[w] - wall.relaxation * std::tanh(error), -wall.max_source,
                             wall.max_source);
    }
  };

  // Stencil of a regular node in its frame: coefficients per edge label, diagonal, cross term.
  struct Stencil {
    std::array<double, 4> a;
    double diag;
    double cross_x;
    double cross_y;
  };
  auto stencil = [&](std::size_t g) {
    const auto& nb = topo.nb[g];
    auto at = [&](int l) { return static_cast<std::size_t>(nb[static_cast<std::size_t>(l)]); };
    const double xs = 0.5 * (x[at(kE)] - x[at(kW)]);
    const double ys = 0.5 * (y[at(kE)] - y[at(kW)]);
    const double xt = 0.5 * (x[at(kN)] - x[at(kS)]);
    const double yt = 0.5 * (y[at(kN)] - y[at(kS)]);
    const double alpha = xt * xt + yt * yt;
    const double beta = xs * xt + ys * yt;
    const double gamma = xs * xs + ys * ys;
    const double xst = 0.25 * (x[at(kNE)] - x[at(kSE)] - x[at(kNW)] + x[at(kSW)]);
    const double yst = 0.25 * (y[at(kNE)] - y[at(kSE)] - y[at(kNW)] + y[at(kSW)]);
    double p = 0.0;
    double q = 0.0;
    for (const WallReach& r : reach[g]) {
      add_source(r.normal, wall_q[r.wall] * r.weight, p, q);
      add_source(r.tangent, wall_p[r.wall] * r.weight, p, q);
    }
    Stencil s{};
    s.a[kE] = alpha * (1.0 + 0.5 * p);
    s.a[kW] = alpha * (1.0 - 0.5 * p);
    s.a[kN] = gamma * (1.0 + 0.5 * q);
    s.a[kS] = gamma * (1.0 - 0.5 * q);
    s.diag = 2.0 * (alpha + gamma);
    s.cross_x = -2.0 * beta * xst;
    s.cross_y = -2.0 * beta * yst;
    return s;
  };

  // Point updates: irregular and sliding nodes.
  std::vector<std::size_t> irregular;
  std::vector<std::size_t> sliding;
  for (std::size_t g = 0; g < n; ++g) {
    if (is(static_cast<std::int64_t>(g), MultiblockNode::Irregular)) {
      irregular.push_back(g);
    } else if (is(static_cast<std::int64_t>(g), MultiblockNode::Sliding)) {
      sliding.push_back(g);
    }
  }
  // Sliding parameters; boundary neighbours that do not slide on the same chain keep the
  // parameter of their projection onto this node's chain.
  std::vector<double> u(n, 0.0);
  struct SlideBounds {
    std::int64_t east;  // sliding neighbour on the same chain, or kNone
    std::int64_t west;
    double east_u;  // fixed parameter when the neighbour does not slide on this chain
    double west_u;
    double lo0;  // within slide_limit of the initial interval on either side
    double hi0;
  };
  std::vector<SlideBounds> bounds(sliding.size());
  auto chain_of = [&](std::size_t g) -> const CurveChain& {
    return *input.chains[input.slides[static_cast<std::size_t>(slide_of[g])].chain];
  };
  for (const std::size_t g : sliding) {
    u[g] = chain_of(g).project_param(Vec3{x[g], y[g], input.z});
  }
  for (std::size_t s = 0; s < sliding.size(); ++s) {
    const std::size_t g = sliding[s];
    const CurveChain& chain = chain_of(g);
    const std::size_t my_chain = input.slides[static_cast<std::size_t>(slide_of[g])].chain;
    SlideBounds b{kNone, kNone, 0.0, 0.0, 0.0, 1.0};
    std::array<double, 2> initial{};
    for (const int side : {kE, kW}) {
      const std::int64_t other = topo.nb[g][static_cast<std::size_t>(side)];
      const auto o = static_cast<std::size_t>(other);
      const bool same = is(other, MultiblockNode::Sliding) &&
                        input.slides[static_cast<std::size_t>(slide_of[o])].chain == my_chain;
      double fixed_u = 0.0;
      if (!same) {
        fixed_u = chain.project_param(Vec3{x[o], y[o], input.z});
        initial[static_cast<std::size_t>(side / 2)] = fixed_u;
      } else {
        initial[static_cast<std::size_t>(side / 2)] = u[o];
      }
      if (side == kE) {
        if (same) {
          b.east = other;
        }
        b.east_u = fixed_u;
      } else {
        if (same) {
          b.west = other;
        }
        b.west_u = fixed_u;
      }
    }
    // Each node keeps to a fraction (slide_limit, default 1/2) of the initial interval on
    // either side, so nodes cannot pile up one after another (e.g. towards a boundary corner
    // where orthogonality is impossible).
    const double a = std::min(initial[0], initial[1]);
    const double c = std::max(initial[0], initial[1]);
    b.lo0 = u[g] - options.slide_limit * (u[g] - a);
    b.hi0 = u[g] + options.slide_limit * (c - u[g]);
    bounds[s] = b;
  }
  auto slide_range = [&](std::size_t s, double& lo, double& hi) {
    const SlideBounds& b = bounds[s];
    double ue = b.east_u;
    double uw = b.west_u;
    if (b.east != kNone) {
      ue = u[static_cast<std::size_t>(b.east)];
    }
    if (b.west != kNone) {
      uw = u[static_cast<std::size_t>(b.west)];
    }
    lo = std::min(ue, uw);
    hi = std::max(ue, uw);
    const double margin = 0.05 * (hi - lo);
    lo = std::max(lo + margin, b.lo0);
    hi = std::min(hi - margin, b.hi0);
    if (hi < lo) {
      hi = lo;
    }
  };
  // Put sliding nodes back on their curves (after an extrapolation). The bounds read the
  // neighbours' parameters, so all new parameters are computed before any is stored.
  std::vector<double> new_u(sliding.size());
  auto reproject = [&] {
    parallel_for(
        0, static_cast<std::int64_t>(sliding.size()),
        [&](std::int64_t s) {
          const std::size_t g = sliding[static_cast<std::size_t>(s)];
          double lo = 0.0;
          double hi = 1.0;
          slide_range(static_cast<std::size_t>(s), lo, hi);
          new_u[static_cast<std::size_t>(s)] =
              chain_of(g).project_param(Vec3{x[g], y[g], input.z}, lo, hi);
        },
        16);
    for (std::size_t s = 0; s < sliding.size(); ++s) {
      const std::size_t g = sliding[s];
      u[g] = new_u[s];
      const Vec3 point = chain_of(g).point_at_param(u[g]);
      x[g] = point.x;
      y[g] = point.y;
    }
  };

  // Line scratch.
  const std::size_t scratch = line_nodes.size() + 2 * lines.size();
  std::vector<double> la(scratch);
  std::vector<double> lb(scratch);
  std::vector<double> lc(scratch);
  std::vector<double> rx(scratch);
  std::vector<double> ry(scratch);
  std::vector<double> line_sum(lines.size(), 0.0);
  std::vector<double> point_sum(irregular.size() + sliding.size(), 0.0);
  std::vector<double> target_x(irregular.size());
  std::vector<double> target_y(irregular.size());
  std::vector<double> target_u(sliding.size());

  auto relax_line = [&](std::size_t l) {
    const Line& line = lines[l];
    const std::size_t o = line.offset + 2 * l;  // scratch offset: size + 2 per line
    const std::size_t m = line.size;
    for (std::size_t k = 0; k < m; ++k) {
      const auto g = static_cast<std::size_t>(line_nodes[line.offset + k]);
      const Stencil s = stencil(g);
      const int prev = line_prev[line.offset + k];
      const int next = line_next[line.offset + k];
      double bx = s.cross_x;
      double by = s.cross_y;
      for (int label = 0; label < 4; ++label) {
        if (label == prev || label == next) {
          continue;
        }
        const auto nb = static_cast<std::size_t>(topo.nb[g][static_cast<std::size_t>(label)]);
        bx += s.a[static_cast<std::size_t>(label)] * x[nb];
        by += s.a[static_cast<std::size_t>(label)] * y[nb];
      }
      la[o + k + 1] = -s.a[static_cast<std::size_t>(prev)];
      lb[o + k + 1] = s.diag;
      lc[o + k + 1] = -s.a[static_cast<std::size_t>(next)];
      rx[o + k + 1] = bx;
      ry[o + k + 1] = by;
    }
    const auto e0 = static_cast<std::size_t>(line.end0);
    const auto e1 = static_cast<std::size_t>(line.end1);
    rx[o] = x[e0];
    ry[o] = y[e0];
    rx[o + m + 1] = x[e1];
    ry[o + m + 1] = y[e1];
    if (!solve_line(&la[o], &lb[o], &lc[o], &rx[o], &ry[o], m + 2)) {
      return 0.0;
    }
    double sum = 0.0;
    for (std::size_t k = 0; k < m; ++k) {
      const auto g = static_cast<std::size_t>(line_nodes[line.offset + k]);
      const double cx = rx[o + k + 1] - x[g];
      const double cy = ry[o + k + 1] - y[g];
      x[g] += options.omega * cx;
      y[g] += options.omega * cy;
      sum += cx * cx + cy * cy;
    }
    return sum;
  };

  // Anderson state: positions and wall sources.
  const std::size_t state_size = 2 * n + 2 * walls.size();
  auto pack = [&](std::vector<double>& state) {
    state.resize(state_size);
    std::copy(x.begin(), x.end(), state.begin());
    std::copy(y.begin(), y.end(), state.begin() + static_cast<std::ptrdiff_t>(n));
    std::copy(wall_p.begin(), wall_p.end(), state.begin() + static_cast<std::ptrdiff_t>(2 * n));
    std::copy(wall_q.begin(), wall_q.end(),
              state.begin() + static_cast<std::ptrdiff_t>(2 * n + walls.size()));
  };
  auto unpack = [&](const std::vector<double>& state) {
    const auto at = [&](std::size_t k) { return state.begin() + static_cast<std::ptrdiff_t>(k); };
    std::copy(at(0), at(n), x.begin());
    std::copy(at(n), at(2 * n), y.begin());
    std::copy(at(2 * n), at(2 * n + walls.size()), wall_p.begin());
    std::copy(at(2 * n + walls.size()), at(state_size), wall_q.begin());
  };
  AndersonAccelerator anderson(state_size, options.acceleration);
  std::vector<double> before;
  std::vector<double> after;

  double updates = 0.0;
  for (const Line& line : lines) {
    updates += static_cast<double>(line.size);
  }
  updates += static_cast<double>(irregular.size() + sliding.size());
  updates = std::max(updates, 1.0);

  double first_residual = 0.0;
  double best_residual = 0.0;
  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    if (progress != nullptr && progress->cancelled()) {
      result.cancelled = true;
      break;
    }
    if (options.acceleration > 0) {
      pack(before);
    }
    update_walls();
    for (const std::vector<std::size_t>& group : by_colour) {
      parallel_for(
          0, static_cast<std::int64_t>(group.size()),
          [&](std::int64_t k) {
            const std::size_t l = group[static_cast<std::size_t>(k)];
            line_sum[l] = relax_line(l);
          },
          4);
    }
    // Irregular nodes (Jacobi): mean of the edge-neighbour average and the parallelogram
    // completions of the incident quads.
    parallel_for(
        0, static_cast<std::int64_t>(irregular.size()),
        [&](std::int64_t k) {
          const std::size_t g = irregular[static_cast<std::size_t>(k)];
          double ex = 0.0;
          double ey = 0.0;
          double px = 0.0;
          double py = 0.0;
          for (const Fan& f : topo.fans[g]) {
            const auto p = static_cast<std::size_t>(f.p);
            const auto d = static_cast<std::size_t>(f.d);
            const auto mm = static_cast<std::size_t>(f.m);
            ex += x[p];
            ey += y[p];
            px += x[p] + x[mm] - x[d];
            py += y[p] + y[mm] - y[d];
          }
          const auto count = static_cast<double>(topo.fans[g].size());
          target_x[static_cast<std::size_t>(k)] = 0.5 * (ex + px) / count;
          target_y[static_cast<std::size_t>(k)] = 0.5 * (ey + py) / count;
        },
        64);
    // Sliding nodes (Jacobi): foot of the grid line leaving the boundary.
    parallel_for(
        0, static_cast<std::int64_t>(sliding.size()),
        [&](std::int64_t s) {
          const std::size_t g = sliding[static_cast<std::size_t>(s)];
          const auto in = static_cast<std::size_t>(topo.nb[g][kN]);
          double lo = 0.0;
          double hi = 1.0;
          slide_range(static_cast<std::size_t>(s), lo, hi);
          target_u[static_cast<std::size_t>(s)] =
              chain_of(g).project_param(Vec3{x[in], y[in], input.z}, lo, hi);
        },
        16);
    for (std::size_t k = 0; k < irregular.size(); ++k) {
      const std::size_t g = irregular[k];
      const double cx = target_x[k] - x[g];
      const double cy = target_y[k] - y[g];
      x[g] += options.omega * cx;
      y[g] += options.omega * cy;
      point_sum[k] = cx * cx + cy * cy;
    }
    for (std::size_t s = 0; s < sliding.size(); ++s) {
      const std::size_t g = sliding[s];
      u[g] += 0.5 * (target_u[s] - u[g]);
      const Vec3 point = chain_of(g).point_at_param(u[g]);
      const double cx = point.x - x[g];
      const double cy = point.y - y[g];
      x[g] = point.x;
      y[g] = point.y;
      point_sum[irregular.size() + s] = cx * cx + cy * cy;
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
      result.warnings.emplace_back("multiblock: iteration diverged; returning the initial grid");
      result.xy = input.xy;
      result.iterations = iteration + 1;
      return result;
    }
    result.residuals.push_back(residual);
    if (progress != nullptr) {
      progress->record(residual);
    }
    if (iteration == 0) {
      first_residual = residual;
      best_residual = residual;
    }
    if (options.acceleration > 0) {
      if (residual > 4.0 * best_residual) {
        anderson.reset();
        best_residual = residual;
      }
      best_residual = std::min(best_residual, residual);
      pack(after);
      if (anderson.step(before, after)) {
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
    if (std::fabs(wall_p[w]) >= walls[w].max_source ||
        std::fabs(wall_q[w]) >= walls[w].max_source) {
      saturated += 1;
    }
  }
  if (saturated > 0) {
    result.warnings.push_back("multiblock: wall control saturated at " + std::to_string(saturated) +
                              " node(s); spacing/orthogonality may not be met (raise max_source "
                              "or lower decay)");
  }
  if (!result.converged && !result.cancelled) {
    result.warnings.push_back("multiblock: not converged after " +
                              std::to_string(result.iterations) + " iterations");
  }
  result.xy.resize(2 * n);
  for (std::size_t g = 0; g < n; ++g) {
    result.xy[2 * g] = x[g];
    result.xy[2 * g + 1] = y[g];
  }
  return result;
}

}  // namespace og
