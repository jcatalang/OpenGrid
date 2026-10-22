#include "opengrid/core/unstructured/inflation.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>

#include "opengrid/core/unstructured/loop_corners.hpp"
#include "opengrid/core/unstructured/size_field.hpp"
#include "opengrid/core/util/predicates.hpp"

namespace og {
namespace {

using Index = std::int64_t;
using P2 = std::array<double, 2>;

std::size_t u(Index i) {
  return static_cast<std::size_t>(i);
}

double dot(const P2& a, const P2& b) {
  return a[0] * b[0] + a[1] * b[1];
}

double distance(const P2& a, const P2& b) {
  return std::hypot(a[0] - b[0], a[1] - b[1]);
}

P2 normalize(const P2& v) {
  const double n = std::hypot(v[0], v[1]);
  if (n == 0.0) {
    return v;
  }
  return {v[0] / n, v[1] / n};
}

P2 left_normal(const P2& a, const P2& b) {
  return normalize({-(b[1] - a[1]), b[0] - a[0]});
}

P2 closest_on_segment(const P2& p, const P2& a, const P2& b) {
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];
  const double len2 = dx * dx + dy * dy;
  double t = 0.0;
  if (len2 > 0.0) {
    t = std::clamp(((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2, 0.0, 1.0);
  }
  return {a[0] + t * dx, a[1] + t * dy};
}

bool on_box(const P2& a, const P2& b, const P2& c) {
  return std::min(a[0], b[0]) <= c[0] && c[0] <= std::max(a[0], b[0]) &&
         std::min(a[1], b[1]) <= c[1] && c[1] <= std::max(a[1], b[1]);
}

// Closed segments ab and cd intersect. With a shared endpoint only a collinear overlap counts.
bool intersect(const P2& a, const P2& b, const P2& c, const P2& d, bool share) {
  if (std::max(a[0], b[0]) < std::min(c[0], d[0]) || std::max(c[0], d[0]) < std::min(a[0], b[0]) ||
      std::max(a[1], b[1]) < std::min(c[1], d[1]) || std::max(c[1], d[1]) < std::min(a[1], b[1])) {
    return false;
  }
  const int o1 = orient2d(a.data(), b.data(), c.data());
  const int o2 = orient2d(a.data(), b.data(), d.data());
  if (share) {
    if (o1 == 0 && o2 == 0) {
      // Collinear with a shared endpoint: overlap if the other ends point the same way.
      const P2* shared = nullptr;
      const P2* x = nullptr;
      const P2* y = nullptr;
      if (a == c) {
        shared = &a;
        x = &b;
        y = &d;
      } else if (a == d) {
        shared = &a;
        x = &b;
        y = &c;
      } else if (b == c) {
        shared = &b;
        x = &a;
        y = &d;
      } else {
        shared = &b;
        x = &a;
        y = &c;
      }
      return ((*x)[0] - (*shared)[0]) * ((*y)[0] - (*shared)[0]) +
                 ((*x)[1] - (*shared)[1]) * ((*y)[1] - (*shared)[1]) >
             0.0;
    }
    return false;
  }
  if (o1 * o2 > 0) {
    return false;
  }
  const int o3 = orient2d(c.data(), d.data(), a.data());
  const int o4 = orient2d(c.data(), d.data(), b.data());
  if (o1 * o2 < 0 && o3 * o4 < 0) {
    return true;
  }
  return (o1 == 0 && on_box(a, b, c)) || (o2 == 0 && on_box(a, b, d)) ||
         (o3 == 0 && on_box(c, d, a)) || (o4 == 0 && on_box(c, d, b));
}

// Uniform bucket grid of segments over a fixed box.
class SegmentGrid {
 public:
  SegmentGrid(double xmin, double ymin, double xmax, double ymax, double cell)
      : x0_(xmin), y0_(ymin), cell_(std::max(cell, 1e-300)) {
    nx_ = std::max(1L, static_cast<long>((xmax - xmin) / cell_) + 1);
    ny_ = std::max(1L, static_cast<long>((ymax - ymin) / cell_) + 1);
    while (nx_ * ny_ > 2'000'000) {
      cell_ *= 2.0;
      nx_ = std::max(1L, static_cast<long>((xmax - xmin) / cell_) + 1);
      ny_ = std::max(1L, static_cast<long>((ymax - ymin) / cell_) + 1);
    }
    cells_.assign(static_cast<std::size_t>(nx_ * ny_), {});
  }

  void clear() {
    for (auto& c : cells_) {
      c.clear();
    }
    seen_.clear();
  }

  void insert(Index id, const P2& a, const P2& b) {
    visit(std::min(a[0], b[0]), std::min(a[1], b[1]), std::max(a[0], b[0]), std::max(a[1], b[1]),
          [&](std::vector<Index>& bucket) { bucket.push_back(id); });
    if (std::cmp_greater_equal(id, seen_.size())) {
      seen_.resize(u(id) + 1, 0);
    }
  }

  template <class F>
  void query(double xa, double ya, double xb, double yb, F&& f) {
    stamp_ += 1;
    visit(xa, ya, xb, yb, [&](std::vector<Index>& bucket) {
      for (const Index id : bucket) {
        if (seen_[u(id)] != stamp_) {
          seen_[u(id)] = stamp_;
          f(id);
        }
      }
    });
  }

 private:
  [[nodiscard]] long cx(double x) const {
    return std::clamp(static_cast<long>(std::floor((x - x0_) / cell_)), 0L, nx_ - 1);
  }
  [[nodiscard]] long cy(double y) const {
    return std::clamp(static_cast<long>(std::floor((y - y0_) / cell_)), 0L, ny_ - 1);
  }
  template <class F>
  void visit(double xa, double ya, double xb, double yb, F&& f) {
    for (long j = cy(ya); j <= cy(yb); ++j) {
      for (long i = cx(xa); i <= cx(xb); ++i) {
        f(cells_[static_cast<std::size_t>(j * nx_ + i)]);
      }
    }
  }

  double x0_;
  double y0_;
  double cell_;
  long nx_ = 1;
  long ny_ = 1;
  std::vector<std::vector<Index>> cells_;
  std::vector<std::uint64_t> seen_;
  std::uint64_t stamp_ = 0;
};

enum class Stop { none, transition, slope, proximity };

// Proximity counts a front segment ahead only when it approaches within 60 degrees of head-on.
constexpr double kHeadOnCos = 0.5;

struct Node {
  Index point = 0;
  bool active = true;
  Stop stop = Stop::none;  // why the column ended before the last layer
  int layers = 0;
  bool preset = false;
  P2 preset_normal{};
  // Per layer:
  P2 normal{};
  P2 n_prev{};  // left normals of the adjacent front segments (zero if none)
  P2 n_next{};
  double scale = 1.0;
  double spacing = 0.0;
  P2 next{};
  double height = 0.0;          // total height of the completed layers
  double origin_spacing = 0.0;  // wall spacing at the wall node the column started from
  double target = 0.0;          // core size at the proposed point
  bool merged = false;          // shares its new point with the previous entry
  Index new_id = 0;
  // An open chain's end column sliding along the boundary (InflationOptions::slide): the
  // boundary vertices it can step through, in order, the next one, and the one proposed.
  bool slider = false;
  std::vector<Index> path{};
  std::size_t path_pos = 0;
  Index target_vertex = -1;
};

struct Chain {
  std::vector<Node> nodes;
  bool closed = false;
};

// A loop of the core boundary: boundary vertices and wall chains in order.
struct Piece {
  bool is_chain = false;
  Index value = 0;  // vertex id or chain index
};

struct Segment {
  P2 a;
  P2 b;
  Index ia;
  Index ib;
  int chain;                                 // -1 for boundary segments
  int pos;                                   // chain position of the first entry
  std::array<std::array<int, 2>, 2> owners;  // (chain, entry) of nodes to stop; chain -1 = none
};

bool cell_ok(const std::vector<P2>& pts, double max_angle_deg) {
  const std::size_t n = pts.size();
  for (std::size_t k = 0; k < n; ++k) {
    const P2& a = pts[(k + n - 1) % n];
    const P2& b = pts[k];
    const P2& c = pts[(k + 1) % n];
    if (orient2d(a.data(), b.data(), c.data()) <= 0) {
      return false;
    }
    const double ux = a[0] - b[0];
    const double uy = a[1] - b[1];
    const double vx = c[0] - b[0];
    const double vy = c[1] - b[1];
    const double angle =
        std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy) * 180.0 / std::numbers::pi;
    if (angle > max_angle_deg) {
      return false;
    }
  }
  return true;
}

class Inflator {
 public:
  Inflator(const TriangulationInput& input, const std::vector<bool>& wall,
           const InflationOptions& options, const SizeFieldOptions& size)
      : input_(input), wall_(wall), options_(options), points_(input.points) {
    // Core size field from the input boundary spacing (as the meshers compute it).
    std::vector<double> sum(input.points.size(), 0.0);
    std::vector<int> count(input.points.size(), 0);
    for (const auto& seg : input.segments) {
      const double length = distance(input.points[u(seg[0])], input.points[u(seg[1])]);
      for (const Index v : seg) {
        sum[u(v)] += length;
        count[u(v)] += 1;
      }
    }
    std::vector<std::array<double, 2>> pts;
    std::vector<double> sizes;
    for (std::size_t k = 0; k < input.points.size(); ++k) {
      if (count[k] > 0) {
        pts.push_back(input.points[k]);
        sizes.push_back(sum[k] / count[k]);
      }
    }
    size_ = std::make_unique<SizeField>(pts, sizes, size);
  }

  InflationResult run(bool even_loops) {
    consumed_.assign(points_.size(), false);
    build_chains();
    for (Chain& chain : chains_) {
      const std::size_t n = chain.nodes.size();
      for (std::size_t j = 0; j < n; ++j) {
        double sum = 0.0;
        int count = 0;
        for (const int direction : {-1, 1}) {
          const Index a = adjacent(chain, j, direction);
          if (a >= 0) {
            sum += distance(points_[u(chain.nodes[j].point)], points_[u(chain.nodes[u(a)].point)]);
            ++count;
          }
        }
        chain.nodes[j].origin_spacing = sum / std::max(count, 1);
      }
    }
    add_fans();
    int max_layers = options_.max_layers;
    if (options_.total_thickness) {
      // Smallest N with first_height * (g^N - 1) / (g - 1) >= total_thickness.
      double sum = 0.0;
      int n = 0;
      double h = options_.first_height;
      while (sum < *options_.total_thickness && n < 100000) {
        sum += h;
        h *= options_.growth_rate;
        ++n;
      }
      max_layers = std::min(max_layers, n);
    }
    double xmin = points_[0][0];
    double xmax = xmin;
    double ymin = points_[0][1];
    double ymax = ymin;
    double mean = 0.0;
    for (const auto& p : points_) {
      xmin = std::min(xmin, p[0]);
      xmax = std::max(xmax, p[0]);
      ymin = std::min(ymin, p[1]);
      ymax = std::max(ymax, p[1]);
    }
    for (const auto& s : input_.segments) {
      mean += distance(points_[u(s[0])], points_[u(s[1])]);
    }
    mean /= static_cast<double>(input_.segments.size());
    SegmentGrid old_grid(xmin, ymin, xmax, ymax, 2.0 * mean);
    SegmentGrid new_grid(xmin, ymin, xmax, ymax, 2.0 * mean);
    int layers_done = 0;
    double h = options_.first_height;
    for (int k = 0; k < max_layers; ++k, h *= options_.growth_rate) {
      if (!propose(k, h)) {
        break;
      }
      merge_converging();
      collect_old_segments(old_grid);
      settle(old_grid, new_grid);
      if (!any_active()) {
        break;
      }
      layer_height_ = h;
      commit();
      layers_done = k + 1;
    }
    InflationResult result;
    result.layers = layers_done;
    for (const Chain& chain : chains_) {
      for (std::size_t j = 0; j < chain.nodes.size(); ++j) {
        if (!chain.closed && (j == 0 || j + 1 == chain.nodes.size())) {
          continue;  // an open chain's ends taper (or slide)
        }
        const Node& node = chain.nodes[j];
        ++result.stops.columns;
        if (node.stop == Stop::transition) {
          ++result.stops.transition;
        } else if (node.stop == Stop::slope) {
          ++result.stops.slope;
        } else if (node.stop == Stop::proximity) {
          ++result.stops.proximity;
        }
      }
    }
    assemble_core(result);
    if (even_loops) {
      make_even(result);
    }
    result.points = std::move(points_);
    result.triangles = std::move(triangles_);
    result.quads = std::move(quads_);
    result.warnings.insert(result.warnings.end(), warnings_.begin(), warnings_.end());
    return result;
  }

 private:
  // ------------------------------------------------------------------ chains and fans
  void build_chains() {
    const std::size_t n = input_.points.size();
    std::vector<Index> next(n, -1);
    std::vector<Index> seg_of(n, -1);
    for (Index s = 0; std::cmp_less(s, input_.segments.size()); ++s) {
      const auto& seg = input_.segments[u(s)];
      if (next[u(seg[0])] >= 0) {
        throw std::invalid_argument("inflation: boundary point " + std::to_string(seg[0]) +
                                    " starts several segments (loops touching at a point)");
      }
      next[u(seg[0])] = seg[1];
      seg_of[u(seg[0])] = s;
    }
    std::vector<bool> visited(n, false);
    for (std::size_t start = 0; start < n; ++start) {
      if (visited[start] || next[start] < 0) {
        continue;
      }
      std::vector<Index> verts;
      auto v = static_cast<Index>(start);
      while (!visited[u(v)]) {
        visited[u(v)] = true;
        verts.push_back(v);
        v = next[u(v)];
      }
      const std::size_t m = verts.size();
      std::vector<bool> is_wall(m);
      bool all = true;
      bool none = true;
      for (std::size_t k = 0; k < m; ++k) {
        is_wall[k] = wall_[u(seg_of[u(verts[k])])];
        all = all && is_wall[k];
        none = none && !is_wall[k];
      }
      std::vector<Piece> pieces;
      if (all) {
        // A closed chain starts as far as possible from its corners (turns above 5 degrees):
        // columns are merged along the chain, not across its seam, so a seam at a corner (a
        // loop starting at one) left the two columns there unmerged, to stop as their cell
        // folded. Without corners (a smooth loop), at its first node.
        std::vector<bool> corner(m);
        bool any_corner = false;
        for (std::size_t k = 0; k < m; ++k) {
          const P2& a = points_[u(verts[(k + m - 1) % m])];
          const P2& b = points_[u(verts[k])];
          const P2& c = points_[u(verts[(k + 1) % m])];
          const P2 d1 = normalize({b[0] - a[0], b[1] - a[1]});
          const P2 d2 = normalize({c[0] - b[0], c[1] - b[1]});
          corner[k] = dot(d1, d2) < std::cos(5.0 * std::numbers::pi / 180.0);
          any_corner = any_corner || corner[k];
        }
        std::size_t seam = 0;
        if (any_corner) {
          // Nodes to the nearest corner, along the loop either way (two sweeps each way).
          std::vector<std::size_t> gap(m, m);
          std::size_t run = m;
          for (std::size_t pass = 0; pass < 2 * m; ++pass) {
            const std::size_t k = pass % m;
            if (corner[k]) {
              run = 0;
            } else if (run < m) {
              ++run;
            }
            gap[k] = std::min(gap[k], run);
          }
          run = m;
          for (std::size_t pass = 0; pass < 2 * m; ++pass) {
            const std::size_t k = (2 * m - 1 - pass) % m;
            if (corner[k]) {
              run = 0;
            } else if (run < m) {
              ++run;
            }
            gap[k] = std::min(gap[k], run);
          }
          for (std::size_t k = 1; k < m; ++k) {
            if (gap[k] > gap[seam]) {
              seam = k;
            }
          }
        }
        Chain chain;
        chain.closed = true;
        for (std::size_t k = 0; k < m; ++k) {
          chain.nodes.push_back(Node{verts[(seam + k) % m]});
        }
        pieces.push_back({true, static_cast<Index>(chains_.size())});
        chains_.push_back(std::move(chain));
      } else {
        std::size_t k0 = 0;
        while (is_wall[k0]) {
          ++k0;
        }
        std::size_t step = 0;
        while (step < m) {
          const std::size_t k = (k0 + step) % m;
          if (!is_wall[k]) {
            pieces.push_back({false, verts[k]});
            ++step;
            continue;
          }
          Chain chain;
          std::size_t run = 0;
          while (is_wall[(k + run) % m]) {
            chain.nodes.push_back(Node{verts[(k + run) % m]});
            ++run;
          }
          chain.nodes.push_back(Node{verts[(k + run) % m]});
          chain.nodes.front().active = false;
          chain.nodes.back().active = false;
          if (options_.slide) {
            // The wall leaves the front junction towards verts[k + 1], arrives at the back one
            // from verts[k + run - 1]; their neighbours on the boundary are verts[k - 1] and
            // verts[k + run + 1].
            if (junction_slides(verts[(k + 1) % m], verts[k], verts[(k + m - 1) % m])) {
              add_slider(chain.nodes.front(), verts, is_wall, k, -1);
            }
            const std::size_t back = (k + run) % m;
            if (junction_slides(verts[(back + 1) % m], verts[back], verts[(back + m - 1) % m])) {
              add_slider(chain.nodes.back(), verts, is_wall, back, 1);
            }
          }
          pieces.push_back({true, static_cast<Index>(chains_.size())});
          chains_.push_back(std::move(chain));
          step += run;
        }
      }
      (void)none;
      loops_.push_back(std::move(pieces));
    }
  }

  // Whether the layers may slide at a junction: the angle inside the domain at `here`, from the
  // direction to the next loop vertex counter-clockwise to the direction to the previous one
  // (the domain lies on the left of the loop), between the slide angle limits.
  [[nodiscard]] bool junction_slides(Index next, Index here, Index previous) const {
    const P2& b = points_[u(here)];
    const P2 to_next = {points_[u(next)][0] - b[0], points_[u(next)][1] - b[1]};
    const P2 to_previous = {points_[u(previous)][0] - b[0], points_[u(previous)][1] - b[1]};
    double angle = std::atan2(to_next[0] * to_previous[1] - to_next[1] * to_previous[0],
                              dot(to_next, to_previous));
    if (angle < 0.0) {
      angle += 2.0 * std::numbers::pi;
    }
    const double degrees = angle * 180.0 / std::numbers::pi;
    return degrees >= options_.slide_min_angle_deg && degrees <= options_.slide_max_angle_deg;
  }

  // The boundary vertices an open chain's end column can slide through: from its junction
  // (verts[junction]) along the non-wall boundary (direction -1: backwards in loop order), up
  // to a corner (a turn above 30 degrees) and never onto the start or end of another wall.
  void add_slider(Node& node, const std::vector<Index>& verts, const std::vector<bool>& is_wall,
                  std::size_t junction, int direction) {
    const std::size_t m = verts.size();
    const double straight = std::cos(30.0 * std::numbers::pi / 180.0);
    std::vector<Index> path;
    std::size_t at = junction;
    P2 previous{};
    for (std::size_t step = 0; step + 1 < m; ++step) {
      std::size_t next = (at + m - 1) % m;
      bool wall_segment = is_wall[next];  // segment verts[next] -> verts[at]
      bool beyond_is_wall = is_wall[(next + m - 1) % m];
      if (direction > 0) {
        next = (at + 1) % m;
        wall_segment = is_wall[at];  // segment verts[at] -> verts[next]
        beyond_is_wall = is_wall[next];
      }
      if (wall_segment) {
        break;
      }
      const P2& a = points_[u(verts[at])];
      const P2& b = points_[u(verts[next])];
      const P2 dir = normalize({b[0] - a[0], b[1] - a[1]});
      if (step > 0 && dot(dir, previous) < straight) {
        break;
      }
      if (beyond_is_wall) {
        break;  // verts[next] is where another wall starts or ends
      }
      path.push_back(verts[next]);
      previous = dir;
      at = next;
    }
    if (!path.empty()) {
      node.slider = true;
      node.active = true;
      node.path = std::move(path);
    }
  }

  // A sliding column's step: its next boundary vertex, if it stands for about the layer
  // height across the front (else, or past its path or a vertex another one took, it stops).
  void propose_slider(Node& node, double h) {
    if (node.path_pos >= node.path.size() || consumed_[u(node.path[node.path_pos])]) {
      node.active = false;
      return;
    }
    const Index vertex = node.path[node.path_pos];
    const P2& here = points_[u(node.point)];
    const P2& there = points_[u(vertex)];
    P2 front = node.n_prev;
    if (front[0] == 0.0 && front[1] == 0.0) {
      front = node.n_next;
    }
    const P2 step = {there[0] - here[0], there[1] - here[1]};
    const double across = std::fabs(dot(step, front));
    if (across < h / options_.slide_ratio || across > h * options_.slide_ratio) {
      node.active = false;
      return;
    }
    node.target =
        std::min((*size_)(there[0], there[1]), std::max(node.spacing, node.origin_spacing));
    if (options_.transition && h >= options_.transition_ratio * node.target) {
      node.active = false;
      node.stop = Stop::transition;
      return;
    }
    node.next = there;
    node.target_vertex = vertex;
  }

  void add_fans() {
    const double threshold = options_.fan_threshold_deg * std::numbers::pi / 180.0;
    const double step = options_.fan_step_deg * std::numbers::pi / 180.0;
    for (Chain& chain : chains_) {
      const std::size_t n = chain.nodes.size();
      std::vector<Node> out;
      for (std::size_t j = 0; j < n; ++j) {
        const Node& node = chain.nodes[j];
        const bool has_prev = chain.closed || j > 0;
        const bool has_next = chain.closed || j + 1 < n;
        if (!node.active || !has_prev || !has_next) {
          out.push_back(node);
          continue;
        }
        const P2& prev = points_[u(chain.nodes[(j + n - 1) % n].point)];
        const P2& here = points_[u(node.point)];
        const P2& next = points_[u(chain.nodes[(j + 1) % n].point)];
        const P2 d1 = {here[0] - prev[0], here[1] - prev[1]};
        const P2 d2 = {next[0] - here[0], next[1] - here[1]};
        const double cross = d1[0] * d2[1] - d1[1] * d2[0];
        const double phi = std::atan2(std::fabs(cross), dot(d1, d2));
        if (cross >= 0.0 || phi <= threshold) {
          out.push_back(node);
          continue;
        }
        // Convex corner: normals from the previous segment's to the next one's.
        const P2 n1 = left_normal(prev, here);
        const P2 n2 = left_normal(here, next);
        const double t1 = std::atan2(n1[1], n1[0]);
        double delta = std::atan2(n2[1], n2[0]) - t1;
        while (delta > std::numbers::pi) {
          delta -= 2.0 * std::numbers::pi;
        }
        while (delta <= -std::numbers::pi) {
          delta += 2.0 * std::numbers::pi;
        }
        const int m = std::max(1, static_cast<int>(std::ceil(phi / step)));
        for (int t = 0; t <= m; ++t) {
          Node fan = node;
          fan.preset = true;
          const double angle = t1 + delta * t / m;
          fan.preset_normal = {std::cos(angle), std::sin(angle)};
          out.push_back(fan);
        }
      }
      chain.nodes = std::move(out);
    }
  }

  // ------------------------------------------------------------------ one layer
  [[nodiscard]] Index neighbour(const Chain& chain, std::size_t j, int direction) const {
    const std::size_t n = chain.nodes.size();
    const Index here = chain.nodes[j].point;
    auto k = static_cast<Index>(j);
    for (std::size_t steps = 0; steps < n; ++steps) {
      k += direction;
      if (chain.closed) {
        k = (k + static_cast<Index>(n)) % static_cast<Index>(n);
      } else if (k < 0 || std::cmp_greater_equal(k, n)) {
        return -1;
      }
      if (chain.nodes[u(k)].point != here) {
        return k;
      }
    }
    return -1;
  }

  [[nodiscard]] Index adjacent(const Chain& chain, std::size_t j, int direction) const {
    const auto n = static_cast<Index>(chain.nodes.size());
    Index k = static_cast<Index>(j) + direction;
    if (chain.closed) {
      return (k + n) % n;
    }
    if (k < 0 || k >= n) {
      return -1;
    }
    return k;
  }

  bool propose(int k, double h) {
    bool any = false;
    for (Chain& chain : chains_) {
      const std::size_t n = chain.nodes.size();
      for (std::size_t j = 0; j < n; ++j) {
        Node& node = chain.nodes[j];
        const P2& here = points_[u(node.point)];
        const Index ip = neighbour(chain, j, -1);
        const Index in = neighbour(chain, j, 1);
        node.n_prev = {0.0, 0.0};
        node.n_next = {0.0, 0.0};
        double spacing = 0.0;
        int count = 0;
        if (ip >= 0) {
          const P2& p = points_[u(chain.nodes[u(ip)].point)];
          node.n_prev = left_normal(p, here);
          spacing += distance(p, here);
          ++count;
        }
        if (in >= 0) {
          const P2& p = points_[u(chain.nodes[u(in)].point)];
          node.n_next = left_normal(here, p);
          spacing += distance(p, here);
          ++count;
        }
        node.spacing = spacing / std::max(count, 1);
        node.normal = normalize({node.n_prev[0] + node.n_next[0], node.n_prev[1] + node.n_next[1]});
        if (std::hypot(node.normal[0], node.normal[1]) < 0.5) {
          node.normal = node.n_prev;
        }
        if (k == 0 && node.preset) {
          node.normal = node.preset_normal;
          node.n_prev = node.preset_normal;
          node.n_next = node.preset_normal;
        }
        if (node.slider && node.active && node.path_pos < node.path.size()) {
          // Along the boundary it slides on (its neighbours' normals are smoothed towards it).
          const P2& there = points_[u(node.path[node.path_pos])];
          node.normal = normalize({there[0] - here[0], there[1] - here[1]});
        }
      }
      // Normal smoothing (Jacobi), keeping every normal visible from both of its segments.
      for (int sweep = 0; sweep < options_.normal_smoothing; ++sweep) {
        std::vector<P2> smoothed(n);
        for (std::size_t j = 0; j < n; ++j) {
          const Node& node = chain.nodes[j];
          smoothed[j] = node.normal;
          if (!node.active || (k == 0 && node.preset)) {
            continue;
          }
          const Index a = adjacent(chain, j, -1);
          const Index b = adjacent(chain, j, 1);
          if (a < 0 || b < 0) {
            continue;
          }
          const P2& na = chain.nodes[u(a)].normal;
          const P2& nb = chain.nodes[u(b)].normal;
          const P2 candidate = normalize(
              {node.normal[0] + 0.25 * (na[0] + nb[0]), node.normal[1] + 0.25 * (na[1] + nb[1])});
          if (dot(candidate, node.n_prev) > 0.3 && dot(candidate, node.n_next) > 0.3) {
            smoothed[j] = candidate;
          }
        }
        for (std::size_t j = 0; j < n; ++j) {
          chain.nodes[j].normal = smoothed[j];
        }
      }
      for (std::size_t j = 0; j < n; ++j) {
        Node& node = chain.nodes[j];
        double worst = 1.0;
        if (node.n_prev[0] != 0.0 || node.n_prev[1] != 0.0) {
          worst = std::min(worst, dot(node.normal, node.n_prev));
        }
        if (node.n_next[0] != 0.0 || node.n_next[1] != 0.0) {
          worst = std::min(worst, dot(node.normal, node.n_next));
        }
        node.scale = 1.0 / std::max(worst, 0.5);
      }
      for (std::size_t j = 0; j < n; ++j) {
        Node& node = chain.nodes[j];
        node.merged = false;
        const P2& here = points_[u(node.point)];
        node.next = here;
        if (!node.active) {
          continue;
        }
        if (node.slider) {
          propose_slider(node, h);
          continue;
        }
        const P2 next = {here[0] + h * node.scale * node.normal[0],
                         here[1] + h * node.scale * node.normal[1]};
        // Isotropic size of the cells around this column: its width along the front (never less
        // than the wall spacing it started from, so fan columns, which start microscopically
        // close together, are not stopped at once), and never more than the core size.
        node.target =
            std::min((*size_)(next[0], next[1]), std::max(node.spacing, node.origin_spacing));
        // Transition: the layer would be as thick as its cells are wide.
        if (options_.transition && h >= options_.transition_ratio * node.target) {
          node.active = false;
          node.stop = Stop::transition;
          continue;
        }
        node.next = next;
      }
      // Slope towards stopped neighbours, repeated until no column stops (a stop can make the
      // next column too steep), so the result does not depend on the processing order.
      bool stopped = true;
      while (stopped) {
        stopped = false;
        for (std::size_t j = 0; j < n; ++j) {
          Node& node = chain.nodes[j];
          if (!node.active) {
            continue;
          }
          const P2& here = points_[u(node.point)];
          for (const int direction : {-1, 1}) {
            const Index a = adjacent(chain, j, direction);
            if (a < 0 || chain.nodes[u(a)].active) {
              continue;
            }
            const Node& other = chain.nodes[u(a)];
            // Gap along the front: the offset to the neighbour's top minus its component along
            // this column's normal (the plain distance would already contain the height
            // difference and never trigger the rule).
            const P2& top = points_[u(other.point)];
            const P2 offset = {top[0] - here[0], top[1] - here[1]};
            const double along = dot(offset, node.normal);
            const double gap =
                std::hypot(offset[0] - along * node.normal[0], offset[1] - along * node.normal[1]);
            // A single-layer step is always allowed (otherwise, once layers are thicker than
            // max_slope * gap, one stop would sweep across the whole front).
            // (Relative tolerance: with equal heights the difference is h up to rounding.)
            if (node.height + h - other.height >
                std::max(options_.max_slope * gap, h) * (1.0 + 1e-9)) {
              node.active = false;
              node.stop = Stop::slope;
              node.next = here;
              stopped = true;
              break;
            }
          }
        }
      }
      for (const Node& node : chain.nodes) {
        any = any || node.active;
      }
    }
    return any;
  }

  // Concave corners: adjacent columns whose new front segment shrinks below half of the old
  // one, or turns back (one would pass the other), are merged; a run of such pairs (both
  // neighbours of a corner column, whose layers have grown taller than the wall spacing) is
  // merged into one node at the run's mean, so both walls' layers meet there.
  void merge_converging() {
    for (Chain& chain : chains_) {
      const std::size_t n = chain.nodes.size();
      std::vector<bool> join(n, false);  // join[j]: node j + 1 merges with node j
      for (std::size_t j = 0; j + 1 < n; ++j) {
        const Node& a = chain.nodes[j];
        const Node& b = chain.nodes[j + 1];
        if (!a.active || !b.active || a.point == b.point || a.slider || b.slider) {
          continue;
        }
        const P2& pa = points_[u(a.point)];
        const P2& pb = points_[u(b.point)];
        const P2 old_segment = {pb[0] - pa[0], pb[1] - pa[1]};
        const P2 new_segment = {b.next[0] - a.next[0], b.next[1] - a.next[1]};
        const double old_length = std::hypot(old_segment[0], old_segment[1]);
        join[j] = std::hypot(new_segment[0], new_segment[1]) < 0.5 * old_length ||
                  dot(old_segment, new_segment) <= 0.0;
      }
      std::size_t first = 0;
      while (first < n) {
        std::size_t last = first;
        while (last + 1 < n && join[last]) {
          ++last;
        }
        if (last > first) {
          P2 mean = {0.0, 0.0};
          for (std::size_t k = first; k <= last; ++k) {
            mean[0] += chain.nodes[k].next[0];
            mean[1] += chain.nodes[k].next[1];
          }
          const auto count = static_cast<double>(last - first + 1);
          mean = {mean[0] / count, mean[1] / count};
          for (std::size_t k = first; k <= last; ++k) {
            chain.nodes[k].next = mean;
            chain.nodes[k].merged = k > first;
          }
        }
        first = last + 1;
      }
    }
  }

  void stop(int c, std::size_t j) {
    Chain& chain = chains_[u(c)];
    Node& node = chain.nodes[j];
    if (!node.active) {
      return;
    }
    node.active = false;
    node.stop = Stop::proximity;
    node.next = points_[u(node.point)];
    // A merged pair stops together.
    if (node.merged && j > 0) {
      node.merged = false;
      stop(c, j - 1);
    }
    if (j + 1 < chain.nodes.size() && chain.nodes[j + 1].merged) {
      chain.nodes[j + 1].merged = false;
      stop(c, j + 1);
    }
  }

  [[nodiscard]] bool any_active() const {
    for (const Chain& chain : chains_) {
      for (const Node& node : chain.nodes) {
        if (node.active) {
          return true;
        }
      }
    }
    return false;
  }

  // Ids of new points (negative) or existing points (>= 0).
  void assign_new_ids() {
    Index next = -1;
    for (Chain& chain : chains_) {
      for (std::size_t j = 0; j < chain.nodes.size(); ++j) {
        Node& node = chain.nodes[j];
        if (!node.active) {
          node.new_id = node.point;
        } else if (node.slider) {
          node.new_id = node.target_vertex;  // an existing boundary vertex
        } else if (node.merged && j > 0) {
          node.new_id = chain.nodes[j - 1].new_id;
        } else {
          node.new_id = next--;
        }
      }
    }
  }

  void collect_old_segments(SegmentGrid& grid) {
    old_.clear();
    grid.clear();
    for (std::size_t s = 0; s < input_.segments.size(); ++s) {
      if (wall_[s]) {
        continue;
      }
      const auto& seg = input_.segments[s];
      old_.push_back(
          {points_[u(seg[0])], points_[u(seg[1])], seg[0], seg[1], -1, 0, {{{-1, 0}, {-1, 0}}}});
    }
    for (int c = 0; std::cmp_less(c, chains_.size()); ++c) {
      const Chain& chain = chains_[u(c)];
      const std::size_t n = chain.nodes.size();
      std::size_t count = n - 1;
      if (chain.closed) {
        count = n;
      }
      for (std::size_t j = 0; j < count; ++j) {
        const Node& a = chain.nodes[j];
        const Node& b = chain.nodes[(j + 1) % n];
        if (a.point == b.point) {
          continue;
        }
        old_.push_back({points_[u(a.point)],
                        points_[u(b.point)],
                        a.point,
                        b.point,
                        c,
                        static_cast<int>(j),
                        {{{-1, 0}, {-1, 0}}}});
      }
    }
    for (Index s = 0; std::cmp_less(s, old_.size()); ++s) {
      grid.insert(s, old_[u(s)].a, old_[u(s)].b);
    }
  }

  void collect_new_segments(SegmentGrid& grid) {
    new_.clear();
    grid.clear();
    for (int c = 0; std::cmp_less(c, chains_.size()); ++c) {
      const Chain& chain = chains_[u(c)];
      const std::size_t n = chain.nodes.size();
      for (std::size_t j = 0; j < n; ++j) {
        const Node& node = chain.nodes[j];
        // (A sliding column's edge is a boundary segment already.)
        if (node.active && !(node.merged && j > 0) && !node.slider) {
          new_.push_back({points_[u(node.point)],
                          node.next,
                          node.point,
                          node.new_id,
                          c,
                          static_cast<int>(j),
                          {{{c, static_cast<int>(j)}, {-1, 0}}}});
        }
      }
      std::size_t count = n - 1;
      if (chain.closed) {
        count = n;
      }
      for (std::size_t j = 0; j < count; ++j) {
        const std::size_t j1 = (j + 1) % n;
        const Node& a = chain.nodes[j];
        const Node& b = chain.nodes[j1];
        if ((!a.active && !b.active) || a.new_id == b.new_id) {
          continue;
        }
        new_.push_back({a.next,
                        b.next,
                        a.new_id,
                        b.new_id,
                        c,
                        static_cast<int>(j),
                        {{{c, static_cast<int>(j)}, {c, static_cast<int>(j1)}}}});
      }
    }
    for (Index s = 0; std::cmp_less(s, new_.size()); ++s) {
      grid.insert(s, new_[u(s)].a, new_[u(s)].b);
    }
  }

  void stop_owners(const Segment& s) {
    for (const auto& [c, j] : s.owners) {
      if (c >= 0) {
        stop(c, u(j));
      }
    }
  }

  [[nodiscard]] bool near_on_chain(const Segment& s, int c, int j) const {
    if (s.chain != c) {
      return false;
    }
    const auto n = static_cast<int>(chains_[u(c)].nodes.size());
    int d = std::abs(s.pos - j);
    if (chains_[u(c)].closed) {
      d = std::min(d, n - d);
    }
    return d <= 2;
  }

  void settle(SegmentGrid& old_grid, SegmentGrid& new_grid) {
    bool changed = true;
    while (changed) {
      changed = false;
      assign_new_ids();
      // Cells.
      for (int c = 0; std::cmp_less(c, chains_.size()); ++c) {
        Chain& chain = chains_[u(c)];
        const std::size_t n = chain.nodes.size();
        std::size_t count = n - 1;
        if (chain.closed) {
          count = n;
        }
        for (std::size_t j = 0; j < count; ++j) {
          const std::size_t j1 = (j + 1) % n;
          const Node& a = chain.nodes[j];
          const Node& b = chain.nodes[j1];
          if (!a.active && !b.active) {
            continue;
          }
          std::vector<P2> pts;
          cell_points(a, b, pts);
          if (pts.size() >= 3 && !cell_ok(pts, options_.max_angle_deg)) {
            stop(c, j);
            stop(c, j1);
            changed = true;
          }
        }
      }
      if (changed) {
        continue;
      }
      assign_new_ids();
      collect_new_segments(new_grid);
      // Crossings with the old fronts and the boundary, and among the new segments.
      for (const Segment& s : new_) {
        bool bad = false;
        const double xa = std::min(s.a[0], s.b[0]);
        const double ya = std::min(s.a[1], s.b[1]);
        const double xb = std::max(s.a[0], s.b[0]);
        const double yb = std::max(s.a[1], s.b[1]);
        old_grid.query(xa, ya, xb, yb, [&](Index id) {
          const Segment& o = old_[u(id)];
          const bool share = s.ia == o.ia || s.ia == o.ib || s.ib == o.ia || s.ib == o.ib;
          bad = bad || intersect(s.a, s.b, o.a, o.b, share);
        });
        new_grid.query(xa, ya, xb, yb, [&](Index id) {
          const Segment& o = new_[u(id)];
          if (&o == &s) {
            return;
          }
          const bool share = s.ia == o.ia || s.ia == o.ib || s.ib == o.ia || s.ib == o.ib;
          if (intersect(s.a, s.b, o.a, o.b, share)) {
            bad = true;
            stop_owners(o);
          }
        });
        if (bad) {
          stop_owners(s);
          changed = true;
        }
      }
      if (changed) {
        continue;
      }
      // Proximity: a facing front or boundary ahead, closer than proximity * spacing.
      for (int c = 0; std::cmp_less(c, chains_.size()); ++c) {
        Chain& chain = chains_[u(c)];
        for (std::size_t j = 0; j < chain.nodes.size(); ++j) {
          Node& node = chain.nodes[j];
          if (!node.active) {
            continue;
          }
          const double r = options_.proximity * node.target;
          bool close = false;
          auto check = [&](const Segment& o) {
            if (close || near_on_chain(o, c, static_cast<int>(j)) || o.ia == node.new_id ||
                o.ib == node.new_id) {
              return;
            }
            const P2 q = closest_on_segment(node.next, o.a, o.b);
            const P2 ahead = {q[0] - node.next[0], q[1] - node.next[1]};
            if (std::hypot(ahead[0], ahead[1]) >= r || dot(ahead, node.normal) <= 0.0) {
              return;
            }
            const P2 facing = left_normal(o.a, o.b);
            // Facing it and approaching head-on (within 60 degrees): the other wall's layers
            // beside it in a concave corner run across its path and meet it there instead.
            if (dot(facing, {-ahead[0], -ahead[1]}) > 0.0 &&
                dot(facing, node.normal) < -kHeadOnCos) {
              close = true;
            }
          };
          old_grid.query(node.next[0] - r, node.next[1] - r, node.next[0] + r, node.next[1] + r,
                         [&](Index id) { check(old_[u(id)]); });
          new_grid.query(node.next[0] - r, node.next[1] - r, node.next[0] + r, node.next[1] + r,
                         [&](Index id) {
                           const Segment& o = new_[u(id)];
                           // Only the new front (not column edges) counts as an obstacle.
                           if (o.owners[1][0] >= 0) {
                             check(o);
                           }
                         });
          if (close) {
            stop(c, j);
            changed = true;
          }
        }
      }
    }
  }

  // Corners of the layer cell under the front segment (a, b): (old a, old b, new b, new a)
  // without repeated points.
  void cell_points(const Node& a, const Node& b, std::vector<P2>& pts) const {
    const std::array<Index, 4> ids = {a.point, b.point, b.new_id, a.new_id};
    const std::array<P2, 4> coords = {points_[u(a.point)], points_[u(b.point)], b.next, a.next};
    pts.clear();
    std::array<Index, 4> kept{};
    std::size_t m = 0;
    for (std::size_t k = 0; k < 4; ++k) {
      if (m > 0 && kept[m - 1] == ids[k]) {
        continue;
      }
      kept[m] = ids[k];
      pts.push_back(coords[k]);
      ++m;
    }
    if (m > 1 && kept[m - 1] == kept[0]) {
      pts.pop_back();
    }
  }

  void commit() {
    for (Chain& chain : chains_) {
      const std::size_t n = chain.nodes.size();
      // New points.
      std::vector<Index> new_index(n);
      for (std::size_t j = 0; j < n; ++j) {
        const Node& node = chain.nodes[j];
        if (!node.active) {
          new_index[j] = node.point;
        } else if (node.slider) {
          new_index[j] = node.target_vertex;
        } else if (node.merged && j > 0) {
          new_index[j] = new_index[j - 1];
        } else {
          new_index[j] = static_cast<Index>(points_.size());
          points_.push_back(node.next);
        }
      }
      // Cells.
      std::size_t count = n - 1;
      if (chain.closed) {
        count = n;
      }
      for (std::size_t j = 0; j < count; ++j) {
        const std::size_t j1 = (j + 1) % n;
        const std::array<Index, 4> ids = {chain.nodes[j].point, chain.nodes[j1].point,
                                          new_index[j1], new_index[j]};
        std::array<Index, 4> kept{};
        std::size_t m = 0;
        for (const Index id : ids) {
          if (m > 0 && kept[m - 1] == id) {
            continue;
          }
          kept[m++] = id;
        }
        if (m > 1 && kept[m - 1] == kept[0]) {
          --m;
        }
        if (m == 4) {
          quads_.push_back(kept);
        } else if (m == 3) {
          triangles_.push_back({kept[0], kept[1], kept[2]});
        }
      }
      // Advance the front.
      std::vector<Node> kept_nodes;
      kept_nodes.reserve(n);
      for (std::size_t j = 0; j < n; ++j) {
        Node node = chain.nodes[j];
        if (node.active) {
          node.layers += 1;
          node.height += layer_height_;
        }
        if (node.active && node.slider) {
          consumed_[u(node.point)] = true;  // a corner of layer cells, off the core loops
          node.path_pos += 1;
        }
        node.point = new_index[j];
        if (node.merged && j > 0) {
          continue;
        }
        node.merged = false;
        node.preset = false;
        kept_nodes.push_back(node);
      }
      chain.nodes = std::move(kept_nodes);
    }
  }

  // ------------------------------------------------------------------ core loops
  void assemble_core(InflationResult& result) {
    for (const auto& pieces : loops_) {
      std::vector<Index> loop;
      for (const Piece& piece : pieces) {
        if (!piece.is_chain) {
          if (!consumed_[u(piece.value)]) {
            loop.push_back(piece.value);
          }
          continue;
        }
        const Chain& chain = chains_[u(piece.value)];
        std::size_t count = chain.nodes.size() - 1;
        if (chain.closed) {
          count = chain.nodes.size();
        }
        for (std::size_t j = 0; j < count; ++j) {
          loop.push_back(chain.nodes[j].point);
        }
      }
      std::vector<Index> clean;
      for (const Index v : loop) {
        if (clean.empty() || clean.back() != v) {
          clean.push_back(v);
        }
      }
      while (clean.size() > 1 && clean.back() == clean.front()) {
        clean.pop_back();
      }
      result.core_loops.push_back(std::move(clean));
    }
  }

  // All-quad cores: the coarse boundary keeps the loop's sharp corners (loop_corners.hpp) and
  // pairs up the nodes between them, so every stretch between consecutive corners (or the whole
  // loop, without corners) needs an even number of segments. An odd stretch loses one node: the
  // shortest pair of adjacent layer nodes on a smooth part of it (both nodes and their neighbours
  // turn less than 20 degrees) is merged.
  void make_even(InflationResult& result) {
    for (auto& loop : result.core_loops) {
      // 1. Odd stretches between corners, while one of them can be fixed.
      for (int attempt = 0; attempt < 64; ++attempt) {
        const std::size_t m = loop.size();
        const std::vector<std::size_t> corners = loop_corners(points_, loop);
        bool merged = false;
        for (std::size_t c = 0; c < corners.size() && !merged; ++c) {
          const std::size_t a = corners[c];
          std::size_t length = (corners[(c + 1) % corners.size()] + m - a) % m;
          if (length == 0) {
            length = m;
          }
          if (length % 2 != 0) {
            merged = merge_flat_pair(loop, a, length);
          }
        }
        if (!merged) {
          break;
        }
      }
      // 2. The loop as a whole must be even (the coarse boundary then pairs by parity).
      if (loop.size() % 2 != 0 && !merge_flat_pair(loop, 0, loop.size()) &&
          !merge_flat_pair(loop, 0, loop.size(), 1e9)) {
        warnings_.emplace_back(
            "inflation: a core loop has an odd number of segments and no pair of layer nodes "
            "could be merged (the all-quad core will be rejected)");
      }
    }
  }

  // Merges the shortest pair of adjacent layer nodes inside loop[start .. start + length) on a
  // smooth part of the loop (the nodes and their neighbours turn less than 20 degrees).
  bool merge_flat_pair(std::vector<Index>& loop, std::size_t start, std::size_t length,
                       double smooth_deg = 20.0) {
    const auto n_input = static_cast<Index>(input_.points.size());
    const double smooth = smooth_deg * std::numbers::pi / 180.0;
    const std::size_t m = loop.size();
    const std::vector<double> turn = loop_turning(points_, loop);
    std::vector<std::pair<double, std::size_t>> candidates;
    for (std::size_t d = 0; d < length; ++d) {
      const std::size_t k = (start + d) % m;
      const Index a = loop[k];
      const Index b = loop[(k + 1) % m];
      bool flat = true;
      for (std::size_t e = 0; e < 4; ++e) {
        flat = flat && turn[(k + m - 1 + e) % m] < smooth;
      }
      if (a >= n_input && b >= n_input && flat) {
        candidates.emplace_back(distance(points_[u(a)], points_[u(b)]), k);
      }
    }
    std::sort(candidates.begin(), candidates.end());
    for (const auto& [length_ab, k] : candidates) {
      if (try_merge(loop[k], loop[(k + 1) % m])) {
        loop.erase(loop.begin() + static_cast<std::ptrdiff_t>((k + 1) % m));
        return true;
      }
    }
    return false;
  }

  bool try_merge(Index a, Index b) {
    const P2 old_a = points_[u(a)];
    points_[u(a)] = {0.5 * (old_a[0] + points_[u(b)][0]), 0.5 * (old_a[1] + points_[u(b)][1])};
    std::vector<std::array<Index, 3>> tris;
    std::vector<std::array<Index, 4>> quads;
    bool ok = true;
    auto rebuild = [&](const Index* v, std::size_t n) {
      std::array<Index, 4> kept{};
      std::size_t m = 0;
      bool touched = false;
      for (std::size_t k = 0; k < n; ++k) {
        Index id = v[k];
        if (id == b) {
          id = a;
        }
        touched = touched || v[k] == a || v[k] == b;
        if (m > 0 && kept[m - 1] == id) {
          continue;
        }
        kept[m++] = id;
      }
      if (m > 1 && kept[m - 1] == kept[0]) {
        --m;
      }
      if (touched && m >= 3) {
        std::vector<P2> pts;
        pts.reserve(m);
        for (std::size_t k = 0; k < m; ++k) {
          pts.push_back(points_[u(kept[k])]);
        }
        ok = ok && cell_ok(pts, 179.0);
      }
      if (m == 4) {
        quads.push_back(kept);
      } else if (m == 3) {
        tris.push_back({kept[0], kept[1], kept[2]});
      }
    };
    for (const auto& t : triangles_) {
      rebuild(t.data(), 3);
    }
    for (const auto& q : quads_) {
      rebuild(q.data(), 4);
    }
    if (!ok) {
      points_[u(a)] = old_a;
      return false;
    }
    triangles_ = std::move(tris);
    quads_ = std::move(quads);
    return true;
  }

  std::unique_ptr<SizeField> size_;
  double layer_height_ = 0.0;
  const TriangulationInput& input_;
  const std::vector<bool>& wall_;
  InflationOptions options_;
  std::vector<P2> points_;
  std::vector<Chain> chains_;
  std::vector<std::vector<Piece>> loops_;
  std::vector<bool> consumed_;  // input vertices a sliding column passed (off the core loops)
  std::vector<Segment> old_;
  std::vector<Segment> new_;
  std::vector<std::array<Index, 3>> triangles_;
  std::vector<std::array<Index, 4>> quads_;
  std::vector<std::string> warnings_;
};

}  // namespace

InflationResult inflate(const TriangulationInput& input, const std::vector<bool>& wall,
                        const InflationOptions& options, const SizeFieldOptions& size,
                        bool even_loops) {
  if (wall.size() != input.segments.size()) {
    throw std::invalid_argument("inflation: one wall flag per boundary segment is required");
  }
  if (!(options.first_height > 0.0) || !(options.growth_rate >= 1.0) || options.max_layers < 0 ||
      (options.total_thickness && !(*options.total_thickness > 0.0)) ||
      !(options.transition_ratio > 0.0) || !(options.fan_step_deg > 0.0) ||
      !(options.proximity >= 0.0) || !(options.max_slope > 0.0) || !(options.slide_ratio > 1.0)) {
    throw std::invalid_argument(
        "inflation: first height and total thickness must be positive, growth rate >= 1, "
        "layers >= 0");
  }
  RefinementOptions check;
  check.refine = false;
  (void)triangulate_domain(input, check);
  Inflator inflator(input, wall, options, size);
  return inflator.run(even_loops);
}

}  // namespace og
