#include "opengrid/core/unstructured/domain_mesher.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "opengrid/core/unstructured/advancing_front.hpp"
#include "opengrid/core/unstructured/loop_corners.hpp"

namespace og {
namespace {

using Index = std::int64_t;

std::size_t u(Index i) {
  return static_cast<std::size_t>(i);
}

std::uint64_t edge_key(Index a, Index b) {
  if (b < a) {
    std::swap(a, b);
  }
  return (static_cast<std::uint64_t>(a) << 32U) ^ static_cast<std::uint64_t>(b);
}

TriangleMesh triangulate(const TriangulationInput& input, const DomainMeshOptions& options) {
  if (options.algorithm == TriangleAlgorithm::advancing_front) {
    return advancing_front(input, options.refinement);
  }
  return triangulate_domain(input, options.refinement);
}

// Coarse boundary for all-quad meshing: every second node of each loop.
struct Coarse {
  TriangulationInput input;
  std::vector<Index> global;                         // coarse boundary point -> input point
  std::unordered_map<std::uint64_t, Index> skipped;  // coarse segment (global ids) -> input point
};

Coarse coarse_boundary(const TriangulationInput& input) {
  const std::size_t n = input.points.size();
  std::vector<Index> next(n, -1);
  for (const auto& s : input.segments) {
    if (next[u(s[0])] >= 0) {
      throw std::invalid_argument("all-quad meshing: boundary point " + std::to_string(s[0]) +
                                  " starts several segments (loops touching at a point)");
    }
    next[u(s[0])] = s[1];
  }
  Coarse coarse;
  std::vector<bool> visited(n, false);
  std::vector<Index> local(n, -1);
  for (std::size_t start = 0; start < n; ++start) {
    if (visited[start] || next[start] < 0) {
      continue;
    }
    std::vector<Index> loop;
    auto v = static_cast<Index>(start);
    while (!visited[u(v)]) {
      visited[u(v)] = true;
      loop.push_back(v);
      v = next[u(v)];
    }
    if (loop.size() % 2 != 0) {
      throw std::invalid_argument(
          "all-quad meshing needs an even number of segments on every boundary loop; the loop "
          "through boundary point " +
          std::to_string(start) + " has " + std::to_string(loop.size()));
    }
    // Coarse nodes: every second node, but sharp corners must be coarse nodes (a coarse chord
    // across a corner folds the subdivided cells next to it). Between consecutive corners the
    // stretch is paired up when its length is even; otherwise (or without corners) the half of
    // the nodes (even or odd positions) holding the sharper corners is kept.
    const std::size_t m = loop.size();
    std::vector<bool> keep(m, false);
    const std::vector<std::size_t> corners = loop_corners(input.points, loop);
    bool paired = !corners.empty();
    for (std::size_t c = 0; c < corners.size() && paired; ++c) {
      const std::size_t a = corners[c];
      const std::size_t b = corners[(c + 1) % corners.size()];
      std::size_t length = (b + m - a) % m;
      if (length == 0) {
        length = m;
      }
      paired = length % 2 == 0;
    }
    if (paired) {
      for (std::size_t c = 0; c < corners.size(); ++c) {
        const std::size_t a = corners[c];
        const std::size_t b = corners[(c + 1) % corners.size()];
        std::size_t length = (b + m - a) % m;
        if (length == 0) {
          length = m;
        }
        for (std::size_t k = 0; k < length; k += 2) {
          keep[(a + k) % m] = true;
        }
      }
    } else {
      const std::vector<double> turn = loop_turning(input.points, loop);
      std::array<double, 2> sharpness = {0.0, 0.0};
      for (std::size_t k = 0; k < m; ++k) {
        sharpness[k % 2] += turn[k] * turn[k];
      }
      std::size_t offset = 0;
      if (sharpness[1] > sharpness[0]) {
        offset = 1;
      }
      for (std::size_t k = offset; k < m; k += 2) {
        keep[k] = true;
      }
    }
    std::vector<std::size_t> kept;
    for (std::size_t k = 0; k < m; ++k) {
      if (keep[k]) {
        kept.push_back(k);
        local[u(loop[k])] = static_cast<Index>(coarse.global.size());
        coarse.global.push_back(loop[k]);
        coarse.input.points.push_back(input.points[u(loop[k])]);
      }
    }
    for (std::size_t c = 0; c < kept.size(); ++c) {
      const std::size_t k = kept[c];
      const Index a = loop[k];
      const Index b = loop[kept[(c + 1) % kept.size()]];
      coarse.input.segments.push_back({local[u(a)], local[u(b)]});
      coarse.skipped[edge_key(a, b)] = loop[(k + 1) % m];
    }
  }
  return coarse;
}

MixedMesh all_quads(const TriangulationInput& input, const DomainMeshOptions& options) {
  const Coarse coarse = coarse_boundary(input);
  DomainMeshOptions coarse_options = options;
  if (coarse_options.refinement.size.max_size) {
    *coarse_options.refinement.size.max_size *= 2.0;
  }
  for (SizeSource& source : coarse_options.refinement.size.sources) {
    source.size *= 2.0;
  }
  coarse_options.refinement.max_points = std::max<Index>(3, options.refinement.max_points / 4);
  const TriangleMesh tri = triangulate(coarse.input, coarse_options);
  const Recombined mixed = recombine(tri.points, tri.triangles, options.recombine);

  MixedMesh out;
  out.rejected = tri.rejected;
  out.warnings = tri.warnings;
  out.points = input.points;
  // Coarse point -> output point: boundary points keep their input index.
  std::vector<Index> map(tri.points.size(), -1);
  for (std::size_t k = 0; k < coarse.global.size(); ++k) {
    map[k] = coarse.global[k];
  }
  for (std::size_t k = coarse.global.size(); k < tri.points.size(); ++k) {
    map[k] = static_cast<Index>(out.points.size());
    out.points.push_back(tri.points[k]);
  }
  std::unordered_map<std::uint64_t, Index> midpoint;
  auto mid = [&](Index ca, Index cb) {
    const Index a = map[u(ca)];
    const Index b = map[u(cb)];
    const std::uint64_t key = edge_key(a, b);
    const auto found = midpoint.find(key);
    if (found != midpoint.end()) {
      return found->second;
    }
    Index id = -1;
    const auto boundary = coarse.skipped.find(key);
    if (boundary != coarse.skipped.end() && std::cmp_less(ca, coarse.global.size()) &&
        std::cmp_less(cb, coarse.global.size())) {
      id = boundary->second;
    } else {
      id = static_cast<Index>(out.points.size());
      out.points.push_back({0.5 * (out.points[u(a)][0] + out.points[u(b)][0]),
                            0.5 * (out.points[u(a)][1] + out.points[u(b)][1])});
    }
    midpoint.emplace(key, id);
    return id;
  };
  auto centre = [&](const Index* corners, int n) {
    double sx = 0.0;
    double sy = 0.0;
    for (int k = 0; k < n; ++k) {
      sx += out.points[u(map[u(corners[k])])][0];
      sy += out.points[u(map[u(corners[k])])][1];
    }
    out.points.push_back({sx / n, sy / n});
    return static_cast<Index>(out.points.size()) - 1;
  };
  for (const auto& t : mixed.triangles) {
    const Index mab = mid(t[0], t[1]);
    const Index mbc = mid(t[1], t[2]);
    const Index mca = mid(t[2], t[0]);
    const Index g = centre(t.data(), 3);
    out.quads.push_back({map[u(t[0])], mab, g, mca});
    out.quads.push_back({map[u(t[1])], mbc, g, mab});
    out.quads.push_back({map[u(t[2])], mca, g, mbc});
  }
  for (const auto& q : mixed.quads) {
    std::array<Index, 4> m{};
    for (std::size_t k = 0; k < 4; ++k) {
      m[k] = mid(q[k], q[(k + 1) % 4]);
    }
    const Index g = centre(q.data(), 4);
    for (std::size_t k = 0; k < 4; ++k) {
      out.quads.push_back({map[u(q[k])], m[k], g, m[(k + 3) % 4]});
    }
  }
  SmoothingOptions settle;
  settle.method = SmoothingMethod::laplacian;
  settle.iterations = 3;
  settle.edge_swap = false;
  smooth_mesh(out.points, out.triangles, out.quads, input.segments, input.points.size(), settle);
  return out;
}

// Layers from the wall segments, then the core region meshed without inflation.
MixedMesh inflated(const TriangulationInput& input, const DomainMeshOptions& options) {
  InflationResult layers = inflate(input, options.wall, options.inflation, options.refinement.size,
                                   options.element_type == ElementType::quads);
  TriangulationInput core;
  std::vector<Index> global;
  for (const auto& loop : layers.core_loops) {
    const auto first = static_cast<Index>(core.points.size());
    const auto n = static_cast<Index>(loop.size());
    for (const Index v : loop) {
      global.push_back(v);
      core.points.push_back(layers.points[u(v)]);
    }
    for (Index k = 0; k < n; ++k) {
      core.segments.push_back({first + k, first + (k + 1) % n});
    }
  }
  DomainMeshOptions core_options = options;
  core_options.wall.clear();
  MixedMesh inner;
  try {
    inner = mesh_domain(core, core_options);
  } catch (const std::invalid_argument& error) {
    throw std::invalid_argument(std::string("core region after inflation: ") + error.what());
  }
  MixedMesh out;
  out.points = std::move(layers.points);
  out.triangles = std::move(layers.triangles);
  out.quads = std::move(layers.quads);
  out.layers = layers.layers;
  out.warnings = std::move(layers.warnings);
  out.warnings.insert(out.warnings.end(), inner.warnings.begin(), inner.warnings.end());
  out.rejected = inner.rejected;
  for (std::size_t k = global.size(); k < inner.points.size(); ++k) {
    global.push_back(static_cast<Index>(out.points.size()));
    out.points.push_back(inner.points[k]);
  }
  for (const auto& t : inner.triangles) {
    out.triangles.push_back({global[u(t[0])], global[u(t[1])], global[u(t[2])]});
  }
  for (const auto& q : inner.quads) {
    out.quads.push_back({global[u(q[0])], global[u(q[1])], global[u(q[2])], global[u(q[3])]});
  }
  // Drop points no cell uses (nodes merged away to even out core loops); input points stay.
  std::vector<Index> map(out.points.size(), -1);
  for (std::size_t k = 0; k < input.points.size(); ++k) {
    map[k] = static_cast<Index>(k);
  }
  auto mark = [&](Index v) { map[u(v)] = 0; };
  for (const auto& t : out.triangles) {
    for (const Index v : t) {
      if (std::cmp_greater_equal(v, input.points.size())) {
        mark(v);
      }
    }
  }
  for (const auto& q : out.quads) {
    for (const Index v : q) {
      if (std::cmp_greater_equal(v, input.points.size())) {
        mark(v);
      }
    }
  }
  std::vector<std::array<double, 2>> points(
      out.points.begin(), out.points.begin() + static_cast<std::ptrdiff_t>(input.points.size()));
  for (std::size_t k = input.points.size(); k < out.points.size(); ++k) {
    if (map[k] == 0) {
      map[k] = static_cast<Index>(points.size());
      points.push_back(out.points[k]);
    }
  }
  out.points = std::move(points);
  for (auto& t : out.triangles) {
    for (Index& v : t) {
      v = map[u(v)];
    }
  }
  for (auto& q : out.quads) {
    for (Index& v : q) {
      v = map[u(v)];
    }
  }
  return out;
}

}  // namespace

MixedMesh mesh_domain(const TriangulationInput& input, const DomainMeshOptions& options) {
  if (std::find(options.wall.begin(), options.wall.end(), true) != options.wall.end()) {
    return inflated(input, options);
  }
  MixedMesh out;
  if (options.element_type == ElementType::quads) {
    // Validate the full boundary first (same messages as the triangle meshers).
    RefinementOptions check = options.refinement;
    check.refine = false;
    (void)triangulate_domain(input, check);
    out = all_quads(input, options);
  } else {
    TriangleMesh tri = triangulate(input, options);
    out.points = std::move(tri.points);
    out.rejected = tri.rejected;
    out.warnings = std::move(tri.warnings);
    if (options.element_type == ElementType::quad_dominant) {
      Recombined mixed = recombine(out.points, tri.triangles, options.recombine);
      out.triangles = std::move(mixed.triangles);
      out.quads = std::move(mixed.quads);
    } else {
      out.triangles = std::move(tri.triangles);
    }
  }
  smooth_mesh(out.points, out.triangles, out.quads, input.segments, input.points.size(),
              options.smoothing);
  // Untangle, and improve cells with a scaled Jacobian below 0.3 (a 25-degree triangle scores
  // 0.49, so well-shaped meshes are left alone).
  const std::int64_t bad = untangle(out.points, out.triangles, out.quads, input.points.size(), 0.3);
  if (bad > 0) {
    out.warnings.push_back(std::to_string(bad) +
                           " cells remain inverted or non-convex after untangling");
  }
  return out;
}

}  // namespace og
