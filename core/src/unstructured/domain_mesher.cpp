#include "opengrid/core/unstructured/domain_mesher.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "opengrid/core/unstructured/advancing_front.hpp"

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
    for (std::size_t k = 0; k < loop.size(); k += 2) {
      local[u(loop[k])] = static_cast<Index>(coarse.global.size());
      coarse.global.push_back(loop[k]);
      coarse.input.points.push_back(input.points[u(loop[k])]);
    }
    for (std::size_t k = 0; k < loop.size(); k += 2) {
      const Index a = loop[k];
      const Index b = loop[(k + 2) % loop.size()];
      coarse.input.segments.push_back({local[u(a)], local[u(b)]});
      coarse.skipped[edge_key(a, b)] = loop[k + 1];
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

}  // namespace

MixedMesh mesh_domain(const TriangulationInput& input, const DomainMeshOptions& options) {
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
  const std::int64_t bad = untangle(out.points, out.triangles, out.quads, input.points.size());
  if (bad > 0) {
    out.warnings.push_back(std::to_string(bad) +
                           " cells remain inverted or non-convex after untangling");
  }
  return out;
}

}  // namespace og
