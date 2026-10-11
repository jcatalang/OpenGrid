#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "opengrid/core/quality/quality_cells.hpp"
#include "opengrid/core/unstructured/domain_mesher.hpp"
#include "opengrid/core/unstructured/recombine.hpp"
#include "opengrid/core/unstructured/smoothing.hpp"

namespace {

using og::MeshIndex;
using Point = std::array<double, 2>;

std::size_t u(MeshIndex i) {
  return static_cast<std::size_t>(i);
}

void add_loop(og::TriangulationInput& in, const std::vector<Point>& pts) {
  const auto first = static_cast<MeshIndex>(in.points.size());
  const auto n = static_cast<MeshIndex>(pts.size());
  for (const auto& p : pts) {
    in.points.push_back(p);
  }
  for (MeshIndex k = 0; k < n; ++k) {
    in.segments.push_back({first + k, first + (k + 1) % n});
  }
}

std::vector<Point> circle(double cx, double cy, double r, int n, bool ccw) {
  std::vector<Point> pts;
  for (int k = 0; k < n; ++k) {
    double t = 2.0 * std::numbers::pi * k / n;
    if (!ccw) {
      t = -t;
    }
    pts.push_back({cx + r * std::cos(t), cy + r * std::sin(t)});
  }
  return pts;
}

double polygon_area(const std::vector<Point>& pts) {
  double a = 0.0;
  for (std::size_t k = 0; k < pts.size(); ++k) {
    const auto& p = pts[k];
    const auto& q = pts[(k + 1) % pts.size()];
    a += 0.5 * (p[0] * q[1] - q[0] * p[1]);
  }
  return a;
}

double total_area(const og::MixedMesh& mesh) {
  const auto q = og::quality_cells(mesh.points, mesh.triangles, mesh.quads);
  double a = 0.0;
  for (const double v : q.area) {
    a += v;
  }
  return a;
}

// Every directed edge at most once; every input segment exactly once, in its direction; no
// interior edge on only one side.
void expect_conforming(const og::TriangulationInput& in, const og::MixedMesh& mesh) {
  std::map<std::pair<MeshIndex, MeshIndex>, int> directed;
  auto add = [&](const MeshIndex* v, int n) {
    for (int k = 0; k < n; ++k) {
      directed[{v[k], v[(k + 1) % n]}] += 1;
    }
  };
  for (const auto& t : mesh.triangles) {
    add(t.data(), 3);
  }
  for (const auto& q : mesh.quads) {
    add(q.data(), 4);
  }
  std::map<std::pair<MeshIndex, MeshIndex>, int> boundary;
  for (const auto& s : in.segments) {
    boundary[{s[0], s[1]}] = 1;
    EXPECT_EQ(directed[std::make_pair(s[0], s[1])], 1) << s[0] << "-" << s[1];
    EXPECT_EQ(directed[std::make_pair(s[1], s[0])], 0);
  }
  for (const auto& [edge, count] : directed) {
    if (count == 0) {
      continue;
    }
    EXPECT_EQ(count, 1);
    if (boundary.count(edge) == 0) {
      EXPECT_EQ(directed[std::make_pair(edge.second, edge.first)], 1) << "hanging edge";
    }
  }
  for (std::size_t k = 0; k < in.points.size(); ++k) {
    EXPECT_EQ(mesh.points[k], in.points[k]);
  }
}

og::TriangulationInput annulus(int n_outer, int n_inner) {
  og::TriangulationInput in;
  add_loop(in, circle(0, 0, 5, n_outer, true));
  add_loop(in, circle(0, 0, 1, n_inner, false));
  return in;
}

double annulus_area(int n_outer, int n_inner) {
  return polygon_area(circle(0, 0, 5, n_outer, true)) -
         polygon_area(circle(0, 0, 1, n_inner, true));
}

}  // namespace

TEST(Recombine, RightTrianglePairsBecomeSquares) {
  constexpr int kN = 8;
  std::vector<Point> pts;
  for (int j = 0; j <= kN; ++j) {
    for (int i = 0; i <= kN; ++i) {
      pts.push_back({static_cast<double>(i), static_cast<double>(j)});
    }
  }
  auto id = [](int i, int j) { return static_cast<MeshIndex>(j) * (kN + 1) + i; };
  std::vector<std::array<MeshIndex, 3>> tris;
  for (int j = 0; j < kN; ++j) {
    for (int i = 0; i < kN; ++i) {
      tris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
      tris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
    }
  }
  const auto out = og::recombine(pts, tris, {});
  EXPECT_TRUE(out.triangles.empty());
  ASSERT_EQ(out.quads.size(), static_cast<std::size_t>(kN * kN));
  for (const auto& q : out.quads) {
    EXPECT_DOUBLE_EQ(og::quad_quality(pts[u(q[0])], pts[u(q[1])], pts[u(q[2])], pts[u(q[3])]), 1.0);
  }
}

TEST(Recombine, QuadDominantAnnulus) {
  const auto in = annulus(80, 32);
  og::DomainMeshOptions options;
  options.element_type = og::ElementType::quad_dominant;
  const auto mesh = og::mesh_domain(in, options);
  expect_conforming(in, mesh);
  EXPECT_NEAR(total_area(mesh), annulus_area(80, 32), 1e-12 * annulus_area(80, 32));
  // Most of the area is in quads, all of them convex.
  const double quad_share = 2.0 * static_cast<double>(mesh.quads.size()) /
                            static_cast<double>(2 * mesh.quads.size() + mesh.triangles.size());
  EXPECT_GE(quad_share, 0.9) << mesh.quads.size() << " quads, " << mesh.triangles.size()
                             << " triangles";
  for (const auto& q : mesh.quads) {
    EXPECT_GT(og::cell_quality(mesh.points, q, 4), 0.0);
  }
  EXPECT_TRUE(mesh.warnings.empty());
}

TEST(AllQuads, AnnulusIsAllQuadConformingAndConvex) {
  for (const auto algorithm :
       {og::TriangleAlgorithm::delaunay, og::TriangleAlgorithm::advancing_front}) {
    const auto in = annulus(80, 32);
    og::DomainMeshOptions options;
    options.algorithm = algorithm;
    options.element_type = og::ElementType::quads;
    const auto mesh = og::mesh_domain(in, options);
    expect_conforming(in, mesh);
    EXPECT_TRUE(mesh.triangles.empty());
    EXPECT_GT(mesh.quads.size(), 100U);
    EXPECT_NEAR(total_area(mesh), annulus_area(80, 32), 1e-12 * annulus_area(80, 32));
    const auto q = og::quality_cells(mesh.points, mesh.triangles, mesh.quads);
    EXPECT_GT(*std::min_element(q.scaled_jacobian.begin(), q.scaled_jacobian.end()), 0.3);
    EXPECT_TRUE(mesh.warnings.empty());
  }
}

TEST(AllQuads, OddLoopIsRejected) {
  const auto in = annulus(81, 32);
  og::DomainMeshOptions options;
  options.element_type = og::ElementType::quads;
  EXPECT_THROW(og::mesh_domain(in, options), std::invalid_argument);
}

TEST(AllQuads, Deterministic) {
  const auto in = annulus(60, 24);
  og::DomainMeshOptions options;
  options.element_type = og::ElementType::quads;
  options.smoothing.method = og::SmoothingMethod::angle_based;
  const auto a = og::mesh_domain(in, options);
  const auto b = og::mesh_domain(in, options);
  EXPECT_EQ(a.points, b.points);
  EXPECT_EQ(a.quads, b.quads);
}

TEST(DomainMesher, TrianglesWithoutSmoothingEqualTheTriangulation) {
  const auto in = annulus(60, 24);
  const auto mesh = og::mesh_domain(in, {});
  const auto tri = og::triangulate_domain(in, {});
  EXPECT_EQ(mesh.points, tri.points);
  EXPECT_EQ(mesh.triangles, tri.triangles);
  EXPECT_TRUE(mesh.quads.empty());
}

TEST(Smoothing, GuardedMethodsImprovePerturbedMesh) {
  for (const auto method : {og::SmoothingMethod::laplacian, og::SmoothingMethod::length_weighted,
                            og::SmoothingMethod::angle_based}) {
    const auto in = annulus(60, 24);
    auto mesh = og::mesh_domain(in, {});
    std::mt19937 rng(3);  // NOLINT(bugprone-random-generator-seed): reproducible test
    std::uniform_real_distribution<double> jitter(-0.04, 0.04);
    for (std::size_t k = in.points.size(); k < mesh.points.size(); ++k) {
      mesh.points[k][0] += jitter(rng);
      mesh.points[k][1] += jitter(rng);
    }
    auto worst = [&]() {
      double w = 2.0;
      for (const auto& t : mesh.triangles) {
        w = std::min(w, og::cell_quality(mesh.points, {t[0], t[1], t[2], -1}, 3));
      }
      return w;
    };
    auto mean_min_angle = [&]() {
      const auto q = og::quality_cells(mesh.points, mesh.triangles, {});
      double s = 0.0;
      for (const double a : q.min_angle) {
        s += a;
      }
      return s / static_cast<double>(q.min_angle.size());
    };
    const double worst_before = worst();
    const double mean_before = mean_min_angle();
    og::SmoothingOptions options;
    options.method = method;
    EXPECT_GT(og::smooth_mesh(mesh.points, mesh.triangles, mesh.quads, in.segments,
                              in.points.size(), options),
              0);
    EXPECT_GE(worst(), worst_before);
    EXPECT_GT(mean_min_angle(), mean_before);
    for (std::size_t k = 0; k < in.points.size(); ++k) {
      EXPECT_EQ(mesh.points[k], in.points[k]);
    }
  }
}

TEST(Untangle, InvertedCellsAreRepaired) {
  const auto in = annulus(60, 24);
  og::DomainMeshOptions options;
  options.element_type = og::ElementType::quads;
  auto mesh = og::mesh_domain(in, options);
  // Throw an interior node far across its neighbours.
  const std::size_t v = in.points.size() + 10;
  mesh.points[v][0] += 0.6;
  mesh.points[v][1] -= 0.4;
  std::int64_t bad = 0;
  for (const auto& q : mesh.quads) {
    if (og::cell_quality(mesh.points, q, 4) <= 0.0) {
      ++bad;
    }
  }
  ASSERT_GT(bad, 0);
  EXPECT_EQ(og::untangle(mesh.points, mesh.triangles, mesh.quads, in.points.size()), 0);
}
