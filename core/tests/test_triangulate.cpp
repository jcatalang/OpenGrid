#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <numbers>
#include <random>
#include <set>
#include <stdexcept>

#include "opengrid/core/unstructured/triangulate.hpp"
#include "opengrid/core/util/predicates.hpp"

using og::MeshIndex;

namespace {

using Int = __int128;

int sign(Int v) {
  if (v > 0) {
    return 1;
  }
  if (v < 0) {
    return -1;
  }
  return 0;
}

// Closed polygon loop appended to input; `ccw` loops keep the domain on the left.
void add_loop(og::TriangulationInput& in, const std::vector<std::array<double, 2>>& pts) {
  const auto first = static_cast<MeshIndex>(in.points.size());
  for (const auto& q : pts) {
    in.points.push_back(q);
  }
  const auto n = static_cast<MeshIndex>(pts.size());
  for (MeshIndex k = 0; k < n; ++k) {
    in.segments.push_back({first + k, first + (k + 1) % n});
  }
}

// Square [x0, x1] x [y0, y1] with m segments per side, counter-clockwise (or clockwise).
std::vector<std::array<double, 2>> square(double x0, double y0, double x1, double y1, int m,
                                          bool ccw) {
  std::vector<std::array<double, 2>> pts;
  pts.reserve(4 * static_cast<std::size_t>(m));
  for (int k = 0; k < m; ++k) {
    pts.push_back({x0 + (x1 - x0) * k / m, y0});
  }
  for (int k = 0; k < m; ++k) {
    pts.push_back({x1, y0 + (y1 - y0) * k / m});
  }
  for (int k = 0; k < m; ++k) {
    pts.push_back({x1 - (x1 - x0) * k / m, y1});
  }
  for (int k = 0; k < m; ++k) {
    pts.push_back({x0, y1 - (y1 - y0) * k / m});
  }
  if (!ccw) {
    std::reverse(pts.begin() + 1, pts.end());
  }
  return pts;
}

double area(const og::TriangleMesh& mesh, std::size_t t) {
  const auto& a = mesh.points[static_cast<std::size_t>(mesh.triangles[t][0])];
  const auto& b = mesh.points[static_cast<std::size_t>(mesh.triangles[t][1])];
  const auto& c = mesh.points[static_cast<std::size_t>(mesh.triangles[t][2])];
  return 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
}

double min_angle_deg(const og::TriangleMesh& mesh, std::size_t t) {
  double worst = 180.0;
  for (int k = 0; k < 3; ++k) {
    const auto& a =
        mesh.points[static_cast<std::size_t>(mesh.triangles[t][static_cast<std::size_t>(k)])];
    const auto& b = mesh.points[static_cast<std::size_t>(
        mesh.triangles[t][static_cast<std::size_t>((k + 1) % 3)])];
    const auto& c = mesh.points[static_cast<std::size_t>(
        mesh.triangles[t][static_cast<std::size_t>((k + 2) % 3)])];
    const double ux = b[0] - a[0];
    const double uy = b[1] - a[1];
    const double vx = c[0] - a[0];
    const double vy = c[1] - a[1];
    worst = std::min(worst, std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy) * 180.0 /
                                std::numbers::pi);
  }
  return worst;
}

// Every input segment is an edge of exactly one triangle (on its left), no triangle is
// inverted, and every interior edge is shared by exactly two triangles.
void expect_conforming(const og::TriangulationInput& in, const og::TriangleMesh& mesh) {
  std::map<std::pair<MeshIndex, MeshIndex>, int> directed;
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    EXPECT_GT(area(mesh, t), 0.0) << "triangle " << t;
    for (int k = 0; k < 3; ++k) {
      const auto a = mesh.triangles[t][static_cast<std::size_t>(k)];
      const auto b = mesh.triangles[t][static_cast<std::size_t>((k + 1) % 3)];
      directed[{a, b}] += 1;
    }
  }
  std::set<std::pair<MeshIndex, MeshIndex>> boundary;
  for (const auto& s : in.segments) {
    EXPECT_EQ(directed[std::make_pair(s[0], s[1])], 1) << s[0] << "-" << s[1];
    EXPECT_EQ(directed.count(std::make_pair(s[1], s[0])), 0U) << s[0] << "-" << s[1];
    boundary.insert({s[0], s[1]});
  }
  for (const auto& [edge, count] : directed) {
    EXPECT_EQ(count, 1);
    if (boundary.count(edge) == 0) {
      EXPECT_EQ(directed.count(std::make_pair(edge.second, edge.first)), 1U)
          << "open interior edge";
    }
  }
  // Boundary points keep their indices and coordinates.
  for (std::size_t i = 0; i < in.points.size(); ++i) {
    EXPECT_EQ(mesh.points[i], in.points[i]);
  }
}

}  // namespace

TEST(Predicates, Orient2dExactOnNearCollinearPoints) {
  // Coordinates X * 2^-40 with integer X: exact orientation from 128-bit integers.
  std::mt19937_64 rng(7);  // NOLINT(bugprone-random-generator-seed): reproducible test
  const long long range = 35'184'372'088'832LL;  // 2^45
  std::uniform_int_distribution<long long> big(-range, range);
  std::uniform_int_distribution<long long> tiny(-3, 3);
  const double unit = std::ldexp(1.0, -40);
  auto point = [unit](long long x, long long y) {
    return std::array<double, 2>{static_cast<double>(x) * unit, static_cast<double>(y) * unit};
  };
  for (int trial = 0; trial < 20000; ++trial) {
    const long long ax = big(rng);
    const long long ay = big(rng);
    const long long dx = big(rng) / 1024;
    const long long dy = big(rng) / 1024;
    const long long cx = ax + 2 * dx + tiny(rng);
    const long long cy = ay + 2 * dy + tiny(rng);
    const long long bx = ax + dx;
    const long long by = ay + dy;
    const auto a = point(ax, ay);
    const auto b = point(bx, by);
    const auto c = point(cx, cy);
    const Int det = static_cast<Int>(ax - cx) * (by - cy) - static_cast<Int>(ay - cy) * (bx - cx);
    ASSERT_EQ(og::orient2d(a.data(), b.data(), c.data()), sign(det)) << trial;
  }
}

TEST(Predicates, IncircleExactOnCocircularPoints) {
  // Lattice points on x^2 + y^2 = 25 are exactly cocircular; tiny shifts decide the side.
  const std::array<double, 2> a = {5, 0};
  const std::array<double, 2> b = {0, 5};
  const std::array<double, 2> c = {-3, -4};
  for (const auto& d : std::vector<std::array<double, 2>>{{3, 4}, {-4, 3}, {4, -3}}) {
    EXPECT_EQ(og::incircle(a.data(), b.data(), c.data(), d.data()), 0);
    const std::array<double, 2> inside = {d[0] * (1.0 - 1e-15), d[1] * (1.0 - 1e-15)};
    const std::array<double, 2> outside = {d[0] * (1.0 + 1e-15), d[1] * (1.0 + 1e-15)};
    EXPECT_EQ(og::incircle(a.data(), b.data(), c.data(), inside.data()), 1);
    EXPECT_EQ(og::incircle(a.data(), b.data(), c.data(), outside.data()), -1);
  }
  const std::array<double, 2> a2 = {1e15 + 5, 0};
  const std::array<double, 2> b2 = {1e15, 5};
  const std::array<double, 2> c2 = {1e15 - 3, -4};
  const std::array<double, 2> d2 = {1e15 + 3, 4};
  EXPECT_EQ(og::incircle(a2.data(), b2.data(), c2.data(), d2.data()), 0);
}

TEST(SizeField, MatchesBruteForce) {
  std::mt19937_64 rng(3);  // NOLINT(bugprone-random-generator-seed): reproducible test
  std::uniform_real_distribution<double> coord(-5.0, 5.0);
  std::uniform_real_distribution<double> h(0.01, 0.5);
  std::vector<std::array<double, 2>> pts;
  std::vector<double> sizes;
  for (int k = 0; k < 500; ++k) {
    pts.push_back({coord(rng), coord(rng)});
    sizes.push_back(h(rng));
  }
  og::SizeFieldOptions options;
  options.growth_rate = 1.15;
  options.max_size = 2.0;
  options.sources.push_back({{0, 0}, {1, 1}, 0.005, 0.3});
  const og::SizeField field(pts, sizes, options);
  for (int k = 0; k < 2000; ++k) {
    const double x = 3.0 * coord(rng);
    const double y = 3.0 * coord(rng);
    double expected = 2.0;
    for (std::size_t j = 0; j < pts.size(); ++j) {
      expected = std::min(expected, sizes[j] + 0.15 * std::hypot(x - pts[j][0], y - pts[j][1]));
    }
    // Source: segment (0,0)-(1,1).
    const double t = std::clamp((x + y) / 2.0, 0.0, 1.0);
    const double ds = std::hypot(x - t, y - t);
    expected = std::min(expected, 0.005 + 0.15 * std::max(0.0, ds - 0.3));
    EXPECT_NEAR(field(x, y), expected, 1e-14);
  }
}

TEST(Triangulate, SquareWithoutRefinement) {
  og::TriangulationInput in;
  add_loop(in, square(0, 0, 1, 1, 10, true));
  og::RefinementOptions options;
  options.refine = false;
  const auto mesh = og::triangulate_domain(in, options);
  expect_conforming(in, mesh);
  EXPECT_EQ(mesh.points.size(), 40U);
  EXPECT_EQ(mesh.triangles.size(), 2U * 40 - 2 - 40);  // T = 2V - 2 - B (all on the boundary)
  double total = 0.0;
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    total += area(mesh, t);
  }
  EXPECT_NEAR(total, 1.0, 1e-14);
}

TEST(Triangulate, RefinedSquareWithHoleQualityAndArea) {
  og::TriangulationInput in;
  add_loop(in, square(0, 0, 4, 4, 40, true));
  add_loop(in, square(1.5, 1.5, 2.5, 2.5, 20, false));  // hole, clockwise
  og::RefinementOptions options;
  options.min_angle_deg = 25.0;
  const auto mesh = og::triangulate_domain(in, options);
  expect_conforming(in, mesh);
  double total = 0.0;
  double worst = 180.0;
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    total += area(mesh, t);
    worst = std::min(worst, min_angle_deg(mesh, t));
  }
  EXPECT_NEAR(total, 16.0 - 1.0, 1e-12 * 15.0);
  EXPECT_GE(worst, 20.7);
  EXPECT_GT(mesh.points.size(), 160U + 80U);  // interior points were added
  for (std::size_t i = 160 + 80; i < mesh.points.size(); ++i) {
    const auto& q = mesh.points[i];
    EXPECT_FALSE(q[0] > 1.5 && q[0] < 2.5 && q[1] > 1.5 && q[1] < 2.5) << "point in the hole";
  }
}

namespace {

og::TriangulationInput l_shape(int coarse, int fine) {
  og::TriangulationInput in;
  std::vector<std::array<double, 2>> pts;
  auto edge = [&pts](double x0, double y0, double x1, double y1, int m) {
    for (int k = 0; k < m; ++k) {
      pts.push_back({x0 + (x1 - x0) * k / m, y0 + (y1 - y0) * k / m});
    }
  };
  edge(0, 0, 2, 0, 2 * coarse);
  edge(2, 0, 2, 1, coarse);
  edge(2, 1, 1, 1, fine);
  edge(1, 1, 1, 2, fine);
  edge(1, 2, 0, 2, coarse);
  edge(0, 2, 0, 0, 2 * coarse);
  add_loop(in, pts);
  return in;
}

}  // namespace

TEST(Triangulate, NonConvexDomainAndGrading) {
  // Consistently spaced boundary (0.02); the interior grades up to the maximum size.
  const auto in = l_shape(50, 50);
  og::RefinementOptions options;
  options.size.growth_rate = 1.2;
  options.size.max_size = 0.2;
  const auto mesh = og::triangulate_domain(in, options);
  expect_conforming(in, mesh);
  double total = 0.0;
  double worst = 180.0;
  double largest = 0.0;
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    total += area(mesh, t);
    largest = std::max(largest, area(mesh, t));
    worst = std::min(worst, min_angle_deg(mesh, t));
  }
  EXPECT_NEAR(total, 3.0, 3e-12);
  EXPECT_GE(worst, 20.7);
  // Graded: the largest triangles are much larger than the boundary triangles (~0.02^2/2).
  EXPECT_GT(largest, 20.0 * 0.0002);
}

TEST(Triangulate, AbruptBoundarySpacingJumpsTerminate) {
  // Tenfold spacing jumps at the re-entrant edges: with fixed boundary nodes no triangulation
  // can have good angles next to the coarse segments there (documented limitation). The mesh
  // must still be conforming, exact in area, terminate, and stay good almost everywhere.
  const auto in = l_shape(10, 100);
  og::RefinementOptions options;
  options.size.growth_rate = 1.1;
  const auto mesh = og::triangulate_domain(in, options);
  expect_conforming(in, mesh);
  double total = 0.0;
  std::size_t bad = 0;
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    total += area(mesh, t);
    if (min_angle_deg(mesh, t) < 20.7) {
      bad += 1;
    }
  }
  EXPECT_NEAR(total, 3.0, 3e-12);
  EXPECT_TRUE(mesh.warnings.empty());
  EXPECT_LT(bad, mesh.triangles.size() / 20);
}

TEST(Triangulate, Deterministic) {
  og::TriangulationInput in;
  add_loop(in, square(0, 0, 3, 1, 30, true));
  add_loop(in, square(1, 0.4, 1.2, 0.6, 8, false));
  const auto a = og::triangulate_domain(in, {});
  const auto b = og::triangulate_domain(in, {});
  EXPECT_EQ(a.points, b.points);
  EXPECT_EQ(a.triangles, b.triangles);
}

TEST(Triangulate, InvalidInput) {
  og::TriangulationInput dup;
  add_loop(dup, {{0, 0}, {1, 0}, {1, 1}, {1, 0}});
  EXPECT_THROW(og::triangulate_domain(dup, {}), std::invalid_argument);
  og::TriangulationInput clockwise;
  add_loop(clockwise, square(0, 0, 1, 1, 4, false));
  EXPECT_THROW(og::triangulate_domain(clockwise, {}), std::invalid_argument);
  og::TriangulationInput crossing;
  add_loop(crossing, {{0, 0}, {2, 0}, {0, 2}, {2, 2}});  // bow tie
  EXPECT_THROW(og::triangulate_domain(crossing, {}), std::invalid_argument);
  og::TriangulationInput bad;
  add_loop(bad, square(0, 0, 1, 1, 4, true));
  bad.segments.push_back({0, 99});
  EXPECT_THROW(og::triangulate_domain(bad, {}), std::invalid_argument);
}

TEST(Triangulate, CollinearBoundaryPointsOnLongEdges) {
  // A segment passing exactly through another boundary vertex must not be "recovered" across it.
  og::TriangulationInput in;
  add_loop(in, {{0, 0}, {4, 0}, {4, 1}, {0, 1}});
  og::RefinementOptions options;
  options.size.max_size = 0.2;
  const auto mesh = og::triangulate_domain(in, options);
  expect_conforming(in, mesh);
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    EXPECT_GT(area(mesh, t), 0.0);
  }
}
