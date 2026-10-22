#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <utility>
#include <vector>

#include "opengrid/core/quality/quality_cells.hpp"
#include "opengrid/core/unstructured/domain_mesher.hpp"
#include "opengrid/core/unstructured/smoothing.hpp"
#include "opengrid/core/util/predicates.hpp"

namespace {

using og::MeshIndex;
using Point = std::array<double, 2>;

std::size_t u(MeshIndex i) {
  return static_cast<std::size_t>(i);
}

struct Case {
  og::TriangulationInput in;
  std::vector<bool> wall;

  void add(const std::vector<Point>& pts, bool is_wall) {
    const auto first = static_cast<MeshIndex>(in.points.size());
    const auto n = static_cast<MeshIndex>(pts.size());
    for (const auto& p : pts) {
      in.points.push_back(p);
    }
    for (MeshIndex k = 0; k < n; ++k) {
      in.segments.push_back({first + k, first + (k + 1) % n});
      wall.push_back(is_wall);
    }
  }
};

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

double naca(double x) {
  return 0.6 * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x -
                0.1036 * x * x * x * x);
}

// NACA 0012 (sharp trailing edge), or with the trailing edge cut at x = cut (blunt), clockwise,
// scaled and shifted.
std::vector<Point> airfoil(int n, double cut, double scale, double dx, double dy, double angle) {
  std::vector<Point> upper;
  std::vector<Point> lower;
  for (int i = 0; i <= n; ++i) {
    const double x = cut * 0.5 * (1.0 - std::cos(std::numbers::pi * i / n));
    upper.push_back({x, naca(x)});
    lower.push_back({x, -naca(x)});
  }
  // Clockwise: leading edge -> upper surface -> trailing edge -> lower surface.
  std::vector<Point> pts;
  for (int i = 0; i <= n; ++i) {
    pts.push_back(upper[u(i)]);
  }
  if (cut < 1.0) {
    const int base = 6;  // points on the blunt base, top to bottom (excluding both ends)
    for (int k = 1; k < base; ++k) {
      const double t = static_cast<double>(k) / base;
      pts.push_back({cut, naca(cut) * (1.0 - 2.0 * t)});
    }
  }
  int start = n;
  if (cut >= 1.0) {
    start = n - 1;  // the sharp trailing edge is shared
  }
  for (int i = start; i >= 1; --i) {
    pts.push_back(lower[u(i)]);
  }
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  for (auto& p : pts) {
    p = {dx + scale * (c * p[0] + s * p[1]), dy + scale * (-s * p[0] + c * p[1])};
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

// Conformity and validity: every directed edge at most once, interior edges matched, boundary
// segments once, input points preserved, every cell positively oriented and convex.
void expect_valid(const Case& c, const og::MixedMesh& mesh) {
  std::map<std::pair<MeshIndex, MeshIndex>, int> directed;
  auto add = [&](const MeshIndex* v, int n) {
    for (int k = 0; k < n; ++k) {
      directed[{v[k], v[(k + 1) % n]}] += 1;
    }
  };
  for (const auto& t : mesh.triangles) {
    add(t.data(), 3);
    EXPECT_GT(og::cell_quality(mesh.points, {t[0], t[1], t[2], -1}, 3), 0.0);
  }
  for (const auto& q : mesh.quads) {
    add(q.data(), 4);
    EXPECT_GT(og::cell_quality(mesh.points, q, 4), 0.0);
  }
  std::map<std::pair<MeshIndex, MeshIndex>, int> boundary;
  for (const auto& s : c.in.segments) {
    boundary[{s[0], s[1]}] = 1;
    EXPECT_EQ(directed[std::make_pair(s[0], s[1])], 1);
  }
  for (const auto& [edge, count] : directed) {
    if (count == 0) {
      continue;
    }
    EXPECT_EQ(count, 1);
    if (boundary.count(edge) == 0) {
      EXPECT_EQ(directed[std::make_pair(edge.second, edge.first)], 1);
    }
  }
  for (std::size_t k = 0; k < c.in.points.size(); ++k) {
    EXPECT_EQ(mesh.points[k], c.in.points[k]);
  }
}

double total_area(const og::MixedMesh& mesh) {
  const auto q = og::quality_cells(mesh.points, mesh.triangles, mesh.quads);
  double a = 0.0;
  for (const double v : q.area) {
    a += v;
  }
  return a;
}

}  // namespace

TEST(Inflation, StraightWallHeightsAndGrowthAreExact) {
  // Channel [0, 4] x [0, 1]; the bottom is a wall, the rest is not. Away from the tapered ends
  // the columns are straight and layer k is exactly h0 g^k thick.
  Case c;
  std::vector<Point> pts;
  pts.reserve(200);
  const int n = 80;
  for (int i = 0; i < n; ++i) {
    pts.push_back({4.0 * i / n, 0.0});
  }
  const std::size_t wall_count = pts.size();
  for (int j = 0; j < 20; ++j) {
    pts.push_back({4.0, j / 20.0});
  }
  for (int i = n; i > 0; --i) {
    pts.push_back({4.0 * i / n, 1.0});
  }
  for (int j = 20; j > 0; --j) {
    pts.push_back({0.0, j / 20.0});
  }
  c.in.points = pts;
  for (std::size_t k = 0; k < pts.size(); ++k) {
    c.in.segments.push_back(
        {static_cast<MeshIndex>(k), static_cast<MeshIndex>((k + 1) % pts.size())});
    c.wall.push_back(k < wall_count);
  }
  og::DomainMeshOptions options;
  options.wall = c.wall;
  options.inflation.first_height = 1e-3;
  options.inflation.growth_rate = 1.25;
  options.inflation.max_layers = 8;
  options.inflation.transition = false;
  const auto mesh = og::mesh_domain(c.in, options);
  expect_valid(c, mesh);
  EXPECT_EQ(mesh.layers, 8);
  EXPECT_NEAR(total_area(mesh), 4.0, 1e-12);
  // Column above x = 2 (wall node 40): heights of its nodes.
  std::vector<double> ys;
  for (std::size_t k = c.in.points.size(); k < mesh.points.size(); ++k) {
    if (std::fabs(mesh.points[k][0] - 2.0) < 1e-12 && mesh.points[k][1] < 0.5) {
      ys.push_back(mesh.points[k][1]);
    }
  }
  std::sort(ys.begin(), ys.end());
  ASSERT_EQ(ys.size(), 8U);
  double previous = 0.0;
  double h = 1e-3;
  for (const double y : ys) {
    EXPECT_NEAR((y - previous) / h, 1.0, 1e-6);
    previous = y;
    h *= 1.25;
  }
}

// Channel [0, 4] x [0, 1] with a wall at the bottom and, on the inlet (x = 0) and outlet
// (x = 4), vertices at the layer heights up to the layers' thickness, then coarser.
Case layered_channel(int layers, double first, double growth, double inlet_deg = 90.0) {
  Case c;
  const int n = 80;
  std::vector<double> ys;  // heights on the side walls, from the bottom
  double y = 0.0;
  double h = first;
  for (int k = 0; k < layers; ++k) {
    y += h;
    ys.push_back(y);
    h *= growth;
  }
  const double top = ys.back();
  for (int j = 1; j <= 12; ++j) {
    ys.push_back(top + (1.0 - top) * j / 13.0);
  }
  std::vector<Point> pts;
  for (int i = 0; i < n; ++i) {
    pts.push_back({4.0 * i / n, 0.0});
  }
  const std::size_t wall_count = pts.size();
  pts.push_back({4.0, 0.0});
  for (const double v : ys) {
    pts.push_back({4.0, v});
  }
  pts.push_back({4.0, 1.0});
  for (int i = n - 1; i > 0; --i) {
    pts.push_back({4.0 * i / n, 1.0});
  }
  if (inlet_deg < 90.0) {  // (the top must stay to the right of the leaning inlet's end)
    while (pts.back()[0] <= 1.0 / std::tan(inlet_deg * std::numbers::pi / 180.0) + 0.05) {
      pts.pop_back();
    }
  }
  // The inlet leans at inlet_deg to the wall; its vertices keep their heights across the
  // layers (so along it they are 1 / sin farther apart).
  double lean = 0.0;  // (exactly: 1 / tan(90 degrees) is 6e-17)
  if (inlet_deg != 90.0) {
    lean = 1.0 / std::tan(inlet_deg * std::numbers::pi / 180.0);
  }
  pts.push_back({lean, 1.0});
  for (auto it = ys.rbegin(); it != ys.rend(); ++it) {
    pts.push_back({lean * *it, *it});
  }
  c.in.points = pts;
  for (std::size_t k = 0; k < pts.size(); ++k) {
    c.in.segments.push_back(
        {static_cast<MeshIndex>(k), static_cast<MeshIndex>((k + 1) % pts.size())});
    c.wall.push_back(k < wall_count);
  }
  return c;
}

TEST(Inflation, SlidingEndsRunTheLayersIntoTheNeighbouringBoundary) {
  const int layers = 8;
  Case c = layered_channel(layers, 1e-3, 1.25);
  og::DomainMeshOptions options;
  options.wall = c.wall;
  options.inflation.first_height = 1e-3;
  options.inflation.growth_rate = 1.25;
  options.inflation.max_layers = layers;
  options.inflation.transition = false;
  // Tapered ends (the default): fewer layer cells than a full stack over the 80 wall segments.
  const auto tapered = og::mesh_domain(c.in, options);
  expect_valid(c, tapered);
  EXPECT_LT(tapered.quads.size(), std::size_t{80 * layers});
  // Sliding ends: every layer runs from the inlet to the outlet, on their own vertices (the
  // input points are all kept, every boundary segment is a cell edge).
  options.inflation.slide = true;
  const auto slid = og::mesh_domain(c.in, options);
  expect_valid(c, slid);
  EXPECT_EQ(slid.layers, layers);
  EXPECT_EQ(slid.quads.size(), std::size_t{80 * layers});
  EXPECT_NEAR(total_area(slid), 4.0, 1e-12);
  // The layer next to the inlet: its cells' left sides are the inlet's segments.
  int on_inlet = 0;
  for (const auto& q : slid.quads) {
    int count = 0;
    for (const MeshIndex v : q) {
      count += static_cast<int>(slid.points[u(v)][0] == 0.0);
    }
    on_inlet += static_cast<int>(count == 2);
  }
  EXPECT_EQ(on_inlet, layers);
  // Vertices too far from the layer heights: the end column stops, the layers taper as before.
  Case coarse = layered_channel(layers, 3e-3, 1.25);
  options.wall = coarse.wall;
  const auto stopped = og::mesh_domain(coarse.in, options);
  expect_valid(coarse, stopped);
  EXPECT_LT(stopped.quads.size(), std::size_t{80 * layers});
}

TEST(Inflation, SlidingNeedsAWallAngleBetween45And135Degrees) {
  // At 60 degrees the end column slides along the leaning inlet (a full layer stack); at 30
  // degrees layers converging into the wedge would be squeezed, so that end tapers.
  for (const auto& [angle, slides] : {std::pair{60.0, true}, std::pair{30.0, false}}) {
    Case c = layered_channel(8, 1e-3, 1.25, angle);
    og::DomainMeshOptions options;
    options.wall = c.wall;
    options.inflation.first_height = 1e-3;
    options.inflation.growth_rate = 1.25;
    options.inflation.max_layers = 8;
    options.inflation.transition = false;
    options.inflation.slide = true;
    const auto mesh = og::mesh_domain(c.in, options);
    expect_valid(c, mesh);
    if (slides) {
      EXPECT_EQ(mesh.quads.size(), std::size_t{80 * 8}) << angle;
    } else {
      EXPECT_LT(mesh.quads.size(), std::size_t{80 * 8}) << angle;
    }
  }
}

TEST(Inflation, ClusteredInletNextToACoarseWallRefinesToAnEnd) {
  // The channel's inlet and outlet clustered at the bottom (1e-3) next to a bottom of 0.05
  // segments: the thin triangles between them cannot be repaired without splitting the bottom
  // segments, and their centroid splits used to cascade until the point limit (hours).
  Case c = layered_channel(8, 1e-3, 1.25);
  og::DomainMeshOptions options;  // no inflation: the bare triangulation
  const auto mesh = og::mesh_domain(c.in, options);
  expect_valid(c, mesh);
  EXPECT_LT(mesh.points.size(), std::size_t{4000});
  EXPECT_NEAR(total_area(mesh), 4.0, 1e-12);
}

TEST(Inflation, CylinderFirstHeightTransitionAndCore) {
  for (const auto element : {og::ElementType::triangles, og::ElementType::quads}) {
    Case c;
    c.add(circle(0, 0, 20, 120, true), false);
    c.add(circle(0, 0, 0.5, 200, false), true);
    og::DomainMeshOptions options;
    options.element_type = element;
    options.wall = c.wall;
    options.inflation.first_height = 1e-4;
    options.inflation.growth_rate = 1.2;
    options.inflation.max_layers = 60;
    options.refinement.size.growth_rate = 1.15;
    const auto mesh = og::mesh_domain(c.in, options);
    expect_valid(c, mesh);
    EXPECT_TRUE(mesh.warnings.empty());
    const double area =
        polygon_area(circle(0, 0, 20, 120, true)) - polygon_area(circle(0, 0, 0.5, 200, true));
    EXPECT_NEAR(total_area(mesh), area, 1e-11 * area);
    // First layer: the first 200 new points lie 1e-4 off the wall polygon (< 2 % curved wall).
    for (std::size_t k = c.in.points.size(); k < c.in.points.size() + 200; ++k) {
      double d = 1e300;
      for (std::size_t s = 120; s < c.in.segments.size(); ++s) {
        const auto& a = c.in.points[u(c.in.segments[s][0])];
        const auto& b = c.in.points[u(c.in.segments[s][1])];
        const double dx = b[0] - a[0];
        const double dy = b[1] - a[1];
        const double t =
            std::clamp(((mesh.points[k][0] - a[0]) * dx + (mesh.points[k][1] - a[1]) * dy) /
                           (dx * dx + dy * dy),
                       0.0, 1.0);
        d = std::min(
            d, std::hypot(mesh.points[k][0] - a[0] - t * dx, mesh.points[k][1] - a[1] - t * dy));
      }
      EXPECT_NEAR(d / 1e-4, 1.0, 0.02);
    }
    // Transition: layers stop at the core size (well before 60 layers).
    EXPECT_GT(mesh.layers, 10);
    EXPECT_LT(mesh.layers, 60);
    if (element == og::ElementType::quads) {
      EXPECT_TRUE(mesh.triangles.empty());
    }
  }
}

TEST(Inflation, SharpAndBluntTrailingEdges) {
  for (const double cut : {1.0, 0.98}) {
    Case c;
    c.add(circle(0.5, 0, 20, 120, true), false);
    c.add(airfoil(100, cut, 1.0, 0.0, 0.0, 0.0), true);
    og::DomainMeshOptions options;
    options.wall = c.wall;
    options.inflation.first_height = 1e-4;
    options.refinement.size.growth_rate = 1.15;
    const auto mesh = og::mesh_domain(c.in, options);
    expect_valid(c, mesh);
    EXPECT_TRUE(mesh.warnings.empty());
    EXPECT_GT(mesh.layers, 15);
    // Fan around the trailing edge: layer nodes behind it (x > TE) close to the wake line.
    int behind = 0;
    for (std::size_t k = c.in.points.size(); k < mesh.points.size(); ++k) {
      if (mesh.points[k][0] > cut + 1e-4 && std::fabs(mesh.points[k][1]) < 2e-3) {
        ++behind;
      }
    }
    EXPECT_GT(behind, 5) << "cut " << cut;
  }
}

TEST(Inflation, TwoElementGapTerminatesOnProximity) {
  // Main element and a flap below/behind its trailing edge with a small gap.
  Case c;
  c.add(circle(0.6, 0, 20, 120, true), false);
  c.add(airfoil(80, 1.0, 1.0, 0.0, 0.0, 0.0), true);
  c.add(airfoil(60, 1.0, 0.3, 0.97, -0.045, 20.0 * std::numbers::pi / 180.0), true);
  og::DomainMeshOptions options;
  options.wall = c.wall;
  options.inflation.first_height = 1e-4;
  options.refinement.size.growth_rate = 1.15;
  const auto mesh = og::mesh_domain(c.in, options);
  expect_valid(c, mesh);
  EXPECT_TRUE(mesh.warnings.empty());
  // The layers of both elements stop short of each other: core triangles fill the gap between
  // the main element's trailing-edge region and the flap.
  int core_in_gap = 0;
  for (const auto& t : mesh.triangles) {
    double x = 0.0;
    double y = 0.0;
    for (const MeshIndex v : t) {
      x += mesh.points[u(v)][0] / 3.0;
      y += mesh.points[u(v)][1] / 3.0;
    }
    if (x > 0.975 && x < 0.985 && y > -0.035 && y < -0.012) {
      ++core_in_gap;
    }
  }
  EXPECT_GT(core_in_gap, 3);
}

TEST(Inflation, LastLayerMatchesCoreSize) {
  // Straight wall with spacing 0.05 and core size 0.05 everywhere: the last layer is within a
  // growth factor (20 %) of the core size.
  Case c;
  std::vector<Point> pts;
  pts.reserve(200);
  for (int i = 0; i < 80; ++i) {
    pts.push_back({0.05 * i, 0.0});
  }
  const std::size_t wall_count = pts.size();
  for (int j = 0; j < 20; ++j) {
    pts.push_back({4.0, 0.05 * j});
  }
  for (int i = 80; i > 0; --i) {
    pts.push_back({0.05 * i, 1.0});
  }
  for (int j = 20; j > 0; --j) {
    pts.push_back({0.0, 0.05 * j});
  }
  c.in.points = pts;
  for (std::size_t k = 0; k < pts.size(); ++k) {
    c.in.segments.push_back(
        {static_cast<MeshIndex>(k), static_cast<MeshIndex>((k + 1) % pts.size())});
    c.wall.push_back(k < wall_count);
  }
  og::DomainMeshOptions options;
  options.wall = c.wall;
  options.inflation.first_height = 1e-4;
  options.inflation.growth_rate = 1.2;
  options.inflation.max_layers = 100;
  options.refinement.size.max_size = 0.05;
  const auto mesh = og::mesh_domain(c.in, options);
  expect_valid(c, mesh);
  std::vector<double> ys;
  for (std::size_t k = c.in.points.size(); k < mesh.points.size(); ++k) {
    if (std::fabs(mesh.points[k][0] - 2.0) < 1e-12 && mesh.points[k][1] < 0.5) {
      ys.push_back(mesh.points[k][1]);
    }
  }
  std::sort(ys.begin(), ys.end());
  ASSERT_GE(ys.size(), 2U);
  EXPECT_EQ(static_cast<int>(ys.size()), mesh.layers);
  const double last = ys.back() - ys[ys.size() - 2];
  EXPECT_GE(last, 0.05 / 1.2 * 0.999);
  EXPECT_LE(last, 0.05);
}

TEST(Inflation, InvalidOptionsAndFlags) {
  Case c;
  c.add(circle(0, 0, 5, 40, true), false);
  c.add(circle(0, 0, 1, 20, false), true);
  og::DomainMeshOptions options;
  options.wall = c.wall;
  options.wall.pop_back();
  EXPECT_THROW(og::mesh_domain(c.in, options), std::invalid_argument);
  options.wall = c.wall;
  options.inflation.first_height = -1.0;
  EXPECT_THROW(og::mesh_domain(c.in, options), std::invalid_argument);
}

TEST(Inflation, Deterministic) {
  Case c;
  c.add(circle(0.5, 0, 20, 120, true), false);
  c.add(airfoil(60, 1.0, 1.0, 0.0, 0.0, 0.0), true);
  og::DomainMeshOptions options;
  options.wall = c.wall;
  options.inflation.first_height = 1e-4;
  const auto a = og::mesh_domain(c.in, options);
  const auto b = og::mesh_domain(c.in, options);
  EXPECT_EQ(a.points, b.points);
  EXPECT_EQ(a.triangles, b.triangles);
  EXPECT_EQ(a.quads, b.quads);
}
