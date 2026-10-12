#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>

#include "curve_fixtures.hpp"
#include "opengrid/core/topo/connector_nodes.hpp"

using og::Vec3;

namespace {

std::shared_ptr<const og::Curve> line(Vec3 a, Vec3 b) {
  return std::make_shared<og::Line>(a, b);
}

std::shared_ptr<const og::Curve> quarter_arc() {
  // From (1, 0) to (0, 1), radius 1.
  return std::make_shared<og::Arc>(Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}, Vec3{1, 0, 0}, 0.0,
                                   std::numbers::pi / 2.0);
}

og::SubConnectorSpec sub(double start, double end, int n, og::DistributionType type,
                         std::optional<double> begin, std::optional<double> finish) {
  og::SubConnectorSpec s;
  s.start = start;
  s.end = end;
  s.n_points = n;
  s.distribution.type = type;
  s.distribution.begin_spacing = begin;
  s.distribution.end_spacing = finish;
  return s;
}

}  // namespace

TEST(CurveChain, LengthAndContinuity) {
  // L-shaped chain: line (2,0)->(1,0) reversed? No: a line from (3,0) to (1,0), then the
  // quarter arc (1,0)->(0,1), then the arc's first half traversed backwards.
  const auto arc = quarter_arc();
  const og::CurveChain chain(
      {{line({3, 0, 0}, {1, 0, 0}), 0.0, 1.0}, {arc, 0.0, 1.0}, {arc, 1.0, 0.5}}, 1e-12);
  const double quarter = std::numbers::pi / 2.0;
  EXPECT_NEAR(chain.length(), 2.0 + quarter + quarter / 2.0, 1e-13);
  EXPECT_LT(og::distance(chain.point_at_length(2.0), {1, 0, 0}), 1e-14);
  EXPECT_LT(og::distance(chain.point_at_length(2.0 + quarter), {0, 1, 0}), 1e-13);
  const Vec3 back = chain.point_at_length(2.0 + quarter + quarter / 4.0);
  EXPECT_NEAR(std::atan2(back.y, back.x), 3.0 * quarter / 4.0, 1e-12);
  EXPECT_LT(og::distance(chain.end(), chain.point_at_length(1e9)), 1e-14);
  // Parametric: u at the junction equals the length fraction.
  EXPECT_LT(og::distance(chain.point_at_param(2.0 / chain.length()), {1, 0, 0}), 1e-13);
  EXPECT_NEAR(chain.length_at_param(0.5), 0.5 * chain.length(), 1e-12);
}

TEST(CurveChain, LocateLength) {
  const auto arc = quarter_arc();
  const og::CurveChain chain({{line({0, 3, 0}, {0, 1, 0}), 0.0, 1.0}, {arc, 1.0, 0.0}}, 1e-12);
  const auto [k0, t0] = chain.locate_length(0.5);
  EXPECT_EQ(k0, 0U);
  EXPECT_NEAR(t0, 0.25, 1e-14);
  // The arc runs backwards: a quarter of its length from its end is t = 0.75.
  const auto [k1, t1] = chain.locate_length(2.0 + std::numbers::pi / 8.0);
  EXPECT_EQ(k1, 1U);
  EXPECT_NEAR(t1, 0.75, 1e-12);
  const double s = 2.3;
  const auto [k2, t2] = chain.locate_length(s);
  EXPECT_LT(og::distance(arc->eval(t2), chain.point_at_length(s)), 1e-15);
  EXPECT_EQ(k2, 1U);
}

TEST(CurveChain, InvalidChains) {
  const auto arc = quarter_arc();
  EXPECT_THROW(og::CurveChain({}, 1e-9), std::invalid_argument);
  EXPECT_THROW(og::CurveChain({{arc, 0.2, 0.2}}, 1e-9), std::invalid_argument);
  EXPECT_THROW(og::CurveChain({{arc, -0.1, 1.0}}, 1e-9), std::invalid_argument);
  try {
    const og::CurveChain chain({{arc, 0.0, 1.0}, {line({0, 1.1, 0}, {0, 3, 0}), 0.0, 1.0}}, 1e-9);
    FAIL() << "expected a gap error";
  } catch (const std::invalid_argument& error) {
    EXPECT_NE(std::string(error.what()).find("segments[1]"), std::string::npos);
  }
}

TEST(ConnectorNodes, EndSpacingsOnCurvedChain) {
  const og::CurveChain chain(
      {{std::make_shared<og::CubicSpline>(og::test::naca0012_upper(101), og::SplineMethod::Clamped,
                                          Vec3{0, 1, 0}, std::nullopt),
        0.0, 1.0}},
      1e-12);
  const double le = 1e-4;
  const double te = 2e-3;
  for (const auto type : {og::DistributionType::Tanh, og::DistributionType::Sine}) {
    const auto nodes =
        og::connector_nodes(chain, og::SpacingMode::ArcLength, {sub(0.0, 1.0, 121, type, le, te)});
    ASSERT_EQ(nodes.points.size(), 121U);
    const auto& s = nodes.arc_positions;
    EXPECT_LT(std::fabs((s[1] - s[0]) - le) / le, 1e-8);
    EXPECT_LT(std::fabs((s[120] - s[119]) - te) / te, 1e-8);
    // Nodes lie on the curve at their arc-length positions.
    for (std::size_t i = 0; i < s.size(); i += 17) {
      EXPECT_LT(og::distance(nodes.points[i], chain.point_at_length(s[i])), 1e-14);
    }
    EXPECT_GT(nodes.max_stretching_ratio, 1.0);
  }
}

TEST(ConnectorNodes, BreakpointsShareNodesAndHonourSpacings) {
  const og::CurveChain chain({{line({0, 0, 0}, {10, 0, 0}), 0.0, 1.0}}, 1e-12);
  const auto nodes =
      og::connector_nodes(chain, og::SpacingMode::ArcLength,
                          {sub(0.0, 0.3, 21, og::DistributionType::Tanh, 0.01, 0.05),
                           sub(0.3, 1.0, 31, og::DistributionType::Geometric, 0.05, std::nullopt)});
  ASSERT_EQ(nodes.points.size(), 51U);
  const auto& s = nodes.arc_positions;
  EXPECT_NEAR(s[20], 3.0, 1e-14);  // breakpoint node index 20
  EXPECT_LT(std::fabs((s[1] - s[0]) - 0.01) / 0.01, 1e-8);
  EXPECT_LT(std::fabs((s[20] - s[19]) - 0.05) / 0.05, 1e-8);
  EXPECT_LT(std::fabs((s[21] - s[20]) - 0.05) / 0.05, 1e-8);
  for (std::size_t i = 1; i < s.size(); ++i) {
    EXPECT_GT(s[i], s[i - 1]);
  }
  EXPECT_TRUE(nodes.warnings.empty());
}

TEST(ConnectorNodes, ParametricModeFollowsCurveParameter) {
  // On a line the parameter is proportional to arc length: both modes agree.
  const og::CurveChain chain({{line({0, 0, 0}, {4, 0, 0}), 0.0, 1.0}}, 1e-12);
  const auto a = og::connector_nodes(chain, og::SpacingMode::ArcLength,
                                     {sub(0.0, 1.0, 11, og::DistributionType::Tanh, 0.05, 0.2)});
  const auto p = og::connector_nodes(chain, og::SpacingMode::Parametric,
                                     {sub(0.0, 1.0, 11, og::DistributionType::Tanh, 0.05, 0.2)});
  for (std::size_t i = 0; i < a.points.size(); ++i) {
    EXPECT_LT(og::distance(a.points[i], p.points[i]), 1e-13);
  }
  // On a polyline with unequal legs (still chord-length parametrised) they agree too; on a
  // NURBS with non-uniform speed they differ.
  const og::CurveChain nurbs({{og::test::rational_cubic(), 0.0, 1.0}}, 1e-12);
  const auto an = og::connector_nodes(
      nurbs, og::SpacingMode::ArcLength,
      {sub(0.0, 1.0, 11, og::DistributionType::Uniform, std::nullopt, std::nullopt)});
  const auto pn = og::connector_nodes(
      nurbs, og::SpacingMode::Parametric,
      {sub(0.0, 1.0, 11, og::DistributionType::Uniform, std::nullopt, std::nullopt)});
  EXPECT_GT(og::distance(an.points[5], pn.points[5]), 1e-3);
  EXPECT_LT(og::distance(an.points.back(), pn.points.back()), 1e-15);
}

TEST(ConnectorNodes, InvalidSubConnectors) {
  const og::CurveChain chain({{line({0, 0, 0}, {1, 0, 0}), 0.0, 1.0}}, 1e-12);
  using og::DistributionType;
  EXPECT_THROW(og::connector_nodes(chain, og::SpacingMode::ArcLength, {}), std::invalid_argument);
  EXPECT_THROW(og::connector_nodes(chain, og::SpacingMode::ArcLength,
                                   {sub(0.0, 0.5, 5, DistributionType::Uniform, {}, {})}),
               std::invalid_argument);
  EXPECT_THROW(og::connector_nodes(chain, og::SpacingMode::ArcLength,
                                   {sub(0.0, 0.5, 5, DistributionType::Uniform, {}, {}),
                                    sub(0.4, 1.0, 5, DistributionType::Uniform, {}, {})}),
               std::invalid_argument);
  EXPECT_THROW(og::connector_nodes(chain, og::SpacingMode::ArcLength,
                                   {sub(0.0, 1.0, 1, DistributionType::Uniform, {}, {})}),
               std::invalid_argument);
}

TEST(ConnectorNodes, InfeasibleSpacingWarnsButStaysValid) {
  const og::CurveChain chain({{line({0, 0, 0}, {1, 0, 0}), 0.0, 1.0}}, 1e-12);
  const auto nodes = og::connector_nodes(chain, og::SpacingMode::ArcLength,
                                         {sub(0.0, 1.0, 11, og::DistributionType::Tanh, 0.7, 0.7)});
  EXPECT_FALSE(nodes.warnings.empty());
  for (std::size_t i = 1; i < nodes.arc_positions.size(); ++i) {
    EXPECT_GT(nodes.arc_positions[i], nodes.arc_positions[i - 1]);
  }
}
