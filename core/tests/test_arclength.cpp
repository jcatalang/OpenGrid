#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "curve_fixtures.hpp"

using og::Vec3;

namespace {

// Independent arc-length reference using only eval(): Richardson-extrapolated chord sums
// (chord-sum error is O(N^-2) per smooth segment, extrapolation removes it to O(N^-4)).
double reference_length(const og::Curve& curve) {
  const std::vector<double> bps = curve.breakpoints();
  long double total = 0.0L;
  for (std::size_t i = 0; i + 1 < bps.size(); ++i) {
    auto chord_sum = [&](int n) {
      long double sum = 0.0L;
      Vec3 prev = curve.eval(bps[i]);
      for (int k = 1; k <= n; ++k) {
        const double t = bps[i] + (bps[i + 1] - bps[i]) * k / n;
        const Vec3 x = curve.eval(t);
        sum += og::distance(x, prev);
        prev = x;
      }
      return sum;
    };
    const long double coarse = chord_sum(4096);
    const long double fine = chord_sum(8192);
    total += (4.0L * fine - coarse) / 3.0L;
  }
  return static_cast<double>(total);
}

}  // namespace

TEST(ArcLength, CircleArcAnalytic) {
  const og::Arc arc({1, 2, 0}, 3.0, {0, 0, 1}, {1, 0, 0}, 0.3, 2.5);
  EXPECT_NEAR(arc.length(), 3.0 * 2.2, 1e-10 * 3.0 * 2.2);
  const og::Arc full({0, 0, 0}, 20.0, {0, 0, 1}, {1, 0, 0}, 0.0, 2.0 * std::numbers::pi);
  const double exact = 40.0 * std::numbers::pi;
  EXPECT_LT(std::fabs(full.length() - exact) / exact, 1e-12);
}

TEST(ArcLength, NurbsCircleAnalytic) {
  const auto circle = og::test::nurbs_circle({0, 0, 0}, 20.0);
  const double exact = 40.0 * std::numbers::pi;
  EXPECT_LT(std::fabs(circle->length() - exact) / exact, 1e-12);
}

TEST(ArcLength, MatchesIndependentReference) {
  for (const auto& [name, curve] : og::test::all_curves()) {
    const double reference = reference_length(*curve);
    EXPECT_LT(std::fabs(curve->length() - reference) / reference, 1e-10) << name;
  }
}

TEST(ArcLength, PartialLengthsAreAdditive) {
  for (const auto& [name, curve] : og::test::all_curves()) {
    const double a = curve->length_between(0.0, 0.37);
    const double b = curve->length_between(0.37, 1.0);
    EXPECT_LT(std::fabs(a + b - curve->length()) / curve->length(), 1e-13) << name;
    EXPECT_DOUBLE_EQ(curve->length_between(0.6, 0.2), -curve->length_between(0.2, 0.6));
  }
}

TEST(ArcLength, Reparametrization) {
  constexpr int kSamples = 200;
  for (const auto& [name, curve] : og::test::all_curves()) {
    const double total = curve->length();
    double previous = -1.0;
    for (int k = 0; k <= kSamples; ++k) {
      const double s = total * k / kSamples;
      const double t = curve->t_at_length(s);
      EXPECT_LT(std::fabs(curve->length_between(0.0, t) - s) / total, 1e-10) << name << " k=" << k;
      EXPECT_GT(t, previous) << name << " k=" << k;
      previous = t;
    }
    EXPECT_EQ(curve->t_at_length(-1.0), 0.0);
    EXPECT_EQ(curve->t_at_length(2.0 * total), 1.0);
  }
}

TEST(ArcLength, ArcReparametrizationIsLinear) {
  const og::Arc arc({0, 0, 0}, 1.5, {0, 0, 1}, {1, 0, 0}, 0.0, 4.0);
  for (int k = 0; k <= 50; ++k) {
    const double frac = k / 50.0;
    EXPECT_NEAR(arc.t_at_length(frac * arc.length()), frac, 1e-12);
  }
}
