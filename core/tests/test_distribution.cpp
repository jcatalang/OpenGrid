#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "opengrid/core/dist/distribution.hpp"

using og::DistributionSpec;
using og::DistributionType;

namespace {

DistributionSpec spec(DistributionType type, std::optional<double> begin, std::optional<double> end,
                      std::optional<double> ratio = std::nullopt) {
  DistributionSpec s;
  s.type = type;
  s.begin_spacing = begin;
  s.end_spacing = end;
  s.growth_ratio = ratio;
  return s;
}

void expect_valid(const std::vector<double>& s, int n_points, const std::string& label) {
  ASSERT_EQ(s.size(), static_cast<std::size_t>(n_points)) << label;
  EXPECT_EQ(s.front(), 0.0) << label;
  EXPECT_EQ(s.back(), 1.0) << label;
  for (std::size_t i = 1; i < s.size(); ++i) {
    ASSERT_TRUE(std::isfinite(s[i])) << label << " i=" << i;
    ASSERT_GT(s[i], s[i - 1]) << label << " i=" << i;
  }
}

double first(const std::vector<double>& s) {
  return s[1] - s[0];
}

double last(const std::vector<double>& s) {
  return s[s.size() - 1] - s[s.size() - 2];
}

}  // namespace

TEST(Distribution, Uniform) {
  const auto d = og::distribute(11, {});
  expect_valid(d.s, 11, "uniform");
  EXPECT_NEAR(d.s[3], 0.3, 1e-16);
  EXPECT_TRUE(d.warnings.empty());
  EXPECT_NEAR(og::max_stretching_ratio(d.s), 1.0, 1e-14);
}

TEST(Distribution, InvalidInputThrows) {
  EXPECT_THROW(og::distribute(1, {}), std::invalid_argument);
  EXPECT_THROW(og::distribute(5, spec(DistributionType::Uniform, 0.1, std::nullopt)),
               std::invalid_argument);
  EXPECT_THROW(og::distribute(5, spec(DistributionType::Tanh, 0.0, std::nullopt)),
               std::invalid_argument);
  EXPECT_THROW(og::distribute(5, spec(DistributionType::Tanh, -1.0, std::nullopt)),
               std::invalid_argument);
  EXPECT_THROW(og::distribute(5, spec(DistributionType::Tanh, NAN, std::nullopt)),
               std::invalid_argument);
  EXPECT_THROW(og::distribute(5, spec(DistributionType::Tanh, 0.1, std::nullopt, 1.2)),
               std::invalid_argument);
  EXPECT_THROW(og::distribute(5, spec(DistributionType::Geometric, 0.1, 0.1)),
               std::invalid_argument);
  EXPECT_THROW(og::distribute(5, spec(DistributionType::Geometric, std::nullopt, std::nullopt)),
               std::invalid_argument);
}

TEST(Distribution, GeometricRatioAndSpacings) {
  for (const int n_points : {3, 11, 101}) {
    for (const double ds : {1e-9, 1e-4, 0.01, 0.2}) {
      const auto begin =
          og::distribute(n_points, spec(DistributionType::Geometric, ds, std::nullopt));
      if (ds * (n_points - 1) >= 1.0) {
        continue;
      }
      expect_valid(begin.s, n_points, "geometric begin");
      EXPECT_TRUE(begin.warnings.empty());
      EXPECT_LT(std::fabs(first(begin.s) - ds) / ds, 1e-12) << n_points << " " << ds;
      // Constant ratio.
      const double r = (begin.s[2] - begin.s[1]) / first(begin.s);
      for (std::size_t i = 2; i + 1 < begin.s.size(); ++i) {
        EXPECT_NEAR((begin.s[i + 1] - begin.s[i]) / (begin.s[i] - begin.s[i - 1]), r, 1e-9 * r);
      }
      if (ds < 1e-7) {
        continue;  // far-end spacings below ~1e-7 are not representable to 1e-8
      }
      const auto end =
          og::distribute(n_points, spec(DistributionType::Geometric, std::nullopt, ds));
      expect_valid(end.s, n_points, "geometric end");
      EXPECT_LT(std::fabs(last(end.s) - ds) / ds, 1e-8) << n_points << " " << ds;
    }
  }
  const auto ratio =
      og::distribute(5, spec(DistributionType::Geometric, std::nullopt, std::nullopt, 2.0));
  EXPECT_NEAR(first(ratio.s), 1.0 / 15.0, 1e-15);
  EXPECT_NEAR(last(ratio.s), 8.0 / 15.0, 1e-15);
}

// Begin spacings down to 1e-12 keep full relative accuracy; end spacings (differences of
// numbers near 1) are tested down to 1e-7 where 1e-8 relative accuracy is representable.
TEST(Distribution, VinokurOneSidedAchievesSpacing) {
  for (const auto type : {DistributionType::Tanh, DistributionType::Sine}) {
    for (const int n_points : {3, 5, 21, 101, 1001}) {
      const int n = n_points - 1;
      for (const double factor : {1e-12, 1e-8, 1e-4, 0.1, 0.5, 0.99, 1.0, 1.01, 1.5, 3.0}) {
        const double ds = factor / n;
        if (ds >= 1.0) {
          continue;
        }
        const std::string label = std::string(og::to_string(type)) +
                                  " n=" + std::to_string(n_points) + " ds=" + std::to_string(ds);
        const auto begin = og::distribute(n_points, spec(type, ds, std::nullopt));
        expect_valid(begin.s, n_points, label);
        if (begin.warnings.empty() || type == DistributionType::Sine) {
          EXPECT_LT(std::fabs(first(begin.s) - ds) / ds, 1e-8) << label;
        }
        if (ds >= 1e-7) {
          const auto end = og::distribute(n_points, spec(type, std::nullopt, ds));
          expect_valid(end.s, n_points, label + " end");
          EXPECT_LT(std::fabs(last(end.s) - ds) / ds, 1e-8) << label << " end";
        }
      }
    }
  }
}

TEST(Distribution, VinokurTwoSidedAchievesSpacings) {
  for (const auto type : {DistributionType::Tanh, DistributionType::Sine}) {
    for (const int n_points : {4, 11, 51, 201, 1001}) {
      const int n = n_points - 1;
      for (const double f0 : {1e-6, 1e-3, 0.1, 0.5, 1.0, 1.4}) {
        for (const double f1 : {1e-4, 0.05, 1.0, 1.3}) {
          if (f1 / n < 1e-7) {
            continue;  // far-end spacing not representable to 1e-8 relative
          }
          const double d0 = f0 / n;
          const double d1 = f1 / n;
          if (d0 + d1 >= 0.9 || (n < 4 && (f0 > 1.0 || f1 > 1.0))) {
            continue;
          }
          const std::string label = std::string(og::to_string(type)) +
                                    " n=" + std::to_string(n_points) + " d0=" + std::to_string(d0) +
                                    " d1=" + std::to_string(d1);
          const auto d = og::distribute(n_points, spec(type, d0, d1));
          expect_valid(d.s, n_points, label);
          if (!d.warnings.empty()) {
            // Only the documented sine -> tanh fallback is acceptable here.
            EXPECT_EQ(type, DistributionType::Sine) << label << ": " << d.warnings.front();
            EXPECT_NE(d.warnings.front().find("used tanh"), std::string::npos) << label;
          }
          EXPECT_LT(std::fabs(first(d.s) - d0) / d0, 1e-8) << label;
          EXPECT_LT(std::fabs(last(d.s) - d1) / d1, 1e-8) << label;
        }
      }
    }
  }
}

TEST(Distribution, ExtremeRequestsFallBackWithWarning) {
  struct Case {
    int n_points;
    DistributionSpec spec;
  };
  const std::vector<Case> cases = {
      {11, spec(DistributionType::Tanh, 2.0, std::nullopt)},       // spacing > length
      {11, spec(DistributionType::Tanh, 0.6, 0.6)},                // sum > 1
      {3, spec(DistributionType::Tanh, 0.1, 0.1)},                 // 2 intervals
      {2, spec(DistributionType::Sine, 0.3, std::nullopt)},        // 1 interval
      {11, spec(DistributionType::Geometric, 5.0, std::nullopt)},  // spacing > length
      {101, spec(DistributionType::Geometric, std::nullopt, std::nullopt, 1e10)},   // ratio
      {101, spec(DistributionType::Geometric, std::nullopt, std::nullopt, 1e-10)},  // underflow
      {6, spec(DistributionType::Sine, 0.9, std::nullopt)},  // too large for sine family
  };
  for (std::size_t k = 0; k < cases.size(); ++k) {
    const auto d = og::distribute(cases[k].n_points, cases[k].spec);
    expect_valid(d.s, cases[k].n_points, "case " + std::to_string(k));
    EXPECT_FALSE(d.warnings.empty()) << "case " << k;
  }
}

TEST(Distribution, GeometricTinyButRepresentableSpacing) {
  // d = 1e-300 with 100 intervals needs r ~ 1070: representable, so no warning.
  const auto d = og::distribute(101, spec(DistributionType::Geometric, 1e-300, std::nullopt));
  expect_valid(d.s, 101, "geometric 1e-300");
  EXPECT_TRUE(d.warnings.empty());
  EXPECT_LT(std::fabs(first(d.s) - 1e-300) / 1e-300, 1e-10);
}

TEST(Distribution, TinySpacingsStayFinite) {
  for (const auto type : {DistributionType::Tanh, DistributionType::Sine}) {
    for (const double ds : {1e-15, 1e-30, 1e-100, 1e-300}) {
      const auto one = og::distribute(41, spec(type, ds, std::nullopt));
      expect_valid(one.s, 41, "one-sided tiny");
      const auto two = og::distribute(41, spec(type, ds, 1e-3));
      expect_valid(two.s, 41, "two-sided tiny");
    }
  }
}

TEST(Distribution, RandomRequestsNeverProduceNaN) {
  std::mt19937_64 rng(12345);  // NOLINT(bugprone-random-generator-seed): reproducible test
  std::uniform_real_distribution<double> log_spacing(-14.0, 0.5);
  std::uniform_int_distribution<int> points(2, 400);
  std::uniform_int_distribution<int> kind(0, 3);
  for (int trial = 0; trial < 3000; ++trial) {
    const int n_points = points(rng);
    const double a = std::pow(10.0, log_spacing(rng));
    const double b = std::pow(10.0, log_spacing(rng));
    DistributionSpec s;
    switch (kind(rng)) {
      case 0:
        s = spec(DistributionType::Tanh, a, b);
        break;
      case 1:
        s = spec(DistributionType::Sine, a, b);
        break;
      case 2:
        s = spec(DistributionType::Tanh, a, std::nullopt);
        break;
      default:
        s = spec(DistributionType::Geometric, std::nullopt, std::nullopt, 0.5 + 2.0 * a);
        break;
    }
    const auto d = og::distribute(n_points, s);
    expect_valid(d.s, n_points, "trial " + std::to_string(trial));
  }
}

TEST(Distribution, Deterministic) {
  const auto a = og::distribute(301, spec(DistributionType::Tanh, 1e-5, 2e-3));
  const auto b = og::distribute(301, spec(DistributionType::Tanh, 1e-5, 2e-3));
  EXPECT_EQ(a.s, b.s);
}

TEST(Distribution, MaxStretchingRatio) {
  EXPECT_DOUBLE_EQ(og::max_stretching_ratio({0.0, 0.1, 0.3, 1.0}), 3.5);
  EXPECT_DOUBLE_EQ(og::max_stretching_ratio({0.0, 1.0}), 1.0);
  const auto g =
      og::distribute(21, spec(DistributionType::Geometric, std::nullopt, std::nullopt, 1.2));
  EXPECT_NEAR(og::max_stretching_ratio(g.s), 1.2, 1e-12);
}

TEST(Distribution, ExplicitPositions) {
  DistributionSpec s;
  s.type = DistributionType::Explicit;
  s.positions = {0.0, 0.1, 0.35, 1.0};
  EXPECT_EQ(og::distribute(4, s).s, s.positions);
  EXPECT_THROW((void)og::distribute(5, s), std::invalid_argument);  // wrong count
  s.positions = {0.0, 0.4, 0.3, 1.0};
  EXPECT_THROW((void)og::distribute(4, s), std::invalid_argument);  // not increasing
  s.positions = {0.0, 0.4, 0.9};
  EXPECT_THROW((void)og::distribute(3, s), std::invalid_argument);  // does not end at 1
  s.positions = {0.0, 0.5, 1.0};
  s.begin_spacing = 0.1;
  EXPECT_THROW((void)og::distribute(3, s), std::invalid_argument);
}
