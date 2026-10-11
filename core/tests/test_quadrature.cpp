#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "opengrid/core/util/gauss_kronrod.hpp"

TEST(GaussKronrod, Polynomial) {
  // K15 is exact for polynomials up to degree 22.
  const double value =
      og::integrate_gk15([](double x) { return std::pow(x, 9) - 3.0 * x; }, -1.0, 2.0);
  EXPECT_NEAR(value, (std::pow(2.0, 10) - 1.0) / 10.0 - 1.5 * (4.0 - 1.0), 1e-12);
}

TEST(GaussKronrod, Sine) {
  EXPECT_NEAR(og::integrate_gk15([](double x) { return std::sin(x); }, 0.0, std::numbers::pi), 2.0,
              1e-14);
}

TEST(GaussKronrod, ReversedBoundsNegate) {
  const double forward = og::integrate_gk15([](double x) { return std::exp(x); }, 0.0, 1.0);
  const double backward = og::integrate_gk15([](double x) { return std::exp(x); }, 1.0, 0.0);
  EXPECT_DOUBLE_EQ(forward, -backward);
  EXPECT_NEAR(forward, std::numbers::e - 1.0, 1e-14);
}

TEST(GaussKronrod, EndpointSingularDerivative) {
  // sqrt(x) has an infinite derivative at 0; adaptivity must still reach 1e-12.
  EXPECT_NEAR(og::integrate_gk15([](double x) { return std::sqrt(x); }, 0.0, 1.0), 2.0 / 3.0,
              1e-12);
}
