// One-dimensional point distributions on the unit interval: s_0 = 0 < s_1 < ... < s_n = 1.
//
// Types:
//  * Uniform.
//  * Geometric: constant growth ratio r between consecutive intervals; given exactly one of the
//    begin spacing, the end spacing or r (file geometric.cpp).
//  * Tanh: Vinokur's tanh / tan stretching, one-sided or two-sided (file vinokur.cpp).
//  * Sine: the sinh / sin family of Vinokur-type stretching functions (file vinokur.cpp).
//  * Explicit: given positions (e.g. the nodes of a connector that was split or joined, so the
//    new connectors keep the old nodes).
// For Tanh and Sine the *discrete* end spacings s_1 - s_0 and s_n - s_{n-1} match the requested
// values (an outer Newton/Brent iteration corrects Vinokur's derivative-based conditions).
// Requests that cannot be met (spacing >= 1, end spacings summing to >= 1, too few points,
// ratios beyond floating-point range, ...) fall back to a valid distribution and add a
// warning; the result is always finite and strictly increasing.
// Spacings are fractions of the total length.
#pragma once

#include <optional>
#include <string>
#include <vector>

namespace og {

enum class DistributionType { Uniform, Geometric, Tanh, Sine, Explicit };

struct DistributionSpec {
  DistributionType type = DistributionType::Uniform;
  std::optional<double> begin_spacing;  // s_1 - s_0
  std::optional<double> end_spacing;    // s_n - s_{n-1}
  std::optional<double> growth_ratio;   // geometric only: Δ_{i+1} / Δ_i
  std::vector<double> positions;        // explicit only: n_points values, 0 = s_0 < ... = 1
};

struct Distribution {
  std::vector<double> s;              // n_points values in [0, 1]
  std::vector<std::string> warnings;  // empty if the request was met
};

// Throws std::invalid_argument for invalid input: n_points < 2, a non-positive or non-finite
// spacing or ratio, spacings given for a Uniform distribution, a growth ratio for Tanh/Sine,
// or a Geometric spec without exactly one of begin spacing, end spacing and growth ratio, or
// Explicit positions that are not n_points strictly increasing values from 0 to 1.
Distribution distribute(int n_points, const DistributionSpec& spec);

// Largest ratio between consecutive intervals of increasing positions (>= 1; 1 if fewer than
// three points).
double max_stretching_ratio(const std::vector<double>& positions);

const char* to_string(DistributionType type);

}  // namespace og
