#include "opengrid/core/topo/connector_nodes.hpp"

#include <stdexcept>

namespace og {
namespace {

std::optional<double> normalised(const std::optional<double>& spacing, double length) {
  if (spacing) {
    return *spacing / length;
  }
  return std::nullopt;
}

}  // namespace

ConnectorNodes connector_nodes(const CurveChain& chain, SpacingMode mode,
                               const std::vector<SubConnectorSpec>& subs) {
  if (subs.empty()) {
    throw std::invalid_argument("connector: at least one sub-connector is required");
  }
  for (std::size_t k = 0; k < subs.size(); ++k) {
    double expected_start = 0.0;
    if (k > 0) {
      expected_start = subs[k - 1].end;
    }
    if (subs[k].start != expected_start || !(subs[k].end > subs[k].start)) {
      throw std::invalid_argument("connector: sub-connectors must tile [0, 1] in order (sub " +
                                  std::to_string(k) + ")");
    }
    if (subs[k].n_points < 2) {
      throw std::invalid_argument("connector: sub-connector " + std::to_string(k) +
                                  " needs at least 2 points");
    }
  }
  if (subs.back().end != 1.0) {
    throw std::invalid_argument("connector: the last sub-connector must end at 1");
  }

  const double total = chain.length();
  ConnectorNodes out;
  for (std::size_t k = 0; k < subs.size(); ++k) {
    const SubConnectorSpec& sub = subs[k];
    const double span = sub.end - sub.start;
    DistributionSpec spec = sub.distribution;
    spec.begin_spacing = normalised(spec.begin_spacing, span * total);
    spec.end_spacing = normalised(spec.end_spacing, span * total);
    const Distribution dist = distribute(sub.n_points, spec);
    for (const std::string& warning : dist.warnings) {
      if (subs.size() > 1) {
        out.warnings.push_back("sub-connector " + std::to_string(k + 1) + ": " + warning);
      } else {
        out.warnings.push_back(warning);
      }
    }
    std::size_t first = 0;
    if (k > 0) {
      first = 1;  // shared breakpoint node
    }
    for (std::size_t i = first; i < dist.s.size(); ++i) {
      double fraction = sub.start + span * dist.s[i];
      if (i + 1 == dist.s.size()) {
        fraction = sub.end;
      }
      if (mode == SpacingMode::ArcLength) {
        out.points.push_back(chain.point_at_length(fraction * total));
        out.arc_positions.push_back(fraction * total);
      } else {
        out.points.push_back(chain.point_at_param(fraction));
        out.arc_positions.push_back(chain.length_at_param(fraction));
      }
    }
  }
  out.points.front() = chain.start();
  out.points.back() = chain.end();
  out.max_stretching_ratio = max_stretching_ratio(out.arc_positions);
  return out;
}

}  // namespace og
