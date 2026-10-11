// Node generation for a connector: a curve chain divided into sub-connectors, each with its
// own number of points and 1D distribution.
//
// Sub-connector k covers the chain fraction [start_k, end_k] (arc-length fractions in
// ArcLength mode, chain-parameter fractions in Parametric mode); consecutive sub-connectors
// share their breakpoint node. Spacings in the specs are absolute lengths; they are divided
// by the sub-connector length (ArcLength) or by (end - start) * L (Parametric, where the
// achieved physical spacing therefore differs from the request on non-uniformly
// parametrised curves).
#pragma once

#include <string>
#include <vector>

#include "opengrid/core/dist/distribution.hpp"
#include "opengrid/core/geom/curve_chain.hpp"

namespace og {

enum class SpacingMode { ArcLength, Parametric };

struct SubConnectorSpec {
  double start = 0.0;
  double end = 1.0;
  int n_points = 2;
  DistributionSpec distribution;  // spacings in absolute length units
};

struct ConnectorNodes {
  std::vector<Vec3> points;
  std::vector<double> arc_positions;  // arc length of every node from the chain start
  std::vector<std::string> warnings;
  double max_stretching_ratio = 1.0;  // over consecutive arc-length intervals
};

// Throws std::invalid_argument if the sub-connectors do not tile [0, 1] in increasing order
// or a sub-connector has fewer than 2 points.
ConnectorNodes connector_nodes(const CurveChain& chain, SpacingMode mode,
                               const std::vector<SubConnectorSpec>& subs);

}  // namespace og
