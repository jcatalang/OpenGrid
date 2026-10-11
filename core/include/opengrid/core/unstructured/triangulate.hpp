// Constrained Delaunay triangulation and Delaunay refinement of a 2D domain bounded by closed
// polygonal loops whose vertices must be preserved exactly (they are connector nodes shared with
// neighbouring domains).
//
// Method:
//  1. Bowyer-Watson incremental insertion of the boundary points inside a large super triangle,
//     with a visibility walk for point location (A. Bowyer, Comput. J. 24 (1981); D. F. Watson,
//     Comput. J. 24 (1981)).
//  2. Segment recovery by edge flips (S. W. Sloan, Comput. Struct. 47 (1993) 441-450), then
//     Lawson flips of unconstrained edges to restore the constrained Delaunay property.
//  3. Classification by flood fill from the oriented boundary segments (domain on the left),
//     never crossing constrained edges; everything else (holes, super triangle) is discarded.
//  4. Delaunay refinement (L. P. Chew 1989; J. Ruppert, J. Algorithms 18 (1995) 548-585):
//     triangles with a radius-edge ratio above 1 / (2 sin(min_angle)) or larger than the size
//     field are split: oversized ones at their circumcenter, badly shaped ones at their
//     off-center (A. Ungor, Comput. Geom. 42 (2009) 109-118). Boundary segments are never split:
//     a new point
//     that lies across a boundary segment or inside the diametral lens of one (the segment
//     subtends more than 120 degrees; J. R. Shewchuk, Comput. Geom. 22 (2002)) is rejected
//     (as Triangle's -Y switch); for oversized triangles without a boundary edge the centroid
//     is tried instead. A triangle's target size is at least the length of its longest boundary
//     edge. The quality guarantee therefore holds away from the boundary; near it, the boundary
//     spacing governs.
// All predicates are exact (util/predicates.hpp); processing order is fixed, so the output is
// deterministic.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "opengrid/core/unstructured/size_field.hpp"

namespace og {

using MeshIndex = std::int64_t;

struct TriangulationInput {
  std::vector<std::array<double, 2>> points;       // boundary nodes, no duplicates
  std::vector<std::array<MeshIndex, 2>> segments;  // oriented: the domain lies on the left
};

struct RefinementOptions {
  bool refine = true;
  double min_angle_deg = 25.0;
  SizeFieldOptions size;
  MeshIndex max_points = 20'000'000;
};

struct TriangleMesh {
  std::vector<std::array<double, 2>> points;  // boundary points first (input order), then Steiner
  std::vector<std::array<MeshIndex, 3>> triangles;  // counter-clockwise
  MeshIndex rejected = 0;  // refinement points rejected for encroaching on the boundary
  std::vector<std::string> warnings;
};

// Throws std::invalid_argument for duplicate points, segments with invalid indices, a boundary
// that is not closed around the domain, or intersecting segments.
TriangleMesh triangulate_domain(const TriangulationInput& input, const RefinementOptions& options);

}  // namespace og
