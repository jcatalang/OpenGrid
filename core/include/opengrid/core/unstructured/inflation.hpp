// 2D boundary-layer inflation by advancing layers from wall segments.
//
// Method (advancing-layer / "advancing normals" family: S. Pirzadeh, AIAA J. 34 (1996) 43-49;
// Y. Kallinderis, S. Ward, AIAA J. 31 (1993) 1440-1447; with local termination of each column):
//  * Wall chains are maximal runs of wall segments along each boundary loop (closed when the
//    whole loop is wall). The domain lies on the left of every segment, so the left normal points
//    into the domain.
//  * Layer k has height h_k = first_height * growth_rate^k. Every front node advances along its
//    marching normal: the bisector of the left normals of its two front segments, smoothed along
//    the front (Jacobi sweeps that keep the normal visible from both segments), scaled by
//    1 / min(n . n_segment) so the layer keeps the height h_k from both adjacent segments. On a
//    straight wall the first layer is exactly first_height high.
//  * Sharp convex corners (turning angle of the normals above fan_threshold_deg, e.g. sharp or
//    blunt trailing edges) get a fan: several columns from the same wall node with normals
//    interpolated in steps of at most fan_step_deg; the first layer there is made of triangles.
//  * Concave corners: adjacent columns whose new front segment shrinks below half of the old one,
//    or turns back (one would pass the other), are merged into one node (the layer cell becomes
//    a triangle); a run of such pairs (both neighbours of a corner column, once the layers are
//    taller than the wall spacing) becomes one node at the run's mean, so the layers of both
//    walls meet there instead of one side stopping.
//  * A column's isotropic size is its width along the front (at least the wall spacing it
//    started from, so fan columns are not stopped at once), capped by the core size field
//    (size_field.hpp: input boundary spacing, growth rate, max size, sources). A column stops
//    for good when: max_layers (or total_thickness) is reached; with transition, when h_k reaches
//    transition_ratio times that size (cells become isotropic, so the last layer matches the
//    size of the core cells next to it);
//    when it would stand more than max_slope times their distance (measured along the front)
//    above a stopped neighbour, and more than one layer (so fronts taper at most at that slope
//    or by single-layer steps instead of forming cliffs); when its new column edge or
//    front segment would cross or touch the old fronts, the boundary or another new segment
//    (exact predicates); when another front or boundary segment that faces it lies ahead closer
//    than proximity times the core size and approaches it within 60 degrees of head-on (two
//    walls of a narrow gap or a sharp wedge; not the other wall's layers beside it in a concave
//    corner, which meet it there); or when one of its cells would be inverted or have an
//    angle above max_angle_deg. Checks are repeated until no column has to stop, so every
//    committed layer is valid.
//  * Open wall chains (wall meeting a non-wall boundary) keep their end nodes fixed (the
//    neighbouring connectors' nodes are shared with other entities); layers taper to zero towards
//    the ends through the slope rule. With slide, each end column instead slides along the
//    neighbouring boundary through its own vertices, one per layer (so no boundary node is added
//    or lost), and the layers run into that boundary: it needs vertices at the layer heights near
//    the wall (the caller distributes them so). A step whose height across the front is not
//    within a factor slide_ratio of the layer's stops that column (the layers then taper from it);
//    it never passes a corner of the boundary (a turn above 30 degrees) or another wall. Only
//    where the wall and the boundary meet at an angle (inside the domain) between
//    slide_min_angle_deg and slide_max_angle_deg: in a sharper wedge the layers would converge
//    into thin distorted cells, and on a wider opening (a wall continuing into a symmetry line)
//    there is nothing to run into, so those ends taper. The vertices passed are corners of the
//    layer cells, no longer of the core loops.
//  * Layer cells are quads; triangles appear at fans, merges and where a column stops next to an
//    advancing one.
// The final fronts and the non-wall boundary form the loops of the core region, which the caller
// meshes with any 2D generator; the core boundary nodes are the front nodes, so the result is
// conforming. With even_loops, adjacent front nodes are merged once on every core loop with an odd
// number of segments (for all-quad cores).
#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "opengrid/core/unstructured/triangulate.hpp"

namespace og {

struct InflationOptions {
  double first_height = 1e-3;
  double growth_rate = 1.2;
  int max_layers = 30;
  std::optional<double> total_thickness;
  bool transition = true;
  double transition_ratio = 1.0;
  double fan_threshold_deg = 60.0;
  double fan_step_deg = 20.0;
  double proximity = 1.0;
  double max_angle_deg = 165.0;
  double max_slope = 0.5;
  int normal_smoothing = 5;
  bool slide = false;        // end columns of open wall chains slide along the boundary
  double slide_ratio = 2.0;  // a slide step's height may differ from the layer's by this factor
  double slide_min_angle_deg = 45.0;   // sharper wall / boundary wedges taper
  double slide_max_angle_deg = 135.0;  // wider openings taper
};

// Wall columns that ended before the last layer, by the rule that stopped them (open chains'
// end columns, which taper or slide, are not counted).
struct InflationStops {
  std::int64_t columns = 0;     // wall columns
  std::int64_t transition = 0;  // a layer as thick as the core cells beside it
  std::int64_t slope = 0;       // too steep next to a stopped column
  std::int64_t proximity = 0;   // near another front or boundary, crossing it, or a bad cell
};

struct InflationResult {
  std::vector<std::array<double, 2>> points;  // input points first, then layer nodes
  std::vector<std::array<std::int64_t, 3>> triangles;
  std::vector<std::array<std::int64_t, 4>> quads;
  std::vector<std::vector<std::int64_t>> core_loops;  // closed, domain on the left
  int layers = 0;                                     // largest number of layers of any column
  InflationStops stops;
  std::vector<std::string> warnings;
};

// wall[s] marks input segment s as a wall. Throws std::invalid_argument for invalid input (as
// triangulate_domain), a wall flag per segment missing, or invalid options.
InflationResult inflate(const TriangulationInput& input, const std::vector<bool>& wall,
                        const InflationOptions& options, const SizeFieldOptions& size,
                        bool even_loops);

}  // namespace og
