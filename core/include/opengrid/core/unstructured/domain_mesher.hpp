// Complete 2D unstructured domain meshing: triangulation, optional quad recombination or all-quad
// subdivision, smoothing, edge swapping and untangling. Boundary nodes are always the input
// points (first, in input order), so neighbouring domains sharing connectors stay conforming.
//
// Element types:
//  * triangles      the triangulation (Delaunay refinement or advancing front);
//  * quad_dominant  the triangulation recombined into quads (recombine.hpp);
//  * quads          all quadrilaterals. Every boundary loop must have an even number of segments
//                   (a necessary condition for any quad mesh of the loop's region). The domain is
//                   meshed with a coarse boundary made of every second node of each loop, at twice
//                   the target size, recombined, and every cell is split (triangle -> 3 quads,
//                   quad -> 4) through edge midpoints and the cell centre; the midpoint of a coarse
//                   boundary edge is the skipped boundary node itself, so the boundary is exactly
//                   the input boundary (the "recombine then subdivide" route of Gmsh's full-quad
//                   option; C. Geuzaine, J.-F. Remacle, Int. J. Numer. Meth. Eng. 79 (2009)). New
//                   nodes are then smoothed (three guarded Laplacian sweeps).
// Afterwards the requested smoothing (smoothing.hpp) runs, and any cell that is still inverted or
// non-convex is untangled; if that fails a warning names the number of such cells.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "opengrid/core/unstructured/recombine.hpp"
#include "opengrid/core/unstructured/smoothing.hpp"
#include "opengrid/core/unstructured/triangulate.hpp"

namespace og {

enum class TriangleAlgorithm { delaunay, advancing_front };
enum class ElementType { triangles, quad_dominant, quads };

struct DomainMeshOptions {
  TriangleAlgorithm algorithm = TriangleAlgorithm::delaunay;
  ElementType element_type = ElementType::triangles;
  RefinementOptions refinement;
  RecombineOptions recombine;
  SmoothingOptions smoothing;
};

struct MixedMesh {
  std::vector<std::array<double, 2>> points;           // boundary points first (input order)
  std::vector<std::array<std::int64_t, 3>> triangles;  // counter-clockwise
  std::vector<std::array<std::int64_t, 4>> quads;      // counter-clockwise
  std::int64_t rejected = 0;
  std::vector<std::string> warnings;
};

// Throws std::invalid_argument for invalid boundaries (as triangulate_domain) and, for all-quad
// meshes, for loops with an odd number of segments or boundary points shared by several loops.
MixedMesh mesh_domain(const TriangulationInput& input, const DomainMeshOptions& options);

}  // namespace og
