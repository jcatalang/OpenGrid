// Advancing-front triangulation of a 2D domain bounded by closed loops whose vertices are
// preserved exactly (same input and output conventions as triangulate_domain).
//
// Method (S. H. Lo, Int. J. Numer. Meth. Eng. 21 (1985) 1403-1426; J. Peraire, M. Vahdati,
// K. Morgan, O. C. Zienkiewicz, J. Comput. Phys. 72 (1987) 449-466):
//  * the front is the set of oriented edges with the unmeshed region on their left, initially the
//    boundary segments; the shortest front edge is processed first (ties by insertion order);
//  * for an edge (a, b) the ideal apex lies on the inward normal at the height of an equilateral
//    triangle whose side is the local size h, clamped to [0.55, 1.3] times the edge length so the
//    size changes gradually from layer to layer;
//  * candidates, in order: existing front nodes within 0.9 h of the ideal apex (by distance),
//    the ideal point, the other nearby front nodes; the first acceptable candidate is taken;
//  * a candidate is acceptable if the triangle is positively oriented, its new edges neither
//    cross nor touch the front (except the front edges it closes), no front node lies inside it,
//    and a new point keeps 0.5 h from front nodes and 0.3 h from front edges; all geometric
//    tests use the exact predicates;
//  * an edge without an acceptable candidate is retried later with a wider search; edges still
//    open at the end are closed by the constrained Delaunay mesher (triangulate_domain) on the
//    remaining polygons, with a warning;
//  * finally: Lawson flips of unconstrained edges, three passes of guarded Laplacian smoothing
//    of the interior nodes (a move is undone if it lowers the smallest angle around the node),
//    and flips again.
// Front nodes and edges are indexed in a uniform bucket grid. Output is deterministic.
#pragma once

#include "opengrid/core/unstructured/triangulate.hpp"

namespace og {

TriangleMesh advancing_front(const TriangulationInput& input, const RefinementOptions& options);

// Lawson flips of all edges not listed as boundary segments until every interior edge is locally
// Delaunay (only convex quadrilaterals are flipped). Returns the number of flips.
std::int64_t flip_to_delaunay(const std::vector<std::array<double, 2>>& points,
                              std::vector<std::array<MeshIndex, 3>>& triangles,
                              const std::vector<std::array<MeshIndex, 2>>& segments);

// Guarded Laplacian smoothing: each node with index >= n_fixed moves to the centroid of its
// neighbours unless that lowers the smallest angle of its incident triangles (or inverts one).
// Returns the number of accepted moves.
std::int64_t smooth_interior(std::vector<std::array<double, 2>>& points,
                             const std::vector<std::array<MeshIndex, 3>>& triangles,
                             std::size_t n_fixed, int passes);

}  // namespace og
