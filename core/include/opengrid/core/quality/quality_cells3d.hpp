// Per-cell quality of unstructured 3D cells (tetrahedra, pyramids, prisms, hexahedra), the same
// definitions as opengrid.model.volumes.volume_quality, and the quality floor every 3D cell must
// meet.
//  * scaled Jacobian: the smallest corner triple product of the three corner edges divided by
//    their lengths, normalised so the ideal cell scores 1 (regular tetrahedron, pyramid with
//    equilateral triangles, right prism on an equilateral triangle with square sides, cube) and
//    capped at 1; <= 0 means an inverted corner. Corners as CELL_CORNERS in volumes.py (a
//    pyramid's apex is not a corner; it enters through the base corners).
//  * dihedral angles: between the two faces at every edge (face normals: outward area vectors,
//    bilinear for quads); min_angle / max_angle in degrees.
//  * equiangle skewness: max((max_angle - a_hi) / (180 - a_hi), (a_lo - min_angle) / a_lo, 0)
//    with the type's ideal smallest / largest dihedral angle (a_lo, a_hi): tetrahedron 70.53 /
//    70.53, pyramid 54.74 / 109.47, prism 60 / 90, hexahedron 90 / 90. 0 is ideal, 1 degenerate.
// Reference: P. M. Knupp, "Achieving finite element mesh quality via optimization of the
// Jacobian matrix norm and associated quantities", Int. J. Numer. Meth. Engng 48 (2000); equiangle
// skewness is the definition CFD solvers' mesh checks use.
// Cell node order: CGNS / VTK (tetrahedron (0, 1, 2) normal towards 3, pyramid base normal
// towards the apex, prism / hexahedron first face normal towards the second).
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace og {

struct CellQuality3D {
  double scaled_jacobian = 1.0;
  double skewness = 0.0;
  double min_angle = 0.0;  // smallest dihedral angle, degrees
  double max_angle = 0.0;  // largest dihedral angle, degrees
};

// The quality every 3D cell must meet: a scaled Jacobian of at least min_jacobian and an
// equiangle skewness of at most max_skewness.
struct QualityFloor {
  double min_jacobian = 0.05;
  double max_skewness = 0.95;
  [[nodiscard]] bool passes(const CellQuality3D& q) const {
    return q.scaled_jacobian >= min_jacobian && q.skewness <= max_skewness;
  }
};

// Signed scaled Jacobian of a cell with n = 4 (tetrahedron), 5 (pyramid), 6 (prism) or 8
// (hexahedron) nodes. Throws std::invalid_argument for other n.
double cell_quality_3d(const std::vector<std::array<double, 3>>& points, const std::int64_t* nodes,
                       int n);

// All metrics of one cell (n as above).
CellQuality3D cell_quality_3d_full(const std::vector<std::array<double, 3>>& points,
                                   const std::int64_t* nodes, int n);

}  // namespace og
