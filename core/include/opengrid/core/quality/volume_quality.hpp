// Per-cell quality metrics of a mixed 3D mesh (tetrahedra, pyramids, prisms, hexahedra).
//
// Method: for every cell, in the order tets, pyramids, prisms, hexes:
//  * scaled Jacobian: the smallest corner det(e1, e2, e3) / (|e1| |e2| |e3|), scaled to 1 for the
//    ideal cell and capped at 1 (tetrahedra: the smallest |corner| with the sign of the volume);
//  * aspect ratio: circumradius / (3 inradius) for tetrahedra, longest / shortest edge otherwise;
//  * smallest / largest dihedral angle between adjacent faces (degrees), and the equiangle
//    skewness against the type's ideal angles (Bhaskaran's definition);
//  * signed volume by the divergence theorem over the outward faces (bilinear quads);
//  * non-orthogonality: largest angle between a face's area vector and the link of the two cell
//    centroids (boundary faces: the face centre);
//  * volume ratio: largest volume ratio with a face neighbour.
// Faces are matched by a counting sort on their smallest node followed by a comparison of the
// sorted node lists within each bucket (linear in the number of faces, deterministic). Per-cell
// work runs in parallel (util/parallel.hpp); the results do not depend on the thread count.
// Reference: the same definitions as quality_cells3d.hpp and the Python reference
// implementation in tests/python (volume_quality_reference).
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace og {

struct CellArray {
  const std::int64_t* nodes = nullptr;  // count x size, row-major
  std::size_t count = 0;
  int size = 4;  // 4, 5, 6 or 8
};

struct VolumeQuality {
  std::vector<double> scaled_jacobian;
  std::vector<double> aspect_ratio;
  std::vector<double> min_angle;
  std::vector<double> max_angle;
  std::vector<double> skewness;
  std::vector<double> volume;
  std::vector<double> non_orthogonality;
  std::vector<double> volume_ratio;
};

// points: n_points x 3. Cell arrays in the order tets, pyramids, prisms, hexes. Node indices
// must be valid (checked by the caller).
VolumeQuality volume_quality(const double* points, std::size_t n_points,
                             const std::vector<CellArray>& cells);

}  // namespace og
