// Transfinite interpolation of a structured block from its six faces.
//
// Method: Boolean sum of the three linear face projectors,
//   x = U + V + W - UV - VW - WU + UVW,
// where U interpolates between the imin and imax faces, UV between the four k-direction edges,
// UVW is trilinear in the eight corners (Gordon & Hall, Int. J. Numer. Meth. Eng. 7 (1973)
// 461-477; Thompson, Warsi & Mastin, Numerical Grid Generation (1985), ch. VIII). The blending
// parameters (xi, eta, zeta) of a node are the arc-length fractions of the grid lines through it
// on the faces, blended bilinearly across the block (the 3D analogue of Soni's 2D choice), so
// opposite faces with equal distributions are reproduced exactly (e.g. an annular sector).
// The faces are reproduced exactly for any parameters.
#pragma once

#include "opengrid/core/structured/grid3d.hpp"

namespace og {

// Throws std::invalid_argument for faces with inconsistent sizes or shared edges / corners
// that do not coincide (tolerance relative to the block size).
Grid3D tfi_block(const BlockFaces& faces);

// The six faces of a block grid (as consumed by tfi_block).
BlockFaces block_faces(const Grid3D& grid);

}  // namespace og
