// Robust geometric predicates for 2D points.
//
// Method: J. R. Shewchuk, "Adaptive Precision Floating-Point Arithmetic and Fast Robust
// Geometric Predicates", Discrete Comput. Geom. 18 (1997) 305-363. Each predicate is first
// evaluated in double precision with a forward error bound (Shewchuk's stage-A bounds); only if
// the result is too close to zero is it recomputed exactly with expansion arithmetic
// (two_sum / two_product / scale_expansion / expansion_sum). Implemented here from the paper
// (the original code is public domain). Requires IEEE double arithmetic without fused
// multiply-add contraction (the build uses -ffp-contract=off).
#pragma once

namespace og {

// +1 if (a, b, c) are in counter-clockwise order, -1 if clockwise, 0 if collinear (exact).
int orient2d(const double* a, const double* b, const double* c);

// +1 if d lies inside the circle through a, b, c (counter-clockwise), -1 outside, 0 on it.
int incircle(const double* a, const double* b, const double* c, const double* d);

}  // namespace og
