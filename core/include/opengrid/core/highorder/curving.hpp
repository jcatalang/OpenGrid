// Extension of boundary-node displacements into a high-order mesh (Phase 11g: curved boundary
// elements).
//
// Method: the high-order nodes moved onto the geometry (curved walls) carry a displacement
// that the rest of the mesh must follow, or the elements next to the wall fold (a sagitta of
// 1e-3 against a first cell of 1e-5). The displacement is extended as the solution of a
// weighted graph Laplacian on the nodes: edges join neighbouring lattice nodes of each element,
// with weights 1 / length^2 (the scaling of the discrete Laplacian), so nodes across thin
// boundary-layer cells are strongly coupled and a layer stack bends with the wall instead of
// being squashed, and the displacement decays over the wall elements' length scale. Fixed
// nodes are Dirichlet values (the moved boundary nodes; the rest of the boundary at zero).
// Free nodes take part only within a region around the moved nodes: a multi-source Dijkstra
// along the edges with each source's distance divided by its radius (normalised distance <=
// 1); the others keep zero displacement. The free nodes' system is symmetric positive
// definite: sparse Cholesky (Eigen SimplicialLDLT, three right-hand sides) up to 20k unknowns,
// conjugate gradients with an incomplete Cholesky preconditioner above (fill-in of 3D
// regions). Weights 1 / length^3 by default: stiffer than the Laplacian's 1 / length^2 across
// thin cells, which halved the loss of scaled Jacobian on a coarse curved wall over thin
// layers (0.20 -> 0.10; 1 / length^4: 0.08 there but worse on a C-grid).
// References: P.-O. Persson, J. Peraire, "Curved mesh generation and mesh refinement using
// Lagrangian solid mechanics", AIAA 2009-949 (the elastic analogy; here its Laplace form);
// Z. Q. Xie, R. Sevilla, O. Hassan, K. Morgan, "The generation of arbitrary order curved meshes
// for 3D finite element analysis", Comput. Mech. 51 (2013) 361-374.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"

namespace og {

struct DisplacementExtensionInput {
  std::vector<Vec3> points;                        // node positions (before the displacement)
  std::vector<std::array<std::int64_t, 2>> edges;  // neighbouring nodes (duplicates allowed)
  std::vector<std::uint8_t> fixed;                 // 1: the displacement is given
  std::vector<Vec3> displacement;                  // given values (fixed nodes; others ignored)
  std::vector<double> radius;  // > 0: a moved node's reach (region seeds); 0 for the others
  double exponent = 3.0;       // edge weights 1 / length^exponent (see the header comment)
};

struct DisplacementExtension {
  std::vector<Vec3> displacement;  // every node: given, solved (region) or zero
  std::int64_t solved = 0;         // free nodes in the region
  std::int64_t seeds = 0;          // moved nodes
  std::int64_t iterations = 0;     // most conjugate-gradient iterations of x, y, z (0: direct)
};

// Throws std::invalid_argument for inconsistent sizes, out-of-range or degenerate edges, and
// std::runtime_error if the factorisation fails.
DisplacementExtension extend_displacement(const DisplacementExtensionInput& input);

}  // namespace og
