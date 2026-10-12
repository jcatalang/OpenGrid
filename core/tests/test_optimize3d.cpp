#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "opengrid/core/unstructured/optimize3d.hpp"

namespace {

using Index = std::int64_t;
using Points = std::vector<std::array<double, 3>>;

TEST(Optimize3D, IdealCellsScoreOne) {
  const double s = std::numbers::sqrt3 / 2.0;
  const Points tet = {{0, 0, 0}, {1, 0, 0}, {0.5, s, 0}, {0.5, s / 3.0, std::sqrt(2.0 / 3.0)}};
  const std::array<Index, 4> t = {0, 1, 2, 3};
  EXPECT_NEAR(og::cell_quality_3d(tet, t.data(), 4), 1.0, 1e-12);
  const std::array<Index, 4> inverted = {0, 2, 1, 3};
  EXPECT_LT(og::cell_quality_3d(tet, inverted.data(), 4), 0.0);
  const Points pyramid = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, std::sqrt(0.5)}};
  const std::array<Index, 5> p = {0, 1, 2, 3, 4};
  EXPECT_NEAR(og::cell_quality_3d(pyramid, p.data(), 5), 1.0, 1e-12);
  const Points prism = {{0, 0, 0}, {1, 0, 0}, {0.5, s, 0}, {0, 0, 1}, {1, 0, 1}, {0.5, s, 1}};
  const std::array<Index, 6> w = {0, 1, 2, 3, 4, 5};
  EXPECT_NEAR(og::cell_quality_3d(prism, w.data(), 6), 1.0, 1e-12);
  const Points cube = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                       {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  const std::array<Index, 8> h = {0, 1, 2, 3, 4, 5, 6, 7};
  EXPECT_NEAR(og::cell_quality_3d(cube, h.data(), 8), 1.0, 1e-12);
}

// The unit cube as 12 tetrahedra joining its centre node to the two triangles of every face.
og::Optimize3DInput cube_around_centre(const std::array<double, 3>& centre) {
  og::Optimize3DInput in;
  in.points = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1},
               {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, centre};
  // Faces with outward normals (counter-clockwise seen from outside).
  const std::vector<std::array<Index, 4>> faces = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                                   {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  for (const auto& f : faces) {
    // A tetrahedron (a, b, c, centre) is positive when (a, b, c) is clockwise seen from outside.
    in.tets.push_back({f[0], f[2], f[1], 8});
    in.tets.push_back({f[0], f[3], f[2], 8});
  }
  in.movable.assign(in.points.size(), false);
  in.movable[8] = true;
  return in;
}

TEST(Optimize3D, UntanglesAnInteriorNode) {
  const og::Optimize3DInput good = cube_around_centre({0.5, 0.5, 0.5});
  for (const auto& t : good.tets) {
    EXPECT_GT(og::cell_quality_3d(good.points, t.data(), 4), 0.0);
  }
  og::Optimize3DInput in = cube_around_centre({1.4, 0.2, 0.7});
  og::Optimize3DOptions opt;
  opt.swap = false;
  const og::Optimize3DResult r = og::optimize_3d(in, opt);
  EXPECT_LT(r.worst_before, 0.0);
  EXPECT_EQ(r.inverted_after, 0);
  EXPECT_GT(r.worst_after, 0.3);
  EXPECT_EQ(r.moved, 1);
  // Fixed nodes stay; the centre ends near the middle.
  for (std::size_t i = 0; i < 8; ++i) {
    EXPECT_EQ(r.points[i], in.points[i]);
  }
  EXPECT_NEAR(r.points[8][0], 0.5, 0.05);
  EXPECT_NEAR(r.points[8][1], 0.5, 0.05);
  EXPECT_NEAR(r.points[8][2], 0.5, 0.05);
  // Nothing to do for a good mesh.
  const og::Optimize3DResult same = og::optimize_3d(good, opt);
  EXPECT_EQ(same.moved, 0);
  EXPECT_EQ(same.points, good.points);
}

// Two tetrahedra on the face (a, b, c) of an equilateral triangle, apexes d (above) and e (below)
// at height h over its centre.
og::Optimize3DInput two_tets(double h) {
  const double s = std::numbers::sqrt3 / 2.0;
  og::Optimize3DInput in;
  in.points = {{0, 0, 0}, {1, 0, 0}, {0.5, s, 0}, {0.5, s / 3.0, h}, {0.5, s / 3.0, -h}};
  in.tets = {{0, 1, 2, 3}, {1, 0, 2, 4}};
  in.movable.assign(in.points.size(), false);
  return in;
}

TEST(Optimize3D, FaceAndEdgeSwaps) {
  og::Optimize3DOptions opt;
  opt.threshold = 0.5;
  opt.move_nodes = false;
  // Flat pair: three tetrahedra around the short edge d-e are better (2-3).
  const og::Optimize3DResult flat = og::optimize_3d(two_tets(0.2), opt);
  EXPECT_EQ(flat.swaps_23, 1);
  EXPECT_EQ(flat.tets.size(), 3U);
  EXPECT_GT(flat.worst_after, flat.worst_before);
  EXPECT_NEAR(flat.worst_after, 0.437, 1e-3);
  // Regular pair: no swap.
  const og::Optimize3DResult regular = og::optimize_3d(two_tets(std::sqrt(2.0 / 3.0)), opt);
  EXPECT_EQ(regular.swaps_23 + regular.swaps_32, 0);
  EXPECT_EQ(regular.tets.size(), 2U);
  // Three tetrahedra around a long interior edge become the regular pair (3-2).
  og::Optimize3DInput three = two_tets(std::sqrt(2.0 / 3.0));
  three.tets = {{0, 1, 4, 3}, {1, 2, 4, 3}, {2, 0, 4, 3}};
  const og::Optimize3DResult merged = og::optimize_3d(three, opt);
  EXPECT_EQ(merged.swaps_32, 1);
  ASSERT_EQ(merged.tets.size(), 2U);
  EXPECT_NEAR(merged.worst_after, 1.0, 1e-9);
  // Fixed tetrahedra are never swapped.
  three.n_fixed_tets = 1;
  const og::Optimize3DResult fixed = og::optimize_3d(three, opt);
  EXPECT_EQ(fixed.swaps_32, 0);
  EXPECT_EQ(fixed.tets, three.tets);
}

TEST(Optimize3D, OtherCellsOnTheEdgeBlockEdgeRemoval) {
  // The three tetrahedra around d-e, and a pyramid (two extra base nodes) that also contains the
  // edge d-e: the edge is not interior to the tetrahedra, so it is not removed.
  og::Optimize3DInput in = two_tets(std::sqrt(2.0 / 3.0));
  in.tets = {{0, 1, 4, 3}, {1, 2, 4, 3}, {2, 0, 4, 3}};
  in.points.push_back({0.5, -1.0, -1.0});
  in.points.push_back({0.5, -1.0, 1.0});
  in.movable.assign(in.points.size(), false);
  in.pyramids = {{5, 6, 3, 4, 0}};
  og::Optimize3DOptions opt;
  opt.threshold = 0.5;
  opt.move_nodes = false;
  const og::Optimize3DResult r = og::optimize_3d(in, opt);
  EXPECT_EQ(r.swaps_32, 0);
  EXPECT_EQ(r.tets.size(), 3U);
}

TEST(Optimize3D, RescueLiftsCellsToTheFloor) {
  // The centre node next to the bottom face: the bottom tetrahedra are flat (valid but far below
  // the floor). With threshold 0 the main passes only untangle, so the rescue does the work.
  og::Optimize3DInput in = cube_around_centre({0.5, 0.5, 0.01});
  og::Optimize3DOptions opt;
  opt.threshold = 0.0;
  opt.swap = false;
  const og::Optimize3DResult r = og::optimize_3d(in, opt);
  EXPECT_GT(r.below_floor_before, 0);
  EXPECT_EQ(r.below_floor_after, 0);
  EXPECT_GE(r.moved, 1);
  for (const auto& t : r.tets) {
    EXPECT_TRUE(opt.floor.passes(og::cell_quality_3d_full(r.points, t.data(), 4)));
  }
}

TEST(Optimize3D, QualityOfCells) {
  // Regular tetrahedron: skewness 0, dihedral 70.53 degrees; cube: 90 degrees everywhere; a flat
  // tetrahedron misses the floor.
  const double s3 = std::numbers::sqrt3 / 2.0;
  const Points tet = {{0, 0, 0}, {1, 0, 0}, {0.5, s3, 0}, {0.5, s3 / 3.0, std::sqrt(2.0 / 3.0)}};
  const std::array<Index, 4> t = {0, 1, 2, 3};
  const og::CellQuality3D q = og::cell_quality_3d_full(tet, t.data(), 4);
  EXPECT_NEAR(q.skewness, 0.0, 1e-12);
  EXPECT_NEAR(q.min_angle, 70.528779365509, 1e-9);
  EXPECT_NEAR(q.max_angle, 70.528779365509, 1e-9);
  const Points cube = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                       {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  const std::array<Index, 8> h = {0, 1, 2, 3, 4, 5, 6, 7};
  const og::CellQuality3D c = og::cell_quality_3d_full(cube, h.data(), 8);
  EXPECT_NEAR(c.skewness, 0.0, 1e-12);
  EXPECT_NEAR(c.min_angle, 90.0, 1e-12);
  const Points flat = {{0, 0, 0}, {1, 0, 0}, {0.5, s3, 0}, {0.5, s3 / 3.0, 0.01}};
  const og::CellQuality3D f = og::cell_quality_3d_full(flat, t.data(), 4);
  EXPECT_FALSE(og::QualityFloor{}.passes(f));
  EXPECT_GT(f.skewness, 0.95);
}

TEST(Optimize3D, InvalidInput) {
  og::Optimize3DInput in = two_tets(0.5);
  in.movable.pop_back();
  EXPECT_THROW((void)og::optimize_3d(in, {}), std::invalid_argument);
  in = two_tets(0.5);
  in.tets.push_back({0, 1, 2, 9});
  EXPECT_THROW((void)og::optimize_3d(in, {}), std::invalid_argument);
  const std::array<Index, 3> three = {0, 1, 2};
  EXPECT_THROW((void)og::cell_quality_3d(in.points, three.data(), 3), std::invalid_argument);
}

}  // namespace
