#include "opengrid/core/highorder/curving.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace og {

namespace {

// Regions up to this many unknowns are solved by sparse Cholesky, larger ones by
// preconditioned conjugate gradients (relative residual kTolerance).
constexpr Eigen::Index kDirectLimit = 20000;
constexpr double kTolerance = 1e-10;
constexpr int kMaxIterations = 5000;

struct Adjacency {
  std::vector<std::int64_t> start;  // n + 1
  std::vector<std::int64_t> node;
  std::vector<double> length;
};

Adjacency build_adjacency(const DisplacementExtensionInput& in) {
  const auto n = static_cast<std::int64_t>(in.points.size());
  std::vector<std::array<std::int64_t, 2>> edges;
  edges.reserve(in.edges.size());
  for (const auto& [a, b] : in.edges) {
    if (a < 0 || b < 0 || a >= n || b >= n || a == b) {
      throw std::invalid_argument("extend_displacement: edge (" + std::to_string(a) + ", " +
                                  std::to_string(b) + ") is out of range or degenerate");
    }
    edges.push_back({std::min(a, b), std::max(a, b)});
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  Adjacency adj;
  adj.start.assign(static_cast<std::size_t>(n) + 1, 0);
  for (const auto& [a, b] : edges) {
    ++adj.start[static_cast<std::size_t>(a) + 1];
    ++adj.start[static_cast<std::size_t>(b) + 1];
  }
  for (std::size_t k = 1; k < adj.start.size(); ++k) {
    adj.start[k] += adj.start[k - 1];
  }
  adj.node.resize(2 * edges.size());
  adj.length.resize(2 * edges.size());
  std::vector<std::int64_t> fill(adj.start.begin(), adj.start.end() - 1);
  for (const auto& [a, b] : edges) {
    const double l =
        norm(in.points[static_cast<std::size_t>(a)] - in.points[static_cast<std::size_t>(b)]);
    if (!(l > 0.0)) {
      throw std::invalid_argument("extend_displacement: nodes " + std::to_string(a) + " and " +
                                  std::to_string(b) + " coincide");
    }
    for (const auto& [from, to] : {std::pair{a, b}, std::pair{b, a}}) {
      const auto slot = static_cast<std::size_t>(fill[static_cast<std::size_t>(from)]++);
      adj.node[slot] = to;
      adj.length[slot] = l;
    }
  }
  return adj;
}

// Normalised distance to the nearest seed: path length / the seed's radius (multi-source
// Dijkstra; ties broken by node index, so the result is deterministic).
std::vector<double> region_distance(const DisplacementExtensionInput& in, const Adjacency& adj,
                                    std::int64_t& seeds) {
  const std::size_t n = in.points.size();
  constexpr double kInf = std::numeric_limits<double>::infinity();
  std::vector<double> best(n, kInf);
  std::vector<double> reach(n, 0.0);  // radius of the seed that reached the node
  using Item = std::tuple<double, std::int64_t>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> queue;
  seeds = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (in.fixed[i] != 0 && in.radius[i] > 0.0 && norm(in.displacement[i]) > 0.0) {
      best[i] = 0.0;
      reach[i] = in.radius[i];
      queue.emplace(0.0, static_cast<std::int64_t>(i));
      ++seeds;
    }
  }
  std::vector<std::uint8_t> done(n, 0);
  while (!queue.empty()) {
    const auto [d, node] = queue.top();
    queue.pop();
    const auto u = static_cast<std::size_t>(node);
    if (done[u] != 0) {
      continue;
    }
    done[u] = 1;
    if (d > 1.0) {
      continue;  // beyond every seed's reach: do not expand
    }
    for (auto k = adj.start[u]; k < adj.start[u + 1]; ++k) {
      const auto v = static_cast<std::size_t>(adj.node[static_cast<std::size_t>(k)]);
      const double candidate = d + adj.length[static_cast<std::size_t>(k)] / reach[u];
      if (candidate < best[v]) {
        best[v] = candidate;
        reach[v] = reach[u];
        queue.emplace(candidate, static_cast<std::int64_t>(v));
      }
    }
  }
  return best;
}

}  // namespace

DisplacementExtension extend_displacement(const DisplacementExtensionInput& in) {
  const std::size_t n = in.points.size();
  if (in.fixed.size() != n || in.displacement.size() != n || in.radius.size() != n) {
    throw std::invalid_argument("extend_displacement: inconsistent sizes");
  }
  const Adjacency adj = build_adjacency(in);
  DisplacementExtension out;
  const std::vector<double> distance = region_distance(in, adj, out.seeds);
  // Unknowns: free nodes in the region, numbered in node order.
  std::vector<std::int64_t> unknown(n, -1);
  std::vector<std::size_t> nodes;
  for (std::size_t i = 0; i < n; ++i) {
    if (in.fixed[i] == 0 && distance[i] <= 1.0) {
      unknown[i] = static_cast<std::int64_t>(nodes.size());
      nodes.push_back(i);
    }
  }
  out.displacement.assign(n, Vec3{});
  for (std::size_t i = 0; i < n; ++i) {
    if (in.fixed[i] != 0) {
      out.displacement[i] = in.displacement[i];
    }
  }
  out.solved = static_cast<std::int64_t>(nodes.size());
  if (nodes.empty()) {
    return out;
  }
  const auto m = static_cast<Eigen::Index>(nodes.size());
  std::vector<Eigen::Triplet<double>> triplets;
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(m, 3);
  for (std::size_t r = 0; r < nodes.size(); ++r) {
    const std::size_t i = nodes[r];
    double diagonal = 0.0;
    for (auto k = adj.start[i]; k < adj.start[i + 1]; ++k) {
      const auto j = static_cast<std::size_t>(adj.node[static_cast<std::size_t>(k)]);
      const double l = adj.length[static_cast<std::size_t>(k)];
      const double w = std::pow(l, -in.exponent);
      diagonal += w;
      if (unknown[j] >= 0) {
        triplets.emplace_back(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(unknown[j]),
                              -w);
      } else if (in.fixed[j] != 0) {
        const Vec3& d = in.displacement[j];
        rhs(static_cast<Eigen::Index>(r), 0) += w * d.x;
        rhs(static_cast<Eigen::Index>(r), 1) += w * d.y;
        rhs(static_cast<Eigen::Index>(r), 2) += w * d.z;
      }
      // Free nodes outside the region keep zero displacement: no right-hand side term.
    }
    triplets.emplace_back(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(r), diagonal);
  }
  Eigen::SparseMatrix<double> matrix(m, m);
  matrix.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::MatrixXd solution;
  if (m <= kDirectLimit) {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(matrix);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("extend_displacement: the factorisation failed");
    }
    solution = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("extend_displacement: the solve failed");
    }
  } else {
    // Large 3D regions: Cholesky fill-in grows too fast; conjugate gradients with an
    // incomplete Cholesky preconditioner instead.
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                             Eigen::IncompleteCholesky<double>>
        solver;
    solver.setTolerance(kTolerance);
    solver.setMaxIterations(kMaxIterations);
    solver.compute(matrix);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("extend_displacement: the preconditioner failed");
    }
    solution = Eigen::MatrixXd::Zero(m, 3);
    for (Eigen::Index column = 0; column < 3; ++column) {
      solution.col(column) = solver.solve(rhs.col(column));
      out.iterations = std::max(out.iterations, static_cast<std::int64_t>(solver.iterations()));
      if (solver.info() != Eigen::Success) {
        throw std::runtime_error("extend_displacement: conjugate gradients did not converge in " +
                                 std::to_string(kMaxIterations) + " iterations");
      }
    }
  }
  for (std::size_t r = 0; r < nodes.size(); ++r) {
    const auto row = static_cast<Eigen::Index>(r);
    out.displacement[nodes[r]] = {solution(row, 0), solution(row, 1), solution(row, 2)};
  }
  return out;
}

}  // namespace og
