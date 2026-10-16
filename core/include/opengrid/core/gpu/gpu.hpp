// GPU compute (OpenCL), optional: device detection, the global mode, and the GPU sweeps of the
// 3D group smoother.
//
// Method: OpenCL 1.2 (Apple, NVIDIA, AMD and Intel GPUs; deprecated on macOS but available).
// The first use selects the GPU with the most compute units, builds the kernels and runs a
// self-test of the float-float arithmetic (a driver that reassociates or flushes it is not
// used). Positions are stored as float-float pairs (hi + lo, ~48 bits) and every other
// quantity is float32: the elliptic sweeps run in delta form, so the arithmetic acts on
// neighbour differences, whose float32 rounding is relative to the local cell size (the
// precision study in tools/gpu: 2e-7 of the cell height at a 1e-6 wall, where plain float32
// positions are 6 % off). No float64 on the device is needed. The residual of float32
// arithmetic levels off near 1e-6 of the first (1M-cell ball: 1.3e-6, within 0.9 % of the
// local cell size of the CPU's converged grid); a GPU run converges at max(tolerance,
// kGpuTolerance). Results are deterministic on one device (no atomics; reductions in a fixed
// order), not bit-identical to the CPU's.
// Status (11d): experimental, off by default. The device sweeps are 2.2x faster than the CPU's
// on the 1M-cell ball, but the iteration only 1.2x (Anderson mixing and the projections of
// sliding nodes stay on the CPU), and the device colouring (Jacobi within a colour) needs more
// iterations: to the CPU's accuracy the GPU run is slower. See docs/PROGRESS.md (11d).
// References: Khronos OpenCL 1.2 specification; T. J. Dekker, "A floating-point technique for
// extending the available precision", Numer. Math. 18 (1971) 224-242 (two-sum); D. E. Knuth,
// TAOCP vol. 2, 4.2.2.
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"

namespace og {

enum class GpuMode : std::uint8_t {
  Off,   // never (the default)
  Auto,  // when a device passed the self-test and the grid has >= kGpuMinNodes nodes
  On,    // whenever a device passed the self-test
};

constexpr std::int64_t kGpuMinNodes = 200000;  // below this the launches cost more than they save
constexpr double kGpuTolerance = 1e-5;         // stopping level above the float32 floor

void set_gpu_mode(GpuMode mode);
[[nodiscard]] GpuMode gpu_mode();

struct GpuInfo {
  bool built = false;      // compiled with OpenCL
  bool available = false;  // a device passed the self-test
  std::string platform;
  std::string device;
  std::string version;
  int compute_units = 0;
  bool float64 = false;  // cl_khr_fp64 (not needed)
  std::string reason;    // why not available
};

// Initialises the runtime on first call (thread-safe).
[[nodiscard]] GpuInfo gpu_info();

// Whether a smoother with `nodes` nodes should use the GPU under the current mode.
[[nodiscard]] bool use_gpu(std::int64_t nodes);

// The interior work of one iteration of the 3D group smoother (structured/
// multiblock_elliptic3d.cpp) on the GPU: frozen coefficients of the regular nodes (27-point
// stencil in each node's frame, wall sources given per node) and the coloured line relaxation
// in delta form. Positions are exchanged with the CPU every iteration (the CPU keeps the
// sliding, irregular and wall-control updates and the Anderson mixing).
struct Multiblock3DTopology {
  std::size_t nodes = 0;
  std::vector<std::int32_t> regular;                  // node ids
  std::vector<std::array<std::int32_t, 27>> stencil;  // per regular node, slot = 13 + i + 3 j + 9 k
  std::vector<std::int32_t> line_offset;              // per line: first entry in line_nodes
  std::vector<std::int32_t> line_size;                // interior nodes
  std::vector<std::int32_t> line_end0;
  std::vector<std::int32_t> line_end1;
  std::vector<std::int32_t> line_nodes;
  std::vector<std::int8_t> line_prev;  // per line node: direction 0..5 of the previous node
  std::vector<std::int8_t> line_next;
  std::vector<std::vector<std::int32_t>> colours;  // lines per colour, in order
};

class Multiblock3DSweeps {
 public:
  // Throws std::runtime_error when no device is available or a buffer cannot be made.
  explicit Multiblock3DSweeps(const Multiblock3DTopology& topology);
  ~Multiblock3DSweeps();
  Multiblock3DSweeps(const Multiblock3DSweeps&) = delete;
  Multiblock3DSweeps& operator=(const Multiblock3DSweeps&) = delete;
  Multiblock3DSweeps(Multiblock3DSweeps&&) = delete;
  Multiblock3DSweeps& operator=(Multiblock3DSweeps&&) = delete;

  // One iteration: x (all nodes) in, wall sources per regular node (in the node's frame axes),
  // relaxation factor; x out (lines relaxed) and the squared update per line.
  void sweep(std::vector<Vec3>& x, const std::vector<std::array<float, 3>>& sources, double omega,
             std::vector<double>& line_sum);

  [[nodiscard]] std::string device() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace og
