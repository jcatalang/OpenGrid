// GPU sweeps of the 3D group smoother (see opengrid/core/gpu/gpu.hpp): the frozen-coefficient
// pass and the coloured line relaxation, in delta form with float-float positions.
//
// Method: one work-item per regular node computes its 27-point coefficients from float-float
// neighbour differences; one work-item per line of a colour assembles and solves its
// tridiagonal system for the update (Thomas, scratch in global memory) and adds it to the
// positions in float-float. The device colours the lines itself, so that lines of one colour
// only share no node (a few colours of thousands of lines: one work-item per line needs that
// many to hide memory latency); a colour reads the positions at its start and writes into a
// copy (double buffering), so lines of a colour see each other's stencil nodes unchanged -
// Jacobi within a colour, Gauss-Seidel across colours (the CPU colours lines that share stencil
// nodes and gets ~20 colours). Positions are converted between float64 (CPU) and float-float
// (device) every iteration on the CPU threads, through mapped buffers (no copy on unified
// memory).
#include <array>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "opengrid/core/gpu/gpu.hpp"
#include "opengrid/core/util/parallel.hpp"

#if defined(OG_WITH_OPENCL)
#include "runtime.hpp"
#endif

namespace og {

#if defined(OG_WITH_OPENCL)

struct Multiblock3DSweeps::Impl {
  const gpu::Runtime* rt = nullptr;
  cl_command_queue queue = nullptr;
  cl_kernel coefficients = nullptr;
  cl_kernel relax = nullptr;
  std::vector<cl_mem> buffers;
  cl_mem regular = nullptr;
  cl_mem stencil = nullptr;
  cl_mem coef_of = nullptr;
  cl_mem src = nullptr;
  std::array<cl_mem, 2> xh{};
  std::array<cl_mem, 2> xl{};
  cl_mem coef = nullptr;
  cl_mem lines = nullptr;
  cl_mem line_offset = nullptr;
  cl_mem line_size = nullptr;
  cl_mem line_nodes = nullptr;
  cl_mem line_prev = nullptr;
  cl_mem line_next = nullptr;
  cl_mem scratch_c = nullptr;
  cl_mem scratch_d = nullptr;
  cl_mem line_sum = nullptr;
  std::size_t nodes = 0;
  std::size_t regular_count = 0;
  std::size_t line_count = 0;
  std::vector<std::pair<int, int>> colours;  // (first, count) into lines

  cl_mem make(std::size_t bytes, const void* data, cl_mem_flags flags) {
    cl_int err = CL_SUCCESS;
    std::size_t size = bytes;
    if (size == 0) {
      size = 4;  // OpenCL refuses empty buffers
    }
    cl_mem_flags all = flags;
    void* host = nullptr;
    if (data != nullptr && bytes > 0) {
      all |= CL_MEM_COPY_HOST_PTR;
      host = const_cast<void*>(data);  // NOLINT(cppcoreguidelines-pro-type-const-cast)
    }
    cl_mem mem = clCreateBuffer(rt->context, all, size, host, &err);
    gpu::check(err, "GPU buffer");
    buffers.push_back(mem);
    return mem;
  }

  ~Impl() {
    for (cl_mem mem : buffers) {
      clReleaseMemObject(mem);
    }
    if (coefficients != nullptr) {
      clReleaseKernel(coefficients);
    }
    if (relax != nullptr) {
      clReleaseKernel(relax);
    }
    if (queue != nullptr) {
      clReleaseCommandQueue(queue);
    }
  }
};

Multiblock3DSweeps::Multiblock3DSweeps(const Multiblock3DTopology& t)
    : impl_(std::make_unique<Impl>()) {
  Impl& m = *impl_;
  m.rt = gpu::runtime();
  if (m.rt == nullptr) {
    throw std::runtime_error("no GPU available: " + gpu_info().reason);
  }
  cl_int err = CL_SUCCESS;
  m.queue = clCreateCommandQueue(m.rt->context, m.rt->device, 0, &err);
  gpu::check(err, "GPU queue");
  m.coefficients = clCreateKernel(m.rt->program, "coefficients", &err);
  gpu::check(err, "GPU kernel coefficients");
  m.relax = clCreateKernel(m.rt->program, "relax_lines", &err);
  gpu::check(err, "GPU kernel relax_lines");
  m.nodes = t.nodes;
  m.regular_count = t.regular.size();
  m.line_count = t.line_offset.size();
  std::vector<std::int32_t> coef_of(t.nodes, -1);
  for (std::size_t r = 0; r < t.regular.size(); ++r) {
    coef_of[static_cast<std::size_t>(t.regular[r])] = static_cast<std::int32_t>(r);
  }
  // Colours for the device: lines conflict when they share an interior node (greedy, in the
  // order of the CPU colours, so it is deterministic).
  std::vector<std::vector<std::int32_t>> lines_at(t.nodes);
  std::vector<std::int32_t> order;
  for (const std::vector<std::int32_t>& colour : t.colours) {
    order.insert(order.end(), colour.begin(), colour.end());
  }
  for (const std::int32_t l : order) {
    const auto first = static_cast<std::size_t>(t.line_offset[static_cast<std::size_t>(l)]);
    const auto size = static_cast<std::size_t>(t.line_size[static_cast<std::size_t>(l)]);
    for (std::size_t q = first; q < first + size; ++q) {
      lines_at[static_cast<std::size_t>(t.line_nodes[q])].push_back(l);
    }
  }
  std::vector<int> colour_of(t.line_offset.size(), -1);
  std::vector<std::vector<std::int32_t>> by_colour;
  for (const std::int32_t l : order) {
    std::vector<std::uint8_t> used(by_colour.size() + 1, 0);
    const auto first = static_cast<std::size_t>(t.line_offset[static_cast<std::size_t>(l)]);
    const auto size = static_cast<std::size_t>(t.line_size[static_cast<std::size_t>(l)]);
    for (std::size_t q = first; q < first + size; ++q) {
      for (const std::int32_t other : lines_at[static_cast<std::size_t>(t.line_nodes[q])]) {
        const int c = colour_of[static_cast<std::size_t>(other)];
        if (other != l && c >= 0) {
          used[static_cast<std::size_t>(c)] = 1;
        }
      }
    }
    std::size_t c = 0;
    while (c < by_colour.size() && used[c] != 0) {
      ++c;
    }
    if (c == by_colour.size()) {
      by_colour.emplace_back();
    }
    by_colour[c].push_back(l);
    colour_of[static_cast<std::size_t>(l)] = static_cast<int>(c);
  }
  std::vector<std::int32_t> lines;
  for (const std::vector<std::int32_t>& colour : by_colour) {
    m.colours.emplace_back(static_cast<int>(lines.size()), static_cast<int>(colour.size()));
    lines.insert(lines.end(), colour.begin(), colour.end());
  }
  const auto i32 = sizeof(std::int32_t);
  m.regular = m.make(t.regular.size() * i32, t.regular.data(), CL_MEM_READ_ONLY);
  m.stencil = m.make(t.stencil.size() * 27 * i32, t.stencil.data(), CL_MEM_READ_ONLY);
  m.coef_of = m.make(coef_of.size() * i32, coef_of.data(), CL_MEM_READ_ONLY);
  m.src = m.make(3 * m.regular_count * sizeof(float), nullptr,
                 CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR);
  for (std::size_t b = 0; b < 2; ++b) {
    m.xh[b] =
        m.make(3 * m.nodes * sizeof(float), nullptr, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR);
    m.xl[b] =
        m.make(3 * m.nodes * sizeof(float), nullptr, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR);
  }
  m.coef = m.make(16 * m.regular_count * sizeof(float), nullptr, CL_MEM_READ_WRITE);
  m.lines = m.make(lines.size() * i32, lines.data(), CL_MEM_READ_ONLY);
  m.line_offset = m.make(t.line_offset.size() * i32, t.line_offset.data(), CL_MEM_READ_ONLY);
  m.line_size = m.make(t.line_size.size() * i32, t.line_size.data(), CL_MEM_READ_ONLY);
  m.line_nodes = m.make(t.line_nodes.size() * i32, t.line_nodes.data(), CL_MEM_READ_ONLY);
  m.line_prev = m.make(t.line_prev.size(), t.line_prev.data(), CL_MEM_READ_ONLY);
  m.line_next = m.make(t.line_next.size(), t.line_next.data(), CL_MEM_READ_ONLY);
  m.scratch_c = m.make(t.line_nodes.size() * sizeof(float), nullptr, CL_MEM_READ_WRITE);
  m.scratch_d = m.make(3 * t.line_nodes.size() * sizeof(float), nullptr, CL_MEM_READ_WRITE);
  m.line_sum =
      m.make(m.line_count * sizeof(float), nullptr, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR);
  // Fixed kernel arguments.
  const auto count = static_cast<cl_int>(m.regular_count);
  gpu::set_arg(m.coefficients, 0, count);
  gpu::set_arg(m.coefficients, 1, m.regular);
  gpu::set_arg(m.coefficients, 2, m.stencil);
  gpu::set_arg(m.coefficients, 3, m.src);
  gpu::set_arg(m.coefficients, 6, m.coef);
  gpu::set_arg(m.relax, 2, m.lines);
  gpu::set_arg(m.relax, 3, m.line_offset);
  gpu::set_arg(m.relax, 4, m.line_size);
  gpu::set_arg(m.relax, 5, m.line_nodes);
  gpu::set_arg(m.relax, 6, m.line_prev);
  gpu::set_arg(m.relax, 7, m.line_next);
  gpu::set_arg(m.relax, 8, m.coef_of);
  gpu::set_arg(m.relax, 9, m.stencil);
  gpu::set_arg(m.relax, 10, m.coef);
  gpu::set_arg(m.relax, 15, m.scratch_c);
  gpu::set_arg(m.relax, 16, m.scratch_d);
  gpu::set_arg(m.relax, 18, m.line_sum);
}

Multiblock3DSweeps::~Multiblock3DSweeps() = default;

std::string Multiblock3DSweeps::device() const {
  return impl_->rt->info.device;
}

void Multiblock3DSweeps::sweep(std::vector<Vec3>& x,
                               const std::vector<std::array<float, 3>>& sources, double omega,
                               std::vector<double>& line_sum) {
  Impl& m = *impl_;
  if (x.size() != m.nodes || sources.size() != m.regular_count) {
    throw std::invalid_argument("GPU sweep: sizes do not match the topology");
  }
  const std::size_t bytes = 3 * m.nodes * sizeof(float);
  cl_int err = CL_SUCCESS;
  auto map = [&](cl_mem mem, cl_map_flags flags, std::size_t size) {
    void* ptr =
        clEnqueueMapBuffer(m.queue, mem, CL_TRUE, flags, 0, size, 0, nullptr, nullptr, &err);
    gpu::check(err, "GPU map");
    return static_cast<float*>(ptr);
  };
  auto unmap = [&](cl_mem mem, void* ptr) {
    gpu::check(clEnqueueUnmapMemObject(m.queue, mem, ptr, 0, nullptr, nullptr), "GPU unmap");
  };
  // Positions in: float64 -> float-float, written straight into the device buffers.
  float* hi = map(m.xh[0], CL_MAP_WRITE_INVALIDATE_REGION, bytes);
  float* lo = map(m.xl[0], CL_MAP_WRITE_INVALIDATE_REGION, bytes);
  parallel_for(
      0, static_cast<std::int64_t>(m.nodes),
      [&](std::int64_t g) {
        const auto k = 3 * static_cast<std::size_t>(g);
        const Vec3& p = x[static_cast<std::size_t>(g)];
        const std::array<double, 3> v = {p.x, p.y, p.z};
        for (std::size_t c = 0; c < 3; ++c) {
          const auto h = static_cast<float>(v[c]);
          hi[k + c] = h;
          lo[k + c] = static_cast<float>(v[c] - static_cast<double>(h));
        }
      },
      16384);
  unmap(m.xh[0], hi);
  unmap(m.xl[0], lo);
  if (m.regular_count > 0) {
    const std::size_t src_bytes = 3 * m.regular_count * sizeof(float);
    float* src = map(m.src, CL_MAP_WRITE_INVALIDATE_REGION, src_bytes);
    std::memcpy(src, sources.data(), src_bytes);
    unmap(m.src, src);
  }
  auto round_up = [](std::size_t value) { return (value + 63) / 64 * 64; };
  if (m.regular_count > 0) {
    gpu::set_arg(m.coefficients, 4, m.xh[0]);
    gpu::set_arg(m.coefficients, 5, m.xl[0]);
    const std::size_t global = round_up(m.regular_count);
    gpu::check(clEnqueueNDRangeKernel(m.queue, m.coefficients, 1, nullptr, &global, nullptr, 0,
                                      nullptr, nullptr),
               "GPU coefficients");
  }
  const auto w = static_cast<float>(omega);
  gpu::set_arg(m.relax, 17, w);
  std::size_t cur = 0;
  for (const auto& [first, count] : m.colours) {
    if (count == 0) {
      continue;
    }
    const std::size_t nxt = 1 - cur;
    gpu::check(clEnqueueCopyBuffer(m.queue, m.xh[cur], m.xh[nxt], 0, 0, bytes, 0, nullptr, nullptr),
               "GPU copy");
    gpu::check(clEnqueueCopyBuffer(m.queue, m.xl[cur], m.xl[nxt], 0, 0, bytes, 0, nullptr, nullptr),
               "GPU copy");
    const cl_int n = count;
    const cl_int f = first;
    gpu::set_arg(m.relax, 0, n);
    gpu::set_arg(m.relax, 1, f);
    gpu::set_arg(m.relax, 11, m.xh[cur]);
    gpu::set_arg(m.relax, 12, m.xl[cur]);
    gpu::set_arg(m.relax, 13, m.xh[nxt]);
    gpu::set_arg(m.relax, 14, m.xl[nxt]);
    const std::size_t global = round_up(static_cast<std::size_t>(count));
    gpu::check(
        clEnqueueNDRangeKernel(m.queue, m.relax, 1, nullptr, &global, nullptr, 0, nullptr, nullptr),
        "GPU line relaxation");
    cur = nxt;
  }
  // Positions out.
  const float* rh = map(m.xh[cur], CL_MAP_READ, bytes);
  const float* rl = map(m.xl[cur], CL_MAP_READ, bytes);
  parallel_for(
      0, static_cast<std::int64_t>(m.nodes),
      [&](std::int64_t g) {
        const auto k = 3 * static_cast<std::size_t>(g);
        x[static_cast<std::size_t>(g)] = {
            static_cast<double>(rh[k]) + static_cast<double>(rl[k]),
            static_cast<double>(rh[k + 1]) + static_cast<double>(rl[k + 1]),
            static_cast<double>(rh[k + 2]) + static_cast<double>(rl[k + 2])};
      },
      16384);
  unmap(m.xh[cur], const_cast<float*>(rh));  // NOLINT(cppcoreguidelines-pro-type-const-cast)
  unmap(m.xl[cur], const_cast<float*>(rl));  // NOLINT(cppcoreguidelines-pro-type-const-cast)
  line_sum.resize(m.line_count);
  if (m.line_count > 0) {
    const float* sums = map(m.line_sum, CL_MAP_READ, m.line_count * sizeof(float));
    for (std::size_t l = 0; l < m.line_count; ++l) {
      line_sum[l] = static_cast<double>(sums[l]);
    }
    unmap(m.line_sum, const_cast<float*>(sums));  // NOLINT(cppcoreguidelines-pro-type-const-cast)
  }
  gpu::check(clFinish(m.queue), "GPU sweep");
}

#else

struct Multiblock3DSweeps::Impl {};

Multiblock3DSweeps::Multiblock3DSweeps(const Multiblock3DTopology& /*topology*/) {
  throw std::runtime_error("built without OpenCL");
}

Multiblock3DSweeps::~Multiblock3DSweeps() = default;

std::string Multiblock3DSweeps::device() const {
  return "";
}

void Multiblock3DSweeps::sweep(std::vector<Vec3>& /*x*/,
                               const std::vector<std::array<float, 3>>& /*sources*/,
                               double /*omega*/, std::vector<double>& /*line_sum*/) {
  throw std::runtime_error("built without OpenCL");
}

#endif

}  // namespace og
