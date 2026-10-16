// GPU sweeps of the 3D group smoother (gpu/gpu.hpp): agreement with the CPU, determinism.
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

#include "opengrid/core/gpu/gpu.hpp"
#include "opengrid/core/structured/multiblock_elliptic3d.hpp"

namespace {

// Two blocks of a sheared box clustered towards k = 0 (first cell 1e-3 of the height).
og::Multiblock3DInput sheared_box(int ni, int nj, int nk) {
  og::Multiblock3DInput input;
  auto id = [&](int i, int j, int k) { return (static_cast<std::int64_t>(k) * nj + j) * ni + i; };
  const double q = 1.25;
  std::vector<double> z(static_cast<std::size_t>(nk), 0.0);
  for (int k = 1; k < nk; ++k) {
    z[static_cast<std::size_t>(k)] = z[static_cast<std::size_t>(k) - 1] + 1e-3 * std::pow(q, k - 1);
  }
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const double s = static_cast<double>(i) / (ni - 1);
        const double t = static_cast<double>(j) / (nj - 1);
        const double r = z[static_cast<std::size_t>(k)] / z.back();
        double x = s;
        double y = t + 0.1 * std::sin(std::numbers::pi * s) * (1.0 - r);
        if (i > 0 && i < ni - 1 && j > 0 && j < nj - 1 && k > 0 && k < nk - 1) {
          x += 0.2 * s * (1.0 - s) * r;  // shear, smoothed out by the sweeps
        }
        input.xyz.insert(input.xyz.end(), {x, y, r});
      }
    }
  }
  input.fixed.assign(input.xyz.size() / 3, 0);
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        if (i == 0 || j == 0 || k == 0 || i == ni - 1 || j == nj - 1 || k == nk - 1) {
          input.fixed[static_cast<std::size_t>(id(i, j, k))] = 1;
        }
      }
    }
  }
  const int half = ni / 2;
  for (const auto& [i0, i1] : {std::pair{0, half}, std::pair{half, ni - 1}}) {
    for (int k = 0; k < nk - 1; ++k) {
      for (int j = 0; j < nj - 1; ++j) {
        for (int i = i0; i < i1; ++i) {
          std::array<std::int64_t, 8> h{};
          for (unsigned c = 0; c < 8; ++c) {
            h[c] = id(i + static_cast<int>(c & 1U), j + static_cast<int>((c >> 1U) & 1U),
                      k + static_cast<int>((c >> 2U) & 1U));
          }
          input.hexes.push_back(h);
        }
      }
    }
  }
  return input;
}

class GpuModeGuard {
 public:
  explicit GpuModeGuard(og::GpuMode mode) : before_(og::gpu_mode()) {
    og::set_gpu_mode(mode);
  }
  ~GpuModeGuard() {
    og::set_gpu_mode(before_);
  }
  GpuModeGuard(const GpuModeGuard&) = delete;
  GpuModeGuard& operator=(const GpuModeGuard&) = delete;
  GpuModeGuard(GpuModeGuard&&) = delete;
  GpuModeGuard& operator=(GpuModeGuard&&) = delete;

 private:
  og::GpuMode before_;
};

}  // namespace

TEST(Gpu, InfoAndModes) {
  const og::GpuInfo info = og::gpu_info();
  if (!info.built) {
    EXPECT_FALSE(info.available);
    EXPECT_FALSE(info.reason.empty());
  }
  const GpuModeGuard off(og::GpuMode::Off);
  EXPECT_FALSE(og::use_gpu(1'000'000'000));
  og::set_gpu_mode(og::GpuMode::Auto);
  EXPECT_FALSE(og::use_gpu(og::kGpuMinNodes - 1));
  EXPECT_EQ(og::use_gpu(og::kGpuMinNodes), info.available);
}

TEST(Gpu, GroupSweepsAgreeWithTheCpuAndRepeat) {
  if (!og::gpu_info().available) {
    GTEST_SKIP() << "no GPU: " << og::gpu_info().reason;
  }
  const og::Multiblock3DInput input = sheared_box(25, 21, 23);
  og::MultiblockOptions options;
  options.tolerance = 1e-10;
  options.max_iterations = 4000;
  og::Multiblock3DResult cpu;
  {
    const GpuModeGuard off(og::GpuMode::Off);
    cpu = og::smooth_multiblock_3d(input, options);
  }
  ASSERT_TRUE(cpu.converged);
  og::Multiblock3DResult gpu;
  og::Multiblock3DResult again;
  {
    const GpuModeGuard on(og::GpuMode::On);
    gpu = og::smooth_multiblock_3d(input, options);
    again = og::smooth_multiblock_3d(input, options);
  }
  ASSERT_TRUE(gpu.converged);
  EXPECT_TRUE(std::any_of(gpu.warnings.begin(), gpu.warnings.end(), [](const std::string& w) {
    return w.find("on the GPU") != std::string::npos;
  }));
  EXPECT_EQ(gpu.xyz, again.xyz);  // deterministic on one device
  // Within a fraction of the local cell size of the converged CPU grid (the GPU stops at a
  // relative residual of 1e-5).
  const std::size_t n = cpu.xyz.size() / 3;
  double worst = 0.0;
  for (const auto& h : input.hexes) {
    for (unsigned c = 0; c < 8; ++c) {
      const auto g = static_cast<std::size_t>(h[c]);
      double cell = 1e300;
      for (const unsigned b : {1U, 2U, 4U}) {
        const auto o = static_cast<std::size_t>(h[c ^ b]);
        cell = std::min(cell, std::hypot(cpu.xyz[3 * g] - cpu.xyz[3 * o],
                                         cpu.xyz[3 * g + 1] - cpu.xyz[3 * o + 1],
                                         cpu.xyz[3 * g + 2] - cpu.xyz[3 * o + 2]));
      }
      const double d =
          std::hypot(gpu.xyz[3 * g] - cpu.xyz[3 * g], gpu.xyz[3 * g + 1] - cpu.xyz[3 * g + 1],
                     gpu.xyz[3 * g + 2] - cpu.xyz[3 * g + 2]);
      worst = std::max(worst, d / cell);
    }
  }
  EXPECT_LT(worst, 0.1) << n << " nodes";
}
