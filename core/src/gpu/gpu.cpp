#include "opengrid/core/gpu/gpu.hpp"

#include <atomic>
#include <cmath>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(OG_WITH_OPENCL)
#include "kernels.hpp"
#include "runtime.hpp"
#endif

namespace og {

namespace {

std::atomic<GpuMode> g_mode{GpuMode::Off};

}  // namespace

void set_gpu_mode(GpuMode mode) {
  g_mode.store(mode);
}

GpuMode gpu_mode() {
  return g_mode.load();
}

#if defined(OG_WITH_OPENCL)

namespace gpu {

void check(cl_int code, const char* what) {
  if (code != CL_SUCCESS) {
    throw std::runtime_error(std::string(what) + ": OpenCL error " + std::to_string(code));
  }
}

namespace {

std::string device_string(cl_device_id device, cl_device_info what) {
  std::size_t size = 0;
  clGetDeviceInfo(device, what, 0, nullptr, &size);
  std::string out(size, '\0');
  clGetDeviceInfo(device, what, size, out.data(), nullptr);
  while (!out.empty() && out.back() == '\0') {
    out.pop_back();
  }
  return out;
}

std::string platform_string(cl_platform_id platform, cl_platform_info what) {
  std::size_t size = 0;
  clGetPlatformInfo(platform, what, 0, nullptr, &size);
  std::string out(size, '\0');
  clGetPlatformInfo(platform, what, size, out.data(), nullptr);
  while (!out.empty() && out.back() == '\0') {
    out.pop_back();
  }
  return out;
}

// Float-float arithmetic as the kernels do it, against float64: a reassociating or flushing
// driver fails here.
void self_test(const Runtime& rt) {
  cl_int err = CL_SUCCESS;
  cl_command_queue queue = clCreateCommandQueue(rt.context, rt.device, 0, &err);
  check(err, "self-test queue");
  constexpr std::size_t kCount = 64;
  std::vector<double> a(std::size_t{6} * kCount);
  std::vector<float> hi(a.size());
  std::vector<float> lo(a.size());
  std::vector<float> d(3 * kCount);
  for (std::size_t k = 0; k < a.size(); ++k) {
    a[k] = 1.0 + 1e-3 * std::sin(static_cast<double>(k)) + 1e-9 * static_cast<double>(k);
    hi[k] = static_cast<float>(a[k]);
    lo[k] = static_cast<float>(a[k] - static_cast<double>(hi[k]));
  }
  for (std::size_t k = 0; k < d.size(); ++k) {
    d[k] = static_cast<float>(1e-7 * std::cos(static_cast<double>(k)));
  }
  auto bytes = [](const std::vector<float>& v) { return v.size() * sizeof(float); };
  cl_mem bh = clCreateBuffer(rt.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bytes(hi),
                             hi.data(), &err);
  check(err, "self-test buffer");
  cl_mem bl = clCreateBuffer(rt.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bytes(lo),
                             lo.data(), &err);
  check(err, "self-test buffer");
  cl_mem bd =
      clCreateBuffer(rt.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, bytes(d), d.data(), &err);
  check(err, "self-test buffer");
  std::vector<float> diff(3 * kCount);
  cl_mem bdiff = clCreateBuffer(rt.context, CL_MEM_WRITE_ONLY, bytes(diff), nullptr, &err);
  check(err, "self-test buffer");
  cl_kernel kernel = clCreateKernel(rt.program, "fftest", &err);
  check(err, "self-test kernel");
  set_arg(kernel, 0, bh);
  set_arg(kernel, 1, bl);
  set_arg(kernel, 2, bd);
  set_arg(kernel, 3, bdiff);
  const std::size_t global = kCount;
  check(clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr),
        "self-test run");
  check(clEnqueueReadBuffer(queue, bh, CL_TRUE, 0, bytes(hi), hi.data(), 0, nullptr, nullptr),
        "self-test read");
  check(clEnqueueReadBuffer(queue, bl, CL_TRUE, 0, bytes(lo), lo.data(), 0, nullptr, nullptr),
        "self-test read");
  check(
      clEnqueueReadBuffer(queue, bdiff, CL_TRUE, 0, bytes(diff), diff.data(), 0, nullptr, nullptr),
      "self-test read");
  double worst_sum = 0.0;
  double worst_diff = 0.0;
  for (std::size_t q = 0; q < kCount; ++q) {
    for (std::size_t c = 0; c < 3; ++c) {
      const std::size_t i = (6 * q) + c;      // node 2q, component c
      const std::size_t j = (6 * q) + 3 + c;  // node 2q + 1
      const std::size_t o = (3 * q) + c;
      const double sum = a[i] + static_cast<double>(d[o]);
      const double got = static_cast<double>(hi[i]) + static_cast<double>(lo[i]);
      worst_sum = std::fmax(worst_sum, std::fabs(got - sum));
      const double want = a[i] - a[j];
      worst_diff = std::fmax(worst_diff, std::fabs(static_cast<double>(diff[o]) - want) /
                                             std::fmax(std::fabs(want), 1e-30));
    }
  }
  clReleaseKernel(kernel);
  clReleaseMemObject(bh);
  clReleaseMemObject(bl);
  clReleaseMemObject(bd);
  clReleaseMemObject(bdiff);
  clReleaseCommandQueue(queue);
  // Float-float sums carry ~48 bits (1e-14 near 1); differences one float32 rounding.
  if (!(worst_sum < 1e-13) || !(worst_diff < 1e-6)) {
    throw std::runtime_error("float-float self-test failed (sum error " +
                             std::to_string(worst_sum) + ", difference error " +
                             std::to_string(worst_diff) + ")");
  }
}

Runtime* create() {
  auto* rt = new Runtime();
  rt->info.built = true;
  try {
    cl_uint platforms = 0;
    if (clGetPlatformIDs(0, nullptr, &platforms) != CL_SUCCESS || platforms == 0) {
      throw std::runtime_error("no OpenCL platform");
    }
    std::vector<cl_platform_id> ids(platforms);
    check(clGetPlatformIDs(platforms, ids.data(), nullptr), "platforms");
    cl_uint best_units = 0;
    for (cl_platform_id platform : ids) {
      cl_uint count = 0;
      if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &count) != CL_SUCCESS) {
        continue;
      }
      std::vector<cl_device_id> devices(count);
      clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, count, devices.data(), nullptr);
      for (cl_device_id device : devices) {
        cl_uint units = 0;
        clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(units), &units, nullptr);
        if (units > best_units) {
          best_units = units;
          rt->device = device;
          rt->info.platform = platform_string(platform, CL_PLATFORM_NAME);
        }
      }
    }
    if (rt->device == nullptr) {
      throw std::runtime_error("no OpenCL GPU device");
    }
    rt->info.device = device_string(rt->device, CL_DEVICE_NAME);
    rt->info.version = device_string(rt->device, CL_DEVICE_VERSION);
    rt->info.compute_units = static_cast<int>(best_units);
    rt->info.float64 =
        device_string(rt->device, CL_DEVICE_EXTENSIONS).find("cl_khr_fp64") != std::string::npos;
    cl_int err = CL_SUCCESS;
    rt->context = clCreateContext(nullptr, 1, &rt->device, nullptr, nullptr, &err);
    check(err, "context");
    const char* source = kKernels;
    rt->program = clCreateProgramWithSource(rt->context, 1, &source, nullptr, &err);
    check(err, "program");
    if (clBuildProgram(rt->program, 1, &rt->device, "-cl-std=CL1.2", nullptr, nullptr) !=
        CL_SUCCESS) {
      std::size_t size = 0;
      clGetProgramBuildInfo(rt->program, rt->device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &size);
      std::string log(size, '\0');
      clGetProgramBuildInfo(rt->program, rt->device, CL_PROGRAM_BUILD_LOG, size, log.data(),
                            nullptr);
      throw std::runtime_error("kernel build failed: " + log.substr(0, 400));
    }
    self_test(*rt);
    rt->info.available = true;
  } catch (const std::exception& error) {
    rt->info.available = false;
    rt->info.reason = error.what();
  }
  return rt;
}

Runtime* instance() {
  static std::once_flag once;
  static Runtime* rt = nullptr;
  std::call_once(once, [] { rt = create(); });
  return rt;
}

}  // namespace

const Runtime* runtime() {
  const Runtime* rt = instance();
  if (!rt->info.available) {
    return nullptr;
  }
  return rt;
}

}  // namespace gpu

GpuInfo gpu_info() {
  return gpu::instance()->info;
}

#else

GpuInfo gpu_info() {
  GpuInfo info;
  info.reason = "built without OpenCL";
  return info;
}

#endif

bool use_gpu(std::int64_t nodes) {
  const GpuMode mode = gpu_mode();
  if (mode == GpuMode::Off) {
    return false;
  }
  if (mode == GpuMode::Auto && nodes < kGpuMinNodes) {
    return false;
  }
  return gpu_info().available;
}

}  // namespace og
