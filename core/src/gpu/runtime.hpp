// OpenCL runtime shared by the GPU algorithms (private to og_core; built with OG_WITH_OPENCL).
#pragma once

#if defined(OG_WITH_OPENCL)

#if defined(__APPLE__)
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include <string>

#include "opengrid/core/gpu/gpu.hpp"

namespace og::gpu {

struct Runtime {
  cl_context context = nullptr;
  cl_device_id device = nullptr;
  cl_program program = nullptr;
  GpuInfo info;
};

// The process-wide runtime (initialised once); nullptr if no device passed the self-test.
const Runtime* runtime();

// Throws std::runtime_error("<what>: OpenCL error <code>") unless code is CL_SUCCESS.
void check(cl_int code, const char* what);

// clSetKernelArg for one value (buffers, scalars).
template <class T>
void set_arg(cl_kernel kernel, cl_uint index, const T& value) {
  check(clSetKernelArg(kernel, index, sizeof(T), static_cast<const void*>(&value)),
        "kernel argument");
}

}  // namespace og::gpu

#endif
