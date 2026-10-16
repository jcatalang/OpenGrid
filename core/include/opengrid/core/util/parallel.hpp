// Deterministic parallel loops.
//
// Method: OpenMP work sharing with dynamic scheduling in chunks (about kChunksPerThread per
// thread: iterations of very different cost, e.g. sliding nodes that project onto NURBS
// surfaces next to plane ones, are balanced; a static partition left 69 % of the threads idle
// in the 3D group smoother) when the library is built with OG_WITH_OPENMP, a plain loop
// otherwise. The result never depends on the number of threads or on the schedule: every
// iteration must write only its own outputs, and reductions are done afterwards in index order
// by the caller (never with `reduction(...)` clauses, whose summation order follows the thread
// partition). Bodies must not throw: validate input before the loop.
#pragma once

#include <cstdint>
#if defined(OG_WITH_OPENMP)
#include <omp.h>
#endif

namespace og {

// Number of threads parallel loops use (1 without OpenMP).
int thread_count();

// Set the number of threads (n >= 1); no effect without OpenMP.
void set_thread_count(int n);

// Whether the library was built with OpenMP.
bool openmp_enabled();

// Calls body(i) for i in [begin, end). Ranges shorter than `grain` run serially (the thread
// start-up costs a few microseconds).
template <class Body>
void parallel_for(std::int64_t begin, std::int64_t end, Body&& body, std::int64_t grain = 1024) {
#if defined(OG_WITH_OPENMP)
  if (end - begin >= grain) {
    constexpr std::int64_t kChunksPerThread = 16;
    const std::int64_t chunk =
        1 + (end - begin) / (kChunksPerThread * static_cast<std::int64_t>(omp_get_max_threads()));
#pragma omp parallel for schedule(dynamic, chunk)
    for (std::int64_t i = begin; i < end; ++i) {
      body(i);
    }
    return;
  }
#endif
  for (std::int64_t i = begin; i < end; ++i) {
    body(i);
  }
}

}  // namespace og
