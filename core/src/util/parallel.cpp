#include "opengrid/core/util/parallel.hpp"

#include <stdexcept>

#if defined(OG_WITH_OPENMP)
#include <omp.h>
#endif

namespace og {

int thread_count() {
#if defined(OG_WITH_OPENMP)
  return omp_get_max_threads();
#else
  return 1;
#endif
}

void set_thread_count(int n) {
  if (n < 1) {
    throw std::invalid_argument("set_thread_count: need at least one thread");
  }
#if defined(OG_WITH_OPENMP)
  omp_set_num_threads(n);
#endif
}

bool openmp_enabled() {
#if defined(OG_WITH_OPENMP)
  return true;
#else
  return false;
#endif
}

}  // namespace og
