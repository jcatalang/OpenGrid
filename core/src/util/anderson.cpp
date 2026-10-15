#include "opengrid/core/util/anderson.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "opengrid/core/util/parallel.hpp"

namespace og {

namespace {

constexpr std::size_t kChunk = 16384;

// body(begin, end) over [0, size) in fixed chunks (the partition never depends on threads).
template <class Body>
void chunked(std::size_t size, Body&& body) {
  const std::size_t chunks = (size + kChunk - 1) / kChunk;
  parallel_for(
      0, static_cast<std::int64_t>(chunks),
      [&](std::int64_t c) {
        const std::size_t begin = static_cast<std::size_t>(c) * kChunk;
        body(begin, std::min(size, begin + kChunk));
      },
      2);
}

}  // namespace

AndersonAccelerator::AndersonAccelerator(std::size_t size, int depth)
    : size_(size), depth_(static_cast<std::size_t>(std::max(depth, 1))) {
  df_.assign(depth_ * size_, 0.0);
  dg_.assign(depth_ * size_, 0.0);
  gram_.assign(depth_ * depth_, 0.0);
  last_f_.assign(size_, 0.0);
  last_g_.assign(size_, 0.0);
  f_.assign(size_, 0.0);
}

void AndersonAccelerator::reset() {
  count_ = 0;
  next_ = 0;
  primed_ = false;
}

double AndersonAccelerator::dot(const double* a, const double* b) const {
  const std::size_t chunks = (size_ + kChunk - 1) / kChunk;
  std::vector<double> partial(chunks, 0.0);
  chunked(size_, [&](std::size_t begin, std::size_t end) {
    double sum = 0.0;
    for (std::size_t k = begin; k < end; ++k) {
      sum += a[k] * b[k];
    }
    partial[begin / kChunk] = sum;
  });
  double total = 0.0;
  for (const double value : partial) {
    total += value;
  }
  return total;
}

bool AndersonAccelerator::step(const std::vector<double>& u, std::vector<double>& g) {
  chunked(size_, [&](std::size_t begin, std::size_t end) {
    for (std::size_t k = begin; k < end; ++k) {
      f_[k] = g[k] - u[k];
    }
  });
  if (primed_) {
    const std::size_t slot = next_;
    double* df = &df_[slot * size_];
    double* dg = &dg_[slot * size_];
    chunked(size_, [&](std::size_t begin, std::size_t end) {
      for (std::size_t k = begin; k < end; ++k) {
        df[k] = f_[k] - last_f_[k];
        dg[k] = g[k] - last_g_[k];
      }
    });
    count_ = std::min(count_ + 1, depth_);
    next_ = (next_ + 1) % depth_;
    for (std::size_t s = 0; s < count_; ++s) {
      const double value = dot(&df_[s * size_], df);
      gram_[slot * depth_ + s] = value;
      gram_[s * depth_ + slot] = value;
    }
  }
  std::copy(f_.begin(), f_.end(), last_f_.begin());
  std::copy(g.begin(), g.end(), last_g_.begin());
  primed_ = true;
  if (count_ == 0) {
    return false;
  }
  std::vector<double> gamma(count_);
  for (std::size_t s = 0; s < count_; ++s) {
    gamma[s] = dot(&df_[s * size_], f_.data());
  }
  if (!solve(gamma)) {
    reset();
    return false;
  }
  chunked(size_, [&](std::size_t begin, std::size_t end) {
    for (std::size_t k = begin; k < end; ++k) {
      double value = g[k];
      for (std::size_t s = 0; s < count_; ++s) {
        value -= gamma[s] * dg_[s * size_ + k];
      }
      g[k] = value;
    }
  });
  return true;
}

bool AndersonAccelerator::solve(std::vector<double>& gamma) const {
  const std::size_t m = count_;
  double scale = 0.0;
  for (std::size_t s = 0; s < m; ++s) {
    scale = std::max(scale, gram_[s * depth_ + s]);
  }
  if (!(scale > 0.0)) {
    return false;
  }
  std::vector<double> a(m * m);
  for (std::size_t r = 0; r < m; ++r) {
    for (std::size_t c = 0; c < m; ++c) {
      a[r * m + c] = gram_[r * depth_ + c];
    }
    a[r * m + r] += 1e-10 * scale;
  }
  for (std::size_t c = 0; c < m; ++c) {
    std::size_t pivot = c;
    for (std::size_t r = c + 1; r < m; ++r) {
      if (std::fabs(a[r * m + c]) > std::fabs(a[pivot * m + c])) {
        pivot = r;
      }
    }
    if (!(std::fabs(a[pivot * m + c]) > 0.0)) {
      return false;
    }
    if (pivot != c) {
      for (std::size_t k = 0; k < m; ++k) {
        std::swap(a[c * m + k], a[pivot * m + k]);
      }
      std::swap(gamma[c], gamma[pivot]);
    }
    for (std::size_t r = c + 1; r < m; ++r) {
      const double w = a[r * m + c] / a[c * m + c];
      for (std::size_t k = c; k < m; ++k) {
        a[r * m + k] -= w * a[c * m + k];
      }
      gamma[r] -= w * gamma[c];
    }
  }
  for (std::size_t c = m; c-- > 0;) {
    double value = gamma[c];
    for (std::size_t k = c + 1; k < m; ++k) {
      value -= a[c * m + k] * gamma[k];
    }
    gamma[c] = value / a[c * m + c];
  }
  return std::all_of(gamma.begin(), gamma.end(), [](double v) { return std::isfinite(v); });
}

}  // namespace og
