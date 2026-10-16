#include "opengrid/core/util/anderson.hpp"

#include <algorithm>
#include <array>
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
  df_.assign(depth_ * size_, 0.0F);
  dg_.assign(depth_ * size_, 0.0F);
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

bool AndersonAccelerator::step(const std::vector<double>& u, std::vector<double>& g) {
  // One pass for the residual, the new difference columns and the saved state; one pass for
  // all dot products (each still summed in index order within fixed chunks, then over the
  // chunks in order: the same values as separate dot products, for any thread count).
  const bool primed = primed_;
  const std::size_t slot = next_;
  float* df = &df_[slot * size_];
  float* dg = &dg_[slot * size_];
  chunked(size_, [&](std::size_t begin, std::size_t end) {
    for (std::size_t k = begin; k < end; ++k) {
      const double fk = g[k] - u[k];
      f_[k] = fk;
      if (primed) {
        df[k] = static_cast<float>(fk - last_f_[k]);
        dg[k] = static_cast<float>(g[k] - last_g_[k]);
      }
      last_f_[k] = fk;
      last_g_[k] = g[k];
    }
  });
  if (primed) {
    count_ = std::min(count_ + 1, depth_);
    next_ = (next_ + 1) % depth_;
  }
  primed_ = true;
  if (count_ == 0) {
    return false;
  }
  const std::size_t m = count_;
  const std::size_t chunks = (size_ + kChunk - 1) / kChunk;
  // Per chunk: m products with the new column (Gram row) and m with the residual.
  std::vector<double> partial(chunks * 2 * m, 0.0);
  chunked(size_, [&](std::size_t begin, std::size_t end) {
    double* out = &partial[(begin / kChunk) * 2 * m];
    for (std::size_t s = 0; s < m; ++s) {
      const float* column = &df_[s * size_];
      // Four independent partial sums (a single accumulator is a latency-bound chain), added
      // in a fixed order: the result is still the same for any thread count.
      std::array<double, 4> with_new{};
      std::array<double, 4> with_f{};
      std::size_t k = begin;
      for (; k + 4 <= end; k += 4) {
        for (std::size_t j = 0; j < 4; ++j) {
          const double value = column[k + j];
          with_new[j] += value * static_cast<double>(df[k + j]);
          with_f[j] += value * f_[k + j];
        }
      }
      for (; k < end; ++k) {
        const double value = column[k];
        with_new[0] += value * static_cast<double>(df[k]);
        with_f[0] += value * f_[k];
      }
      out[2 * s] = (with_new[0] + with_new[1]) + (with_new[2] + with_new[3]);
      out[2 * s + 1] = (with_f[0] + with_f[1]) + (with_f[2] + with_f[3]);
    }
  });
  std::vector<double> gamma(m, 0.0);
  std::vector<double> row(m, 0.0);
  for (std::size_t c = 0; c < chunks; ++c) {
    for (std::size_t s = 0; s < m; ++s) {
      row[s] += partial[(c * m + s) * 2];
      gamma[s] += partial[(c * m + s) * 2 + 1];
    }
  }
  if (primed) {
    for (std::size_t s = 0; s < m; ++s) {
      gram_[slot * depth_ + s] = row[s];
      gram_[s * depth_ + slot] = row[s];
    }
  }
  if (!solve(gamma)) {
    reset();
    return false;
  }
  // Column by column (each unknown still sees the columns in the same order, so the values
  // are those of the unknown-by-unknown loop; this form streams and vectorises).
  chunked(size_, [&](std::size_t begin, std::size_t end) {
    for (std::size_t s = 0; s < m; ++s) {
      const double factor = gamma[s];
      const float* column = &dg_[s * size_];
      for (std::size_t k = begin; k < end; ++k) {
        g[k] -= factor * static_cast<double>(column[k]);
      }
    }
  });
  return true;
}

bool AndersonAccelerator::guarded_step(const std::vector<double>& u, std::vector<double>& g,
                                       double residual) {
  if (!started_ && std::isfinite(residual)) {
    started_ = true;
    best_ = residual;
    lowest_ = residual;
    mark_ = residual;
    best_state_ = g;
  }
  if (!enabled_) {
    return false;
  }
  if (residual < 0.9 * mark_) {
    mark_ = residual;
    since_mark_ = 0;
  } else {
    since_mark_ += 1;
  }
  if (!started_) {
    return false;  // no finite residual yet: nothing to restore
  }
  if (!std::isfinite(residual) || residual > kBlowup * lowest_ || since_mark_ > kStall) {
    restarts_ += 1;
    reset();
    g = best_state_;
    best_ = lowest_;
    mark_ = lowest_;
    since_mark_ = 0;
    if (restarts_ >= kMaxRestarts) {
      enabled_ = false;
    }
    return true;
  }
  if (residual < lowest_) {
    lowest_ = residual;
    best_state_ = g;
  }
  // Restart the history when the residual grows well above its best value.
  if (residual > 4.0 * best_) {
    reset();
    best_ = residual;
  }
  best_ = std::min(best_, residual);
  return step(u, g);
}

std::string AndersonAccelerator::restart_note() const {
  if (restarts_ == 0) {
    return "";
  }
  std::string note = "acceleration diverged or stalled " + std::to_string(restarts_) +
                     " time(s); restarted from the best iterate";
  if (!enabled_) {
    note += " and continued without acceleration";
  }
  return note;
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
