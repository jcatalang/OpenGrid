// Zeroed scratch arrays for hot inner loops: on the stack up to kStack doubles, else on the
// heap (geometry evaluation runs millions of times inside the smoothers' projections, and a
// heap allocation per call cost half of their busy time).
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace og {

class Scratch {
 public:
  static constexpr std::size_t kStack = 256;

  explicit Scratch(std::size_t n) {
    if (n <= kStack) {
      data_ = stack_.data();
    } else {
      heap_.assign(n, 0.0);
      data_ = heap_.data();
    }
    std::fill(data_, data_ + n, 0.0);
  }
  Scratch(const Scratch&) = delete;
  Scratch& operator=(const Scratch&) = delete;
  Scratch(Scratch&&) = delete;
  Scratch& operator=(Scratch&&) = delete;
  ~Scratch() = default;

  [[nodiscard]] double* data() {
    return data_;
  }

 private:
  std::array<double, kStack> stack_;
  std::vector<double> heap_;
  double* data_ = nullptr;
};

}  // namespace og
