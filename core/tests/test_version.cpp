#include <gtest/gtest.h>

#include "opengrid/core/version.hpp"

TEST(Version, IsNotEmpty) {
  EXPECT_FALSE(og::version().empty());
}

TEST(Add, SmallIntegers) {
  EXPECT_EQ(og::add(2, 3), 5);
  EXPECT_EQ(og::add(-7, 7), 0);
}

TEST(Add, SixtyFourBit) {
  const std::int64_t two_pow_40 = 1'099'511'627'776;
  const std::int64_t two_pow_41 = 2'199'023'255'552;
  EXPECT_EQ(og::add(two_pow_40, two_pow_40), two_pow_41);
}
