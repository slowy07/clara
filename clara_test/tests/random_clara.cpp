#include <algorithm>
#include <climits>

#include "../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_rand_integer, AllTests) {
  bigint a = 42, b = a;
  EXPECT_EQ(a, clara::rand(a, b));

  a = 1;
  b = 5;
  for (int i = 0; i < 1000; ++i) {
    int n = clara::rand(a, b);
    EXPECT_GE(n, a);
    EXPECT_LE(n, b);
  }
}

TEST(clara_rand_min_max, AllTests) {
  bigint a = INT_MIN, b = INT_MAX;
  for (int i = 0; i < 1000; ++i) {
    int n = clara::rand(a, b);
    EXPECT_GE(n, a);
    EXPECT_LE(n, b);
  }
}

TEST(clara_rand_double, AllTest) {
  idx N = 1000;
  double a = 0, b = 1;
  for (idx i = 0; i < N; ++i) {
    double n = clara::rand(a, b);
    EXPECT_GE(n, a);
    EXPECT_LE(n, b);
  }

  N = 10000;
  a = -10, b = 10;
  double average = 0;
  for (idx i = 0; i < N; ++i) {
    double n = clara::rand(a, b);
    EXPECT_GE(n, a);
    EXPECT_LE(n, b);
    average += n;
  }
  average /= N;
  EXPECT_NEAR(0, average, 2e-1);
}
