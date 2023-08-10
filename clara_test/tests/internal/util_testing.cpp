#include <vector>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;
using namespace clara::internal;

TEST(clara_util_test, N2MultiIdx) {
  const idx dims[] = {3, 4, 5};
  idx result[3];

  n2multiidx(10, 3, dims, result);
  EXPECT_EQ(result[0], 0);
  EXPECT_EQ(result[1], 2);
  EXPECT_EQ(result[2], 0);

  n2multiidx(24, 3, dims, result);
  EXPECT_EQ(result[0], 1);
  EXPECT_EQ(result[1], 0);
  EXPECT_EQ(result[2], 4);
}

TEST(UtilTest, MultiIdx2N) {
  const idx dims[] = {3, 4, 5};
  idx midx[] = {2, 3, 1};

  idx n = multiidx2n(midx, 3, dims);
  EXPECT_EQ(n, 56);

  idx midx2[] = {0, 3, 4};
  idx n2 = multiidx2n(midx2, 3, dims);
  EXPECT_EQ(n2, 19);
}
