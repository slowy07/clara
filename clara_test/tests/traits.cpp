#include <cmath>
#include <vector>
#include <list>
#include "../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_is_complex, AllTests) {
  EXPECT_TRUE(clara::is_complex<std::complex<double>>::value);
  EXPECT_TRUE(clara::is_complex<std::complex<int>>::value);
}

TEST(clara_is_matrix_expression, AllTests) {
  dmat A, B;
  int x{}, y{}, z{};
  
  EXPECT_TRUE(clara::is_matrix_expression<decltype(3 * A)>::value);
  EXPECT_TRUE(clara::is_matrix_expression<decltype(A + B)>::value);
  EXPECT_FALSE(clara::is_matrix_expression<decltype(x + y * z)>::value);
}
