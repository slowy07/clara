#include <cmath>
#include <complex>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_States_zero, AllTests) {
  idx n = 1;
  EXPECT_NEAR(0, norm(clara::st.z0 - clara::st.zero(n)), 1e-7);
  
  n = 2;
  EXPECT_NEAR(0, norm(kron(st.z0, st.z0) - clara::st.zero(n)), 1e-7);
}
