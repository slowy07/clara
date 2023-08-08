#include <cmath>
#include <complex>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

const double TOLERANCE = 1e-9;

bool complexVectorsAlmostEqual(const clara::ket& v1, const clara::ket& v2, double tol) {
  if (v1.size() != v2.size()) {
    return false;
  }

  for (clara::idx i = 0; i < v1.size(); ++i) {
    if (std::abs(v1[i].real() - v2[i].real()) > tol ||
        std::abs(v1[i].imag() - v2[i].imag()) > tol) {
      return false;
    }
  }
  return true;
}

TEST(clara_states_testing, SingletonInstance) {
  const States& states = States::get_instance();
  const States& states2 = States::get_instance();
  EXPECT_EQ(&states, &states2);
}

TEST(clara_states_testing, OtherTest) {
  idx n = 1;
  EXPECT_NEAR(0, norm(clara::st.z0 - clara::st.zero(n)), 1e-7);

  n = 2;
  EXPECT_NEAR(0, norm(kron(st.z0, st.z0) - clara::st.zero(n)), 1e-7);
}
