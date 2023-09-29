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

TEST(clara_states_testing, MesState) {
  const States& states = States::get_instance();
  ket mes_state = states.mes();

  ket expected_mes_state = ket::Zero(4);
  expected_mes_state(0) = expected_mes_state(3) = 1.0 / std::sqrt(2.0);

  EXPECT_EQ(mes_state, expected_mes_state);
}

TEST(clara_states_testing, TestWStates) {
  const States& states = States::get_instance();
  ket w_state = states.W;
  ket expected_w_state = ket::Zero(8);
  expected_w_state(1) = expected_w_state(2) = expected_w_state(4) = 1.0 / std::sqrt(3.0);
  EXPECT_EQ(w_state, expected_w_state);
}

TEST(clara_states_testing, TestTwoQubit) {
  const States& states = States::get_instance();
  ket plus_state = states.plus(2);
  ket expected_plus_state = ket::Zero(4);
  std::complex<double> value = 0.5;
  expected_plus_state(0) = expected_plus_state(1) = expected_plus_state(2) =
      expected_plus_state(3) = value;

  EXPECT_EQ(plus_state, expected_plus_state);
}

TEST(clara_states_testing, OtherTest) {
  idx n = 1;
  EXPECT_NEAR(0, norm(clara::st.z0 - clara::st.zero(n)), 1e-7);

  n = 2;
  EXPECT_NEAR(0, norm(kron(st.z0, st.z0) - clara::st.zero(n)), 1e-7);
}
