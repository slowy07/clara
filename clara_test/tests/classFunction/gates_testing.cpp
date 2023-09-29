#include <cmath>
#include <complex>
#include <vector>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

const double tolerance = 1e-9;

TEST(clara_gates_test, Initialization) {
  const Gates& gates = Gates::get_instance();
  EXPECT_EQ(gates.Id2, cmat::Identity(2, 2));
}

TEST(clara_gates_test, RotationGates) {
  const Gates& gates = Gates::get_instance();
  cmat rx_gate = gates.RX(0.5);
  cmat expected_rx_gate(2, 2);
  expected_rx_gate << std::cos(0.25), -1.0_i * std::sin(0.25), -1.0_i * std::sin(0.25),
      std::cos(0.25);
  EXPECT_TRUE(rx_gate.isApprox(expected_rx_gate, tolerance));
}

bool complexMatricesEqual(const cmat& A, const cmat& B, double tolerance) {
  return (A - B).norm() <= tolerance;
}

TEST(clara_gates_test, ZdGate) {
  const Gates& gates = Gates::get_instance();
  cmat expectedZ2 = cmat::Zero(2, 2);
  expectedZ2(0, 0) = 1.0;
  expectedZ2(1, 1) = -1.0;
  EXPECT_TRUE(complexMatricesEqual(gates.Zd(2), expectedZ2, tolerance));
}

TEST(clara_Gates_control, Qubits) {
  cmat CTRL1 = gt.CTRL(gt.X, {0}, {1}, 2);
  EXPECT_EQ(CTRL1, gt.CNOT);

  cmat CTRL2 = gt.CTRL(gt.X, {1}, {0}, 2);
  EXPECT_EQ(CTRL2, gt.CNOTba);

  cmat CTRL3 = gt.CTRL(gt.X, {0, 1}, {2}, 3);
  EXPECT_EQ(CTRL3, gt.TOF);
  CTRL3 = gt.CTRL(gt.X, {0, 1}, {2}, 3, 2);
  EXPECT_EQ(CTRL3, gt.TOF);

  cmat U = randU(2);
  cmat CTRL4 = gt.CTRL(U, {0, 2}, {1}, 3);
  ket psi1 = mket({0, 0, 1});
  ket res1 = mket({0, 0, 1});
  EXPECT_NEAR(0, norm(CTRL4 * psi1 - res1), 1e-7);

  ket psi2 = mket({1, 1, 1});
  ket res2 = kron(st.z1, U * st.z1, st.z1);
  EXPECT_NEAR(0, norm(CTRL4 * psi2 - res2), 1e-7);
}

TEST(clara_Gates_Xd, AllTest) {
  for (idx D = 1; D < 10; ++D) {
    cmat Xd = gt.Xd(D);
    for (idx i = 0; i < D; ++i) {
      ket psi = mket({i}, D);
      ket res = mket({(i + 1) % D}, D);
      EXPECT_NEAR(0, norm(res - Xd * psi), 1e-7);
    }
  }
}

TEST(clara_Gates_Zd, AllTest) {
  for (idx D = 1; D < 10; ++D) {
    cmat Zd = gt.Zd(D);
    cplx oD = omega(D);
    for (idx i = 0; i < D; ++i) {
      ket psi = mket({i}, D);
      ket res = std::pow(oD, i) * psi;
      EXPECT_NEAR(0, norm(res - Zd * psi), 1e-7);
    }
  }
}
