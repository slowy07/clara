#include <random>
#include <vector>

#include "../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_entanglement_negativity_test, AllTest) {
  cmat rho = randrho(3);
  cmat sigma = randrho(3);
  EXPECT_NEAR(0, clara::negativity(kron(rho, sigma), {3, 3}), 1e-7);
}

TEST(clara_schmidtprobs_qubits, random_2x2_product_state) {
  idx dimension = 2;
  cmat UA = randU(dimension);
  cmat UB = randU(dimension);
  ket psi = kron(UA, UB) * mket({0, 0});
  std::vector<double> result_vect = clara::schmidtprobs(psi);

  dyn_col_vect<double> result =
      Eigen::Map<dyn_col_vect<double>>(result_vect.data(), result_vect.size());
  dyn_col_vect<double> expected = dyn_col_vect<double>::Zero(dimension);
  expected(0) = 1;
  EXPECT_NEAR(0, norm(result - expected), 1e-7);
}

TEST(clara_schmidtprobs_qubits, random_3x3_product_state) {
  idx dimension = 3;
  cmat UA = randU(dimension);
  cmat UB = randU(dimension);
  ket psi = kron(UA, UB) * mket({1, 1}, dimension);
  std::vector<double> result_vect = clara::schmidtprobs(psi, dimension);
  dyn_col_vect<double> result =
      Eigen::Map<dyn_col_vect<double>>(result_vect.data(), result_vect.size());
  dyn_col_vect<double> expected = dyn_col_vect<double>::Zero(dimension);
  expected(0) = 1;
  EXPECT_NEAR(0, norm(result - expected), 1e-7);
}

TEST(clara_schmidtprobs_qubits, random_degenerate_1x1_product_state) {
  idx dA = 1, dB = 1, D = dA * dB, minD = std::min(dA, dB);
  cmat UA = randU(dA);
  cmat UB = randU(dB);
  ket psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
  dyn_col_vect<double> result = clara::schmidtcoeffs(psi, {dA, dB});
  dyn_col_vect<double> expected(minD);
  expected << 1;
  EXPECT_NEAR(0, norm(result - expected), 1e-7);
}

TEST(clara_schmidtprobs_qubits, random_degenerated_3x1_product_state) {
  idx dA = 3, dB = 1, D = dA * dB, minD = std::min(dA, dB);
  cmat UA = randU(dA);
  cmat UB = randU(dB);
  ket psi = kron(UA, UB) * mket({0, 0}, {dA, dB});
  dyn_col_vect<double> result = clara::schmidtcoeffs(psi, {dA, dB});
  dyn_col_vect<double> expected = dyn_col_vect<double>(minD);
  expected << 1;
  EXPECT_NEAR(0, norm(result - expected), 1e-7);
}

TEST(clara_schmidtprobs_qubits, random_2x2_state_with_fixed_schmidt_coeff) {
  idx dimension = 2;
  double c0 = 0.8, c1 = 0.6;
  cmat UA = randU(dimension);
  cmat UB = randU(dimension);
  ket psi = kron(UA, UB) * (c0 * st.zero(2) + c1 * st.one(2));
  std::vector<double> result_vect = clara::schmidtprobs(psi);
  dyn_col_vect<double> result =
      Eigen::Map<dyn_col_vect<double>>(result_vect.data(), result_vect.size());
  dyn_col_vect<double> expected = dyn_col_vect<double>::Zero(dimension);
  expected << c0 * c0, c1 * c1;
  EXPECT_NEAR(0, norm(result - expected), 1e-7);
}

TEST(clara_schmidtprobs_qubits, random_3x3_state_with_fixed_schmidt_coeff) {
  idx dimension = 3;
  double c0 = 0.8, c1 = 0.5;
  double c2 = std::sqrt(1 - c0 * c0 - c1 * c1);
  cmat UA = randU(dimension);
  cmat UB = randU(dimension);
  ket psi = kron(UA, UB) * (c0 * mket({0, 0}, dimension) + c1 * mket({1, 1}, dimension) +
                            c2 * mket({2, 2}, dimension));
  std::vector<double> result_vect = clara::schmidtprobs(psi, dimension);
  dyn_col_vect<double> result =
      Eigen::Map<dyn_col_vect<double>>(result_vect.data(), result_vect.size());
  dyn_col_vect<double> expected = dyn_col_vect<double>::Zero(dimension);
  expected << c0 * c0, c1 * c1, c2 * c2;
  EXPECT_NEAR(0, norm(result - expected), 1e-7);
}
