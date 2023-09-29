#include <complex>
#include <cstdlib>
#include <numeric>
#include <vector>

#include "../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_apply, AllTests) {
  // Eigen::Matrix2cd state;
  //
  // state << std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0),
  //     std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 1.0);
  // Eigen::MatrixX2cd gate;
  // gate << std::complex<double>(0.0, 1.0), std::complex<double>(1.0, 0.0),
  //     std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 1.0);
  // std::vector<idx> subsys = {0};
  // std::vector<idx> dims = {2};
  // dyn_mat<std::complex<double>> result = apply(state, gate, subsys, dims);
  // EXPECT_EQ(result(0, 0), std::complex<double>(0.0, 1.0));
  ket psi = 1_ket;

  ket resultX = clara::apply(psi, gt.X, {0}, std::vector<idx>({2}));
  EXPECT_EQ(0_ket, resultX);
}

TEST(clara_apply, CustomGateOnQubit) {}

TEST(clara_applyCTRL, EmptyControl) {
  // 3 qubits
  std::vector<idx> dims{2, 2, 2, 2};
  idx D = prod(dims);

  // apply control
  std::vector<idx> ctrl{};
  // target
  std::vector<idx> target{1, 0, 3};

  // dimension of the target susbystem
  idx Dtarget = 1;
  for (idx i = 0; i < target.size(); ++i)
    Dtarget *= dims[target[i]];

  // random n qudit pure state
  ket psi = randket(D);

  // corresponding density matrix
  cmat rho = psi * adjoint(psi);
  cmat U = randU(Dtarget);

  // applyCTRL on pure state
  ket A = applyCTRL(psi, U, ctrl, target, dims);

  // applyCtrl on density matrix
  cmat B = applyCTRL(rho, U, ctrl, target, dims);

  cmat result_psi = A * adjoint(A);
  cmat result_rho = B;

  double res = norm(result_psi - result_rho);
  EXPECT_NEAR(0, res, 1e-7);
}
