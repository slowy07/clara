#include <numeric>
#include <vector>
#include "gtest/gtest.h"
#include "../../include/clara.h"

using namespace clara;

TEST(clara_applyCTRL, NoEmptyControl) {
  std::vector<idx> dims{2, 2, 2, 2};
  idx D = prod(dims);
  
  std::vector<idx> ctrl{2, 0};
  std::vector<idx> target{1, 3};
  
  // dimension of the target subsystem
  idx Dtarget = 1;
  for (idx i = 0; i < target.size(); ++i)
    Dtarget *= dims[target[i]];

  // random n qudit pure state
  ket psi = randket(D);

  // the corresponding density matrix
  cmat rho = psi * adjoint(psi);
  // some random unitary on the target
  cmat U = randU(Dtarget);

  // applyCTRL on pure state
  ket A = applyCTRL(psi, U, ctrl, target, dims);

  // applyCTRL on density matrix
  cmat B = applyCTRL(rho, U, ctrl, target, dims);

  cmat result_psi = A * adjoint(A);
  cmat result_rho = B;
  
  double res = norm(result_psi - result_rho);
  EXPECT_NEAR(0, res, 1e-7);
}

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
