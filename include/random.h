#ifndef RANDOM_H_
#define RANDOM_H_

#include "constants.h"
#include "stat.h"
#include "types.h"

namespace clara {
// random double matrix with entries in uniform
inline Eigen::MatrixXd rand(size_t rows, size_t cols) {
  return Eigen::MatrixXd::Random(rows, cols);
}

// random double square matrix with entries in uniform
inline Eigen::MatrixXd rand(size_t rows) { return rand(rows, rows); }

// random double matrix with entries in normal
inline Eigen::MatrixXd randn(size_t rows, size_t cols) {
  stat::NormalDistribution nd;
  Eigen::MatrixXd A(rows, cols);

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      A(i, j) = nd.sample();
  return A;
}

// random suare matrix with entries in normal
inline Eigen::MatrixXd randn(size_t rows) { return randn(rows, rows); }

// random unitary matrix
inline Eigen::MatrixXcd rand_unitary(size_t size) {
  Eigen::MatrixXcd X(size, size);
  X.real() = 1. / sqrt(2) * randn(size);
  X.imag() = 1. / sqrt(2) * randn(size);
  Eigen::HouseholderQR<Eigen::MatrixXcd> qr(X);

  Eigen::MatrixXcd Q = qr.householderQ();
  /*
   * phase correction so that the resultation matrix
   * is uniformly distributed according to haar
   * measure
   */
  Eigen::VectorXcd phases = (rand(size, 1)).template cast<types::cplx>();
  for (size_t i = 0; i < (size_t)phases.rows(); i++)
    phases(i) = std::exp((types::cplx)(2 * ct::pi * ct::ii * phases(i)));
  Q = Q * phases.asDiagonal();
  return Q;
}

}  // namespace clara

#endif  // !RANDOM_H_
