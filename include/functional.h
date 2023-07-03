#ifndef FUNCTIONAL_H_
#define FUNCTIONAL_H_

#include <cmath>
#include <eigen3/Eigen/QR>
#include <stdexcept>

#include "internal.h"
#include "types.h"

namespace clara {
template <typename Derived>
Eigen::MatrixXcd funm(const Eigen::MatrixBase<Derived> &A, types::cplx (*f)(const types::cplx &)) {
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: funm matrix must be square");
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A.template cast<types::cplx>());
  Eigen::MatrixXcd evects = es.eigenvectors();
  Eigen::MatrixXcd evals = es.eigenvalues();
  for (int i = 0; i < evals.rows(); i++)
    evals(i) = (*f)(evals(i));

  Eigen::MatrixXcd evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

template <typename Derived>
Eigen::MatrixXcd expm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::exp);
}

// matrix square root
template <typename Derived>
Eigen::MatrixXcd sqrtm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::sqrt);
}
// matrix sin
template <typename Derived>
Eigen::MatrixXcd sinm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::sin);
}

// matrix cos
template <typename Derived>
Eigen::MatrixXcd cosm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::cos);
}

// matrix hyperbolic sin
template <typename Derived>
Eigen::MatrixXcd sinhm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::sinh);
}

// matrix hyperbolic cos
template <typename Derived>
Eigen::MatrixXcd coshm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::cosh);
}

// matrix absolute value
template <typename Derived>
Eigen::MatrixXcd absm(const Eigen::MatrixBase<Derived> &A) {
  return funm(adjoint(A) * A, [](const types::cplx &x) -> types::cplx { return std::sqrt(x); });
}

}  // namespace clara

#endif  // !FUNCTIONAL_H_
