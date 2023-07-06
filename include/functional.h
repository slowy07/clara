#ifndef FUNCTIONAL_H_
#define FUNCTIONAL_H_

#include <cmath>
#include <stdexcept>
#include "internal.h"
#include "types.h"

namespace clara {

// compute f(A), where (*f) is the function pointer
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

// apply f(A) component wise, where (*F) is function pointer
template <typename FunctionInputType, typename FunctionOutputType, typename MatrixInputType>
Eigen::Matrix<FunctionInputType, Eigen::Dynamic, Eigen::Dynamic> fun(
    const Eigen::MatrixBase<MatrixInputType> &A,
    FunctionOutputType (*f)(const FunctionInputType &)) {
  /*
    * the type of A is MatrixInputType
    * the function is the form FunctionOutputType f(const FunctionInputType &)
    * the output is an Eigen::Matrix of the type FunctionOutputType
    * 
    * the MatrixInputType is an general automatically deduced
    * if (*f) is not overloaded, then FunctionInputType and FunctionOutputType are
    * also automatically deduced
    */
  Eigen::Matrix<FunctionOutputType, Eigen::Dynamic, Eigen::Dynamic> result(A.rows(), A.cols());
  for (size_t i = 0; i < A.rows(); i++)
    for (size_t j = 0; j < A.cols(); j++)
      result(i, j) = (*f)(A(i,j));
  return result;
}

template <typename Derived>
Eigen::MatrixXcd expm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::exp);
}

// matrix logarithm
template<typename Derived>
Eigen::MatrixXcd logm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::log);
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

// matrix absolute value
template <typename Derived>
Eigen::MatrixXcd absm(const Eigen::MatrixBase<Derived> &A) {
  return funm(adjoint(A) * A, [](const types::cplx &x) -> types::cplx { return std::sqrt(x); });
}

}  // namespace clara

#endif  // !FUNCTIONAL_H_
