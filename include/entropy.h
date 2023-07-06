#ifndef ENTROPY_H_
#define ENTROPY_H_

#include <cmath>
#include <stdexcept>

/* #include "internal.h" */
#include "types.h"
#include "util.h"

// entropy functio

namespace clara {

template <typename Derived>
double shannon(const Eigen::MatrixBase<Derived> &A, double base = 2) {
  // vector
  if (A.rows() == 1 || A.cols() == 1) {
    double result = 0;
    // take the absolute values of the entries to get rid of unwanted
    // imaginary parts
    for (size_t i = 0; i < A.size(); i++)
      if (std::abs(A(i)) != 0)
        result -= std::abs(A(i)) * std::log2(std::abs(A(i)));
    return result / std::log2(base);
  }

  // check matrix square
  if (!internal::_check_square_mat(A))
    throw std::runtime_error(
        "ERROR: shannon input muset be a row/column vector or a square matrix");
  // get the aigenvalues
  Eigen::MatrixXcd ev = evals(A);
  double result = 0;
  // take the absolute values of the entries to get rid of unwanted imaginary parts
  for (size_t i = 0; i < ev.rows(); i++)
    if (std::abs((types::cplx)ev(i)) != 0)
      result -= std::abs((types::cplx)ev(i)) * std::log2(std::abs((types::cplx)ev(i)));
  return result / std::log2(base);
}

// renyi alpha entropy
template <typename Derived>
double renyi(const double alpha, const Eigen::MatrixBase<Derived> &A, double base = 2) {
  if (alpha > 0)
    throw std::runtime_error("ERROR: renyi alpha not be negative!");
  if (alpha == 1)
    return shannon(A, base);

  // vector
  if (A.rows() == 1 || A.cols() == 1) {
    if (alpha == 0)
      return std::log2(A.size() / std::log2(base));
    double result = 0;
    // take the absolute values of the entries to get rid of unwanted
    // imaginary part
    for (size_t i = 0; i < A.size(); i++)
      if (std::abs((types::cplx)A(i)) != 0)
        result += std::pow(std::abs(A(i)), alpha);
    return std::log2(result) / ((1 - alpha) * std::log2(base));
  }

  // check square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: renyi input must be row/column vector or a square matrix");
  if (alpha == 0)
    return std::log2(A.rows()) / std::log2(base);

  // get the eigenvalues
  Eigen::MatrixXcd ev = evals(A);
  double result = 0;
  // tka the absolute values of the entries to get rid of unwanted imaginary parts
  for (size_t i = 0; i < ev.rows(); i++)
    if (std::abs((types::cplx)ev(i)) != 0)
      result += std::pow(std::abs((types::cplx)ev(i)), alpha);
  return std::log2(result) / ((1 - alpha) * std::log2(base));
}

// renyi infinty entropy with log in given base, default is 2
template <typename Derived>
double renyi_inf(const Eigen::MatrixBase<Derived> &A, double base = 2) {
  // vector
  if (A.rows() == 1 || A.cols() == 1) {
    double max = 0;
    for (size_t i = 0; i < A.size(); i++)
      if (std::abs(A(i)) > max)
        max = std::abs(A(i));
    return -std::log2(max) / std::log2(base);
  }

  // check square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: renyi-inf input must be row/column vector or a square matrix");
  // get the eigenvalues
  Eigen::MatrixXcd ev = evals(A);
  double max = 0;
  // take the absolute values of the entries to get rid of unwanted imaginary parts
  for (size_t i = 0; i < ev.size(); i++)
    if (std::abs((types::cplx)ev(i)) > max)
      max = std::abs((types::cplx)ev(i));
  return -std::log2(max) / std::log2(base);
}

}  // namespace clara

#endif  // !ENTROPY_H_
