#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <ios>
#include <iterator>
#include <memory>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "classFunction/exception.h"
#include "constants.h"
#include "internal/util.h"
#include "traits.h"
#include "types.h"

namespace clara {
// eigen function wrappers
template <typename Derived>
dyn_mat<typename Derived::Scalar> transpose(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::transpose()", Exception::Type::ZERO_SIZE);
  return rA.transpose();
}

// complex conjugate
template <typename Derived>
dyn_mat<typename Derived::Scalar> conjugate(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::conjugate()", Exception::Type::ZERO_SIZE);
  return rA.conjugate();
}

// ajdoint
template <typename Derived>
dyn_mat<typename Derived::Scalar> adjoint(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::ajdoint()", Exception::Type::ZERO_SIZE);
  return rA.adjoint();
}

// inverse
template <typename Derived>
dyn_mat<typename Derived::Scalar> inverse(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::inverse()", Exception::Type::ZERO_SIZE);
  return rA.inverse();
}

// trace
template <typename Derived>
typename Derived::Scalar trace(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::trace()", Exception::Type::ZERO_SIZE);
  return rA.trace();
}

// determinant
template <typename Derived>
typename Derived::Scalar det(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::det()", Exception::Type::ZERO_SIZE);
  return rA.determinant();
}

/**
 * logarithm of the determinant
 * usefull when the determinant overflows/underflows
 */
template <typename Derived>
typename Derived::Scalar logdet(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::logdet()", Exception::Type::ZERO_SIZE);

  // check square matrix
  if (!internal::check_square_mat(rA))
    throw Exception("clara::logdet()", Exception::Type::MATRIX_NOT_SQUARE);

  Eigen::PartialPivLU<dyn_mat<typename Derived::Scalar>> lu(rA);
  dyn_mat<typename Derived::Scalar> U = lu.matrixLU().template triangularView<Eigen::Upper>();
  typename Derived::Scalar result = std::log(U(0, 0));

  for (idx i = 1; i < static_cast<idx>(rA.rows()); i++)
    result += std::log(U(i, i));
  return result;
}

// element-wise sum of A
template <typename Derived>
typename Derived::Scalar sum(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::sum()", Exception::Type::ZERO_SIZE);
  return rA.sum();
}

// element-wise product of A
template <typename Derived>
typename Derived::Scalar prod(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::prod()", Exception::Type::ZERO_SIZE);
  return rA.prod();
}

// fronbenius form
template <typename Derived>
double norm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::norm()", Exception::Type::ZERO_SIZE);
  // convert matrix to complex then return its norm
  return (rA.template cast<cplx>()).norm();
}

// full eigen decomposition
template <typename Derived>
std::pair<dyn_col_vect<cplx>, cmat> eig(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::eig()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::eig()", Exception::Type::MATRIX_NOT_SQUARE);

  Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());
  return std::make_pair(es.eigenvalues(), es.eigenvectors());
}

// brief eigenvalues
template <typename Derived>
dyn_col_vect<cplx> evals(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::evals()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::evals()", Exception::Type::MATRIX_NOT_SQUARE);

  return eig(rA).first;
}

template <typename Derived>
cmat evects(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::evects()", Exception::Type::MATRIX_NOT_SQUARE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::evects()", Exception::Type::MATRIX_NOT_SQUARE);
  Eigen::ComplexEigenSolver<cmat> es(rA.template cas<cplx>());

  return eig(rA).second;
}

template <typename Derived>
std::pair<dyn_col_vect<double>, cmat> heig(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::heig()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::heig()", Exception::Type::MATRIX_NOT_SQUARE);
  Eigen::SelfAdjointEigenSolver<cmat> es(rA.template cast<cplx>());
  return std::make_pair(es.eigenvalues(), es.eigenvectors());
}

template <typename Derived>
dyn_col_vect<double> hevals(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_square_mat(rA))
    throw Exception("clara::hevals()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::hevals()", Exception::Type::MATRIX_NOT_SQUARE);
  return heig(rA).first;
}

// hermitian eigenvectors
template <typename Derived>
cmat hevects(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::hevects()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::hevects()", Exception::Type::MATRIX_NOT_SQUARE);
  return heig(rA).second;
}

// full singular value decomposition
template <typename Derived>
std::tuple<cmat, dyn_col_vect<double>, cmat> svd(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::svd()", Exception::Type::ZERO_SIZE);
  Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(
      rA, Eigen::DecompositionOptions::ComputeFullU | Eigen::DecompositionOptions::ComputeFullV);
  return std::make_tuple(sv.matrix(), sv.singularValues(), sv.matrixV());
}

/**
 * singular values
 */
template <typename Derived>
dyn_col_vect<double> svals(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::svals()", Exception::Type::ZERO_SIZE);
  Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(rA);
  return sv.singularValues();
}

// left singular vectors
template <typename Derived>
cmat svdU(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.dervied();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::svdU()", Exception::Type::ZERO_SIZE);
  Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(rA,
                                                         Eigen::DecompositionOptions::ComputeFullU);
  return sv.matrixU();
}

// right singular vectors
template <typename Derived>
cmat svdV(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::svdV()", Exception::Type::ZERO_SIZE);
  Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(rA,
                                                         Eigen::DecompositionOptions::ComputeFullV);
  return sv.matrix();
}

// functional calculus
template <typename Derived>
cmat funm(const Eigen::MatrixBase<Derived>& A, cplx (*f)(const cplx&)) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::funm()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::funm()", Exception::Type::MATRIX_NOT_SQUARE);
  Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());
  cmat evects = es.eigenvectors();
  cmat evals = es.eigenvalues();
  for (idx i = 0; i < static_cast<idx>(evals.rows()); ++i)
    evals(i) = (*f)(evals(i));
  cmat evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

// matrix square root
template <typename Derived>
cmat sqrtm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::sqrtm()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::sqrtm()", Exception::Type::MATRIX_NOT_SQUARE);
  return funm(rA, &std::sqrt);
}

// matrix absolute value
template <typename Derived>
cmat absm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::absm()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::absm()", Exception::Type::MATRIX_NOT_SQUARE);
  return sqrtm(adjoint(rA) * rA);
}

// matrix exponential
template <typename Derived>
cmat expm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::expm()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::expm()", Exception::Type::MATRIX_NOT_SQUARE);
  return fnum(rA, &std::exp);
}

// matrix logarithm
template <typename Derived>
cmat logm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::logm()", Exception::Type::ZERO_SIZE);
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw Exception("clara::logm()", Exception::Type::MATRIX_NOT_SQUARE);
  return funm(rA, std::log);
}

// matrix sin
template <typename Derived>
cmat sinm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::sinm()", Exception::Type::ZERO_SIZE);
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw Exception("clara::sinm()", Exception::Type::MATRIX_NOT_SQUARE);
  return funm(rA, &std::sin);
}

// matrix cos
template <typename Derived>
cmat cosm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::cosm()", Exception::Type::ZERO_SIZE);
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw Exception("clara::cosm()", Exception::Type::MATRIX_NOT_SQUARE);
  return funm(rA, &std::cos);
}

// matrix power
// clara::pow()
template <typename Derived>
cmat spectralpow(const Eigen::MatrixBase<Derived>& A, const cplx z) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::spectralpow()", Exception::Type::ZERO_SIZE);
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw Exception("clara::spectralpow()", Exception::Type::MATRIX_NOT_SQUARE);

  // define A^0 = Id
  if (real(z) == 0 && imag(z) == 0)
    return cmat::Identity(rA.rows(), rA.rows());

  Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());
  cmat evects = es.eigenvalues();
  cmat evals = es.eigenvalues();
  for (idx i = 0; i < static_cast<idx>(evals.rows()); ++i)
    evals(i) = std::pow(evals(i), z);
  cmat evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

// fast matrix power based on the SQUARE-AND-MULTIPLY algorithm
template <typename Derived>
dyn_mat<typename Derived::Scalar> powm(const Eigen::MatrixBase<Derived>& A, idx n) {
  // check zero-size
  if (!internal::check_nonzero_size(A))
    throw Exception("clara::powm()", Exception::Type::ZERO_SIZE);
  // check square matrix
  if (!internal::check_square_mat(A))
    throw Exception("clara::powm()", Exception::Type::MATRIX_NOT_SQUARE);
  if (n == 1)
    return A;
  dyn_mat<typename Derived::Scalar> result =
      dyn_mat<typename Derived::Scalar>::Identity(A.rows(), A.rows());
  if (n == 0)
    return result;

  dyn_mat<typename Derived::Scalar> cA = A.derived();

  // fast matrix power
  for (; n > 0; n /= 2) {
    if (n % 2)
      result = (result * cA).eval();
    cA = (cA * cA).eval();
  }
  return result;
}

// schatten matrix norm
template <typename Derived>
double schatten(const Eigen::MatrixBase<Derived>& A, double p) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::schatten()", Exception::Type::ZERO_SIZE);
  if (p < 1)
    throw Exception("clara::schatten()", Exception::Type::OUT_OF_RANGE);

  if (p == clara::inifinity)
    return svals(rA)(0);
  return std::pow(trace(powm(absm(rA), p)).real(), 1. / p);
}

template <typename OutputScalar, typename Derived>
dyn_mat<OutputScalar> cwise(const Eigen::MatrixBase<Derived>& A,
                            OutputScalar (*f)(const typename Derived::Scalar&)) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::cwise()", Exception::Type::ZERO_SIZE);
  dyn_mat<OutputScalar> result(rA.rows(), rA.cols());

#ifndef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // !WITH_OPENMP_
  for (idx j = 0; j < static_cast<idx>(rA.cols()); ++j)
    for (idx i = 0; i < static_cast<idx>(rA.rows()); ++i)
      result(i, j) = (*f)(rA(i, j));
  return result;
}

/**
 * kronecker product of multiple matrices, preserve return type
 * variadic template
 */
template <typename T>
dyn_mat<typename T::Scalar> kron(const T& head) {
  return head;
}

/**
 * kronecker product
 * return kronecker product of all input parameters,
 * evaluated from left to right, as a dynamic matrix
 * over the same scalar field as its arguments
 */
template <typename T, typename... Args>
dyn_mat<typename T::Scalar> kron(const T& head, const Args&... tail) {
  return internal::kron2(head, kron(tail...));
}

/**
 * kronecker product
 * return kronecker product of all elements in As,
 * evalueted from left to right, as a dynamic matrix over
 * the same scalar field as its arguments
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> kron(const std::vector<Derived>& As) {
  if (As.size() == 0)
    throw Exception("clara::kron()", Exception::Type::ZERO_SIZE);

  for (auto&& it : As)
    if (!internal::check_nonzero_size(it))
      throw Exception("clara::kron()", Exception::Type::ZERO_SIZE);

  dyn_mat<typename Derived::Scalar> result = As[0].derived();
  for (idx i = 1; i < As.size(); ++i) {
    result = kron(result, As[i]);
  }
  return result;
}

/**
 * kronecker prodcut of a list matrices, preserve, return type
 * educe the template parameters from intilizer_list
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> kron(const std::initializer_list<Derived>& As) {
  return kron(std::vector<Derived>(As));
}

/**
 * kronecker product
 * return kronecker product of A with itself n items as dynamic
 * matrix over the same scalar field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> kronpow(const Eigen::MatrixBase<Derived>& A, idx n) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::kronpow()", Exception::Type::ZERO_SIZE);
  if (n == 0)
    throw Exception("clara::kronpow()", Exception::Type::OUT_OF_RANGE);
  std::vector<dyn_mat<typename Derived::Scalar>> As(n, rA);
  return kron(As);
}

/**
 * direct sum of multiple matrices, preserve return type
 * variadic template
 */
template <typename T>
dyn_mat<typename T::Scalar> dirsum(const T& head) {
  return head;
}

/**
 * direct sum
 * return sum of all inputs parameters,
 * evaluatefrom left to right, as dynamic matrix
 * over the same scalar field as its arguments
 */
template <typename T, typename... Args>
dyn_mat<typename T::Scalar> dirsum(const T& head, const Args&... tail) {
  return internal::dirsum2(head, dirsum(tail...));
}

/**
 * direct sum
 * return direct sum of all elements in As,
 * evaluated from left to right, as a dynamic matrix
 * over the same scalar field as its arguments
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> dirsum(const std::vector<Derived>& As) {
  if (As.size() == 0)
    throw Exception("clara::dirsum()", Exception::Type::ZERO_SIZE);
  for (auto&& it : As)
    if (!internal::check_nonzero_size(it))
      throw Exception("clara::dirsum()", Exception::Type::ZERO_SIZE);
  idx total_rows = 0, total_cols = 0;
  for (idx i = 0; i < As.size(); ++i) {
    total_rows += static_cast<idx>(As[i].rows());
    total_cols += static_cast<idx>(As[i].cols());
  }
  dyn_mat<typename Derived::Scalar> result =
      dyn_mat<typename Derived::Scalar>::Zero(total_rows, total_cols);
  idx cur_row = 0, cur_col = 0;
  for (idx i = 0; i < As.size(); i++) {
    result.block(cur_row, cur_col, As[i].rows(), As[i].cols()) = As[i];
    cur_row += static_cast<idx>(As[i].rows());
    cur_col += static_cast<idx>(As[i].cols());
  }
  return result;
}

/**
 * direct sum of a list of matrices, preserve return type
 * deduce the template parameters from initializer_list
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> dirsum(const std::initializer_list<Derived>& As) {
  return dirsum(std::vector<Derived>(As));
}

/**
 * direct sum power
 * return direct sum of A with itself n times, as a dynamics
 * matrix over the same scalar field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> dirsumpow(const Eigen::MatrixBase<Derived>& A, idx n) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::dirsumpow()", Exception::Type::ZERO_SIZE);
  if (n == 0)
    throw Exception("clara::dirsumpow()", Exception::Type::OUT_OF_RANGE);
  std::vector<dyn_mat<typename Derived::Scalar>> As(n, rA);
  return dirsum(As);
}

/**
 * reshape
 * use column-major order when reshaping
 * return reshaped matrix with rows and cols columns,
 * as a dynamic matrix over the same scalar field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> reshape(const Eigen::MatrixBase<Derived>& A, idx rows, idx cols) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  idx Arows = static_cast<idx>(rA.rows());
  idx Acols = static_cast<idx>(rA.cols());

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::reshape()", Exception::Type::ZERO_SIZE);
  if (Arows * Acols != rows * cols)
    throw Exception("clara::reshape()", Exception::Type::DIMS_MISMATCH_MATRIX);
  return Eigen::Map<dyn_mat<typename Derived::Scalar>>(
      const_cast<typename Derived::Scalar*>(rA.data()), rows, cols);
}

// comulator
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> comm(const Eigen::MatrixBase<Derived1>& A,
                                        const Eigen::MatrixBase<Derived2>& B) {
  const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
  const dyn_mat<typename Derived2::Scalar>& rB = B.derived();
  if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
    throw Exception("clara::comm()", Exception::Type::TYPE_MISMATCH);
  if (!internal::check_nonzero_size(rA) || !internal::check_nonzero_size(rB))
    throw Exception("clara::comm()", Exception::Type::ZERO_SIZE);

  // check square matrices
  if (!internal::check_square_mat(rA) || !internal::check_square_mat(rB))
    throw Exception("clara::comm()", Exception::Type::MATRIX_NOT_SQUARE);

  // check equal dimension
  if (rA.rows() != rB.rows())
    throw Exception("clara::comm()", Exception::Type::DIMS_NOT_EQUAL);
  return rA * rB - rB * rA;
}

/**
 * anti cumulator
 * return anti-commutator \f$AB + BA\f$ as a dynamic matrix
 * over the same scalar field as
 */
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> anticomm(const Eigen::MatrixBase<Derived1>& A,
                                            const Eigen::MatrixBase<Derived2>& B) {
  const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
  const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

  if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
    throw Exception("clara::anticomm()", Exception::Type::TYPE_MISMATCH);

  if (!internal::check_nonzero_size(rA) || !internal::check_nonzero_size(rB))
    throw Exception("clara::anticomm()", Exception::Type::ZERO_SIZE);

  // check square matrix
  if (internal::check_square_mat(rA) || !internal::check_square_mat(rB))
    throw Exception("clara::anticomm()", Exception::Type::MATRIX_NOT_SQUARE);

  // check equal dimension
  if (rA.rows() != rB.rows())
    throw Exception("clara::anticomm()", Exception::Type::DIMS_NOT_EQUAL);
  return rA * rB + rB * rA;
}

template <typename Derived>
dyn_mat<typename Derived::Scalar> prj(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::prJ()", Exception::Type::ZERO_SIZE);

  // check column vector
  if (!internal::check_cvector(rA))
    throw Exception("clara::prj()", Exception::Type::MATRIX_NOT_CVECTOR);
  double normA = norm(rA);
  if (normA > eps)
    return rA * adjoint(rA) / (normA * normA);
  else
    return dyn_mat<typename Derived::Scalar>::Zero(rA.rows(), rA.rows());
}

// gram-shmidt orthogonalization
template <typename Derived>
dyn_mat<typename Derived::Scalar> grams(const std::vector<Derived>& As) {
  if (!internal::check_nonzero_size(As))
    throw Exception("clara::grams()", Exception::Type::ZERO_SIZE);

  for (auto&& it : As)
    if (!internal::check_nonzero_size(it))
      throw Exception("clara::grams()", Exception::Type::ZERO_SIZE);

  // check that As[0] is a column vector
  if (!internal::check_cvector(As[0]))
    throw Exception("clara::grams()", Exception::Type::MATRIX_NOT_CVECTOR);

  // check that all rest match the size of the first vector
  for (auto&& it : As)
    if (it.rows() != As[0].rows() || it.cols() != 1)
      throw Exception("clara::grams()", Exception::Type::DIMS_NOT_EQUAL);

  dyn_mat<typename Derived::Scalar> cut =
      dyn_mat<typename Derived::Scalar>::Identity(As[0].rows(), As[0].rows());
  dyn_mat<typename Derived::Scalar> vi = dyn_mat<typename Derived::Scalar>::Zero(As[0].rows(), 1);

  std::vector<dyn_mat<typename Derived::Scalar>> outvecs;
  // unf the first non-zero vector in the list
  idx pos = 0;
  for (pos = 0; pos < As.size(); ++pos) {
    if (norm(As[pos]) > eps) {
      outvecs.push_back(As[pos]);
      break;
    }
  }

  // start the process
  for (idx i = pos + 1; i < As.size(); ++i) {
    cut -= prj(outvecs[i - 1 - pos]);
    vi = cut * As[i];
    outvecs.push_back(vi);
  }

  dyn_mat<typename Derived::Scalar> result(As[0].rows(), outvecs.size());
  idx cnt = 0;
  for (auto&& it : outvecs) {
    double normA = norm(it);
    if (normA > eps) {
      result.col(cnt) = it / normA;
      cnt++;
    }
  }
  return result.block(0, 0, As[0].rows(), cnt);
}

/**
 * deduce template parameter from initializer_list
 * return gram-shmidt vectors of As as column of a dynamic matrix
 * over the same scala field as its arguments
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> grams(const std::initializer_list<Derived>& As) {
  return grams(std::vector<Derived>(As));
}

/**
 * gram-shmidt orthogonalization
 * return gram-shmidt vectors of columns of A, as columns
 * of a dynamic matrix over the same scalar field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> grams(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::grams()", Exception::Type::ZERO_SIZE);
  std::vector<dyn_mat<typename Derived::Scalar>> input;
  for (idx i = 0; i < static_cast<idx>(rA.cols()); ++i)
    input.push_back(rA.cols(i));
  return grams<dyn_mat<typename Derived::Scalar>>(input);
}

/**
 * non negative integer index to multi-index
 * rturn multi index of the same size as dims
 */
inline std::vector<idx> n2multiidx(idx n, const std::vector<idx>& dims) {
  if (!internal::check_dims(dims))
    throw Exception("clara::n2multiidx()", Exception::Type::DIMS_INVALID);
  if (n >= std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                           std::multiplies<idx>()))
    throw Exception("clara::n2multiidx()", Exception::Type::OUT_OF_RANGE);
  idx result[2 * maxn];
  internal::n2multiidx(n, dims.size(), dims.data(), result);
  return std::vector<idx>(result, result + dims.size());
}

/**
 * multi-index to non-negative integer index
 * return non-negative integer index
 */
inline idx multiidx2n(const std::vector<idx>& midx, const std::vector<idx>& dims) {
  if (!internal::check_dims(dims))
    throw Exception("clara::multiidx2n()", Exception::Type::DIMS_INVALID);
  for (idx i = 0; i < dims.size(); ++i)
    if (midx[i] >= dims[i])
      throw Exception("clara::multiidx2n()", Exception::Type::OUT_OF_RANGE);
  return internal::multiidx2n(midx.data(), dims.size(), dims.data());
}

/**
 * multi-partite qudit ket
 * construct the multi-partite qudit ket \f$|\mathrm{mask}\rangle\f$,
 * where mask is a std::vector of non-negative integers.
 * each element in mask has to be smaller than the corresponding element
 * in dims
 */
inline ket mket(const std::vector<idx>& mask, const std::vector<idx>& dims) {
  idx N = mask.size();
  idx D = std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                          std::multiplies<idx>());

  if (N == 0)
    throw Exception("clara::mket()", Exception::Type::ZERO_SIZE);
  // check valid dims
  if (!internal::check_dims(dims))
    throw Exception("clara::mket()", Exception::Type::DIMS_INVALID);
  // check mask and dims have the same size
  if (mask.size() != dims.size())
    throw Exception("clara::mket()", Exception::Type::SUBSYS_MISMATCH_DIMS);
  // check mask is a valid vector
  for (idx i = 0; i < N; ++i)
    if (mask[i] >= dims[i])
      throw Exception("clara::mket()", Exception::Type::SUBSYS_MISMATCH_DIMS);

  ket result = ket::Zero(D);
  idx pos = multiidx2n(mask, dims);
  result(pos) = 1;
  return result;
}

/**
 * multi-partite qudit ket
 * construct the multi-partite qudit ket \f$|\mathrm{mask}\rangle\f$,
 * all subsystem having equal dimension d
 * mask is a std::vector of non-negative integers, and each element in masuk to be strictly smaller
 * than d
 */
inline ket mket(const std::vector<idx>& mask, idx d = 2) {
  idx N = mask.size();
  idx D = static_cast<idx>(std::llround(std::pow(d, N)));

  // check zero-size
  if (N == 0)
    throw Exception("clara::mket()", Exception::Type::ZERO_SIZE);
  // check valid dims
  if (d == 0)
    throw Exception("clara::mket()", Exception::Type::DIMS_INVALID);

  for (idx i = 0; i < N; ++i)
    if (mask[i] >= d)
      throw Exception("clara::mket()", Exception::Type::SUBSYS_MISMATCH_DIMS);
  ket result = ket::Zero(D);
  std::vector<idx> dims(N, d);
  idx pos = multiidx2n(mask, dims);
  result(pos) = 1;

  return result;
}

/**
 * projector onto multi-partite qudit ket
 * construct the projector into the multi-partite qudit ket
 * \f$|\mathrm{mask}\rangle\f$,
 * where mask is a std::vector of non-negative integers
 * each element in mask to be smaller than the corresponding element in
 * dims
 */
inline cmat mprj(const std::vector<idx>& mask, const std::vector<idx>& dims) {
  idx N = mask.size();
  idx D = std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                          std::multiplies<idx>());

  // check zero-size
  if (N == 0)
    throw Exception("clara::mprj()", Exception::Type::ZERO_SIZE);
  // check valid dims
  if (!internal::check_dims(dims))
    throw Exception("clara::mprj()", Exception::Type::DIMS_INVALID);
  // check mask and dims have the same size
  if (mask.size() != dims.size())
    throw Exception("clara::mprj()", Exception::Type::SUBSYS_MISMATCH_DIMS);
  // check mask is valid vector
  for (idx i = 0; i < N; ++i)
    if (mask[i] >= dims[i])
      throw Exception("clara::mprj()", Exception::Type::SUBSYS_MISMATCH_DIMS);
  cmat result = cmat::Zero(D, D);
  idx pos = multiidx2n(mask, dims);
  result(pos, pos) = 1;
  return result;
}

/**
 * projector into multi-partite qudit ket
 * return projector into multi-partite qudit state vector
 * as complex dynamic matrix
 */
inline cmat mprj(const std::vector<idx>& mask, idx d = 2) {
  idx N = mask.size();
  idx D = static_cast<idx>(std::llround(std::pow(d, N)));

  // check zero size
  if (N == 0)
    throw Exception("clara::mprj()", Exception::Type::ZERO_SIZE);

  // check valid dims
  if (d == 0)
    throw Exception("clara::mprj()", Exception::Type::DIMS_INVALID);

  // check mask is a valid vector
  for (idx i = 0; i < N; ++i)
    if (mask[i] >= d)
      throw Exception("clara::mprj()", Exception::Type::SUBSYS_MISMATCH_DIMS);

  cmat result = cmat::Zero(D, D);
  std::vector<idx> dims(N, d);
  idx pos = multiidx2n(mask, dims);
  result(pos, pos) = 1;
  return result;
}

template <typename InputIterator>
std::vector<double> abssq(InputIterator first, InputIterator last) {
  std::vector<double> weights(std::distance(first, last));
  std::transform(first, last, std::begin(weights), [](cplx z) -> double { return std::norm(z); });
  return weights;
}

/**
 * computes the absolute values square of an STL-like container
 * real vector consitingof container's absolute values squared
 */
template <typename Container>
std::vector<double> abssq(const Container& c,
                          typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  /**
   * std::enable_if to SFINAE out Eigen expression
   * that will otherwise match, instead of matching
   * the overload
   */
  return abssq(std::begin(c), std::end(c));
}

template <typename Derived>
std::vector<double> abssq(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::abssq()", Exception::Type::ZERO_SIZE);
  return abssq(rA.data(), rA.data() + rA.size());
}

/**
 * element-wise sum of an STL-like range
 * returen element-wise of the range, as a scalar
 * over the same scalar field as the range
 */
template <typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type sum(InputIterator first,
                                                             InputIterator last) {
  using value_type = typename std::iterator_traits<InputIterator>::value_type;
  return std::accumulate(first, last, static_cast<value_type>(0));
}

/**
 * elemen-wise sum of the element of an STL-like container
 * as scalar over the same scalar field as the container
 */
template <typename Container>
typename Container::value_type sum(
    const Container& c, typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  return sum(std::begin(c), std::end(c));
}

/**
 * element-wise product of an STL-like range
 * return element-wise product of the range, as a
 * scalar over the the same scalar field as the range
 */
template <typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type prod(InputIterator first,
                                                              InputIterator last) {
  using value_type = typename std::iterator_traits<InputIterator>::value_type;
  return std::accumulate(first, last, static_cast<value_type>(1), std::multiplies<value_type>());
}

/**
 * element-wise product elements of an STL-like container
 * return element-wise product of the element of the container,
 * as a scalar over the same scalar field as the container
 */
template <typename Container>
typename Container::value_type prod(
    const Container& c, typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  return prod(std::begin(c), std::end(c));
}

/**
 * finds the pure state representation of a matrix proprtional to a projector into a pure state
 * return the unique non-zero eigenvectors of A, as dynamic column vector over the same scalar field
 * as A
 */
template <typename Derived>
dyn_col_vect<typename Derived::Scalar> rho2pure(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::rho2pure()", Exception::Type::ZERO_SIZE);
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw Exception("clara::rho2pure*()", Exception::Type::MATRIX_NOT_SQUARE);

  dyn_col_vect<double> tmp_evals = hevals(rA);
  cmat tmp_evects = hevects(rA);
  dyn_col_vect<typename Derived::Scalar> result =
      dyn_col_vect<typename Derived::Scalar>::Zero(rA.rows());
  for (idx k = 0; k < static_cast<idx>(rA.rows()); ++k) {
    if (std::abs(tmp_evals(k)) > eps) {
      result = tmp_evects.col(k);
      break;
    }
  }
  return result;
}

/**
 * construct the complement of susbsystem vector
 * return a complement of subsys with respect to the set
 * \f$\{0, 1, \ldots, N - 1\}\f$
 */
template <typename T>
std::vector<T> complement(std::vector<T> subsys, idx N) {
  if (N < subsys.size())
    throw Exception("clara::complement()", Exception::Type::OUT_OF_RANGE);
  std::vector<T> all(N);
  std::vector<T> subsys_bar(N - subsys.size());

  std::iota(std::begin(all), std::end(all), 0);
  std::sort(std::begin(subsys), std::end(subsys));
  std::set_difference(std::begin(all), std::end(all), std::begin(subsys), std::end(subsys),
                      std::begin(subsys_bar));
  return subsys_bar;
}

template <typename Derived>
std::vector<double> rho2bloch(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check qubit matrix
  if (!internal::check_qubit_matrix(rA))
    throw Exception("clara::rho2bloch()", Exception::Type::NOT_QUBIT_MATRIX);

  std::vector<double> result(3);
  cmat X(2, 2), Y(2, 2), Z(2, 2);
  X << 0, 1, 1, 0;
  Y << 0, -1_i, 1_i, 0;
  Z << 1, 0, 0, -1;
  result[0] = std::real(trace(rA * X));
  result[1] = std::real(trace(rA * Y));
  result[2] = std::real(trace(rA * Z));
  return result;
}

/**
 * compute the density matrix corresponding to the 3-dimensional real bloch vector r
 * and return qubit density matrix
 */
inline cmat bloch2rho(const std::vector<double>& r) {
  if (r.size() != 3)
    throw Exception("clara::bloch2rho", "r is not a 3-dimensional vector!");

  cmat X(2, 2), Y(2, 2), Z(2, 2), Id2(2, 2);
  X << 0, 1, 1, 0;
  Y << 0, -1_i, 1_i, 0;
  Z << 1, 0, 0, -1;
  Id2 << 1, 0, 0, 1;

  return (Id2 + r[0] * X + r[1] * Y + r[2] * Z) / 2.;
}

}  // namespace clara

#endif  // !FUNCTIONS_H_
