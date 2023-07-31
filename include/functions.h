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
#include <string>
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
/**
 * @brief transpose the input matrix
 * transposing a matrix swaps its rows and colsumns, effectively turning rows into columns
 * and columns into rows
 *
 * @param A the input matrix to be transposed
 * @return then transposed matrix of 'A'
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 2);
 * inputMatrix << 1, 2,
 *                3, 4,
 *                5, 6;
 *
 * // transpose the input matrix
 * Eigen::MatrixXd transposedMatrix = tranpose(inputMatrix);
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> transpose(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::transpose()");
  return rA.transpose();
}

/**
 * @brief compute the complex conjugate of the input matrix
 * the complex conjugate of a matrix is obtained by taking the complex conjugate
 * of each element in the matrix. for realmatrices, the complex conjugate operation
 * has no effect, as the imaginary part of each element is zero
 *
 * @param A the input matrix for which the complex conjugate is computed
 * @return the complex conjugate of the input matrix 'A'
 *
 * @exception exception::ZeroSize thrown if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXcd inputMatrix(2,2);
 * inputMatrix << std::complex<double>(1, 2), std::complex<double>(3,4), std::complex<double>(5, 6);
 *
 * // conjugate the complex conjugate of the input matrix
 * Eigen::MatrixXcd conjugateMatrix = conjugateMatrix(inputMatrix);
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> conjugate(const Eigen::MatrixBase<Derived>& A) {
  // check if the input matrix has a non-zero size
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::conjugate()");
  // return the complex conjutate of the matrix
  return rA.conjugate();
}

/**
 * @brief compute the adjoint (conjugate transpose) of the input matrix
 * the adjoint of matrix is obtained by taking the complex conjugate of each
 * element in the matrix and transposing it. for real matrices, the adjoint
 * operation is equivalent to the transpose operation
 *
 * @param A the input matrix for which the djoint is computed
 * @return the adjoint (conjugate transpose) of the input matrix 'A'
 *
 * @exception exception::ZeroSize thrown if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXcd inputMatrix(2, 2);
 * input matrix << std::complex<double(1, 2), std::complex(3, 4);
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> adjoint(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::adjoint()");
  return rA.adjoint();
}

/**
 * @brief compute the inverse of the input matrix
 * the inverse of a square matrix 'A' is another square matrix 'B' and that the
 * product of 'A' and 'B' is Identity matrix. in mathematrical terms, if 'A' is
 * asquare matrix and 'B' is its inverse, then 'A * B = B * A = I', where 'I' is
 * the identity matrix
 *
 * @param A the input matrix fro which the inverse is computed
 * @return the inverse of the input matrix 'A'
 *
 * @exception exception::ZeroSize if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(2, 2);
 * inputMatrix << 2, 1,
 *                4, 3;
 *
 * // compute the inverse of the input matrix
 * Eigen::MatrixXd inverseMatrix = invverse(inputMatrix);
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> inverse(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::inverse()");
  return rA.inverse();
}

/**
 * @brief compute the trace of the input matrix
 * the trace of square matrix is the sum of the elements on its main as diagonal.
 * in mathematical terms, if 'A' is a square matrix, then the trace 'Tr(A)'is
 * defined as the sum 'A(i, i)' for all 'i', where 'A(i, i)' is the element of
 * matrix 'A' at the i-th row and i-th column
 *
 * @param A the input matrix for which the trace is computed
 * @param the trace of the input matrix 'A'
 *
 * @exception exception::ZeroSize thrown if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(2, 2);
 * inputMatrix << 2, 1,
 *                4, 3;
 *
 * // compute the inverse of the input matrix
 * Eigen::MatrixXd inverseMatrix = inverse(inputMatrix);
 */
template <typename Derived>
typename Derived::Scalar trace(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::trace()");
  return rA.trace();
}

/**
 * @brief compute the determinant of the input matrix
 * the determinant of the square matrix is a scalar value that represent the
 * signed volume of the paralleliped spanned by the column vectors of the matrix
 *
 * @param A input matrix for wwhich the determinanat is computed
 * @return the determinant of the input matrix 'A'
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                4, 5, 6,
 *                7, 8, 9;
 *
 * // calculate the determinant of the input matrix
 * double determinantValue = det(inputMatrix);
 */
template <typename Derived>
typename Derived::Scalar det(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::det()");
  return rA.determinant();
}

/**
 * @brief compute the determinant of the input matrix
 *
 * the determinant of a square matrix is a scalar value that represent
 * the signed volume of the paralleliped spanned by the column vectors of the matrix
 *
 * @param A the input matrix for which the determinant is computed
 * @return the determinant of the input Matrix 'A'
 *
 * @exception exception::ZeroSize thrown if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                4, 5, 6,
 *                7, 8, 9;
 *
 * // compute the determinant of the input matrix
 * double determinantValue = det(inputMatrix);
 */
template <typename Derived>
typename Derived::Scalar logdet(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::logdet()");

  // check square matrix
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::logdet()");

  Eigen::PartialPivLU<dyn_mat<typename Derived::Scalar>> lu(rA);
  dyn_mat<typename Derived::Scalar> U = lu.matrixLU().template triangularView<Eigen::Upper>();
  typename Derived::Scalar result = std::log(U(0, 0));

  for (idx i = 1; i < static_cast<idx>(rA.rows()); i++)
    result += std::log(U(i, i));
  return result;
}

/**
 * @brief compute the element-wise sum of the input matrix
 *
 * @param A the input matrix for which the element-wise sum is computed
 * @return the sum of all elements in the input matrix 'A'
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                4, 5, 6,
 *                7, 8, 9;
 * // compute the element-wise sum of the input matrix
 * dpuble sumValue = sum(inputMatrix);
 */
template <typename Derived>
typename Derived::Scalar sum(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::sum()");
  return rA.sum();
}

/**
 * @brief compute the element-wise product of the input matrix
 *
 * @param A the input matrix which the element-wise product product is computed
 * @param the input the product of all elements in the input matrix 'A'
 *
 * @exception exception::ZeroSize thrown if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3
 *                4, 5, 6,
 *                7, 8, 9;
 * // compute the element-wise the product of the input matrix
 * double productValue= prod(inputMatrix);
 */
template <typename Derived>
typename Derived::Scalar prod(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::prod()");
  return rA.prod();
}

/**
 * @brief compute the fronbenius norm of the input matrix
 *
 * @param A the input matrix for which the fronbenius norm is computed
 * @return the fronbenius norm of the input matrix 'A'
 *
 * @exception exception::ZeroSize thrown if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3
 *                4, 5, 6,
 *                7, 8, 9;
 *
 * // compute the fronbenius norm of the input matrix
 * double normValue = nrom(inputMatrix);
 */
template <typename Derived>
double norm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::norm()");
  // convert matrix to complex then return its norm
  return (rA.template cast<cplx>()).norm();
}

/**
 * @brief compute the full eigen decomposition of the input square matrix
 *
 * @param A the input square matrix for which the eigen decomposition is computed
 * @param A std::pair containing eigenvalues (as a dynamic column vector)
 *
 * @exception exception::ZeroSize thrown if input matrix 'A' has zero size
 * @exception exception::MatrixNotSquare thrown if the input matrix 'A' is not square
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                4, 5, 6,
 *                7, 8, 9;
 * // compute the full eigen decomposition of input matrix
 * auto eigenDecomp = eig(inputMatrix);
 * Eigen::vectorXd eigenvalues = eigenDecomp.first.real();
 * Eigen::MatrixXd eigenvectors = eigenDecomp.second;
 */
template <typename Derived>
std::pair<dyn_col_vect<cplx>, cmat> eig(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::eig()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::eig()");

  Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());
  return std::make_pair(es.eigenvalues(), es.eigenvectors());
}

/**
 * @brief compute the eigenvalues of the input square matrix
 *
 * @param A the input square matrix for which the eigenvalues are computed
 * @param A dynamic column vector containing the eigenvalues
 *
 * @exception exception::ZeroSize thrown if the input matrix 'A' has zero size
 * @exeption exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                4, 5, 6,
 *                7, 8, 9;
 * // compute the eigenvalues of the input matrix
 * Eigen::vectorXd eigenvalues = evals(inputMatrix).real();
 * std::cout << "eigen value " << eigenvalues << std::endl;
 */
template <typename Derived>
dyn_col_vect<cplx> evals(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::evals()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::evals()");

  return eig(rA).first;
}

/**
 * @brief compute the eigenvectors of the input square matrix
 *
 * @param A the input matrix for which the eigenvectors are computed
 * @param A complex matrix containing the eigenvectors as columns
 *
 * @exception exception::ZeroSize thrown if the input matrix 'A' has zero size
 * @exception exception::MatrixNotSquare thrown if the input matrix 'A' is not square
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3
 *                4, 5, 6,
 *                7, 8, 9;
 * // compute the eigenvectors of the input matrix
 * Eigen::MatrixXcd eigenvectors = evects(inputMatrix)
 * std::cout << "eigenvectors:\n" << eigenvectors << std::endl;
 */
template <typename Derived>
cmat evects(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::evects()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::evects()");

  Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());

  return es.eigenvectors();
}

/**
 * @brief compute the eigenvalues and eigenvectors of the input hermitian matrix
 *
 * this function computes the eigenvalues and eigenvectors of the input hermitian marix 'A'
 * the function ensure that input matrix is square and has a non-zero before performing computation
 * the eigenvalues are returned as real vector, and the eigenvectors are returned as a complex
 * matrix each eigenvectors as columns
 *
 * @tparam derived the derived type of the input matrix 'A'
 * @param A the input hermitian matrix for which the eigenvalues and eigenvectors are computed
 * @return A pair of a real vector containing the eigenvalues and a complex matrix containing the
 *            eigenvectors as colsumns
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                2, 5, 6,
 *                3, 6, 9;
 * // compute the eigenvalues and eigenvectors of the input hermitian matrix
 * auto result = heig(inputMatrix);
 * Eigen::vectorXd eigenvalues = result.first;
 * Eigen::MatrixXcd eigenvectors = result.second;
 *
 * std::cout << "eigenvalues: \n" << eigenvalues << std::endl;
 * std::cout << "eigenvectors: \n" << eigenvectors << std::endl;
 */
template <typename Derived>
std::pair<dyn_col_vect<double>, cmat> heig(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::heig()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::heig()");
  // compute the eigenvalues and eigenvectors using 'SelfAdjointEigenSolver'
  // and return the result as a pair of real vector and complex matrix
  Eigen::SelfAdjointEigenSolver<cmat> es(rA.template cast<cplx>());
  return std::make_pair(es.eigenvalues(), es.eigenvectors());
}

/**
 * @brief compute the eigenvalues of the input hermitian matrix
 *
 * this function computes the eigenvalues of the input matrix 'A'
 * this function ensure that the input matrix is square and has non-zero before performing
 * computation the eigenvalues are returned as a real vector
 *
 * @tparam derived the derived type of the input matrix 'A'
 * @param A the input hermitian matrix for which for which the eigenvalues are computed
 * @retrun a real vector containing the eigenvalues
 *
 * @exception exception::ZeroSize if the input matrix 'A' zero size
 * @exception exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                2, 5, 6,
 *                3, 6, 9;
 * // compute the eigenvalues of the input hermitian matrix
 * Eigen::vectorXd eigenvalues = hevals(inputMatrix);
 * std::cout << "eigenvalues:\n" << eigenvalues << std::endl;
 */
template <typename Derived>
dyn_col_vect<double> hevals(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_square_mat(rA))
    throw exception::ZeroSize("clara::hevals()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::hevals()");
  // compute the eigenvalues using the 'heigh'
  return heig(rA).first;
}

/**
 * @brief compute the eigenvectors of the input hermitian matrix
 * this function computes the eigenvectors of the input hermitian matrix 'A'
 * the function ensures that input matrix is square and has a non-zero before performing
 * the computation. the eigenvectors returned as a complex matrix
 *
 * @tparam derived the derived type of the input matrix 'A'
 * @param A the input hermitian matrix for which the eigenvectors are computed
 * @return A complex matrix containing the eigenvectors as columns
 *
 * @exception exception::ZeroSize if the input matrix 'A' has zero size
 * @exception exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                2, 5, 6,
 *                3, 6, 9;
 * // compute the eigenvectors of the input hermitian matrix
 * Eigen::MatrixXcd eigenvectors = hevects(inputMatrix);
 * std::vector << "eigenvectors:\n" << eigenvectors << std::endl;
 */
template <typename Derived>
cmat hevects(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::hevects()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::hevects()");
  return heig(rA).second;
}

/**
 * @brief compute the Singular value decomposition (SVD) of the input matrix
 * this function computes the singluar value decomposition of the input matrix 'A'
 * the function ensure that input matrix has non-zero size before performing the computation
 * the SVD is returned as three matrices 'U', 'S', and 'V, such that 'A = U * S * V.conjugate()'
 *
 * @tparam Derived derived type of the input matrix 'A'
 * @param A the input matrix for which the SVD is computed
 * @return A tuple containing three matrices: 'U', 'S', and 'V', representing the SVD of 'A'
 *
 * @exception exception::ZeroSize if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 2);
 * inputMatrix << 1, 2
 *                3, 4,
 *                5, 6;
 * // compute the singular value decomposition (SVD) of the input matrix
 * auto [U, S, V] = svd(inputMatrix);
 *
 * std::cout <<"matrix U:\n" << U << std::endl;
 * std::cout <<"singular values:\n" << S << std::endl;
 * std::cout << "matrix V:\n" << V << std::endl;
 */
template <typename Derived>
std::tuple<cmat, dyn_col_vect<double>, cmat> svd(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::svd()");
  // compute the singular decomposition (SVD) using jacobi SVD and return U, S, and V matrices as
  // tuple
  Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(
      rA, Eigen::DecompositionOptions::ComputeFullU | Eigen::DecompositionOptions::ComputeFullV);
  return std::make_tuple(sv.matrix(), sv.singularValues(), sv.matrixV());
}

/**
 * @brief compute the singular values of the inpute matrix
 * this function computes the singular values of the input matrix 'A'
 * the function ensure that the input matrix has non-zero size before performing the computation
 * the singular value are returned as a column vector
 *
 * @tparam derived the derived type of the input matrix 'A'
 * @param A the input matrix for which the singluar values are computed
 * @return a column vector containing the singular values of 'A'
 *
 * @exception exception::ZeroSize if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                4, 5, 6,
 *                7, 8, 9;
 * // compute the singluar values of the input matrix
 * auto singularValues = svals(inputMatrix);
 * std::cout << "singular values:\n" << singularValues << std::endl;
 */
template <typename Derived>
dyn_col_vect<double> svals(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::svals()");
  // compute the singular using jacobiSVD and return as a column vector
  Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(rA);
  return sv.singularValues();
}

/**
 * @brief compute the left singular vectors of the input matrix
 * this function computes the left singular vectors of the input matrix 'A'
 * the function ensure that the input matrix has a non-zero size before performing the computation.
 * the left singular vectors are returned as a complex matrix
 *
 * @tparam derived derived type of the input matrix 'A'
 * @param A input matrix for which the left singular vectors are computed
 * @return A complex matrix containing the left singular vectors of 'A'
 *
 * @exception exception::ZeroSize if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                4, 5, 6,
 *                7, 8, 9;
 *
 * // compute the singular vectors of the input matrix
 * auto leftSingularVectors = svdU(inputMatrix);
 * std::cout << "left singular vectors:\n" << leftSingularVectors << std::endl;
 */
template <typename Derived>
cmat svdU(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.dervied();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::svdU()");
  Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(rA,
                                                         Eigen::DecompositionOptions::ComputeFullU);
  return sv.matrixU();
}

/**
 * @brief compute the right singular vectors of the input matrix
 * this function computes the right singular vectors of the input matrix 'A'.
 * the function ensure that the input matrix has a non-zero size before performing the computation
 * the right singular vectors are returned as a complex matrix
 *
 * @tparam Derived the derived type of the input matrix 'A'
 * @param A the input matrix for which the right singular vectors are computed
 * @return A complex matrix containing the right singular vectors of 'A'
 *
 * @exception exception::ZeroSize if the input matrix 'A' has zero size
 *
 * @example
 * Eigen::MatrixXd inputMatrix(3, 3);
 * inputMatrix << 1, 2, 3,
 *                4, 5, 6,
 *                7, 8, 9;
 * // compute the right singular vectors of the input matrix
 * auto rightSingularVectors = svdV(inputMatrix);
 * std::cout << "right singular vectors:\n" << rightSingularVectors << std::endl;
 */

template <typename Derived>
cmat svdV(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::svdV()");
  Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(rA,
                                                         Eigen::DecompositionOptions::ComputeFullV);
  return sv.matrix();
}

/**
 * @brief functional calculus of a complex matrix
 * this function applies a given complex-valued function 'f' element-wise to the eigenvalues of the
 * input matrix 'A'. the function then recronstruct the matrix using eigenvectors and the
 * transformed eigenvalues. the function ensure that the input matrix is square and has non-zero
 * size before performing the computation.
 *
 * @tparam Derived the derived type of the input matrix 'A'
 * @param A the input matrix to which the functional calculus is applied
 * @param f A pointer to a complex-valued function that is applied element-wise to the eigen values
 * of 'A'
 * @return A complex matrix obtained by applying the functional calculus to A
 *
 * @exception exception::ZeroSize if the input matrix 'A' has zero-size
 * @exception exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::Matrix2cd inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 * // defined a complex-valued function to be applied element-wise to the eigenvalues
 * auto func = [](const cplx& x) {return std::exp(x); };
 * // apply the functional calculus to the input matrix
 * auto resultMatrix = funm(inputMatrix, func);
 * std::cout << "result matrix:\n" << resultMatrix << std::endl;
 */
template <typename Derived>
cmat funm(const Eigen::MatrixBase<Derived>& A, cplx (*f)(const cplx&)) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::funm()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::funm()");
  Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());
  cmat evects = es.eigenvectors();
  cmat evals = es.eigenvalues();
  for (idx i = 0; i < static_cast<idx>(evals.rows()); ++i)
    evals(i) = (*f)(evals(i));
  cmat evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

/**
 * @brief matrix square root of a complex matrix
 * this function computes the square root of a complex matrix 'A' applying the functional calculus
 * it uses the standard library function 'std::sqrt' as the complex-valued function to compute the
 * square root element-wise. the function ensure that the input matrix square and has a non-zero
 * before performing the computation
 *
 * @tparam Dervied the derived type of the input matrix 'A'
 * @param A the input matrix for which the square root is computed
 * @return the square root of the input matrix 'A'
 *
 * @throws exception::ZeroSize if the input matrix 'A' has zero-size
 * @throws exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::Matrix2cd inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 * // compute the square root of the input matrix
 * auto sqrtMatrix = sqrtm(inputMatrix);
 * std::cout << "square root matrix:\n" << sqrtMatrix << std::endl;
 */
template <typename Derived>
cmat sqrtm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::sqrtm()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::sqrtm()");
  // compute the square root using functional calculus with std::sqrt
  // as the complex-valued function
  return funm(rA, &std::sqrt);
}

/**
 * @brief absolute value of a complex matrix
 * this function computes the absolute value of a complex matrix 'A' using the matrix square root.
 * it first computes the hermitian product of the adjoint of 'A' with 'A' and then takes the square
 * root of the result. the function ensure that the input mtarix is square and has a non-zero size
 * before performing the computation.
 *
 * @tparam Derived the derived type of the input matrix 'A'
 * @param A the input matrix for which the absolute value computed
 * @return the absolute value of the input matrix 'A'
 *
 * @throws exception::ZeroSize if the input matrix 'A' has zero-size
 * @throws exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::Matrix2cd inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 *
 * // compute the absolute value the input matrix
 * auto absMatrix = absm(inputMatrix);
 * std::cout << "absolute value matrix:\n" << absMatrix << std::endl;
 */
template <typename Derived>
cmat absm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::absm()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::absm()");
  // compute the absolute value using matrix square root of adjoint(A) * A;
  return sqrtm(adjoint(rA) * rA);
}

/**
 * @brief matrix exponential
 * this function computes the matrix exponential of square matrix 'A' using the functional
 * calculus method. it applies the exponential function to the matrix 'A', which is computed
 * using its eigenvalues and eigenvectors. the function ensures that the input matrix is square
 * and has non-zero size before performing the computation
 *
 * @tparam Derived the derived type of the input matrix 'A'
 * @param A the input matrix for which the matrix exponential is computed
 * @return the matrix exponential of the input matrix 'A'
 *
 * @throws exception::ZeroSize if the input matrix 'A' has zero-size
 * @throws exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::Matrix2cd inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 * // compute the matrix exponential of the input matrix
 * auto expMatrix = expm(inputMatrix);
 * std::cout << "matrix exponential:\n" << expMatrix << std::endl;
 */
template <typename Derived>
cmat expm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::expm()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::expm()");
  // compute the matrix exponential using the functional calculus method
  // exp function
  return fnum(rA, &std::exp);
}

/**
 * @brief matrix logarithm
 * this function compute the matrix logarithm of a square matrix 'A' using functional calculus
 * method it applies the natural logarithm function to the matrix 'A' which is computed using its
 * eigenvalues and eigenvectors the functions ensures that the input matrix square and has non-zero
 * size before performing the computation
 *
 * @tparam Derived the derived type of the input matrix 'A'
 * @param A the input matrix for which the matrix logarithm is computed
 * @return the matrix logarithm of the input matrix 'A'
 *
 * @throws exception::ZeroSize if the input matrix 'A'
 * @throws exception::MatrixNotSquare if the input matrix 'A' not square
 *
 * @example
 * Eigen::Matrix2cd inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 *
 * // compute the matrix logarithm of the input matrix
 * auto logMatrix = logm(inputMatrix);
 * std::cout << "matrix logarithm:\n" << logMatrix << std::endl;
 */
template <typename Derived>
cmat logm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::logm()");
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::logm()");
  // compute the matrix logarithm using functional calculus method natural
  // logarithm function
  return funm(rA, std::log);
}

/**
 * @brief matrix sine function
 * this function computes the matrix 'A' using the functional calculus method.
 * it applies the sine function to the matrix 'A' which is computed using its eigenvalues
 * and eigenvectors. the function ensures that the input matrix is square and has a non-zero size
 * before performing the computation
 *
 * @tparam Derived the derived type of the input matrix 'A'
 * @paran A the input matrix for the matrix sine is computed
 * @return matrix sine of the input matrix 'A'
 *
 * @throws exception::ZeroSize if the input matrix 'A' has zero size.
 * @throws exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::Matrix2cd inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 *
 * // compute the matrix sine of the input matrix
 * auto sinMatrix = sinm(inputMatrix);
 * std::cout << "matrix sine:\n" << sinMatrix << std::endl;
 */
template <typename Derived>
cmat sinm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::sinm()");
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::sinm()");
  return funm(rA, &std::sin);
}

/**
 * @brief matrix cosine function
 * this function computes the matrix cosine of a square matrix 'A' using the functional calculus
 * method it applies the cosine function to the matrix 'A' which is computed using its eigenvalues
 * and eigenvectors the function ensures that the input matrix is square and has a non-zero size
 * before performing the computation
 *
 * @tparam Derived the derived of the input matrix 'A'
 * @param the input matrix for which the matrix cosine is computed
 * @return the matri cosine of input matrix 'A'
 *
 * @throws exception::ZeroSize if the input matrix 'A' has zero size
 * @throws exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::Matrix2cd inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 *
 * auto cosMatrix = cosm(inputMatrix);
 * std::cout << "Matrix cosinue:\n" << cosMatrix << std::endl;
 */
template <typename Derived>
cmat cosm(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::cosm()");
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::cosm()");
  // compute the matrix cosinye using the functional calculus method with cosinue functions
  return funm(rA, &std::cos);
}

/**
 * @brief power of using spectral decomposition
 * this function computes the matrix power of a square matrix 'A' using the spectral decomposition
 * method. it computes A^z, where 'z' is complex number. by applying power to the eigenvalues and
 * recronstructing the matrix using the eigenvectors the functions ensures that input matrix is
 * square and has a non-zero size before performing the computation the function also handles the
 * spacial case when z is equal to 0, in which case result is the identity matrix
 *
 * @tparam Derived the derived tpe of the input matrix A
 * @param A the input matrix for which the matrix power is computed
 * @param z the complex number representing the power
 * @return the matrix power A^z of the input matrix 'A'
 *
 * @throws exception::ZeroSize if the input matrix 'A' has zero size
 * @throws exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::Matrix2cd inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 *
 * std::complex<double> power(2.0, -1.0)
 * // compute the matrix power A^z of the input matrix
 * auto result = spectralpow(inputMatrix, power);
 * std::cout << "matrix pow:\n" << result << std::endl;
 */
template <typename Derived>
cmat spectralpow(const Eigen::MatrixBase<Derived>& A, const cplx z) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::spectralpow()");
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::spectralpow()");

  // define A^0 = Id
  if (real(z) == 0 && imag(z) == 0)
    return cmat::Identity(rA.rows(), rA.rows());

  // compute the matrix power using spectral decomposition
  Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());
  cmat evects = es.eigenvalues();
  cmat evals = es.eigenvalues();
  for (idx i = 0; i < static_cast<idx>(evals.rows()); ++i)
    evals(i) = std::pow(evals(i), z);
  cmat evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

/**
 * @brief fast matrix power based on SQUARE-AND-MULTIPLY algorithm
 * this function computes the matrix power A^n of a square matrix 'A' using SQUARE-AND-MULTIPLY
 * algorithm. the function effectly computes the power of the matrix by repeatedly squaring 'A'
 * and multiplying the result.
 * the function ensures that the input matrix square and has non-zero size before performing the
 * computation
 *
 * @tparam Dervied the derived type of the input matrix 'A'
 * @param A the input matrix for which the matrix power is computed
 * @param n the exponent to which the matrix 'A' is raised
 * @return the matrix power A^n the input matrix 'A'
 *
 * @throws exception::ZeroSize if the input matrix 'A' has zero size
 * @throws exception::MatrixNotSquare if the input matrix 'A' is not square
 *
 * @example
 * Eigen::Matrix2d inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 * int power = 3;
 *
 * // compute the matrix power A^n of the input matrix
 * auto result = powm(inputMatrix, power);
 * std::cout << "matrix power:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> powm(const Eigen::MatrixBase<Derived>& A, idx n) {
  // check zero-size
  if (!internal::check_nonzero_size(A))
    throw exception::ZeroSize("clara::powm()");
  // check square matrix
  if (!internal::check_square_mat(A))
    throw exception::MatrixNotSquare("clara::powm()");
  if (n == 1)
    return A;
  dyn_mat<typename Derived::Scalar> result =
      dyn_mat<typename Derived::Scalar>::Identity(A.rows(), A.rows());
  if (n == 0)
    return result;

  dyn_mat<typename Derived::Scalar> cA = A.derived();

  // fast matrix power using SQUARE-AND-MULTIPLY algorithm
  for (; n > 0; n /= 2) {
    if (n % 2)
      result = (result * cA).eval();
    cA = (cA * cA).eval();
  }
  return result;
}

/**
 * @brief schatten matrix norm
 * this function computes the schatten matrix nrom of a matrix 'A' with a given exponent 'p'.
 * the schatten matrix norm is defined as the p-norm of the vector singular value of 'A'
 * its a generalization of the matrix frobenius norm (p = 2) and is well-defined for any position
 * value of 'p the function ensures that the input matrix is non-zero and has valid exponent 'p'
 * before performing the computation
 *
 * @tparam Derived the derived tpe of the input matrix 'A'
 * @param A the input matrix for which the schatten matrix norm is computed
 * @param p is the exponent of the schatten norm, where p >= 1
 * @return the schatten matrix norm of the input matrix 'A' with the given exponent 'p'
 *
 * @throws exception::ZeroSize if the input matrix 'A' has zero size
 * @throws exception::OutOfRange if the exponent 'p' less than 1
 *
 * @example
 * Eigen::Matrix2d inputMatrix;
 * inputMatrix << 1, 2, 3, 4;
 *
 * double p = 1.5;
 * // compute the schatten matrix nrom of the input matrix with exponent 'p'
 * double result = schatten(inputMatrix, p);
 *
 * std::cout << "schatten matrix norm "<< result << std::endl;
 */
template <typename Derived>
double schatten(const Eigen::MatrixBase<Derived>& A, double p) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::schatten()");
  if (p < 1)
    throw exception::OutOfRange("clara::schatten()");

  const dyn_col_vect<double> sv = svals(rA);
  double result = 0;
  for (idx i = 0; i < static_cast<idx>(sv.rows()); ++i)
    result += std::pow(sv[i], p);
  return std::pow(result, 1. / p);
}

template <typename OutputScalar, typename Derived>
dyn_mat<OutputScalar> cwise(const Eigen::MatrixBase<Derived>& A,
                            OutputScalar (*f)(const typename Derived::Scalar&)) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::cwise()");
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
 * @brief kronecker product of multiple matrices, prserving the return type
 * this function compute the kronecker product of multiple matrices using a variadic template
 * the kronecker product of two matrices A and B is denoted as A \tensorprod B result in
 * block matrix. ot preserve the return ype and ensures the input matrices have compatible
 * dimensions
 *
 * @tparam T the first matrix type
 * @param head the first matrix in the krocnker product sequence
 * @return the kronecker product of the input matrices, preseving the return type
 *
 * @example
 * Eigen::<atrix2d A;
 * A << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d B;
 * B << 5, 6, 7, 8;
 * Eigen::Matrix4d C = clara::kron(A, B);
 * std::cout << "kroncker product of A and B:\n" << C << std::endl;
 */
template <typename T>
dyn_mat<typename T::Scalar> kron(const T& head) {
  return head;
}

/**
 * @brief kronecker product of multiple matrices, preserving the return type
 * this function computes the kronecker product of multiple matrices using a variadic template.
 * the kroncker product of two matrices A and B is denoted as A \tensorprod B and result in
 * block matrix. it preserves the return type and ensures that the input matrices have compatible
 * dimension
 *
 * @tparam T the matrix type
 * @tparam Args the types of additional matrices of the compute the Kroncker Product
 * @param head the first matrix in the kronecker product sequence
 * @param tail additional matrices to compute the kronecker product
 * @return the kronecker product of the input matrices, preserving the return type
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d B;
 * B << 5, 6, 7, 8;
 *
 * Eigen::Matrix2d C;
 * C << 9, 10, 11, 12;
 *
 * Eigen::Matrix8d result = clara::kron(A, B, C)
 * std::cout << "kronecker product of A, B, C: \n" << result << std::endl;
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
    throw exception::ZeroSize("clara::kron()");

  for (auto&& it : As)
    if (!internal::check_nonzero_size(it))
      throw exception::ZeroSize("clara::kron()");

  dyn_mat<typename Derived::Scalar> result = As[0].derived();
  for (idx i = 1; i < As.size(); ++i) {
    result = kron(result, As[i]);
  }
  return result;
}

/**
 * @brief kronecker prodcut of a list of matrices, preserving the return type
 * this function computes the kronecker product of a list of matrices provided in an
 * initializer_list. the kronecker product of two matrices A and B is denoted as A \tensorprod B and
 * result in block matrix. it preserves the return type and ensures that the input matrices have
 * compatible the dimension
 *
 * @tparam Derived the type of matrices initializer_list
 * @param As the list of matrices to compute the kronecker product.
 * @return the kronecker of the matrices in the initializer_list, preserving the return type
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * @Eigen::Matrix2d B;
 * B << 5, 6, 7, 8;
 *
 * Eigen::Matrix2d C;
 * C << 9, 10, 11, 12;
 *
 * Eigen::Matrix8d result = clara::kron({A, B, C});
 * std::cout << "kronecker product of A, B, and C:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> kron(const std::initializer_list<Derived>& As) {
  return kron(std::vector<Derived>(As));
}

/**
 * @brief kronecker product of matrix with itself 'n' times is a dynamic matrix over the same
 * scalar field
 *
 * this function computes the kroncker product power of a matrix 'A' with itself 'n' times. the
 * kronekcer priduct of a matrix A with itself 'n' times is denoted as A \tensorprod \ A \tensorprod
 * ... (n times), and it result in a block matrix with 'n' factors of 'A' in the product
 *
 * @tparam Dervied the type of the input matrix 'A'
 * @param A the matrix to compute the kronecker product power with
 * @param n the number of times to compute the kroncker product with 'A'
 * @return the kronecker product power 'A' with itself 'n' times as dynamic matrix over the same
 *            acalar the field as 'A'
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * int n = 3;
 * Eigen::Matrix8d result = clara::kronpow(A, n);
 * std::cout << "kronecker product power A with itself 3 times:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> kronpow(const Eigen::MatrixBase<Derived>& A, idx n) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::kronpow()");
  if (n == 0)
    throw exception::OutOfRange("clara::kronpow()");
  std::vector<dyn_mat<typename Derived::Scalar>> As(n, rA);
  return kron(As);
}

/**
 * @brief direct sum of multiple matrices, preserving the return type
 * this function computes the direct sum of mulitiple matrices, preserving the return type
 * the direct sum of two matrices A and B, denoted as A \tensorprod B in a block matrix where A and
 * B placed diagonally
 *
 * @tparam T the type of the first matrix
 * @param head the first matrix to start the direct sum operation
 * @return the direct sum of multiple matrices, preserving the return type
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d B;
 * B << 5, 6, 7, 8;
 *
 * Eigen::Matrix2d C;
 * C << 9, 10, 11, 12;
 */
template <typename T>
dyn_mat<typename T::Scalar> dirsum(const T& head) {
  return head;
}

/**
 * @brief direct sum of multiple matrices, preserving the return type
 *
 * this function computes the direct sum multipe matries, preserving the reutn type
 * the direct sum of two matrices A and B, denoted as A \tensorprod B, result in block matrix where
 * A and B place diagonally. the direct sum of multiple matrices
 * A, B, C, ..., denoted as A \tensorprod B \tensorprod C ..., is the
 * concatention of their block diagonal matrics
 *
 * @tparam T the first matrix
 * @param head the first matrix to start the direct sum operation
 * @return the direct sum of multiple matrices, preserving the return type
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d B;
 * B << 5, 6, 7, 8;
 *
 * Eigen::Matrix2d C;
 * C << 9, 10, 11, 12;
 *
 * Eigen::Matrix6d result = clara::dirsum(A, B, C);
 * std::cout << "direct sum of A, B, and C:\n" << result << std::endl;
 */
template <typename T, typename... Args>
dyn_mat<typename T::Scalar> dirsum(const T& head, const Args&... tail) {
  return internal::dirsum2(head, dirsum(tail...));
}

/**
 * @brief direct sum of matrices
 *
 * this function computes the direct of all elements in the input vector 'As'
 * evaluated from left to right, as dynamic matrix over the same scalar field as its arguments.
 * the direct sum of two matrices A and B, denoted as A \tensorprod B, result in block
 * matrix where A and B are placed diagonally. The direct sum of multiple matrices A, B, C ...,
 *
 * @tparam Derived the type of the matrice in the vector
 * @param As a vector containing matrices to perform the direct sum operation
 * @return the direct sum of all elements in the vector 'As'
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d B;
 * B << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d C;
 * C << 9, 10, 11, 12;
 *
 * std::vector<Eigen::Matrix2d> matrices = {A, B, C};
 * Eigen::Matrix6d result = clara::dirsum(matrices);
 *
 * std::cout << "direct of A, B, and C:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> dirsum(const std::vector<Derived>& As) {
  if (As.size() == 0)
    throw exception::ZeroSize("clara::dirsum()");
  for (auto&& it : As)
    if (!internal::check_nonzero_size(it))
      throw exception::ZeroSize("clara::dirsum()");
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
 * @brief direct sum of matrices
 * this function computes the direct sum of all elements in the nput initializer_list 'As'
 * evaluated from left to right.
 * as dynamic matrix over the same scalar field as its arguments. the direct sum of two matrices
 * denoted as A \tensorprod B, result in a block matrix where A and B are placed diagonally
 * the direct sum of multiplying matrices A, B, C ..., denoted as A \tensorprod B \tensorprod C
 * \tensorprod ..., is the concatention of their block
 *
 * @tparam Dervied the type of the matrices in the initializer_list
 * @param As an initializer_list containing matrices to pefrom the direct sum operation
 * @return the direct sum of all elements in the input initializer_list 'As'
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d B;
 * B << 5, 6, 7, 8;
 *
 * Eigen::Matrix2d C;
 * C << 9, 10, 11, 12;
 * std::cout << "direct sum of A, B, C:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> dirsum(const std::initializer_list<Derived>& As) {
  return dirsum(std::vector<Derived>(As));
}

/**
 * @brief direct sum power of a matrix
 * this function computes the direct sum of the input matrix A with itself 'n' times
 * as a dynamic matrix over the same scalar field as A.
 * the direct sum power of a matrix A with itself n times, denoted as A \tensorprod A (n times),
 * result in a block diagonal matrix where A appears on the main diagonal `n` times
 *
 * @tparam Derived the type of the input matrix.
 * @param A the input matrix to perform the direct sum power operation on
 * @param n the number of times the input matrix A is directly summed with itself
 * @return the direct sum of the input matrix A with itself 'n' times
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * idx n = 3;
 * auto result = clara::dirsumpow(A, n);
 *
 * std::cout << "direct sum power A (n = 3):\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> dirsumpow(const Eigen::MatrixBase<Derived>& A, idx n) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::dirsumpow()");
  if (n == 0)
    throw exception::OutOfRange("clara::dirsumpow()");
  std::vector<dyn_mat<typename Derived::Scalar>> As(n, rA);
  return dirsum(As);
}

/**
 * @brief reshape a matrix
 *
 * this function reshape the input matrix A to a new shape defined by the number of rows and columns
 * the reshape operation preserve the order of elements in the original matrix and uses oclumn-major
 * order.
 * @tparam Derived the type of the input matrix
 * @param A the input matrix to reshaped
 * @param rows the number of reows in the reshape matrix
 * @param cols the number of columns in the reshape matrix
 * @return the reshape matrix with the specified number of rows and columns
 *
 * @example
 * Eigen::Matrix3d A;
 * Eigen << 1, 2, 3, 4, 5, 6, 7, 8, 9;
 *
 * idx rows = 2;
 * idx cols = 4;
 * auto result = clara::reshape(A, rows, cols);
 * std::cout << "reshaped matrix:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> reshape(const Eigen::MatrixBase<Derived>& A, idx rows, idx cols) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  idx Arows = static_cast<idx>(rA.rows());
  idx Acols = static_cast<idx>(rA.cols());

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::reshape()");
  if (Arows * Acols != rows * cols)
    throw exception::DimsMismatchMatrix("clara::reshape()");
  return Eigen::Map<dyn_mat<typename Derived::Scalar>>(
      const_cast<typename Derived::Scalar*>(rA.data()), rows, cols);
}

/**
 * @brief calculate the comutator of two matrices A and B
 * the comutator of two square matricesA and B is given by [A, B] = AB - BA
 * this function returns the comutator of the input matrices A and B
 *
 * @tparam Derived1 the type of the first input matrix A
 * @tparam Dertived2 the type of the second input matrix B
 * @param A the first input matrix A
 * @param B the second input matrix B
 * @return the comutator matrices A and B
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d B;
 * B << 5, 6, 7, 8;
 * auto result = clara::comm(A, B);
 *
 * std::cout << "comutator of A and B:\n" << result << std::endl;
 */
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> comm(const Eigen::MatrixBase<Derived1>& A,
                                        const Eigen::MatrixBase<Derived2>& B) {
  const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
  const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

  if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
    throw exception::TypeMismatch("clara::comm()");
  if (!internal::check_nonzero_size(rA) || !internal::check_nonzero_size(rB))
    throw exception::ZeroSize("clara::comm()");

  // check square matrices
  if (!internal::check_square_mat(rA) || !internal::check_square_mat(rB))
    throw exception::MatrixNotSquare("clara::comm()");

  // check equal dimension
  if (rA.rows() != rB.rows())
    throw exception::DimsNotEqual("clara::comm()");
  return rA * rB - rB * rA;
}

/**
 * @brief calculate the anti comutator of two matrices A and B
 * the antri-comutator of two matrices A and B is given by {A, B} = AB + BA
 * this function returns the anti-commutator of the input matrices A and B
 *
 * @tparam Derived1 the type of the first input matrix A
 * @tparam Derived2 the type of the second input matrix B
 * @param A the first input matrix A
 * @param B the second input matrix B
 * @return the anti-commutator of matrices A and B
 *
 * @example
 * Eigen::Matrix2d A;
 * A << 1, 2, 3, 4;
 *
 * Eigen::Matrix2d B;
 * B << 5, 6, 7, 8;
 *
 * auto result = clara::anticomm(A, B);
 * std::cout << "anti-commutator of A and B: \n" << rsult << std::endl;
 */
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> anticomm(const Eigen::MatrixBase<Derived1>& A,
                                            const Eigen::MatrixBase<Derived2>& B) {
  const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
  const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

  if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
    throw exception::TypeMismatch("clara::anticomm()");

  if (!internal::check_nonzero_size(rA) || !internal::check_nonzero_size(rB))
    throw exception::ZeroSize("clara::anticomm()");

  // check square matrix
  if (internal::check_square_mat(rA) || !internal::check_square_mat(rB))
    throw exception::MatrixNotSquare("clara::anticomm()");

  // check equal dimension
  if (rA.rows() != rB.rows())
    throw exception::DimsNotEqual("clara::anticomm()");
  return rA * rB + rB * rA;
}

/**
 * @brief calculate the projection matrix of a complex column vector
 * the projection matrix of a complex column vector A is given by Prj(A) = A * adjoint / ||A||^2
 * where ||A|| is the euclidean norm of A. if A is a zero vector, the function
 * return zero matrix
 *
 * @tparam Dervied the type of the input column vector A
 * @param A the input complex column vector A
 * @return the projection matrix Prj(A) = A * adjoint(A) / ||A||^2
 *
 * @example
 * Eigen::VectorXcd A(3);
 * A << std::complex<double>(1.0, 2.0),
 *      std::complex<double>(3.0, 4.0),
 *      std::complex<double>(5.0, 6.0);
 * auto result = clara::prj(A);
 * std::cout << "projection matrix:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> prj(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::prj()");

  // check column vector
  if (!internal::check_cvector(rA))
    throw exception::MatrixNotCvector("clara::prj()");
  double normA = norm(rA);
  if (normA > eps)
    return rA * adjoint(rA) / (normA * normA);
  else
    return dyn_mat<typename Derived::Scalar>::Zero(rA.rows(), rA.rows());
}

/**
 * @brief perform Gram-Schmidt orthogonalization on a list of column vectors
 * given a list of non-zero complex column vectors, this function perform
 * Gram-Schmidt orthogonalization to obtain an orthonormal set of vectors. the
 * funcion returns a matrix where each column represent an orthonormal vector from the
 * original list. if some vectors in the list are linearly dependent, they will not
 * be included in the output matrix
 *
 * @tparam Derived the type of the input column vectors
 * @param As a list of non-zero complex column vector
 * @return A matrix with orthonormal vectors obtained from Gram-Schmidt orthogonalization
 *
 * @example
 * Eigen::VectorXcd A(3), B(3), C(3);
 * A << std::complex<double>(1.0, 2.0),
 *      std::complex<double>(3.0, 4.0),
 *      std::complex<double>(5.0, 6.0);
 * B << std::complex<double>(-1.0, 1.0),
 *      std::complex<double>(2.0, -2.0),
 *      std::complex<double>(-3.0, 3.0);
 * C << std::complex<double>(1.0, -1.0),
 *      std::complex<double>(2.0, -2.0),
 *      std::complex<double>(3.0, -3.0);
 *
 * std::vector<Eigen::VectorXcd> vectors = {A, B, C};
 * auto result = clara::grams(vectors);
 * std::cout << "orthonormal vectors:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> grams(const std::vector<Derived>& As) {
  if (!internal::check_nonzero_size(As))
    throw exception::ZeroSize("clara::grams()");

  for (auto&& it : As)
    if (!internal::check_nonzero_size(it))
      throw exception::ZeroSize("clara::grams()");

  // check that As[0] is a column vector
  if (!internal::check_cvector(As[0]))
    throw exception::MatrixNotCvector("clara::grams()");

  // check that all rest match the size of the first vector
  for (auto&& it : As)
    if (it.rows() != As[0].rows() || it.cols() != 1)
      throw exception::DimsNotEqual("clara::grams()");

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
 * @brief perform Gram-Schmidt orthogonalization on a list column vectors.
 * given a list of non-zero column vectors, this function pefroms Gram-Schmidt
 * orthogonalization to obtrain an orthonormal set of vectors. the function returns
 * a matrix where where each column represent an orthonormal vector from original
 * list. if some vectors in the list are linearly dependent, they will not be included
 * in the output matrix.
 *
 * @tparam Derived the type of the input column vectors
 * @param As A list of non-zero complex column vectors
 * @return a mwtrix with orthonormal vectors obrained from gram-shmidt orthogonalization
 *
 * @example
 * Eigen::VectorXcd A(3), B(3), C(3);
 * A << std::complex<double>(1.0, 2.0),
 *      std::complex<double>(3.0, 4.0),
 *      std::complex<double>(5.0, 6.0);
 * B << std::complex<double>(-1.0, 1.0),
 *      std::complex<double>(2.0, -2.0),
 *      std::complex<double>(-3.0, 3.0);
 * C << std::complex<double>(1.0, -1.0),
 *      std::complex<double>(2.0, -2.0),
 *      std::complex<double>(3.0, -3.0);
 *
 * auto result = clara::grams({A, B, C});
 * std::cout << "orthonormal vectors:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> grams(const std::initializer_list<Derived>& As) {
  return grams(std::vector<Derived>(As));
}

/**
 * @brief perform Gram-Schmidt orthogonalization on the column of a matrix
 * given matrix with non-zero complex columns, this function perform Gram-Schmidt
 * orthogonalization on the columns to obtain an orthonormal set of vectors. the function
 * returns a matrix where each column represent orthonormal vector from the original matrix.
 * if some columns in the matrix are linearly dependent, they will not be included in
 * the output matrix
 *
 * @tparam Derived the type of the input matrix
 * @param A matrix with non-zero complex columns
 * @return A matrix with orthonormal vectors obtained from Gram-Schmidt orthogonalization
 *
 * @example
 * Eigen::MatrixXcd A(3, 3);
 * A << std::complex<double>(1.0, 2.0), std::complex<double>(-1.0, 1.0), std::complex<double>(1.0,
 -1.0),
 *      std::complex<double>(3.0, 4.0), std::complex<double>(2.0, -2.0), std::complex<double>(2.0,
 -2.0),
 *      std::complex<double>(5.0, 6.0), std::complex<double>(-3.0, 3.0), std::complex<double>(3.0,
 -3.0);
 * auto result = clara::grams(A);
 * std::cout << "orthonormal vectors:\n" << result << std::endl;
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> grams(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::grams()");
  // prepare a vector to store individual column vectors of the input matrix
  std::vector<dyn_mat<typename Derived::Scalar>> input;
  for (idx i = 0; i < static_cast<idx>(rA.cols()); ++i)
    input.push_back(rA.cols(i));
  // use the overloaded 'grams' functions that takes std::vector as input
  return grams<dyn_mat<typename Derived::Scalar>>(input);
}

/**
 * @brief convert a non-negative integer index to multi-index
 * given a non-negative integer index 'n' and vector 'dims' representing the dimension
 * of a multidimensional array, this function returns the corresponding multi-index.
 * the multi-index is a vector with same size as 'dims' , containing indices that repersent
 * the position of 'n' in the multidimensional array
 *
 * @param n the non-negative integer index to be converted to a multi-index
 * @param dims A vector repersenting the dimension of a multidimensional array
 * @return the multi-index representing the position of 'n' in the multidimensional array.
 *
 * @throws exception::DimsInvalid if the dimension vector 'dims' is invalid
 * @throws exception::OutOfRange if the index 'n' is greater than the or equal to the
 *                            tpta; number of elements in multi index
 * @example
 * std::vector<idx> dims = {3, 4, 2};
 * idx n = 7;
 * std::vector<idx> result = clara::n2multiidx(n, dims);
 *
 * // print result
 * std::cout << "mult index for index: ";
 * for (idx i : result) {
 *    std::cout << i << " ";
 * }
 */
inline std::vector<idx> n2multiidx(idx n, const std::vector<idx>& dims) {
  // check if the dimension vector 'dims' is valid (non-empty and containing non-zero element)
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::n2multiidx()");
  if (n >= std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                           std::multiplies<idx>()))
    throw exception::OutOfRange("clara::n2multiidx()");
  idx result[2 * maxn];
  internal::n2multiidx(n, dims.size(), dims.data(), result);
  // create a vector from resulting array and return it as the multi-index
  return std::vector<idx>(result, result + dims.size());
}

/**
 * @brief convert a multi-index to a non-negative integer index
 * given a multi-index 'midx' and a vector 'dims' repersenting the dimension of
 * multidimensional array, this function returns the corresponding non-negative integer
 * index. the multi-index is a vector containing indices that repersent the position
 * of a specific element in the multidimensional array. the non-negative integer index is the unique
 * identifier of the element within the array, obtained by mapping the multi-index to a linear
 * index.
 *
 * @param midx the multi index representing the position of an element in the multidimensional array
 * @param dims a vector represeinting the dimension of the multidimensional array
 * @return the non-negative integer index of the element
 *
 * @throws exception::DimsInvalid if the dimension vectors 'dims' is invalid
 * @throws exception::OutOfRange if any element of the multi-index 'midx' exceeds the corresponding
 * dimension dims
 *
 * @example
 * std::vector<idx> dims = {3, 4, 2};
 * std::vector<idx> midx = {1, 3, 1};
 *
 * idx n = clara::multiidx2n(midx, dims);
 * std::cout << "non negative integer index for multi index {1, 3, 1}: " << n << std::endl;
 */
inline idx multiidx2n(const std::vector<idx>& midx, const std::vector<idx>& dims) {
  // check if the dimension vector 'dims' is valid
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::multiidx2n()");
  // check if any alement of the multi-index 'midx' exceeds the corresponding
  // dimension in 'dims'
  for (idx i = 0; i < dims.size(); ++i)
    if (midx[i] >= dims[i])
      throw exception::OutOfRange("clara::multiidx2n()");
  // convert the multi-index to a non-negative integer index using internal helper
  // function
  return internal::multiidx2n(midx.data(), dims.size(), dims.data());
}

/**
 * @brief a multi-partite qudit ket presented by the multi-index 'mask' and the dimension of each
 * subsystem the dimension of each subsystem 'dims', this function constructs the corresponding
 * qudit ket in hilbert space. the multi-index 'mask' is a vector containing non-negative integers,
 * where each element repersent is a vector containing non-negative integers, where each element
 * repersent the basis state index for the corresponding subsystem. the dimensions 'dims' is a
 * vector representing the number of basis states in each subsystem. the resulting qudit ket
 * repersented as a column vectro in hilbert space
 *
 * @throws clara::ZeroSize if the size of 'mask' or 'dims' is zero
 * @throws clara::DimsInvalid if 'dims' contains invalid values
 * @throws clara::SubsysMismatchdims if 'mask' and 'dims' have different size, or if any element in
 * 'mask' exceeds the corresponding dimension in 'dims'
 *
 * @example
 * std::vector<idx> dims = {2, 3, 2};
 * std::vector<idx> mask = {1, 0, 1};
 * ket qudit_ket = clara::mket(mask, dims);
 * std::cout << "multi-partite qudit ket for mask {1, 0, 1}:\n" << qudit_ket << std::endl;
 */
inline ket mket(const std::vector<idx>& mask, const std::vector<idx>& dims) {
  // get the number of subsystem and the total dimension of the hilbert space
  idx N = mask.size();
  idx D = std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                          std::multiplies<idx>());

  if (N == 0)
    throw exception::ZeroSize("clrara::mket()");
  // check valid dims
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::mket()");
  // check mask and dims have the same size
  if (mask.size() != dims.size())
    throw exception::SubsysMismatchdims("clara::mket()");
  // check mask is a valid vector
  for (idx i = 0; i < N; ++i)
    if (mask[i] >= dims[i])
      throw exception::SubsysMismatchdims("clara::mket()");

  // create a ket vector with zero entries and set the position corresponding to the multi-index
  // 'mask' to 1
  ket result = ket::Zero(D);
  idx pos = multiidx2n(mask, dims);
  result(pos) = 1;
  return result;
}

/**
 * @brief construct a multi-partite qudit ket with equal dimension for all subsystem
 * given a multi-partite qudit ket repersented by the multi-index 'mask' and a common
 * dimension 'd', this function construct the corresponding qudit ket in the hilbert space.
 * the multi-index 'mask' is a vector containing non-negative integers, where each element
 * represents the basis state index for the corresponding subsystem. the common dimension 'd'
 * represent the number of basis states in each subsystem, and all subsystem have the same dimension
 * 'd'. the resulting qudit ket is represented as a column vector in the hilbert space
 *
 * @param mask the multi-index repersenting the basis state indices for each subsystem
 * @param d the common dimension of all subsystem. default is 2
 * @return the multi-partite qudit ket as a column vector in the hilbert space
 *
 * @throws exception::ZeroSize if the size of 'mask' is zero
 * @throws exception::DimsInvalid if 'd' is zero
 * @throws exception::SubsysMismatchdims if any element in 'mask' is greater than or equal to 'd'
 *
 * @example
 * // multi-index repersenting basis state indices for each subsystem
 * std::vector<idx> mask = {1, 0, 1};
 * // common dimension for all subsystem
 * idx d = 3;
 * ket qudit = clara::mket(mask, d);
 * std::cout << "multi-partite qudit ket for mask {1, 0, 1} and dimension 3: \n" << qudit_ket;
 */
inline ket mket(const std::vector<idx>& mask, idx d = 2) {
  // get the number of subsystem and the total dimension of the hilbert space
  idx N = mask.size();
  idx D = static_cast<idx>(std::llround(std::pow(d, N)));

  // check zero-size
  if (N == 0)
    throw exception::ZeroSize("clara::mket()");
  // check valid dims
  if (d == 0)
    throw exception::DimsInvalid("clara::mket()");

  for (idx i = 0; i < N; ++i)
    if (mask[i] >= d)
      throw exception::SubsysMismatchdims("clara::mket()");

  // create a ket vector with zero entries and set the position corresponding
  // to the multi-index 'mask' to 1
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
    throw exception::ZeroSize("clara::mprj()");
  // check valid dims
  if (!internal::check_dims(dims))
    throw exception::ZeroSize("clara::mprj()");
  // check mask and dims have the same size
  if (mask.size() != dims.size())
    throw exception::SubsysMismatchdims("clara::mprj()");
  // check mask is valid vector
  for (idx i = 0; i < N; ++i)
    if (mask[i] >= dims[i])
      throw exception::SubsysMismatchdims("clara::mprj()");
  cmat result = cmat::Zero(D, D);
  idx pos = multiidx2n(mask, dims);
  result(pos, pos) = 1;
  return result;
}

/**
 * @brief construct the projector onto multi-partite qudit ket
 * given a multi-partite qudit ket repersented the mutli-index 'mask'
 * and a vector of dimension 'dims', this function construct the corresponding porjector
 * onto the multi-partite qudit ket in the hilbert space. the multi-index 'mask' is
 * a vector containing non-negative integers, where each element repersents the basis
 * state index for the corresponding subsystem. the vector 'dims' contains the dimension of all
 * subsystem, and each element in 'mask' must be strictly smaller than the corresponding element
 * in 'dims'. the resulting projector is represented as a complex-valued matrix
 * in the hilber space
 *
 * @param mask the multi-index representing the basis state indices for each subsystem
 * @param dims the dimension of all subsystem
 * @return the projector onto the multi-partite qudit ket as a complex-valued matrix
 *
 * @throws clara::ZeroSize if the size of 'mask' is zero
 * @throws clara::ZeroSize if the size of 'dims' is zero
 * @throws clara::SubsysMismatchdims if the size of 'mask' and 'dims' do not matcj
 *
 * @example
 * // multi-index representing basis state indices for each subsystem
 * std::vector<idx> mask = {1, 0, 1}
 * // dimension of each subsystem
 * std::vector<idx> dims = {2, 3, 2};
 * cmat projector = clara::mprj(mask, dims);
 *
 * std::cout << "projector into multi-partite qudit ket for mask {1, 0, 1} and dims {2, 3, 2}: \n"
 *          << projector << std::endl;
 */
inline cmat mprj(const std::vector<idx>& mask, idx d = 2) {
  idx N = mask.size();
  idx D = static_cast<idx>(std::llround(std::pow(d, N)));

  // check zero size
  if (N == 0)
    throw exception::ZeroSize("clara::mprj()");

  // check valid dims
  if (d == 0)
    throw exception::DimsInvalid("clara::mprj()");

  // check mask is a valid vector
  for (idx i = 0; i < N; ++i)
    if (mask[i] >= d)
      throw exception::SubsysMismatchdims("clara::mprj()");

  // create a complex-valued matrix with zero entries and set the position
  // corresponding to the multi-index 'mask'
  cmat result = cmat::Zero(D, D);
  std::vector<idx> dims(N, d);
  idx pos = multiidx2n(mask, dims);
  result(pos, pos) = 1;
  return result;
}

/**
 * @brief calculate the absolute square of complex numbers in the range [first, last]
 * given a range of complex numbers specified by iterator, first and last, this function
 * calculate the absolute square (modulus squared) of each complex number and returns a vector
 * containing the result. the input range should contain complex number, and the output vector
 * will contain the corresponding absolute square as double values
 *
 * @tparam InputIterator iterator type for the input range of complex number
 * @param first iterator to the beginning of the input range
 * @param last iterator to the end of the last input range
 * @retun a vector containing the absolute squares of complex numbers in the input range
 *
 * @example
 * std::vector<cplx> complexNumber = {cplx(1.0, 2.0), cplx(-3.0, 4.0), cplx(0.0, -1.0)};
 * std::vector<double> absoluteSquare = abssq(std::begin(complexNumber), std::end(complexNumber));
 *
 * for (double absSq : absoluteSquare) {
 *  std::cout << "absolute square " << absSq << std::endl;
 * }
 */
template <typename InputIterator>
std::vector<double> abssq(InputIterator first, InputIterator last) {
  // calculate the number of complex numbers in the input range
  std::vector<double> weights(std::distance(first, last));
  // calculate the absolute to store the absolute the complex number
  std::transform(first, last, std::begin(weights), [](cplx z) -> double { return std::norm(z); });
  return weights;
}

/**
 * @brief calculate the absolute square of element in an STL-like container
 * given an STL-like container 'c' containing elements that support the 'std::begin()' and
 * 'std::end()' functions (e.g std::vector, std::list, etc), this function calculates, the absolute
 * square (modulus squared) of each element in the container and returns a vector containing the
 * result the input container should contain elements that can be square, and the output vector will
 * container the corresponding absolute square as double values
 *
 * @tparam container the type of the input STL-like container
 * @param c the STL-like container to calculate the absolute square for
 * @return vector containing the absolute square of the elements in the input container
 *
 * @example
 * std::vector<double> values = {1.0, -3.0, 2.5, -4.5};
 * std::vector<double> absoluteSquare = abssq(values);
 *
 * for (double absSq : absoluteSquare) {
 *  std::cout << "absolute square: " << absSq << std::endl;
 * }
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
    throw exception::ZeroSize("clara::abssq()");
  return abssq(rA.data(), rA.data() + rA.size());
}

/**
 * @brief calculate the absolute of elements in an Eigen matrix
 * given an Eigen matrix 'A', this function calculates the absolute square (modulus squared)
 * of each element in the matrix and returns a vector containing the result. the input matrix
 * should have elements that can be squared, and the output vector will contain the corresponding
 * absolute squares as double values
 *
 * @tparam Derived the type of the Eigen matrix
 * @param A the eigen matrix to calculate the absolute square for
 * @return vector containing the absolute square of the elements in the input matrix
 *
 * @example
 * Eigen::Matrix2d matrix;
 * matrix << 1.0, -3.0, 2.5, -4.5;
 * std::vector<double> absoluteSquares = abssq(matrix);
 *
 * for (double absSq : absoluteSquares) {
 *  std::cout << "absolute square: " << absSq << sd::endl;
 * }
 */
template <typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type sum(InputIterator first,
                                                             InputIterator last) {
  using value_type = typename std::iterator_traits<InputIterator>::value_type;
  return std::accumulate(first, last, static_cast<value_type>(0));
}

/**
 * @brief calculate the element wise sum of elements in an STL-like containers
 * given an STL-like container 'c', this function calculates the element-wise sum
 * of all elements in the container and returns the result as a scalar over the same scalar
 * field as container's element. the input container sougld have elements that can be summed,
 * and the output will be of the same type as the container's elements.
 *
 * @tparam container the type of the STL-like container
 * @param c the STL-like container to calculate the element-wise sum for
 * @return the element-wise sum of all elements in the container
 *
 * @example
 * std::vector<int> numbers = {1, 2, 3, 4, 5};
 *
 * int sumResult = sum(numbers);
 * std::cout << "sum: " << sumResult << std::endl;
 */
template <typename Container>
typename Container::value_type sum(
    const Container& c, typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  return sum(std::begin(c), std::end(c));
}

/**
 * @brief calculate the element-wise sum of elements in an STL-like container
 * given an STL-like range specified by the iterators 'first' and 'last', this function
 * calculate the element-wise prodcut of all elements in the range and return the result
 * as a scalar over the same scalar field as the range's element. the input range should
 * have elements that can be multiplied, and the output will be of the same type as range
 * elements.
 *
 * @tparam InputIterator type of the iterator for the STL-like range
 * @param first iterator pointing to the first element of the range
 * @param last iterator pointing to the position after the last element of the range
 * @return the element-wise product of elements in the range
 *
 * @example
 * std::vector<int> number = {1, 2, 3, 4, 5};
 *
 * int producttresult = prod(number.begin(), number.end());
 * std::cout << "product: " << producttresult << std::endl;
 */
template <typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type prod(InputIterator first,
                                                              InputIterator last) {
  using value_type = typename std::iterator_traits<InputIterator>::value_type;
  // using the element-wise product using std::accumulate and std::multiplies
  // starting from 1 to act as an identity element multiplication
  return std::accumulate(first, last, static_cast<value_type>(1), std::multiplies<value_type>());
}

/**
 * @brief calculate the element-wise product of an STL-like container
 * given an STL-like container and returns the result as a scalar over the same scalar
 * field as the container's elements. the container's lement should be able to be multiplied
 * , and the output will be of the same type as the container's elements
 *
 * @tparam container the type of the STL-like container
 * @param c the STL-like container for which the element-wise product is calculated
 * @return the element-wise prodcut of all elements in containers
 *
 * @examplea
 * std::vector<int> numbers = {1, 2, 3, 4, 5};
 * int productResult = prod(numbers);
 * std::cout << "Product: " << productResult << std::endl;
 */
template <typename Container>
typename Container::value_type prod(
    const Container& c, typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  // call the prod function with iterator of the container to calculate the result
  return prod(std::begin(c), std::end(c));
}

/**
 * @brief convert a density matrix (hermitian matrix) to pure state vector (ket)
 * given a density matrix 'A' (hermitian matrix), this function calculate the corresponding
 * pure state vector (ket) representing the most significant eigenvector of the matrix. it
 * returns the pure state vector as a dynamic column vector over the same scalar field vector
 * as a dynamic column vector over same scalar field as the input matrix 'A'
 *
 * @tparam Derived the matrix of type of the input 'A'
 * @param A the density matrix (hermitian matrix) to convert to a pure state vector
 * @return the pure state vector (ket) representing the most significant eigenvectors of 'A'
 *
 * @example
 * // 3x3 density matrix (hermitian matrix)
 * Eigen::Matrix3cd rho
 *
 * Eigen::Vector3cd pureState = rho2pure(rho);
 * std::cout << "pure state vector:\n" << pureState << std::endl;
 */
template <typename Derived>
dyn_col_vect<typename Derived::Scalar> rho2pure(const Eigen::MatrixBase<Derived>& A) {
  // obtain a const reference to the derived matrix type from the input 'A'
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::rho2pure()");
  // check square matrix
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::rho2pure()");

  // calculate the eigenvalues and eigenvectors of the hermitian matrix 'A'
  dyn_col_vect<double> tmp_evals = hevals(rA);
  cmat tmp_evects = hevects(rA);
  dyn_col_vect<typename Derived::Scalar> result =
      dyn_col_vect<typename Derived::Scalar>::Zero(rA.rows());

  // find the most significant eigenvectors by looking for the first non-zero eigenvalues
  for (idx k = 0; k < static_cast<idx>(rA.rows()); ++k) {
    if (std::abs(tmp_evals(k)) > eps) {
      result = tmp_evects.col(k);
      break;
    }
  }

  // return the pure state vector representing the most significant eigenvectors of 'A'
  return result;
}

/**
 * @brief construct the complement of a subsystem vector
 * given a subsystem vector 'subsys' and the total number of elements 'N', this function
 * construct and returns the complement of 'subsys' with respect to the set {0, 1 ..., N - 1}.
 * the complemment is a vector containing elements from {0, 1, ..., N - 1} that are not present
 * in 'subsys'
 *
 * @tparam T the data of the elements in the 'subsys' vector
 * @oaram subsys subsystem vector whose complement is to be computed
 * @param N total number of element in the set {0, 1, ..., N - 1}
 * @return the complement of 'subsys' vector
 *
 * @example
 * std::vector<int> subsystem = {1, 3};
 * int totalElements = 5;
 * std::vector<int> complement = complement(subsystem, totalElements);
 */
template <typename T>
std::vector<T> complement(std::vector<T> subsys, idx N) {
  if (N < subsys.size())
    throw exception::OutOfRange("clara::complement()");
  std::vector<T> all(N);
  std::vector<T> subsys_bar(N - subsys.size());

  // generate a Vector containing elements from {0, 1, ..., N - 1}
  std::iota(std::begin(all), std::end(all), 0);
  // sort the 'subsys' vector for set_difference operation
  std::sort(std::begin(subsys), std::end(subsys));
  // compute the complement of 'subsys' using set_difference
  std::set_difference(std::begin(all), std::end(all), std::begin(subsys), std::end(subsys),
                      std::begin(subsys_bar));
  return subsys_bar;
}

/**
 * @brief convert a qubit density matrix to bloch vector
 * given a qubit density matrix 'A' this function calculates the corresponding bloch vector
 * representation. It returns a vector containing the three component (X, Y, and Z) of the bloch
 * vector. the input matrix 'A' must be a 2x2 qubit density matrix
 *
 * @tparam Derived the matrix type of the input 'A'
 * @param A the qubit density matrix to convert the bloch vector
 * @return the bloch vector representation of the qubit density matrix
 *
 * @example
 * Eigen::Matrix2cd rho;
 *
 * std::vector<double> blochVector = rho2bloch(rho);
 * std::cout << "bloch vector (X, Y, Z) " << blochVector[0] << ", "
 *           << blochVector[1] << ", " << blochVector[2] << std::endl;
 */
template <typename Derived>
std::vector<double> rho2bloch(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check qubit matrix
  if (!internal::check_qubit_matrix(rA))
    throw exception::NotQubitMatrix("clara::rho2bloch()");

  std::vector<double> result(3);
  cmat X(2, 2), Y(2, 2), Z(2, 2);
  X << 0, 1, 1, 0;
  Y << 0, -1_i, 1_i, 0;
  Z << 1, 0, 0, -1;

  // compute the X, Y, Z components of the bloch vector
  result[0] = std::real(trace(rA * X));
  result[1] = std::real(trace(rA * Y));
  result[2] = std::real(trace(rA * Z));
  return result;
}

/**
 * @brief compute the density matrix corresponding to the 3-dimensional real bloch vector 'r'
 * given an 3-dimensional real bloch vector 'r', this function calculates the corresponding qubit
 * density matrix. the density matrix represent quantum state of a qubit in a 2-dimensional
 * hilbert-space. the function returns a 2x2 complex matrix, which the qubit density matrix
 * associated with the bloch vector 'r'
 *
 * @param r the 3-dimensional real bloch vector representing a quantum state
 * @return the qubit density matrix corresponding to the bloch vector 'r'
 *
 * @example
 * std::vector<double> blochVector = {0.1, 0.2, 0.3};
 * cmat densityMatrix = bloch2rho(blochVector);
 */
inline cmat bloch2rho(const std::vector<double>& r) {
  // check if the input vector 'r' is 3-dimensional
  if (r.size() != 3)
    throw exception::CustomException("clara::bloch2rho", "r is not a 3-dimensional vector!");

  cmat X(2, 2), Y(2, 2), Z(2, 2), Id2(2, 2);

  // pauli matrices
  X << 0, 1, 1, 0;
  Y << 0, -1_i, 1_i, 0;
  Z << 1, 0, 0, -1;
  Id2 << 1, 0, 0, 1;

  // calculate the qubit density matrix using the bloch vector 'r'
  // rho = (Id2 + r[0] * X + r[1] * Y + r[2] * Z) / 2
  return (Id2 + r[0] * X + r[1] * Y + r[2] * Z) / 2.;
}

/**
 * @brief multi-partite qubit ket user-defined literal
 * this user-defined literal construct the multi partite qubit ket \f$|\mathrm{Bits}\rangle\f$
 * where 'Bits' is a binary repersentation of the quantum state. the literal allow constructing
 * qubit kets directly using binary string repersentation. for example 010_ket repersent the
 * qubit ket $|010\rangle$
 *
 * @tparam bits the binary representation of the quantum state, provided as a string of '0' and '1'
 * @return multi-partite qubit corresponding to the binary repersentation
 *
 * @example
 * // usage the user-defined literal to construct qubit ket |010⟩
 * ket qubitKet = 010_ket;
 */
template <char... Bits>
ket operator"" _ket() {
  constexpr idx n = sizeof...(Bits);
  constexpr char bits[n + 1] = {Bits..., '\0'};
  clara::ket q = clara::ket::Zero(std::pow(2, n));
  idx pos = 0;

  // check if the input binary repersentation bits contains only
  // '0' and '1'
  for (idx i = 0; i < n; ++i) {
    if (bits[i] != '0' && bits[i] != '1')
      throw exception::OutOfRange(R"xxx(clara::operator ""_ket())xxx");
  }

  // convert the binary representation to an index 'pos' in decimal (base-10) format
  pos = std::stoi(bits, nullptr, 2);
  q(pos) = 1;
  return q;
}

/**
 * @brief multi-partite qubit bra user-defined literal
 * this user defined literal construct the multi-partite qubit bra $\langle\mathrm{Bits}|\f$,
 * where 'Bits' is a binary representation of the quantum state. the literal allows constructing
 * qubit bra directly using binary string representation.
 *
 * @tparam Bits the binary representation of the quantum state, provided string of '0' and '1'
 * @return multi-partite qubit bra corresponding to the binary repersentation, as a complex dynamic
 *          row vector
 
 * @example
 * bra qubitBra = 010_bra;
 */
template <char... Bits>
bra operator"" _bra() {
  constexpr idx n = sizeof...(Bits);
  constexpr char bits[n + 1] = {Bits..., '\0'};
  clara::bra q = clara::ket::Zero(std::pow(2, n));
  idx pos = 2;

  // check valid multi-partite qubit state
  for (idx i = 0; i < n; ++i) {
    if (bits[i] != '0' && bits[i] != '1')
      throw exception::OutOfRange(R"xxx(qpp::operator "" _bra())xxx");
  }

  // convert the binary representation to an index 'pos' in decimal (base-10) format
  pos = std::stoi(bits, nullptr, 2);
  q(pos) = 1;

  return q;
}

/**
 * @brief multi-partite qubit projector user-defined literal
 * this user-defined literal coonstruct the multi-partite qubit projector
 * where 'Bits' is a binary repersentation of quantum state. the literal allows constructing
 * qubit projectors directly using a binary string representation.
 *
 * @tparam Bits the binary representation of the quantum state, provided as string of '0' and '1'
 * @return multi-partite qubit projector corresponding to the binary representation as a complex
 *                        dynamic matrix
 *
 * @example
 * cmat qubitProjector = 101_prj
 */
template <char... Bits>
cmat operator"" _prj() {
  constexpr idx n = sizeof...(Bits);
  constexpr char bits[n + 1] = {Bits..., '\0'};

  // check valid multi-partite qubit state
  for (idx i = 0; i < n; ++i) {
    if (bits[i] != '0' && bits[i] != '1')
      throw exception::OutOfRange(R"xxx(qpp::operator "" _prj())xxx");
  }

  // using the kronecker product to construct the qubit projector
  return kron(operator""_ket<Bits...>(), operator""_bra<Bits...>());
}

}  // namespace clara

#endif  // !FUNCTIONS_H_
