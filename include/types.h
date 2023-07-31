#ifndef TYPES_H_
#define TYPES_H_

#include <complex>
#include <cstddef>
#include <eigen3/Eigen/Dense>

namespace clara {

/**
 * @brief alias for the size type for indexing and dimension
 */
using idx = std::size_t;

/**
 * @brief  alias for signed integer type with a larger size, suitable for lager integers
 */
using bigint = long long int;

/**
 * @brief alias for an unsigned integer type with a larger size, suitable for large positive
 * integers
 */
using ubigint = unsigned long long int;

/**
 * @brief alias for complex number type using double precision
 */
using cplx = std::complex<double>;

/**
 * @brief alias for a column vector of complex nubmers using Eigen library
 */
using ket = Eigen::VectorXcd;

/**
 * @brief alias for a row vector of complex numbers using Eigen library
 */
using bra = Eigen::RowVectorXcd;

/**
 * @brief alias for matrix of complex number using Eigen library
 */
using cmat = Eigen::MatrixXcd;

/**
 * @brief alias template for dynamically sized matrices with elements of tyype 'Scalar'
 */
using dmat = Eigen::MatrixXd;

/**
 * @brief alias template for dynamically sized column vectors with elements of type 'Scalar'
 */
template <typename Scalar>
using dyn_mat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * @brief alias template for dynamically sized column vectors with elements of type 'Scalar'
 */
template <typename Scalar>
using dyn_col_vect = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

/**
 * @brief alias template for dynamically sized row vectors with elements of type 'Scalar'
 */
template <typename Scalar>
using dyn_row_vect = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

}  // namespace clara

#endif  // !TYPES_H_
