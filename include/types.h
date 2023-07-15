#ifndef TYPES_H_
#define TYPES_H_

#include <complex>
#include <cstddef>
#include <eigen3/Eigen/Dense>

namespace clara {

using idx = std::size_t;
using bigint = long long int;
using ubigint = unsigned long long int;
using cplx = std::complex<double>;
using ket = Eigen::VectorXcd;
using bra = Eigen::RowVectorXcd;
using cmat = Eigen::MatrixXcd;
using dmat = Eigen::MatrixXd;

template <typename Scalar>
using dyn_mat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Scalar>
using dyn_col_vect = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

template <typename Scalar>
using dyn_row_vect = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

}  // namespace clara

#endif  // !TYPES_H_
