#ifndef ENTANGLEMENT_H_
#define ENTANGLEMENT_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <ios>
#include <iterator>

#include "classFunction/exception.h"
#include "classFunction/gates.h"
#include "functions.h"
#include "internal/util.h"
#include "operations.h"
#include "types.h"
namespace clara {
/**
 * @brief schidmit coefficients of bi-partite pure state A
 * @note the sum square of the schidmit coefficients equals 1
 */
template <typename Derived>
dyn_col_vect<double> schmidcoeffs(const Eigen::MatrixBase<Derived>& A,
                                  const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::schmidcoeffs()", Exception::Type::ZERO_SIZE);
  if (dims.size() != 2)
    throw Exception("clara::schmidcoeffs()", Exception::Type::NOT_BIPARTITE);
  if (!internal::check_cvector(rA))
    throw Exception("clara::schmidcoeffs()", Exception::Type::MATRIX_NOT_CVECTOR);
  if (!internal::check_dims_match_mat(dims, rA))
    throw Exception("clara::schmidcoeffs()", Exception::Type::DIMS_MISMATCH_MATRIX);

  return svals(transpose(reshape(rA, dims[1], dims[0])));
}

/**
 * @brief schmidt coefficients of the bi-partite pure state A
 * @note the sum of the squares of schmidt coefficients equals
 * @return schmidt coefficients of A, as a real dynamic column vector
 */
template <typename Derived>
dyn_col_vect<double> schmidtcoeffs(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(A))
    throw Exception("clara::schmidtcoeffs()", Exception::Type::ZERO_SIZE);
  if (d == 0)
    throw Exception("clara::schmidtcoeffs()", Exception::Type::DIMS_INVALID);
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);

  return schmidtcoeffs(A, dims);
}

/**
 * @brief schmidt basis on alice slide
 * @return unitary matrix \f$ U \f$ whose columns represent
 * the schmidt basis vectors on alice slide
 */
template <typename Derived>
cmat schmidtA(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::schmidtU()", Exception::Type::ZERO_SIZE);
  if (dims.size() != 2)
    throw Exception("clara::schmidtU()", Exception::Type::NOT_BIPARTITE);
  if (!internal::check_cvector(rA))
    throw Exception("clara::schmidtU()", Exception::Type::MATRIX_NOT_CVECTOR);
  if (!internal::check_dims_match_mat(dims, rA))
    throw Exception("clara::schmidtU()", Exception::Type::DIMS_MISMATCH_MATRIX);

  return svdU(transpose(reshape(rA, dims[1], dims[0])));
}

/**
 * @brief schmidt basis on alice side
 * @return unitary matrix \f$ U \f$ whose columns represent
 * the schmidt basis vectors on alice side
 */
template <typename Derived>
cmat schmidtA(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(A))
    throw Exception("clara::schmidtA()", Exception::Type::ZERO_SIZE);
  if (d == 0)
    throw Exception("clara::schmidtA()", Exception::Type::DIMS_INVALID);

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return schmidtA(A, dims);
}

/**
 * @brief schmidt basis on bob side
 * @return unitary matrix ``V`` whose columns represent
 * the schmidt basis vectors on bob side
 */
template <typename Derived>
cmat schmidtB(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::schmidtV()", Exception::Type::ZERO_SIZE);
  // check bi-partite
  if (dims.size() != 2)
    throw Exception("clara::schmidtV()", Exception::Type::NOT_BIPARTITE);
  if (!internal::check_cvector(rA))
    throw Exception("clara::schmidtV()", Exception::Type::MATRIX_NOT_CVECTOR);

  if (!internal::check_dims_match_mat(dims, rA))
    throw Exception("clara::schmidtV()", Exception::Type::DIMS_MISMATCH_MATRIX);

  return svdV(transpose(reshape(conjugate(rA), dims[1], dims[0])));
}

/**
 * @brief schmidt basis on bob side
 * @return unitary matrix V whose columns repersent the schmidt basis
 * vectors on bob side
 */
template <typename Derived>
cmat schmidtB(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(A))
    throw Exception("clara::schmidtB()", Exception::Type::ZERO_SIZE);
  if (d == 0)
    throw Exception("clara::schmidtB()", Exception::Type::DIMS_INVALID);

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return schmidtB(A, dims);
}

/**
 * @brief schmidt probabilities of the bi-partite pure state A
 * define as the square of the schmidt coefficients
 * the sum the schmidt probabilities equals 1.
 * @return real vector consistring of the schmidt probabilities.
 */
template <typename Derived>
std::vector<double> schmidtprobs(const Eigen::MatrixBase<Derived>& A,
                                 const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::schmidtprobs()", Exception::Type::ZERO_SIZE);
  if (dims.size() != 2)
    throw Exception("clara::schmidtprobs()", Exception::Type::NOT_BIPARTITE);

  if (!internal::check_cvector(rA))
    throw Exception("clara::schmidtprobs()", Exception::Type::MATRIX_NOT_CVECTOR);

  if (!internal::check_dims_match_mat(dims, rA))
    throw Exception("clara::schmidtprobs()", Exception::Type::DIMS_MISMATCH_MATRIX);

  std::vector<double> result;
  dyn_col_vect<double> scf = schmidtcoeffs(rA, dims);
  for (idx i = 0; i < static_cast<idx>(scf.rows()); ++i)
    result.push_back(std::pow(scf(i), 2));
}

/**
 * @brief schmidt probabilities of the bi-partite pure state A
 * defined as the square of the schmidt coefficients probabilities equals 1
 * @return real vector consiting of the schmidt probabilities of A
 */
template <typename Derived>
std::vector<double> schmidtprobs(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar> rA = A.derived();

  if (!internal::check_nonzero_size(A))
    throw Exception("clara::schmidtprobs()", Exception::Type::ZERO_SIZE);
  if (d == 0)
    throw Exception("clara::schmidtprobs()", Exception::Type::DIMS_INVALID);

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return schmidtprobs(A, dims);
}

/**
 * @brief entanglement of the bi-partite pure state A
 * defined as the von neuman entropy of the reduced density matrix
 * of one of the subsystem
 * @return entanglement with logarithm in base
 */
template <typename Derived>
double entanglement(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::entanglement()", Exception::Type::ZERO_SIZE);
  if (dims.size() != 2)
    throw Exception("clara::entanglement()", Exception::Type::NOT_BIPARTITE);
  if (!internal::check_cvector(rA))
    throw Exception("clara::entanglement()", Exception::Type::MATRIX_NOT_CVECTOR);
  // check matching dimensions
  if (!internal::check_dims_match_mat(dims, rA))
    throw Exception("clara::entanglement()", Exception::Type::DIMS_MISMATCH_MATRIX);
  return entropy(schmidtprobs(rA, dims));
}

/**
 * @brief entanglement of bi-partite pure state A
 * defined as the von-neumann entropy of the reduced density matrix
 * of one the subsystem
 * @return entanglement with logarithm in base 2
 */
template <typename Derived>
double entanglement(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(A))
    throw Exception("clara::entanglement()", Exception::Type::ZERO_SIZE);
  if (d == 0)
    throw Exception("clara::entanglement()", Exception::Type::DIMS_INVALID);

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return entanglement(A, dims);
}

/**
 * @brief G-concurrence of the bi-partite pure state A
 * @return G-concurrence
 */
template <typename Derived>
double gconcurrence(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::gconcurrence()", Exception::Type::ZERO_SIZE);
  if (!internal::check_cvector(rA))
    throw Exception("clara::gconcurrence()", Exception::Type::MATRIX_NOT_CVECTOR);
  idx d = internal::get_dim_subsystem(static_cast<idx>(rA.rows()), 2);

  if (d * d != static_cast<idx>(rA.rows()))
    throw Exception("clara::gconcurrence()", Exception::Type::DIMS_NOT_EQUAL);
  return d * std::abs(std::exp(2. / d * logdet(reshape(rA, d, d))));
}

/**
 * @brief negativity of the bi-partite mixed state A
 * @return negativity
 */
template <typename Derived>
double negativity(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::negativity()", Exception::Type::ZERO_SIZE);
  if (dims.size() != 2)
    throw Exception("clara::negativity()", Exception::Type::NOT_BIPARTITE);
  // check square matrix vector
  if (!internal::check_square_mat(rA))
    throw Exception("clara::negativity()", Exception::Type::MATRIX_NOT_SQUARE);

  if (!internal::check_dims_match_mat(dims, rA))
    throw Exception("clara::negativity()", Exception::Type::DIMS_MISMATCH_MATRIX);
  return (schatten(ptranspose(rA, {0}, dims), 1) - 1.) / 2.;
}

/**
 * @brief negativity of the bi-partite mixed state A
 * @return negativity
 */
template <typename Derived>
double negetivity(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(A))
    throw Exception("clara::negativity()", Exception::Type::ZERO_SIZE);
  if (d == 0)
    throw Exception("clara::negativity()", Exception::Type::DIMS_INVALID);

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return negativity(A, dims);
}

/**
 * @brief logarithmic negativity bi-partite mixed state A
 * @return logarithmic negativity with the logarithm in base 2
 */
template <typename Derived>
double lognegativity(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::lognegativity()", Exception::Type::ZERO_SIZE);
  if (dims.size() != 2)
    throw Exception("clara::lognegativity()", Exception::Type::NOT_BIPARTITE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::lognegativity()", Exception::Type::MATRIX_NOT_SQUARE);
  if (!internal::check_dims_match_mat(dims, rA))
    throw Exception("clara::lognegativity()", Exception::Type::DIMS_MISMATCH_MATRIX);
  return std::log2(2 * negativity(rA, dims) + 1);
}

/**
 * @brief logarithmic negativity of the bo-partite mixed state A
 * @return logarithmic negativity, with the logarithmic in base 2
 */
template <typename Derived>
double lognegativity(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(A))
    throw Exception("clara::lognegativity()", Exception::Type::ZERO_SIZE);
  if (d == 0)
    throw Exception("clara::lognegativity()", Exception::Type::DIMS_INVALID);

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return lognegativity(A, dims);
}

/**
 * @brief wootters concurrence of the bi-partite qubit mixed state A
 * @return wotters concurrence
 */
template <typename Derived>
double concurrence(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::concurrence()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::concurrence()", Exception::Type::MATRIX_NOT_SQUARE);
  if (rA.rows() != 4)
    throw Exception("clara::concurrence()", Exception::Type::NOT_QUBIT_SUBSYS);
  cmat sigmaY = Gates::get_instance().Y;
  dyn_col_vect<double> lambdas =
      evals(rA * kron(sigmaY, sigmaY) * conjugate(rA) * kron(sigmaY, sigmaY)).real();
  std::vector<double> lambdas_sorted(lambdas.data(), lambdas.data() + lambdas.size());

  std::sort(std::begin(lambdas_sorted), std::end(lambdas_sorted), std::greater<double>());
  std::transform(std::begin(lambdas_sorted), std::end(lambdas_sorted), std::begin(lambdas_sorted),
                 [](double elem) { return std::sqrt(std::abs(elem)); });

  return std::max(0.,
                  lambdas_sorted[0] - lambdas_sorted[1] - lambdas_sorted[2] - lambdas_sorted[3]);
}

}  // namespace clara

#endif  // !ENTANGLEMENT_H_
