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
    throw exception::ZeroSize("clara::schmidtcoeffs()");
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::schmidtcoeffs()");
  if (!internal::check_cvector(rA))
    throw exception::MatrixNotSquareNotCvector("clara::schmidtcoeffs()");
  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchCvector("clara::schmidtcoeffs()");

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
    throw exception::ZeroSize("clara::schmidtcoeffs()");
  if (d < 2)
    throw exception::DimsInvalid("clara::schmidtcoeffs()");
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
    throw exception::ZeroSize("clara::schmidtA()");
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::schmidtA()");
  if (!internal::check_cvector(rA))
    throw exception::MatrixNotCvector("clara::schmidtA()");
  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchCvector("clara::schmidtA()");

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
    throw exception::ZeroSize("clara::schmidtA()");
  if (d < 2)
    throw exception::DimsInvalid("clara::schmidtA()");

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
    throw exception::ZeroSize("clara::schmidtB()");
  // check bi-partite
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::schmidtB()");
  if (!internal::check_cvector(rA))
    throw exception::MatrixNotCvector("clara::schmidtB()");

  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchCvector("clara::schmidtB()");

  // by default returns U_B*
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
    throw exception::ZeroSize("clara::schmidtB()");
  if (d < 0)
    throw exception::DimsInvalid("clara::schmidtB()");

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  // local dimension vector
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
    throw exception::ZeroSize("clara::schmidtprobs()");
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::schmidtprobs()");

  if (!internal::check_cvector(rA))
    throw exception::MatrixNotCvector("clara::schmidtprobs()");

  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchCvector("clara::schmidtprobs()");

  std::vector<double> result;
  dyn_col_vect<double> scf = schmidtcoeffs(rA, dims);
  for (idx i = 0; i < static_cast<idx>(scf.rows()); ++i)
    result.push_back(std::pow(scf(i), 2));
  return result;
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
    throw exception::ZeroSize("clara::schmidtprobs()");
  if (d < 2)
    throw exception::DimsInvalid("clara:schmidtprobs()");

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
    throw exception::ZeroSize("clara::entanglement()");
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::entanglement()");
  if (!internal::check_cvector(rA))
    throw exception::MatrixNotCvector("clara::entanglement()");
  // check matching dimensions
  if (!internal::check_dims_match_cvect(dims, rA))
    throw exception::DimsMismatchCvector("clara::entanglement()");

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
    throw exception::ZeroSize("clara::entanglement()");
  if (d < 2)
    throw exception::DimsInvalid("clara::entanglement()");

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
    throw exception::ZeroSize("clara::gconcurrence()");
  if (!internal::check_cvector(rA))
    throw exception::MatrixNotCvector("clara::gconcurrence()");
  idx d = internal::get_dim_subsystem(static_cast<idx>(rA.rows()), 2);

  if (d * d != static_cast<idx>(rA.rows()))
    throw exception::DimsNotEqual("clara::gconcurrence()");
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
    throw exception::ZeroSize("clara::negativity()");
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::negativity()");
  // check square matrix vector
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::negativity()");

  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchMatrix("clara::negativity()");
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
    throw exception::ZeroSize("clara::negativity()");
  if (d < 2)
    throw exception::DimsInvalid("clara::negativity()");

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
    throw exception::ZeroSize("clara::lognegativity()");
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::lognegativity()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::lognegativity()");
  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchMatrix("clara::lognegativity()");
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
    throw exception::ZeroSize("clara::lognegativity()");
  if (d < 0)
    throw exception::DimsInvalid("clara::lognegativity()");

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
    throw exception::ZeroSize("clara::concurrence()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::concurrence()");
  if (rA.rows() != 4)
    throw exception::NotQubitSubsys("clara::concurrence()");
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
