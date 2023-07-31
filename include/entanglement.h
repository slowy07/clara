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
 * @brief calculate the schmidt coefficients of bipartite pure state 'A'
 * @note the sum of square of schmidt coefficients equal 1
 * @tparam derived the matrix expression type
 * @param A the eigen matrix of matrix expression representing the bipartitepure state
 * @param dims vector containing the dimension of the tow subsystem
 *              the size of the vector should be 2, and dims[0] should be the dimension
 *              of the first dimension of the second subsystem
 * @return dyn_col_vect<double> A column vector containing the schmidt coefficients
 *                              of pure state 'A'
 *
 * @exception exception::ZeroSize throw if 'A' has zero size
 * @exception exception::NotBipartite thrown if 'dims' does not contain excatly 2 dimension
 * @exception exception::DimsMismatchCvector thrown if dimension in 'dims' do not match
 *                                          the dimension of 'A'
 * @example
 * Eigen::MatrixXd bipartiteState(3, 4);
 *
 * // define the dimension of the two subsystem
 * std::vector<idx> dimensions = {3, 4};
 *
 * // calculat the schmidt coefficients of the bipartite pure state
 * dyn_col_vect<double> schmidCoeffs = schmidtcoeffs(bipartiteState, dimensions);
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

  // calculate the singlar value of the transposed and reshaped matrix
  return svals(transpose(reshape(rA, dims[1], dims[0])));
}

/**
 * @brief calculate the schmidt coefficients of a bipartite pure state 'A'
 * @note the sum of the square of schmidt coefficients equals 1
 * @tparam derived the matrix expression type
 * @param A the eigen matrix or matrix expression representing the bipartite pure state
 * @param d the dimension of each subsystem. default value is 2
 * @return dyn_col_vect<double> A column vector containing schmidt coefficients of the
 *                              pure state 'A'
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::DimsInvalid thrown if 'd' is less than 2
 * @exception exception::NotBipartite thrown if the dimension 'd' result in a non bipartite
 *
 * @example
 * Eigen::MatrixXd bipartiteState(4, 4);
 *
 * // calculate the schmidt coefficients of the bipartite pure state with a default dimensio
 * // of 3
 * dyn_col_vect<double> schmidtCoeffs = schmidtcoeffs(bipartiteState);
 *
 * // calculate the schmidt coefficients of the bipartite purestate
 * // with a specified dimension of 3
 * dyn_col_vect<double> schmidCoeffs3 = schmidtcoeffs(bipartiteState, 3);
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
 * @brief calculate the schmidt basis on alices's side for a bipartite pure state 'A'
 * @tparam derived the matrix expression type
 * @param dims A vector containing the dimension of the tow subsystem
 *            the size of the vector should be 2, and dims[0] shuold be the dimension of
 *            the first subsystem, and dims[1] should be the dimension of the second subsystem
 * @return cmat the unitary matrix 'U' whose column represent the schmidt basis vectors on alice
 * slide
 *
 * @exception exception::ZeroSize thrown if 'A' is zero size
 * @exception exception::NotBipartite thrown is 'dims' does not conain excatly 2 dimension
 * @exception exception::MatrixNotCvector thrown if 'A' is not column vector
 * @exception exception::DimsMismatchCvector thrown if the dimension in 'dims' do not
 *                                        match dimension of 'A'
 * @example
 * Eigen::MatrixXd bipartiteState(3, 4);
 *
 * // define the dimension of the two subystem
 * std::vector<idx> dimension = {3, 4};
 *
 * // calculate the schmidt basis on alice side for the bipartite pure state
 * cmat schmidtMatrix = schmidtA(bipartiteState, dimension);
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
 * @brief calculate the schmidt basis on Alice's side for a bipartite pure 'A'
 * @tparam derived the matrix expression type
 * @param A eigen matrix or matrix expression representing the bipartite pure state
 * @param d the dimension fo each subystem default is 2
 * @return cmat the unitary matrix `u` whose column represent the schmidt basis vector on alice's
 * side.
 *
 * @exception exception::ZeroSize thrown is 'A' has zero size
 * @exception exception::DimsInvalid throw if 'd' is less than 2
 * @exception exception::NotBipartite thrown if the dimension 'd' result in non-bipartite state
 *
 * @example
 * Eigen::MatrixXd bipartiteState(4, 4);
 *
 * // calculate the schmidt basis on Alice's side with default dimension of 2 for each subystem
 * cmat schmidtMatrix = schmidtA(bipartiteState);
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
 * @brief calculate the schmidt basis on bob's side for a bipartite pure state 'A'
 * @tparam the eigen matrix or matrix expression represnting the bipartite pure state
 * @param dims a vector containing the dimension of the two subsystem
 *              the size of the vector should be 2, and dims[0] should be the dimension of
 *              the first subsystem, and dims[1] should be the dimension of the second
 *              subsystem.
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::NotBipartite thrown if 'dims' does not contain excatly 2 dimension
 * @exception exception::MatrixNotCvector thrown if 'A' is not column vector
 * @exception exception::DimsMismatchCvector thrown if the dimensions in 'dims' do not match of 'A'
 *
 * @example
 * Eigen::MatrixXd bipartiteState(3, 4);
 *
 * // define the dimension of the two subystem
 * std::vector<idx> dimensions = {3, 4};
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
 * @brief calculate the schmidt basis bob's side for a bipartite pure state 'A'
 * @tparam derived the matrix expression type
 * @param d the dimension of each subsystem expression representing the bipartite pure state
 * @param d the dimension of each subsystem. default value value is 2
 *
 * @exception exception::ZeroSize thrown if 'A' has hero size
 * @exception exception::DimsInvalid thrown if 'd' is less than 0
 * @exception exception::NotBipartite thrown if the dimension 'd' result in non bipartite-state
 *
 * @example
 * Eigen::MatrixXd bipartiteState(4, 4);
 *
 * // calculate the schmidt basis on bob's side with a specified dimension of
 * // of 3 for each subsystem
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
 * @brief calculate the squared schmidt coefficients (schmidt probabilites) for a bipartite
 *        state 'A'
 * @tparam derived the marix expression type
 * @param A the eigen matrix or matrix expression representing the bipartite pure state
 * @param dims A vector containing the dimension of the two subsystem
 *            the size of the vector should be 2, and dims[0] should be the dimension of the
 *            first subsystem, and dims[1] should be the dimension of the second subsystem
 * @return std::vector<double> A vector containing the squared schmidt coefficients
 *
 * @example
 * Eigen::MatrixXd bipartiteState(3, 4);
 *
 * // Define the dimension of the two subsystem
 * std::vector<idx> dimension = {3, 4};
 *
 * // caculate the squared schmidt coefficients for the bipartite pure state
 * std::vector<double> schmidtProbs = schmidtprobs(bipartiteState, dimensions);
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
 * @brief calculate the square schmidt probabilites for a bipartite pure state 'A'
 *
 * @tparam derived matrix expression type
 * @param A the eigen matrix or matrix expression representing the bipartite pure state
 * @param d the dimension of each subsystem. default value is 2
 * @return std::vector<double> A vector containgin the squared schmidt probabilites
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::DimsInvalid thrown if 'd' is less than 2
 *
 * @example
 * Eigen::MatrixXd bipartiteState(4, 4);
 *
 * // calculate the squared schmidt probabilites with a default dimension of 2 for each
 * // subsystem
 * std::vector<double> schmidtProbs = schmidtprobs(bipartiteState);
 *
 * // calculate the squared schmidt probabilites with a specified dimension fo 3 for each
 * // subystem
 * std::vector<double> schmidtProbs3 = schmidtprobs(bipartiteState, 3);
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
 * @brief calculate the logarithmic negativity of a bipartite mixed state 'A'
 * @tparam derived the matrix expression type
 * @param A the eigen matrix or matrix expression representing the bipartite mixed state
 * @param d the dimension of each subsystem. default value is 2
 * @return double the logarithmic negativity with the logarithm in base 2
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::DimsInvalid thrown if 'd' less than 0
 *
 * @example
 * Eigen::MatrixXd bipartiteMixedState(4, 4);
 *
 * // calculate the logarithmic negativity for the bipartite mixed state with a default
 * double logNegativityValue = lognegativity(bipartiteMixedState);
 *
 * // calculate the logarithmic negativity for the bipartite mixed state with a specific
 * // subystem
 * double logNegativityValue3 = lognegativity(bipartiteMixedState, 3);
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
 * @brief calculate the entanglement of a bipartite pure state 'A'
 * @tparam Derived the matrix expression type
 * @param A the eigen matrix or matrix expression representing the bipartite pure state
 * @param d the diension of each subsystem, default is 2
 * @return double the entanglement with logarithm in base 2
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::DimsInvalid thrown if 'd' is less than 2
 *
 * @example
 * Eigen::MatrixXd bipartiteState(4, 4);
 *
 * // calculate the entanglement with a defualt dimension of 2 for each subsystem
 * double entanglementValue = entanglement(bipartiteState);
 *
 * // calculate the entanglement with a specified dimension of
 * double entanglementValue = entanglement(bipartiteState, 3);
 *
 * // calculate the entanglement with a specified dimension of 3 for each subsystem
 * double entanglementValue3 = entanglement(bipartiteState, 3);
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
 * @brief calculated the G-concurrence of a bipartite pure state 'A'
 * @tparam derived the matrix expression type
 * @param A the eigen matrix or matrixx expression representing the bipartite pure state.
 * @return double the G-concurrence
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::MatrixNotCvector thrown if 'A' is not a column
 * @exception exception::DimsNotEqual thrown if the number
 *
 * @example
 * Eigen::MatrixXd bipartiteState(4, 4);
 *
 * // calculate the G-concurrence for the bipartite pure state
 * double gConccurenceValue = gconcurrence(bipartiteState);
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
 * @brief calculate the negativity of a bipartite mixed state 'A'.
 * @param A the eigen matrix or matrix expression representing the bipartite mixed state
 * @param dims the dimension of the subystem as a vector [dim_subys1, dims_subsys2]
 * @return double the negativity of the bipartite mixed state
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::NotBipartite thrown if the dims vector does not have exactly two elements
 * @exception exception::MatrixNotSquare thrown if 'A' is not a square matrix
 * @exception exception::DimsMismatchMatrix thrown if the dimension
 *                                specified by 'dims' do not match the matrix
 *
 * @example
 * Eigen::MatrixXd bipartiteMixedState(4, 4);
 *
 * std::vector<idx> dims = {2, 2};
 *
 * // calculate the negativity for the bipartite mixed state
 * double negativityValue = negativity(bipartiteMixedState, dims);
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
 * @brief calculate the negativity of a bipartite mixed state 'A'
 * @tparam Derived the matrix expression representing the bipartite
 * @param A the eigen matrix or matrix expression representing the bipartite mixed state
 * @param d the dimension of each subsystem. default value is 2
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::DimsInvalid thrown if 'd' is less than 2
 *
 * @example
 * Eigen::MatrixXd bipartiteMixedState(4, 4);
 *
 * // calculate the negativity for the bipartite mixed state
 * double negativityValue = negativity(bipartiteMixedState);
 *
 * // calculate the negativity for the bipartite mixed state with specified dimension of 3 for each
 * subystem double negativityValue3 = negativity(bipartiteMixedState, 3);
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
 * @brief calculate the logarithmic negativity of a bipartite mixed state 'A'
 * @tparam Derived the matrix expression type
 * @param A the eigen matrix expression representing the bipartite mixed state.
 * @param dims the dimension of the two subsystem as vector [dims_subsys1, dims_subsys2]
 *
 * @exception exception::ZeroSize thrown if 'A' has zeros
 * @exception exception::NotBipartite thrown 'dims' vector does not have exactly two elements
 * @exception exception::MatrixNotSquare thrown if 'A' is not a square matrix
 * @exception exception::DimsMismatchMatrix thrown if the dimension specified  by 'dims'
 *                                   not match matrix 'A'
 * @example
 * Eigen::MatrixXd bipartiteMixedState(4, 4);
 *
 * std::vector<idx> dims = {2, 2};
 *
 * // calculate the logarithmic negativity for the bipartite mixed state
 * double logNegativityValue = lognegativity(bipartiteMixedState, dims);
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
 * @brief calculate the logarithmic negativity of a bipartite mixed state 'A'
 * @param A the Eigen matrix or matrix expression representing the bipartite mixed state
 * @param d the dimension of each subsystem, default value is 2.
 * @return double the logarithmic negativity with the logarithmic in base 2
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::DimsInvalid thrown if 'd' is less than 0
 *
 * @example
 * Eigen::MatrixXd bipartiteMixedState(4, 4);
 *
 * // calculate the logarithmic negativityfor the bipartite mixed state
 * // with mixed state with default dimension of 2 for each
 * double logNegativityValue = lognegativity
 *
 * // calculate the logarithmic negativity for the bipartite mixed state with specified
 * // dimension of 3 for each subsystem
 * double logNegativityValue3 = lognegativity(bipartiteMixedState, 3);
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
 * @brief calculate the wotters concurrence of a bipartite mixed state 'A'
 * @tparam derived the matrix expression type
 * @param A the eigen matrix or matrix expression representing the bipartite qubit mixed state
 * @return double the wooters concurrence
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::MatrixNotSquare thrown if 'A' is not a square matrixqq
 * @exception exception::NotQubitSubsys thrown if 'A' does not represent a qubit subsystem
 *
 * @example
 * Eigen::Matrix2cd bipartiteMixedState;
 *
 * // calculate the wooters concurrence for the bipartite qubit  mixed state
 * double concurrenceValue = concurrence(bipartiteMixedState);
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
