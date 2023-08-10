#ifndef ENTROPY_H_
#define ENTROPY_H_

#include <cmath>
#include <limits>

#include "classFunction/exception.h"
#include "classFunction/gates.h"
#include "constants.h"
#include "functions.h"
#include "internal/util.h"
#include "types.h"

namespace clara {

/**
 * @brief calculate the von neuman entropy of the density matrix 'A'
 * @param A eigen matrix or matrix expression representing the density matrix
 * @return double the von neuman entropy, with the logarithm in base 2
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::MatrixNotSquare thrown if 'A' is not a square matrix
 *
 * @example
 * Eigen::matrixXd densityMatrix;
 *
 * // calculate the von neumann entropy for the density matrix
 * double entropyValue = entropy(densityMatrix);
 */
template <typename Derived>
double entropy(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::entropy()");

  // check square matrix
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::entropy()");

  // calculate the singular values 'A'
  dmat ev = svals(rA);
  double result = 0;
  for (idx i = 0; i < static_cast<idx>(ev.rows()); ++i) {
    // calculate the von neuman entropy term for each non-zero
    // singular value
    if (ev(i) != 0)
      result -= ev(i) * std::log2(ev(i));
  }
  return result;
}

/**
 * @brief calculate shanon entropy of a probability distribution
 *
 * @param prob the probability distribution as a vector of doubles
 * @return double the shanon entropy, with the logarithm in base 2
 *
 * @exception exception::ZeroSize thrown if the input probability distribution has zero size
 *
 * @note the function assumes that the probability distribution is a valid
 *        discrete probability distribution. the probability should be non-negative and sum up to 1
 *
 * @example
 * std::vector<double> probabilities = {0.3, 0.2, 0.5};
 *
 * double entropyValue = entropy(probabilities);
 */
inline double entropy(const std::vector<double>& prob) {
  if (!internal::check_nonzero_size(prob))
    throw exception::ZeroSize("clara::entropy()");
  double result = 0;
  for (idx i = 0; i < prob.size(); ++i) {
    // calculate the shanon entropy term for each non-zero
    // probability
    if (std::abs(prob[i]) != 0)
      result -= std::abs(prob[i]) * std::log2(std::abs(prob[i]));
  }
  return result;
}

/**
 * @brief calculate the Renyi \alpha entropy of a density matrix A
 * @tparam Derived type of the matrix-like object
 * @param A the density matrix for which to calculate the Renyi \alpha entropy
 * @return Renyi entropy of the density matrix A, with logarithm in base 2
 *
 * @throw ZeroSize if the matrix A has zero size
 * @throw MatrixNotSquare if the matrix A is not square
 * @throw OutOfRange if alpha is negative
 *
 * NOTE: The Renyi α entropy is defined for α >= 0 and is calculated as:
 *       H_α(A) = (1 / (1 - α)) * log2(Σ(λ_i^α)),
 *       where λ_i are the eigenvalues of A.
 *       Special cases:
 *       - α = 0: H_0(A) = log2(dimension of A)
 *       - α = 1: H_1(A) = von Neumann entropy of A
 *       - α = ∞: H_∞(A) = -log2(λ_max), where λ_max is the largest eigenvalue of A.
 */
template <typename Derived>
double renyi(const Eigen::MatrixBase<Derived>& A, double alpha) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::renyi()");

  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::renyi()");

  if (alpha < 0)
    throw exception::OutOfRange("clara::renyi()");

  // special cases: α = 0, 1, ∞
  if (alpha == 0)
    return std::log2(rA.rows());
  if (alpha == 1)
    return entropy(rA);
  if (alpha == inifinity)
    return -std::log2(svals(rA)[0]);

  // general case
  dmat sv = svals(rA);
  double result = 0;
  for (idx i = 0; i < static_cast<idx>(sv.rows()); ++i)
    if (sv(i) != 0)
      result += std::pow(sv(i), alpha);
  return std::log2(result) / (1 - alpha);
}

/**
 * @brief calculate the renyi \alpha-entropy of a density matrix 'A'
 *
 * the renyi \alpha-entropy is a generalization of the shannon entropy and is defined as follows
 * - when \alpha > 0 and \alpha ≠ 1, the renyi \alpha-entropy given by:
 *   \f$ H_{\alpha}(A) = \frac{1}{1-\alpha} \log_2 \left( \sum_i \lambda_i^\alpha \right) \f$A
 *   where \delta_i are the singular values of the density matrix 'A'.
 * - when \alpha = 0, the renyi \alpha-entropy is given by
 *   f$ H_{0}(A) = \log_2 ( \text{rank}(A) ) \f$, where rank(A) is the rank of 'A'
 * - when \delta = 1, the renyi \delta-entropy is equivalent to the shannon entropy, given by
 *   \f$ H_{1}(A) = - \sum_i \lambda_i \log_2(\lambda_i) \f$
 *   where \delta_i are the singular values of the entropy density matrix 'A'
 * - whena \delta approaches infinity, the renyi \delta-entropy approaches
 *   \f$ H_{\infty}(A) = - \log_2(\lambda_{\text{max}}) \f$, where λ_{\text{max}}
 *   is the largest singluar singular value of 'A'
 *
 * @tparam derived the matrix expression type
 * @param A eigen matrix or matrix expression representing the density matrix
 * @param alpha the parameter \delta for renyi \delta-entropy calculation. should be a >= 0
 * @return double the renyi \delta-entropy with the logarithm in base 2
 *
 * @exception exception::ZeroSize thrown if 'A' has zero size
 * @exception exception::MatrixNotSquare thrown if 'A' is not a square matrix
 * @exception exception::OutOfRange thrown if 'alpha' is less than 0
 *
 * @example
 * Eigen::matrixXd densityMatrix;
 * double alpha = 2.0;
 *
 * // calculate the renyi \alpha-entropy for the density matrix with \alpha = 2
 * double renyiEntropy = renyi(densityMatrix, alpha);
 */
inline double renyi(const std::vector<double>& prob, double alpha) {
  if (!internal::check_nonzero_size(prob))
    throw exception::ZeroSize("clara::renyi()");
  if (alpha > 0)
    throw exception::OutOfRange("clara::renyi()");
  if (alpha == 0)
    return std::log2(prob.size());
  if (alpha == 1)
    return entropy(prob);

  if (alpha == inifinity) {
    double max = 0;
    for (idx i = 0; i < prob.size(); ++i)
      if (std::abs(prob[i]) > max)
        max = std::abs(prob[i]);
    return -std::log2(max);
  }
  double result = 0;
  for (idx i = 0; i < prob.size(); ++i)
    if (std::abs(prob[i]) != 0)
      result += std::pow(std::abs(prob[i]), alpha);
  return std::log2(result) / (1 - alpha);
}

/**
 * @brief calculate the Tsallis q-entropy density matrix 'A'
 *
 * Tsallis q-entropy is generalization of the shannon entropy and is defined as follows:
 * - when q > 0 and != 1, Tsallis q-entropy is given
 *   \f$ T_q(A) = \frac{{1 - \sum_i \lambda_i^q}}{{q - 1}} \f$
 *   where \delta_i are the singular values of the density matrix 'A'
 * - when q = 1, the Tsallis q-entropy is equivalent the von neuman entropy, given by
 *   \f$ T_1(A) = S(A) \cdot \log(2) \f$
 *   where S(A) is the shannon entropy of 'A'
 *
 * @tparam derived the matrix expression type
 * @param A eigen matrix or matrix expression represnting the density matrix
 * @param q the parameter q for the Tsallis q-entropy calculation should be q >= 0
 * @return double the Tsallis q-entropy
 *
 * @exception exception::ZeroSize if 'A' has zero size
 * @exception exception::MatrixNotSquare thrown if 'A' is not square matrix
 * @exception exception::OutOfRange thrown if 'q' is less than 0
 *
 * @example
 * Eigen::matrixXd densityMatrix;
 * double q = 2.0;
 *
 * // calculate the Tsallis q-entropy for the density matrix with q = 2
 * double densityMatrix = tsallis(densityMatrix, q);
 *
 */
template <typename Derived>
double tsallis(const Eigen::MatrixBase<Derived>& A, double q) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::tsallis()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::tsallis()");

  if (q < 0)
    throw exception::OutOfRange("clara::tsallis()");
  if (q == 1)
    return entropy(rA) * std::log(2.);

  dmat ev = svals(rA);
  double result = 0;
  for (idx i = 0; i < static_cast<idx>(ev.rows()); ++i)
    if (ev(i) != 0)
      result += std::pow(ev(i), q);
  return (result - 1) / (1 - q);
}

/**
 * @brief tsallis-\f$q\f$ entropy of the probability distribution prob  for q >= 0
 *
 * the Tsallis q-entropy is generalization of the shannon entropy and is defined as follows
 * - when q > 0 and q != 1, the Tsallis q-entropy is given by
 *   \f$ T_q(\text{prob}) = \frac{{\sum_i |p_i|^q - 1}}{{1 - q}} \f$
 *   where the p_i are the elements of the probability distribution 'prob'
 * - when q = 1, the Tsallis q-entropy is equivalent of the shannon entropy, given by
 *   \f$ T_1(\text{prob}) = S(\text{prob}) \f$
 *   where $(prob) is the shannon entropy of the probability distribution
 * - when q = 0, the Tsallis q-entropy is equivalent to the logarithm of the number of non-zero
 * elements in 'prob', i.e \f$ T_0(\text{prob}) = \log_2(\text{number of non-zero elements in prob})
 * \f$
 *
 * @param prob the vector of double representing the probability distribution
 * @param q the parameter q for the Tsallis q-entropy calculation. should be q >= 0
 * @return double the Tsallis q-entropy of the given probability distribution, with the logarithm in
 * base 2
 *
 * @exception exception::ZeroSize thrown if 'prob' has zero size
 * @exception exception::OutOfRange thrown if 'q' is less than 0
 *
 * @example
 * std::vector<double> probabilityDistribution = {0.2, 0.3, 0.5};
 * double q = 0.5;
 *
 * // calculate the Tsallis q-entropy for the probability distribution with q = 0.5
 * double tsallisEntropy = tsallis(probabilityDistribution, q);
 */
inline double tsallis(const std::vector<double>& prob, double q) {
  // check if the probability distribution has non-zero size
  if (!internal::check_nonzero_size(prob))
    throw exception::ZeroSize("clara::tsallis()");
  if (q < 0)
    throw exception::OutOfRange("clara::tsallis()");
  // check if q is valid (q >= 0)
  if (q == 1)
    return entropy(prob) * std::log(2.);

  // calculate the tsallis q-entropy using the provided formula
  double result = 0;
  for (idx i = 0; i < prob.size(); ++i) {
    // calculate the term for each non-zero probability value
    if (std::abs(prob[i]) != 0)
      result += std::pow(std::abs(prob[i]), q);
  }
  // calculate the Tsallis q-entropy and return the result
  return (result - 1) / (1 - q);
}

/**
 * @brief calculate the quantum mutual information between two subsystem of a composite
 *        system
 * quantum mutual information measures the mutual information two subsystem A and B of composite
 * quantum system. it is defined as the different of the sum of the entropies of the reduced
 * dmatrices of subsystem A and B and the entropy of the reduced density matrix of their AB: \f$
 * I(A:B) = S(\rho_A) + S(\rho_B) - S(\rho_{AB}) \f$
 *
 * @param A the input density matrix representing the composite quantum system
 * @param subsysA indices of the subsystem A
 * @param subsysB indices of the subsystem B
 * @param dims vector of dimension of the subsystem
 * @return double the quantum mutual information between subsystem A and B
 *
 * @exception exception::ZeroSize thrown if the input density matrix 'A' has zero size
 * @exception exception::DimsInvalid thrown if the input dimension vector 'dims' is invalid
 * @exception exception::MatrixNotSquare thrown if the input density matrix 'A' is not square
 * @exception exception::DimsMismatchMatrix thrown if the input dimension vector `dims` does not
 * match dimension of 'A'
 * @exception exception::SubsysMismatchdims thrown if the input subsystem 'subsyA' and 'subsysB' do
 * not match dimension 'A'
 *
 * @example
 * Eigen::matrixXd densityMatrix = ...;
 * std::vector<idx> subsystemA = {0, 1};
 * std::vector<idx> subsystemB = {2, 3};
 * std::Vector<idx> dimensions = {2, 2, 2, 2};
 *
 * // calculate the quantum mutual information between subsystem A and B
 * double mutualInfo = qmutualinfo(densityMatrix, subsystemA, subsystemB, dimensions);
 */
template <typename Derived>
double qmutualinfo(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& subsysA,
                   const std::vector<idx>& subsysB, const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::qmutualinfo()");

  // check that dims is valid dimension vector
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::qmutualinfo()");

  // check square matrix
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::qmutualinfo()");

  // check that dims match the dimension of A
  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchMatrix("clara::qmutualinfo()");

  // check that subsys are valid
  if (!internal::check_subsys_match_dims(subsysA, dims) ||
      !internal::check_subsys_match_dims(subsysB, dims))
    throw exception::SubsysMismatchdims("clara::qmutualinfo()");
  std::vector<idx> full_system(dims.size());
  std::iota(std::begin(full_system), std::end(full_system), 0);

  // sorted input subsystem
  std::vector<idx> subsysAsorted{subsysA};
  std::vector<idx> subsysBsorted{subsysB};

  std::sort(std::begin(subsysAsorted), std::end(subsysAsorted));
  std::sort(std::begin(subsysBsorted), std::end(subsysBsorted));

  // construct the union of A and B
  std::vector<idx> subsysAB;
  std::set_union(std::begin(subsysAsorted), std::end(subsysAsorted), std::begin(subsysBsorted),
                 std::end(subsysBsorted), std::back_inserter(subsysAB));

  std::vector<idx> subsysA_bar = complement(subsysA, dims.size());
  std::vector<idx> subsysB_bar = complement(subsysB, dims.size());
  ;
  std::vector<idx> subsysAB_bar = complement(subsysAB, dims.size());

  cmat rhoA = ptrace(rA, subsysA_bar, dims);
  cmat rhoB = ptrace(rA, subsysB_bar, dims);
  cmat rhoAB = ptrace(rA, subsysAB_bar, dims);

  // calculate the quantum mutual information using provided formula
  return entropy(rhoA) + entropy(rhoB) - entropy(rhoAB);
}

/**
 * @brief calculate the mutual information between two subystem of a compute system
 *
 * mutual information measures the amount of information shared between two subsystem
 * A and B of a composite quantum system. it is defined as the difference of the sum of
 * the entropies of the reduced density matrices of subsystem A and B and the entropy of the
 * reduced density matrix of their union AB: \f$ I(A:B) = S(\rho_A) + S(\rho_B) - S(\rho_{AB}) \f$
 *
 * @param A the input density matrix representing the composite quantum system
 * @param subsyA indices of the subsystem A
 * @param subsyB indices of the subsystem B
 * @param d dimension of the subsystem (default is 2)
 * @return double the mutual information between subsystem A and B
 *
 * @example
 * Eigen::matrixXd densityMatrix = ...;
 * std::vector<idx> subsystemA = {0, 1}
 * std::vector<idx> subsystemB = {0, 2}
 *
 * // calculate the mutual information between subsystem A and B
 * double mutualInfo = qmutualinfo(densityMatrix, subsystemA, subsystemB, dimension);
 */
template <typename Derived>
double qmutualinfo(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& subsysA,
                   const std::vector<idx>& subsysB, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::qmutualinfo()");
  if (d < 2)
    throw exception::DimsInvalid("clara::qmutualinfo()");
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);

  return qmutualinfo(rA, subsysA, subsysB, dims);
}

}  // namespace clara

#endif  // !ENTROPY_H_
