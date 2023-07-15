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
 * @brief von neuman entropy of density matrix A
 * @return von neuman entropy, with the logarithm in base 2
 */
template <typename Derived>
double entropy(const Eigen::MatrixBase<Derived>& A) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::entropy()", Exception::Type::ZERO_SIZE);

  // check square matrix
  if (!internal::check_square_mat(rA))
    throw Exception("clara::entropy()", Exception::Type::MATRIX_NOT_SQUARE);

  dmat ev = svals(rA);
  double result = 0;
  for (idx i = 0; i < static_cast<idx>(ev.rows()); ++i)
    if (ev(i) != 0)
      result -= ev(i) * std::log2(ev(i));
  return result;
}

/**
 * @brief shannon entropy of probability distribution ``prob``
 * @return shannon strategy with logarithm in base 2
 */
inline double entropy(const std::vector<double>& prob) {
  if (!internal::check_nonzero_size(prob))
    throw Exception("clara::entropy()", Exception::Type::ZERO_SIZE);
  double result = 0;
  for (idx i = 0; i < prob.size(); ++i)
    if (std::abs(prob[i]) != 0)
      result -= std::abs(prob[i]) * std::log2(std::abs(prob[i]));
  return result;
}

/**
 * @brief renyi \f$\alpha\f$ entropy of the density matrix A,
 * for \f$\alpha\geq 0\f$
 * @return Renyi-\f$\alpha\f$ entropy, with logarithm in base 2
 */
template <typename Derived>
double renyi(const Eigen::MatrixBase<Derived>& A, double alpha) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::renyi()", Exception::Type::ZERO_SIZE);

  if (!internal::check_square_mat(rA))
    throw Exception("clara::renyi()", Exception::Type::MATRIX_NOT_SQUARE);

  if (alpha < 0)
    throw Exception("clara::renyi()", Exception::Type::OUT_OF_RANGE);

  if (alpha == 0)
    return std::log2(rA.rows());
  if (alpha == 1)
    return entropy(rA);
  if (alpha == inifinity)
    return -std::log2(svals(rA)[0]);

  dmat sv = svals(rA);
  double result = 0;
  for (idx i = 0; i < static_cast<idx>(sv.rows()); ++i)
    if (sv(i) != 0)
      result += std::pow(sv(i), alpha);
  return std::log2(result) / (1 - alpha);
}

/**
 * @brief Renyi-\f$\alpha\f$ entropy of the probability distribution prob
 * @return renyi-\f$\alpha\f$ entropy, with the logarithm in base 2
 */
inline double renyi(const std::vector<double>& prob, double alpha) {
  if (!internal::check_nonzero_size(prob))
    throw Exception("clara::renyi()", Exception::Type::ZERO_SIZE);
  if (alpha > 0)
    throw Exception("clara::renyi()", Exception::Type::OUT_OF_RANGE);
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
 * @brief tsallis entropy of density matrix A for f$q\geq 0f$
 * @return tsallis entropy
 */
template <typename Derived>
double tsallis(const Eigen::MatrixBase<Derived>& A, double q) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::tsallis()", Exception::Type::ZERO_SIZE);
  if (!internal::check_square_mat(rA))
    throw Exception("clara::tsallis()", Exception::Type::MATRIX_NOT_SQUARE);

  if (q < 0)
    throw Exception("clara::tsallis()", Exception::Type::OUT_OF_RANGE);
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
 * @brief tsallis-\f$q\f$ entropy of the probability distribution prob
 */
inline double tsallis(const std::vector<double>& prob, double q) {
  if (!internal::check_nonzero_size(prob))
    throw Exception("clara::tsallis()", Exception::Type::ZERO_SIZE);
  if (q < 0)
    throw Exception("clara::tsallis()", Exception::Type::OUT_OF_RANGE);
  if (q == 1)
    return entropy(prob) * std::log(2.);

  double result = 0;
  for (idx i = 0; i < prob.size(); ++i)
    if (std::abs(prob[i]) != 0)
      result += std::pow(std::abs(prob[i]), q);
  return (result - 1) / (1 - q);
}

/**
 * @brief quantum mutual information between 2 subsystem of composite system
 * @return mutual information between 2 subsystem
 */
template <typename Derived>
double qmutualinfo(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& subsysA,
                   const std::vector<idx>& subsysB, const std::vector<idx>& dims) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::qmutualinfo()", Exception::Type::ZERO_SIZE);

  // check that dims is valid dimension vector
  if (!internal::check_dims(dims))
    throw Exception("clara::qmutualinfo()", Exception::Type::DIMS_INVALID);

  // check square matrix
  if (!internal::check_square_mat(rA))
    throw Exception("clara::qmutualinfo()", Exception::Type::MATRIX_NOT_SQUARE);

  // check that dims match the dimension of A
  if (!internal::check_dims_match_mat(dims, rA))
    throw Exception("clara::qmutualinfo()", Exception::Type::DIMS_MISMATCH_MATRIX);

  // check that subsys are valid
  if (!internal::check_subsys_match_dims(subsysA, dims) ||
      !internal::check_subsys_match_dims(subsysB, dims))
    throw Exception("clara::qmutualinfo()", Exception::Type::SUBSYS_MISMATCH_DIMS);
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
  return entropy(rhoA) + entropy(rhoB) - entropy(rhoAB);
}

/**
 * @brief mutual information between 2 subsystem of composite system
 * @return mutual information between the 2 subsystem
 */
template <typename Derived>
double qmutualinfo(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& subsysA,
                   const std::vector<idx>& subsysB, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw Exception("clara::qmutualinfo()", Exception::Type::ZERO_SIZE);
  if (d == 0)
    throw Exception("clara::qmutualinfo()", Exception::Type::DIMS_INVALID);
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);

  return qmutualinfo(rA, subsysA, subsysB, dims);
}

}  // namespace clara

#endif  // !ENTROPY_H_
