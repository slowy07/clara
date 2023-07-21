#ifndef INSTRUMENTS_H_
#define INSTRUMENTS_H_

#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <random>
#include <tuple>
#include <vector>

#include "classFunction/exception.h"
#include "classFunction/gates.h"
#include "classFunction/random_devices.h"
#include "functions.h"
#include "internal/util.h"
#include "operations.h"
#include "types.h"

namespace clara {

/**
 * @brief generalized inner product
 * @return the inner product \f$\langle \phi_{subsys}|\psi\rangle\f$,
 * as scalar or column vector over the remaining hilbert space
 */
template <typename Derived>
dyn_col_vect<typename Derived::Scalar> ip(const Eigen::MatrixBase<Derived>& phi,
                                          const Eigen::MatrixBase<Derived>& psi,
                                          const std::vector<idx>& subsys,
                                          const std::vector<idx>& dims) {
  const dyn_col_vect<typename Derived::Scalar>& rphi = phi.derived();
  const dyn_col_vect<typename Derived::Scalar>& rpsi = psi.derived();

  if (!internal::check_nonzero_size(rphi))
    throw exception::ZeroSize("clara::ip()");
  if (!internal::check_nonzero_size(rpsi))
    throw exception::ZeroSize("clara::ip()");

  // check column vector
  if (!internal::check_cvector(rphi))
    throw exception::MatrixNotCvector("clara::ip()");
  // check column vector
  if (!internal::check_cvector(rpsi))
    throw exception::MatrixNotCvector("clara::ip()");
  // check dims is a valid dimenion vector
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::ip()");
  // check that subsys are valid
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::MatrixMismatchSubsys("clara::ip()");
  // check that dims match state vector psi
  if (!internal::check_dims_match_cvect(dims, rpsi))
    throw exception::DimsMismatchCvector("clara::ip()");

  std::vector<idx> subsys_dims(subsys.size());
  for (idx i = 0; i < subsys.size(); ++i)
    subsys_dims[i] = dims[subsys[i]];
  if (!internal::check_dims_match_cvect(subsys_dims, rphi))
    throw exception::DimsMismatchCvector("clara::ip()");

  idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));
  idx D = static_cast<idx>(rpsi.rows());
  idx Dsubsys_bar = D / Dsubsys;

  idx N = dims.size();
  idx Nsubsys = subsys.size();
  idx Nsubsys_bar = N - Nsubsys;

  idx Cdims[maxn];
  idx Csubsys[maxn];
  idx Cdimssubsys[maxn];
  idx Csubsys_bar[maxn];
  idx Cdimssubsys_bar[maxn];

  std::vector<idx> subsys_bar = complement(subsys, N);
  std::copy(std::begin(subsys_bar), std::end(subsys_bar), std::begin(Cdimssubsys_bar));

  for (idx i = 0; i < N; ++i) {
    Cdims[i] = dims[i];
  }
  for (idx i = 0; i < Nsubsys; ++i) {
    Csubsys[i] = subsys[i];
    Cdimssubsys[i] = dims[subsys[i]];
  }
  for (idx i = 0; i < Nsubsys_bar; ++i) {
    Cdimssubsys_bar[i] = dims[subsys_bar[i]];
  }

  auto worker = [=](idx b) noexcept -> typename Derived::Scalar {
    idx Cmidxrow[maxn];
    idx Cmidxrowsubsys[maxn];
    idx Cmidxcolsubsys_bar[maxn];

    internal::n2multiidx(b, Nsubsys_bar, Cdimssubsys_bar, Cmidxcolsubsys_bar);
    // write in the global row multi-index
    for (idx k = 0; k < Nsubsys_bar; ++k) {
      Cmidxrow[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
    }
    typename Derived::Scalar result = 0;
    for (idx a = 0; a < Dsubsys; ++a) {
      internal::n2multiidx(a, Nsubsys, Cdimssubsys, Cmidxrowsubsys);
      for (idx k = 0; k < Nsubsys; ++k) {
        Cmidxrow[Csubsys[k]] = Cmidxrowsubsys[k];
      }

      idx i = internal::multiidx2n(Cmidxrow, N, Cdims);
      result += std::conj(rphi(a)) * rpsi(i);
    }
    return result;
  };
  dyn_col_vect<typename Derived::Scalar> result(Dsubsys_bar);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
  for (idx m = 0; m < Dsubsys_bar; ++m)
    result(m) = worker(m);
  return result;
}

/**
 * @brief generalized inner product
 * @return the inner product \f$\langle \phi_{subsys}|\psi\rangle\f$, as a scalar
 * or column vector over the remaining hilbert space
 */
template <typename Derived>
dyn_col_vect<typename Derived::Scalar> ip(const Eigen::MatrixBase<Derived>& phi,
                                          const Eigen::MatrixBase<Derived>& psi,
                                          const std::vector<idx>& subsys, idx d = 2) {
  const dyn_col_vect<typename Derived::Scalar>& rphi = phi.derived();
  const dyn_col_vect<typename Derived::Scalar>& rpsi = psi.derived();

  if (!internal::check_nonzero_size(rpsi))
    throw exception::ZeroSize("clara::ip()");
  if (d < 0)
    throw exception::DimsInvalid("clara::ip()");

  idx N = internal::get_num_subsys(static_cast<idx>(rpsi.rows()), d);
  std::vector<idx> dims(N, d);
  return ip(phi, psi, subsys, dims);
}

/**
 * @brief measures the state A using the set of krays operator Ks
 * @return tuple of
 *   - result of the measurement
 *   - vector of outcome proabilities
 *   - vector of post-measureemnt normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const std::vector<cmat>& Ks) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");

  if (Ks.size() == 0)
    throw exception::ZeroSize("clara::measure()");
  if (!internal::check_square_mat(Ks[0]))
    throw exception::MatrixNotSquare("clara::measure()");
  if (Ks[0].rows() != rA.rows())
    throw exception::MatrixNotSquare("clara::measure()");
  for (auto&& it : Ks)
    if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
      throw exception::DimsNotEqual("clara::measure()");

  // proabilities
  std::vector<double> prob(Ks.size());
  // resulting states
  std::vector<cmat> outstates(Ks.size());

  if (internal::check_square_mat(rA)) {
    for (idx i = 0; i < Ks.size(); ++i) {
      outstates[i] = cmat::Zero(rA.rows(), rA.rows());
      cmat tmp = Ks[i] * rA * adjoint(Ks[i]);
      prob[i] = std::abs(trace(tmp));
      if (prob[i] > eps)
        outstates[i] = tmp / prob[i];
    }
  } else if (internal::check_cvector(rA)) {
    for (idx i = 0; i < Ks.size(); ++i) {
      outstates[i] = ket::Zero(rA.rows());
      ket tmp = Ks[i] * rA;
      prob[i] = std::pow(norm(tmp), 2);
      if (prob[i] > eps)
        outstates[i] = tmp / std::sqrt(prob[i]);
    }
  } else
    throw exception::MatrixNotSquareNotCvector("clara::measure()");

  // sample from the probability discrete_distribution
  std::discrete_distribution<idx> dd(std::begin(prob), std::end(prob));
  idx result = dd(RandomDevices::get_instance().get_prng());

  return std::make_tuple(result, prob, outstates);
}

/**
 * @brief measure the state A in the orthonormal basis specified
 * by the unitary matrix U
 * @return tuple:
 *   - result of the measurement
 *   - vector of outcome proabilities
 *   - vector of post measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
    const Eigen::MatrixBase<Derived>& A, const std::initializer_list<cmat>& Ks) {
  return measure(A, std::vector<cmat>(Ks));
}

/**
 * @brief measure the state A in the orthonormal basis specified by unitary matrix
 * @return tuple:
 *   - result of the measuremnt
 *   - vector of outcome probabilites
 *   - vector of post-measureemnt normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const cmat& U) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");

  // check the unitary basis matrix U
  if (!internal::check_nonzero_size(U))
    throw exception::ZeroSize("clara::measure()");
  if (!internal::check_square_mat(U))
    throw exception::MatrixNotSquare("clara::measure()");
  if (U.rows() != rA.rows())
    throw exception::DimsMismatchMatrix("clara::measure()");

  std::vector<cmat> Ks(U.rows());
  for (idx i = 0; i < static_cast<idx>(U.rows()); ++i)
    Ks[i] = U.col(i) * adjoint(U.col(i));

  return measure(rA, Ks);
}

/**
 * @brief measure the part subsys of the multi-partite state vector
 * or density matrix A using the set of kraus operator Ks
 * @note the dimension of all Ks must match the dimension of subsys
 * the measurement is destructive
 * @return tuple of:
 *   - result of the measurement
 *   - vector of outcome proabilities
 *   - vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const std::vector<cmat>& Ks,
                                                                const std::vector<idx>& subsys,
                                                                const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");
  // check that dimension is valid
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::measure()");
  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchMatrix("clara::mesaure()");

  // check subsy is valid w.r.t dims
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::measure()");
  std::vector<idx> subsys_dims(subsys.size());
  for (idx i = 0; i < subsys.size(); ++i)
    subsys_dims[i] = dims[subsys_dims[i]];

  idx D = prod(std::begin(dims), std::end(dims));
  idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));
  idx Dsubsys_bar = D / Dsubsys;

  // check kraus operator
  if (Ks.size() == 0)
    throw exception::ZeroSize("clara::measure()");
  if (!internal::check_square_mat(Ks[0]))
    throw exception::MatrixNotSquare("clara::measure()");
  if (Dsubsys != static_cast<idx>(Ks[0].rows()))
    throw exception::DimsMismatchMatrix("clara::measure()");
  for (auto&& it : Ks)
    if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
      throw exception::DimsNotEqual("clara::measure()");

  // proabilities
  std::vector<double> prob(Ks.size());
  // resulting states
  std::vector<cmat> outstates(Ks.size(), cmat::Zero(Dsubsys_bar, Dsubsys_bar));

  if (internal::check_square_mat(rA)) {
    for (idx i = 0; i < Ks.size(); ++i) {
      cmat tmp = apply(rA, Ks[i], subsys, dims);
      tmp = ptrace(tmp, subsys, dims);
      prob[i] = std::abs(trace(tmp));
      if (prob[i] > eps) {
        /* normalized output state
         * corresponding the measurement result i */
        outstates[i] = tmp / prob[i];
      }
    }
  } else if (internal::check_cvector(rA)) {
    for (idx i = 0; i < Ks.size(); ++i) {
      ket tmp = apply(rA, Ks[i], subsys, dims);
      prob[i] = std::pow(norm(tmp), 2);
      if (prob[i] > eps) {
        tmp /= std::sqrt(prob[i]);
        outstates[i] = ptrace(tmp, subsys, dims);
      }
    }
  } else
    throw exception::MatrixNotSquareNotCvector("clara::measure()");

  std::discrete_distribution<idx> dd(std::begin(prob), std::end(prob));
  idx result = dd(RandomDevices::get_instance().get_prng());
  return std::make_tuple(result, prob, outstates);
}

/**
 * @brief measure the part subsys of the multi-partite state vector or density matrix A
 * using the set of kraus operator
 * @note the dimensoin of all Ks must matech the dimenion of subsys the measurement is
 * destructive
 * @return tuple of:
 *   - result of measurement
 *   - vector of outcome proabilities
 *   - vector of post-measureemnt normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
    const Eigen::MatrixBase<Derived>& A, const std::initializer_list<cmat>& Ks,
    const std::vector<idx>& subsys, const std::vector<idx>& dims) {
  return measure(A, std::vector<cmat>(Ks), subsys, dims);
}

/**
 * @brief measure the part of subsy of the multi-partite state vector or density matrix A
 * using the set of kraus operator Ks
 * @note the dimension of all Ks must match the dimension of subsys
 * the measurement is destructive, the measured susbsystem are traced away
 * @return tuple of:
 *   - result of measurement
 *   - vector of outcome probabilites
 *   - vector post-measureemnt normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const std::vector<cmat>& Ks,
                                                                const std::vector<idx>& subsys,
                                                                idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");
  if (d < 2)
    exception::DimsInvalid("clara::measure()");
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return measure(rA, Ks, subsys, dims);
}

/**
 * @brief measure the part of subsys of the multi-partite state vector or density matrix A
 * using the set of kraus operator Ks
 * @note the dimension of all Ks must match the dimension of subsys
 * the measurement is destructive
 * @return tuple of:
 *   - result the measurement
 *   - vector outcome probabilites
 *   - vector of post measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
    const Eigen::MatrixBase<Derived>& A, const std::initializer_list<cmat>& Ks,
    const std::vector<idx>& subsys, idx d = 2) {
  return measure(A, std::vector<cmat>(Ks), subsys, d);
}

/**
 * @brief measures the part subsy of the multi-partite state vector or density matrix A
 * in the orthonormal basis rank-1 povm specified by matrix V
 * @note the dimension V must match the dimension of subsy. the measurement is destructive
 * @return tuple of:
 *   - result of measurement
 *   - vector outcome probabilites
 *   - vector of post-measurement nomralized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const cmat& V,
                                                                const std::vector<idx>& subsys,
                                                                const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::measure()");
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::measure()");
  std::vector<idx> subsys_dims(subsys.size());
  for (idx i = 0; i < subsys.size(); ++i)
    subsys_dims[i] = dims[subsys[i]];
  idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));

  if (!internal::check_nonzero_size(V))
    throw exception::ZeroSize("clara::measure()");
  if (Dsubsys != static_cast<idx>(V.rows()))
    throw exception::DimsMismatchMatrix("clara::measure()");

  idx M = static_cast<idx>(V.cols());
  if (internal::check_cvector(rA)) {
    const ket& rpsi = A.derived();
    if (!internal::check_dims_match_cvect(dims, rA))
      throw exception::DimsMismatchCvector("clara::measure()");

    std::vector<double> prob(M);
    std::vector<cmat> outstates(M);

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
    for (idx i = 0; i < M; ++i)
      outstates[i] = ip(static_cast<const ket&>(V.col(i)), rpsi, subsys, dims);
    for (idx i = 0; i < M; ++i) {
      double tmp = norm(outstates[i]);
      prob[i] = tmp * tmp;
      if (prob[i] > eps) {
        outstates[i] /= tmp;
      }
    }
    // sample from the probability distribution
    std::discrete_distribution<idx> dd(std::begin(prob), std::end(prob));
    idx result = dd(RandomDevices::get_instance().get_prng());
    return std::make_tuple(result, prob, outstates);
  } else if (internal::check_square_mat(rA)) {
    if (!internal::check_dims_match_mat(dims, rA))
      throw exception::DimsMismatchMatrix("clara::measure()");
    std::vector<cmat> Ks(M);
    for (idx i = 0; i < M; ++i)
      Ks[i] = V.col(i) * adjoint(V.col(i));
    return measure(rA, Ks, subsys, dims);
  }
  throw exception::MatrixNotSquareNotCvector("clara::measure()");
}

/**
 * @brief measure the part subsys of the multi-partite state vector or density matrix A
 * in the orthonormal basis or rank-1 POVM specified by the matrix V
 * @note the dimension V must match dimension of subsys. the measurement is destructive
 * @return tuple of:
 *   - result of the measurement
 *   - vector outcome proabilities
 *   - vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const cmat& V,
                                                                const std::vector<idx>& subsys,
                                                                idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EigenvaluesReturnType& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");
  if (d < 2)
    throw exception::DimsInvalid("clara::measure()");
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return measure(rA, V, subsys, dims);
}

/**
 * @brief sequentially measure the part subsys
 * of the multi-partite state vector or density matrix A in the computational basis
 * @return tuple of:
 *   - vector of outcome result of the measurement (ordered in increasing order
 *     with smallest index)
 *   - outcome probability
 *   - post-measurement normalized state
 */
template <typename Derived>
std::tuple<std::vector<idx>, double, cmat> measure_seq(const Eigen::MatrixBase<Derived>& A,
                                                       std::vector<idx> subsys,
                                                       std::vector<idx> dims) {
  dyn_mat<typename Derived::Scalar> cA = A.derived();

  if (!internal::check_nonzero_size(cA))
    throw exception::ZeroSize("clara::measure_seq()");
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::measure_seq()");
  if (!internal::check_square_mat(cA)) {
    if (!internal::check_dims_match_mat(dims, cA))
      throw exception::DimsMismatchMatrix("clara::measure_seq()");
  } else if (!internal::check_cvector(cA)) {
    if (!internal::check_dims_match_cvect(dims, cA))
      throw exception::DimsMismatchMatrix("clara::measure_seq()");
  } else
    throw exception::MatrixNotSquareNotCvector("clara::measure_seq");

  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::measure_seq()");
  std::vector<idx> result;
  double prob = 1;

  std::sort(std::begin(subsys), std::end(subsys), std::greater<idx>{});

  while (subsys.size() > 0) {
    auto tmp = measure(cA, Gates::get_instance().Id(dims[subsys[0]]), {subsys[0], dims});
    result.push_back(std::get<0>(tmp));
  }
}

/**
 * @brief sequentially measure the part subsys of multi-partite state vector
 * density matwrix A in the computational basis
 * @return tuple of:
 *   - vector of outcome result of the measurement (ordered increasing oredre with
 *     smallest index)
 *   - outcome probability
 *   - post-measurement normalized state
 */
template <typename Derived>
std::tuple<std::vector<idx>, double, cmat> measure_seq(const Eigen::MatrixBase<Derived>& A,
                                                       std::vector<idx> subsys, idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure_seq()");
  if (d < 2)
    throw exception::DimsInvalid("clara::measure_seq()");
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return measure_seq(rA, subsys, dims);
}

}  // namespace clara

#endif  // !INSTRUMENTS_H_
