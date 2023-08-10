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
 * @brief generalized inner product between two quantum states
 * the function computes the inner product \f$\langle \phi_{\text{subsys}}|\psi\rangle\f$
 * between two quantum states \f$\phi\f$ and \f$\psi\f$ on a multi-partite quantum system.
 * the function allows for the subsystem to have different dimension. it returns the inner
 * product as a scalar or column vector over the remaining hilbert stapce after tracing out
 * the subsystem
 *
 * @tparam Derived the type of the input quantum state, can be any matrix expression that Eigen
 *          can handle
 * @param phi the first quantum state, repersented as a column vector
 * @param psi the second quantum state, repersented as column vector
 * @param subsys vector specifying the subsystem for which the inner product is computed
 * @param dims vector  specifying the dimension of the overall quantum system.
 * @return the inner product  \f$\langle \phi_{\text{subsys}}|\psi\rangle\f$ as scalar or
 *              column vector over the remaining hilbert space
 *
 * @example
 * // represent the qubit state |01>
 * cmat phi = 01_ket;
 * // represent the qubit state |10>
 * cmat psi = 10_ket;
 * // compute the inner product only for subsystem (the first qubit)
 * std::vector<idx> subsys = {0};
 * // the quantum system consist of two qubits
 * std::vector<idx> dims = {2, 2};
 * cmat result = ip(phi, psi, subsys, dims);
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

  // calculate various dimension and size required for the computation
  idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));
  idx D = static_cast<idx>(rpsi.rows());
  idx Dsubsys_bar = D / Dsubsys;

  idx N = dims.size();
  idx Nsubsys = subsys.size();
  idx Nsubsys_bar = N - Nsubsys;

  // pre-allocate arrays for multi-index conversion
  idx Cdims[maxn];
  idx Csubsys[maxn];
  idx Cdimssubsys[maxn];
  idx Csubsys_bar[maxn];
  idx Cdimssubsys_bar[maxn];

  // calculate complement of the subsystem indices
  std::vector<idx> subsys_bar = complement(subsys, N);
  std::copy(std::begin(subsys_bar), std::end(subsys_bar), std::begin(Cdimssubsys_bar));

  // copy dimension and subsystem indices to the pre-allocated arrays for multi-index conversion
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

  // lambda function for the worker threads to calculate the result in parallel
  auto worker = [=](idx b) noexcept -> typename Derived::Scalar {
    idx Cmidxrow[maxn];
    idx Cmidxrowsubsys[maxn];
    idx Cmidxcolsubsys_bar[maxn];

    internal::n2multiidx(b, Nsubsys_bar, Cdimssubsys_bar, Cmidxcolsubsys_bar);
    // write in the global row multi-index of the global row vector by inserting the
    // multi-index of the remaining hilbert space
    for (idx k = 0; k < Nsubsys_bar; ++k) {
      Cmidxrow[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
    }
    typename Derived::Scalar result = 0;

    // iterate over all possible multi-indices the subsystem and calculate the inner product
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

  // allocate the result vector and parallize the computation using OpenMP
  dyn_col_vect<typename Derived::Scalar> result(Dsubsys_bar);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
  for (idx m = 0; m < Dsubsys_bar; ++m)
    result(m) = worker(m);
  return result;
}

/**
 * @brief generalized inner product between two quantum states
 * this function computes the inner product \f$\langle \phi_{\text{subsys}}|\psi\rangle\f$
 * between two quantum state \f$\langle \phi_{\text{subsys}}|\psi\rangle\f$ on multi
 * partite quantum system. the function allows for the subsystem to have different dimension.
 * it returns the inner product as a scalar or dolumn vector over the remaining hilbert space
 * after tracing out the subsystem.
 *
 * @tparam Derive the type of the input quantum states. can be any matrix expression that Eigen
 *                can handle
 * @param phi first quantum state, represented as a column vector
 * @param psi second quantum state, represented as column vector
 * @param subsys  vector specifying the subsystem for which the inner product is computed
 * @param d the dimension of each subsystem. default value is 2 (qubit dimension)
 * @return  the inner product \f$\langle \phi_{\text{subsys}}|\psi\rangle\f$ as scalar or column
 *        vector over the remaining the hilbert space
 *
 * @example
 * // usage of the ip function to compute the inner product between two qubit states
 *
 * // represent the qubit state |01>
 * cmat phi = 01_ket;
 * // represent the qubit state |10>
 * cmat psi = 10_ket;
 *
 * // compute the inner product only for subsystem 0 (the first qubit)
 * idx d = 2;
 * cmat result = ip(phi, psi, subsys, d);
 */
template <typename Derived>
dyn_col_vect<typename Derived::Scalar> ip(const Eigen::MatrixBase<Derived>& phi,
                                          const Eigen::MatrixBase<Derived>& psi,
                                          const std::vector<idx>& subsys, idx d = 2) {
  const dyn_col_vect<typename Derived::Scalar>& rphi = phi.derived();
  const dyn_col_vect<typename Derived::Scalar>& rpsi = psi.derived();

  // check for zero-size input matrices
  if (!internal::check_nonzero_size(rpsi))
    throw exception::ZeroSize("clara::ip()");
  // check that dimension 'd' is valid
  if (d < 0)
    throw exception::DimsInvalid("clara::ip()");

  // calculate the number of subsystem 'N' for the given quantum state 'psi' and dimension 'd'
  idx N = internal::get_num_subsys(static_cast<idx>(rpsi.rows()), d);
  std::vector<idx> dims(N, d);
  // call the main 'ip' function with the computed 'dims'
  return ip(phi, psi, subsys, dims);
}

/**
 * @brief measures the quantum state 'A' using a set of kraus operators 'Ks'
 * this function perform quantum measurement on the input quantum state 'A' using a set of kraus
 * operators 'Ks' it returns a tuple containing the following information:
 *  - the index of the outcome of the measurement
 *  - the vector of outcome probabilities for each kraus operator
 *  - the vector of post measurement normalized states for each kraus operator
 *  - std::vector<cmat> the vector of post-measurement normalized states for each kraus operator
 *
 * @example
 * // usage of the measure function to perform a quantum measurement on a qubit state
 *
 * // represent the density matrix |01><01|
 * cmat rho = 01_prj
 * // measurement operators for outcomes |02><02| and |12><12|
 * std::tuple<idx, std::vector<double>, std::vector<cmat>> measurement = measure(rho, ks);
 * idx outcome_index = std::get<0>(measurement);
 * std::vector<double> outcome_probs = std::get<i1>(measurement);
 * std::vector<cmat> post_measurement_states = std::get<2>(measurement);
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const std::vector<cmat>& Ks) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check for zero-size input matrices
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");

  // check that all kraus operators are square and have the same dimension
  if (Ks.size() == 0)
    throw exception::ZeroSize("clara::measure()");
  if (!internal::check_square_mat(Ks[0]))
    throw exception::MatrixNotSquare("clara::measure()");
  if (Ks[0].rows() != rA.rows())
    throw exception::MatrixNotSquare("clara::measure()");
  for (auto&& it : Ks)
    if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
      throw exception::DimsNotEqual("clara::measure()");

  // initialize vector to store outcome probabilites and post-measurement states
  std::vector<double> prob(Ks.size());
  std::vector<cmat> outstates(Ks.size());

  // perform the measurement based on the type of quantum state 'A'
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

  // return the measurement result in a tuple
  return std::make_tuple(result, prob, outstates);
}

/**
 * @brief measures the quantum state 'A' in the orthonormal basis specified by the unitary matrix U
 * this function perform quantum measurement on the input quantum state 'A' using the orthonormal
 * basis sepcified by the unitary matrix 'U'. it returns a tuple containing the following
 * information:
 *  - the index of the outcome measurement
 *  - the vector of outcome probabilites for each measurement basis state
 *  - the vector of post-measurement normalized states for each measurement basis state
 *
 * @tparam Derived the type of the input quantum state. can be any matrix expression that eigen can
 * handle
 * @param A the input quantum state, represented as a square matrix (for density matrix) or column
 * vector (for pure state)
 * @param Ks the list of measurement basis state represented as a std::initializer_list of complex
 * matrices
 * @return a tuple containing:
 *  - 'idx': the index of the outcome of the measurement
 *  - 'std::vector<double>': the vector outcome probabilites for each measurement basis state
 *  - 'std::vector<cmat>': the vector post-measurement normalized states for each measurement basis
 * state
 *
 * @example
 * // usage of the measrue function to perform a quantum measurement on a qubit state
 *
 * // represent the density matrix |01> <01|
 * cmat rho = 01_prj;
 * // unitary matrix representing measurement bassis
 * std::tuple<idx, std::vector<double>, std::vector<cmat>> measurement = measure(rho, {02_prj,
 * 12_prj}); idx outcome_index = std::get<0>(measurement); std::vector<double> outcome_probs =
 * std::get<1>(measurement); std::vector<cmat> post_measurement_states = std::get<2>(measurement);
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
    const Eigen::MatrixBase<Derived>& A, const std::initializer_list<cmat>& Ks) {
  // call the main 'measure' function with the input measurement basis states as a std::vector
  return measure(A, std::vector<cmat>(Ks));
}

/**
 * @brief measure the quantum state 'A' in the orthonormal basis specified by the unitary matrix 'U'
 * this function perform quantum measurement on the input quantum state 'A' using the orthonormal
 * basis specified by the unitary matirx 'U', it returns a tuple containing the following
 * information
 *  - the index of the outcome of the measurement
 *  - the vector of outcome probabilites for each measurement basis state.
 *  - the vector of post-measurement normalized states for each measurement basis state
 *
 * @tparam Derived the type of the input quantum state. can be any matrix expression that Eigen can
 * handle
 * @param A the input quantum state, repersented as a square matrix (for density) or a column vector
 * @param U the unitary matri representing the measurement basis. it must have the same number of
 * rows 'A'
 * @return A tuple containing:
 *  - 'idx': the index of the outcome of the measurement
 *  - 'std::vector<double>': the vector of outcome probabilites for each measurement basis state
 *  - 'std::vector<cmat>': the vector of post-measurement normalized states for each measurement
 * basis state
 * @example
 * // usage of the measure function to perform a quantum measurement on a qubit state
 *
 * // represent the density matrix |01> <01|
 * cmat rho = 01_prj;
 * // unitary matrix representing the measurement basis
 * cmat U = cmat::Zero(4, 4);
 * // define the unitary matrix 'U' representing the measurement basis
 * // the columns of 'U' are the orthonormal measurement basis state
 * // 'U.col(0)' represents the first basis state, 'U.col(1)' represents the second basis state, and
 * so on std::tuple<idx, std::vector<double>, std::vector<cmat>> measurement = measure(rho, U); idx
 * outcome_index = std::get<0>(measurement); std::vector<double> outcome_probs =
 * std::get<1>(measurement); std::vector<cmat> post_measurement_states = std::get<2>(measurement);
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

  // calculate the measurement basis operators from the unitary matrix U
  std::vector<cmat> Ks(U.rows());
  for (idx i = 0; i < static_cast<idx>(U.rows()); ++i)
    Ks[i] = U.col(i) * adjoint(U.col(i));

  // perform the measurement using the measurement basis operators
  return measure(rA, Ks);
}

/**
 * @brief measure the part of the multi-partite state vector or density matrix 'A'
 *        corresponding to the susbsystem specified by 'subsys', using the set of kraus oeprators
 * 'Ks'
 *
 * this function perform a destructive measurement on the part of the multi-partite state 'A'
 * corresponding to the subsystem specified by 'subsys', using the set of kraus operator 'Ks'
 * it returns a tuple containing the following information:
 *  - the index of the outcome of the measurement
 *  - the vector of outcome probabilites for each measurement basis state.
 *  - the vector of post-measurement normalized states for each outcome, corresponding
 *    to each measure basis state
 *
 * @tparam Derived the type of the input quantum state. can be any matrix epression that eigen can
 * handle
 * @param A the input multi-partite state, repersented as a square matrix (for density matrices) or
 * column (for pure state)
 * @param Ks the set of kraus operators representing the measurement basis for the specified
 * subsystem
 * @param subsys the indices of the subsystem for which the measurement is performed
 * @param dims the dimension of the individiual subsystem, as a vector of integers
 * @return a tuple containing:
 *    - 'idx': the index of the outcome measurement
 *    - 'std::vector<double>': the vector of outcome probabilites for each measurement basis state
 *    - 'atd::vector<cmat>': the vector post-measurement normalized state for each outcome,
 *                          corresponding to each measurement basis state.
 *
 * @example
 * // usage of the measure function to perform a measurement on a multi-partite quantum state
 *
 * // represent the density matrix |01> <01|
 * cmat rho = 01_prj;
 * // set of kraus operators representing the measurement basis
 * std::vector<cmat> Ks
 * // define the kraus operators 'Ks' representing the measurement basis for the subsystem
 * // 'Ks[0]' represent the first measurement operator, 'Ks[1]' represent measurement operator
 * measure(rho, Ks, {0, 1}, {2, 2});
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
 * @brief measure the part of the multi-partite state vector or density matrix 'a'
 *        corresponding to the subsystem specified by 'subsys', using the set of kraus operators
 *        'ks'
 * this function performs a destructive measurement on the part of the multi-partite state 'a'
 * corresponding to the subystem specified by 'subsys', using the set of kraus operators 'ks'
 * it returns a tuple containing the following information:
 *  - the index of the outcome of the measurement
 *  - the vector of outcome probabilites for each measurement basis state
 *  - the vector of post-measurement normalized states for each outcome, corresponding to each
 *    measurement basis state.
 *
 * @tparam derived type of the input quantum state. can be any matrix expression that eigen can
 * handle
 * @param a the input multi-partite state, represented as a square matrix (for density matrices) or
 * a column vector (for pure state)
 * @param ks the set kraus operators representing the measurement basis for the specified subsystem
 * @param dims the dimension of the subsystem for which the measurement is performed
 * @param dims the dimension of the individiual subystem, as a vector of integers
 * @return a tuple containing
 *  - 'idx': the index of the outcome of the measurement
 *  - 'std::vector<double>': the vector of outcome probabilites for each measurement basis state
 *  - 'std::vector<cmat>': the vector of post-measurement normalized states for each outcome,
 * corresponding to each measurement basis state
 * @example
 * // usage of the measure function to perform a measurement on a multi-partite quantum state
 *
 * // represent the density matrix |01> <01|
 * cmat rho = 01_prj;
 * // define the kraus operator 'k0' and 'k1' representing the measurement basis for the subsystem
 * cmat k0 = 00_prj;
 * cmat k1 = 11_prj;
 * measure(rho, {k0, k1}, {0, 1}, {2, 2});
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
    const Eigen::MatrixBase<Derived>& A, const std::initializer_list<cmat>& Ks,
    const std::vector<idx>& subsys, const std::vector<idx>& dims) {
  return measure(A, std::vector<cmat>(Ks), subsys, dims);
}

/**
 * @brief measure the part of the multi-partite state vector or density matrix 'A'
 *        corresponding to the subsystem specified by 'subsys', using the set of kraus
 *        operators 'Ks'
 * this function performs a destructive measurement on the part of the multi-partite state 'A'
 * corresponding to the subystem specified by 'subsys', using the set of kraus operators 'Ks'
 * it return a tuple containing the following information:
 *  - the index of the outcome of the measurement
 *  - the vector of outcome probabilities for each measurement basis state
 *  - the vector of post-measurement normalized states for each outcome, corresponding to each
 * measurement basis state
 *
 * @tparam Derived the type of the input quantum state. can be any matrix expression that Eigen can
 * handle
 * @param A the input multi-partite state, repersented as a square matrix (for density matrices) or
 * a column vector
 * @param Ks the set of kraus operator representing the measurement basis for the specified
 * subsystem
 * @param subsys the indices of the subsystem for which the measurement is performed
 * @Param d the dimension of each individiual subsystem (default: 2 for qubits)
 * @return a tuple containing:
 *  - 'idx': the index of the outcome of the measurement
 *  - 'std::vector<double>': the vector of outcome probabilities afor each measurement basis state
 *  - 'std::vectpr<cmat>': the vector of post-measurement normalized states for each outcome,
 *                          corresponding to each measurement basis state
 *
 * @note the dimension 'd' must be greater than or equal to 2, the dimension of all kraus operator
 * 'Ks' must match the dimension of the specified 'subsys'. if 'd' is not provided, its default to 2
 *        which assumes qubits (2-dimensional subsystem)
 *
 * @example
 * // usage of the measure function to perform a measurement on a multi-partite quantum state
 * cmat rho = 01_prj;
 * cmat K0 = 00_prj;
 * cmat K1 = 11_prj;
 * measure(rho, {K0, K1}, {0, 1}, 2);
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const std::vector<cmat>& Ks,
                                                                const std::vector<idx>& subsys,
                                                                idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();
  // check if the input state matrix is not empty
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");
  // check if the dimension 'd' is valid
  if (d < 2)
    exception::DimsInvalid("clara::measure()");
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  // create a vector 'dims' containing the dimension 'd' for each individiual subsystem
  std::vector<idx> dims(N, d);
  // call the measure function with the calculated 'dims'
  return measure(rA, Ks, subsys, dims);
}

/**
 * @brief measure the part of the multi-partite state vector or density or density matrix 'A'
 *        corresponding to the subsystem sepcified by 'subsys', using the set of kraus operators
 * 'Ks'
 *
 * this function perform a destructive measurement on the part of the multi-partite state 'A'
 * corresponding to the subystem specified by 'subsys', using the set of kraus operatos 'Ks'
 * it returns a tuple containing following information
 *  - the index of the outcome of the measurement
 *  - the index of outcome probabilities for each measurement basis state
 *  - the vector of post-measurement normalized state for each outcome, corresponding
 *    to each measurement basis state
 *
 * @tparam Derived the type of the input quantum state. can be any matrix expression that Eigen can
 * handle
 * @param A the input multi-partite state, represented as a square matrix (for density matrices) or
 * a column vector
 * @param Ks the set of kraus operators representing the measurement basis for the specified
 * subsystem
 * @param d the dimension of each individiual subsystem (default: 2 for qubits)
 * @return a tuple containing
 *  - 'idx': the index of the outcome of the measurement
 *  - 'std::vector<double>': the vector of outcome probabilities for each measurement basis state
 *  - 'std::vector<cmat>': the vector of post-measurement normalized states for each outcome
 *                          corresponding to each measurement basis state
 * @note the dimension of all kraus 'Ks' must match the dimension of the specified 'subsys',
 *        if 'd' is not provided, it default to 2, which assumes qubits (2-dimensional subsystem)
 *
 * @example
 * // usage of the measure function to perform a measurement on a multi-partite quantum state
 *
 * // represent the density matrix |01><01|
 * cmat rho = 01_prj;
 * // represent the projector |00><00|
 * cmat K0 = 00_prj;
 * // repersente the projector |11><11|
 * cmat K1 = 11_prj;
 *
 * // measurement subsystem with indices 0 and 1, each of dimension 2 (qubits)
 * measure(rho, {K0, K1}, {0, 1}, 2);
 * cmat rho = 01_prj;
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
    const Eigen::MatrixBase<Derived>& A, const std::initializer_list<cmat>& Ks,
    const std::vector<idx>& subsys, idx d = 2) {
  return measure(A, std::vector<cmat>(Ks), subsys, d);
}

/**
 * @brief measure the part of the multi-partite state vector or density matrix 'A'
 *        corresponding to the subsystem specified by 'subsys' in the orthonormal basis
 *        rank-1 POV (positive operator-valued) measurement specified by matrix 'V'
 *
 * this function performs a destructive measurement on the part of the multi-partite state 'A'
 * corresponding to the subsystem specified by 'subsys', using the orthonormal basis rank-1 POV
 * measurement specified by the matrix 'V'. it returns a tuple containing the following information
 *  - the index of the outcome of measurement
 *  - the vector of outcome probabilites for each measurement basis state
 *  - the vector of post-measurement normalized states for each outcome, corresponding to each
 *    measurement basis state
 *
 * @tparam Derived the type of the input quantum state. can be any matrix expression that Eigen can
 * handle
 * @param A the input multi-partite state, repersented as a square matrix (for density matrix) or a
 * column column vector (for pure state)
 * @param V  the matrix repersenting the orthonormal basis rank-1 POV measurement
 * @param subsys the indices of the subsystem for hich the measurement is performed
 * @param dims the vector of dimension of each individiual subsystem in the multi-partite state
 * @return A tuple containing:
 *  - 'idx': the index of the outcome of the measurement
 *  - 'std::vector<double>': the vector of outcome probabilites fro each measurement basis state
 *  - 'std::vector<cmat>': the vector of post-measurement normalized state
 *
 * @note the dimension of matrix 'V' must match the dimension of the specified 'subsys'. the
 * dimension of all subsystem in 'dims' must match the respective dimension of the multi-partite
 * state. the measurement is destructive
 *
 * @example
 * // usage of the measure function to perform measurement on multi-partite quantum state
 *
 * // represent the density matrix |01><01|
 * cmat rho = 01_prj;
 * // define the orthonormal basis rank-1 POV matrix 'V' representing the measurement basis for the
 * // subsystem
 * cmat V = 01_ket;
 * // measure the subsystem with indices 0 and 1, each of dimension 2 (qubits)
 * measure(rhi, V, {0, 1}, {2, 2});
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const cmat& V,
                                                                const std::vector<idx>& subsys,
                                                                const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  // check if the input state amtrix is not empty
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");
  // check if the dimension 'dims' are valid
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::measure()");
  // check if the specified subsystem 'subsys' match the dimension dims
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::measure()");

  // calculate the dimension of the specified subsystem
  std::vector<idx> subsys_dims(subsys.size());
  for (idx i = 0; i < subsys.size(); ++i)
    subsys_dims[i] = dims[subsys[i]];

  // calculate the total dimension of the specified subsystem
  idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));

  // check if the POVM 'V' is not empty and matches the dimension of the specified
  // subsystem
  if (!internal::check_nonzero_size(V))
    throw exception::ZeroSize("clara::measure()");
  if (Dsubsys != static_cast<idx>(V.rows()))
    throw exception::DimsMismatchMatrix("clara::measure()");

  // calcualte the nu,ber of measurement outcomes 'M' based on the columns of the POVM 'V'
  idx M = static_cast<idx>(V.cols());

  // perform measurement based on the type of type input state 'A'
  if (internal::check_cvector(rA)) {
    // if 'A' is a pure state (column vector), perform measurement using the
    // provided POVM 'V'
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
 * @brief measure the part of the multi-partite state vector or density matrix 'A' in the
 * orthonormal basis specified by the rank-1 positive operator-valued measure (POVM) repersented in
 * 'V'
 *
 * this function perform a destructive measurement on the part of the multi-partite state 'A'
 * corresponding to the subsystem specified by 'subsys' using the rank-1 POVM 'V'
 * it returns a tuple containing the following information:
 *  - the inde of the outcomes of the measurement
 *  - the vector of outcome probabilites for each measurement basis state
 *  - the vector of post-measurement normalized state for each outcome, corresponding to each
 * measurement basis state
 *
 * @tparam Derived the type of the input quantum state, Can be any matrix expression that Eigen can
 * handle
 * @param A the input multi-partite state, represented as a square matrix (for density matrix) or
 * column column vector (for pure state)
 * @param subsys the indices of the subsystem for which the measurement is performed
 * @param dims the dimension of the subsystem in the multi-partite state
 * @return a tuple containing
 *  - 'idx': the index of the outcome of the measurement
 *  - 'std::vector<double>': the vector of outcome probabilites for each measurement basis state
 *  - 'std::vector<cmat>': the vector post-measurement normalized state for each outcome,
 * corresponding to each measurement basis state.
 *
 * @note the dimension of the POVM 'V' must match the dimension of the specified subsystem 'subsys'
 *      the input state matrix 'A' must not be empty
 *      the dimension 'dims' must be valid and match the dimension of 'A'
 *      the measurement outcome index is randomly sampled according to the outcome probabilites
 *      the function perform different measurement strategies based on whether 'A' is a pure state
 * or a density matrix (square matrix)
 *
 * @example
 * // usage of the measure function to perform a measurement on a multi-partite quantum state
 *
 * // repersent the density matrix |01> <01|
 * cmat rho = 01_prj;
 * // define  the rank-1 POVM 'V' representing the measurement basis for the subsystem
 * cmat V(2, 1);
 * V << 1, 1;
 * // measure the subystem with indices 0 and 1, each of dimension 2 (qubits)
 * measure(rho, V, {0, 1}, {2, 2});
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(const Eigen::MatrixBase<Derived>& A,
                                                                const cmat& V,
                                                                const std::vector<idx>& subsys,
                                                                idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EigenvaluesReturnType& rA = A.derived();

  // check if the input state matrix is not empty
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure()");

  // check if the dimension 'd' is valid
  if (d < 2)
    throw exception::DimsInvalid("clara::measure()");

  // calculate the number of subsystem 'N' based on the input state matrix and dimension 'd'
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  // create a vector containing dimension 'd' for all subsystem
  std::vector<idx> dims(N, d);
  // perform measurement on the specified subsystem using the orthonormal basis or POVM 'V'
  return measure(rA, V, subsys, dims);
}

/**
 * @brief sequentially measure the pasrt of the multi-partite state vector or density matrix 'A'
 *        in the computational basis
 *
 * this function perform sequential measurement on the subsystem specified by 'subsys' of the
 * multi-partite state 'A' using the computational basis. it returns a tuple containing the
 * following information
 *  - a vector containing the outcome of the measurement, oredered in increasing order with the
 * smallest index first
 *  - the overall outcome probability of the sequential measurement.
 *
 * @tparam Derived the type of the input quantum state. can be any matrix expression that Eigen can
 * handle
 * @param A the input multi-partite state, repersented as square matrix (for density matrix) or a
 * column vector (for pure state)
 * @param subsys the indices of the subsystem to be sequentially measured
 * @param dims the dimension of the subsystem
 * @return A tuple containing
 *  - 'std::vector<idx>': A vector containing the outcome result of the measurement
 *  - 'double': the overall outcome probability of the sequential measurement
 *  - 'cmat': the post-measurement normalized state obtained after the sequential measurement
 *
 * @note the input state matrix 'A' must not be empty. the dimension of the subsystem in 'dims'
 * souhld mactch the dimension of the corresponding subsystem specified 'subsys'. the function
 * perform sequential measurement in the computational basis, starting from the subsystem with the
 * largest index. after each measurement, the post-measurement normalized state is update for the
 * next measurement
 *
 * @example
 * // usage of the measure_seq function to perform sequential measurement on a multi-partite quantum
 * state
 *
 * // represent the density matrix |01> <01|
 * cmat rho = kron(01_ket, 01_bra);
 * // specify the dimension of the subystem
 * std::vector<idx> subsys = {0, 1}
 * std::vector<idx> dims = {2, 2};
 * measure_seq(rho, subsys, dims);
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

  // check if the input state matrix is square or a column vector, and if the dimension match the
  // input matrix
  if (!internal::check_square_mat(cA)) {
    if (!internal::check_dims_match_mat(dims, cA))
      throw exception::DimsMismatchMatrix("clara::measure_seq()");
  } else if (!internal::check_cvector(cA)) {
    if (!internal::check_dims_match_cvect(dims, cA))
      throw exception::DimsMismatchMatrix("clara::measure_seq()");
  } else
    throw exception::MatrixNotSquareNotCvector("clara::measure_seq");

  // check if the specified subsystem match the dimension
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::measure_seq()");
  std::vector<idx> result;
  double prob = 1;

  // sort the subsystem in decreasing order to start from the subsystem with the largest index
  std::sort(std::begin(subsys), std::end(subsys), std::greater<idx>{});

  // perform sequential measurement on the specified subsystem using the computational
  // basis
  while (subsys.size() > 0) {
    auto tmp = measure(cA, Gates::get_instance().Id(dims[subsys[0]]), {subsys[0], dims});
    result.push_back(std::get<0>(tmp));
    prob *= std::get<1>(tmp)[std::get<0>(tmp)];
    cA = std::get<2>(tmp)[std::get<0>(tmp)];
    dims.erase(std::next(std::begin(dims), subsys[0]));
    subsys.erase(std::begin(subsys));
  }

  std::reverse(std::begin(result), std::end(result));
  return std::make_tuple(result, prob, cA);
}

/**
 * @brief sequentially measures the specified subsystem of the multi-partite state
 *        vector or density matrix 'A' in the computational matrix
 *
 * this function perform sequential measurement on the subsystem specified by 'subsys' of the
 * multi-partite state 'A' using computational basis, it returns a tuple containing the
 * following:
 *  - a vector containing the outcome result of the measruements, ordered in increasing order with
 *    the ssmallest index first
 *  - the overall outcome probability of the sequential measurement
 *  - the post-measurement normalized state obtained after the sequential measruements
 *
 * @tparam Derived the type of the input quantum state, can be any matrix expression that Eigen can
 * hamdle
 * @param A the input multi-partite state, repersetend as a square matrix
 *        (for density  matrix) or a column vector (for pure states)
 * @param subsys the indices of the subsystem to be sequentially measured
 * @param d the dimension of the subsystem in the computational basis. default is 2
 * @return a tuple containing:
 *  - 'std::vector<idx>': A vector containing the outcome result of the measurement
 *  - 'double': the overall outcome probability of the sequential measurement
 *  - 'cmat': the post-measurement normalized state obtained after the sequential measurement
 *
 * @note the input state matrix 'A' must not be empty. the function perform sequential measruements
 *       in the computationalbsis, starting from the subsystem with the largest index. after each
 *       measurement, the post-measurement normalized state is update for the next measruement
 *       the parameter 'q' speficies the dimension of the subsystem in the computational basis.
 *       it must be greater than or equal to 2
 *
 * @example
 * // usage of the measure_seq function to perform sequential measurement on a multi-partite quantum
 * state
 *
 * // represente the density matrix |01> <01|
 * cmat rho = kron(01_ket, 01_bra);
 * // specify the subsystem to be sequentially measured
 * std::vector<idx> subsys = {0, 1}
 * // perform sequential measurement on the subsystem
 * measure_seq(rho, subsys);
 */
template <typename Derived>
std::tuple<std::vector<idx>, double, cmat> measure_seq(const Eigen::MatrixBase<Derived>& A,
                                                       std::vector<idx> subsys, idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  // check the input state matrix is not emtpy
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::measure_seq()");
  // check if the dimension 'd' is valid
  if (d < 2)
    throw exception::DimsInvalid("clara::measure_seq()");
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);

  // perform the sequential measurement using the measure_seq function with the calculated
  // dimension 'dims'
  std::vector<idx> dims(N, d);
  return measure_seq(rA, subsys, dims);
}

}  // namespace clara

#endif  // !INSTRUMENTS_H_
