// Copyright (c) 2023 arfy slowy
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef CLASSFUNCTION_GATES_H_
#define CLASSFUNCTION_GATES_H_

#include <complex.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>
#include <optional>
#include <vector>

#include "../functions.h"
#include "../internal/classFunction/singleton.h"
#include "../internal/util.h"
#include "../number_theory.h"
#include "exception.h"

namespace clara {

/**
 * @class clara::Gates
 * @brief Singleton
 *
 * this class is singleton and provides a collection of quantum gates and operation commonly
 * used in quantum computing. it contains gates such as Pauli-X, Pauli-Y, pauli-Z, Hadamard,
 * S, T, CNOT, CZ, SWAP, TOFFOLI, Fredkin, and more. the class also inclueds methods to generate
 * rotation gates, generalized Z gates, Fourier transform gates, and Identity gates for qudits.
 * additionally, it provides methods to construct multi-parity multiple controlled gates. the gates
 * class ensure that gate are initialized only once and canbe accessed globally using the singleton
 * pattern
 */
class Gates final : public internal::Singleton<const Gates> {
  friend class internal::Singleton<const Gates>;

 public:
  // define various quantum gates as member variables
  cmat Id2{cmat::Identity(2, 2)};
  // hadamard gate
  cmat H{cmat::Zero(2, 2)};
  // pauli sigma X gate
  cmat X{cmat::Zero(2, 2)};
  // pauli sigma Y gate
  cmat Y{cmat::Zero(2, 2)};
  // pauli sigma Z gate
  cmat Z{cmat::Zero(2, 2)};
  /// S gate
  cmat S{cmat::Zero(2, 2)};
  // T gate
  cmat T{cmat::Zero(2, 2)};

  // two qubit gates

  // controlled-not controll target gate
  cmat CNOT{cmat::Identity(4, 4)};
  // controlled phase gate
  cmat CZ{cmat::Identity(4, 4)};
  // controlled not gate control gate
  cmat CNOTba{cmat::Zero(4, 4)};
  // swap gate
  cmat SWAP{cmat::Identity(4, 4)};

  // toffoli gate
  cmat TOF{cmat::Identity(8, 8)};
  // fredkin gate
  cmat FRED{cmat::Identity(8, 8)};

 private:
  /**
   * @brief Construct for the gates class
   *
   * the constructor is private and can only be accessed the Singleton pattern. it initializes
   * the quantum gates with their corresponding matrices
   */
  Gates() {
    // initialize quantum gates with their corresponding matrices
    H << 1 / std::sqrt(2.), 1 / std::sqrt(2.), 1 / std::sqrt(2.), -1 / std::sqrt(2.);
    X << 0, 1, 1, 0;
    Z << 1, 0, 0, -1;
    Y << 0, -1_i, 1_i, 0;
    S << 1, 0, 0, 1_i;
    T << 1, 0, 0, std::exp(1_i * clara::pi / 4.0);
    CNOT.block(2, 2, 2, 2) = X;
    CNOTba(0, 0) = 1;
    CNOTba(1, 3) = 1;
    CNOTba(2, 2) = 1;
    CNOTba(3, 1) = 1;
    CZ(3, 3) = -1;

    SWAP.block(1, 1, 2, 2) = X;
    TOF.block(6, 6, 2, 2) = X;
    FRED.block(4, 4, 4, 4) = SWAP;
  }
  // DECLARE constructor private to prefent direct instantion
  ~Gates() = default;

 public:
  /**
   * @brief qubit rotation about a 3-dimensional real unit vector (Rn) gate
   *
   * @param theta the angle of rtation in radians
   * @param n A 3-dimensional vector representing the unit vector about which to rotate
   * @return the rotation gate as 2x2 complex matrix
   *
   * NOTE: the rn Gate rotates the qubit about the specified 3-dimensional real unit vector
   * by the angle theta. the function returns the rotation matrix
   */
  cmat Rn(double theta, const std::vector<double>& n) const {
    if (n.size() != 3) {
      throw exception::CustomException("clara::Gates::Rn()", "n is not a 3-dimensional vector!");
    }
    cmat result(2, 2);
    result =
        std::cos(theta / 2) * Id2 - 1_i * std::sin(theta / 2) * (n[0] * X + n[1] * Y + n[2] * Z);
    return result;
  }

  /**
   * @brief the 2x2 rotation matrix around the X-axis in quantum circuit
   * @param theta the angle of rotation in radians
   * @return the 2x2 rotation amtrix corresponding to the RX gate
   */
  cmat RX(double theta) const {
    // call the Rn function to generate the rotation matrix for the RX gate
    // the axis of rotation is {1, 0, 0} which corresponds to the X-axis
    return Rn(theta, {1, 0, 0});
  }

  /**
   * @brief calculates the 2x2 rotation matrix around the Y axis in quamtum circuit
   * @param theta the angle of rotation in radians
   * @return the 2x2 rotation matrix corresponding to the RY gate
   */
  cmat RY(double theta) const {
    // call the Rn function to generate the rotation matrix for the RY gate
    // the axis of rotation is {0, 1, 0}, which corresponds to the Y-axis
    return Rn(theta, {0, 1, 0});
  }

  /**
   * @brief calculate the 2x2 rotation amtrix around the Z-axis in a quantum circuit
   * @param theta the angle of rotation in radians
   * @return the 2x2 rotation matrix corresponding to the RZ gate
   */
  cmat RZ(double theta) const {
    // call the Rn function to generate the rotation matrix for the RZ gate
    // the axis of rotation is {0, 0, 1} which corresponds to the Z-axis
    return Rn(theta, {0, 0, 1});
  }

  /**
   * @brief generalized Z (Zd) gate for qudits
   *
   * @param D the dimension of the qudit (default 2)
   * @return the Zd gate as DxD complex matrix
   *
   * NOTE: the Zd gate is a generalization of the pauli-Z gate (Z) for qudit (quantum system with
   * dimension D) its defined as \f$ Z = \sum_{j=0}^{D-1} \exp(2\pi \mathrm{i} j/D) |j\rangle\langle
   * j| \f$. the function returns the Zd gate matrix for the specified dimension D
   */
  cmat Zd(idx D = 2) const {
    if (D == 0) {
      throw exception::DimsInvalid("clara::Gates::Zd()");
    }
    cmat result = cmat::Zero(D, D);
    for (idx i = 0; i < D; ++i) {
      result(i, i) = std::pow(omega(D), static_cast<double>(i));
    }
    return result;
  }

  /**
   * @brief generate the SWAP gate for a given dimension
   * the SWAP gate exchange the position of two qubits
   * @param D dimension of the qubit
   * @return the SWAP gate matrix of size DxDxDxD
   *
   * @throw DimsInvalid exception if D is invalid
   */
  cmat SWAPd(idx D = 2) const {
    if (D == 0)
      throw exception::DimsInvalid("clara::Gates::SWAPd()");
    cmat result = cmat::Zero(D * D, D * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif
    for (idx j = 0; j < D; ++j)
      for (idx i = 0; i < D; ++i)
        result(D * i + j, i + D * j) = 1;
    return result;
  }

  /**
   * @brief Fourier transform (Fd) gate for qudit
   *
   * @param D the dimension of the qudit (default: 2)
   * @return the Fd gate as a DxD complex matrix
   *
   * NOTE: the Fd gate is the fourier gate for qudits (quantum system with dimension D).
   * it is defined as \f$ F = \sum_{j,k=0}^{D-1} \exp(2\pi \mathrm{i} jk/D) |j\rangle\langle k| \f$.
   * the function returns the Fd gate matrix for the specified dimension D
   */
  cmat Fd(idx D = 2) const {
    if (D == 0)
      throw exception::DimsInvalid("clara::Gates::Fd()");
    cmat result(D, D);
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // DEBUG
    for (idx j = 0; j < D; ++j) {
      for (idx i = 0; i < D; ++i) {
        result(i, j) = 1 / std::sqrt(D) * std::pow(omega(D), static_cast<double>(i * j));
      }
    }
    return result;
  }

  /**
   * @brief genrate MODMUL gate for a given value a, N, and n
   * the MODMUL gate implements modular multiplication
   * @param a the integer value a
   * @param N the modulis value N
   * @param n the number of qubits (size of the gate matrix)
   * @return the modmul gate matrix of size 2^n x 2^n
   *
   * @throws OutOfRange exception if a, N or n is out of valid range
   */
  [[clara::parallel]] cmat MODMUL(idx a, idx N, idx n) const {
    assert(gcd(a, N) == 1);
    if (N < 3 || a >= N) {
      throw exception::OutOfRange("clara::Gates::MODMUL()");
    }
    if (n < static_cast<idx>(std::ceil(std::log2(N)))) {
      throw exception::OutOfRange("clara::Gtes::MODMUL()");
    }

    // calculate the dimension of the gate matrix
    idx D = static_cast<idx>(std::llround(std::pow(2, n)));
    cmat result = cmat::Zero(D, D);

#ifdef CLARA_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // WITH_OPENMP_
    // poplulate the MODMUL gate matrix using a loop
    for (idx j = 0; j < N; ++j) {
      for (idx i = 0; i < N; ++i) {
        if (static_cast<idx>(modmul(static_cast<bigint>(j), static_cast<bigint>(a),
                                    static_cast<bigint>(N))) == i) {
          result(i, j) = 1;
        }
      }
    }

#ifdef CLARA_OPENMP_
#pragma omp parallel for
#endif  // WITH_OPENMP_
    // set diagonal elements of the gate matrix for remaining indices
    for (idx i = N; i < D; ++i) {
      result(i, i) = 1;
    }
    return result;
  }

  /**
   * @brief generalized X (Xd) gate for qudit
   *
   * @param D the dimension of the qudit (default 2)
   * @return the xd Xd gate as DxD complex matrix
   *
   * NOTE: the Xd gate is a generalization of the Pauli gate (x) for qudit (quantum system with
   * dimension D) it is defined as \f$ X = \sum_{j=0}^{D-1} |j\oplus 1\rangle\langle j| \f$, where
   * \f$ \oplus \f$ denotes addition modulo D. the Xd gate acts as the raising operatos f$
   * X|j\rangle = |j\oplus 1\rangle \f$ the function return the Xd gate for the specified dimension
   * D
   */
  cmat Xd(idx D = 2) const {
    if (D == 0) {
      throw exception::DimsInvalid("clara::Gates::Xd()");
    }

    if (D == 2) {
      return X;
    }
    return Fd(D).inverse() * Zd(D) * Fd(D);
  }

  /**
   * @brief identity gate (id) for qudits
   *
   * @tparam Derived the type of the of the return matrix (default: complex matrix)
   * @param D the dimension of the qudit (default: 2)
   * @return the Id gate as DxD matrix of the specified matrix for qudits (qunatum system with
   * dimension D)
   *
   * NOTE: the template parameter 'Derived' can be explicitly specified to change the
   * return type for the default complex matrix
   */
  template <typename Derived = cmat>
  Derived Id(idx D = 2) const {
    if (D == 0) {
      throw exception::DimsInvalid("clara::Gates::Id()");
    }
    return Derived::Identity(D, D);
  }

  /**
   * @brief embed a small unitary matrix into a larger hilbert space a quantum gate
   *
   * this function applies quantum gate `A` specific subsystem (`target`) within
   * composite system describe by `dims`, resulting matrix acts as identity on all
   * other subsystem
   *
   * @tparam Derived type derived from Eigen::MatrixBase
   * @param A gate matrix to embed, must be square and match the dimension of the target
   * @param target indices of subsystem on which the gate acts
   * @param dims vector of dimension for each subsystem in the full system
   * @return embedded gate matrix of size D x D where D = product(dims)
   */
  template <typename Derived>
  dyn_mat<typename Derived::Scalar> GATE(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& target,
                                         const std::vector<idx>& dims) const {
    // get derived type (for access row, cols)
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    if (!internal::check_nonzero_size(rA)) {
      throw exception::ZeroSize("clara::Gates::GATE()", "A");
    }

    if (!internal::check_square_mat(rA)) {
      throw exception::MatrixNotSquare("clara::Gates::GATE()", "A");
    }

    if (target.empty()) {
      throw exception::ZeroSize("clara::Gates::GATE()", "target");
    }

    if (!internal::check_dims(dims)) {
      throw exception::DimsInvalid("clara::Gates::GATE()", "dims");
    }

    if (!internal::check_subsys_match_dims(target, dims)) {
      throw exception::SubsysMismatchDims("clara::Gates::GATE()", "dims/target");
    }

    // compute total dimension of the gate's subspace (product of target subsystem)
    using Index = typename dyn_mat<typename Derived::Scalar>::Index;

    idx DA = 1;
    for (idx elem : target) {
      DA *= dims[elem];
    }

    // make sure input matrix matching the expected gate dimension
    if (rA.rows() != static_cast<Index>(DA)) {
      throw exception::MatrixMismatchSubsys("clara::Gates::GATE()", "A/dims/target");
    }

    // internal variable using for indexing
    idx Cdims[internal::maxn];     // full system dimension
    idx midx_row[internal::maxn];  // multi-index for row
    idx midx_col[internal::maxn];  // multi-index for column

    idx CdimsA[internal::maxn];     // dimension of target subsystem
    idx midxA_row[internal::maxn];  // row indices for gate subspace
    idx midxA_col[internal::maxn];  // column indices for gate subspace

    idx CdimsA_bar[internal::maxn];   // dimension of non-target subsystem
    idx Csubsys_bar[internal::maxn];  // indices of non-target subsystem
    idx midx_bar[internal::maxn];     // multi-index for non-target subspace

    idx n = dims.size();                   // total number subsystem
    idx n_gate = target.size();            // number of subsystem acted upon
    idx n_subsys_bar = n - target.size();  // number of untouched subsystem

    // get indices subsystem not acted upon (complement or target)
    std::vector<idx> subsys_bar = complement(target, n);

    // total dimension of entire system
    idx D = prod(dims);
    idx Dsubsys_bar = 1;

    for (idx elem : subsys_bar) {
      Dsubsys_bar *= dims[elem];
    }

    // copy untargeted indices and dimension for easier access
    std::copy(subsys_bar.begin(), subsys_bar.end(), std::begin(Csubsys_bar));
    for (idx k = 0; k < n_subsys_bar; ++k) {
      CdimsA_bar[k] = dims[subsys_bar[k]];
      midx_bar[k] = 0;  // initialize to zero
    }

    // initialize dimension and indices for target subsystem
    for (idx k = 0; k < n_subsys_bar; ++k) {
      CdimsA_bar[k] = dims[subsys_bar[k]];
      midx_bar[k] = 0;
    }

    for (idx k = 0; k < n_gate; ++k) {
      midxA_row[k] = midxA_col[k] = 0;
      CdimsA[k] = dims[target[k]];
    }

    // start with an identity matrix of full system size
    dyn_mat<typename Derived::Scalar> result = dyn_mat<typename Derived::Scalar>::Identity(D, D);

    // loop over all configuration of the untouched subsystem
    for (idx i = 0; i < Dsubsys_bar; ++i) {
      internal::n2multiidx(i, n_subsys_bar, CdimsA_bar, midx_bar);
      // loop over rows of the gate matrix
      for (idx a = 0; a < DA; ++a) {
        internal::n2multiidx(a, n_gate, CdimsA, midxA_row);

        // set row indices in the full system
        for (idx k = 0; k < n_gate; ++k) {
          midx_row[target[k]] = midxA_row[k];
        }

        // set bar indices (untouched subsystem) in both row and cols
        for (idx k = 0; k < n_subsys_bar; ++k) {
          midx_row[Csubsys_bar[k]] = midx_col[Csubsys_bar[k]] = midx_bar[k];
        }

        // loop over all columns of the gate matrix
        for (idx b = 0; b < DA; ++b) {
          internal::n2multiidx(b, n_gate, CdimsA, midxA_col);

          // set column indices in the full system
          for (idx k = 0; k < n_gate; ++k) {
            midx_col[target[k]] = midxA_col[k];
          }

          // compute linear indices and assign gate value
          result(internal::multiidx2n(midx_row, n, Cdims),
                 internal::multiidx2n(midx_col, n, Cdims)) = rA(a, b);
        }
      }
    }

    return result;
  }

  /**
   * @brief overloaded version GATE for system with uniform subystem dimension
   *
   * this function provide a simplified interface to the full `GATE(...)` function
   * by assuming that all subystem in the system have the same dimension `d`
   *
   * @tparam Derived type from Eigen::MatrixBase
   * @param A gate matrix to embed, must be square and match dimension of the target
   * @param target indices of subystem on which the gate acts
   * @param n total_number total number of subsystem in composite system
   * @param d dimension of each subsystem (default: 2 for qubits)
   * @return embedded gate matrix of size d^n x d^n
   */
  template <typename Derived>
  dyn_mat<typename Derived::Scalar> GATE(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& target, idx n, idx d = 2) const {
    if (d == 0) {
      throw exception::DimsInvalid("clara::Gates::GATE()", "d");
    }

    // forward to the general version GATE that acc vector of dimension
    return GATE(A, target, std::vector<idx>(n, d));
  }

  /**
   * @brief generate the multi-paritite multi-controlled in matrix form
   *
   * @tparam Derived the type of the input matrix A and the return matrix
   * @param A the gate matrix to be controlled
   * @param ctrl list of control qubit indices (0-bates) for the gate
   * @param subsy list of subsystem qubit indices for the gate
   * @param B the total number qubits
   * @param d the local dimension of each qubit
   * @return the controlled gate matrix as complex matrix of the specified mtype
   *
   * NOTE: `CTRL` function generate the multi-paritite multi-controlled gate in matrix form
   * the function takes an input matrix `A` (the gate to be controlled) and lists of control and
   * subsystem qubit indices. it returns the controlled gate matrix as a complex matrix of the
   * specfied type
   */
  template <typename Derived>
  dyn_mat<typename Derived::Scalar> CTRL(
      const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& ctrl,
      const std::vector<idx>& target, idx n, idx d = 2,
      std::optional<std::vector<idx>> shift = std::nullopt) const {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // check matrix zero size
    if (!internal::check_nonzero_size(rA)) {
      throw exception::ZeroSize("clara::Gates::CTRL()");
    }

    // check square matrix
    if (!internal::check_square_mat(rA)) {
      throw exception::MatrixNotSquare("clara::Gates::CTRL()");
    }

    // check list zero size
    if (ctrl.empty()) {
      throw exception::ZeroSize("clara::Gates::CTRL()");
    }

    // check out of range
    if (n == 0) {
      throw exception::OutOfRange("clara::Gates::CTRL()");
    }
    // check valid local dimension
    if (d == 0) {
      throw exception::OutOfRange("clara::Gates::CTRL()");
    }

    std::vector<idx> ctrlgate = ctrl;
    std::sort(ctrlgate.begin(), ctrlgate.end());
    ctrlgate.insert(ctrlgate.end(), target.begin(), target.end());

    for (idx elem_ctrl : ctrl) {
      for (idx elem_target : target) {
        if (elem_ctrl == elem_target) {
          throw exception::OutOfRange("clara::Gates::CTRL()", "ctrl/target");
        }
      }
    }

    std::vector<idx> dims(n, d);
    if (!internal::check_subsys_match_dims(ctrlgate, dims)) {
      throw exception::SubsysMismatchDims("clara::Gates::CTRL()", "ctrl/dims");
    }

    idx DA = rA.rows();
    if (DA != internal::safe_pow(d, target.size())) {
      throw exception::MatrixMismatchSubsys("clara::Gates::CTRL()", "A/d/target");
    }

    if (shift.has_value()) {
      for (idx& elem : shift.value()) {
        if (elem >= d) {
          throw exception::OutOfRange("clara::Gates::CTRL()", "shift");
        }
        elem = d - elem;
      }
    }

    if (!shift.has_value()) {
      shift = std::vector<idx>(ctrl.size(), 0);
    }

    idx D = prod(dims);
    idx Dctrl = internal::safe_pow(d, ctrl.size());
    idx ctrl_size = ctrl.size();

    dyn_mat<typename Derived::Scalar> result = dyn_mat<typename Derived::Scalar>::Zero(D, D);

    std::vector<idx> ctrl_bar = complement(ctrlgate, n);
    std::vector<idx> ctrlgate_bar = complement(ctrlgate, n);
    idx Dctrlagate_bar = 1;
    for (idx elem : ctrlgate_bar) {
      Dctrlagate_bar *= dims[elem];
    }

    dyn_mat<typename Derived::Scalar> Id_ctrlgate_bar =
        dyn_mat<typename Derived::Scalar>::Identity(Dctrlagate_bar, Dctrlagate_bar);

    dyn_mat<typename Derived::Scalar> prj_bar =
        dyn_mat<typename Derived::Scalar>::Identity(Dctrl, Dctrl);

    for (idx k = 0; k < d; ++k) {
      std::vector<idx> ctrl_shift = shift.value();

      std::transform(ctrl_shift.begin(), ctrl_shift.end(), ctrl_shift.begin(),
                     [k, d](idx elem) { return (elem + k) % d; });

      idx pos = multiidx2n(ctrl_shift, std::vector<idx>(ctrl_size, d));

      dyn_mat<typename Derived::Scalar> prj_mat =
          dyn_mat<typename Derived::Scalar>::Zero(Dctrl, Dctrl);

      prj_bar -= prj_mat;

      dyn_mat<typename Derived::Scalar> Ak = powm(rA, k);
      dyn_mat<typename Derived::Scalar> gate = kron(prj_mat, Ak);

      result += GATE(gate, ctrlgate, n, d);
    }

    result +=
        GATE(kron(prj_bar, dyn_mat<typename Derived::Scalar>::Identity(DA, DA)), ctrlgate, n, d);

    return result;
  }

  /**
   * @brief expand the matrix A to a larger tensor product space
   *
   * @tparam Derived the type of the input matrix A and the return matrix
   * @param A the input matrix to be expanded
   * @param pos the position of the tensor product space to which a should expanded
   * @param dims the dimension of the tebsor product space
   * @return the expanded matrix as a complex matrix of the specified type
   *
   * NOTE: `expandout` function expand the matrix A to a larger tensor product space specified by
   * the dimension in the `pos` parameter ditermines the position at which the expansion should be
   * done the function returns the expanded matrix as a complex matrix of the specfied type. the
   * input matrix A should be square, and the dimension at position in `dims` should be valid fro a
   * tensor product space
   *
   * @throws clara::exception::ZeroSize if the input matrix A has zero size
   * @throws clara::exception::DimsInvalid the dimension in `dims` are invalid
   * @throws clara::exception::MatrixNotSquare if the input matrix A is not square
   * @throws clara::exception::OutOfRange if the `pos` parameter is out of range for the dimension
   * in `dims`
   * @throws clara::exception::DimsMismatchMatrix if the dimension of A does not match dimension at
   * position `pos in `dims`
   *
   * `
   */
  template <typename Derived>
  dyn_mat<typename Derived::Scalar> expandout(const Eigen::MatrixBase<Derived>& A, idx pos,
                                              const std::vector<idx>& dims) const {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();
    if (!internal::check_nonzero_size(rA)) {
      throw exception::ZeroSize("clara::Gates::expandout()");
    }
    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims)) {
      throw exception::DimsInvalid("clara::Gates::expandout()");
    }
    // check square matrix
    if (!internal::check_square_mat(rA)) {
      throw exception::MatrixNotSquare("clara::Gates::expandout()");
    }

    // check that position is valid
    if (pos > dims.size() - 1) {
      throw exception::OutOfRange("clara::Gates::expandout()");
    }
    // check that dims[pos] match dimension of A
    if (static_cast<idx>(rA.rows()) != dims[pos]) {
      throw exception::DimsMismatchMatrix("clara::Gates::expandout()");
    }

    return GATE(A, {pos}, dims);
  }

  /**
   * @brief expand the matri A to larger tensor product space using a list of dimensional
   *
   * @tparam Derived the type of the input matrix and the return matrix
   * @param A the input matri to expanded
   * @param pos the position of the tensor product space to which A should be expanded
   * @param dims the dimensions of the tensor product space as an initializer list
   * @return the expanded matrix a complex matrix of the specfied type
   *
   * NOTE: the `expandout` function expands the matrix A to larger tensor product space specfied by
   * the number of subsystem `N`. the `pos` parameter determines the position at which the expansion
   * should be done. the `d` parameter specifies the local dimension of each subsystem (default is 2
   * for qubits) the function returns the expanded matrix as a complex matrix of the specfied
   */
  template <typename Derived>
  dyn_mat<typename Derived::Scalar> expandout(const Eigen::MatrixBase<Derived>& A, idx pos,
                                              const std::initializer_list<idx>& dims) const {
    return expandout(A, pos, std::vector<idx>(dims));
  }

  /**
   * @brief expand a given matrix into a larger multi-parite system
   * @param A the input matrix to be expanded
   * @param pos the position (index) where the matrix A will be inserted
   * @param N the total number of parties
   * @param d the dimension of each subsystem (default is 2)
   * @return expanded matrix in the multi-paritite system
   *
   * @throws ZeroSize exception if the input matrix matrix has zero size
   * @throws DimsInvalid exception if d is invalid
   */
  template <typename Derived>
  dyn_mat<typename Derived::Scalar> expandout(const Eigen::MatrixBase<Derived>& A, idx pos, idx n,
                                              idx d = 2) const {
    if (!internal::check_nonzero_size(A)) {
      throw exception::ZeroSize("clara::Gates::expandout()");
    }
    // check valid dims
    if (d == 0) {
      throw exception::DimsInvalid("clara::Gates::expandout()");
    }

    return expandout(A, pos, std::vector<idx>(n, d));
  }

  /**
   * @brief return the name of a quantum gate based on its matrix representation
   *
   * utility compare the given unitary matrix `U` againts known standard quantum gates
   * if a match is found, the corresponding gate name is return as an `std::optional<std::string>`
   * if not match is found, an empty optional is returned
   *
   * @param U complex matrix representing a quantum gate
   * @return std::optional<std::string> name of the matchied gate
   */
  std::optional<std::string> get_name(const cmat& U) const {
    // validate input matrix non-empty
    if (!internal::check_nonzero_size(U)) {
      throw exception::ZeroSize("clara::Gates::gate_name()", "U");
    }

    if (!internal::check_square_mat(U)) {
      return {};  // return empty if matrix not square
    }

    // get dimension of the matrix (number of row / columns)
    const idx D = static_cast<idx>(U.rows());

    // handle switch-case for handle different gate dimensions
    switch (D) {
      case 2:
        // single-qubit gate (dimension 2x2)
        if (U == Id2) {
          return "Id2";  // identity gate
        } else if (U == H) {
          return "H";  // hadamard gate
        } else if (U == X) {
          return "X";  // pauli-x gate
        } else if (U == Y) {
          return "Y";  // pauli-y gate
        } else if (U == Z) {
          return "Z";  // pauli-z gate
        } else if (U == S) {
          return "S";  // phase gate
        } else if (U == adjoint(S)) {
          return "S+";  // adjoint of phase gate
        } else if (U == T) {
          return "T";  // pi / 8 gate
        } else if (U == adjoint(T)) {
          return "T+";  // adjoint of T gate
        } else {
          return {};
        }
      case 4:
        // two-qubit gate (dimension 4x4)
        if (U == CNOT) {
          return "CNOT";  // controlled-not gate
        } else if (U == CZ) {
          return "CZ";  // controlled-z gate
        } else if (U == CNOTba) {
          return "CNOTba";  // iverse cnot (target -> control)
        } else if (U == SWAP) {
          return "SWAP";  // swap gate
        } else {
          return {};
        }
      case 8:
        // three-qubit gate (dimension 8x8)
        if (U == TOF) {
          return "TOF";  // toffoli gate
        } else if (U == FRED) {
          return "FRED";  // fredkin gate
        } else {
          return {};
        }

      default:
        return {};
    }
  }
};

}  // namespace clara

#endif  // !CLASSFUNCTION_GATES_H_
