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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>

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
    if (n.size() != 3)
      throw exception::CustomException("clara::Gates::Rn()", "n is not a 3-dimensional vector!");
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
    if (D == 0)
      throw exception::DimsInvalid("clara::Gates::Zd()");
    cmat result = cmat::Zero(D, D);
    for (idx i = 0; i < D; ++i)
      result(i, i) = std::pow(omega(D), static_cast<double>(i));
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
    for (idx j = 0; j < D; ++j)
      for (idx i = 0; i < D; ++i)
        result(i, j) = 1 / std::sqrt(D) * std::pow(omega(D), static_cast<double>(i * j));
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
  cmat MODMUL(idx a, idx N, idx n) const {
#ifndef DEBUG
    assert(gcd(a, N) == 1);
#endif  // !DEBUG
    if (N < 3 || a >= N) {
      throw exception::OutOfRange("clara::Gates::MODMUL()");
    }
    if (n < static_cast<idx>(std::ceil(std::log2(N)))) {
      throw exception::OutOfRange("clara::Gtes::MODMUL()");
    }

    // calculate the dimension of the gate matrix
    idx D = static_cast<idx>(std::llround(std::pow(2, n)));
    cmat result = cmat::Zero(D, D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // WITH_OPENMP_
    // poplulate the MODMUL gate matrix using a loop
    for (idx j = 0; j < N; ++j)
      for (idx i = 0; i < N; ++i)
        if (static_cast<idx>(modmul(j, a, N)) == i)
          result(i, j) = 1;

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // WITH_OPENMP_
    // set diagonal elements of the gate matrix for remaining indices
    for (idx i = N; i < D; ++i)
      result(i, i) = 1;
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
    if (D == 0)
      throw exception::DimsInvalid("clara::Gates::Xd()");
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
  template <typename Derived = Eigen::MatrixXcd>
  Derived Id(idx D = 2) const {
    if (D == 0)
      throw exception::DimsInvalid("clara::Gates::Id()");
    return Derived::Identity(D, D);
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
  dyn_mat<typename Derived::Scalar> CTRL(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& ctrl,
                                         const std::vector<idx>& subsys, idx N, idx d = 2) const {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // check matrix zero size
    if (!internal::check_nonzero_size(rA))
      throw exception::ZeroSize("clara::Gates::CTRL()");

    // check square matrix
    if (!internal::check_square_mat(rA))
      throw exception::MatrixNotSquare("clara::Gates::CTRL()");

    // check list zero size
    if (ctrl.size() == 0)
      throw exception::ZeroSize("clara::Gates::CTRL()");
    if (subsys.size() == 0)
      throw exception::ZeroSize("clara::Gates::CTRL()");

    // check out of range
    if (N == 0)
      throw exception::OutOfRange("clara::Gates::CTRL()");
    // check valid local dimension
    if (d == 0)
      throw exception::OutOfRange("clara::Gates::CTRL()");

    // control gate subsystem
    std::vector<idx> ctrlgate = ctrl;
    ctrlgate.insert(std::end(ctrlgate), std::begin(subsys), std::end(subsys));
    std::sort(std::begin(ctrlgate), std::end(ctrlgate));

    std::vector<idx> dims(N, d);

    /**
     * check that control + gats subystem is valid
     * with respect to local dimensions
     */
    if (!internal::check_subsys_match_dims(ctrlgate, dims))
      throw exception::SubsysMismatchdims("clara::Gates::CTRL()");

    // check that subsys list match the dimension of the matrix
    using Index = typename dyn_mat<typename Derived::Scalar>::Index;
    if (A.rows() != static_cast<Index>(std::llround(std::pow(d, subsys.size()))))
      throw exception::DimsMismatchMatrix("clara::Gates::CTRL()");

    idx Cdims[maxn];
    idx midx_row[maxn];
    idx midx_col[maxn];

    idx CdimsA[maxn];
    idx midxA_row[maxn];
    idx midxA_col[maxn];

    idx Cdims_bar[maxn];
    idx Csubsys_bar[maxn];
    idx midx_bar[maxn];

    idx Ngate = subsys.size();
    idx Nctrl = ctrl.size();
    idx Nsubsys_bar = N - ctrlgate.size();
    idx D = static_cast<idx>(std::llround(std::pow(d, N)));
    idx DA = static_cast<idx>(rA.rows());
    idx Dsubsys_bar = static_cast<idx>(std::llround(std::pow(d, Nsubsys_bar)));

    // compute the complementary subystem of control gate w.r.t
    std::vector<idx> subsys_bar = complement(ctrlgate, N);
    std::copy(std::begin(subsys_bar), std::end(subsys_bar), std::begin(Csubsys_bar));

    for (idx k = 0; k < N; ++k) {
      midx_row[k] = midx_col[k] = 0;
      Cdims[k] = d;
    }

    for (idx k = 0; k < Nsubsys_bar; ++k) {
      Cdims_bar[k] = d;
      midx_bar[k] = 0;
    }

    for (idx k = 0; k < Ngate; ++k) {
      midxA_row[k] = midxA_col[k] = 0;
      CdimsA[k] = d;
    }

    dyn_mat<typename Derived::Scalar> result = dyn_mat<typename Derived::Scalar>::Identity(D, D);
    dyn_mat<typename Derived::Scalar> Ak;

    // run over complement indexes
    for (idx i = 0; i < Dsubsys_bar; ++i) {
      internal::n2multiidx(i, Nsubsys_bar, Cdims_bar, midx_bar);
      for (idx k = 0; k < d; ++k) {
        Ak = powm(rA, k);
        // run ober the subsys row multi-index
        for (idx a = 0; a < DA; ++a) {
          // get the subsys row multi-index
          internal::n2multiidx(a, Ngate, CdimsA, midxA_row);
          /**
           * construct the result row multi-index
           * first the ctrl part (equal for both row and colum)
           */
          for (idx c = 0; c < Nctrl; ++c)
            midx_row[ctrl[c]] = midx_col[ctrl[c]] = k;
          // then the complement part
          for (idx c = 0; c < Nsubsys_bar; ++c)
            midx_row[Csubsys_bar[c]] = midx_col[Csubsys_bar[c]] = midx_bar[c];
          // then the subsys part
          for (idx c = 0; c < Ngate; ++c)
            midx_row[subsys[c]] = midxA_row[c];
          // run over the subsys column multi-index
          for (idx b = 0; b < DA; ++b) {
            internal::n2multiidx(b, Ngate, CdimsA, midxA_col);
            // construct the result column multi-index
            for (idx c = 0; c < Ngate; ++c)
              midx_col[subsys[c]] = midxA_col[c];
            // write the value
            result(internal::multiidx2n(midx_row, N, Cdims),
                   internal::multiidx2n(midx_col, N, Cdims)) = Ak(a, b);
          }
        }
      }
    }
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
    if (!internal::check_nonzero_size(rA))
      throw exception::ZeroSize("clara::Gates::expandout()");
    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
      throw exception::DimsInvalid("clara::Gates::expandout()");
    // check square matrix
    if (!internal::check_square_mat(rA))
      throw exception::MatrixNotSquare("clara::Gates::expandout()");

    // check that position is valid
    if (pos > dims.size() - 1)
      throw exception::OutOfRange("clara::Gates::expandout()");
    // check that dims[pos] match dimension of A
    if (static_cast<idx>(rA.rows()) != dims[pos])
      throw exception::DimsMismatchMatrix("clara::Gates::expandout()");

    idx D = std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                            std::multiplies<idx>());
    dyn_mat<typename Derived::Scalar> result = dyn_mat<typename Derived::Scalar>::Identity(D, D);
    idx Cdims[maxn];
    idx midx_row[maxn];
    idx midx_col[maxn];

    for (idx k = 0; k < dims.size(); k++) {
      midx_row[k] = midx_col[k] = 0;
      Cdims[k] = dims[k];
    }

    // run over the main diagonal multi-index
    for (idx i = 0; i < D; ++i) {
      internal::n2multiidx(i, dims.size(), Cdims, midx_row);
      // get column multi_index
      internal::n2multiidx(i, dims.size(), Cdims, midx_col);
      // run over the gate row multi-index
      for (idx a = 0; a < static_cast<idx>(rA.rows()); ++a) {
        midx_row[pos] = a;
        // run over the gate column multi-index
        for (idx b = 0; b < static_cast<idx>(rA.cols()); ++b) {
          // construct the total column multi-index
          midx_col[pos] = b;
          result(internal::multiidx2n(midx_row, dims.size(), Cdims),
                 internal::multiidx2n(midx_col, dims.size(), Cdims)) = rA(a, b);
        }
      }
    }
    return result;
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
    return this->expandout(A, pos, std::vector<idx>(dims));
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
  dyn_mat<typename Derived::Scalar> expandout(const Eigen::MatrixBase<Derived>& A, idx pos, idx N,
                                              idx d = 2) const {
    if (!internal::check_nonzero_size(A))
      throw exception::ZeroSize("clara::Gates::expandout()");
    // check valid dims
    if (d == 0)
      throw exception::DimsInvalid("clara::Gates::expandout()");

    std::vector<idx> dims(N, d);
    return this->expandout(A, pos, dims);
  }
};

}  // namespace clara

#endif  // !CLASSFUNCTION_GATES_H_
