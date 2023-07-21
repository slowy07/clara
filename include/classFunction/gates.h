#ifndef CLASSFUNCTION_GATES_H_
#define CLASSFUNCTION_GATES_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>

#include "../functions.h"
#include "../internal/classFunction/singleton.h"
#include "../internal/util.h"
#include "exception.h"

namespace clara {

/**
 * @class clara::Gates
 * @brief const singleton class that implements most commonly used gates
 */
class Gates final : public internal::Singleton<const Gates> {
  friend class internal::Singleton<const Gates>;

 public:
  cmat Id2{cmat::Identity(2, 2)};
  cmat H{cmat::Zero(2, 2)};
  cmat X{cmat::Zero(2, 2)};
  cmat Y{cmat::Zero(2, 2)};
  cmat Z{cmat::Zero(2, 2)};
  cmat S{cmat::Zero(2, 2)};
  cmat T{cmat::Zero(2, 2)};

  cmat CNOT{cmat::Identity(4, 4)};
  cmat CZ{cmat::Identity(4, 4)};
  cmat CNOTba{cmat::Zero(4, 4)};
  cmat SWAP{cmat::Identity(4, 4)};

  cmat TOF{cmat::Identity(8, 8)};
  cmat FRED{cmat::Identity(8, 8)};

 private:
  Gates() {
    H << 1 / std::sqrt(2.), 1 / std::sqrt(2.), 1 / std::sqrt(2.), -1 / std::sqrt(2.);
    X << 0, 1, 1, 0;
    Z << 1, 0, 0, -1;
    Y << 0, -1_i, 1_i, 0;
    S << 1, 0, 0, 1_i;
    T << 1, 0, 0, std::exp(1_i * pi / 4.0);
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
  /*
   * @brief default constructor
   */
  ~Gates() = default;

 public:
  /**
   * @brief Qubit roration of theta about the 3 dimensional real unit vector n
   * @return rotation gate
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
   * @brief generalized Z gates for qudits
   * @note define as \f$ Z = \sum_{j=0}^{D-1} \exp(2\pi \mathrm{i} j/D) |j\rangle\langle j| \f$
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
   * @brief fourier transform gate for qudits
   * @note define as \f$ F = \sum_{j,k=0}^{D-1} \exp(2\pi \mathrm{i} jk/D) |j\rangle\langle k| \f$
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
   * @brief generalized X gate for qudits
   * @note define as f$ X = \sum_{j=0}^{D-1} |j\oplus 1\rangle\langle j| \f$, ie raising operator
   * \f$ X|j\rangle = |j\oplus 1\rangle\f$
   */
  cmat Xd(idx D = 2) const {
    if (D == 0)
      throw exception::DimsInvalid("clara::Gates::Xd()");
    return Fd(D).inverse() * Zd(D) * Fd(D);
  }

  /**
   * @brief inditity gate
   * @note can change the return type from complex matrix (default) by explicitly
   * specifying the template paramter
   */
  template <typename Derived = Eigen::MatrixXcd>
  Derived Id(idx D = 2) const {
    if (D == 0)
      throw exception::DimsInvalid("clara::Gates::Id()");
    return Derived::Identity(D, D);
  }

  /**
   * @brief generate the multi-paritite multiple contrilled A gate in matrix form
   * @note the dimension of the gate A must match the dimension of subsys
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
   * @brief expands out clara::kron()
   * expand out A as a matrix in a multi-paritite system. faster than using clara::kron()
   * @return tensor product
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
   * @brief expand out
   * @note the std::initializer_list overload exists because otherwise, in the
   * degenerate case when dims only one element, the one element list is implicitly
   * converted to the element's underlying type
   */
  template <typename Derived>
  dyn_mat<typename Derived::Scalar> expandout(const Eigen::MatrixBase<Derived>& A, idx pos,
                                              const std::initializer_list<idx>& dims) const {
    return this->expandout(A, pos, std::vector<idx>(dims));
  }

  /**
   * @brief expands out
   * @note expands out A as a matrix in a multi-paritite system
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
