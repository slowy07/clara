#ifndef OPERATIONS_H_
#define OPERATIONS_H_

#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <ios>
#include <iterator>
#include <tuple>
#include <utility>

#include "classFunction/exception.h"
#include "constants.h"
#include "functions.h"
#include "internal/util.h"
#include "types.h"

namespace clara {

template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> applyCTRL(const Eigen::MatrixBase<Derived1>& state,
                                             const Eigen::MatrixBase<Derived2>& A,
                                             const std::vector<idx>& ctrl,
                                             const std::vector<idx>& subsys,
                                             const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate = state.derived();
  const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

  if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
    throw exception::TypeMismatch("clara::applyCTRL()");

  // check zero size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::applyCTRL()");
  if (!internal::check_nonzero_size(rstate))
    throw exception::ZeroSize("clara::applyCTRL()");

  // check square matrix for the gate
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::applyCTRL()");

  // check that all control subsystem have the same dimension
  idx d = ctrl.size() > 0 ? dims[ctrl[0]] : 1;
  for (idx i = 1; i < ctrl.size(); ++i)
    if (dims[ctrl[i]] != d)
      throw exception::DimsNotEqual("clara::applyCTRL()");
  // check that dimension is valid
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::applyCTRL()");

  // check subsys
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::applyCTRL()");

  // check that gate matches the dimension of the subsys
  std::vector<idx> subsys_dims(subsys.size());
  for (idx i = 0; i < subsys.size(); ++i)
    subsys_dims[i] = dims[subsys[i]];
  if (!internal::check_dims_match_mat(subsys_dims, rA))
    throw exception::MatrixMismatchSubsys("clara::applyCTRL()");

  std::vector<idx> ctrlgate = ctrl;
  ctrlgate.insert(std::end(ctrlgate), std::begin(subsys), std::end(subsys));
  std::sort(std::begin(ctrlgate), std::end(ctrlgate));

  if (!internal::check_subsys_match_dims(ctrlgate, dims))
    throw exception::SubsysMismatchdims("clara::applyCTRL()");

  // contruct the table of A^i
  std::vector<dyn_mat<typename Derived1::Scalar>> Ai;
  std::vector<dyn_mat<typename Derived1::Scalar>> Aidagger;
  for (idx i = 0; i < std::max(d, static_cast<idx>(2)); ++i) {
    Ai.push_back(powm(rA, i));
    Aidagger.push_back(powm(adjoint(rA), i));
  }
  
  // total dimension
  idx D = static_cast<idx>(rstate.rows());
  // total number of subsystem
  idx n = dims.size();
  // number of ctrl subsystem
  idx ctrlsize = ctrl.size();
  // number of ctrl + gate subsystem
  idx ctrlgatesize = ctrlgate.size();
  idx subsyssize = subsys.size();
  // dimension of ctrl subsystem
  idx Dctrl = static_cast<idx>(std::llround(std::pow(d, ctrlsize)));
  // dimension of gate subsystem
  idx DA = static_cast<idx>(rA.rows());

  // local dimension
  idx Cdims[maxn];
  idx CdimsA[maxn];
  idx CdimsCTRL[maxn];
  idx CdimsCTRLA_bar[maxn];

  // compute the complementary subsystem
  std::vector<idx> ctrlgate_bar = complement(ctrlgate, n);
  // number of subsystem that are complementary to the ctrlgate
  idx ctrlgate_barsize = ctrlgate_bar.size();
  // dimension of the rest
  idx DCTRLA_bar = 1;
  for (idx i = 0; i < ctrlgate_barsize; ++i)
    DCTRLA_bar *= dims[ctrlgate_bar[i]];

  for (idx k = 0; k < n; ++k)
    Cdims[k] = dims[k];
  for (idx k = 0; k < subsyssize; ++k)
    CdimsA[k] = dims[subsys[k]];
  for (idx k = 0; k < ctrlsize; ++k)
    CdimsCTRL[k] = d;
  for (idx k = 0; k < ctrlgate_barsize; ++k)
    CdimsCTRLA_bar[k] = dims[ctrlgate_bar[k]];

  /**
   * worker computes the coefficient and index for ket case
   * used in #pragma omp parallel fro collapse
   */
  auto coeff_idx_ket = [&](const idx i_, const idx m_,
                           const idx r_) noexcept -> std::pair<typename Derived1::Scalar, idx> {
    idx indx = 0;
    typename Derived1::Scalar coeff = 0;

    idx Cmidx[maxn];
    idx CmidxA[maxn];
    idx CdmixCTRLA_bar[maxn];

    // compute the index
    for (idx k = 0; k < ctrlsize; ++k) {
      Cmidx[ctrl[k]] = i_;
    }

    internal::n2multiidx(r_, n - ctrlgatesize, CdimsCTRLA_bar, CdmixCTRLA_bar);
    for (idx k = 0; k < n - ctrlgatesize; ++k) {
      Cmidx[ctrlgate_bar[k]] = CdmixCTRLA_bar[k];
    }

    // set the A part
    internal::n2multiidx(m_, subsyssize, CdimsA, CmidxA);
    for (idx k = 0; k < subsyssize; ++k) {
      Cmidx[subsys[k]] = CmidxA[k];
    }

    // get the total index
    indx = internal::multiidx2n(Cmidx, n, Cdims);

    // compute the coefficient
    for (idx n_ = 0; n_ < DA; ++n_) {
      internal::n2multiidx(n_, subsyssize, CdimsA, CmidxA);
      for (idx k = 0; k < subsyssize; ++k) {
        Cmidx[subsys[k]] = CmidxA[k];
      }
      coeff += Ai[i_](m_, n_) * rstate(internal::multiidx2n(Cmidx, n, Cdims));
    }
    return std::make_pair(coeff, indx);
  };

  /**
   * worker computes the coefficient and the index
   * for the density matrix case used in #pragma omp parallel for collapse
   */
  auto coeff_idx_rho =
      [&](const idx i1_, const idx m1_, const idx r1_, const idx i2_, const idx m2_,
          const idx r2_) noexcept -> std::tuple<typename Derived1::Scalar, idx, idx> {
    idx idxrow = 0;
    idx idxcol = 0;
    typename Derived1::Scalar coeff = 0, lhs = 1, rhs = 1;

    idx Cmidxrow[maxn];
    idx Cmidxcol[maxn];
    idx CmidxArow[maxn];
    idx CmidxAcol[maxn];
    idx CmidxCTRLrow[maxn];
    idx CmidxCTRLcol[maxn];
    idx CmidxCTRLA_barrow[maxn];
    idx CmidxCTRLA_barcol[maxn];

    // compute the ket/bra indexes

    // set the CTRL part
    internal::n2multiidx(i1_, ctrlsize, CdimsCTRL, CmidxCTRLrow);
    internal::n2multiidx(i2_, ctrlsize, CdimsCTRL, CmidxCTRLcol);

    for (idx k = 0; k < ctrlsize; ++k) {
      Cmidxrow[ctrl[k]] = CmidxCTRLrow[k];
      Cmidxcol[ctrl[k]] = CmidxCTRLcol[k];
    }

    // set the rest
    internal::n2multiidx(r1_, n - ctrlgatesize, CdimsCTRLA_bar, CmidxCTRLA_barrow);
    internal::n2multiidx(r2_, n - ctrlgatesize, CdimsCTRLA_bar, CmidxCTRLA_barcol);
    for (idx k = 0; k < n - ctrlgatesize; ++k) {
      Cmidxrow[ctrlgate_bar[k]] = CmidxCTRLA_barrow[k];
      Cmidxcol[ctrlgate_bar[k]] = CmidxCTRLA_barcol[k];
    }

    // set the A part
    internal::n2multiidx(m1_, subsyssize, CdimsA, CmidxArow);
    internal::n2multiidx(m2_, subsyssize, CdimsA, CmidxAcol);
    for (idx k = 0; k < subsys.size(); ++k) {
      Cmidxrow[subsys[k]] = CmidxArow[k];
      Cmidxcol[subsys[k]] = CmidxAcol[k];
    }
    idxrow = internal::multiidx2n(Cmidxrow, n, Cdims);
    idxcol = internal::multiidx2n(Cmidxcol, n, Cdims);

    bool all_ctrl_rows_equal = true;
    bool all_ctrl_cols_equal = true;

    idx first_ctrl_row, first_ctrl_col;
    if (ctrlsize > 0) {
      first_ctrl_row = CmidxCTRLrow[0];
      first_ctrl_col = CmidxCTRLcol[0];
    } else {
      first_ctrl_row = first_ctrl_col = 1;
    }

    for (idx k = 1; k < ctrlsize; ++k) {
      if (CmidxCTRLrow[k] != first_ctrl_row) {
        all_ctrl_rows_equal = false;
        break;
      }
    }
    for (idx k = 1; k < ctrlsize; ++k) {
      if (CmidxCTRLcol[k] != first_ctrl_col) {
        all_ctrl_cols_equal = false;
        break;
      }
    }

    // at least one control activated, compute the coefficient
    for (idx n1_ = 0; n1_ < DA; ++n1_) {
      internal::n2multiidx(n1_, subsyssize, CdimsA, CmidxArow);
      for (idx k = 0; k < subsyssize; ++k) {
        Cmidxrow[subsys[k]] = CmidxArow[k];
      }
      idx idxrowtmp = internal::multiidx2n(Cmidxrow, n, Cdims);

      if (all_ctrl_rows_equal) {
        lhs = Ai[first_ctrl_row](m1_, n1_);
      } else {
        lhs = (m1_ == n1_) ? 1 : 0;  // identity matrix
      }

      for (idx n2_ = 0; n2_ < DA; ++n2_) {
        internal::n2multiidx(n2_, subsyssize, CdimsA, CmidxAcol);
        for (idx k = 0; k < subsyssize; ++k) {
          Cmidxcol[subsys[k]] = CmidxAcol[k];
        }

        if (all_ctrl_cols_equal) {
          rhs = Aidagger[first_ctrl_col](n2_, m2_);
        } else {
          rhs = (n2_ == m2_) ? 1 : 0;  // identity matrix
        }

        idx idxcoltmp = internal::multiidx2n(Cmidxcol, n, Cdims);

        coeff += lhs * rstate(idxrowtmp, idxcoltmp) * rhs;
      }
    }

    return std::make_tuple(coeff, idxrow, idxcol);
  };

  if (internal::check_cvector(rstate)) {
    if (!internal::check_dims_match_cvect(dims, rstate))
      throw exception::DimsMismatchCvector("clara::applyCTRL()");
    if (D == 1)
      return rstate;
    dyn_mat<typename Derived1::Scalar> result = rstate;

#ifndef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // !WITH_OPENMP_
    for (idx m = 0; m < DA; ++m)
      for (idx r = 0; r < DCTRLA_bar; ++r) {
        if (ctrlsize == 0) {
          result(coeff_idx_ket(1, m, r).second) = coeff_idx_ket(1, m, r).first;
        } else
          for (idx i = 0; i < d; i++) {
            result(coeff_idx_ket(i, m, r).second) = coeff_idx_ket(i, m, r).first;
          }
      }
    return result;
  } else if (internal::check_square_mat(rstate)) {
    // density matrix
    if (!internal::check_dims_match_mat(dims, rstate))
      throw exception::DimsMismatchMatrix("clara::applyCTRL()");

    if (D == 1)
      return rstate;
    dyn_mat<typename Derived1::Scalar> result = rstate;

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(4)
#endif  // DEBUG
    for (idx m1 = 0; m1 < DA; ++m1)
      for (idx r1 = 0; r1 < DCTRLA_bar; ++r1)
        for (idx m2 = 0; m2 < DA; ++m2)
          for (idx r2 = 0; r2 < DCTRLA_bar; ++r2)
            if (ctrlsize == 0) {
              auto coeff_idxes = coeff_idx_rho(1, m1, r1, 1, m2, r2);
              result(std::get<1>(coeff_idxes), std::get<2>(coeff_idxes)) = std::get<0>(coeff_idxes);
            } else {
              for (idx i1 = 0; i1 < Dctrl; ++i1)
                for (idx i2 = 0; i2 < Dctrl; ++i2) {
                  auto coeff_idxes = coeff_idx_rho(i1, m2, r1, i2, m2, r2);
                  result(std::get<1>(coeff_idxes), std::get<2>(coeff_idxes)) =
                      std::get<0>(coeff_idxes);
                }
            }
    return result;
  } else
    throw exception::MatrixNotSquareNotCvector("clara::applyCTRL()");
}

/**
 * @brief applies the controlled-gate A to the part subsys
 * of the multi partite state vector or density matrix state
 * @note the dimension of the gate A must match the dimension
 * of subsys
 */
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> applyCTRL(const Eigen::MatrixBase<Derived1>& state,
                                             const Eigen::MatrixBase<Derived2>& A,
                                             const std::vector<idx>& ctrl,
                                             const std::vector<idx>& subsys, idx d = 2) {
  const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate = state.derived();
  const dyn_mat<typename Derived1::Scalar>& rA = A.derived();

  // check zero size
  if (!internal::check_nonzero_size(rstate))
    throw exception::ZeroSize("clara::applyCTRL()");

  if (d == 0)
    throw exception::DimsInvalid("clara::applyCTRL()");

  idx N = internal::get_num_subsys(static_cast<idx>(rstate.rows()), d);
  std::vector<idx> dims(N, d);
  return applyCTRL(rstate, rA, ctrl, subsys, dims);
}

/**
 * @brief applies the gate A to the part subys
 * of the multi partite state vector or density matrix
 * @note the dimension of the gate A must match the dimension
 * of a subsys
 * @return gate A applied to the part subsys of state
 */
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> apply(const Eigen::MatrixBase<Derived1>& state,
                                         const Eigen::MatrixBase<Derived2>& A,
                                         const std::vector<idx>& subsys,
                                         const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate = state.derived();
  const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

  // check types
  if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
    throw exception::TypeMismatch("clara::apply()");

  // check zero-size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::apply()");

  // check zero size
  if (!internal::check_nonzero_size(rstate))
    throw exception::ZeroSize("clara::apply()");

  // check square mat
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::apply()");

  // check that dimensional is valid
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::apply()");

  // check subsys is valid w.r.t dims
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::apply()");

  // check that gate matches the dimension of the subsys
  std::vector<idx> subsys_dims(subsys.size());
  for (idx i = 0; i < subsys.size(); ++i)
    subsys_dims[i] = dims[subsys[i]];
  if (!internal::check_dims_match_mat(subsys_dims, rA))
    throw exception::MatrixMismatchSubsys("clara::apply()");

  if (internal::check_cvector(rstate)) {
    if (!internal::check_dims_match_cvect(dims, rstate))
      throw exception::DimsMismatchCvector("clara::apply()");
    return applyCTRL(rstate, rA, {}, subsys, dims);
  } else if (internal::check_square_mat(rstate)) {
    if (!internal::check_dims_match_mat(dims, rstate))
      throw exception::DimsMismatchMatrix("clara::apply()");
    return applyCTRL(rstate, rA, {}, subsys, dims);
  } else
    throw exception::MatrixNotSquareNotCvector("clara::apply()");
}

/**
 * @brief applies the gate A to the part subsys
 * of the multi-partite vector or density matrix state
 * @note the dimension of the gate must match the dimension
 * subsys
 */
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> apply(const Eigen::MatrixBase<Derived1>& state,
                                         const Eigen::MatrixBase<Derived2>& A,
                                         const std::vector<idx>& subsys, idx d = 2) {
  const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate = state.derived();
  const dyn_mat<typename Derived1::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rstate))
    throw exception::ZeroSize("clara::apply()");

  if (d < 2)
    throw exception::DimsInvalid("clara::apply()");

  idx N = internal::get_num_subsys(static_cast<idx>(rstate.rows()), d);
  std::vector<idx> dims(N, d);

  return apply(rstate, rA, subsys, dims);
}

template <typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks) {
  const cmat& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::apply()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::apply()");
  if (Ks.size() == 0)
    throw exception::ZeroSize("clara::apply()");
  if (!internal::check_square_mat(Ks[0]))
    throw exception::MatrixNotSquare("clara::apply()");
  if (Ks[0].rows() != rA.rows())
    throw exception::DimsMismatchMatrix("clara::apply()");

  for (auto&& it : Ks)
    if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
      throw exception::DimsNotEqual("clara::apply()");
  cmat result = cmat::Zero(rA.rows(), rA.rows());

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
  for (idx i = 0; i < Ks.size(); ++i) {
#ifdef WITH_OPENMP_
#pragma omp critical
#endif  // DEBUG
    { result += Ks[i] * rA * adjoint(Ks[i]); }
  }
  return result;
}

/**
 * @brief applies the specified by the set of kraus operators Ks to
 * part subsys of the multi-partite density matrix A
 * @return output density matrix after the action of the channel
 */
template <typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& subsys, const std::vector<idx>& dims) {
  const cmat& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::apply()");
  if (!internal::check_square_mat(rA))
    throw exception::MatrixNotSquare("clara::apply()");
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::apply()");
  if (!internal::check_dims_match_mat(dims, rA))
    throw exception::DimsMismatchMatrix("clara::apply()");

  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::apply()");

  std::vector<idx> subsys_dims(subsys.size());
  for (idx i = 0; i < subsys.size(); ++i)
    subsys_dims[i] = dims[subsys[i]];

  if (Ks.size() == 0)
    throw exception::ZeroSize("clara::apply()");
  if (!internal::check_square_mat(Ks[0]))
    throw exception::MatrixNotSquare("clara::apply()");
  if (!internal::check_dims_match_mat(subsys_dims, Ks[0]))
    throw exception::MatrixMismatchSubsys("clara::apply()");

  for (auto&& it : Ks)
    if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
      throw exception::DimsNotEqual("clara::apply()");

  cmat result = cmat::Zero(rA.rows(), rA.rows());
  for (idx i = 0; i < Ks.size(); ++i)
    result += apply(rA, Ks[i], subsys, dims);
  return result;
}

/**
 * @brief applies the channel specified the set of kraus operators ks to the part
 * subsys of the multi density matrix A
 * @return output density matrix after the action of the channel
 */
template <typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& subsys, idx d = 2) {
  const cmat& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::apply()");
  if (d < 2)
    throw exception::DimsInvalid("clara::apply()");
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);

  return apply(rA, Ks, subsys, dims);
}

/**
 * @brief superoperator matrix
 * construct the superoperator matrix of the channel specified by the rest of kraus operators Ks in
 * the standard operator basis \f$\{|i\rangle\langle j|\}\f$ ordered in lexicographical order
 * @return superoperator matrix
 */
inline cmat kraus2super(const std::vector<cmat>& Ks) {
  if (Ks.size() == 0)
    throw exception::ZeroSize("clara::kraus2super()");
  if (!internal::check_nonzero_size(Ks[0]))
    throw exception::ZeroSize("clara::kraus2super()");
  if (!internal::check_square_mat(Ks[0]))
    throw exception::MatrixNotSquare("kraus2super()");

  for (auto&& it : Ks)
    if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
      throw exception::DimsNotEqual("clara::kraus2super()");

  idx D = static_cast<idx>(Ks[0].rows());
  cmat result(D * D, D * D);
  cmat MN = cmat::Zero(D, D);
  bra A = bra::Zero(D);
  ket B = ket::Zero(D);
  cmat EMN = cmat::Zero(D, D);

#ifdef WITH_OPENMP_
#pragma omp parallel
#endif  // DEBUG
  for (idx m = 0; m < D; ++m) {
    for (idx n = 0; n < D; ++n) {
#ifdef WITH_OPENMP_
#pragma omp critical
#endif  // DEBUG
      {
        MN(m, n) = 1;
        for (idx i = 0; i < Ks.size(); ++i)
          EMN += Ks[i] * MN * adjoint(Ks[i]);
        MN(m, n) = 0;
        for (idx a = 0; a < D; ++a) {
          A(a) = 1;
          for (idx b = 0; b < D; ++b) {
            B(b) = 1;
            result(a * D + b, m * D + n) = static_cast<cmat>(A * EMN * B).value();
            B(b) = 0;
          }
          A(a) = 0;
        }
        EMN = cmat::Zero(D, D);
      }
    }
  }
  return result;
}

/**
 * @brief choi matrix
 * construct the choi matrix of the channel specified by the set of kraus
 * operators Ks in the standard operators basisa \f$\{|i\rangle\langle j|\}\f$
 * ordered in lexicographical order
 * @return choi matrix
 */
inline cmat kraus2choi(const std::vector<cmat>& Ks) {
  if (Ks.size() == 0)
    throw exception::ZeroSize("clara::kraus2choi()");
  if (!internal::check_nonzero_size(Ks[0]))
    throw exception::ZeroSize("clara::kraus2choi()");
  if (!internal::check_square_mat(Ks[0]))
    throw exception::MatrixNotSquare("clara::kraus2choi()");

  for (auto&& it : Ks)
    if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
      throw exception::DimsNotEqual("clara::kraus2choi()");

  idx D = static_cast<idx>(Ks[0].rows());

  /**
   * construct the D x D \sum |jj> vector
   */
  cmat MES = cmat::Zero(D * D, 1);
  for (idx a = 0; a < D; ++a)
    MES(a * D + a) = 1;

  cmat Omega = MES * adjoint(MES);
  cmat result = cmat::Zero(D * D, D * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
  for (idx i = 0; i < Ks.size(); ++i) {
#ifdef WITH_OPENMP_
#pragma omp critical
#endif  // DEBUG
    {
      result +=
          kron(cmat::Identity(D, D), Ks[i]) * Omega * adjoint(kron(cmat::Identity(D, D), Ks[i]));
    }
  }
  return result;
}

/**
 * @brief othogonal kraus operators from choi matrix
 * extract a set of orthogonal kraus operatro from the choi matrix
 * @return set of orthogonal kraus operators
 */
inline std::vector<cmat> choi2kraus(const cmat& A) {
  if (!internal::check_nonzero_size(A))
    throw exception::ZeroSize("clara::choi2kraus()");
  if (!internal::check_square_mat(A))
    throw exception::MatrixNotSquare("clara::choi2kraus()");
  idx D = internal::get_dim_subsystem(A.rows(), 2);
  if (D * D != static_cast<idx>(A.rows()))
    throw exception::DimsInvalid("clara::choi2kraus()");

  dmat ev = hevals(A);
  cmat evec = hevects(A);
  std::vector<cmat> result;

  for (idx i = 0; i < D * D; ++i) {
    if (std::abs(ev(i)) > eps)
      result.push_back(std::sqrt(std::abs(ev(i))) * reshape(evec.col(i), D, D));
  }
  return result;
}

/**
 * @brief convert choi matrix to superoperator matrix
 * @return superoperator matrix
 */
inline cmat choi2super(const cmat& A) {
  if (!internal::check_nonzero_size(A))
    throw exception::ZeroSize("clara::choi2super()");
  if (!internal::check_square_mat(A))
    throw exception::MatrixNotSquare("clara::choi2super()");

  idx D = internal::get_dim_subsystem(static_cast<idx>(A.rows()), 2);
  if (D * D != static_cast<idx>(A.rows()))
    throw exception::DimsInvalid("clara::choi2super()");

  cmat result(D * D, D * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(4)
#endif  // DEBUG
  for (idx a = 0; a < D; ++a)
    for (idx b = 0; b < D; ++b)
      for (idx m = 0; m < D; ++m)
        for (idx n = 0; n < D; ++n)
          result(a * D + b, m * D + n) = A(m * D + a, n * D + b);

  return result;
}

/**
 * @brief convert the superoperator matrix to choi matrix
 * @return choi matrix
 */
inline cmat super2choi(const cmat& A) {
  if (!internal::check_nonzero_size(A))
    throw exception::ZeroSize("clara::super2choi()");
  if (!internal::check_square_mat(A))
    throw exception::MatrixNotSquare("clara::super2choi()");
  idx D = internal::get_dim_subsystem(static_cast<idx>(A.rows()), 2);

  if (D * D != static_cast<idx>(A.rows()))
    throw exception::DimsInvalid("clara::super2choi()");

  cmat result(D * D, D * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(4)
#endif  // DEBUG
  for (idx a = 0; a < D; ++a)
    for (idx b = 0; b < D; ++b)
      for (idx m = 0; m < D; ++m)
        for (idx n = 0; n < D; ++n)
          result(m * D + a, n * D + b) = A(a * D + b, m * D + n);
  return result;
}

/**
 * @brief partial trace
 * partial trace over the first subsystem of bi-partite state vector or density matrix
 * @return partial trace \f$Tr_{A}(\cdot)\f$ over the first subsytem \f$A\f$
 * in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix over the same field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace1(const Eigen::MatrixBase<Derived>& A,
                                          const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::ptrace1()");
  // check that dims is a vlid dimension vector
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::ptrace1()");
  // check dims has only 2 elements
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::ptrace1()");
  idx DA = dims[0];
  idx DB = dims[1];
  dyn_mat<typename Derived::Scalar> result = dyn_mat<typename Derived::Scalar>::Zero(DB, DB);

  if (internal::check_cvector(rA)) {
    if (!internal::check_dims_match_cvect(dims, rA))
      throw exception::DimsMismatchCvector("clara::ptrace1()");

    auto worker = [&](idx i, idx j) noexcept -> typename Derived::Scalar {
      typename Derived::Scalar sum = 0;
      for (idx m = 0; m < DA; ++m)
        sum += rA(m * DB + i) * std::conj(rA(m * DB + j));
      return sum;
    };
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // DEBUG
    for (idx j = 0; j < DB; ++j)
      for (idx i = 0; i < DB; ++i)
        result(i, j) = worker(i, j);
    return result;
  } else if (internal::check_square_mat(rA)) {
    if (!internal::check_dims_match_mat(dims, rA))
      throw exception::DimsMismatchMatrix("clara::ptrace1()");
    auto worker = [=](idx i, idx j) noexcept -> typename Derived::Scalar {
      typename Derived::Scalar sum = 0;
      for (idx m = 0; m < DA; ++m)
        sum += rA(m * DB + i, m * DB + j);
      return sum;
    };
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // DEBUG
    for (idx j = 0; j < DB; ++j)
      for (idx i = 0; i < DB; ++i)
        result(i, j) = worker(i, j);
    return result;
  } else
    throw exception::MatrixNotSquareNotCvector("clara::ptrace1()");
}

/**
 * @brief partial trace
 * partial trace over the first subsystem of bi-partite state vector or density matrix
 * @return partial trace \f$Tr_{A}(\cdot)\f$ over the first subsytem \f$A\f$ in a
 * bi-partite system \f$A\otimes B\f$, as a dynamic matrix over the same scalar.
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace1(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();
  // check zero size
  if (!internal::check_nonzero_size(A))
    throw exception::ZeroSize("clara::ptrace1()");
  if (d == 0)
    throw exception::DimsInvalid("clara::ptrace1()");

  // local dimanesion vector
  std::vector<idx> dims(2, d);
  return ptrace1(A < dims);
}

/**
 * @brief partial trace
 * partial trace over the second subsystem of bi-partite state vector or density matrix
 * @return patrtial trace \f$Tr_{B}(\cdot)\f$ over the second subsytem \f$B\f$
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace2(const Eigen::MatrixBase<Derived>& A,
                                          const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::ptrace2()");
  // check that dims is valid dimension vector
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::ptrace2()");
  // check dims has only 2 elements
  if (dims.size() != 2)
    throw exception::NotBipartite("clara::ptrace2()");
  idx DA = dims[0];
  idx DB = dims[1];

  dyn_mat<typename Derived::Scalar> result = dyn_mat<typename Derived::Scalar>::Zero(DA, DA);

  if (!internal::check_cvector(rA)) {
    if (!internal::check_dims_match_cvect(dims, rA))
      throw exception::DimsMismatchCvector("clara::ptrace2()");

    auto worker = [=](idx i, idx j) noexcept -> typename Derived::Scalar {
      typename Derived::Scalar sum = 0;
      for (idx m = 0; m < DB; ++m)
        sum += rA(i * DB + m) * std::conj(rA(j * DB + m));
      return sum;
    };
#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // DEBUG
    for (idx j = 0; j < DA; ++j)
      for (idx i = 0; i < DA; ++i)
        result(i, j) = worker(i, j);
    return result;
  } else if (internal::check_square_mat(rA)) {
    if (!internal::check_dims_match_mat(dims, rA))
      throw exception::DimsMismatchMatrix("clara::ptrace2()");

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // DEBUG
    for (idx j = 0; j < DA; ++j)
      for (idx i = 0; i < DA; ++i)
        result(i, j) = trace(rA.block(i * DB, j * DB, DB, DB));
    return result;
  } else
    throw exception::MatrixNotSquareNotCvector("clara::ptrace1()");
}

/**
 * @brief partial trace
 * partial trace over the second subsystem
 * of bi-partite state vector or density matrix
 * @return partial trace \f$Tr_{B}(\cdot)\f$ over the second subsytem \f$B\f$
 * in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix over the same scalar field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace2(const Eigen::MatrixBase<Derived>& A, idx d = 2) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(A))
    throw exception::ZeroSize("clara::ptrace2()");
  if (d == 0)
    throw exception::DimsInvalid("clara::ptrace2()");
  std::vector<idx> dims(2, d);
  return ptrace2(A, dims);
}

/**
 * @brief partial trace
 * partial trace of the multi-partite state vector or density matrix
 * over list of subsystem
 * @return partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsytems \a subsys
 * in multi-partit system, as a dynamic matrix
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& subsys,
                                         const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::ptrace()");

  // check that dims is a vlid dimension vector
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::ptrace()");

  // check that subsys are valid
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::ptrace()");

  idx D = static_cast<idx>(rA.rows());
  idx N = dims.size();
  idx Nsubsys = subsys.size();
  idx Nsubsys_bar = N - Nsubsys;
  idx Dsubsys = 1;
  for (idx i = 0; i < Nsubsys; ++i)
    Dsubsys *= dims[subsys[i]];
  idx Dsubsys_bar = D / Dsubsys;

  idx Cdims[maxn];
  idx Csubsys[maxn];
  idx Cdimssubsys[maxn];
  idx Csubsys_bar[maxn];
  idx Cdimssubsys_bar[maxn];

  idx Cmidxcolsubsys_bar[maxn];

  std::vector<idx> subsys_bar = complement(subsys, N);
  std::copy(std::begin(subsys_bar), std::end(subsys_bar), std::begin(Csubsys_bar));

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

  dyn_mat<typename Derived::Scalar> result =
      dyn_mat<typename Derived::Scalar>(Dsubsys_bar, Dsubsys_bar);

  if (internal::check_cvector(rA)) {
    if (!internal::check_dims_match_cvect(dims, rA))
      throw exception::DimsMismatchCvector("clara::ptrace()");
    if (subsys.size() == dims.size()) {
      result(0, 0) = (adjoint(rA) * rA).value();
      return result;
    }
    if (subsys.size() == 0)
      return rA * adjoint(rA);

    auto worker = [=, &Cmidxcolsubsys_bar](idx i) noexcept -> typename Derived::Scalar {
      // use static allocation for speed
      idx Cmidxrow[maxn];
      idx Cmidxcol[maxn];
      idx Cmidxrowsubsys_bar[maxn];
      idx Cmidxsubsys[maxn];

      // get the row multi-indexes of the complement
      internal::n2multiidx(i, Nsubsys_bar, Cdimssubsys_bar, Cmidxrowsubsys_bar);

      for (idx k = 0; k < Nsubsys_bar; ++k) {
        Cmidxrow[Csubsys_bar[k]] = Cmidxrowsubsys_bar[k];
        Cmidxcol[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
      }
      typename Derived::Scalar sm = 0;
      for (idx a = 0; a < Dsubsys; ++a) {
        // get the multi-index over which we do the summation
        internal::n2multiidx(a, Nsubsys, Cdimssubsys, Cmidxsubsys);
        for (idx k = 0; k < Nsubsys; ++k)
          Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]] = Cmidxsubsys[k];

        sm += rA(internal::multiidx2n(Cmidxrow, N, Cdims)) *
              std::conj(rA(internal::multiidx2n(Cmidxcol, N, Cdims)));
      }
      return sm;
    };
    for (idx j = 0; j < Dsubsys_bar; ++j) {
      internal::n2multiidx(j, Nsubsys_bar, Cdimssubsys_bar, Cmidxcolsubsys_bar);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
      for (idx i = 0; i < Dsubsys_bar; ++i) {
        result(i, j) = worker(i);
      }
    }
    return result;
  } else if (internal::check_square_mat(rA)) {
    if (!internal::check_dims_match_mat(dims, rA))
      throw exception::DimsMismatchMatrix("clara::ptrace()");

    if (subsys.size() == dims.size()) {
      result(0, 0) = rA.trace();
      return result;
    }
    if (subsys.size() == 0)
      return rA;

    auto worker = [=, &Cmidxcolsubsys_bar](idx i) noexcept -> typename Derived::Scalar {
      idx Cmidxrow[maxn];
      idx Cmidxcol[maxn];
      idx Cmidxrowsubsys_bar[maxn];
      idx Cmidxsubsys[maxn];

      internal::n2multiidx(i, Nsubsys_bar, Cdimssubsys_bar, Cmidxrowsubsys_bar);
      for (idx k = 0; k < Nsubsys_bar; ++k) {
        Cmidxrow[Csubsys_bar[k]] = Cmidxrowsubsys_bar[k];
        Cmidxcol[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
      }
      typename Derived::Scalar sm = 0;
      for (idx a = 0; a < Dsubsys; ++a) {
        internal::n2multiidx(a, Nsubsys, Cdimssubsys, Cmidxsubsys);
        for (idx k = 0; k < Nsubsys; ++k)
          Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]] = Cmidxsubsys[k];
        sm +=
            rA(internal::multiidx2n(Cmidxrow, N, Cdims), internal::multiidx2n(Cmidxcol, N, Cdims));
      }
      return sm;
    };

    for (idx j = 0; j < Dsubsys_bar; ++j) {
      internal::n2multiidx(j, Nsubsys_bar, Cdimssubsys_bar, Cmidxcolsubsys_bar);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
      for (idx i = 0; i < Dsubsys_bar; ++i) {
        result(i, j) = worker(i);
      }
    }
    return result;
  } else
    throw exception::MatrixNotSquareNotCvector("clara::ptrace()");
}

/**
 * @brief atrace of the multi-partite state vector or density matrix
 * over the list of subsystem
 * @return partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsytems \a subsys
 * in a multi-partite system, as a dynamic matrix over the same
 * scalar field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& subsys, idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::ptrace()");
  if (d < 2)
    throw exception::DimsInvalid("clara::ptrace()");
  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return ptrace(rA, subsys, dims);
}

/**
 * @brief partial transpose
 * partial transpose of the multi-partite state vector or density matrix over
 * a list of subsystem
 * @return partial transpose \f$(\cdot)^{T_{subsys}}\f$ over the subsystem subsys in a
 * multi-partite system, as a dynamic matrix over the same scalar field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptranspose(const Eigen::MatrixBase<Derived>& A,
                                             const std::vector<idx>& subsys,
                                             const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::ptranspose()");
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::ptranspose()");
  if (!internal::check_subsys_match_dims(subsys, dims))
    throw exception::SubsysMismatchdims("clara::ptranspose()");

  idx D = static_cast<idx>(rA.rows());
  idx N = dims.size();
  idx Nsubsys = subsys.size();
  idx Cdims[maxn];
  idx Cmidxcol[maxn];
  idx Csubsys[maxn];

  for (idx i = 0; i < N; ++i)
    Cdims[i] = dims[i];
  for (idx i = 0; i < Nsubsys; ++i)
    Csubsys[i] = subsys[i];

  dyn_mat<typename Derived::Scalar> result(D, D);

  if (internal::check_cvector(rA)) {
    if (!internal::check_dims_match_cvect(dims, rA))
      throw exception::DimsMismatchCvector("clara::ptranspose()");
    if (subsys.size() == dims.size())
      return (rA * adjoint(rA)).transpose();

    if (subsys.size() == 0)
      return rA * adjoint(rA);

    auto worker = [=, &Cmidxcol](idx i) noexcept -> typename Derived::Scalar {
      idx midxcoltmp[maxn];
      idx midxrow[maxn];

      for (idx k = 0; k < N; ++k)
        midxcoltmp[k] = Cmidxcol[k];

      internal::n2multiidx(i, N, Cdims, midxrow);

      return rA(internal::multiidx2n(midxrow, N, Cdims)) *
             std::conj(rA(internal::multiidx2n(midxcoltmp, N, Cdims)));
    };
    for (idx j = 0; j < D; ++j) {
      internal::n2multiidx(j, N, Cdims, Cmidxcol);

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
      for (idx i = 0; i < D; ++i)
        result(i, j) = worker(i);
    }
    return result;
  } else if (internal::check_square_mat(rA)) {
    if (!internal::check_dims_match_mat(dims, rA))
      throw exception::DimsMismatchCvector("clara::ptranspose()");
    if (subsys.size() == dims.size())
      return (rA * adjoiint(rA)).transpose();
    if (subsys.size() == 0)
      return rA * adjoint(rA);

    auto worker = [=, &Cmidxcol](idx i) noexcept -> typename Derived::Scalar {
      idx midxcoltmp[maxn];
      idx midxrow[maxn];

      for (idx k = 0; k < N; ++k)
        midxcoltmp[k] = Cmidxcol[k];
      internal::n2multiidx(i, N, Cdims, midxrow);

      for (idx k = 0; k < Nsubsys; ++k)
        std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);

      // write the result
      return rA(internal::multiidx2n(midxrow, N, Cdims),
                internal::multiidx2n(midxcoltmp, N, Cdims));
    };

    for (idx j = 0; j < D; ++j) {
      // compute the column multi-index
      internal::n2multiidx(j, N, Cdims, Cmidxcol);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
      for (idx i = 0; i < D; ++i)
        result(i, j) = worker(i);
    }
    return result;
  } else
    throw exception::MatrixNotSquareNotCvector("clara::ptranspose()");
}

/**
 * @brief partitial transpose
 * partial transpose of multi-partite state vector or density matrix
 * over a list of subsystem
 * @return partial transpose \f$(\cdot)^{T_{subsys}}\f$ over the subsystem
 * in a multi-partite system, as dynamic matrix over the same scalar field A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptranspose(const Eigen::MatrixBase<Derived>& A,
                                             const std::vector<idx>& subsys, idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::ptranspose()");

  if (d < 2)
    throw exception::DimsInvalid("clara::ptranspose()");

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);
  return ptranspose(rA, subsys, dims);
}

/**
 * @brief susbsystem permutation
 * permute the subsystem of state vector or density matrix. the qubit
 * perm[i] is permuted to the location to i
 * @return permuted system, as a dynamic matrix over the same scalar field as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> syspermute(const Eigen::MatrixBase<Derived>& A,
                                             const std::vector<idx>& perm,
                                             const std::vector<idx>& dims) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::syspermute()");
  // check that dims is a vlid dimension vector
  if (!internal::check_dims(dims))
    throw exception::DimsInvalid("clara::syspermute()");
  // check that we a valid permutation
  if (!internal::check_perm(perm))
    throw exception::PermInvalid("clara::syspermute()");
  // check that permutation match dimension
  if (perm.size() != dims.size())
    throw exception::PermMismatchDims("clara::syspermute()");

  idx D = static_cast<idx>(rA.rows());
  idx N = dims.size();

  dyn_mat<typename Derived::Scalar> result;

  if (internal::check_cvector(rA)) {
    idx Cdims[maxn];
    idx Cperm[maxn];

    // check dims match of the dimension rA
    if (!internal::check_dims_match_cvect(dims, rA))
      throw exception::DimsMismatchCvector("clara::syspermute()");

    // copy dims in Cdims and perm in Cperm
    for (idx i = 0; i < N; ++i) {
      Cdims[i] = dims[i];
      Cperm[i] = perm[i];
    }
    result.resize(D, 1);

    auto worker = [&Cdims, &Cperm, N](idx i) noexcept -> idx {
      /**
       * use static allocation for speeding
       * double the size for matrices reshaped as vectors
       */
      idx midx[maxn];
      idx midxtmp[maxn];
      idx perdims[maxn];

      // compute the multi-index
      internal::n2multiidx(i, N, Cdims, midx);
      for (idx k = 0; k < N; ++k) {
        perdims[k] = Cdims[Cperm[k]];
        midxtmp[k] = midx[Cperm[k]];
      }
      return internal::multiidx2n(midxtmp, N, perdims);
    };
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
    for (idx i = 0; i < D; ++i)
      result(worker(i)) = rA(i);
    return result;
  } else if (internal::check_square_mat(rA)) {
    idx Cdims[2 * maxn];
    idx Cperm[2 * maxn];

    // check that dims match the dimnesion of rA
    if (!internal::check_dims_match_mat(dims, rA))
      throw exception::DimsMismatchMatrix("clara::syspermute()");
    for (idx i = 0; i < N; ++i) {
      Cdims[i] = dims[i];
      Cdims[i + N] = dims[i];
      Cperm[i] = perm[i];
      Cperm[i + N] = perm[i] + N;
    }
    result.resiz(D * D, 1);
    // map A to column vector
    dyn_mat<typename Derived::Scalar> vectA = Eigen::Map<dyn_mat<typename Derived::Scalar>>(
        const_cast<typename Derived::Scalar*>(rA.data()), D * D, 1);

    auto worker = [&Cdims, &Cperm, N](idx i) noexcept -> idx {
      /**
       * use static allocation for speed
       * double the size for matrices reshaped vectros
       */
      idx midx[2 * maxn];
      idx midxtmp[2 * maxn];
      idx permdims[2 * maxn];

      internal::n2multiidx(i, 2 * N, Cdims, midx);

      for (idx k = 0; k < 2 * N; ++k) {
        permdims[k] = Cdims[Cperm[k]];
        midxtmp[k] = midx[Cperm[k]];
      }
      return internal::multiidx2n(midxtmp, 2 * N, permdims);
    };
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif  // DEBUG
    for (idx i = 0; i < D * D; ++i)
      result(worker(i)) = rA(i);
    return reshape(result, D, D);
  } else
    throw exception::MatrixNotSquareNotCvector("clara::syspermute()");
}

/**
 * @brief subsystem permutation
 * permutate subsystem of a state vector or density matrix
 * the quibt perm[i] is permuted to the location
 * @return permuted system, as dynamic matrix over the same scalar fiedl as A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> syspermute(const Eigen::MatrixBase<Derived>& A,
                                             const std::vector<idx>& perm, idx d = 2) {
  const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::syspermute()");

  // check valid dims
  if (d < 2)
    throw exception::DimsInvalid("clara::syspermute()");

  idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
  std::vector<idx> dims(N, d);

  return syspermute(rA, perm, dims);
}

}  // namespace clara

#endif  // !OPERATIONS_H_
