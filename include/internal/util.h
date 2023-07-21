#ifndef INTERNAL_UTIL_H_
#define INTERNAL_UTIL_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <type_traits>

#include "../constants.h"
#include "../types.h"

#if (__GNUC__ && !__clang__)
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif  // (__GNUC__ && !__clang__)

namespace clara {
namespace internal {
inline void n2multiidx(idx n, idx numdims, const idx* const dims, idx* result) noexcept {
#ifndef NDEBUG
  if (numdims > 0) {
    idx D = 1;
    for (idx i = 0; i < numdims; ++i)
      D *= dims[i];
    assert(n < D);
  }
#endif  // !NDEBUG
  // no error check in release version to improve speed
  for (idx i = 0; i < numdims; ++i) {
    result[numdims - i - 1] = n % (dims[numdims - i - 1]);
    n /= (dims[numdims - i - 1]);
  }
}

#if (__GNUC__ && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
/**
 * multiindex to integer
 * standard lexicographical order
 */
inline idx multiidx2n(const idx* const midx, idx numdims, const idx* const dims) {
#ifndef NDEBUG
  assert(numdims > 0);
#endif  // !NDEBUG
  idx part_prod[2 * maxn];
  idx result = 0;
  part_prod[numdims - 1] = 1;
  for (idx i = 1; i < numdims; i++) {
    part_prod[numdims - i - 1] = part_prod[numdims - i] * dims[numdims - i];
    result += midx[numdims - i - 1] * part_prod[numdims - i - 1];
  }
  return result + midx[numdims - 1];
}

#if (__GNUC__ && !__clang__)
#pragma GCC diagnostic pop
#endif

// check square matrix
template <typename Derived>
bool check_square_mat(const Eigen::MatrixBase<Derived>& A) {
  return A.rows() == A.cols();
}

// check whether intpu is a vector or not
template <typename Derived>
bool check_vector(const Eigen::MatrixBase<Derived>& A) {
  return A.rows() == 1 || A.cols() == 1;
}

// check wheter input is a column vector or not
template <typename Derived>
bool check_rvector(const Eigen::MatrixBase<Derived>& A) {
  return A.rows() == 1;
}

// check wheter input is a column vector or not
template <typename Derived>
bool check_cvector(const Eigen::MatrixBase<Derived>& A) {
  return A.cols() == 1;
}

// check non-zero size of object that support size() function
template <typename T>
bool check_nonzero_size(const T& x) noexcept {
  return x.size() != 0;
}

// check that all size match
template <typename T1, typename T2>
bool check_matching_sizes(const T1& lhs, const T2& rhs) noexcept {
  return lhs.size() == rhs.size();
}

// check that dims is a valid dimension vector
inline bool check_dims(const std::vector<idx>& dims) {
  if (dims.size() == 0)
    return false;
  return std::find_if(std::begin(dims), std::end(dims), [dims](idx i) -> bool {
           if (i == 0)
             return true;
           else
             return false;
         }) == std::end(dims);
}

/**
 * check that valid dims match the dimension
 * of valid (non-zero sized) square matrix
 */
template <typename Derived>
bool check_dims_match_mat(const std::vector<idx>& dims, const Eigen::MatrixBase<Derived>& A) {
#ifndef NDEBUG
  assert(dims.size() > 0);
  assert(A.rows() == A.cols());
#endif  // !NDEBUG
  idx proddim = std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                                std::multiplies<idx>());
  return proddim == static_cast<idx>(A.rows());
}

// check that valid dims match the dimension of valid column vector
template <typename Derived>
bool check_dims_match_cvect(const std::vector<idx>& dims, const Eigen::MatrixBase<Derived>& A) {
#ifndef NDEBUG
  assert(dims.size() > 0);
  assert(A.rows() > 0);
  assert(A.cols() == 1);
#endif  // !NDEBUG
  idx proddim = std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                                std::multiplies<idx>());
  return proddim == static_cast<idx>(A.rows());
}

// check that valid dims mtch the dimension of valid row vector
template <typename Derived>
bool check_dims_match_rvect(const std::vector<idx>& dims, const Eigen::MatrixBase<Derived>& A) {
#ifndef NDEBUG
  assert(dims.size() > 0);
  assert(A.cols() > 0);
  assert(A.rows() == 1);
#endif  // !NDEBUG
  idx proddim = std::accumulate(std::begin(dims), std::end(dims), static_cast<idx>(1),
                                std::multiplies<idx>());
  return proddim == static_cast<idx>(A.cols());
}

// check that all elements in valid dims equal to dim
inline bool check_eq_dims(const std::vector<idx>& dims, idx dim) noexcept {
#ifndef NDEBUG
  assert(dims.size() > 0);
#endif  // !NDEBUG

  for (idx i : dims)
    if (i != dim)
      return false;
  return true;
}

// check that subsys is valid with respect to valid dims
inline bool check_subsys_match_dims(const std::vector<idx>& subsys, const std::vector<idx>& dims) {
  if (subsys.size() > dims.size())
    return false;
  // sort the subsys
  std::vector<idx> subsyssort = subsys;
  std::sort(std::begin(subsyssort), std::end(subsyssort));

  // check duplicates
  if (std::unique(std::begin(subsyssort), std::end(subsyssort)) != std::end(subsyssort))
    return false;
  // check range of subsystem
  return std::find_if(std::begin(subsyssort), std::end(subsyssort), [dims](idx i) -> bool {
           return i > dims.size() - 1;
         }) == std::end(subsyssort);
}

// check matrix is 2 x 2
template <typename Derived>
bool check_qubit_matrix(const Eigen::MatrixBase<Derived>& A) noexcept {
  return A.rows() == 2 && A.cols() == 2;
}

// check column vector is 2 x 1
template <typename Derived>
bool check_qubit_cvector(const Eigen::MatrixBase<Derived>& A) noexcept {
  return A.rows() == 2 && A.cols() == 1;
}

// check row vector 1 x 2
template <typename Derived>
bool check_qubit_rvector(const Eigen::MatrixBase<Derived>& A) noexcept {
  return A.rows() == 1 && A.cols() == 2;
}

// check row vector is 1 x 2 or 2 x 1
template <typename Derived>
bool check_qubit_vector(const Eigen::MatrixBase<Derived>& A) noexcept {
  return (A.rows() == 1 && A.cols() == 2) || (A.rows() == 2 && A.cols() == 1);
}

// check valid permutation
inline bool check_perm(const std::vector<idx>& perm) {
  if (perm.size() == 0)
    return false;
  std::vector<idx> ordered(perm.size());
  std::iota(std::begin(ordered), std::end(ordered), 0);

  return std::is_permutation(std::begin(ordered), std::end(ordered), std::begin(perm));
}

/**
 * kronecker product of 2 matrices, preserve return type
 * internal function for the variadic tempalte function wrapper kron()
 */
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> kron2(const Eigen::MatrixBase<Derived1>& A,
                                         const Eigen::MatrixBase<Derived2>& B) {
  const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
  const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

  if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
    throw exception::TypeMismatch("clara::kron()");

  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::kron()");

  // check zero size
  if (!internal::check_nonzero_size(rB))
    throw exception::ZeroSize("clara::kron()");

  idx Acols = static_cast<idx>(rA.cols());
  idx Arows = static_cast<idx>(rA.rows());
  idx Bcols = static_cast<idx>(rB.cols());
  idx Brows = static_cast<idx>(rB.rows());

  dyn_mat<typename Derived1::Scalar> result;
  result.resize(Arows * Brows, Acols * Bcols);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif  // DEBUG
  for (idx j = 0; j < Acols; ++j)
    for (idx i = 0; i < Arows; ++i)
      result.block(i * Brows, j * Bcols, Brows, Bcols) = rA(i, j) * rB;
  return result;
}

template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> dirsum2(const Eigen::MatrixBase<Derived1>& A,
                                           const Eigen::MatrixBase<Derived2>& B) {
  const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
  const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

  if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
    throw exception::ZeroSize("clara::dirsum()");
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::dirsum()");
  if (!internal::check_nonzero_size(rB))
    throw exception::ZeroSize("clara::dirsum()");
  idx Acols = static_cast<idx>(rA.cols());
  idx Arows = static_cast<idx>(rA.rows());
  idx Bcols = static_cast<idx>(rB.cols());
  idx Brows = static_cast<idx>(rB.rows());

  dyn_mat<typename Derived1::Scalar> result =
      dyn_mat<typename Derived1::Scalar>::Zero(Arows + Brows, Acols + Bcols);
  result.block(0, 0, Arows, Acols) = rA;
  result.block(Arows, Acols, Brows, Bcols) = rB;
  return result;
}

// extract variadic template argyment pack into std::vector
template <typename T>
void variadict_vector_emplace(std::vector<T>&) {}

template <typename T, typename First, typename... Args>
void variadic_vector_emplace(std::vector<T>& v, First&& first, Args&&... args) {
  v.emplace_back(std::forward<First>(first));
  variadic_vector_emplace(v, std::forward<Args>(args)...);
}

/**
 * return the number of subystem (each subsystem assumed of the same
 * dimension d) from an object (ket/bra/density matrix) of size sz
 */
inline idx get_num_subsys(idx sz, idx d) {
#ifndef NDEBUG
  assert(sz > 0);
  assert(d > 1);
#endif  // !NDEBUG
  return static_cast<idx>(std::llround(std::log2(sz) / std::log2(d)));
}

/**
 * return the dimension of a subsystem (each subsystem assumed of the
 * dimension d) from an object (ket/bra/density matrix) of size sz
 * consiting of N subsystem
 */
inline idx get_dim_subsystem(idx sz, idx N) {
  #ifndef NDEBUG
  assert(N > 0);
  assert(sz > 0);
  #endif // !NDEBUG
  if (N == 2)
    return static_cast<idx>(std::llround(std::sqrt(sz)));
  return static_cast<idx>(std::llround(std::pow(sz, 1. / N)));
}

// implementation detasil for pretty formating
struct Display_Impl_ {
  template <typename T>
  // T must support rows(), cols(), operator()(idx, idx) const
  std::ostream& display_impl_(const T& A, std::ostream& os, double chop = clara::chop) const {
    std::ostringstream ostr;
    ostr.copyfmt(os);
    std::vector<std::string> vstr;
    std::string strA;

    for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
      for (idx j = 0; j < static_cast<idx>(A.cols()); ++j) {
        strA.clear();
        ostr.clear();
        ostr.str(std::string{});

        // convert to complex
        double re = static_cast<cplx>(A(i, j)).real();
        double im = static_cast<cplx>(A(i, j)).imag();
        if (std::abs(re) < chop && std::abs(im) < chop) {
          ostr << "0 ";
          vstr.push_back(ostr.str());
        } else if (std::abs(re) < chop) {
          ostr << im;
          vstr.push_back(ostr.str() + "i");
        } else if (std::abs(im) < chop) {
          ostr << re;
          vstr.push_back(ostr.str() + " ");
        } else {
          ostr << re;
          strA = ostr.str();
          strA += (im > 0 ? " + " : " - ");
          ostr.clear();
          ostr.str(std::string());
          ostr << std::abs(im);
          strA += ostr.str();
          strA += "i";
          vstr.push_back(strA);
        }
      }
    }
    // determine the maximum length of the entries in each column
    std::vector<idx> maxlengthcols(A.cols(), 0);
    for (idx i = 0; i < static_cast<idx>(A.rows()); i++)
      for (idx j = 0; j < static_cast<idx>(A.cols()); ++j)
        if (vstr[i * A.cols() + j].size() > maxlengthcols[j])
          maxlengthcols[j] = vstr[i * A.cols() + j].size();

    for (idx i = 0; i < static_cast<idx>(A.rows()); i++) {
      os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right << vstr[i * A.cols()];
      for (idx j = 1; j < static_cast<idx>(A.cols()); ++j)
        os << std::setw(static_cast<int>(maxlengthcols[j] + 2)) << std::right
           << vstr[i * A.cols() + j];

      if (i < static_cast<idx>(A.rows()) - 1)
        os << std::endl;
    }
    return os;
  }
};

}  // namespace internal
}  // namespace clara
#endif  // !INTERNAL_UTIL_H_
