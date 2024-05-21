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
/**
 * @brief convert a linear index to multi-dimensional indices for a given set of dimension
 * @param n the linear index
 * @param numdims the number of dimension
 * @param dims pointer to the array of dimension
 * @param result pointer to the array to store the multi-dimensional indices
 * NOTE: this function is used to convert a linear index to multi-dimensional indices corresponding
 *        to a multi-dimensional array with the given dimensions. the result array will store the
 *        multi-dimensional indices in reverse order, with the first index corresponding to the last
 *        dimension and so on. for example, for a 3x3x3 array, the linear index 5 would be converted
 * to multidimensional indices [2, 1, 0]
 */
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
 * @brief convert multi-dimensional indices to a linear indexx
 * @param midx pointer to the array of multi-dimensional indces
 * @param numdims the number of dimensions
 * @param dims pointer to the array of dimensions
 * @return the linear index corresponding the multi-dimensional
 *
 * NOTE: this function is used to convert multi-dimensional indices to a linear index for a
 * multi-dimensional mutli-dimensional array with the given dimensions. the midx array array is
 * expected to store the multi-dimensional indices in reverse order, with the first index
 * corresponding to the last dimensions and so on. for eample a 3x3x3 array, the multi-dimensional
 * indices [2, 1, 0] would be converted to the linear index 5
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

/**
 * @brief check if a matrix is square
 * @tparam derived Eigen matrix type
 * @param A eigen matrix to check
 * @return True if the matrix is square, false otherwise
 */
template <typename Derived>
bool check_square_mat(const Eigen::MatrixBase<Derived>& A) {
  return A.rows() == A.cols();
}

/**
 * @brief check if a matrix is a vector
 * @tparam Derived Eigen matrix type
 * @param A Eigen matrix to check
 * @return True if the matrix is a vector, false otherwise
 */
template <typename Derived>
bool check_vector(const Eigen::MatrixBase<Derived>& A) {
  return A.rows() == 1 || A.cols() == 1;
}

/**
 * @brief check if a matrix is a row vector
 * @tparam Derived Eigen matrix type
 * @param A Eigen matrix to check
 * @return true if type matrix is a row vector, false otherwise
 */
template <typename Derived>
bool check_rvector(const Eigen::MatrixBase<Derived>& A) {
  return A.rows() == 1;
}

/**
 * @brief check if a matrix is a column vector
 * @tparam Derived Eigen matrix type
 * @param A eigen matrix to check
 * @return true if the matrix is a column vector, false otherwise
 */
template <typename Derived>
bool check_cvector(const Eigen::MatrixBase<Derived>& A) {
  return A.cols() == 1;
}

/**
 * @brief check if an object has a non-zero size
 * @tparam T type of the object
 * @param x object to check
 * @return true if the object has a non-zero size, false otherwise
 */
template <typename T>
bool check_nonzero_size(const T& x) noexcept {
  return x.size() != 0;
}

/**
 * @brief check if two objects have matching size
 * @tparam T1 type of the first object
 * @tparam T2 type of the second object
 * @param lhs the first object to compare
 * @param rhs the second object to compare
 * @return true if the size matrch, false otherwise
 */
template <typename T1, typename T2>
bool check_matching_sizes(const T1& lhs, const T2& rhs) noexcept {
  return lhs.size() == rhs.size();
}

/**
 * @brief check if a vector represent valid dimension
 * @param dims vector of dimension to check
 * @return true if dimension are valid, false otherwise
 */
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
 * @brief check if valid dims match the dimension of a valid (non-zero sized) square matrix
 * @param dims the vector of dimension
 * @param A the square matrix
 * @return true if valid dims match the dimension of the matrix, false otherwise
 *
 * NOTE: this function check whether the product of all dimensions in dims is equal to the number
 *        of rows/columns in the square matrix A, the dimensions in dims are assumed to be far valid
 * square matrix of the same dimension
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

/**
 * @brief check if valid dimension match dimension of a valid column vector
 * @tparam Derived Eigen matrix type
 * @param dims vector of dimensions to check
 * @param A eigen matrix (column vector) to compare dimension
 * @return true if dimension match false otherwise
 *
 * NOTE: this function check if the product of dimension in `dims` matches the number
 * of rows in the column vector `A`.
 */
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

/**
 * @brief check if valid dimensions mtach the dimension of valid row vector
 * @tparam Derived Eigen matrix type
 * @param dims vector of dimension to checks
 * @param A Eigen matrix to compare dimension
 * @return true if dimension match, false otherwise
 *
 * NOTE: function check if the product of dimension in `dims` matches the number of columns in the
 * row vector `A`.
 */
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

/**
 * @brief check if all elements in valid dimension vector are equal to a given dimension
 * @param dims vector of dimensions to check
 * @param dim dimension to compare
 * @return true if all dimension in `dims` are equal to `dim`, false otherwise
 */
inline bool check_eq_dims(const std::vector<idx>& dims, idx dim) noexcept {
#ifndef NDEBUG
  assert(dims.size() > 0);
#endif  // !NDEBUG

  for (idx i : dims)
    if (i != dim)
      return false;
  return true;
}

/**
 * @brief check if subsystem indices are valid with respect to given diemension
 * @param subsys vector of subsystem indices to check
 * @param dims vector a valid dimension
 * @param dims vector of valid dimensions
 * @return true if subsystem indices are valid, false otherwise
 *
 * NOTE: the function check if each index in `subsys` is within the valid range of dimension
 * specified by `dims`
 */
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

/**
 * @brief check if matrix is 2x2 matrix
 * @tparam Derived Eigen matrix type
 * @param A Eigen matrix to check
 * @return true if the matrix is 2x2 matrix, false otherwise
 */
template <typename Derived>
bool check_qubit_matrix(const Eigen::MatrixBase<Derived>& A) noexcept {
  return A.rows() == 2 && A.cols() == 2;
}

/**
 * @brief check if collumn vector is a 2x1 vector
 * @tparam Derived Eigen matrix type
 * @param A Eigen column vector to check
 * @return true if the column vector is a 2x1 vector, false otherwise
 */
template <typename Derived>
bool check_qubit_cvector(const Eigen::MatrixBase<Derived>& A) noexcept {
  return A.rows() == 2 && A.cols() == 1;
}

/**
 * @brief check if row vector is a 1x2 vector
 * @tparam Derived Eigen matrix type
 * @param A eigen row vector to check
 * @return true if the row vector is a 1x2 vector, false otherwise
 */
template <typename Derived>
bool check_qubit_rvector(const Eigen::MatrixBase<Derived>& A) noexcept {
  return A.rows() == 1 && A.cols() == 2;
}

/**
 * @brief check if vector is a 1x2 or 2x1 vector
 * @tparam Derived Eigen matrix type
 * @param A Eigen vector to check
 * @return true if the vector is a 1x2 or 2x1 vector, false otherwise
 */
template <typename Derived>
bool check_qubit_vector(const Eigen::MatrixBase<Derived>& A) noexcept {
  return (A.rows() == 1 && A.cols() == 2) || (A.rows() == 2 && A.cols() == 1);
}

/**
 * @brief check if a given vector represent a valid permutation
 * @param perm vector repersenting a permutation
 * @return true if the vector is a valid permutation, false otherwise
 *
 * NOTE: this function check if the input vector contains a valid permutation of indices from 0 to
 * (size - 1)
 */
inline bool check_perm(const std::vector<idx>& perm) {
  if (perm.size() == 0)
    return false;
  std::vector<idx> ordered(perm.size());
  std::iota(std::begin(ordered), std::end(ordered), 0);

  return std::is_permutation(std::begin(ordered), std::end(ordered), std::begin(perm));
}

/**
 * @brief kronecker product of two matrices
 * @tparam Derived1 type of the first matrix
 * @tparam Derived2 type of the second matrix
 * @param A the first matrix
 * @param B the second matrix
 * @return the kronecker product of A and B
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

/**
 * @brief Direct sum of two matrices
 * @tparam Derived1 type of the first matrix
 * @tparam Derived2 type
 * @param A the first matrix
 * @param B the second matrix
 * @return the direct sum of A and B
 */
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

/**
 * @brief extract a variadic template argument pack and emplaces them into a std::vector
 * @tparam T type of the elements to be emplaced in the vector
 * @param v the vector to emplace elements
 * @param first the first element to emplace
 * @param args the rest of the elements to emplace
 *
 * NOTE: this function recursively emplace the provide elements into the vector
 */
template <typename T>
void variadict_vector_emplace(std::vector<T>&) {}
template <typename T, typename First, typename... Args>
void variadic_vector_emplace(std::vector<T>& v, First&& first, Args&&... args) {
  v.emplace_back(std::forward<First>(first));
  variadic_vector_emplace(v, std::forward<Args>(args)...);
}

/**
 * @brief get the number of subsystem from an object of size sz
 * @param sz the size of the object
 * @param d the dimension of each subsystem
 * @return the number of subsystem
 *
 * NOTE: this function calculates the number of subsystem in an boject (ket/bra/density matrix) of
 *       size, where each subsystem has the same dimension d
 */
inline idx get_num_subsys(idx sz, idx d) {
#ifndef NDEBUG
  assert(sz > 0);
  assert(d > 1);
#endif  // !NDEBUG
  return static_cast<idx>(std::llround(std::log2(sz) / std::log2(d)));
}

/**
 * @brief get the dimension of a sbusystem from an object of size sz consiting of N subsystem
 * @param sz the size of the object
 * @param N the number of subsystem
 * @return the dimension of each subystem
 *
 * NOTE: this function calculates the dimension of each subsystem in an object (ket/bra/density
 * matrix) of size sz, consiting of N subsystem. the dimension in the object are assumed to be the
 * same
 */
inline idx get_dim_subsystem(idx sz, idx N) {
#ifndef NDEBUG
  assert(N > 0);
  assert(sz > 0);
#endif  // !NDEBUG
  if (N == 2)
    return static_cast<idx>(std::llround(std::sqrt(sz)));
  return static_cast<idx>(std::llround(std::pow(sz, 1. / N)));
}

/**
 * @brief implementation of displaying matrix-line object with formatting options
 *
 * NOTE: this struct provides a methods to display a matrix-like object with customizable formatting
 */
struct Display_Impl_ {
  /**
   * @brief display a matrix-like object with customizable formatting
   * @tparam T type of the matrix-like object
   * @param A the matrix-like object to be displayed
   * @param os the output stream to which the formatted matrix will be written
   * @param chop the threshold for chopping small values
   * @return the modified output stream
   */
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
