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

#ifndef INPUT_OUTPUT_H_
#define INPUT_OUTPUT_H_

#include <complex>
#include <cstddef>
#include <fstream>
#include <istream>
#include <iterator>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "classFunction/exception.h"
#include "constants.h"
#include "internal/classFunction/iomanip.h"
#include "internal/util.h"
#include "options.h"
#include "traits.h"
#include "types.h"

namespace clara {

/**
 * @brief factory function creating display manipulator for scalar value
 *
 * this function enable custom format of scalar value
 *
 * @tparam Scalar type of the scalar value
 * @param scalar the scalar value to be display
 * @param opts optional formatting options
 * @return IOManipScalar<Scalar> object ready for streaming
 *
 * example:
 * ```
 * double x = 1e-12;
 * std::cout << disp(x);
 * std::cout << disp(x, opts);
 */
template <typename Scalar, typename std::enable_if_t<std::is_arithmetic_v<Scalar>>* = nullptr>
inline internal::IOManipScalar<Scalar> disp(Scalar scalar, IOManipScalarOpts opts = {}) {
  return internal::IOManipScalar<Scalar>{scalar, opts};
}

/**
 * @brief factory function for creating display manipulator for complex number
 *
 * support feature like
 * - supressing near-zero real or imag part using a threshold
 * - customizing delimiter
 * design for use in generic ouput system where fine control over complex number
 * formatting is needed
 *
 * @tparam T underlying floating-point type of the complex number
 * @param z complex number to be displayed
 * @param opts optional formatting option
 * @return IOManipScalar<std::complex<T>> object ready for streaming
 *
 * example:
 * ```
 * std::complex<double> z{0.0, -1.0};
 * std::cout << disp(z);
 * std::cout << disp(z, opts);
 */
template <typename T>
inline internal::IOManipScalar<std::complex<T>> disp(std::complex<T> z,
                                                     IOManipComplexOpts opts = {}) {
  return internal::IOManipScalar<std::complex<T>>{z, opts};
}

/**
 * @brief factory function for creating display manipulator for range elements
 *
 * wrapping sequence defined by two iterator
 *
 * design for use with any input iterator pair, making it suitable for
 * - standard container
 * - raw array
 * - custom memory buffer
 *
 * @tparam InputIterator type of the input iterator
 * @param first iterator pointing to the beginning of the range
 * @param last iterator pointing to the end of the range
 * @param opts optional formatting options
 * @return IOManipRange<InputIterator> object ready for streaming
 *
 * example:
 * ```
 * std::vector<double> vektor = {1.0, 0.000001, 3.0};
 * std::cout << disp(vektor.begin(), vektor.end());
 * std::cout << disp(vektor.begin(), vektor.end(), opts);
 * ```
 */
template <typename InputIterator>
internal::IOManipRange<InputIterator> disp(InputIterator first, InputIterator last,
                                           IOManipRangeOpts opts = {}) {
  return internal::IOManipRange<InputIterator>{first, last, opts};
}

/**
 * @brief factory function for creating a display manipulator for a range of elmeents
 *
 * wraps squence defined by two interator into `IOManipRange` object
 * which allow customizable formatting when printing to an output stream
 *
 * support feature
 * - custom delimiter around the entire range
 * - raw array
 * - custom memory buffer
 *
 * @tparam InputIterator type of the input iterator
 * @param first iterator pointing to the begining of the range
 * @param last iterator pointing to the end of the range
 * @param opts optional formatting options
 * @return an IOManipRange<InputIterator> object ready for streaming
 *
 * example
 * ```
 * std::vector<double> vektor = {1.0, 0.000001, 3.0};
 * std::cout << disp(vektor.begin(), vektor.end());
 * std::cout << disp(vektor.begin(), vektor.end(), opts);
 * ```
 */
template <typename Container>
internal::IOManipRange<typename Container::const_iterator> disp(
    const Container& c, IOManipContainerOpts opts = {},
    typename std::enable_if_t<is_iterable_v<Container>>* = nullptr,
    typename std::enable_if_t<!is_matrix_expression_v<Container>>* = nullptr) {
  return internal::IOManipRange<typename Container::const_iterator>{std::begin(c), std::end(c),
                                                                    opts};
}

/**
 * @brief factory afunction for creating display manipulator for pointer-barsed sequence
 *
 * wrap raw pointer its length into an `IOManipPointer` object, which allows
 * customizable formatting when printing the data to an output stream
 *
 * support feature like:
 * - custom delimiters around the entire squence
 * - custom separator between element
 * - suppression of near-zero value using threshold
 *
 * @tparam PointerType underlying value type
 * @param P pointer to the start of the data to display
 * @param N number of element to display
 * @param opts optional formatting options
 * @return an IOManipPointer<PointerType> object ready for streaming
 */
template <typename PointerType>
internal::IOManipPointer<PointerType> disp(const PointerType* p, idx N,
                                           IOManipPointerOpts opts = {}) {
  return internal::IOManipPointer<PointerType>{p, N, opts};
}

template <typename Derived>
internal::IOManipEigen disp(const Eigen::MatrixBase<Derived>& A, IOManipEigenOpts opts = {}) {
  return internal::IOManipEigen(A, opts);
}

template <typename Scalar>
internal::IOManipDirac<Scalar> disp(const clara::dirac_t<Scalar>& A, IOManipDiracOpts opts = {}) {
  return internal::IOManipDirac<Scalar>{A, opts};
}

/**
 * @brief serializes an Eigen matrix to a text-based output stream
 *
 * this function write the eigen matrix to the specified output stream custom format:
 * - first line contains the matrix dimension
 * - subsequent lines contain each row, with element separated by space
 * - complex number are written as (real, imag)
 *
 * @tparam derived the derived type of the eigen matrix expression
 * @param A eigen matrix to serialize
 * @param os output stream
 * @throws exception::ZeroSize if the matrix has zero size
 * @throws std::runtime_error if the output stream is not in good state
 */
template <typename Derived>
void save(const Eigen::MatrixBase<Derived>& A, std::ostream& os) {
  const clara::dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA)) {
    throw exception::ZeroSize("clara::save()", "A");
  }

  if (!os.good()) {
    throw std::runtime_error("clara::save(): error writting ouput stream");
  }

  // extract matrix dimension
  idx rows = static_cast<idx>(rA.rows());
  idx cols = static_cast<idx>(rA.cols());
  os << rows << " " << cols << '\n';

  bool is_cplx = clara::is_complex_v<typename Derived::Scalar>;

  // loop throught each row of the matrix
  for (idx i = 0; i < rows; ++i) {
    std::string sep;
    // loop through each column in the current row
    for (idx j = 0; j < cols; ++j) {
      if (is_cplx) {
        // if the matrix contains complex number, wrap them in parentheses
        os << sep << '(' << internal::real2text(std::real(rA(i, j)));
        os << ',' << internal::real2text(std::imag(rA(i, j))) << ')';
      } else {
        os << sep << internal::real2text(rA(i, j));
      }
      sep = " ";
    }
    os << '\n';
  }
}

/**
 * @brief loads a complex-value matrix from an input stream
 *
 * this function parses text-based representation of complex matrix
 *
 * design as counterpart to the `save(...)` function, this enable deterministic
 * serialization and deserialization of complex-valued eigen matrices
 *
 * @tparam derived eigen matrix type with complex scalar
 * @param is input stream to read
 * @param dummy SFINAE constraint to enable only for complex type
 * @return complex matrix recronstucted from the stream
 * @throws std::runtime_error if stream is invalid or data malformed
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> load(
    std::istream& is, std::enable_if_t<is_complex_v<typename Derived::Scalar>>* = nullptr) {
  if (!is.good()) {
    throw std::runtime_error("clara::load(): error opening input stream");
  }

  // read matrix dimension
  idx rows, cols;
  is >> rows >> cols;

  // create matrix of appropriate size
  clara::dyn_mat<typename Derived::Scalar> A(rows, cols);
  // temporary variables for parsing complex number
  char skip;
  decltype(std::declval<typename Derived::Scalar>().real()) re, im;

  // read each element row by row
  for (idx i = 0; i < rows; ++i) {
    for (idx j = 0; j < cols; ++j) {
      is >> skip >> re >> skip >> im >> skip;
      A(i, j) = typename Derived::Scalar{re, im};
    }
  }
  return A;
}

/**
 * @brief load real-valued matrix from an input stream
 *
 * design as counterpart to the `save(...)` function, this enables deterministic
 * serialization and deserialization of real-valued eigen matrices
 *
 * @tparam Derived the derived Eigen matrix with real scalar
 * @param is input stream to read from
 * @param dummy SFINAE constraint to eanble only fro real type
 * @return real-valued matrix recronstucted from the stream
 * @throws std::runtime_error if stream is invalid or data malformed
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> load(
    std::istream& is, std::enable_if_t<!clara::is_complex_v<typename Derived::Scalar>>* = nullptr) {
  if (!is.good()) {
    throw std::runtime_error("clara:load(): error opening input stream");
  }

  // read matrix dimension
  idx rows, cols;
  is >> rows >> cols;

  // create matrix of appropriate size
  clara::dyn_mat<typename Derived::Scalar> A(rows, cols);
  // read each element row by row
  for (clara::idx i = 0; i < rows; ++i) {
    for (idx j = 0; j < cols; ++j) {
      is >> A(i, j);
    }
  }
  return A;
}

namespace stalestuff {

/**
 * @brief serialization an eigen matrix to binary output stream
 *
 * this useful for
 * - fast serialization / deserialization
 * - efficient storage of large matrices
 * - interfacing with system requiring binary input / output
 *
 * @tparam Derived the derived type of the eigen matrix expression
 * @param A eigen matrix to serialize
 * @param os output stream to write to
 * @throws expcetion::ZeroSize if the matrix has zero size
 * @throws std::runtime_error if the output stream is not a good state
 */
template <typename Derived>
void save(const Eigen::MatrixBase<Derived>& A, std::ostream& os) {
  // get the concrete derived matrix type
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  if (!internal::check_nonzero_size(rA)) {
    throw exception::ZeroSize("clara::stalestuff::save()", "A");
  }

  if (!os.good()) {
    throw std::runtime_error("clara::stalestuff::save(): error writing output stream");
  }

  // write human-readable header to identify the file content
  const std::string header_ = "TYPE::Eigen:Matrix";
  os.write(header_.c_str(), static_cast<std::ptrdiff_t>(header_.length()));

  // extract matrix dimension
  idx rows = static_cast<idx>(rA.rows());
  idx cols = static_cast<idx>(rA.cols());

  // write the number of row and column in binary format
  os.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
  os.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
  // write raw binary data of the matrix
  // this assume the matrix is stored in column-major order
  os.write(reinterpret_cast<const char*>(rA.data()),
           sizeof(typename Derived::Scalar) * rows * cols);
}

/**
 * @brief loads an eigen matrix from binary input stream
 *
 *
 * - fixed ASCII header string: `"TYPE::Eigen::Matrix`
 * - binary-encoded number of rows and columns
 * - raw binary data of all matrix elements
 *
 * @tparam Derived the derived type of the Eigen matrix expression
 * @param is input stream to read from
 * @return A matrix recronstucted from the binary stream
 * @throws std::runtime_error if stream is invalid or file is corrupt
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> load(std::istream& is) {
  if (!is.good()) {
    throw std::runtime_error("clara::stalestuff::load(): error opening input stream");
  }

  // defined expected header string used during save
  const std::string header_ = "TYPE::Eigen::Matrix";
  // allocate a buffer to read the header from the stream
  std::unique_ptr<char[]> fheader_{new char[header_.length()]};
  if (std::string(fheader_.get(), header_.length()) != header_) {
    throw std::runtime_error("clara::stalestuff::load(): input stream is corrupt");
  }

  // read matrix dimension
  idx rows, cols;
  // read row count from binary stream
  is.read(reinterpret_cast<char*>(&rows), sizeof(rows));
  // read column count from binary stream
  is.read(reinterpret_cast<char*>(&cols), sizeof(cols));

  // create dynamic sized matrix of appropriate scalar type
  dyn_mat<typename Derived::Scalar> A(rows, cols);
  is.read(reinterpret_cast<char*>(A.data()), sizeof(typename Derived::Scalar) * rows * cols);

  return A;
}

}  // namespace stalestuff

}  // namespace clara

#endif  // !INPUT_OUTPUT_H_
