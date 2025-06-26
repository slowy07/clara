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

#ifndef CLARA_TRAITS_H_
#define CLARA_TRAITS_H_

#include <charconv>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <type_traits>

namespace clara {

// supressing -Weffc++ warning for GCC version 4.8
// diagnostic push/pop block only applies to GCC 4.8 and not Clang or ICC
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4)
        // && (__GNUC_MINOR__ == 8)

// primary template for is_iterable - which default to false_type
// this trait will be check whether a type T can be used in range-based looping
template <typename T, typename = void>
struct is_iterable : std::false_type {};
// restoring previous diagnostic settings after specialization definition
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)

// specialization of is_iterable - and become true_type if T are support
// - begin() and end() methods
// - dereferencing result of begin
// using SFINAE via std::void_t to conditionally enable the specialization
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4)
        // && (__GNUC_MINOR__ == 8)
template <typename T>
struct is_iterable<T, std::void_t<decltype(std::declval<T>().begin()),  // check for the begin
                                  decltype(std::declval<T>().end()),    // check for end
                                  decltype(*(std::declval<T>().begin()))>> : std::true_type {};

// trait to be check if type derive from Eigen::MatrixBase<>
// using to detect Eigen expression like Matrix, array, product, dll
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4)
        // && (__GNUC_MINOR__ == 8)

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic "-Weffc++"
#endif  // defined (__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ ==
        // 4) && (__GNUC_MINOR__ == 8)
template <typename Derived>
struct is_matrix_expression
    : std::is_base_of<Eigen::MatrixBase<std::decay_t<Derived>>, std::decay_t<Derived>> {};

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4)
        // && (__GNUC_MINOR__ == 8)

template <typename T>
inline constexpr bool is_matrix_expression_v = is_matrix_expression<T>::value;

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4)
        // && (__GNUC_MINOR__ == 8)
template <typename T>
struct is_complex : std::false_type {};
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4)
        // && (__GNUC_MINOR__ == 8)

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
template<typename T>
struct is_complex<std::complex<T>> : std::true_type {};
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic pop
#endif // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && (__GNUC_MINOR__ == 8)

template <typename T>
inline constexpr bool is_complex_v = is_complex<T>::value;

namespace internal {
// helper alias which to get evaluating type of an Eigen expression
// eval_t<Derived> give the concrete matrix type result from eval the expression
template <typename Derived>
using eval_t = std::decay_t<typename Eigen::MatrixBase<Derived>::EvalReturnType>;
}  // namespace internal

// checking at compile-time if matrix expression has exactly one row
// useful for enforcing dimension constraint on vectors
template <typename Derived>
bool constexpr is_bra_v() {
  return (internal::eval_t<Derived>::RowsAtCompileTime == 1);
}

// check at compile-time if matrix expression has exactly one column
template <typename Derived>
bool constexpr is_ket_v() {
  return (internal::eval_t<Derived>::ColsAtCompileTime == 1);
}

}  // namespace clara

#endif  // !CLARA_TRAITS_H_
