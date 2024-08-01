// Copyright (c) 2023 arfy slowy

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE

#ifndef CLARA_TRAITS_HPP_
#define CLARA_TRAITS_HPP_

#include <Eigen/Eigen>
#include <complex>
#include <type_traits>
namespace clara {
template <typename...>
struct make_void {
  typedef void type;
};

template <typename... Ts>
using to_void = typename make_void<Ts...>::type;

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 0) && !__clang__)

template <typename T, typename = void>
struct is_iterable : std::false_type {};

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
template <typename T>
struct is_iterable<T,
                   to_void<decltype(std::declval<T>().begin()), decltype(std::declval<T>().end()),
                           decltype(*(std::declval<T>().begin()))>> : std::true_type {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
template <typename Derived>
struct is_matrix_expression : std::is_base_of<Eigen::MatrixBase<typename std::decay<Derived>::type>,
                                              typename std::decay<Derived>::type> {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
template <typename T>
struct is_complex : std::false_type {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};
#if ((__GNUC__) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif  // ((__GNUC__) && (__GNUC_MINOR__ == 8) && !__clang__)

}  // namespace clara

#endif  // !CLARA_TRAITS_HPP_
