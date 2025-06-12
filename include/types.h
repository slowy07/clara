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

#ifndef CLARA_TYPES_H_
#define CLARA_TYPES_H_

#include <complex>
#include <cstddef>
#include <eigen3/Eigen/Dense>
#include <tuple>
#include <type_traits>

#include "traits.h"

namespace clara {

template <typename T>
inline constexpr bool is_iterable_v = std::is_convertible_v<is_iterable<T>, bool>;

#if defined(CLARA_IDX_DEFAULT)
using idx = std::size_t;
#elif defined(CLARA_IDX_SHORT)
using idx = short int;
#elif defined(CLARA_IDX_LONG)
using idx = long int;
#elif defined(CLARA_IDX_LONG_LONG)
using idx = long long int;
#elif defined(CLARA_IDX_USHORT)
using idx = unsigned short int;
#elif defined(CLARA_IDX_UINT)
using idx = unsigned int;
#elif defined(CLARA_IDX_ULONG)
using idx = unsigned long int;
#elif defined(CLARA_IDX_ULONG_LONG)
using idx = unsigned long long int;
#else
using idx = std::size_t;
#endif
static_assert(std::is_integral_v<idx>, "the type of idx must be integral");
static_assert(sizeof(idx) > 1, "the type must be at least 2 bytes long");

#if defined(CLARA_BIGINT_DEFAULT)
using bigint = long long int;
#elif defined(CLARA_BIGINT_SHORT)
using bigint = short int;
#elif defined(CLARA_BIGINT_INT)
using bigint = long int;
#elif defined(CLARA_BIGINT_LONG_LONG)
using bigint = long long int;
#else
using bigint = long long int;
#endif  // defined(CLARA_BIGINT_DEFAULT)
static_assert(std::is_integral_v<bigint>, "the type must be integral");
static_assert(std::is_integral_v<bigint>, "the type must be signed");
static_assert(sizeof(bigint) > 1, "the type must be at least 2 bytes long");

#if defined(CLARA_FP_DEFAULT)
using realT = double;
#elif defined(CLARA_FP_FLOAT)
using realT = float;
#elif defined(CLARA_FP_DOUBLE)
using realT = double;
#elif defined(CLARA_FP_LONG_DOUBLE)
using realT = long double;
#else
using realT = double;
#endif  // defined(CLARA_FP_DEFAULT)
static_assert(std::is_floating_point_v<realT>, "the type must be floating-point");

using ubigint = std::make_unsigned_t<bigint>;
using cplx = std::complex<realT>;

template <typename Scalar>
using dyn_mat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Scalar>
using dyn_col_vect = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

template <typename Scalar>
using dyn_row_vect = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

using ket = dyn_col_vect<cplx>;
using bra = dyn_row_vect<cplx>;
using cmat = dyn_mat<cplx>;
using rmat = dyn_mat<realT>;

template <typename Scalar>
struct dirac_t {
  std::vector<idx> dims_rows{};
  std::vector<idx> dims_cols{};
  std::vector<std::pair<Scalar, std::vector<idx>>> states{};

  bool operator==(const dirac_t& rhs) const {
    return std::tie(dims_rows, dims_cols, states) ==
           std::tie(rhs.dims_rows, rhs.dims_cols, rhs.states);
  }

  bool operator!=(const dirac_t& rhs) const { return !(*this == rhs); }
};

// representing quantum ram strucure, which is essentially a vector of indices
using quantum_ram = std::vector<idx>;

/**
 * @brief expression template type alias for eigen matrices
 *
 * @tparam Derived the eigen expression type being evaluated
 *
 * @detailes
 *  - scalar type is extracted from the evaluated form of the derived expression
 *  - matrix dimension are preserved at compile-time if known, otherwise to dynamic
 */
template <typename Derived>
using expr_t =
    Eigen::Matrix<typename internal::eval_t<Derived>::scalar,
                  internal::eval_t<Derived>::RowsAtCompileTime == 1 ? 1 : Eigen::Dynamic,
                  internal::eval_t<Derived>::ColsAtCompileTime == 1 ? 1 : Eigen::Dynamic,
                  internal::eval_t<Derived>::Options,
                  internal::eval_t<Derived>::RowsAtCompileTime == 1 ? 1 : Eigen::Dynamic,
                  internal::eval_t<Derived>::ColsAtCompileTime == 1 ? 1 : Eigen::Dynamic>;

/**
 * @brief utility class to create lambda overload sets
 *
 * CRTP-style helper class merges multiple callable objects
 *
 * example:
 * ```
 * auto vist = overloaded {
 *    [](int i) {
 *        std:cout << "int: " << i;
 *    },
 *    [](double d) {
 *        std::cout << "double: " << d;
 *    }
 * };
 * ```
 */
template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;  // import all operator() overload from base type
};

/**
 * @brief deduction guid for `overloaded<Ts...>`
 *
 * this allow the compile to deduce the template argument for `overloaded`
 * when it is constructed with a braced list of lambdas or functors
 */
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

}  // namespace clara

#endif  // !CLARA_TYPES_H_
