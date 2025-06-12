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

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <limits>

#include "classFunction/exception.h"
#include "types.h"

namespace clara {

inline namespace literals {
inline constexpr cplx operator"" _i(unsigned long long int x) noexcept {
  return {0., static_cast<double>(x)};
}

inline constexpr cplx operator"" _i(long double x) noexcept { return {0., static_cast<double>(x)}; }

inline constexpr std::complex<float> operator"" _if(unsigned long long int x) noexcept {
  return {0., static_cast<float>(x)};
}

inline constexpr std::complex<float> operator"" _if(long double x) noexcept {
  return {0., static_cast<float>(x)};
}

}  // namespace literals

/**
 * @brief used to decide whether a number or epression in double precision
 * os zero or not
 */
constexpr realT eps = 1e-12;
/**
 * @brief maxium number of allowed qubit
 * used internally to allocate array on the stack
 */
constexpr idx maxn = 64;
constexpr realT pi = 3.141592653589793238462643383279502884;

/**
 * @brief representing the mathematical constant e
 */
constexpr realT ee = 2.718281828459045235360287471352662497;

/**
 * @brief the constant infinity representing the maximum value of a double
 *        this is used to represent an unbound value in certain calculations
 */
constexpr double inifinity = std::numeric_limits<realT>::infinity();

/**
 * @brief function to calculate the complex number omega, used in certain
 *        quamtum computation
 * @param D input value 'D' used to calculate the omega value
 * @return std::complex<double> the complex number omga calculated as exp(2.0 * pi * 1_1 /
 * static_cast<double>D)
 * @exception exception::OutOfRange thrown if 'D' is zero
 */
inline cplx omega(idx D) {
  if (D == 0)
    throw exception::OutOfRange("clara::omega()");
  return exp(2.0 * pi * 1_i / static_cast<double>(D));
}

enum {
  RES = 0,   // measurement
  PROB = 1,  // probability
  ST = 2,    // output state
};

}  // namespace clara

#endif  // ! CONSTANTS_H_
