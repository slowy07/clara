#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <limits>

#include "classFunction/exception.h"
#include "types.h"

namespace clara {
inline constexpr cplx operator"" _i(unsigned long long int x) noexcept {
  return {0., static_cast<double>(x)};
}

inline constexpr cplx operator"" _i(long double x) noexcept { return {0., static_cast<double>(x)}; }

/**
 * @brief double chop used in clara::disp() for setting to zero numbers
 * that have their absolute value smaller than clara:chop
 */
constexpr double chop = 1e-10;
/**
 * @brief used to decide whether a number or epression in double precision
 * os zero or not
 */
constexpr double eps = 1e-12;
/**
 * @brief maxium number of allowed qubit
 * used internally to allocate array on the stack
 */
constexpr idx maxn = 64;
constexpr double pi = 3.141592653589793238462643383279502884;

/**
 * @brief representing the mathematical constant e
 */
constexpr double ee = 2.718281828459045235360287471352662497;

/**
 * @brief the constant infinity representing the maximum value of a double
 *        this is used to represent an unbound value in certain calculations
 */
constexpr double inifinity = std::numeric_limits<double>::max();

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

}  // namespace clara

#endif  // ! CONSTANTS_H_
