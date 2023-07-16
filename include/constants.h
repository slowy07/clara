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
constexpr double ee = 2.718281828459045235360287471352662497;
constexpr double inifinity = std::numeric_limits<double>::max();

inline cplx omega(idx D) {
  if (D == 0)
    throw Exception::Type::OUT_OF_RANGE;
  return exp(2.0 * pi * 1_i / static_cast<double>(D));
}

}  // namespace clara

#endif  // ! CONSTANTS_H_
