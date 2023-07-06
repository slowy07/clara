#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "types.h"

namespace clara {
namespace ct {

// for setting to zero everything that is smaller in absolute magniteude
const double chop = 1e-10;

// imaginary i (square root of -1)
const types::cplx ii(0, 1);
// pi notation
const double pi = 3.141592653589793238462643383279502884;
// base natural log
const double ee = 2.718281828459045235360287471352662497;

inline types::cplx omega(size_t D) { return exp(2.0 * pi * ii / static_cast<double>(D)); }

}  // namespace ct
}  // namespace clara

#endif  // !CONSTANTS_H_
