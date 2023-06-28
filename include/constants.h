#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "types.h"

namespace clara {
namespace ct {
const types::cplx ii(0, 1);
const double pi = 3.141592653589793238462643383279502884;
const double ee = 2.718281828459045235360287471352662497;

inline types::cplx omega(size_t D) {
  return exp(2.0 * pi * ii / static_cast<double>(D));
}

}
}

#endif // !CONSTANTS_H_

