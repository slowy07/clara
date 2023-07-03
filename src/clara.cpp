#include <cstdlib>
#include <ctime>
#include "../include/clara.h"

namespace clara {
// random number initialization

// make random device generator
std::random_device stat::_rd;
// make the mt19937 generator
std::mt19937 stat::_rng(_rd());

// initialization gate
namespace gt {
// various matrices
Eigen::MatrixXcd H, Id2, X, Y, Z, S, T;
Eigen::MatrixXcd CNOT, CP;
Eigen::MatrixXcd TOF(8, 8);
}

// initialize the library
int _init() {
  // intialize gate
  gt::_init_gates();
  // seed the standard random number generator
  std::srand(static_cast<unsigned int>(stat::_rd()));
  return 0;
}
}
