#include <cstdlib>
#include <ctime>
#include "../include/clara.h"
#include "../include/gates.h"
#include "../include/stat.h"

namespace clara {

// make the random_device generator visible
std::random_device stat::_rd;
// make the mt19937 generator visible
std::mt19937 stat::_rng(_rd());

int _init() {
  // initialize the gates
  gt::_init_gates();
  // seed the standard random number generator
  std::srand(static_cast<unsigned int>(stat::_rd()));
  return 0;
}
}
