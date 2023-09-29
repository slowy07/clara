#ifndef CLASSFUNCTION_INIT_H_
#define CLASSFUNCTION_INIT_H_

#include <iomanip>
#include <iostream>

#include "../internal/classFunction/singleton.h"

namespace clara {

/**
 * @class Init
 * @brief singleton class for initializing the library
 *
 * the init class is a singleton class that provides initialization for the library
 * it ensures that certain configuration are set up before using other components
 * of the library
 */
class Init final : public internal::Singleton<const Init> {
  friend class internal::Singleton<const Init>;

 private:
  Init() {
    // std::cout << std::fixed;
    // std::cout << std::setprecision(4);
  }
  ~Init() {}
};
}  // namespace clara

#endif  // !CLASSFUNCTION_INIT_H_
