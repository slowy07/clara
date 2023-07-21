#ifndef CLASSFUNCTION_INIT_H_
#define CLASSFUNCTION_INIT_H_

#include <iomanip>
#include <iostream>

#include "../internal/classFunction/singleton.h"

namespace clara {

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
