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
