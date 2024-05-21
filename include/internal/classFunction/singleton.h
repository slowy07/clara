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

#ifndef INTERNAL_CLASSFUNCTION_SINGLETON_H_
#define INTERNAL_CLASSFUNCTION_SINGLETON_H_

#include <type_traits>
namespace clara {
namespace internal {

/**
 * @class Singleton
 * @brief Implementation of a singleton pattern
 * @tparam T the type of the class that needs to be a singleton. the class should have
 *         a private constructor and descructor
 */
template <typename T>
class Singleton {
 protected:
  Singleton() noexcept = default;
  Singleton(const Singleton&) = delete;
  Singleton& operator=(const Singleton&) = delete;
  virtual ~Singleton() = default;

 public:
  /**
   * @brief get the instance of the singleton class
   * @return reference to the singleton instance
   * NOTE: the first call to this function will construct the singleton instance
   *       subsequemt calls will return the same instance.
   */
  static T& get_instance() noexcept(std::is_nothrow_constructible<T>::value) {
    static T instance;
    return instance;
  }
#ifndef NO_THREAD_LOCAL_
  /**
   * @brief get the thread-local instance of the singleton class
   * @return reference to the thread-local singleton instance
   * NOTE: the first class to this function in each thread will construct a separate
   *       instance for that thread. subsequent  calls in the same thread will return
   *       the same thread-local instance
   */
  static T& get_thread_local_instance() noexcept(std::is_nothrow_constructible<T>::value) {
    thread_local static T instance;
    return instance;
  }
#endif  // !NO_THREAD_LOCAL_
};

}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_CLASSFUNCTION_SINGLETON_H_
