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
 * @brief template singleton class to enforce single instance of type T
 *
 * this class provides a thread-safe way to manage a singleton object
 * it support both global (shared across threads) and thread-local instances,
 * depending on whether NO_THREAD_LOCAL is defined
 */
template <typename T>
class Singleton {
 public:
  /**
   * @brief delete copy constructor and assignment operator
   *
   * prevent copying of the singleton instance to make sure only one exists
   */
  Singleton(const Singleton&) = delete;
  Singleton& operator=(const Singleton&) = delete;

  /**
   * get global (non-thread-local) singleton instance
   *
   * usign static local variable which is initialized once and shared globally
   * thread safe in C++11 and later due guaranteed initialization order
   */
  static T& get_no_thread_local_instance() noexcept(std::is_nothrow_constructible_v<T>) {
    static T instance{};
    return instance;
  }

#ifndef NO_THREAD_LOCAL
  /**
   * @brief get thread-local singleton instance
   *
   * each thread gets its own unique instance of T
   * avoids synchronization overhead and race condition between threads
   */
  static T& get_thread_local_instance() noexcept(std::is_nothrow_constructible_v<T>) {
    thread_local static T instance{};  // one instance per thread
    return instance;
  }
#endif  // !NO_THREAD_LOCAL

  /**
   * @brief get the appropriate singleton instance based on configuration
   *
   * if NO_THREAD_LOCAL is defined, returns global instance
   * otherwise, return thread-local instance
   */
  static T& get_instance() noexcept(std::is_nothrow_constructible_v<T>) {
#ifdef NO_THREAD_LOCAL
    return get_no_thread_local_instance();  // using non-thread-local version
#else
    return get_thread_local_instance();  // used thread-local version
#endif  // NO_THREAD_LOCAL
  }

 protected:
  /**
   * @brief protected default constructor
   *
   * prevents external instatiation while allowing derived classes to construct
   * marked noexcept to support environments where exceptions are disabled
   */
  Singleton() noexcept = default;

  /**
   * @brief virtual destructor for safe inheritance
   *
   * make sure correct destruction wheren deleting through a base class pointer
   * pure virtual destructor make this class abstract
   */
  virtual ~Singleton() = default;
};

}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_CLASSFUNCTION_SINGLETON_H_
