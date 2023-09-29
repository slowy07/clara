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
