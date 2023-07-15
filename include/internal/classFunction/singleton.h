#ifndef INTERNAL_CLASSFUNCTION_SINGLETON_H_
#define INTERNAL_CLASSFUNCTION_SINGLETON_H_

#include <type_traits>
namespace clara {
namespace internal {

/**
* to implement a singleton, derive class from clara::internal::Singleton,
* make clara::internal::Singleton of class, then declare the constructor
* and desctructor of class as private. to get an instance, use the static
* member function clara::internal::Singleton::get_instance()
* (clara::internal::Singleton::get_thread_local_instance()), which return
* a reference (thread_local reference) to newly created singleton
*/
template <typename T>
class Singleton {
 protected:
  Singleton() noexcept = default;
  Singleton(const Singleton&) = delete;
  Singleton& operator=(const Singleton&) = delete;
  virtual ~Singleton() = default;

 public:
  static T& get_instance() noexcept(std::is_nothrow_constructible<T>::value) {
    static T instance;
    return instance;
  }
#ifndef NO_THREAD_LOCAL_
  static T& get_thread_local_instance() noexcept(std::is_nothrow_constructible<T>::value) {
    thread_local static T instance;
    return instance;
  }
#endif  // !NO_THREAD_LOCAL_
};

}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_CLASSFUNCTION_SINGLETON_H_
