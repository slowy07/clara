#ifndef INTERNAL_CLASSFUNCTION_VECTOR_PROXY_H_
#define INTERNAL_CLASSFUNCTION_VECTOR_PROXY_H_

#include <cstddef>
#include <numeric>
#include <type_traits>
#include <vector>

#include "../../classFunction/exception.h"

namespace clara {
namespace internal {

/**
 * @class VectorProxy
 * @brief provide proxy interface to access elements of a vector through an indirection layer
 *
 * this class allows accessing elements of a base vector using a label array
 *
 * @tparam T type element stored underlying vector
 * @tparam is_const boolean flag indicating whether access should be read-only
 */
template <class T, bool is_const>
class VectorProxy {
  // define VecType based on whether the proxy is const or not
  using VecType = std::conditional_t<is_const, const std::vector<T>, std::vector<T>>;
  // reference to the underlying data vector
  VecType& data_{};
  // label array that maps proxy indices to actual vector indices
  std::vector<std::size_t> label_{};

 public:
  /**
   * @brief construct new VectorProxy with given data and explicit label mapping
   *
   * this constructor binds the proxy to a specific data vector and provides a custom index mapping
   *
   * @param data reference to the original vector to be accessed via proxy
   * @param label vector of indices into `data`, defining the logical order
   */
  VectorProxy(std::vector<T>& data, std::vector<std::size_t> label)
      : data_{data}, label_{std::move(label)} {}

  /**
   * @brief construct new VectorProxy with default labeling
   *
   * if no label provided, it default to squence from 0 to size - 1
   * effectively providing direct access to all elements in the original order
   *
   * @param data reference to the original vector to be accessed via proxy
   */
  VectorProxy(std::vector<T>& data) : data_{data} {
    label_.resize(data_.size());
    std::iota(label_.begin(), label_.end(), 0);
  }

  /**
   * @brief const version of operator[] access element at index `i`
   *
   * bounds check access to the proxied element, throw if the index is out of range
   * or refers to an invalid position in the underlying data
   *
   * @param i logical index used to acess the proxied data
   * @return const reference to the element at mapped index
   */
  const T& operator[](std::size_t i) const {
    if (i + 1 > label_.size()) {
      throw exception::OutOfRange("clara::internal::VectorProxy::operator[]() const", "i");
    }

    if (label_[i] + 1 > data_.size()) {
      throw exception::OutOfRange("clara::internal::VectorProxy::operator[]() const", "i");
    }
    return data_[label_[i]];
  }

  /**
   * @brief non-const version of operator[] access element at index `i`
   *
   * only availabel when `is_const == false`, provide mutable access to the element at the mapped
   * index after performing bounds checks
   *
   * @tparam B internal template parameter to enable SFINAE condition
   * @tparam std::enable_if_t<!B> make sure overload only exist for non-cost proxied
   * @param i logical index used to acess the proxied data
   * @return reference to the element at mapped index
   */
  template <bool B = is_const, typename = std::enable_if_t<!B>>
  T& operator[](std::size_t i) {
    if (i + 1 > label_.size()) {
      throw exception::OutOfRange("clara::internal::VectorProxy::operator[]()", "i");
    }

    if (label_[i] + 1 > data_.size()) {
      throw exception::OutOfRange("clara::internal::VectorProxy::operator[]()", "i");
    }
    return data_[label_[i]];
  }
};
}  // namespace internal
}  // namespace clara

#endif  // !ENOLA_CLASSFUNCTION_VECTOR_PROXY_H_
