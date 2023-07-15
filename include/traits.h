#ifndef TRAITS_H_
#define TRAITS_H_

#include <complex>
#include <eigen3/Eigen/Dense>
#include <type_traits>
namespace clara {

template <typename...>
using to_void = void;

template <typename T, typename = void>
struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, to_void<decltype(std::declval<T>().begin()),
                              decltype(std::declval<T>().end()), typename T::value_type>>
    : std::true_type {};

template <typename Derived>
struct is_matrix_expression : std::false_type {};

template <typename Derived>
struct is_matrix_expression<typename Eigen::MatrixBase<Derived>> : std::true_type {};

template<typename T>
struct is_complex: std::false_type {};

template<typename T>
struct is_complex<std::complex<T>> : std::true_type {};

}  // namespace clara

#endif  // !TRAITS_H_
