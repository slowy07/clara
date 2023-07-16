#ifndef TRAITS_H_
#define TRAITS_H_

#include <complex>
#include <eigen3/Eigen/Dense>
#include <type_traits>
namespace clara {

template <typename... Ts>
struct make_void {
  typedef void type;
};
template <typename... Ts>
using to_void = typename make_void<Ts...>::type;

/**
 * @brief check whether T is compatible with STL-like iterable container
 * provide the constant member \a value which is equal to \a true,
 * if \a T compabitle with an iterable container.
 */
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
template <typename T, typename = void>
struct is_iterable : std::false_type {};

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

/**
 * @brief check whether T is compatible with an STL-like iterable container
 * specialization for STL-like iterable containers
 */
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
template <typename T>
struct is_iterable<T, to_void<decltype(std::declval<T>().begin()),
                              decltype(std::declval<T>().end()), typename T::value_type>>
    : std::true_type {};

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

/**
 * @brief check whether the type is an Eigen matrix expression
 * provides the constant member value which is equal to true,
 * if the type is an Eigen matrix expression of type Eigen::MatrixBase<Derived>,
 */
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
template <typename Derived>
struct is_matrix_expression : std::is_base_of<Eigen::MatrixBase<typename std::decay<Derived>::type>,
                                              typename std::decay<Derived>::type> {};

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

/**
* @brief check whether the type is complex type
* provide the constant member value which is equal to true,
* if the type is a complex type.
*/
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
template <typename T>
struct is_complex : std::false_type {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

/**
* @brief check whether the type is complex number type,
* specialization for complex type
*/
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragam GCC diagnostic ignored "-Weffc++"
#endif // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
 #endif // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)

}  // namespace clara

#endif  // !TRAITS_H_
