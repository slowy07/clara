#ifndef CLARA_TYPES_H_
#define CLARA_TYPES_H_

#include <complex>
#include <cstddef>
#include <eigen3/Eigen/Dense>
#include <tuple>
#include <type_traits>
#include "traits.h"

namespace clara {

template<typename T>
inline constexpr bool is_iterable_v = std::is_convertible_v<is_iterable<T>, bool>;

#if defined(CLARA_IDX_DEFAULT)
using idx = std::size_t;
#elif defined(CLARA_IDX_SHORT)
using idx = short int;
#elif defined(CLARA_IDX_LONG)
using idx = long int;
#elif defined(CLARA_IDX_LONG_LONG)
using idx = long long int;
#elif defined(CLARA_IDX_USHORT)
using idx = unsigned short int;
#elif defined(CLARA_IDX_UINT)
using idx = unsigned int;
#elif defined(CLARA_IDX_ULONG)
using idx = unsigned long int;
#elif defined(CLARA_IDX_ULONG_LONG)
using idx = unsigned long long int;
#else
using idx = std::size_t;

static_assert(std::is_integral_v<idx>, "the type of idx must be integral");
static_assert(sizeof(idx) > 1, "the type must be at least 2 bytes long");

#if defined(CLARA_BIGINT_DEFAULT)
using bigint = long long int;
#elif defined(CLARA_BIGINT_SHORT)
using bigint = short int;
#elif defined(CLARA_BIGINT_INT)
using bigint = long int;
#elif defined(CLARA_BIGINT_LONG_LONG)
using bigint = long long int;
#else
using bigint = long long int;
#endif  // defined(CLARA_BIGINT_DEFAULT)
static_assert(std::is_integral_v<bigint>, "the type must be integral");
static_assert(std::is_integral_v<bigint>, "the type must be signed");
static_assert(sizeof(bigint) > 1, "the type must be at least 2 bytes long");

#if defined(CLARA_FP_DEFAULT)
using realT = double;
#elif defined(CLARA_FP_FLOAT)
using realT = float;
#elif defined(CLARA_FP_DOUBLE)
using realT = double;
#elif defined(CLARA_FP_LONG_DOUBLE)
using realT = long double;
#else
using realT = double;
#endif  // defined(CLARA_FP_DEFAULT)
static_assert(std::is_floating_point_v<realT>, "the type must be floating-point");

using ubigint = std::make_unsigned_t<bigint>;
using cplx = std::complex<realT>;

template <typename Scalar>
using dyn_mat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Scalar>
using dyn_col_vect = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

template <typename Scalar>
using dyn_row_vect = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

using ket = dyn_col_vect<cplx>;
using bra = dyn_row_vect<cplx>;
using cmat = dyn_mat<cplx>;
using rmat = dyn_mat<realT>;

#endif  // defined(CLARA_IDX_DEFAULT)

template <typename Scalar>
struct dirac_t {
  std::vector<idx> dims_rows{};
  std::vector<idx> dims_cols{};
  std::vector<std::pair<Scalar, std::vector<idx>>> states{};

  bool operator==(const dirac_t& rhs) const {
    return std::tie(dims_rows, dims_cols, states) ==
           std::tie(rhs.dims_rows, rhs.dims_cols, rhs.states);
  }

  bool operator!=(const dirac_t& rhs) const { return !(*this == rhs); }
};

}  // namespace clara

#endif  // !CLARA_TYPES_H_
