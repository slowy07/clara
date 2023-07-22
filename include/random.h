#ifndef RANDOM_H_
#define RANDOM_H_

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <random>
#include <vector>

#include "classFunction/exception.h"
#include "classFunction/random_devices.h"
#include "functions.h"
#include "input_output.h"
#include "types.h"

namespace clara {
/**
 * @brief genrate random real number uniformly distributed in the interval [a, b]
 * @return random real number (double) uniformly distributed in [a, b]
 */
inline double rand(double a, double b) {
  if (a >= b)
    throw exception::OutOfRange("clara::rand()");
  std::uniform_real_distribution<> ud(a, b);

#ifdef NO_THREAD_LOCAL_
  return ud(RandomDevices::get_instance().get_prng());
#else
  return ud(RandomDevices::get_thread_local_instance().get_prng());
#endif  // DEBUG
}

/**
 * @brief generate a random big integer uniformly distributed in the interval [a, b]
 * @return random big integer uniformly distributed in the interval
 */
inline bigint rand(bigint a, bigint b) {
  if (a > b)
    throw exception::OutOfRange("clara::rand()");
  std::uniform_int_distribution<bigint> uid(a, b);
#ifdef NO_THREAD_LOCAL_
  return uid(RandomDevices::get_instance().get_prng());
#else
  return uid(RandomDevices::get_thread_local_instance().get_prng());
#endif  // DEBUG
}

/**
 * @brief generate random index (idx) uniformly distribyted in the interval [a,b]
 * @return random index uniformly distributed in the interval[a, b]
 */
inline idx randidx(idx a = std::numeric_limits<idx>::min(),
                   idx b = std::numeric_limits<idx>::max()) {
  if (a > b)
    throw exception::OutOfRange("clara::randidx()");
  std::uniform_int_distribution<idx> uid(a, b);

#ifdef NO_THREAD_LOCAL_
  return uid(RandomDevices::get_instance().get_prng());
#else
  return uid(RandomDevices::get_thread_local_instance().get_prng());
#endif  // DEBUG
}

/**
 * @brief generate a random matrix with entries uniformly distributed
 * in the interval [a, b]
 * if complex, then both real and imaginary parts are uniformly distributed
 * in [a, b]
 */
template <typename Derived>
Derived rand(idx rows, idx cols, double a = 0, double b = 1) {
  (void)rows;
  (void)cols;
  (void)a;
  (void)b;
  throw exception::UndefinedType("clara::rand()");
}

/**
 * @brief generate a random real matrix with entries uniformly
 * distributed in the interval [a, b], the template parameter cannot
 * be automatically deduced and must be explicitly provied
 * @return random real matrix
 */
template <>
inline dmat rand(idx rows, idx cols, double a, double b) {
  if (rows == 0 || cols == 0)
    throw exception::ZeroSize("clara::rand()");
  if (a >= b)
    throw exception::OutOfRange("clara::rand()");
  return dmat::Zero(rows, cols).unaryExpr([a, b](double) { return rand(a, b); });
}

/**
 * @brief generate a random complex matrix with entries (both and real and imaginary) uniformly
 * distributed in the interval [a, b], specialized for complex matrices
 * the template parameter cannot be automatically deduced and must be explicitly
 * provided
 * @return random complex matrix
 */

template <>
inline cmat rand(idx rows, idx cols, double a, double b) {
  if (rows == 0 || cols == 0)
    throw exception::ZeroSize("clara::rand()");
  if (a >= b)
    throw exception::OutOfRange("clara::rand()");

  return rand<dmat>(rows, cols, a, b).cast<cplx>() +
         1_i * rand<dmat>(rows, cols, a, b).cast<cplx>();
}

/**
 * @brief generate a random matrix with entries normally distributed in N(mean, sigma)
 * if complex, then both real and imaginary parts are normally distributed in N(mean, sigma)
 */
template <typename Derived>
Derived randn(idx rows, idx cols, double mean = 0, double sigma = 1) {
  (void)rows;
  (void)cols;
  (void)mean;
  (void)sigma;
  throw exception::UndefinedType("clara::randn()");
}

/**
 * @brief generate a random real matrix with entries normally
 * distributed in N(mean, sigma) specialization for double matrices
 * this template parameter cannot be automatically deduced and
 * must be explicitly provied
 */
template <>
inline dmat randn(idx rows, idx cols, double mean, double sigma) {
  if (rows == 0 || cols == 0)
    throw exception::ZeroSize("clara::randn()");
  std::normal_distribution<> nd(mean, sigma);
  return dmat::Zero(rows, cols).unaryExpr([&nd](double) {
#ifdef NO_THREAD_LOCAL_
    return nd(RandomDevices::get_instance().get_prng());
#else
    return nd(RandomDevices::get_thread_local_instance().get_prng());
#endif  // DEBUG
  });
}

/**
 * @brief generate random complex matrix with entries (both real and imaginary)
 * normally distributer in N(mean, sigma)
 * @return random complex matrix
 */
template <>
inline cmat randn(idx rows, idx cols, double mean, double sigma) {
  if (rows == 0 || cols == 0)
    throw exception::ZeroSize("clara::randn()");
  return randn<dmat>(rows, cols, mean, sigma).cast<cplx>() +
         1_i * randn<dmat>(rows, cols, mean, sigma).cast<cplx>();
}

/**
 * @brief generate random real number (double) normally distributed in N(mean, sigma)
 * @return random real number normally distributed in N(mean, sigma)
 */
inline double randn(double mean = 0, double sigma = 1) {
  std::normal_distribution<> nd(mean, sigma);
#ifdef NO_THREAD_LOCAL_
  return nd(RandomDevices::get_instance().get_prng());
#else
  return nd(RandomDevices::get_thread_local_instance().get_prng());
#endif  // DEBUG
}

inline cmat randU(idx D) {
  if (D == 0)
    throw exception::DimsInvalid("clara::randU()");

  cmat X = 1 / std::sqrt(2.) * randn<cmat>(D, D);
  Eigen::HouseholderQR<cmat> qr(X);

  cmat Q = qr.householderQ();
  Eigen::VectorXcd phases = (rand<dmat>(D, 1)).cast<cplx>();
  for (idx i = 0; i < static_cast<idx>(phases.rows()); ++i)
    phases(i) = std::exp(2 * pi * 1_i * phases(i));
  Q = Q * phases.asDiagonal();
  return Q;
}

/**
 * @brief generate random isometry matrix
 * @return random isometry matrix
 */
inline cmat randV(idx Din, idx Dout) {
  if (Din == 0 || Dout == 0 || Din > Dout)
    throw exception::DimsInvalid("clara::randV()");
  return randU(Dout).block(0, 0, Dout, Din);
}

/**
 * @brief generate set random of kraus operator
 * @return set N kraus operators satisfying the closure condition
 */
inline std::vector<cmat> randkraus(idx N, idx D) {
  if (N == 0)
    throw exception::OutOfRange("clara::randkraus()");
  if (D == 0)
    throw exception::DimsInvalid("clara::randkraus()");

  std::vector<cmat> result(N);
  for (idx i = 0; i < N; ++i)
    result[i] = cmat::Zero(D, D);
  cmat Fk(D, D);
  cmat U = randU(N * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(3)
#endif  // DEBUG
  for (idx k = 0; k < N; ++k)
    for (idx a = 0; a < D; ++a)
      for (idx b = 0; b < D; ++b)
        result[k](a, b) = U(a * N + k, b * N);
  return result;
}

/**
 * @brief generate a random hermitian matrix
 * @return random hermitian matrix
 */
inline cmat randH(idx D) {
  if (D == 0)
    throw exception::DimsInvalid("clara::randH()");
  cmat H = 2 * rand<cmat>(D, D) - (1. + 1_i) * cmat::Ones(D, D);
  return H + adjoint(H);
}

/**
 * @brief generate random normalized ket (pure state vector)
 * @return random normalized ket
 */
inline ket randket(idx D) {
  if (D == 0)
    throw exception::DimsInvalid("clara::randket()");
  ket kt = randn<cmat>(D, 1);
  return kt / norm(kt);
}

/**
 * @brief generate a random density matrix
 * @return random density matrix
 */
inline cmat randrho(idx D) {
  if (D == 0)
    throw exception::DimsInvalid("clara::randrho()");
  cmat result = 10 * randH(D);
  result = result * adjoint(result);
  return result / trace(result);
}

/**
 * @brief generate a random uniformly distributed permutation
 * uses knuth shuffle method (as implemented by std::shuffle)
 * so that all permutation are equally probable
 * @return random permutation of size N
 */
inline std::vector<idx> randperm(idx N) {
  if (N == 0)
    throw exception::PermInvalid("clara::randperm()");
  std::vector<idx> result(N);

  // fill increasing oreder
  std::iota(std::begin(result), std::end(result), 0);
#ifdef NO_THREAD_LOCAL_
  std::shuffle(std::begin(result), std::end(result), RandomDevices::get_instance().get_prng());
#else
  std::shuffle(std::begin(result), std::end(result),
               RandomDevices::get_thread_local_instance().get_prng());
#endif  // NO_THREAD_LOCAL_
  return result;
}

/**
 * @brief generate random probability vector uniformly distributed over the
 * probability simplex
 * @return random probability vector
 */
inline std::vector<double> randprob(idx N) {
  if (N == 0)
    throw exception::OutOfRange("clara::randprob()");
  std::vector<double> result(N);

  // generate
  std::exponential_distribution<> ed(1);
  for (idx i = 0; i < N; ++i) {
#ifdef NO_THREAD_LOCAL_
    result[i] = ed(clara::RandomDevices::get_instance().get_prng());
#else
    result[i] = ed(clara::RandomDevices::get_thread_local_instance().get_prng());
#endif  // NO_THREAD_LOCAL_
  }
  // normalize
  double sumprob = sum(result);
  for (idx i = 0; i < N; ++i)
    result[i] /= sumprob;
  return result;
}

}  // namespace clara

#endif  // !RANDOM_H_
