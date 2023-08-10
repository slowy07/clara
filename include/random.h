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
 * @brief generate a random real number uniformly distributed in the internal [a, b]
 *
 * this function generate a random real number uniformly distributed in the interval a and b
 * it uses std::uniform_real_distribution from cpp standard library to achieve uniform distribution
 * the function takes two parameters a and b, representing the lower and upper bound of the
 * interval it throws an exception if 'a' is greater than equal to 'b'
 *
 * @param a the lower bound of the interval for the random real number
 * @param b the upper bound of the interval for the random real number
 * @return return a random number uniformly distributed in a b
 *
 * @note the function check if 'a' is greater than or equal to 'b' and throws and exception if so.
 *       it uses std::uniform_real_distribution to generate a random real number in the specified
 * interval if the macro NO_THREAD_LOCAL_ is defined, it uses the global
 * RandomDevices::get_instance().get_prng() to obtain the random number generator. then returns the
 * generated random real number
 *
 * @example
 * // usage of rand functio to generate a random real number in the interval [0.0, 1.0]
 * double lower_bound = 0.0;
 * double upper_bound = 1.0;
 * double result = rand(lower_bound, upper_bound)
 * */
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
 * @brief generate a random big integer uniformly sistributed in the interval [a, b]
 * this function generates a random big integer uniformly distributed in the interval [a, b]
 * it uses std::uniform_int_distribution from standard library to achieve uniform distribution
 * the function takes two parameters 'a' and 'b', representing the lower and upper bounds of the
 * interval it throws an exception if 'a' is gerater than 'b'
 *
 * @param a the lower bound of the interval for the random big integer
 * @param b the upper bound of the interval for the random big integer
 *
 * @note the function check if 'a' is greater than 'b' and throws an exception if so
 *       it uses std::uniform_int_distribution to generate random big integer in the specified
 * inteval. if macro NO_THREAD_LOCAL_ is defined, it uses the global
 * RandomDevices::get_instance().get_prng() to obtain the random number generator, otherwise, it
 * uses the thread local RandomDevices::get_thread_local_instance().ger_prng() to obtain the random
 * number generatro the function then returns the gnerated random big integer
 *
 * @example
 * // usage of rand function to generate a random big integer in the interval [0, 100]
 * bigint lower_bound = 0;
 * bigint upper_bound = 100;
 * bigint result = rand(lower_bound, upper_bound);
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
 * @brief generate a random index (idx) uniformly in distribted in the interval [a, b]
 * this function generate a random index (idx) uniformly distributed in the interval [a, b]
 * it uses std::uniform_int_distribution from standard library to achieve uniform distribution
 * the function takes two parameters a and b. representing the lower and upper bound of the
 * interval it throws exception if 'a' is greater than b
 *
 * @param a the lower bound of the interval for the random index
 * @param b the upper bound of the interval for the random index
 * @return random uniformly distributed in the interval [a, b]
 *
 * @note the function check if 'a' is greater than 'b' and throws an exception if so
 *       it uses std::uniform_int_distribution to generate a random index in the specified interval
 *       if the macro macro NO_THREAD_LOCAL_ is defined, it uses the global
 * RandomDevices::get_instance().get_prng() to obtain the random number generator, otherwise it uses
 * the thread-local
 *
 * @eample
 * // usage of randidx function to generate a random index in the interval [0, 10]
 * idx lower_bound = 0;
 * idx upper_bound = 10;
 * idx result = randidx(lower_bound, upper_bound);
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
 * @brief generate random matrix with entries uniformly distributed in the inteval [a, b]
 * this function generate a random matrix with entries uniformly dostributed in the interval [a, b]
 * the matrixx can be of any type that is derived from the eigen library
 * the function takes four parameters rows, cols, a , b
 * 'rows' and 'cols' represent the number of rows and columns of the generated matrix, respectively
 * 'a' and 'b' represent the lower and upper bounds of the interval for the random entries
 * if the matrix is complex, both the real and imaginary parts of the entries are uniformly
 * distributed in the inteval [a, b]
 *
 * @tparam Derived the type of the matrix, which must be derived from the eigen library
 * @param rows the number of rows of the generated matrix
 * @param cols the number of columns of the generated matrix
 * @param a the lower bound of the interval for the random entries
 * @parma b the upper bound of the interval for the random matrix
 *
 * @throws exception::UndefinedType if the function is called with a type that is not derived from
 *                                  the eigen library
 *
 * @example
 * // using of the rand function to generate a 3x3 random matrix with entries in the interval [-1,
 * 1] int rows = 3; int cols = 3; Eigen::MatrixXd random_matrix = rand<Eigen::MatrixXd>(rows, cols,
 * lower_bound, upper_bound);
 *
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
 * @brief generate random real matrix with entries uniformly distributed in the interval [a, b]
 * this function generate a random real matrix with entries uniformly distributed in the interval
 * The template parameter 'dmat' must be explicitly provided since it cannot be automatically
 * deduced the function takes four parameters row, cols, a and b 'rows' and 'cols' represent the
 * number of rows and columns of the generated matrix, respectively 'a' and 'b' represent the lower
 * and upper bounds of the interval for the random entries
 *
 * @param rows the number of rows of the generated matrix
 * @param cols the number of columns the generated matrix
 * @param a lower bound of the inteval for the random entries
 * @param b upper bound of the interval for the random entries
 * @return a random real matrix with entries uniformly distributed in the interval [a, b]
 *
 * @throws exception::ZeroSize if either 'rows' or 'cols' is zero
 * @throws exception::OutOfRange if 'a' is greater than or equal to 'b'
 *
 * @example
 * // usage of the rand function to generate 3x3 random real matrix with entries interval [-1, 1]
 * int rows = 3;
 * int cols = 3;
 * double lower_bound = -1;
 * double upper_bound = 1;
 * Eigen::MatrixXd random_matrix = rand<Eigen::MatrixXd>(rows, cols, lower_bound, upper_bound);
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
 * @brief generate a random real matrix with entries uniformly distributed in the interval [a, b]
 *
 * this function is specialized for complex matrices (cmat) and uses the 'rand' function to generate
 * random real matrices (dmat) with entries uniformly distributed in the interval [a, b]. The
 * template parameter 'cmat' must be explicitly provided since it cannot be automatically deduced.
 * the function takes four parameters: 'rows', 'cols', 'a', and 'b'.
 * 'rows' and 'cols' represent the number of rows and columns of the generated complex matrix,
 * respectively. 'a' and 'b' represent the lower and upper bounds of the interval for the random
 * real and imaginary entries.
 *
 * @param rows the number of rows of the generated matrix
 * @param cols the number of columns of the generated matrix
 * @param a lower bound of the interval for the random entries
 * @param b upper bpund of the interval for the random entries
 * @return random real matrix with entries uniformly distributed in the interval
 *
 * @throws exception::ZeroSize if either 'rows' or 'cols' is zero
 * @throws exception::OutOfRange if 'a' is greater than or equal to 'b'
 *
 * @example
 * // usage of the rand function to generated a 3x3 random real matrix with entries interval [-1, 1]
 * int rows = 3;
 * int cols = 3;
 * double lower_bound = 1;
 * double upper_bound = 1;
 * Eigen::MatrixXd random_complex_matrix = rand<Eigen::MatrixXd>(rows, cols, lower_bound,
 * upper_bound);
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
 * @brief generate a random matrix with entries normally distribution in N(mean, sigma);
 *
 * this function is template function that generate random matrices with entries normally
 * distributed in N(mean, sigma). the template parameter 'Derived' represents the type of the
 * generated matrix and must be explicitly provided as Eigen's matrix type the function takes four
 * parameters rows, cols, mean, and sigma. rows and cols repersent the number of rows and columns of
 * the generated matrix, respectively. 'mean' repersent the mean (average) value of the normal
 * distribution, and 'sigma' represent the standard deviation of the normal distribution
 *
 * @tparam Derived the type of the generated matrix, which must be explicitly provided as Eigen's
 * matrix type
 * @param rows the number of rows of generated matrix
 * @param cols the number of columns of the gnerated matrix
 * @param mean the mean (average) value of the normal distribution
 * @param sigma the standard deviation of the normal distribution
 * @return random matrix with entries normally distributed in N(mean, sigma)
 *
 * @throw exception::UndefinedType if the template parameter 'Derived' is not provided or is not a
 *                                 valid Eigen Matrix type
 *
 * @example
 * // usage of the randn function to generate a 3x3 random real matrix with entries normally
 * // distributed
 * int rows = 3;
 * int cols = 3;
 * double mean = 0;
 * double sigma = 1;
 * Eigen::MatrixXd random_real_matrix = randn<Eigen::MatrixXd>(rows, cols, mean, sigma);
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
 * @brief generate a random real matrix with entries normally distributed in N(mean, sigma)
 *
 * this is specialization of the 'randn' function for generating random real matrices with
 * entries normally distributed in N(mean, sigma)
 * the function takes four parameters: 'rows', 'cols', 'mean', and 'sigma',
 * 'rows' and 'cols' represent the number of rows and columns of the generated matrix, respectively
 * 'mean' represents the mean (average) value of the normal distribution, and 'sigma' represents the
 * standard deviation of the normal distribution
 *
 * @throws exception::ZeroSize if either 'rows' or 'cols' is zero, indicating an invalid matrix size
 *
 * @example
 * // usage of the randn function to generate a 3x3 random real matrix with entries normally
 * // distributed with mean 0 and standard deviation 1
 * int rows = 3;
 * int cols = 3;
 * double mean = 0;
 * double sigma = 1;
 * Eigen::MatrixXd random_real_matrix = randn(rows, cols, mean, sigma);
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
 * @brief generate a random complex matrix with entries (both real and imaginary parts)
 *        anormally distributed in N(mean, sigma)
 *
 * this is specialization of the 'randn' function for generating random comple matrices with
 * entries normally distributed in N(mean, sigma) the function takes four parameters:
 * rows, cols, mean, sigma. rows and cols repersent tje number of rows and columns of the
 * generated matrix respectively. 'mean' repersent the mean (average) value of the normal
 * distribution and 'sigma' represent the standard deviation of the normal distribution
 *
 * @param rows the number of the generated matrix
 * @param cols the number of the columns of the generated matrix
 * @param mean the mean (average) value of the normal distribution
 * @param sigma the standard deviation of the normal distribution
 * @return random complex matrix with entries normally distributed in N(mean, sigma)
 *
 * @throws exception::ZeroSize if either 'rows' or 'cols' is zero, indicating an invalid matrix size
 *
 * @example
 * // usage of the randn function to generate 3x3 random complex matrix with entries normally
 * distributed
 * // with mean 0 and standard deviation 1
 * int rows = 3;
 * int cols = 3;
 * double mean = 0;
 * double sigma = 1;
 * Eigen::MatrixXcd random_complex_matrix = randn(rows, cols, mean, sigma);
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
 * @brief generate a random unitary matrix of size DxD using Haar measure
 *
 * this function generates a random unitary matrix of size DxD using the Haar measure
 * which is natural measure for the group of unitary matrices. the generated matrix
 * is uniformly distributed over the unitary group
 *
 * @param D the size of the square matrix
 * @return random unitary matrix of size DxD
 *
 * @throws exception::DimsInvalid if 'D' is zero, indicating an invalid matrix ZeroSize
 *
 * @example
 * int D = 2;
 * Eigen::MatrixXcd random_unitary_matrix = randU(D);
 */
inline cmat randV(idx Din, idx Dout) {
  if (Din == 0 || Dout == 0 || Din > Dout)
    throw exception::DimsInvalid("clara::randV()");
  return randU(Dout).block(0, 0, Dout, Din);
}

/**
 * @brief generate a set of random kraus operators
 *
 * this function generates a set of N kraus operators satfisying the closure condition for a quantum
 * operation on a D-dimensional quantum system, the generated kraus operators are represented as a
 * vector of complex matrices
 *
 * @param N the number of krause operators
 * @param D the size of the square kraus operators
 * @return vector containing N random kraus operators, each of size DxD
 *
 * @throws exception::OutOfRange if 'N' is zero, indicating an invalid number of Kraus operators
 * @throws exception::DimsInvalid if 'D' is zero, indicating an invalid matrix size
 *
 * @example
 * // usage of the randkraus function to generate a set of 2 random 2x2 kraus operators
 * int N = 2;
 * int D = 2;
 * std::vector<cmat> random_kraus_operators = randkraus(N, D);
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
 * @brief generate random hermitian matrix
 *
 * this function generates a random hermitian matrix of size 'D'. the generated matrix is hermitian
 * meaning that is equal to its own conjugate transpose. the matrix is constructed by first
 * generating a random complex matrix 'H' of size 'D x D' with entries uniformly distributed in the
 * interval [0, 1], then, its is transformed into hermitian matrix using formula H = H + H^†, where
 * H^† is the conjugate transpose of 'H'
 *
 * @param D the size of the hermitian matrix, its must be a positive inetger
 * @return a random hermitian matrix of size 'DxD'
 *
 * @throws exception::DimsInvalid if 'D' is zero, indicating an invalid input size
 *
 * @example
 * // usage of randH function to generate random hermitian matrix of size 3x3
 * idx D = 3;
 * cmat randomHermitianMatrix = randH(D);
 */
inline cmat randH(idx D) {
  if (D == 0)
    throw exception::DimsInvalid("clara::randH()");
  // generate random complex matrix of size D x D with entries uniformly distributed in [0, 1]
  cmat H = 2 * rand<cmat>(D, D) - (1. + 1_i) * cmat::Ones(D, D);
  // make the matrix hermitian by adding its conjugate transpose
  return H + adjoint(H);
}

/**
 * @brief generated a random quantum ket (column vector)
 *
 * this function generate a random quantum ket (column vector) of size Dx1, representing a quantum
 * state the elements of the vector are generated with a normal distribution N(mean, sigma)m and the
 * resulting ket is normalized to have a unit norm
 *
 * @param D the dimension of the quantum ket (column vector)
 * @return Dx1 random quantum ket (rnomalized column vector)
 *
 * @throws exception::DimsInvalid if 'D' is zero, indicating an invalid vector size
 *
 * @example
 * // usage of the randket function to generate a random 3-dimensional quantim ket
 * int D = 3;
 * ket random_quantum_ket = randket(D);
 */
inline ket randket(idx D) {
  if (D == 0)
    throw exception::DimsInvalid("clara::randket()");
  ket kt = randn<cmat>(D, 1);
  return kt / norm(kt);
}

/**
 * @brief generate a random quantum density matrix
 *
 * this function generate a random quantum density matrix of size Dx1, representing a quantum state
 * the eelement of the density matrix are generated by first generating a random hermitian matrix
 * and its adjoint (conjugate transpose), and finally normalized to have a unit trace, ensuring it
 * represente a valid quantum state
 *
 * @param D the dimension of the quantum state
 * @return DxD random quantum density matrix
 *
 * @throws exception::DimsInvalid if 'D' is zero, indicating invalid matrix size
 *
 * @example
 * // usage of the randrho function to generate a random 3-dimensional quantum density matrix
 * int D = 3;
 * cmat random_density_matrix = randrho(3);
 */
inline cmat randrho(idx D) {
  if (D == 0)
    throw exception::DimsInvalid("clara::randrho()");
  cmat result = 10 * randH(D);
  result = result * adjoint(result);
  return result / trace(result);
}

/**
 * @brief generate a random uniformly distribution permutation
 *
 * this function generates a random permutation of size N using the Knuth shuffle method, which
 * ensures that all permutation are equally probable. the function first create a vector of size N
 * with elements increasing order. then, it use std::shuffle algorithm to randomly shuffle the
 * elements of the vector, generating a random permutation
 *
 * @param N the size of the permutation to be generated
 * @return vector representing a random permutation of size N
 *
 * @throws exception::PermInvalid if 'N' is zero, indicating an invalid permutation size
 *
 * @example
 * // usage of the randperm function to generate a random permutation of size 5
 * int N = 5;
 * std::vector<idx> random_permutation = randperm(N);
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
 * @brief generate random probability vector uniformly distributed over the probability simplex
 *
 * this function generates a random probability vector of size N, where the elements are uniformly
 * distributed over the probability simple. the probability simplex is the set of all probability
 * vector whose elements are non-negative and sum up to 1. the function first generates N random
 * numbers from an exponential distribution with a rate parameter of 1 these random numbers are then
 * normalized to ensure that they sum up to 1, thus forming a valid probability vector
 *
 * @param N the size of the probability vector to be generated
 * @return vector representing a random probability vector of size N
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
