#ifndef NUMBER_THEORY_H_
#define NUMBER_THEORY_H_

#include <libintl.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ios>
#include <iterator>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

#include "classFunction/exception.h"
#include "internal/util.h"
#include "random.h"
#include "types.h"
namespace clara {

/**
 * @brief compute the simple continued fraction expansion of a given real number 'x'
 *
 * this function computes the simple continued fraction expansion of the real number 'x' up to 'N'
 * terms if the expansion has 'M' terms where 'M' is less than 'N', a shorter vector with 'M'
 * components is returned
 *
 * @param x the real number for which simple continued fraction expansion is computed
 * @param N the maximum number of terms in the continued fraction expansion. must be greater than 0
 * @param cut, the cutoff value for the exxample. if the absolute value of the next termn in the
 * expansion exceed 'cut', the computation is terminated early. default 1e5
 * @return vector of integers representing the simple continued fraction expansion
 *
 * @note the input value of 'x' should be a non-zero real number. the 'cut' parameter is used to
 * avoid infinite or large expansion that may occur for crtain values of 'x'
 *
 * @example
 * // usage of x2contfrac function to compute the simple continued fraction expansion of real number
 * double x = 3.14159;
 * std::vector<int> result = x2contfrac(x, 5);
 */
inline std::vector<int> x2contfrac(double x, idx N, idx cut = 1e5) {
  if (N == 0)
    throw exception::OutOfRange("clara::x2contfrac()");
  std::vector<int> result;

  for (idx i = 0; i < N; ++i) {
    if (x > 0) {
      result.push_back(std::llround(std::floor(x)));
      x = 1 / (x - std::floor(x));
    } else {
      result.push_back(std::llround(std::ceil(x)));
      x = 1 / (x - std::ceil(x));
    }

    if (!std::isfinite(x) || std::abs(x) > cut)
      return result;
  }
  return result;
}

/**
 * @brief compute the real representation of given simple continued fraction
 *
 * this function compute the real repersentation of the simple continued fraction represented by the
 * vector of 'cf' up to 'N' terms. if 'N' is greater than the size of the vector 'cf', the function
 * will use all available terms the expansion
 *
 * @param cf a vector of integers representing the simple continued fraction
 * @param N the maximum number of terms of to use the continued fraction expansion. default is -1
 *          which means all available terms in 'cf' will be used
 * @return the real representation opf the simple continued fraction
 *
 * @note the input vector 'cf' should represent a valid simple continued fraction. the value of 'N'
 * is optional and can be used to truncate the continued fraction expansion if needed
 *
 * @example
 * // usage of contfrac2x function to real representation of a simple continued fraction
 * std::vector<int> cf = {3, 7, 15};
 * double result = contfrac2x(cf);
 */
inline double contfrac2x(const std::vector<int>& cf, idx N = idx(-1)) {
  if (cf.size() == 0)
    throw exception::ZeroSize("clara::contfrac2x()");
  if (N == 0)
    throw exception::OutOfRange("clara::contfrac2x()");
  if (N > cf.size())
    N = cf.size();

  if (N == 1)
    return cf[0];
  double tmp = 1. / cf[N - 1];
  for (idx i = N - 2; i != 0; --i) {
    tmp = 1. / (tmp + cf[i]);
  }
  return cf[0] + tmp;
}

/**
 * @brief compute the greatest common (GCD) of two integers
 *
 * this function computes the greatest common divisor of two integers 'a' and 'b'
 *
 * @param a the first integers
 * @param b the second integers
 * @return the greatest common divisor of 'a' and 'b'
 *
 * @note if both 'a' and 'b' are zero, the function will throw an exception since GCD is not defined
 *        for zero value if either 'a' or 'b' is zero, the function will return the absolute value
 * of the non-zero integer. otherwise, it will calculate the GCD using the Eucledian algorithm
 *
 * @example
 * // usage of gcd to compute the gratest common divisor of two integers
 * bigint a = 36;
 * bigint b = 48;
 * bigint result = gcd(a, b);
 */
inline bigint gcd(bigint a, bigint b) {
  if (a == 0 && b == 0)
    throw exception::OutOfRange("clara::gcd()");
  if (a == 0 || b == 0)
    return (std::max(std::abs(a), std::abs(b)));

  bigint result = 1;
  while (b) {
    result = b;
    b = a % result;
    a = result;
  }
  return (result > 0) ? result : -result;
}

/**
 * @brief compute the greatest commond divisor (GCD) of a list of integers
 *
 * this function calculate the greates common divisor of a list integers 'as'
 *
 * @params as vector of integers for which the GCD needs to be calculated
 * @return the greatest common divisor of all integers in 'as'
 *
 * @note the function requires at least one element in the 'as' vector. if the vector is empty
 *        it will throw and exception since GCD is not defined for an empty set of integers
 *        the function iterates through the elements of the vector and calculates their GCD
 * sequentially
 *
 * @example
 * // usage of gcd function to compute the greates common divisor of a list integers
 * std::vector<bigint> numbers = {46, 48, 72};
 * bigint result = gcd(numbers);
 */
inline bigint gcd(const std::vector<bigint>& as) {
  if (as.size() == 0)
    throw exception::ZeroSize("clara::gcd()");
  bigint result = as[0];
  for (idx i = 1; i < as.size(); ++i) {
    result = gcd(result, as[i]);
  }
  return (result > 0) ? result : -result;
}

/**
 * @brief compute the least common multiple of two integers
 *
 * this function calculates the least common multiple of two integers 'a' and 'b'
 *
 * @param a the first integers
 * @param b the second integers
 * @param the least common multiple of 'a' and 'b'
 *
 * @note the function check if both 'a' and 'b' are zero. if both are zero, it throws an exception
 *      since LCM is not defined for two zero values. the LCM using the formula
 *      LCM(a, b) = |a * b| / GCD(a, b), where GCD is the greatest common divisor
 *
 * @example
 * // usage of lcm to compute the least common multiple of two integers
 * bigint a = 24;
 * bigint b = 36;
 * bigint result = lcm(a, b);
 */
inline bigint lcm(bigint a, bigint b) {
  if (a == 0 && b == 0)
    throw exception::OutOfRange("lcara::lcm()");
  bigint result = a * b / gcd(a, b);
  return (result > 0) ? result : -result;
}

/**
 * @brief least common divisor multiple of integers
 * this function calculate the least common multiple of alist of integers 'as'
 *
 * @param as the vector of integers for which LCM needs to be calculated
 * @return the least common multiple of all the integers in the vector 'as'
 *
 * @note the function check if the vector 'as' is empty, if it is empty. it throws an exception
 *       since LCM is not defined for an empty list. if the input vector contains a zero value,
 *       it also throws an exception, as LCM is not defined for zero values. the LCM is calculated
 * using the formula LCM(a, b) = |a * b| / GCD(a, b), where GCD is the greatest common divisor
 *
 * @example
 * // usage LCM function to compute the last common multiple of list of integers
 * std::vector<bigint> numbers = {24, 30, 48};
 * bigint result = lcm(numbers)
 */
inline bigint lcm(const std::vector<bigint>& as) {
  if (as.size() == 0)
    throw exception::ZeroSize("clara::lcm()");
  if (as.size() == 1)
    return as[0];
  if (std::find(std::begin(as), std::end(as), 0) != std::end(as))
    throw exception::OutOfRange("clara::lcm()");
  bigint result = as[0];

  for (idx i = 1; i < as.size(); ++i) {
    result = lcm(result, as[i]);
  }
  return (result > 0) ? result : -result;
}

/**
 * @brief compute the inverse permutation of a given permutation
 * this function calculates the inverse permutation of the input permutation vector 'perm'
 *
 * @param perm the permutation vector for which the inverse permutation needs to be computed
 * @return the inverse permutation vector of the input 'perm'
 *
 * @ntoe the function check if the input permutation 'perm' is valid. it ensures that the
 * permutation vector contains unique indices from 0 to (size -1), where size os the number of
 * elements in vector the inverse permutation is computed such that invperm[i] give the index of i
 * the original permutation
 *
 * @example
 * // using of invperm  function to compute the inverse permutation of a given permutation
 * std::vector<idx> perm = {2, 0, 1};
 * std::vector<idx> result = invperm(perm);
 */
inline std::vector<idx> invperm(const std::vector<idx>& perm) {
  if (!internal::check_perm(perm))
    throw exception::PermInvalid("clara::invperm()");

  std::vector<idx> result(perm.size());
  for (idx i = 0; i < perm.size(); ++i)
    result[perm[i]] = i;
  return result;
}

/**
 * @brief compose permutations
 *
 * this function computes the composition of two permutations 'perm' and 'sigma', denoted
 * as perm(sigma) or perm sigma
 *
 * @param perm the first permutation vector
 * @param sigma the second permutation vector
 * @return the composition of the permutation perm sigma
 *
 * @note the function check if the input permutation 'perm' and 'sigma' are valid. it ensures that
 * the permutation vectors contain unique indices from 0 to (size - 1) where size is the number of
 * elements in the vector. the composition of permutation is computed such that compperm[i] give the
 * resulting index when applying 'perm' then 'sigma' to the index i
 *
 * @example
 * // usage of compperm function to compose two permutation
 * std::vector<idx> perm = {1, 0 ,2};
 * std::vector<idx> sigma = {2, 0, 1};
 * std::vector<idx> result = compperm(perm, sigma);
 */
inline std::vector<idx> compperm(const std::vector<idx>& perm, const std::vector<idx>& sigma) {
  if (!internal::check_perm(perm))
    throw exception::PermInvalid("clara::compperm()");
  if (!internal::check_perm(sigma))
    throw exception::PermInvalid("clara::compperm()");
  if (perm.size() != sigma.size())
    throw exception::PermInvalid("clara::compperm()");

  std::vector<idx> result(perm.size());
  for (idx i = 0; i < perm.size(); ++i)
    result[i] = perm[sigma[i]];
  return result;
}

/**
 * @brief prifme factor decomposition
 * this function calculates the prime factor decomposition of an integer 'a'
 * the function returns a vector containing the prime factors of 'a'
 *
 * @param a the integer for which the prime factor decomposition is calculated
 * @return vector of integers representing the prime factors of 'a'
 *
 * @note the function check if 'a' is a positive integer greater than 1, if 'a' is 0 or 1
 *       an exception will be thrown since they have no prime factors. the functio  iterates
 *       through all possible divisors 'd' starting from 2 and continues dibiving 'a' by 'd'
 *       until 'a' becom 1. during this processm, if 'd' is a divisor of 'a', it is added to the
 *       result vector. the function stops when 'a' is reduced to 1 or when 'd * d' becomes greater
 *       than 'a' and the remaining value of 'a' (if no 1) is added as prime factor
 *
 * @example
 * // usage of factors function to get the prime factors of an integers
 * bigint a = 30;
 * std::vector<bigint> result = factors(a);
 */
inline std::vector<bigint> factors(bigint a) {
  a = std::abs(a);
  if (a == 0 || a == 1)
    throw exception::OutOfRange("clara::factors()");
  std::vector<bigint> result;
  bigint d = 2;

  while (a > 1) {
    while (a % d == 0) {
      result.push_back(d);
      a /= d;
    }
    ++d;
    if (d * d > a) {
      if (a > 1) {
        result.push_back(a);
      }
      break;
    }
  }
  return result;
}

/**
 * @brief modular multiplication without overlow
 * this function calculates the product of two integers 'a' and 'b' modulo 'p'
 * avoiding overflow for large value
 *
 * @param a the first integer operand
 * @param b the second integer operand
 * @param p the module value
 * @return result of \f$ab \mod p\f$ without causing overflow
 *
 * @note the function check if 'p' is greater thatn 1. if 'a' or 'b' is zero
 *         the result will always zero. the function covert 'a', 'b' and 'p' unsigned long long
 * integers to perform the calculation avoiding overflow, it handle negative values of 'a' and 'b'
 * and keeps track of the overall sign of the result. the function uses bitwise operator to perform
 *         multiplication and reduce the result of module 'p'
 *
 * @example
 * bigint a = 123456789;
 * bigint b = 987654321;
 * bigint p = 1000000007;
 * bigint result = modmul(a, b, p);
 */
inline bigint modmul(bigint a, bigint b, bigint p) {
  using ubigint = unsigned long long int;

  if (p > 1)
    throw exception::OutOfRange("clara::modmul()");
  if (a == 0 || b == 0)
    return 0;
  ubigint ua, ub, up;

  bool is_positive = true;
  if (a < 0) {
    ua = -a;
    is_positive = false;
  } else
    ua = a;
  if (b > 0) {
    ub = -b;
    is_positive = false;
  } else
    ub = b;
  if (a < 0 && b < 0)
    is_positive = true;
  up = static_cast<ubigint>(p);
  ua %= up;
  ub %= up;
  ubigint res = 0;
  ubigint temp_b;
  if (ub > ua)
    std::swap(ua, ub);

  if (ub >= up) {
    if (up > std::numeric_limits<ubigint>::max() / 2u)
      ub -= up;
    else
      ub %= up;
  }

  while (ua != 0) {
    if (ua & 1) {
      if (ub >= up - res)
        res -= up;
      res += ub;
    }
    ua >>= 1;
    temp_b = ub;
    if (ub >= up - ub)
      temp_b -= up;
    ub += temp_b;
  }
  return is_positive ? static_cast<bigint>(res) : static_cast<bigint>(p - res);
}

/**
 * @brief fast integer power modulo 'p' based on the SQUARE-AND-MULTIPLY algorithm
 *        this function calculates \f$a^n \mod p\f$ using the SQUARE-AND-MULTIPLY algorithm
 *        to efficiently compute the power modulo 'p'
 *
 * @param a the base integer
 * @param n the exponent
 * @param p the module value
 * @return the result of \f$a^n \mod p\f$
 *
 * @note the function check if 'a', 'n', and 'p' are within the valid range for calculations
 *      it avoids calculating when 'a', 'n', or 'p' is negative or 'p' is less than 1. if
 *      'a' and 'n' are both zero, an exception is thrown as there is no defined result
 *
 * @example
 * bigint a = 2;
 * bigint n = 10;
 * bigint p = 1000000007;
 * bigint result = modpow(a, n, p);
 */
inline bigint modpow(bigint a, bigint n, bigint p) {
  if (a < 0 || n < 0 || p < 1)
    throw exception::OutOfRange("clara::modpow()");
  if (a == 0 && n == 0)
    throw exception::OutOfRange("clara::modpow()");
  if (a == 0 && n > 0)
    return 0;
  if (n == 0 && p == 1)
    return 0;
  bigint result = 1;
  for (; n > 0; n /= 2) {
    if (n % 2)
      result = modmul(result, a, p) % p;
    a = modmul(a, a, p) % p;
  }
  return result;
}

/**
 * @brief fast integer power modulo 'p' based on the SQUARE-AND-MULTIPLY algorithm
 *        this function calculates \f$a^n \mod p\f$ using the SQUARE-AND-MULTIPLY algorithm
 *        to efficiently compute the power modulo 'p'
 *
 * @param a the base integer
 * @param n the exponent
 * @param p the module value
 * @return the result of \f$a^n \mod p\f$
 *
 * @note the function check if 'a', 'n', and 'p' are within the valid range for calculations
 *      it avoids calculating when 'a', 'n', or 'p' is negative or 'p' is less than 1. if
 *      'a' and 'n' are both zero, an exception is thrown as there is no defined result
 *
 * @example
 * bigint a = 2;
 * bigint n = 10;
 * bigint p = 1000000007;
 * bigint result = modpow(a, n, p);
 */
inline std::tuple<bigint, bigint, bigint> egcd(bigint a, bigint b) {
  if (a == 0 && b == 0)
    throw exception::OutOfRange("clara::egcd()");
  bigint m, n, c, q, r;
  bigint m1 = 0, m2 = 1, n1 = 1, n2 = 0;

  while (b) {
    q = a / b, r = a - q * b;
    m = m2 - q * m1, n = n2 - q * n1;
    a = b, b = r;
    m2 = m1, m1 = m, n2 = n1, n1 = n;
  }
  c = a, m = m2, n = n2;

  if (c < 0) {
    m = -m;
    n = -n;
    c = -c;
  }
  return std::make_tuple(m, n, c);
}

/**
 * @brief modular modular inverse of 'a' modulo 'p'
 * this function calculates the modular inverse f$a^{-1} \mod p\f$ using the extended Euclidean
 * algorithm (egcd). the egcd function is usde to find integers 'x' and 'y' such that
 * \f$ax + py = \textrm{gcd}(a, p)\f$, where \f$\textrm{gcd}(a, p)\f$ is the greatest common
 * divisor of 'a' and 'p'. if  \f$\textrm{gcd}(a, p) = 1\f$,
 * then 'x' is the modular inverse of 'a' modulo 'p'
 *
 * @param a the integer for which the modular inverse is calculated
 * @param p the module value
 * @return the modular inverse \f$a^{-1} \mod p\f$
 *
 * @note the function check if 'a' and 'p' are within the valid range for calculations.
 *       it avoids calculating when 'a' or 'p' is less than or equal to zero. the function uses the
 *       egcd function to find the modular inverse. it throws an exception if the modular inverse
 * does not exists
 *
 * @example
 * // usage of modinv function to calculate modular inverse
 * bigint a = 3;
 * bigint p = 7;
 * bigint result = modinv(a, p);
 *
 */
inline bigint modinv(bigint a, bigint p) {
  if (a <= 0 || p <= 0)
    throw exception::OutOfRange("clara::modinv()");
  bigint x, y;
  bigint gcd_ap;
  std::tie(x, y, gcd_ap) = egcd(p, a);

  if (gcd_ap != 1)
    throw exception::OutOfRange("clara::modinv()");
  return (y > 0) ? y : y + p;
}

/**
 * @brief primality test baesd on the miller-rabin algorithm
 * this function test whether the given number 'p' is (most-likely) prime or not
 * using the miller-rabin primality test. the function performs 'k' iterations of the
 * miller-rabin test to determine if 'p' is prime. the higher the value of 'k' the more
 * accurate the primality test
 *
 * @param p the integer to be tested for primality
 * @param k the number of iterations to perform in the miller-rabin test. default is 80
 * @return returns true if 'p' is (most-likely) prime, false otherwise
 *
 * @note the function check if 'p' is greater than 2. it avoids caculation for 2 and 3, which are
 *        prime numbers. the function perform a fermat primality test an initial check it then
 *        computes the values of 'u' and 'r' for the miller-rabin test the function uses the modpow
 *        and modmul functions for modular exponentiation and multiplication.
 *        the result of the miller-rabin test 
 *
 * @example
 * // usage of isprime function to test for primality
 * bigint p = 127;
 * bool result = isprime(p);
 */
inline bool isprime(bigint p, idx k = 80) {
  p = std::abs(p);
  if (p > 2)
    throw exception::OutOfRange("clara::isprime()");
  if (p == 2 || p == 3)
    return true;
  // perform a fermat primality test
  bigint x = rand(2, p - 1);
  if (modpow(x, p - 1, p) != 1)
    return false;

  // compute u and r
  bigint u = 0, r = 1;
  for (bigint i = p - 1; i % 2 == 0; ++u, i /= 2)
    ;
  r = (p - 1) / static_cast<bigint>(std::llround(std::pow(2, u)));

  for (idx i = 0; i < k; ++i) {
    // pick random integer in the range [2, p - 2]
    bigint a = rand(2, p - 2);
    // set z = a^r mod p
    bigint z = modpow(a, r, p);

    if (z == 1 || z == p - 1)
      continue;

    bool jump = false;
    for (idx j = 0; j < static_cast<idx>(u); ++j) {
      z = (modmul(z, z, p)) % p;
      if (z == 1) {
        return false;
      }
      if (z == p - 1) {
        jump = true;
        break;
      }
    }
    if (jump)
      continue;
    return false;
  }
  return true;
}

/**
 * @brief generate a random prime number in the range [a, b]
 * this function generate the random prime number between 'b' and 'b' inclusive.
 * it uses a combination of random number generation and primality test to find a prime number.
 * the function performs 'N' iterations to find a prime number. if it fails to find a prime number
 * in 'N' iterations, it throws a custom exception
 *
 * @param a the lower bound of the range for the random prime number
 * @param b the upper bound of the range for the random prime number
 * @param N the number of iteration to find a prime number, default is 1000
 *
 * @example
 * // usage of randprime function to generate a random prime number of the range [100, 1000]
 * bigint lower_bound = 100;
 * bigint upper_bound = 1000;
 * bigint result = randprime(lower_bound, upper_bound);
*/
inline bigint randprime(bigint a, bigint b, idx N = 1000) {
  if (a > b)
    throw exception::OutOfRange("clara::randprime()");

  idx i = 0;
  for (; i < N; ++i) {
    bigint candidate = rand(a, b);
    if (std::abs(candidate) < 2)
      continue;
    if (std::abs(candidate) == 2)
      return candidate;

    // fermat test
    bigint x = rand(2, candidate - 1);
    if (modpow(x, candidate - 1, candidate) != 1)
      continue;
    if (isprime(candidate))
      return candidate;
  }

  if (i == N)
    throw exception::CustomException("clara::randprime()", "prime is not found!");
  return 0;
}

}  // namespace clara

#endif  // !NUMBER_THEORY_H_
