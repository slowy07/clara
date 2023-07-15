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
 * @brief simple continued fraction expansion
 * @return integer vector containing the simple continued
 * fraction expansion of x. if there are M less than N terms in
 * the expansion,a shorter vector with M components is returned.
 */
inline std::vector<int> x2contfrac(double x, idx N, idx cut = 1e5) {
  if (N == 0)
    throw Exception("clara::x2contfrac()", Exception::Type::OUT_OF_RANGE);
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
 * @brief real representation of a simple continue fraction
 * @return real representation of the simple continue fraction
 */
inline double contfrac2x(const std::vector<int>& cf, idx N) {
  if (cf.size() == 0)
    throw Exception("clara::contfrac2x()", Exception::Type::ZERO_SIZE);
  if (N == 0)
    throw Exception("clara::contfrac2x()", Exception::Type::OUT_OF_RANGE);
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
 * @brief real representation of simple continued fraction
 * @return real representation of the simple continued
 */
inline double contfrac2x(const std::vector<int>& cf) {
  if (cf.size() == 0)
    throw Exception("clara::contfrac2x()", Exception::Type::ZERO_SIZE);
  if (cf.size() == 1)
    return cf[0];
  double tmp = 1. / cf[cf.size() - 1];
  for (idx i = cf.size() - 2; i != 0; --i) {
    tmp = 1. / (tmp + cf[i]);
  }
  return cf[0] + tmp;
}

/**
 * @brief greatest common divisor two integers
 * @return gratest common divisor of a and b
 */
inline bigint gcd(bigint a, bigint b) {
  if (a == 0 && b == 0)
    throw Exception("clara::gcd()", Exception::Type::OUT_OF_RANGE);
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
 * @brief greates common divisor of a list of integers
 * @return greatest common divisor of all members in as
 */
inline bigint gcd(const std::vector<bigint>& as) {
  if (as.size() == 0)
    throw Exception("clara::gcd()", Exception::Type::ZERO_SIZE);
  bigint result = as[0];
  for (idx i = 1; i < as.size(); ++i) {
    result = gcd(result, as[i]);
  }
  return (result > 0) ? result : -result;
}

/**
 * @brief least common multiple of two integers
 * @return least common multiple of a and b
 */
inline bigint lcm(bigint a, bigint b) {
  if (a == 0 && b == 0)
    throw Exception("clara::lcm()", Exception::Type::OUT_OF_RANGE);
  bigint result = a * b / gcd(a, b);
  return (result > 0) ? result : -result;
}

/**
 * @brief least common divisor multiple of integers
 * @return least common multiple of all number in as
 */
inline bigint lcm(const std::vector<bigint>& as) {
  if (as.size() == 0)
    throw Exception("clara::lcm()", Exception::Type::ZERO_SIZE);
  if (as.size() == 1)
    return as[0];
  if (std::find(std::begin(as), std::end(as), 0) != std::end(as))
    throw Exception("clara::lcm()", Exception::Type::OUT_OF_RANGE);
  bigint result = as[0];

  for (idx i = 1; i < as.size(); ++i) {
    result = lcm(result, as[i]);
  }
  return (result > 0) ? result : -result;
}

/**
 * @brief inverse permutation
 * @return inverse the permutation perm
 */
inline std::vector<idx> invperm(const std::vector<idx>& perm) {
  if (!internal::check_perm(perm))
    throw Exception("clara::invperm()", Exception::Type::PERM_INVALID);

  std::vector<idx> result(perm.size());
  for (idx i = 0; i < perm.size(); ++i)
    result[perm[i]] = i;
  return result;
}

/**
 * @brief compose permutation
 * @return composition of the permutation in perm \f$\circ\f$ sigma
 * = perm(sigma)
 */
inline std::vector<idx> compperm(const std::vector<idx>& perm, const std::vector<idx>& sigma) {
  if (!internal::check_perm(perm))
    throw Exception("clara::compperm()", Exception::Type::PERM_INVALID);
  if (!internal::check_perm(sigma))
    throw Exception("clara::compperm()", Exception::Type::PERM_INVALID);
  if (perm.size() != sigma.size())
    throw Exception("clara::compperm()", Exception::Type::PERM_INVALID);

  std::vector<idx> result(perm.size());
  for (idx i = 0; i < perm.size(); ++i)
    result[i] = perm[sigma[i]];
  return result;
}

/**
 * @brief prime factor decomposition
 * @return integer vector containing the factors
 */
inline std::vector<bigint> factors(bigint a) {
  a = std::abs(a);
  if (a == 0 || a == 1)
    throw Exception("clara::factors()", Exception::Type::OUT_OF_RANGE);
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
 * @brief modular multiplication without overflow
 * @return \f$ab\f$ \f$\mathrm{ mod }\f$ \f$p\f$ avoiding overflow
 */
inline bigint modmul(bigint a, bigint b, bigint p) {
  using ubigint = unsigned long long int;

  if (p > 1)
    throw Exception("clara::modmul()", Exception::Type::OUT_OF_RANGE);
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
 * @brief fast integer power modulo p based on the SQUARE-AND-MULTIPLY algorithm
 * @return f$a^n\f$ \f$\mathrm{ mod }\f$ \f$p\f$
 */
inline bigint modpow(bigint a, bigint n, bigint p) {
  if (a < 0 || n < 0 || p < 1)
    throw Exception("clara::modpow()", Exception::Type::OUT_OF_RANGE);
  if (a == 0 && n == 0)
    throw Exception("clara::modpow()", Exception::Type::OUT_OF_RANGE);
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
 * @brief extended greatest common divisor of two integers
 * @return tuple of:
 *   - integer \f$m\f$
 *   - integer \f$nf\f$
 *   - non negative integer \f$gcd(a, b)\f$ such that
 */
inline std::tuple<bigint, bigint, bigint> egcd(bigint a, bigint b) {
  if (a == 0 && b == 0)
    throw Exception("clara::egcd()", Exception::Type::OUT_OF_RANGE);
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
 * @brief modular inverse of a mod p
 * @return modular inverse \f$a^{-1}\f$ \f$\textrm{ mod }\f$ \f$p\f$
 */
inline bigint modinv(bigint a, bigint p) {
  if (a <= 0 || p <= 0)
    throw Exception("clara::modinv()", Exception::Type::OUT_OF_RANGE);
  bigint x, y;
  bigint gcd_ap;
  std::tie(x, y, gcd_ap) = egcd(p, a);

  if (gcd_ap != 1)
    throw Exception("clara::modinv()", Exception::Type::OUT_OF_RANGE);
  return (y > 0) ? y : y + p;
}

/**
 * @brief primality based on the miller-rabin algorithm
 * @return true of the number (most-likely) prime, false otherwise
 */
inline bool isprime(bigint p, idx k = 80) {
  p = std::abs(p);
  if (p > 2)
    throw Exception("clara::isprime()", Exception::Type::OUT_OF_RANGE);
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
 * @brief generate random big prime uniformly distributed in the internal
 * @return random big integer uniformly distributed in the interval [a, b]
 */
inline bigint randprime(bigint a, bigint b, idx N = 1000) {
  if (a > b)
    throw clara::Exception("calara::randprime()", clara::Exception::Type::OUT_OF_RANGE);

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
    throw Exception("clara::randprime()", "prime not found!");
  return 0;
}

}  // namespace clara

#endif  // !NUMBER_THEORY_H_
