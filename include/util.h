#ifndef UTIL_H_
#define UTIL_H_

#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "types.h"
#include "constants.h"

namespace clara {
// kronecker product
types::cmat kron(const types::cmat &A, const types::cmat &B);
types::cmat kron_list(const std::vector<types::cmat> & list);
types::cmat kron_n(const types::cmat &A, int n);
types::cmat ptrace2(const types::cmat &AB, const std::vector<size_t> dims);
types::cmat syspermute(const types::cmat &A, const std::vector<size_t> &dims, const std::vector<size_t> perm);
types::cmat ptrace(const types::cmat &A, const std::vector<size_t> &dims, const std::vector<size_t> &subsys);

types::cmat mat_pow(const types::cmat &A, const types::cplx z);
// matrix functional calculus
types::cmat mat_f(const types::cmat &A, types::cplx (*)(const types::cplx &));
// exponential 
types::cmat mat_exp(const types::cmat &A);
// random matrix with entries uniform [0, 1]
types::cmat rand(const size_t rows, const size_t cols);
types::cmat rand(const size_t rows);
// random sqrt matrix
types::cmat randn(const size_t rows, const size_t cols);
types::cmat randn(const size_t rows);

// random unitary matrix
types::cmat rand_unitary(const size_t size);

// display a complex eigen::matrix (types::cmat)
void disp(const types::cmat &A, std::ostream& os = std::cout, unsigned int precision = 4, double eps = 1e-16);
// save matrix in text file
void save_text(const types::cmat &A, const std::string& fname, size_t precision = 16);
// load matrix from text file
types::cmat load_text(const std::string& fname);
// save metrix to a binary file in double precision
void save(const types::cmat &A, const std::string& fname);

// load matrix from binary file
types::cmat load(const std::string& fname);

// reshape the columns of A and returns a cmat with m rows and n columns
types::cmat reshape(const types::cmat& A, size_t rows, size_t cols);

template<typename T>
inline void print_container(const T& x) {
  for (typename T::const_iterator it = x.begin(); it != x.end(); it++)
    std::cout<<*it<<" ";
  std::cout<<std::endl;
}

// inline index to multi-index
inline std::vector<size_t> n2multiidx(const size_t &n, const std::vector<size_t> &dims) {
  size_t numdims = dims.size();
  std::vector<size_t> result(numdims, 0);
  int _n = n;
  size_t maxn = 1;
  for (size_t i = 0; i < numdims; ++i)
    maxn *= dims[i];
  if (n > maxn - 1)
    throw std::runtime_error("number too large, out of bounds!");
  size_t tmp = 0;
  for (size_t i = 0; i < numdims; i++) {
    tmp = _n % static_cast<int>(dims[numdims - i - 1]);
    result[numdims - i - 1] = tmp;
    _n = _n / static_cast<int>(dims[numdims - i - 1]);
  }
  return result;
}

// multi index to integer index
inline size_t multiidx2n(const std::vector<size_t> &midx, const std::vector<size_t> &dims) {
  size_t numdims = dims.size();
  std::vector<size_t> part_prod(numdims, 1);
  size_t result = 0;
  for (size_t i = 0; i < numdims; i++)
    if (midx[i] >= dims[i])
      throw std::runtime_error(
        "sub index exceeds corresponding dimension"
      );

  for (size_t j = 0; j < numdims; j++)
    if (j == 0)
      part_prod[numdims - 1] = 1;
    else
      part_prod[numdims - j - 1] = part_prod[numdims - j] * dims[numdims - j];

  for (size_t i = 0; i < numdims; i++)
    result += midx[i] * part_prod[i];
  return result;
}

inline types::cmat transpose(const types::cmat& A) {
  return A.transpose();
}

// conjugate
inline types::cmat conjugate(const types::cmat& A) {
  return A.conjugate();
}

// adjoint
inline types::cmat adjoint(const types::cmat& A) {
  return A.adjoint();
}

// trace
inline types::cplx trace(const types::cmat& A) {
  return A.trace();
}

inline double norm(const types::cmat& A) {
  return A.norm();
}

}

#endif // UTIL_H_
