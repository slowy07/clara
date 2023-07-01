#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "types.h"
#include "constants.h"

namespace clara {


inline types::cmat transpose(const types::cmat& A) {
  return A.transpose();
}

inline types::cmat conjugate(const types::cmat& A) {
  return A.conjugate();
}

inline types::cmat adjoint(const types::cmat& A) {
  return A.adjoint();
}

inline types::cplx trace(const types::cmat& A) {
  return A.trace();
}

inline double norm(const types::cmat& A) {
  return A.norm();
}

inline types::cvect evals(const types::cmat& A) {
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A);
  return es.eigenvalues();
}

inline types::cmat evects(const types::cmat &A) {
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A);
  return es.eigenvectors();
}

// kronecker product
types::cmat kron(const types::cmat &A, const types::cmat &B);
types::cmat kron_list(const std::vector<types::cmat> & list);
types::cmat kron_pow(const types::cmat &A, int n);

void _n2multiidx(const size_t n, const size_t numdims, const size_t *dims, size_t *result);
// multi index to integer index
size_t _multiidx2n(const size_t *midx, const size_t numdims, const size_t *dims);

// partial trace over subsys
types::cmat ptrace2(const types::cmat &AB, const std::vector<size_t> dims);

// permute the subsystem in cmat
types::cmat syspermute(const types::cmat &A, const std::vector<size_t> perm, const std::vector<size_t> &dims);

// partial trace
types::cmat ptrace(const types::cmat &A, const std::vector<size_t> &subsys, const std::vector<size_t> &dims);

// matrix power A^z
types::cmat mat_pow(const types::cmat &A, const types::cplx z);

/*
* matrix exponential calculus
* compute f(A), where f() is the function pointer
*/
types::cmat mat_f(const types::cmat &A, types::cplx (*)(const types::cplx &));

// matrix exponential
types::cmat mat_exp(const types::cmat &A);

// random matrix with entries in uniform
types::cmat rand(const size_t rows, const size_t cols);

// random square matrix with entries in uniform
types::cmat rand(const size_t rows);

// random matrix with entries in normal
types::cmat randn(const size_t rows, const size_t cols);

// random square matrix entries in normal
types::cmat randn(const size_t rows);

// random unitary matrix
types::cmat rand_unitary(const size_t size);

// display the complex eigen matrix
void disp(const types::cmat &A, std::ostream& os = std::cout, unsigned int precision = 4, double eps = 1e-16);
// display a complex vector in frienly form
void disp(const types::cvect &v, std::ostream& os = std::cout, unsigned int precision = 4, double eps = 1e-16);
// display a complex number in  friendly form
void disp(const types::cplx &c, std::ostream& os = std::cout, unsigned int precision = 4, double eps = 1e-16);

// display an integer vector
void disp(const types::ivect &v, std::ostream& os = std::cout);

// save save matrix to binaru file in double precision
void save(const types::ivect &v, std::ostream& os = std::cout);
// load to binary file in double precision
types::cmat load(const std::string& fname);

/*
* reshape the columns of A and return a cmat with rows and n columns
* use column-major order
*/
types::cmat reshape(const types::cmat& A, size_t rows, size_t cols);

// inline template
template<typename T>
inline void print_container(const T& x) {
  for (typename T::const_iterator it = x.begin(); it != x.end(); it++)
    std::cout<<*it<<" ";
  std::cout<<std::endl;
}

// used inside the #pragma omp parallel for syspermute
void _syspermute_worker(const size_t numdims, const size_t *cdims, const size_t *cperm, const size_t i, const size_t j, size_t &iperm, size_t &jprem, const types::cmat &A, types::cmat &result);

}


#endif // UTIL_H_
