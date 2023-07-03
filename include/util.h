#ifndef UTIL_H_
#define UTIL_H_

#include <eigen3/Eigen/QR>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "types.h"
#include "constants.h"
#include "internal.h"
#include "stat.h"

namespace clara {


template<typename Derived>
Derived transpose(const Eigen::MatrixBase<Derived>& A) {
  return A.transpose();
}

template<typename Derived>
Derived conjugate(const Eigen::MatrixBase<Derived>& A) {
  return A.conjugate();
}

template<typename Derived>
Derived adjoint(const Eigen::MatrixBase<Derived>& A) {
  return A.adjoint();
}

// trace, perserve return type
template<typename Derived>
typename Derived::Scalar trace(const Eigen::MatrixBase<Derived>& A) {
  return A.trace();
}

// absolute values component-wise, return double
template<typename Derived>
Derived abs(const Eigen::MatrixBase<Derived>& A) {
  Derived result = Derived::Zero(A.rows(), A.cols());
  for (typename Derived::Index i = 0; i < A.rows(); i++)
    for (typename Derived::Index j = 0; j < A.cols(); j++)
      result(i, j) = std::abs(A(i, j));
  return result;
}

// trace-norm
template<typename Derived>
typename Eigen::MatrixXd::Scalar norm(const Eigen::MatrixBase<Derived>& A) {
  // convert matrix to complex then return its norm
  return (A.template cast<types::cplx>()).norm();
}

// eigenvalues
template<typename Derived>
Eigen::MatrixXcd evals(const Eigen::MatrixBase<Derived>& A) {
  return (A.template cast<types::cplx>()).eigenvalues();
}

// eigenvectors
template<typename Derived>
Eigen::MatrixXcd evects(const Eigen::MatrixBase<Derived>& A) {
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A.template cast<types::cplx>());
  return es.eigenvectors();
}

// kronecker product
template<typename Derived>
Derived kron(const Eigen::MatrixBase<Derived> &A, const Eigen::MatrixBase<Derived> &B) {
  int Acols = A.cols();
  int Arows = A.rows();
  int Bcols = B.cols();
  int Brows = B.rows();
  
  Derived result;
  result.resize(Arows * Brows, Acols * Bcols);

  for (int i = 0; i < Arows; i++)
    for (int j = 0; j < Acols; j++)
      result.block(i * Brows, j * Bcols, Brows, Bcols) = A(i, j) * B;
  return result;
}

/**
 * kronecker prodoct of alist of matris, preserve return type
 * <Derived> is force to be a mterix by invocation of kron inside
 * the function
 */
template<typename Derived>
Derived kron_list(const std::vector<Derived> &list) {
  Derived result = list[0];
  for (size_t i = 1; i < list.size(); i++)
    result = kron(result, list[i]);
  return result;
}

/**
 * kronecker product of matrix with itself of $n$ times, preserve
 * return type
 */
template<typename Derived>
Derived kron_pow(const Eigen::MatrixBase<Derived> &A, size_t n) {
  std::vector<Derived> list;
  for (size_t i = 0; i < n; i++)
    list.push_back(A);
  return kron_list(list);
}

// matrix power A^z
template<typename Derived>
Eigen::MatrixXcd mpower(const Eigen::MatrixBase<Derived> &A, const types::cplx z) {
  // check square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: mpower matrix must be square");

  // define A^0
  if (real(z) == 0 && imag(z) == 0) {
    Eigen::MatrixXcd result(A.rows(), A.rows());
    result.setIdentity();
    return result;
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A.template cast<types::cplx>());
  Eigen::MatrixXcd evects = es.eigenvectors();
  Eigen::MatrixXcd evals = es.eigenvalues();
  for (int i = 0; i < evals.rows(); i++)
    evals(i) = std::pow(static_cast<types::cplx>(evals(i)), static_cast<types::cplx>(z));
  Eigen::MatrixXcd evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

/**
 * itneger matrix power, preserve return type
 * explictly multiply the matrix with itself
 */
template<typename Derived>
Derived mpower(const Eigen::MatrixBase<Derived> &A, const size_t n) {
  // check square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: mpower matrix must be square");
  Derived result = A;
  if (n == 0)
    return result.setIdentity();
  for (size_t i = 1; i < n; i++)
    result *= A;
  return result;
}

/**
 * matrix functional calculus
 *  f(A), where (*f) is the function pointer
 */
template<typename Derived>
Eigen::MatrixXcd funm(const Eigen::MatrixBase<Derived> &A, types::cplx (*f)(const types::cplx)) {
  // check square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: matrix must be square");
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A.template cast<types::cplx>());
  Eigen::MatrixXcd evects = es.eigenvectors();
  Eigen::MatrixXcd evals = es.eigenvalues();
  for (int i = 0; i < evals.rows(); i++)
    evals(i) = (*f)(evals(i));
  Eigen::MatrixXcd evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

// matrix exponential
template<typename Derived>
Eigen::MatrixXcd expm(const Eigen::MatrixBase<Derived> &A) {
  return funm(A, std::exp);
}

// random double matrix with entries in normal
inline Eigen::MatrixXd randn(const size_t rows, const size_t cols) {
  stat::NormalDistribution nd;
  Eigen::MatrixXd A(rows, cols);
  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      A(i, j) = nd.sample();
  return A;
}

// random square matrix with entries in nornal
inline Eigen::MatrixXd randn(const size_t rows) {
  return randn(rows, rows);
}

// random unitary matrix
inline Eigen::MatrixXcd rand_unitary(const size_t size) {
  Eigen::MatrixXcd x(size, size);
  
  x.real() = 1. / sqrt(2) * randn(size);
  x.imag() = 1. / sqrt(2) * randn(size);
  Eigen::HouseholderQR<Eigen::MatrixXcd> qr(x);
  Eigen::MatrixXcd Q = qr.householderQ();
  return Q;
}

/**
 * reshape the columns of A and returns a matrix with m rows and n columns
 * ise column major order
 */
template<typename Derived>
Derived reshape(const Eigen::MatrixBase<Derived>& A, size_t rows, size_t cols) {
  size_t rowsA = A.rows();
  size_t colsA = A.cols();

  if (rowsA * colsA != rows * cols)
    throw std::runtime_error("ERROR: dimension mismatch");
  return Eigen::Map<Derived>(static_cast<Derived>(A).data(), rows, cols);
}

// permutes the subsystem in a matrix
template<typename Derived>
Derived syspermute(const Eigen::MatrixBase<Derived> &A, const std::vector<size_t> perm, const std::vector<size_t> &dims)  {
  // check square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: syspermute matrix must be square");
  // check that dims is a valid dimension vector
  if (!internal::_check_dims(dims))
    throw std::runtime_error("ERROR: syspermute invalid dimension vector");
  // check that the size of the permutation is ok
  if (!internal::_check_perm(perm, dims))
    throw std::runtime_error("ERROR: syspermute invalid permutation size");

  size_t dim = static_cast<size_t>(A.rows());
  const size_t numdims = dims.size();
  size_t *cdims = new size_t[numdims];
  size_t *cperm = new size_t[numdims];
  Derived result(dim, dim);

  // copy dims in cdims and perm in cperm
  for (size_t i = 0; i < numdims; i++) {
    cdims[i] = dims[i];
    cperm[i] = perm[i];
  }
  size_t iperm = 0;
  size_t jperm = 0;
  for (size_t i = 0; i < dim; i++)
  #pragma omp parallel for
    for (size_t j = 0; j < dim; j++)
      internal::_syspermute_worker(numdims, cdims, cperm, i, j, iperm, jperm, A, result);
  delete[] cdims;
  delete[] cperm;
  return result;
}

// partial trace over subsystem B in a D_A x D_B system
template<typename Derived>
Derived ptrace2(const Eigen::MatrixBase<Derived> &A, const std::vector<size_t> dims) {
  // check square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: ptrace2 matrix must be square");
  // check dims has only 2 elements
  if (dims.size() != 2)
    throw std::runtime_error("ERROR: ptrace2 must habe only 2 dimension");
  // check that dim as valid dimension vector
  if (!internal::_check_dims(dims))
    throw std::runtime_error("ERROR: dimension vector does not match the dimension of the matrix");

  size_t DA = dims[0];
  size_t DB = dims[1];
  Derived result = Derived::Zero(DA, DA);
  
  for (size_t i = 0; i < DA; i++)
  #pragma omp parallel for
    for (size_t j = 0; j < DA; j++) {
      result(i, j) = trace(
        static_cast<Derived>(A.block(i * DB, j * DB, DB, DB))
      );
    }
  return result;
}

// partial trace
template<typename Derived>
Derived ptrace(const Eigen::MatrixBase<Derived> &A, const std::vector<size_t> &subsys, const std::vector<size_t> &dims) {
  // check square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: ptrace matrix must be a square");
  // check that dims is a valid dimension vector
  if (!internal::_check_dims(dims))
    throw std::runtime_error("ERROR: ptrace invalid dimension vector");
  // check that dims match dimension of A
  if (!internal::_check_dims_match_mat(dims, A))
    throw std::runtime_error("ERROR: ptrace dimension vector does not match the dimension of the matrix");
  if (!internal::_check_subsys(subsys, dims))
    throw std::runtime_error("ERROR: ptrace invalid subsystem");

  size_t dim = static_cast<size_t>(A.rows());
   // number of subsystems trace out
  size_t numsubsys = subsys.size();
   // total number of subsystems
  size_t numdims = dims.size();
  // the permutation vector
  std::vector<size_t> perm(numdims, 0);
  // the permuted dimensions
  std::vector<size_t> permdims;

  Derived result;
  
  // the total dimension of the taced-out subsystem
  size_t dimsubsys = 1;
  for (size_t i = 0; i < numsubsys; i++)
    dimsubsys *= dims[subsys[i]];
  std::vector<size_t> sizeAB;
	sizeAB.push_back(dim / dimsubsys);
	sizeAB.push_back(dimsubsys);

  // construct the permutation that bring the traced out to the end
  size_t cnt0 = 0;
  size_t cnt1 = 0;
  for (size_t i = 0; i < numdims; i++) {
    // finding that belong the subsys
    if (std::find(subsys.begin(), subsys.end(), i) != subsys.end()) {
      perm[numdims - numsubsys + cnt0] = i;
      cnt0++;
    } else {
      perm[cnt1] = i;
      cnt1++;
    }
  }
  return ptrace2(syspermute(A, perm, dims), sizeAB);
}

// partial transpose
template<typename Derived>
Derived ptranspose(const Eigen::MatrixBase<Derived> &A, const std::vector<size_t>& subsys, const std::vector<size_t>& dims) {
  // check the square matrix
  if (!internal::_check_square_mat(A))
    throw std::runtime_error("ERROR: ptranspose matrix must be square");
  // check that dims is a valid dimension vector
  if (!internal::_check_dims(dims))
    throw std::runtime_error("ERROR: ptranspose invalid dimension vector");
  // check that dims match the dimension of A
  if (!internal::_check_dims_match_mat(dims, A))
    throw std::runtime_error("ERROR: ptranspose dimension vectore does not match the dimension of the matrix");
  if (!internal::_check_subsys(subsys, dims))
    throw std::runtime_error("ERROR: ptranspose invalid subsystem");

  size_t dim = static_cast<size_t>(A.rows());
  const size_t numdims = dims.size();
  const size_t numsubsys = subsys.size();
  size_t *cdims = new size_t[numdims];
	size_t *midxrow = new size_t[numdims];
	size_t *csubsys = new size_t[numsubsys];

  Derived result = A;

  // copy dims in cdims and subsys in csubsys
  for (size_t i = 0; i < numdims; i++)
    cdims[i] = dims[i];
  for (size_t i = 0; i < numsubsys; i++)
    csubsys[i] = subsys[i];

  size_t iperm = 0;
  size_t jperm = 0;
  for (size_t i = 0; i < dim; i++) {
    internal::_n2multiidx(i, numdims, cdims, midxrow);
#pragma omp parallel for
    for (size_t j = 0; j < dim; j++)
      internal::_ptranspose_worker(midxrow, numdims, numsubsys, cdims, csubsys, i, j, iperm, jperm, A, result);
  }
  
  delete[] midxrow;
  delete[] cdims;
  delete[] csubsys;
  return result;
}

// random matrix with entries in uniform
template<typename Derived>
Derived rand(const size_t rows, const size_t cols) {
  return Derived::Random(rows, cols);
}

// random square matrix with entries in uniform
template<typename Derived>
Derived rand(const size_t rows) {
  return rand<Derived>(rows, rows);
}

// save matrix to a binary file in double precision
template<typename Derived>
void save(const Eigen::MatrixBase<Derived> & A, const std::string& fname) {
  std::fstream fout;
  fout.open(fname.c_str(), std::ios::out | std::ios::binary);
  if (fout.fail()) {
    throw std::runtime_error(
      "ERROR: writting output file \"" + std::string(fname) + "\" error!"
    );

    // write header to file
    const char _header[] = "TYPE::Eigen::Matrix";
    fout.write(_header, sizeof(_header));
    
    size_t rows = A.cols();
    size_t cols = A.cols();
    fout.write((char*) &rows, sizeof(rows));
    fout.write((char*) &cols, sizeof(cols));
    for (size_t i = 0; i < rows; i++)
      for (size_t j = 0; j < cols; j++)
        fout.write((char*) &A(i, j), sizeof(A(i, j)));
    fout.close();
  }
}

// load matrix from binary
template<typename Derived>
Derived load(const std::string& fname) {
  std::fstream fin;
  fin.open(fname.c_str(), std::ios::in | std::ios::binary);

  if (fin.fail()) {
    throw std::runtime_error(
      "ERROR: opening input file \"" + std::string(fname) + "\" error!"
    );

    // write header to file
    const char _header[] = "TYPE::Eigen::Matrix";
    // include the zero-terminating string
    char *_fheader = new char[sizeof(_header)];

    // read the header
    fin.read(_fheader, sizeof(_fheader));
    if (strcmp(_fheader, _header)) {
      delete[] _fheader;
      throw std::runtime_error(
        "ERROR: load input file \"" + std::string(fname) + "\" is corrupted!"
      );
    }
    delete[] _fheader;

    size_t rows, cols;
    fin.read((char*) &rows, sizeof(rows));
    fin.read((char*) &cols, sizeof(cols));
    
    Derived A(rows, cols);
    for (size_t i = 0; i < rows; i++)
      for (size_t j = 0; j < cols; j++)
        fin.read((char*) &A(i, j), sizeof(A(i, j)));
    fin.close();
    return A;
  }
}

// display a an Eigen::MatrixX in friendly form
template<typename Derived>
void disp(const Eigen::MatrixBase<Derived> &A, std::ostream& os = std::cout, unsigned int precision = 4) {
  std::cout << "typeid: " << typeid(Derived).name() << std::endl;
  os << std::setprecision(precision) << std::fixed << A;
}

template<>
inline void disp(const Eigen::MatrixBase<Eigen::MatrixXcd> &A, std::ostream& os, unsigned int precision) {
  if (A.rows() * A.cols() == 0) {
    os << "empty [" << A.rows() <<  " x " << A.cols() << "] matrixx";
    os << std::endl;
    return;
  };
  // everything smaller than eps is displayed as 0
  double eps = std::numeric_limits<double>::epsilon();
  std::ostringstream ostr;
  std::vector<std::string> vstr;
  std::string strtmp;

  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      strtmp.clear();
      ostr.clear();
      ostr.str(std::string());
      
      double re = static_cast<types::cplx>(A(i, j)).real();
      double im = static_cast<types::cplx>(A(i, j)).imag();

      if (std::abs(re) < eps && std::abs(im) < eps) {
        vstr.push_back("0 ");
      } else if (std::abs(re) < eps) {
        ostr << std::setprecision(precision) << std::fixed << im;
        vstr.push_back(ostr.str() + "i");
      } else if (std::abs(im) < eps) {
        ostr << std::setprecision(precision) << std::fixed << re;
        vstr.push_back(ostr.str() + " ");
      } else {
        ostr << std::setprecision(precision) << std::fixed << re;
        strtmp = ostr.str();
        strtmp += (im > 0 ? " + " : " - ");
        ostr.clear();
        ostr.str(std::string());
        ostr << std::setprecision(precision) << std::fixed << std::abs(im);
        strtmp += ostr.str();
        strtmp += "i";
        vstr.push_back(strtmp);
      }
    }
  }
  // determine the maximum length of the entries in each column
  std::vector<size_t> maxlengthcols(A.cols(), 0);
  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      if (vstr[i * A.cols() + j].size() > maxlengthcols[j])
        maxlengthcols[j] = vstr[i * A.cols() + j].size();
  
  // display
  for (int i = 0; i < A.rows(); i++) {
    os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right << vstr[i * A.cols()];
    for (int j = 1; j < A.cols(); j++)
      os << std::setw(static_cast<int>(maxlengthcols[j] + 2)) << std::right << vstr[i * A.cols() + j];
    if (i < A.rows() - 1)
      os << std::endl;
  }
}

// display complex number in friendly form
inline void disp(const types::cplx c, std::ostream& os = std::cout, unsigned int precision = 4) {
  Eigen::MatrixXcd tmp(1, 1);
  tmp(0, 0) = c;
  disp(tmp, os, precision);
}

}
#endif // UTIL_H_
