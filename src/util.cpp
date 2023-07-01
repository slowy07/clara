#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include "../include/util.h"
#include "../include/types.h"
#include "../include/stat.h"
#include "../include/constants.h"

namespace clara {

types::cmat kron(const types::cmat &A, const types::cmat &B) {
  int Acols = A.cols();
  int Arows = A.rows();
  int Bcols = B.cols();
  int Brows = B.cols();
  
  types::cmat result(Arows * Brows, Acols * Bcols);

  for (int i = 0; i < Arows; i++)
    for (int j = 0; j < Acols; ++j)
      result.block(i * Brows, j * Bcols, Brows, Bcols) = A(i, j) * B;
  return result;
}

// kronecker product a list of matrices
types::cmat kron_list(const std::vector<types::cmat>& list) {
  types::cmat result;
  result = list[0];
  for (unsigned int i = 1; i < list.size(); i++)
    result = kron(result, (types::cmat) list[i]);
  return result;
}

// kronecker product a matrix with itself $n$ times
types::cmat kron_pow(const types::cmat &A, size_t n) {
  std::vector<types::cmat> list;
  for (size_t i = 0; i < n; i++)
    list.push_back(A);
  return kron_list(list);
}

// integer index to multi index
void _n2multiidx(const size_t n, const size_t numdims, const size_t *dims, size_t *result) {
  size_t maxn = 1;
  for (size_t i = 0; i < numdims; i++)
    maxn *= dims[i];
  if (n > maxn - 1)
    throw std::runtime_error("number to large, out of bound");

  size_t _n = n;
  for (size_t i = 0; i < numdims; i++) {
    result[numdims - i - 1] = _n % static_cast<int>(dims[numdims - i - 1]);
    _n = _n / static_cast<int>(dims[numdims - i - 1]);
  }
}

size_t _multiidx2n(const size_t *midx, size_t numdims, const size_t *dims) {
  for (size_t i = 0; i < numdims; i++)
    if (midx[i] >= dims[i])
      throw std::runtime_error("sub index exceeds corresponding dimension");
  size_t *part_prod = new size_t[numdims];
  part_prod[numdims - 1] = 1;
  for (size_t j = 1; j < numdims; j++)
    part_prod[numdims - j - 1] = part_prod[numdims - j] * dims[numdims - j];
  size_t result = 0;
  for (size_t i = 0; i < numdims; i++)
    result += midx[i] * part_prod[i];
  delete[] part_prod;
  return result;
}

// partial tace over subsystem
types::cmat ptrace2(const types::cmat &AB, const std::vector<size_t> dims) {
  size_t D = static_cast<size_t>(AB.rows());
  if (D != static_cast<size_t>(AB.cols()))
    throw std::runtime_error("matrix must be square");
  size_t DA = dims[0];
  size_t DB = dims[1];
  if (DA * DB != D)
    throw std::runtime_error(
      "product of partial dimension must equal the dimension of the matrix"
    );
  types::cmat result = types::cmat::Zero(DA, DA);
  for (size_t i = 0; i < DA; i++)
    #pragma omp parallel for
    for (size_t j = 0; j < DA; j++) {
      result(i, j) = trace(static_cast<types::cmat>(AB.block(i * DB, j * DB, DB, DB)));
    }
  return result;
}

// permute the subsys in cmat
types::cmat syspermute(const types::cmat &A, const std::vector<size_t> perm, const std::vector<size_t> &dims) {
  size_t dimA = static_cast<size_t>(A.rows());
  if (dimA != static_cast<size_t>(A.cols()))
    throw std::runtime_error("matrix must be square");

  const size_t numdims = dims.size();
  if (numdims != perm.size())
    throw std::runtime_error("invalid permutation size");
  
  // check consistent dimension
  size_t dimtotal = 1;
  for (size_t i = 0; i < numdims; i++)
      dimtotal *= dims[i];
  if (dimtotal != dimA)
    throw std::runtime_error("dimension list does not match the dimension of the matrix");

  // check permutation is valid
  std::vector<size_t> sort_perm  = perm;
  std::sort(sort_perm.begin(), sort_perm.end());
  for (size_t i = 0; i < numdims; i++) {
    if (sort_perm[i] != i) {
      throw std::runtime_error("not valid permutation");
    }
  }

  size_t tmp = 1;
  for (size_t i = 0; i < numdims; i++)
    tmp *= dims[i];
  if (tmp != dimA)
    throw std::runtime_error("dimension mismatch");

  types::cmat result(dimA, dimA);
  size_t *cdims = new size_t[numdims];
  size_t *cperm = new size_t[numdims];
  
  for (size_t i = 0; i < numdims; i++) {
    cdims[i] = dims[i];
    cperm[i] = perm[i];
  }
  size_t iperm = 0;
  size_t jperm = 0;

  for (size_t i = 0; i < dimA; i++)
    #pragma omp parallel for
    for (size_t j = 0; j < dimA; j++)
      _syspermute_worker(numdims, cdims, cperm, i, j, iperm, jperm, A, result);
  delete[] cdims;
  delete[] cperm;

  return result;
}

void _syspermute_worker(const size_t numdims, const size_t *cdims, const size_t *cperm, const size_t i, const size_t j, size_t &iperm, size_t &jperm, const types::cmat &A, types::cmat &result) {
  size_t *midxrow = new size_t[numdims];
	size_t *midxcol = new size_t[numdims];
	size_t *midxrowtmp = new size_t[numdims];
	size_t *midxcoltmp = new size_t[numdims];
	size_t *permdims = new size_t[numdims];

  for (size_t i = 0; i < numdims; i++)
    permdims[i] = cdims[cperm[i]];

  _n2multiidx(i, numdims, cdims, midxrow);
  _n2multiidx(j, numdims, cdims, midxcol);

  for (size_t k = 0; k < numdims; k++) {
    midxrowtmp[k] = midxrow[cperm[k]];
		midxcoltmp[k] = midxcol[cperm[k]];
  }

  iperm = _multiidx2n(midxrowtmp, numdims, permdims);
  jperm = _multiidx2n(midxcoltmp, numdims, permdims);
  result(iperm, jperm) = A(i, j);

  delete[] midxrow;
  delete[] midxcol;
  delete[] midxrowtmp;
  delete[] midxcoltmp;
  delete[] permdims;
}

// partial trace
types::cmat ptrace(const types::cmat &A, const std::vector<size_t> &subsys, const std::vector<size_t> &dims) {
  types::cmat result;
  std::vector<size_t> perdims;
  std::vector<size_t> subsyssort = subsys;

  // sort the subsys
  std::sort(subsyssort.begin(), subsyssort.end());

  size_t numsubsys = subsyssort.size();
  size_t numdims = dims.size();
  std::vector<size_t> perm(numdims, 0);

  // total dimension of the trace subsys
  size_t dimsubsys = 1;
  for (size_t i = 0; i < numsubsys; i++)
    dimsubsys *= dims[subsyssort[i]];

  size_t dimtotal = 1;
  for (size_t i = 0; i < numdims; i++)
    dimtotal *= dims[i];

  std::vector<size_t> sizeAB;
  sizeAB.push_back(dimtotal / dimsubsys);
  sizeAB.push_back(dimsubsys);

  size_t cnt0 = 0;
  size_t cnt1 = 0;
  for (size_t i = 0; i < numdims; i++) {
    if (std::find(subsyssort.begin(), subsyssort.end(), i) != subsyssort.end()) {
      perm[numdims - numsubsys + cnt0] = i;
      cnt0++;
    } else {
      perm[cnt1] = i;
      cnt1++;
    }
  }
  
  return ptrace2(syspermute(A, perm, dims), sizeAB);
}

types::cmat mat_pow(const types::cmat &A, const types::cplx z) {
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A);
  types::cmat evects = es.eigenvectors();
  types::cmat evals = es.eigenvalues();
  for (int i = 0; i < evals.rows(); i++)
    evals(i) = std::pow(static_cast<types::cplx>(evals(i)), static_cast<types::cplx>(z));
  types::cmat evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

/*
* matrix functional calculus
* compute f(A), where f(*f) is function pointer
*/
types::cmat mat_f(const types::cmat &A, types::cplx (*f)(const types::cplx &)) {
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A);
  types::cmat evects = es.eigenvectors();
  types::cmat evals = es.eigenvalues();
  for (int i = 0; i < evals.rows(); i++)
    evals(i) = (*f)(evals(i));
  types::cmat evalsdiag = evals.asDiagonal();
  return evects * evalsdiag * evects.inverse();
}

types::cmat mat_exp(const types::cmat &A) {
  return mat_f(A, std::exp);
}

// random matrix with entries in uniform
types::cmat rand(const size_t rows, const size_t cols) {
  return Eigen::MatrixXcd::Random(rows, cols);
}

// random square matrix with entries in uniform
types::cmat rand(const size_t rows) {
  return rand(rows, rows);
}

// random matrix with entries
types::cmat randn(const size_t rows, const size_t cols) {
  stat::NormalDistribution nd;
  types::cmat A(rows, cols);
  double re, im;

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++) {
      re = nd.sample();
      im = nd.sample();
      A(i, j) = re + ct::ii * im;
    }
  return A;
}

types::cmat randn(const size_t rows) {
  return randn(rows, rows);
}

// random unitary matrix
types::cmat rand_unitary(const size_t size) {
  types::cmat H = randn(size);
  H = (H + adjoint(H)) / 2;
  return mat_exp(static_cast<types::cmat>(ct::ii * H));
}

void disp(const types::cmat &A, std::ostream& os, unsigned int precision, double eps) {
  if (A.rows() * A.cols() == 0) {
    os << "empty [" << A.rows() << " x " << A.cols() << "] matrix" << std::endl;
    return;
  }
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

  std::vector<size_t> maxlengthcols(A.cols(), 0);
  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      if (vstr[i * A.cols() + j].size() > maxlengthcols[j])
        maxlengthcols[j] = vstr[i * A.cols() + j].size();

  for (int i = 0; i < A.rows(); i++) {
    os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right;
    for (int j = 0; j < A.cols(); j++)
      os << std::setw(static_cast<int>(maxlengthcols[j] + 2)) << std::right << vstr[i * A.cols() + j];
    if (i < A.rows() - 1)
      os << std::endl;
  }
}

void disp(const types::cvect &v, std::ostream& os, unsigned int precision, double eps) {
  for (size_t i = 0; i < static_cast<size_t>(v.size() - 1); i++) {
    disp(v(i), os, precision, eps);
    std::cout<<" ";
  };
  disp(v(v.size() - 1), os, precision, eps);
}


void disp(const types::cplx &c, std::ostream& os, unsigned int precision, double eps) {
  types::cmat tmp(1, 1);
  tmp(0, 0) = c;
  disp(tmp, os, precision, eps);
}


// save matrix in text file
void save_text(const types::cmat & A, const std::string& fname, size_t precision) {
  size_t rows = A.rows();
  size_t cols = A.cols();
  std::ofstream fout(fname.c_str());
  if (!fout.is_open()) {
    throw std::runtime_error("error writing output file \"" + fname + "\"!");
  }

  fout << rows << " " << cols << std::endl;
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++)
      fout << std::setprecision(precision) << A(i, j);
    fout << std::endl;
  }
  fout.close();
}

// save matrix to a binary file in double precision
void save(const types::cmat & A, const std::string& fname) {
  std::fstream fout;
  fout.open(fname.c_str(), std::ios::out | std::ios::binary);
  if (fout.fail()) {
    throw std::runtime_error("error writing output file \"" + std::string(fname) + "\"!");
  }
  const char _header[] = "TYPES::CMAT";
  fout.write(_header, sizeof(_header));

  size_t rows = A.rows();
  size_t cols = A.cols();
  fout.write((char*) &rows, sizeof(rows));
  fout.write((char*) &cols, sizeof(cols));
  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      fout.write((char*) &A(i, j), sizeof(A(i, j)));
  fout.close();
}

// load matrix from binary file
types::cmat load(const std::string& fname) {
  std::fstream fin;
  fin.open(fname.c_str(), std::ios::in | std::ios::binary);
  if (fin.fail()) {
    throw std::runtime_error("error opening input file \"" + std::string(fname) + "\"!");
  }
  const char _header[] = "TYPES::CMAT";
  char *_fheader = new char[sizeof(_header)];

  // read the header
  fin.read(_fheader, sizeof(_header));
  if (strcmp(_fheader, _header)) {
    delete[] _fheader;
    throw std::runtime_error("input file \"" + std::string(fname) + "\" corrupted!");
  }
  delete[] _fheader;
  size_t rows, cols;
  fin.read((char*) &rows, sizeof(rows));
  fin.read((char*) &cols, sizeof(cols));

  types::cmat A(rows, cols);
  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      fin.read((char*) &A(i, j), sizeof(A(i, j)));
  fin.close();
  return A;
}

types::cmat reshape(const types::cmat& A, size_t rows, size_t cols) {
  size_t rowsA = A.rows();
  size_t colsA = A.cols();
  if (rowsA * colsA != rows * cols)
    throw std::runtime_error("dimension mismatch, cannot reshape!");
  Eigen::MatrixXd realA = A.real();
  Eigen::MatrixXd imagA = A.imag();

  realA = Eigen::Map<Eigen::MatrixXd>(realA.data(), rows, cols);
  imagA = Eigen::Map<Eigen::MatrixXd>(imagA.data(), rows, cols);
  return realA.cast<types::cplx>() + ct::ii * imagA.cast<types::cplx>();
}


}
