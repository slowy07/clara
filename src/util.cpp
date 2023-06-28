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
types::cmat kronn(const types::cmat &A, size_t n) {
  std::vector<types::cmat> list;
  for (size_t i = 0; i < n; i++)
    list.push_back(A);
  return kron_list(list);
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
  types::cmat Aij = types::cmat::Zero(DB, DB);
  types::cmat result = types::cmat::Zero(DA, DA);
  for (size_t i = 0; i < DA; i++)
    for (size_t j = 0; j < DA; j++) {
      Aij = AB.block(i * DB, j * DB, DB, DB);
      result(i, j) = trace(Aij);
    }
  return result;
}

// permute the subsys in cmat
types::cmat syspermute(const types::cmat &A, const std::vector<size_t> &dims, const std::vector<size_t> perm) {
  size_t dimA = static_cast<size_t>(A.rows());
  if (dimA != static_cast<size_t>(A.cols()))
    throw std::runtime_error("matrix must be square");

  size_t numdims = dims.size();
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
  std::vector<size_t> midxrow(numdims);
  std::vector<size_t> midxcol(numdims);
	std::vector<size_t> midxrowtmp(numdims);
	std::vector<size_t> midxcoltmp(numdims);
	std::vector<size_t> permdims(numdims);
  size_t iperm = 0;
  size_t jperm = 0;
  
  for (size_t i = 0; i < numdims; i++)
    permdims[i] = dims[perm[i]];

  for (size_t i = 0; i < dimA; i++)
    for (size_t j = 0; j < dimA; j++) {
      midxrow = n2multiidx(i, dims);
			midxrowtmp = midxrow;
			midxcol = n2multiidx(j, dims);
			midxcoltmp = midxcol;
      for (size_t k = 0; k < numdims; k++) {
        midxrowtmp[k] = midxrow[perm[k]];
        midxcoltmp[k] = midxcol[perm[k]];
      }
      iperm = multiidx2n(midxrowtmp, permdims);
      jperm = multiidx2n(midxcoltmp, permdims);
      result(iperm, jperm) = A(i, j);
    }
  return result;
}

// partial trace
types::cmat ptrace(const types::cmat &A, const std::vector<size_t> &dims, const std::vector<size_t> &subsys) {
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
  result = syspermute(A, dims, perm);
  result = ptrace2(result, sizeAB);
  
  return result;
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
    os << "Empty ["<< A.rows() << " x " << A.cols() << "] matrix." << std::endl;
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
      } else if (std::abs(im) < eps) {
        ostr << std::setprecision(precision) << std::fixed << re;
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

  // determine the maxium length of the entries in each column
  std::vector<size_t> maxlengthcols(A.cols(), 0);
  
  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      if (vstr[i * A.cols() + j].size() > maxlengthcols[j])
        maxlengthcols[j] = vstr[i * A.cols() + j].size();

  for (int i = 0; i < A.rows(); i++) {
    // display first column
    os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right << vstr[i * A.cols()];
    for (int j = 1; j < A.cols(); j++)
      os << std::setw(static_cast<int>(maxlengthcols[j] + 2)) << std::right << vstr[i * A.cols() + j];
    os << std::endl;
  }
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

// load matrix from text file
types::cmat load_text(const std::string& fname) {
  std::ifstream fin(fname.c_str());
  if (!fin.is_open()) {
    throw std::runtime_error("error opening input file \"" + fname  + "\"!");
  }
  size_t rows, cols;

  fin >> rows;
  fin >> cols;
  types::cmat A(rows, cols);
  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      fin >> A(i, j);
  fin.close();
  return A;
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
