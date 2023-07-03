#ifndef GATES_H_
#define GATES_H_

#include "types.h"
#include "constants.h"
#include "util.h"

namespace clara {
namespace gt {

// one qubit gates

// hadamard matrix
extern Eigen::MatrixXcd H;
// identity matrix
extern Eigen::MatrixXcd Id2;
// x, y, z, s t gate
extern Eigen::MatrixXcd X;
extern Eigen::MatrixXcd Y;
extern Eigen::MatrixXcd Z;
extern Eigen::MatrixXcd S;
extern Eigen::MatrixXcd T;

// two qubit gates
extern Eigen::MatrixXcd CNOT;
extern Eigen::MatrixXcd CP;

// three qubits gates
extern types::cmat TOF;

// initialize gates
inline void _init_gates() {
  // initialize the constant and gates
  H = Id2 = X = Y = Z = S = T = Eigen::MatrixXcd::Zero(2,2);
  CNOT = CP = Eigen::MatrixXcd::Zero(4,4);
  TOF = Eigen::MatrixXcd::Zero(8,8);

  H << 1 / sqrt(2), 1 / sqrt(2), 1 / sqrt(2), -1 / sqrt(2);
  Id2 << 1, 0, 0, 1;
  X << 0, 1, 1, 0;
	Z << 1, 0, 0, -1;
	Y(0, 1) = -ct::ii;
	Y(1, 0) = ct::ii;
	S(0, 0) = 1;
	S(1, 1) = ct::ii;
	T(0, 0) = 1;
	T(1, 1) = exp(ct::ii * ct::pi / 4.0);
	CNOT(0, 0) = 1;
	CNOT(1, 1) = 1;
	CNOT(2, 3) = 1;
	CNOT(3, 2) = 1;
	CP(0, 0) = 1;
	CP(1, 1) = 1;
	CP(2, 2) = 1;
	CP(3, 3) = -1;
	TOF(0, 0) = 1;
	TOF(1, 1) = 1;
	TOF(2, 2) = 1;
	TOF(3, 3) = 1;
	TOF(4, 4) = 1;
	TOF(5, 5) = 1;
	TOF(6, 7) = 1;
	TOF(7, 6) = 1;
}

// gates with variable dimension

// one qubit gates
inline Eigen::MatrixXcd Rtheta(double theta) {
  const std::complex<double> i(0.0, 1.0);
  Eigen::MatrixXcd result(2, 2);
  result << 1.0, 0.0, 0.0, std::exp(i * theta);
  return result;
}

// two qubit gates
inline Eigen::MatrixXcd CU(const Eigen::MatrixXcd &U) {
  Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(4, 4);
  result(0, 0) = 1;
  result(1, 1) = 1;
  result.block(2, 2, 2, 2) = U;
  return result;
}

// one quDit gates
inline Eigen::MatrixXcd Zd(size_t D) {
  Eigen::MatrixXcd result(D, D);
  for (size_t i = 0; i < D; i++)
    result(i, i) = pow(ct::omega(D), i);
  return result;
}

inline Eigen::MatrixXcd Fd(size_t D) {
  Eigen::MatrixXcd result(D, D);
  result.setZero();
  double sqrtD = std::sqrt(static_cast<double>(D));
  for (size_t i = 0; i < D; i++) {
    for (size_t j = 0; j < D; j++)
      result(i, j) = 1.0 / sqrtD * std::pow(ct::omega(D), static_cast<double>(i * j));
  }
  return result;
}

inline Eigen::MatrixXcd Xd(size_t D) {
  return Fd(D) * Zd(D) * Fd(D).inverse();
}

// two qudit gates
inline Eigen::MatrixXcd CUd(const Eigen::MatrixXcd &U) {
  // retrieves the dimension from the size of U
  size_t D = U.cols();
  Eigen::MatrixXcd result(D * D, D * D);
  result = Eigen::MatrixXcd::Zero(D * D, D * D);
  Eigen::MatrixXcd tmp(D, D);
  tmp = Eigen::MatrixXcd::Zero(D, D);

  for (size_t i = 0; i < D; i++) {
    if (i > 0)
      tmp(i - 1, i - 1) = 0;
    tmp(i, i) = 1;
    result += kron(tmp, mpower(U, i));
  }
  return result;
}

}
}

#endif // !GATES_H_
