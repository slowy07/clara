#ifndef GATES_H_
#define GATES_H_

#include "types.h"

namespace clara {
namespace gt {
extern types::cmat2 H;
extern types::cmat2 Id2;
extern types::cmat2 X;
extern types::cmat2 Y;
extern types::cmat2 Z;
extern types::cmat2 T;
types::cmat2 Rtheta(double theta);

// two qubit gates
extern types::cmat4 CNOT;
extern types::cmat4 CP;
// controlled-U, for arbitrary U
extern types::cmat4 CU(const types::cmat2 &);

// three qubits gates
extern types::cmat TOF;

// generalized Z gate
extern types::cmat Zd(size_t);
// generalized fourier gate
extern types::cmat Fd(size_t);
// generalized X gate
extern types::cmat Xd(size_t);

extern types::cmat CUd(const types::cmat &);

int _init_gates();
}
}

#endif // !GATES_H_
