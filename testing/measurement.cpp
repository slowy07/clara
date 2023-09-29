#include <iostream>
#include <tuple>

#include "../include/clara.h"

using namespace clara;
int main() {
  // create quantum state psi initialize to the |00⟩ sate
  ket psi = 00_ket;
  // create a CNOT gate U by tensor product of hadamard (H) and identity gates
  cmat U = gt.CNOT * kron(gt.H, gt.Id2);
  // apply the gate U to the quantum state psi to produce new state
  ket result = U * psi;

  std::cout << "producing bell state:\n";
  std::cout << disp(result) << std::endl;
  
  // apply the X gate to the result state at qubit 1
  result = apply(result, gt.X, {1});
  std::cout << "producing bell state:\n";
  std::cout << disp(result) << std::endl;

  // measure the result state using the hadamard gate on qubit 0
  auto measured = measure(result, gt.H, {0});
  std::cout << "measurement result: " << std::get<0>(measured) << std::endl;
  std::cout << "probabilities: ";
  std::cout << disp(std::get<1>(measured), ", ") << std::endl;
  // print the result state for each ppossible measure outcome
  std::cout << "result state:\n";
  for (auto&& it : std::get<2>(measured))
    std::cout << disp(it) << "\n\n";
}
