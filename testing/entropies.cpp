#include <iostream>

#include "../include/clara.h"

using namespace clara;

int main() {
  // create complex matrix `rho` using the predefined st.pb00 state
  cmat rho = st.pb00;

  // compute the partial trace over subsystem B and store it in `rhoA`
  cmat rhoA = ptrace(rho, {1});
  // display the original state matrix
  std::cout << "State: " << std::endl << disp(rho) << std::endl;
  // display the partial trace result
  std::cout << "Partial Trace over B:" << std::endl << disp(rhoA) << std::endl;

  // calculate and display von Neumann entropy of the partial trace result
  std::cout << "von-neumann entropy: " << entropy(rhoA) << std::endl;
  // calculate and display Tsallis-1 entropy for the given state matrix `rho`
  std::cout << "Tsallis-1 entropy: " << tsallis(rho, 1) << std::endl;
  // calculate and display Tsallis-2 entropy for the given state matrix `rho`
  std::cout << "Tsallis-2 entropy: " << tsallis(rho, 2) << std::endl;
  // calculate and display quantum mutual information between subsystem A and B
  std::cout << "Quantum mutual information Between A and B: " << qmutualinfo(rho, {0}, {1})
            << std::endl;
}
