#include <iostream>

#include "../include/clara.h"

using namespace clara;

// define constant for kraus operators
const cmat PZ0 = st.pz0;
const cmat PZ1 = st.pz1;

void ChannelOperations() {
  // create an initial quantum state
  cmat rho = st.pb00;
  std::cout << "initial state " << std::endl;
  std::cout << disp(rho) << std::endl;

  // perform partial transpose of first subsystem
  cmat rhoTA = ptranspose(rho, {0});
  std::cout << "eigenvalues of the partial transpose of bell-0 state are " << std::endl;
  std::cout << disp(transpose(hevals(rhoTA))) << std::endl;

  // set up measurement channel with 2 kraus operators
  std::cout << "measurement channel with 2 kraus operators" << std::endl;
  std::vector<cmat> Ks{PZ0, PZ1};
  std::cout << disp(Ks[0]) << "\nand\n" << disp(Ks[1]) << std::endl;

  // compute the superoperator matrix of the channel
  std::cout << "superoperator matrix of channel: " << std::endl
            << disp(kraus2super(Ks)) << std::endl;

  // partially trace down the second subsystem
  std::cout << "choi matrix of the channels: " << std::endl;
  std::cout << disp(kraus2choi(Ks)) << std::endl;

  // apply the measurement channel into the first subsystem
  cmat rhoOut = apply(rho, Ks, {0});
  std::cout << "after applying the measurement channel on the first qubit" << std::endl;
  std::cout << disp(rhoOut);

  // partial trace down the second subsystem
  cmat rhoA = ptrace(rhoOut, {1});
  std::cout << "after partially tracing down the second subsystem " << std::endl
            << disp(rhoA) << std::endl;

  // compute the von-neumann entropy of the resulting state
  double entropies = entropy(rhoA);
  std::cout << "entropy: " << entropies << std::endl;
}

int main() { ChannelOperations(); }
