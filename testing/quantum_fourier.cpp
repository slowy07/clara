#include <omp.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "../include/clara.h"

using namespace clara;

int main(int argc, char **argv) {
  // check if the correct number of command line arguments is provided
  if (argc != 3) {
    std::cerr << "ERROR: please specify number of core, qubit" << std::endl;
    exit(EXIT_FAILURE);
  }

  // number of cpu core to use
  idx num_core = std::stoi(argv[1]);
  // number of qubits in the quantum system
  idx n = std::stoi(argv[2]);
  // set the number openMP threads to use
  omp_set_num_threads(num_core);

  // create a vector to store qubit indices
  std::vector<idx> qubits(n);
  // create a vector to store qubit indices
  ket psi = mket(qubits);
  // create a ket representing the intial quantum state
  ket result = psi;

  // measure the execution time usign the timer class
  Timer<> t;
  for (idx i = 0; i < n; ++i) {
    result = apply(result, gt.H, {i});
    for (idx j = 2; j <= n - i; ++j) {
      cmat Rj(2, 2);
      Rj << 1, 0, 0, omega(std::pow(2, j));
      // apply controlled-phase gate with control on qubit i and target on qubit i + j - 1
      result = applyCTRL(result, Rj, {i + j - 1}, {i});
    }
  }

  // apply SWAP gates to exchange qubit states
  for (idx i = 0; i < n / 2; ++i) {
    result = apply(result, gt.SWAP, {i, n - i - 1});
  }

  // print the number of cores, number of qubits and the execution time
  std::cout << "core: " << num_core << ", qubits" << n << ", time" << t.toc() << " seconds"
            << std::endl;
}
