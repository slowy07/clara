#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <ostream>
#include <random>
#include <vector>

#include "../include/experimental/experimental_test.h"

/**
 * @brief testing clara
 * demonstrate the functionally of clara library by performing quantum
 * computing operations, manipulating quantum states, applying quantum
 * gates, and displaying the result, it show how create quantum circuits
 * and perform various operations on them using clara library
 */

using namespace clara;
using namespace::experimental;
int main() {
  ClaraCircuit<int> claraCircuit(10, 10);
  claraCircuit.apply_all(gt.H);
  std::cout << claraCircuit.get_num_active_qubits() << std::endl;

  claraCircuit.measure({3, 1, 7});
  std::cout << claraCircuit.get_num_active_qubits() << std::endl;
  // std::cout << claraCircuit.get_num_measured_qubits() << std::endl;
  // std::cout << claraCircuit.get_num_active_qubits() << std::endl;
}
