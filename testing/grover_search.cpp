/**
 * @brief grove search algorithm
 * a quantum algorithm for searching an unsorted
 * database with N entries in O(N^(1/2)) time and using
 * O(logN) stroage space
 */

#include <cmath>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "../include/clara.h"

using namespace clara;

int main() {
  // represent the number of qubit used in algorithm
  clara::idx n = 4;
  std::cout << "grover on n " << n << "qubits" << std::endl;

  /**
   * local dimension is vector with `n` elements, each
   * elements representing the dimension of a qubit in the
   * quantum system, in this case, all qubits have 2 dimension
   * (standard quantum bit or qubit)
   */
  std::vector<clara::idx> dims(n, 2);
  /**
   * `N` calculated 2^n, representing the size of the quantum
   * state  or the database to be searched
   */
  clara::idx N = std::round(std::pow(2, n));
  std::cout << "database size " << N << std::endl;

  /**
   * random element `marked` from the database is selected
   * for demonstration purposes, the element is printed,
   * showing its binary representation using the function
   * clara::n2multiidx
   */
  clara::idx marked = clara::randidx(0, N - 1);
  std::cout << "mark state " << marked << "->";
  std::cout << clara::disp(clara::n2multiidx(marked, dims), " ") << std::endl;

  /**
   * computational state `psi` is initialized to `|0>^n`, meaning
   * all qubits are in the `|0>` state
   */
  clara::ket psi = clara::mket(clara::n2multiidx(0, dims));

  /**
   * hadamard gate `H` is applied to each qubit `n` times
   * using function `clara::kronpow`, where `gt.H`
   * is the hadamard gate, the storead in `psi`, representing
   * the superposition of all possible states
   */
  psi = (clara::kronpow(clara::gt.H, n) * psi).eval();
  /**
   * diffusion operator `G` is calculated by subtract the identity
   * matrix `gt.Id(N)` from the projection of `psi` onto itsself
   * (2 * prj(psi)). the diffusion operator amplifies the
   * amplitude of the marked element in the superposition of all
   * possible sates
   */
  clara::cmat G = 2 * clara::prj(psi) - clara::gt.Id(N);

  /**
   * number of queries is set based on the square root of the database
   * size `N`
   */
  clara::idx nqueries = std::ceil(clara::pi / 4 * std::sqrt(N));
  /**
   * the grover algorithm is executed by applying the oracle operation
   * (inverting the amplitude of the marked element) followed the
   * diffusion operator in a loop for `nqueries` times. this amplifies
   * the amplitude of the marked state in the superposition
   */
  std::cout << "run " << nqueries << " query" << std::endl;
  for (clara::idx i = 0; i < nqueries; ++i) {
    // apply oracle first, no aliasing
    psi(marked) = -psi(marked);
    // then the deiffusion operator
    psi = (G * psi).eval();
  }

  /**
   * the probability of the marked state is calculated and printed,
   * representing the probability of finding the marked element in
   * the superposition, the probability of all states in the superposition
   * are printed using the `clara::disp()` function
   */
  auto measured = clara::measure(psi, gt.Id(N));
  std::cout << "probability of marked state: " << std::get<1>(measured)[marked] << std::endl;
  std::cout << "probability of all results: ";
  std::cout << clara::disp(std::get<1>(measured), ", ") << std::endl;

  // sample
  std::cout << "sample " << std::endl;
  idx result = std::get<0>(measured);
  if (result == marked)
    std::cout << "correct result obrained ";
  else
    std::cout << "not yet obtained";
  std::cout << result << " -> ";
  std::cout << clara::disp(n2multiidx(result, dims), " ") << std::endl;
}
