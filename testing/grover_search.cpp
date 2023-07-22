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
  idx n = 4;
  std::cout << "grover on n " << n << "qubits" << std::endl;

  // local dimension
  std::vector<idx> dims(n, 2);
  idx N = std::round(std::pow(2, n));
  std::cout << "database size " << N << std::endl;

  // mark an element randomly
  idx marked = clara::randidx(0, N - 1);
  std::cout << "mark state " << marked << "->";
  std::cout << clara::disp(clara::n2multiidx(marked, dims), " ") << std::endl;

  // computational |0>^times n
  clara::ket psi = clara::mket(n2multiidx(0, dims));

  // apply H^times n, no aliasing
  psi = (kronpow(gt.H, n) * psi).eval();
  // diffusion operator
  cmat G = 2 * prj(psi) - gt.Id(N);

  // number queries
  idx nqueries = std::ceil(pi / 4 * std::sqrt(N));
  std::cout << "run " << nqueries << " query" << std::endl;
  for (idx i = 0; i < nqueries; ++i) {
    // apply oracle first, no aliasing
    psi(marked) = -psi(marked);
    // then the deiffusion operator
    psi = (G * psi).eval();
  }

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
