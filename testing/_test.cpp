#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <ostream>
#include <random>
#include <vector>

#include "../include/clara.h"
#include "../include/experimental/experimental_test.h"

int main() {
  using namespace clara;
  std::cout << "testing testing \n";

  const idx bits = 70;
  experimental::Bit_circuit b{bits};

  const idx trials = 20;

  // number of trials
  b.rand();
  auto c = b;

  std::random_device rd;
  std::mt19937 gen{rd()};
  std::vector<std::vector<idx>> indices(trials);

  for (idx i = 0; i < trials; ++i) {
    std::vector<idx> v(bits);
    std::iota(v.begin(), v.end(), 0);
    std::shuffle(v.begin(), v.end(), gen);
    std::vector<idx> tof(v.data(), v.data() + 3);
    indices[i] = tof;
  }

  for (idx i = 0; i < trials; ++i) {
    std::cout << "first: ";
    for (auto&& elem : indices[i])
      std::cout << elem << " ";
    std::cout << std::endl;
    b.TOF(indices[i]);
  }

  for (idx i = trials; i-- > 0;) {
    std::cout << "second: ";
    for (auto&& elem : indices[i])
      std::cout << elem << " ";
    std::cout << std::endl;
    b.TOF(indices[i]);
  }

  std::cout << "initial: " << b << std::endl;
  std::cout << "final: " << c << std::endl;
  std::cout << "hamming weight: " << b.count() << std::endl;

  std::cout << b.gate_count.NOT << " " << b.gate_count.X << " " << b.gate_count.TOF << std::endl;
  std::cout << (b == c) << std::endl;
  std::cout << (b != c) << std::endl;

  experimental::Dynamic_bitset bb(9);
  bb.set(1).set(3).set(8);
  std::cout << bb << std::endl;

  std::cout << "info: " << std::endl;
  std::cout << bb.to_string('o', 'X') << std::endl;

  experimental::Dynamic_bitset vlad(20);
  std::cout << vlad << std::endl;

  std::vector<unsigned int> vv(20);
  for (auto& elem : vv) {
    std::cout << elem;
  }
  std::cout << std::endl;

  ket x = (10_ket + 01_ket) / std::sqrt(2);
  std::cout << disp(x) << std::endl;

  bra y = (10_bra + 01_bra) / std::sqrt(2);
  std::cout << disp(x) << std::endl;

  cmat z = 110_prj;
  std::cout << disp(z) << std::endl;
}
