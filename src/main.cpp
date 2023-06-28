#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include "../include/clara.h"

using namespace std;
using namespace clara;
using namespace clara::types;

int main() {
  _init();

  stat::UniformRealDistribution a(-2, 2);
  cout << endl << a.sample() << endl;

  stat::NormalDistribution b;
  cout << endl << b.sample() << endl;
  
  cmat A = cmat::Random(2, 2);
  cout << endl;
  disp(A);

  cout << endl;
  int n = 41;
  cout << "the " << n << " root of unity is: " << ct::omega(n) << endl;
}
