#include <iostream>
#include "../include/clara.h"

using namespace clara;

int main() {
  std::cout << "testing quantum state |0> state: \n";
  std::cout <<  disp(st.z0) << std::endl;
}
