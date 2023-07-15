#include <iostream>
#include "../include/clara.h"
#include "../include/experimental/experimental_test.h"

using namespace clara;

int main() {
  std::cout << "testing testing \n";
  PRINTLN("testing dbg message");
  ERRORLN("testing error debug message");
}
