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

  vector<size_t> dims={2, 3, 4, 5, 6, 7};
  vector<size_t> subsys={1, 2, 4};
  size_t dim = 5040;
  
  size_t cnt = 0;
  cmat A(dim, dim);
  for (size_t i = 0; i<dim; i++)
    for (size_t j = 0; j<dim; j++)
      A(i, j) = cnt++;
  disp(ptrace(A, dims, subsys));
  cout<<endl<<norm(ptrace(A, dims, subsys)) << endl;
}
