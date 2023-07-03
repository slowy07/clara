#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include "../include/clara.h"
#include "../include/internal.h"

using namespace std;
using namespace clara;
using namespace clara::types;

int main() {
  _init();
  // Eigen::MatrixXd _a = Eigen::MatrixXd::Random(4, 4);
  // save(_a, "../clara/_a");
  //
  // Eigen::MatrixXd a = load<Eigen::MatrixXd>("../clara/_a");
  // std::vector<size_t> subsys = {2, 2};
  // std::vector<size_t> perm = {1, 0};
  //
  // cout << "error in norm difference load/save" << norm(_a - a) << endl;
  // 
  // disp(ptrace2(a, {2, 2}));
  // cout<<endl<<endl;
  // disp(ptrace(a, {0}, {2, 2}));
  // cout <<endl<<endl;
  
  imat kt(3, 1);
  kt << 0, 0, 1;

  imat bt(1, 3);
  bt << 0, 1, 0;
  
  disp(kron(kt, bt).template cast<double>());
  cout <<endl<<endl;

  disp(bt * kt);
  cout<<endl<<endl;
  
  size_t dim = 10;
  cout << "generate random unitary" << endl;
  cmat u = rand_unitary(dim);
  cout << "done generating random unitary";
  disp(u);
  cout<<endl;
  double normdiff = norm((cmat) (u * adjoint(u) - cmat::Identity(dim, dim)));
  cout << "norm difference" << normdiff << endl;
  disp(normdiff, std::cout, 18);
  if (normdiff > std::numeric_limits<double>::epsilon())
    cout << "yes";
  else
    cout << "no";
  cout << endl << endl;
  cout <<"the eigen values of are: "<<endl;
  disp(transpose(evals(u)));
  cout << endl << endl;
  
  cout << "the absolute values of the eigen of are: " << endl;
  disp(abs(transpose(evals(u))));
  cout<<endl<<endl;

  imat im(2, 2);
  im << 1, 2, 3, 4;
  disp(funm(adjoint(im) *im, [](const cplx x) -> cplx {return std::sqrt(x);}));
}
