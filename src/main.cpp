#include <iostream>

#include "../include/clara.h"

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
  cout << endl << endl << "Clara testing finish..." << endl;
  
}
