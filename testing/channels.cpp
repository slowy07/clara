#include <complex>
#include <iostream>
#include <valarray>
#include <vector>

#include "../include/clara.h"

using namespace clara;

int main() {
  idx nk = 5;
  idx D = 3;

  std::cout << "generating random channel with " << nk << " kraus operators on a " << D
            << " dimensional space\n";
  std::vector<cmat> Ks = randkraus(nk, D);

  // random input state
  cmat rho_in = randrho(D);
  // output state
  cmat rho_out = clara::apply(rho_in, Ks);

  std::cout << ">> computing choi matrix ... \n";
  cmat choim = kraus2choi(Ks);
  std::cout << ">> choi matrix:\n" << disp(choim) << std::endl;

  std::cout << ">> the eigenvalues of the choi matrix are: \n"
            << disp(transpose(hevals(choim))) << "\n";

  std::cout << "their sum is: " << sum(hevals(choim)) << "\n";

  std::vector<cmat> Kperps = choi2kraus(choim);
  std::cout << "kraus rank of channel : " << Kperps.size() << "\n";

  cmat rho_out1 = clara::apply(rho_in, Kperps);
  // verification norm difference
  std::cout << "norm difference ouput states: " << norm(rho_out1 - rho_out) << "\n";

  std::cout << "superoperator matrix:\n";
  cmat smat = kraus2super(Ks);
  std::cout << disp(smat) << "\n";

  std::cout << "the eigenvalues of the superoperator matrix :\n";
  dyn_col_vect<cplx> evalsupop = evals(smat);
  std::cout << disp(transpose(evalsupop)) << "\n";

  std::cout << "their absolute values are: \n";
  for (idx i = 0; i < (idx)evalsupop.size(); ++i)
    std::cout << std::abs(evalsupop(i)) << " ";

  // verification
  std::cout << "\n norm difference for the superoperator action: ";
  cmat rho_out2 = transpose(reshape(smat * reshape(transpose(rho_in), D * D, 1), D, D));
  std::cout << norm(rho_out - rho_out2) << "\n";
}
