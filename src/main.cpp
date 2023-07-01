#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include "../include/clara.h"

using namespace std;
using namespace clara;
using namespace clara::types;

int main(int /*argc*/, char **/*argv*/) {
  _init();

  size_t n = 3;
  size_t nout = n - 1;
  
  size_t dim = pow(2, n);

  ivect dims(n, 2);

  ivect subsys(nout);
	for (size_t i = 0; i < nout; i++) {
		subsys(i)=i;
  }

  ivect perm(n - 1);
  for(size_t i = 1; i < n; i++)
    perm(i - 1) = i;
  perm(n - 1) = 0;


  cmat A = randn(dim);

  cvect v(3);
  v<<ct::ii,2,3;
  disp(v);
  cout<<endl<<endl;
  cout<<v<<endl<<endl;

  ivect a(2);
  a<<1,2;

  cmat c = randn(3);
  disp(mat_pow(c, 1));
  cout<<endl<<endl;
  disp(static_cast<cmat>(transpose(v) * v));
  
  // cplx vv=trace(c);
}
