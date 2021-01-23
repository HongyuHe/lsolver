#include "SeldonLib.hxx"
#include <assert.h>

using namespace Seldon;

int main()
{
  Vector<double> v(2);
  Matrix<double> A(4, 4);
  
  v(0) = 1.0;

  A(2, 3) = 2.0;

  TinyVector<double, 2> x(0.3, 0.5);
  TinyVector<double, 2> y(0.7, 0.9), z;

  z = x + y;

  DISP(z);
  
  return 0;
}
