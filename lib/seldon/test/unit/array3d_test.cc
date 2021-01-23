#include "SeldonLib.hxx"

using namespace Seldon;

namespace std
{
  inline bool isnan(const complex<double>& x)
  {
    if (isnan(real(x)))
      return true;
    
    if (isnan(imag(x)))
      return true;
    
    return false;
  }
}

Real_wp threshold;

template<class T>
void FillRand(Vector<T> & x)
{
  x.FillRand();
  Mlt(1e-9, x);
}

template<class T>
void FillRand(Vector<complex<T> > & x)
{
  for (int i = 0; i < x.GetM(); i++)
    x(i) = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T>
void FillRand(Array3D<T> & A)
{
  A.FillRand();
  Mlt(1e-9, A);
}

template<class T>
void FillRand(Array3D<complex<T> > & A)
{
  for (int i = 0; i < A.GetLength1(); i++)
    for (int j = 0; j < A.GetLength2(); j++)
      for (int k = 0; k < A.GetLength3(); k++)
        A(i, j, k) = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T>
void GetRand(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T> & x)
{
  int j = rand()%3;
  if (j == 0)
    x = complex<T>(0, rand())/T(RAND_MAX);
  else if (j == 1)
    x = complex<T>(rand(), 0)/T(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T>
bool EqualArray(const Array3D<T>& A, const Array3D<T>& C)
{
  if ( (A.GetLength1() != C.GetLength1())
       || (A.GetLength2() != C.GetLength2())
       || (A.GetLength3() != C.GetLength3()) )
    {
      return false;
    }
  
  for (int i = 0; i < A.GetLength1(); i++)
    for (int j = 0; j < A.GetLength2(); j++)
      for (int k = 0; k < A.GetLength3(); k++)
        if (isnan(A(i, j, k)) || isnan(C(i, j, k))
            || (abs(A(i, j, k) - C(i, j, k)) > threshold))
          return false;
  
  return true;
}

template<class T>
void CheckArray3D(Array3D<T>& A)
{
  int m = 4, n = 5, p  =3;
  A.Reallocate(m, n, p);
  
  if ( (A.GetLength1() != m) || (A.GetLength2() != n) || (A.GetLength3() != p)
       || (A.GetSize() != m*n*p) || (A.GetDataSize() != m*n*p) )
    {
      cout << "GetLength1, GetLength2, GetLength3, GetSize incorrect" << endl;
      abort();
    }
  
  Array3D<T> B(m, n, p);
  if ( (B.GetLength1() != m) || (B.GetLength2() != n) || (B.GetLength3() != p) )
    {
      cout << "Constructor incorrect" << endl;
      abort();
    }
  
  FillRand(A);
  Array3D<T> C(A);
  if (!EqualArray(A, C))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  T zero, x;
  SetComplexZero(zero);
  GetRand(x);
  C.Fill(x);
  for (int i = 0; i < A.GetLength1(); i++)
    for (int j = 0; j < A.GetLength2(); j++)
      for (int k = 0; k < A.GetLength3(); k++)
	if ((abs(C(i, j, k)-x) > threshold) || isnan(C(i, j, k)))
          {
            cout << "Fill incorrect " << endl;
            abort();
          }
  
  C.Fill(zero);
  C = A;
  if (!EqualArray(A, C))
    {
      cout << "Operator = incorrect" << endl;
      abort();
    }
  
  Vector<T> y(m*n*p);
  FillRand(y);
  for (int i = 0; i < A.GetLength1(); i++)
    for (int j = 0; j < A.GetLength2(); j++)
      for (int k = 0; k < A.GetLength3(); k++)
        A(i, j, k) = y(n*p*i + p*j + k);
  
  for (int i = 0; i < A.GetLength1(); i++)
    for (int j = 0; j < A.GetLength2(); j++)
      for (int k = 0; k < A.GetLength3(); k++)
        if ((A(i, j, k) != y(n*p*i + p*j + k)) || isnan(A(i, j, k)))
          {
            cout << "Operator () incorrect" << endl;
            abort();
          }
  
  A.Write("toto.dat");
  C.Fill(zero);
  C.Read("toto.dat");
  
  if (!EqualArray(A, C))
    {
      cout << "Read/Write incorrect" << endl;
      abort();
    }
  
  C.Clear();
  if ( (C.GetLength1() != 0) || (C.GetLength2() != 0) || (C.GetLength3() != 0)
       || (C.GetSize() != 0) || (C.GetDataSize() != 0) )
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
}

int main(int argc, char** argv)
{
  threshold = 1e-12;

  srand(0);
  
  {
    Array3D<Real_wp> A;
    CheckArray3D(A);
  }

  {
    Array3D<Complex_wp> A;
    CheckArray3D(A);
  }
  
  std::remove("toto.dat");
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
