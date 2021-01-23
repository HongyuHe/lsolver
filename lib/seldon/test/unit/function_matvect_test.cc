#include "SeldonLib.hxx"

using namespace Seldon;

typedef Vector<Real_wp> VectReal_wp;

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

template<bool f>
class GhostIf
{
};

Real_wp threshold;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz, GhostIf<true>& sparse_form,
                          GhostIf<false>& triang_form, bool low)
{
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  A.Reallocate(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRandNumber(x);
      A.Set(i, j, x);
    }
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz, GhostIf<false>& sparse_form,
                          GhostIf<false>& triang_form, bool low)
{
  A.Reallocate(m, n);
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
        GetRandNumber(x);
        A.Set(i, j, x);
      }
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz, GhostIf<false>& sparse_form,
                          GhostIf<true>& triang_form, bool low)
{
  A.Reallocate(m, n);
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  if (low)
    for (int i = 0; i < m; i++)
      for (int j = 0; j <= i; j++)
        {
          GetRandNumber(x);
          A.Set(i, j, x);
        }
  else
    for (int i = 0; i < m; i++)
      for (int j = i; j < n; j++)
        {
          GetRandNumber(x);
          A.Set(i, j, x);
        }
}


template<class T, class Prop, class Storage, class Allocator, class T1, class T2>
void MltTest(const Matrix<T, Prop, Storage, Allocator> & A,
             const Vector<T1> & X, Vector<T2>& Y)
{
  T2 zero, val;
  SetComplexZero(zero);
  Y.Fill(zero);
  for (int i = 0; i < Y.GetM(); i++)
    {
      val = zero;
      for (int j = 0; j < X.GetM(); j++)
        val += A(i, j)*X(j);
      
      Y(i) = val;
    }
}

template<class T, class Prop, class Storage, class Allocator, class T1, class T2>
void MltTransTest(const Matrix<T, Prop, Storage, Allocator> & A,
                  const Vector<T1> & X, Vector<T2>& Y)
{
  T2 zero, val;
  SetComplexZero(zero);
  Y.Fill(zero);
  for (int i = 0; i < Y.GetM(); i++)
    {
      val = zero;
      for (int j = 0; j < X.GetM(); j++)
        val += A(j, i)*X(j);
      
      Y(i) = val;
    }
}

template<class T, class Prop, class Storage, class Allocator, class T1, class T2>
void MltConjTransTest(const Matrix<T, Prop, Storage, Allocator> & A,
                      const Vector<T1> & X, Vector<T2>& Y)
{
  T2 zero, val;
  SetComplexZero(zero);
  Y.Fill(zero);
  for (int i = 0; i < Y.GetM(); i++)
    {
      val = zero;
      for (int j = 0; j < X.GetM(); j++)
        val += conjugate(A(j, i))*X(j);
      
      Y(i) = val;
    }
}

template<class T1, class Prop, class Storage, class Allocator,
	 class T2, class T3>
void SorTest(const Matrix<T1, Prop, Storage, Allocator>& A,
	     Vector<T2>& X, const Vector<T2>& B,
	     const T3& omega, int iter, int type)
{
  // forward sweep
  T2 val, zero, one;
  SetComplexZero(zero);
  SetComplexOne(one);
  int n = A.GetM();
  if (type%2 == 0)
    for (int num = 0; num < iter; num++)
      for (int i = 0; i < n; i++)
	{
	  val = B(i);
	  for (int j = 0; j < n; j++)
	    if (i != j)
	      val -= A(i, j)*X(j);
	  
	  X(i) = (one - omega) * X(i) + omega * val / A(i, i);
	}
  
  if (type%3 == 0)
    for (int num = 0; num < iter; num++)
      for (int i = n-1; i >= 0; i--)
	{
	  val = B(i);
	  for (int j = 0; j < n; j++)
	    if (i != j)
	      val -= A(i, j)*X(j);
	  
	  X(i) = (one - omega) * X(i) + omega * val / A(i, i);
	}
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y, Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > eps) || isnan(abs(x(i)-y(i))))
      {
        DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void CheckGaussLU(Matrix<T, Prop, Storage, Allocator>& A, 
                  Matrix<T, Prop, Storage, Allocator>& B, 
                  GhostIf<true>& herm_form, Vector<T>& x, Vector<T>& y, Vector<T>& z)
{
}

template<class T, class Prop, class Storage, class Allocator>
void CheckGaussLU(Matrix<T, Prop, Storage, Allocator>& A, 
                  Matrix<T, Prop, Storage, Allocator>& B, 
                  GhostIf<false>& herm_form, Vector<T>& x, Vector<T>& y, Vector<T>& z)
{
}

template<class T, class Storage, class Allocator>
void CheckGaussLU(Matrix<T, General, Storage, Allocator>& A, 
                  Matrix<T, General, Storage, Allocator>& B, 
                  GhostIf<false>& herm_form, Vector<T>& x, Vector<T>& y, Vector<T>& z)
{
      Gauss(A, y);
      
      if (!EqualVector(x, y, 1000*threshold))
        {
          cout << "Gauss incorrect" << endl;
          abort();
        }
      
      A = B;
      y = z;      
      GetLU(A);
      SolveLU(A, y);
      if (!EqualVector(x, y, 1000*threshold))
        {
          cout << "GetLU/SolveLU incorrect" << endl;
          abort();
        }
}

template<class T, class Prop, class Storage, class Allocator, bool herm>
void CheckSolveMatrix(Matrix<T, Prop, Storage, Allocator>& A, 
                      GhostIf<false>& sparse_form, GhostIf<herm>& herm_form,
                      GhostIf<false>& triang_form)
{
  int m = 72, n = m, nnz = 300;
  
  Vector<T> x, y, z;
  T zero;
  SetComplexZero(zero);
  
  GenerateRandomMatrix(A, m, n, nnz, sparse_form, triang_form, true);
  GenerateRandomVector(x, n);
  GenerateRandomVector(y, m);

  if (herm)
    {
      T val;
      for (int i = 0; i < A.GetM(); i++)
        {
          SetComplexReal(realpart(A(i, i)), val);
          A.Set(i, i, val);
        }
    }
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  Mlt(A, x, y);
  z = y;

  CheckGaussLU(A, B, herm_form, x, y, z);
  
  A = B;
  for (int i = 0; i < A.GetM(); i++)
    {
      T val;
      SetComplexReal(3, val);
      A.Get(i, i) += val;
    }
  
  Mlt(A, x, z);
  
  x.Fill(zero);
  typedef typename ClassComplexType<T>::Treal Treal;
  Treal omega = Treal(1)/2;
  SorTest(A, x, z, omega, 3, 0);

  y.Fill(zero);
  SOR(A, y, z, omega, 3, 0);

  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }

  x.Fill(zero);
  SorTest(A, x, z, omega, 3, 2);

  y.Fill(zero);
  SOR(A, y, z, omega, 3, 2);

  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }

  x.Fill(zero);
  SorTest(A, x, z, omega, 3, 3);

  y.Fill(zero);
  SOR(A, y, z, omega, 3, 3);

  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }

  omega = Treal(1);
  x.Fill(zero);
  SorTest(A, x, z, omega, 3, 0);

  y.Fill(zero);
  GaussSeidel(A, y, z, 3, 0);

  if (!EqualVector(x, y, 100*threshold))
    {
      cout << "GaussSeidel incorrect" << endl;
      abort();
    }

  x.Fill(zero);
  SorTest(A, x, z, omega, 3, 2);

  y.Fill(zero);
  GaussSeidel(A, y, z, 3, 2);

  if (!EqualVector(x, y, 100*threshold))
    {
      cout << "GaussSeidel incorrect" << endl;
      abort();
    }

  x.Fill(zero);
  SorTest(A, x, z, omega, 3, 3);

  y.Fill(zero);
  GaussSeidel(A, y, z, 3, 3);

  if (!EqualVector(x, y, 100*threshold))
    {
      cout << "GaussSeidel incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class Storage, class Allocator, bool herm>
void CheckSolveMatrix(Matrix<T, Prop, Storage, Allocator>& A, 
                      GhostIf<false>& sparse_form, GhostIf<herm>& herm_form,
                      GhostIf<true>& triang_form)
{
  for (int i = 0; i < A.GetM(); i++)
    {
      T sum;
      SetComplexZero(sum);
      for (int j = 0; j < A.GetM(); j++)
	if (j != i)
	  sum += abs(A(i, j));
      
      A.Set(i, i, abs(A(i, i)) + sum);
    }
  
  // testing Mlt of triangular systems (Blas interface)
  int n = A.GetM();
  Vector<T> X, Y, X0;
  GenerateRandomVector(X, n);
  
  X0 = X;
  Y = X;
  MltTest(A, X, Y);
  Mlt(A, X);
  
  if (!EqualVector(X, Y))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Solve(A, X);

  if (!EqualVector(X, X0, 1e3*threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  X = X0;
  Y = X;
  MltTest(A, X, Y);
  Mlt(SeldonNoTrans, SeldonNonUnit, A, X);
  
  if (!EqualVector(X, Y))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }

  Solve(SeldonNoTrans, SeldonNonUnit, A, X);

  if (!EqualVector(X, X0, 1e3*threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
    

  X = X0;
  Y = X;
  MltTransTest(A, X, Y);
  Mlt(SeldonTrans, SeldonNonUnit, A, X);
  
  if (!EqualVector(X, Y))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Solve(SeldonTrans, SeldonNonUnit, A, X);

  if (!EqualVector(X, X0, 1000*threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  if (IsComplexMatrix(A))
    {
      X = X0;
      Y = X;
      MltConjTransTest(A, X, Y);
      Mlt(SeldonConjTrans, SeldonNonUnit, A, X);
      
      if (!EqualVector(X, Y))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}

      Solve(SeldonConjTrans, SeldonNonUnit, A, X);
      
      if (!EqualVector(X, X0, 1e3*threshold))
	{
	  cout << "Solve incorrect" << endl;
	  abort();
	}
    }

  T one;
  SetComplexOne(one);
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  for (int i = 0; i < n; i++)
    B.Set(i, i, one);

  X = X0;
  Y = X;
  MltTest(B, X, Y);
  Mlt(SeldonNoTrans, SeldonUnit, A, X);
  
  if (!EqualVector(X, Y))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }

  Solve(SeldonNoTrans, SeldonUnit, A, X);
  
  if (!EqualVector(X, X0, 1000*threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  X = X0;
  Y = X;
  MltTransTest(B, X, Y);
  Mlt(SeldonTrans, SeldonUnit, A, X);
  
  if (!EqualVector(X, Y))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }

  Solve(SeldonTrans, SeldonUnit, A, X);

  if (!EqualVector(X, X0, 1e3*threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  if (IsComplexMatrix(A))
    {
      X = X0;
      Y = X;
      MltConjTransTest(B, X, Y);
      Mlt(SeldonConjTrans, SeldonUnit, A, X);
      
      if (!EqualVector(X, Y))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}

      Solve(SeldonConjTrans, SeldonUnit, A, X);
      
      if (!EqualVector(X, X0, 1000*threshold))
	{
	  cout << "Solve incorrect" << endl;
	  abort();
	}
    }
}

template<class T, class Prop, class Storage, class Allocator>
void CheckSolveMatrix(Matrix<T, Prop, Storage, Allocator>& A, 
                      GhostIf<true>& sparse_form, GhostIf<false>& herm_form,
                      GhostIf<false>& triang_form)
{
  // not tested in this file
}

template<class T, class Prop, class Storage, class Allocator, bool triang>
void CheckRealMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                     GhostIf<triang>& triang_form, bool low)
{
  int m = 70, n = 36, nnz = 200;
  if (IsSymmetricMatrix(A) || triang)
    {
      m = 53;
      n = 53;
    }
  
  T one, zero;
  SetComplexZero(zero);
  SetComplexOne(one);
  
  GhostIf<Storage::Sparse> sparse_form;
  GenerateRandomMatrix(A, m, n, nnz, sparse_form, triang_form, low);
  
  VectReal_wp x, y, z;
  GenerateRandomVector(x, n);
  GenerateRandomVector(y, m);
  GenerateRandomVector(z, m);
  
  // testing y = A*x
  Mlt(A, x, y);
  MltTest(A, x, z);
  
  if (!EqualVector(y, z))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  y.Fill(zero);
  Mlt(SeldonNoTrans, A, x, y);
  
  if (!EqualVector(y, z))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  y.Fill(zero);
  Mlt(one, A, x, y);
  
  if (!EqualVector(y, z))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  y.Fill(zero);
  T alpha(2.3);
  Mlt(alpha, A, x, y);
  
  Mlt(alpha, z);
  if (!EqualVector(y, z))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Mlt(T(1.0)/alpha, z);

  // testing y = A*x when x and y are complex  
  Vector<complex<T> > xc, yc, zc;
  complex<T> czero(0, 0), cone(1, 0);
  
  GenerateRandomVector(xc, n);
  GenerateRandomVector(yc, m);
  GenerateRandomVector(zc, m);
  
  Mlt(A, xc, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      DISP(xc); DISP(yc);
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  yc.Fill(czero);
  Mlt(SeldonNoTrans, A, xc, yc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  yc.Fill(czero);
  Mlt(cone, A, xc, yc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  yc.Fill(czero);
  complex<T> alphac(2.3, 0.7);
  Mlt(alphac, A, xc, yc);
  
  Mlt(alphac, zc);
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Mlt(cone/alphac, zc);
  yc.Fill(czero);
  
  // testing y = A^T*x  
  GenerateRandomVector(x, m);
  GenerateRandomVector(y, n);
  GenerateRandomVector(z, n);
  GenerateRandomVector(zc, n);

  Mlt(SeldonTrans, A, x, y);
  MltTransTest(A, x, z);
  if (!EqualVector(y, z))
    {
      cout << "Mlt(SeldonTrans) incorrect" << endl;
      abort();
    }

  GenerateRandomVector(xc, m);
  GenerateRandomVector(yc, n);
  
  Mlt(SeldonTrans, A, xc, yc);
  MltTransTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt(SeldonTrans) incorrect" << endl;
      abort();
    }
  
  // testing MltAdd
  GenerateRandomVector(x, n);
  GenerateRandomVector(y, m);
  GenerateRandomVector(z, m);
  
  MltAdd(one, A, x, zero, y);
  MltTest(A, x, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-one, A, x, zero, y);
  MltTest(A, x, z);
  Mlt(-one, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, A, x, zero, y);
  MltTest(A, x, z);
  Mlt(alpha, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Vector<T> w;
  GenerateRandomVector(w, m);
  T beta(0.83);
  
  Copy(w, y);
  MltAdd(alpha, A, x, beta, y);
  
  MltTest(A, x, z);
  Mlt(alpha, z);
  Add(beta, w, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(xc, n);
  GenerateRandomVector(yc, m);
  GenerateRandomVector(zc, m);
  
  MltAdd(cone, A, xc, czero, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(one, A, xc, zero, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-one, A, xc, zero, yc);
  MltTest(A, xc, zc);
  Mlt(-one, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alphac, A, xc, czero, yc);
  MltTest(A, xc, zc);
  Mlt(alphac, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, A, xc, zero, yc);
  MltTest(A, xc, zc);
  Mlt(alpha, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Vector<complex<T> > wc;
  GenerateRandomVector(wc, m);
  complex<T> betac(1.17, 2.8);
  
  Copy(wc, yc);
  MltAdd(alphac, A, xc, betac, yc);
  
  MltTest(A, xc, zc);
  Mlt(alphac, zc);
  Add(betac, wc, zc);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Copy(wc, yc);
  MltAdd(alpha, A, xc, beta, yc);
  
  MltTest(A, xc, zc);
  Mlt(alpha, zc);
  Add(beta, wc, zc);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  
  
  // testing MltAdd with SeldonNoTrans
  GenerateRandomVector(x, n);
  GenerateRandomVector(y, m);
  GenerateRandomVector(z, m);
  
  MltAdd(one, SeldonNoTrans, A, x, zero, y);
  MltTest(A, x, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-one, SeldonNoTrans, A, x, zero, y);
  MltTest(A, x, z);
  Mlt(-one, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, SeldonNoTrans, A, x, zero, y);
  MltTest(A, x, z);
  Mlt(alpha, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(w, m);
  
  Copy(w, y);
  MltAdd(alpha, SeldonNoTrans, A, x, beta, y);
  
  MltTest(A, x, z);
  Mlt(alpha, z);
  Add(beta, w, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(xc, n);
  GenerateRandomVector(yc, m);
  GenerateRandomVector(zc, m);
  
  MltAdd(cone, SeldonNoTrans, A, xc, czero, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(one, SeldonNoTrans, A, xc, zero, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-one, SeldonNoTrans, A, xc, zero, yc);
  MltTest(A, xc, zc);
  Mlt(-one, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alphac, SeldonNoTrans, A, xc, czero, yc);
  MltTest(A, xc, zc);
  Mlt(alphac, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, SeldonNoTrans, A, xc, zero, yc);
  MltTest(A, xc, zc);
  Mlt(alpha, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(wc, m);
  
  Copy(wc, yc);
  MltAdd(alphac, SeldonNoTrans, A, xc, betac, yc);
  
  MltTest(A, xc, zc);
  Mlt(alphac, zc);
  Add(betac, wc, zc);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Copy(wc, yc);
  MltAdd(alpha, SeldonNoTrans, A, xc, beta, yc);
  
  MltTest(A, xc, zc);
  Mlt(alpha, zc);
  Add(beta, wc, zc);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  // testing MltAdd with transpose
  GenerateRandomVector(x, m);
  GenerateRandomVector(y, n);
  GenerateRandomVector(z, n);
  
  MltAdd(one, SeldonTrans, A, x, zero, y);
  MltTransTest(A, x, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-one, SeldonTrans, A, x, zero, y);
  MltTransTest(A, x, z);
  Mlt(-one, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, SeldonTrans, A, x, zero, y);
  MltTransTest(A, x, z);
  Mlt(alpha, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(w, n);
  
  Copy(w, y);
  MltAdd(alpha, SeldonTrans, A, x, beta, y);
  
  MltTransTest(A, x, z);
  Mlt(alpha, z);
  Add(beta, w, z);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(xc, m);
  GenerateRandomVector(yc, n);
  GenerateRandomVector(zc, n);
  
  MltAdd(cone, SeldonTrans, A, xc, czero, yc);
  MltTransTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(one, SeldonTrans, A, xc, zero, yc);
  MltTransTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-one, SeldonTrans, A, xc, zero, yc);
  MltTransTest(A, xc, zc);
  Mlt(-one, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alphac, SeldonTrans, A, xc, czero, yc);
  MltTransTest(A, xc, zc);
  Mlt(alphac, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, SeldonTrans, A, xc, zero, yc);
  MltTransTest(A, xc, zc);
  Mlt(alpha, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(wc, n);
  
  Copy(wc, yc);
  MltAdd(alphac, SeldonTrans, A, xc, betac, yc);
  
  MltTransTest(A, xc, zc);
  Mlt(alphac, zc);
  Add(betac, wc, zc);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Copy(wc, yc);
  MltAdd(alpha, SeldonTrans, A, xc, beta, yc);
  
  MltTransTest(A, xc, zc);
  Mlt(alpha, zc);
  Add(beta, wc, zc);
  
  if (!EqualVector(y, z))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  
  
  GhostIf<false> herm_form;
  CheckSolveMatrix(A, sparse_form, herm_form, triang_form);
}

template<class Tcplx, class Prop, class Storage, class Allocator, bool herm, bool triang>
void CheckComplexMatrix(Matrix<Tcplx, Prop, Storage, Allocator>& A,
                        GhostIf<herm>& herm_form, GhostIf<triang>& triang_form,
                        bool low)
{
  typedef typename ClassComplexType<Tcplx>::Treal T;
  
  int m = 70, n = 36, nnz = 200;
  if ( IsSymmetricMatrix(A) || herm || triang)
    {
      m = 55;
      n = 55;
    }
  
  T one, zero;
  SetComplexZero(zero);
  SetComplexOne(one);

  complex<T> cone, czero;
  SetComplexZero(czero);
  SetComplexOne(cone);
  
  GhostIf<Storage::Sparse> sparse_form;
  GenerateRandomMatrix(A, m, n, nnz, sparse_form, triang_form, low);
  
  if (herm)
    for (int i = 0; i < m; i++)
      A.Set(i, i, complex<T>(real(A(i, i)), 0));
  
  // testing y = A*x when x and y are complex  
  Vector<complex<T> > xc, yc, zc;
  
  GenerateRandomVector(xc, n);
  GenerateRandomVector(yc, m);
  GenerateRandomVector(zc, m);

  Mlt(A, xc, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  yc.Fill(czero);
  Mlt(SeldonNoTrans, A, xc, yc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  yc.Fill(czero);
  Mlt(cone, A, xc, yc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  complex<T> alphac(2.3, 0.7);
  T alpha(2.3);
  
  yc.Fill(czero);  
  Mlt(alphac, A, xc, yc);
  
  Mlt(alphac, zc);
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }

  Mlt(cone/alphac, zc);

  // testing y = A^T*x  
  GenerateRandomVector(xc, m);
  GenerateRandomVector(yc, n);
  GenerateRandomVector(zc, n);
  
  Mlt(SeldonTrans, A, xc, yc);
  MltTransTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "Mlt(SeldonTrans) incorrect" << endl;
      abort();
    }
  
  // testing MltAdd
  GenerateRandomVector(xc, n);
  GenerateRandomVector(yc, m);
  GenerateRandomVector(zc, m);
  
  MltAdd(cone, A, xc, czero, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(one, A, xc, zero, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-cone, A, xc, czero, yc);
  MltTest(A, xc, zc);
  Mlt(-cone, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alphac, A, xc, czero, yc);
  MltTest(A, xc, zc);
  Mlt(alphac, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Vector<complex<T> > wc;
  GenerateRandomVector(wc, m);
  complex<T> betac(1.17, 2.8);
  T beta(3.2);
  
  Copy(wc, yc);
  MltAdd(alphac, A, xc, betac, yc);
  
  MltTest(A, xc, zc);
  Mlt(alphac, zc);
  Add(betac, wc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Copy(wc, yc);
  MltAdd(alpha, A, xc, beta, yc);
  
  MltTest(A, xc, zc);
  Mlt(alpha, zc);
  Add(beta, wc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  
  
  // testing MltAdd with SeldonNoTrans
  GenerateRandomVector(xc, n);
  GenerateRandomVector(yc, m);
  GenerateRandomVector(zc, m);
  
  MltAdd(cone, SeldonNoTrans, A, xc, czero, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(one, SeldonNoTrans, A, xc, zero, yc);
  MltTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-cone, SeldonNoTrans, A, xc, czero, yc);
  MltTest(A, xc, zc);
  Mlt(-one, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alphac, SeldonNoTrans, A, xc, czero, yc);
  MltTest(A, xc, zc);
  Mlt(alphac, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, SeldonNoTrans, A, xc, zero, yc);
  MltTest(A, xc, zc);
  Mlt(alpha, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(wc, m);
  
  Copy(wc, yc);
  MltAdd(alphac, SeldonNoTrans, A, xc, betac, yc);
  
  MltTest(A, xc, zc);
  Mlt(alphac, zc);
  Add(betac, wc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Copy(wc, yc);
  MltAdd(alpha, SeldonNoTrans, A, xc, beta, yc);
  
  MltTest(A, xc, zc);
  Mlt(alpha, zc);
  Add(beta, wc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  // testing MltAdd with transpose
  GenerateRandomVector(xc, m);
  GenerateRandomVector(yc, n);
  GenerateRandomVector(zc, n);
  
  MltAdd(cone, SeldonTrans, A, xc, czero, yc);
  MltTransTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(one, SeldonTrans, A, xc, zero, yc);
  MltTransTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-cone, SeldonTrans, A, xc, czero, yc);
  MltTransTest(A, xc, zc);
  Mlt(-one, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alphac, SeldonTrans, A, xc, czero, yc);
  MltTransTest(A, xc, zc);
  Mlt(alphac, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, SeldonTrans, A, xc, zero, yc);
  MltTransTest(A, xc, zc);
  Mlt(alpha, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(wc, n);
  
  Copy(wc, yc);
  MltAdd(alphac, SeldonTrans, A, xc, betac, yc);
  
  MltTransTest(A, xc, zc);
  Mlt(alphac, zc);
  Add(betac, wc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Copy(wc, yc);
  MltAdd(alpha, SeldonTrans, A, xc, beta, yc);
  
  MltTransTest(A, xc, zc);
  Mlt(alpha, zc);
  Add(beta, wc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  
  
  // testing MltAdd with transpose
  GenerateRandomVector(xc, m);
  GenerateRandomVector(yc, n);
  GenerateRandomVector(zc, n);
  
  MltAdd(cone, SeldonConjTrans, A, xc, czero, yc);
  MltConjTransTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(one, SeldonConjTrans, A, xc, zero, yc);
  MltConjTransTest(A, xc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(-cone, SeldonConjTrans, A, xc, czero, yc);
  MltConjTransTest(A, xc, zc);
  Mlt(-one, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alphac, SeldonConjTrans, A, xc, czero, yc);
  MltConjTransTest(A, xc, zc);
  Mlt(alphac, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  MltAdd(alpha, SeldonConjTrans, A, xc, zero, yc);
  MltConjTransTest(A, xc, zc);
  Mlt(alpha, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  GenerateRandomVector(wc, n);
  
  Copy(wc, yc);
  MltAdd(alphac, SeldonConjTrans, A, xc, betac, yc);
  
  MltConjTransTest(A, xc, zc);
  Mlt(alphac, zc);
  Add(betac, wc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  

  Copy(wc, yc);
  MltAdd(alpha, SeldonConjTrans, A, xc, beta, yc);
  
  MltConjTransTest(A, xc, zc);
  Mlt(alpha, zc);
  Add(beta, wc, zc);
  
  if (!EqualVector(yc, zc))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }  
  
  CheckSolveMatrix(A, sparse_form, herm_form, triang_form);
}

int main(int argc, char** argv)
{
  threshold = 1e-12;
  
  //srand(time(NULL));

  GhostIf<false> triang;
  
  {
    Matrix<Real_wp, General, RowMajor> A;
    CheckRealMatrix(A, triang, false);
  }
  
  {
    Matrix<Real_wp, General, ColMajor> A;
    CheckRealMatrix(A, triang, false);
  }
  
  {
    Matrix<Real_wp, Symmetric, RowSymPacked> A;
    CheckRealMatrix(A, triang, false);
  }
  
#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSymPacked> A;
    CheckRealMatrix(A, triang, false);
  }
#endif
  
  {
    Matrix<Real_wp, Symmetric, RowSym> A;
    CheckRealMatrix(A, triang, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSym> A;
    CheckRealMatrix(A, triang, false);
  }
#endif
  
  {
    Matrix<Complex_wp, Hermitian, RowHermPacked> A;
    GhostIf<true> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }
  
#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Hermitian, ColHermPacked> A;
    GhostIf<true> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }
#endif

 {
    Matrix<Complex_wp, Hermitian, RowHerm> A;
    GhostIf<true> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
 {
    Matrix<Complex_wp, Hermitian, ColHerm> A;
    GhostIf<true> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }
#endif
  
  {
    Matrix<Real_wp, General, RowSparse> A;
    CheckRealMatrix(A, triang, false);
  }
  
  {
    Matrix<Complex_wp, General, RowSparse> A;
    GhostIf<false> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }
  
  {
    Matrix<Real_wp, General, ColSparse> A;
    CheckRealMatrix(A, triang, false);
  }

  {
    Matrix<Complex_wp, General, ColSparse> A;
    GhostIf<false> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }

  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckRealMatrix(A, triang, false);
  }

  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    GhostIf<false> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ArrayColSparse> A;
    CheckRealMatrix(A, triang, false);
  }

  {
    Matrix<Complex_wp, General, ArrayColSparse> A;
    GhostIf<false> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }
#endif

  {
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    CheckRealMatrix(A, triang, false);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymSparse> A;
    GhostIf<false> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSymSparse> A;
    CheckRealMatrix(A, triang, false);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSymSparse> A;
    GhostIf<false> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }
#endif
  
  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckRealMatrix(A, triang, false);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    GhostIf<false> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A;
    CheckRealMatrix(A, triang, false);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayColSymSparse> A;
    GhostIf<false> herm_form;
    CheckComplexMatrix(A, herm_form, triang, false);
  }
#endif
  
  GhostIf<true> triang2;
  GhostIf<false> non_herm;
  
  {
    Matrix<Real_wp, General, RowLoTriang> A;
    CheckRealMatrix(A, triang2, true);
  }

  {
    Matrix<Complex_wp, General, RowLoTriang> A;
    CheckComplexMatrix(A, non_herm, triang2, true);
  }

  {
    Matrix<Real_wp, General, RowLoTriangPacked> A;
    CheckRealMatrix(A, triang2, true);
  }

  {
    Matrix<Complex_wp, General, RowLoTriangPacked> A;
    CheckComplexMatrix(A, non_herm, triang2, true);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColLoTriang> A;
    CheckRealMatrix(A, triang2, true);
  }

  {
    Matrix<Complex_wp, General, ColLoTriang> A;
    CheckComplexMatrix(A, non_herm, triang2, true);
  }

  {
    Matrix<Real_wp, General, ColLoTriangPacked> A;
    CheckRealMatrix(A, triang2, true);
  }

  {
    Matrix<Complex_wp, General, ColLoTriangPacked> A;
    CheckComplexMatrix(A, non_herm, triang2, true);
  }
#endif

  {
    Matrix<Real_wp, General, RowUpTriang> A;
    CheckRealMatrix(A, triang2, false);
  }

  {
    Matrix<Complex_wp, General, RowUpTriang> A;
    CheckComplexMatrix(A, non_herm, triang2, false);
  }

  {
    Matrix<Real_wp, General, RowUpTriangPacked> A;
    CheckRealMatrix(A, triang2, false);
  }

  {
    Matrix<Complex_wp, General, RowUpTriangPacked> A;
    CheckComplexMatrix(A, non_herm, triang2, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColUpTriang> A;
    CheckRealMatrix(A, triang2, false);
  }

  {
    Matrix<Complex_wp, General, ColUpTriang> A;
    CheckComplexMatrix(A, non_herm, triang2, false);
  }

  {
    Matrix<Real_wp, General, ColUpTriangPacked> A;
    CheckRealMatrix(A, triang2, false);
  }

  {
    Matrix<Complex_wp, General, ColUpTriangPacked> A;
    CheckComplexMatrix(A, non_herm, triang2, false);
  }
#endif  

  {
    Matrix<Complex_wp, General, RowComplexSparse> A;
    CheckComplexMatrix(A, non_herm, triang, false);
  }
  
#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, General, ColComplexSparse> A;
    CheckComplexMatrix(A, non_herm, triang, false);
  }
#endif

  {
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    CheckComplexMatrix(A, non_herm, triang, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A;
    CheckComplexMatrix(A, non_herm, triang, false);
  }
#endif

  {
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    CheckComplexMatrix(A, non_herm, triang, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, General, ArrayColComplexSparse> A;
    CheckComplexMatrix(A, non_herm, triang, false);
  }
#endif

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;
    CheckComplexMatrix(A, non_herm, triang, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A;
    CheckComplexMatrix(A, non_herm, triang, false);
  }
#endif
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
