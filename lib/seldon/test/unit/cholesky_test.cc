// Copyright (C) 2010 Marc Duruflé
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


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

template<class T>
bool CheckVector(Vector<T>& x)
{
  bool test = true;
  T icplx, coef;
  SetComplexReal(x.GetM(), icplx);
  SetComplexOne(coef); coef /= icplx;
  for (int i = 0; i < x.GetM(); i++)
    {
      SetComplexReal(i, icplx);
      if ((abs(x(i) - icplx)*coef > threshold) || isnan(x(i)))
        test = false;
    }

  return test;
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y,
		 Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > eps) || isnan(x(i)) || isnan(y(i)))
      return false;
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void CheckDenseCholesky(Matrix<T, Prop, Storage, Allocator>& A)
{
  Matrix<T, Prop, Storage, Allocator> Adense;
  Matrix<T, Symmetric, ArrayRowSymSparse> B;
  Vector<T> x, y, b;
  
  B.ReadText("test/matrix/MhSparse.dat");
    
  Copy(B, Adense);
  
  GetCholesky(Adense);
  
  x.Reallocate(B.GetN());
  b.Reallocate(B.GetN());
  y.Reallocate(B.GetN());
  x.Fill();
  Mlt(B, x, b);
  x = b;
  SolveCholesky(SeldonNoTrans, Adense, x);
  SolveCholesky(SeldonTrans, Adense, x);
  if (!CheckVector(x))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  y.Fill();
  SolveCholesky(SeldonNoTrans, Adense, y);
  SolveCholesky(SeldonTrans, Adense, y);

  x = y;
  MltCholesky(SeldonTrans, Adense, x);
  MltCholesky(SeldonNoTrans, Adense, x);

  if (!CheckVector(x))
    {
      cout << "MltCholesky incorrect" << endl;
      abort();
    }
  
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomTriangular(Matrix<T, Prop, Storage, Allocator>& A,
			      int m, int n, bool low)
{
  T one, zero;
  SetComplexOne(one);
  SetComplexZero(zero);

  A.Reallocate(m, n);
  A.Fill(zero);
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
  
  for (int i = 0; i < m; i++)
    A.Get(i, i) = one;
  
}

template<class T, class Prop, class Storage, class Allocator>
void CheckHermitianCholesky(Matrix<T, Prop, Storage, Allocator>& A)
{
  T one, zero;
  SetComplexOne(one);
  SetComplexZero(zero);

  int n = 5;
  Matrix<T, General, RowMajor> L, B(n, n), C(n, n), invL(n, n);
  GenerateRandomTriangular(L, n, n, true);
  invL = L;
  GetInverse(invL);
  
  B.Fill(zero);
  MltAdd(one, SeldonNoTrans, L, SeldonConjTrans, L, zero, B);
  
  Matrix<T, Prop, Storage, Allocator> Adense(n, n);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      Adense.Val(i, j) = B(i, j);
  
  GetCholesky(Adense);
  
  C.Fill(zero);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      C(i, j) = Adense.Val(i, j);
  
  Vector<T> x(n), b(n), y(n);
  GenerateRandomVector(x, n);
  y = x;
  Mlt(B, x, b);

  x = b;
  SolveCholesky(SeldonNoTrans, Adense, x);
  Mlt(SeldonNoTrans, invL, b, y);
  if (!EqualVector(x, y))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  x = b;
  SolveCholesky(SeldonConjTrans, Adense, x);
  Mlt(SeldonConjTrans, invL, b, y);
  if (!EqualVector(x, y))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  x = b;
  SolveCholesky(SeldonTrans, Adense, x);
  Mlt(SeldonTrans, invL, b, y);
  if (!EqualVector(x, y))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  x = y;
  MltCholesky(SeldonConjTrans, Adense, x);
  Mlt(SeldonConjTrans, L, y, b);

  if (!EqualVector(x, b))
    {
      cout << "MltCholesky incorrect" << endl;
      abort();
    }
  
  x = y;
  MltCholesky(SeldonNoTrans, Adense, x);
  Mlt(SeldonNoTrans, L, y, b);

  if (!EqualVector(x, b))
    {
      cout << "MltCholesky incorrect" << endl;
      abort();
    }
  
  x = y;
  MltCholesky(SeldonTrans, Adense, x);
  Mlt(SeldonTrans, L, y, b);
  
  if (!EqualVector(x, b))
    {
      cout << "MltCholesky incorrect" << endl;
      abort();
    }
  
}


template<class T, class Prop, class Storage, class Allocator>
void CheckSparseCholesky(Matrix<T, Prop, Storage, Allocator>& A)
{
  // symmetric matrix is read in a file
  A.ReadText("test/matrix/MhSparse.dat");
  
  // creation of a right hand side b = A*[0;1;...;n-1]
  Vector<T> x(A.GetM()), b(A.GetM()), y(A.GetM());
  x.Fill();
  Mlt(A, x, b);
  
  /*********************
   * ArrayRowSymSparse *
   *********************/
  
  // testing Cholesky factorisation using Seldon function
  GetCholesky(A);
  
  x = b;
  SolveCholesky(SeldonNoTrans, A, x);
  SolveCholesky(SeldonTrans, A, x);

  if (!CheckVector(x))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  y.Fill();
  SolveCholesky(SeldonNoTrans, A, y);
  SolveCholesky(SeldonTrans, A, y);

  x = y;
  MltCholesky(SeldonTrans, A, x);
  MltCholesky(SeldonNoTrans, A, x);

  if (!CheckVector(x))
    {
      cout << "MltCholesky incorrect" << endl;
      abort();
    }

#ifdef SELDON_WITH_CHOLMOD
  
  /***********
   * Cholmod *
   ***********/

  // testing Cholesky factorisation using Cholmod
  MatrixCholmod mat_chol;
  mat_chol.ShowMessages();
			
  A.ReadText("test/matrix/MhSparse.dat");

  GetCholesky(A, mat_chol);

  x = b;
  SolveCholesky(SeldonNoTrans, mat_chol, x);
  SolveCholesky(SeldonTrans, mat_chol, x);

  if (!CheckVector(x))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  y.Fill();
  SolveCholesky(SeldonNoTrans, mat_chol, y);
  SolveCholesky(SeldonTrans, mat_chol, y);
  
  x = y;
  MltCholesky(SeldonTrans, mat_chol, x);
  MltCholesky(SeldonNoTrans, mat_chol, x);

  if (!CheckVector(x))
    {
      cout << "MltCholesky incorrect" << endl;
      abort();
    }
  
  mat_chol.Clear();
#endif


  // testing SparseCholeskySolver
  SparseCholeskySolver<T> mat_lu;

  A.ReadText("test/matrix/MhSparse.dat");

  mat_lu.Factorize(A);

  x = b;
  mat_lu.Solve(SeldonNoTrans, x);
  mat_lu.Solve(SeldonTrans, x);

  if (!CheckVector(x))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  y.Fill();
  mat_lu.Solve(SeldonNoTrans, y);
  mat_lu.Solve(SeldonTrans, y);
  
  x = y;
  mat_lu.Mlt(SeldonTrans, x);
  mat_lu.Mlt(SeldonNoTrans, x);

  if (!CheckVector(x))
    {
      cout << "MltCholesky incorrect" << endl;
      abort();
    }

}

int main(int argc, char** argv)
{
  InitSeldon(argc, argv);
  
  threshold = 1e-12;
  DISP(threshold);
  
  {
    Matrix<Real_wp, Symmetric, RowSymPacked> A;
    CheckDenseCholesky(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSymPacked> A;
    CheckDenseCholesky(A);
  }
#endif

  {
    Matrix<Real_wp, Symmetric, RowSym> A;
    CheckDenseCholesky(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSym> A;
    CheckDenseCholesky(A);
  }
#endif
    
  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckSparseCholesky(A);
  }

  {
    Matrix<Complex_wp, Hermitian, RowHermPacked> A;
    CheckHermitianCholesky(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Hermitian, ColHermPacked> A;
    CheckHermitianCholesky(A);
  }
#endif

  {
    Matrix<Complex_wp, Hermitian, RowHerm> A;
    CheckHermitianCholesky(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Hermitian, ColHerm> A;
    CheckHermitianCholesky(A);
  }
#endif

  cout << "All tests passed successfully" << endl;
    
  return FinalizeSeldon();
}
