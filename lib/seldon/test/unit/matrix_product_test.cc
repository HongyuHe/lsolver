#include "SeldonLib.hxx"

using namespace Seldon;

typedef Vector<double> VectReal_wp;

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
                          int m, int n, int nnz, GhostIf<true>& sparse_form)
{
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  A.Clear();
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
                          int m, int n, int nnz, GhostIf<false>& sparse_form)
{
  A.Reallocate(m, n);
  T x;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
        GetRandNumber(x);
        A.Set(i, j, x);
      }
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomTriangular(Matrix<T, Prop, Storage, Allocator>& A,
			      int m, int n, bool lower)
{
  A.Reallocate(m, n);
  T x;
  if (lower)
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
    {
      T sum;
      SetComplexZero(sum);
      for (int j = 0; j < n; j++)
	if (j != i)
	  sum += abs(A(i, j));
      
      A.Get(i, i) += sum;
    }
  
}

template<class T1, class Prop1, class Storage1, class Allocator1,
         class T2, class Prop2, class Storage2, class Allocator2>
bool EqualMatrix(const Matrix<T1, Prop1, Storage1, Allocator1>& A,
                 const Matrix<T2, Prop2, Storage2, Allocator2>& B)
{
  if ( (A.GetM() != B.GetM())  || (A.GetN() != B.GetN()) )
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j) - B(i, j))))
        {
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
          return false;
        }
  
  return true;
}

template<class T1, class Prop1, class Storage1, class Allocator1,
         class T2, class Prop2, class Storage2, class Allocator2,
         class T3, class Prop3, class Storage3, class Allocator3>
void MltTest(const Matrix<T1, Prop1, Storage1, Allocator1>& A,
             const Matrix<T2, Prop2, Storage2, Allocator2>& B,
             Matrix<T3, Prop3, Storage3, Allocator3>& C)
{
  C.Fill(0);
  T3 val;
  for (int i = 0; i < C.GetM(); i++)
    for (int j = 0; j < C.GetN(); j++)
      {
        SetComplexZero(val);
        for (int k = 0; k < A.GetN(); k++)
          val += A(i, k)*B(k, j);
        
        C.Set(i, j, val);
      }
}

template<class T, class Prop, class Storage, class Allocator>
void Display(const Matrix<T, Prop, Storage, Allocator>& A)
{
  
}

template<class T, class Prop, class Allocator>
void Display(const Matrix<T, Prop, RowSparse, Allocator>& A)
{
  DISP(A.GetM()); DISP(A.GetN());
  cout << "Ptr  = " << endl;
  for (int i = 0; i <= A.GetM(); i++)
    cout << A.GetPtr()[i]+1 << " ";

  cout << endl << endl << "Ind = " << endl;
  for (int i = 0; i < A.GetDataSize(); i++)
    cout << A.GetInd()[i]+1 << " ";

  cout << endl << endl << "Val = " << endl;
  for (int i = 0; i < A.GetDataSize(); i++)
    cout << A.GetData()[i] << " ";
  
  cout << endl << endl;
}

// testing real sparse matrices
template<class T, class Prop, class Storage, class Allocator, bool sparse>
void CheckRealMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                     GhostIf<sparse>& sparse_form, bool all_test)
{
  int m = 64, n = 54, k = 59, nnz = 300;
  //int m = 8, n = 6, k = 7, nnz = 15;
  Matrix<T, Prop, Storage, Allocator> B, At, Bt, C, D, C2, D2;
  GenerateRandomMatrix(A, m, n, nnz, sparse_form);
  GenerateRandomMatrix(B, n, k, nnz-4, sparse_form);
  GenerateRandomMatrix(C, m, k, nnz-12, sparse_form);
  
  D = C;
  C2 = C;
  Mlt(A, B, C);
  MltTest(A, B, C2);
  if (!EqualMatrix(C, C2))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Real_wp alpha, beta;
  GetRandNumber(alpha);
  GetRandNumber(beta);
  Mlt(alpha, A, B, C);
  Mlt(alpha, C2);
  if (!EqualMatrix(C, C2))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  MltTest(A, B, C2);

  C = D;
  MltAdd(alpha, A, B, beta, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j)) > threshold)
          || isnan(abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j))))
        {
          cout << "MltAdd incorrect" << endl;
          abort();
        }
  
  Transpose(A, At); Transpose(B, Bt);

  C = D;
  MltAdd(alpha, SeldonNoTrans, A, SeldonNoTrans, B, beta, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j)) > threshold)
          || isnan(abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j))))
        {
          cout << "MltAdd incorrect" << endl;
          abort();
        }
  
  if (!sparse || all_test)
    {
      C = D;
      MltAdd(alpha, SeldonTrans, At, SeldonNoTrans, B, beta, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }
    }
  
  C = D;
  MltAdd(alpha, SeldonNoTrans, A, SeldonTrans, Bt, beta, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j)) > threshold)
          || isnan(abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j))))
        {
          cout << "MltAdd incorrect" << endl;
          abort();
        }
  

  if (!sparse || all_test)
    {
      C = D;
      MltAdd(alpha, SeldonTrans, At, SeldonTrans, Bt, beta, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - beta*D(i, j) - alpha*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }
    }
  
  
}

// checking complex sparse matrix
template<class T, class Prop, class Storage, class Allocator, bool sparse>
void CheckComplexMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                        GhostIf<sparse>& sparse_form, bool all_test)
{
  int m = 64, n = 54, k = 59, nnz = 300;
  //int m = 4, n = 4, k = 4, nnz = 10;
  Matrix<T, Prop, Storage, Allocator> B, At, Bt, Act, Bct, C, D, C2, D2;
  GenerateRandomMatrix(A, m, n, nnz, sparse_form);
  GenerateRandomMatrix(B, n, k, nnz-12, sparse_form);
  GenerateRandomMatrix(C, m, k, nnz+11, sparse_form);
  
  D = C;
  C2 = C;
  Mlt(A, B, C);
  MltTest(A, B, C2);
  if (!EqualMatrix(C, C2))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Complex_wp alphac, betac;
  GetRandNumber(alphac);
  GetRandNumber(betac);
  Mlt(alphac, A, B, C);
  Mlt(alphac, C2);
  if (!EqualMatrix(C, C2))
    {
      DISP(alphac); DISP(A); DISP(B); DISP(C); DISP(C2);
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Complex_wp one(1, 0), zero(0, 0);
  MltAdd(one, A, B, zero, C);
  MltTest(A, B, C2);
  if (!EqualMatrix(C, C2))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }
  
  C = D;
  MltAdd(alphac, A, B, betac, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
          || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
        {
          cout << "MltAdd incorrect" << endl;
          abort();
        }
  
  // testing with alpha, beta complexes
  C = D;
  MltAdd(alphac, A, B, betac, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
          || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
        {
          cout << "MltAdd incorrect" << endl;
          abort();
        }
  
  C = D;
  MltAdd(alphac, SeldonNoTrans, A, SeldonNoTrans, B, betac, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
          || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
        {
          cout << "MltAdd incorrect" << endl;
          abort();
        }

  Transpose(A, At);
  TransposeConj(A, Act);
  Transpose(B, Bt);
  TransposeConj(B, Bct);
  if (!sparse || all_test)
    {
      C = D;
      MltAdd(alphac, SeldonTrans, At, SeldonNoTrans, B, betac, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }
    }
  
  C = D;
  MltAdd(alphac, SeldonNoTrans, A, SeldonTrans, Bt, betac, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
          || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
        {
          cout << "MltAdd incorrect" << endl;
          abort();
        }
  

  if (!sparse || all_test)
    {
      C = D;
      MltAdd(alphac, SeldonTrans, At, SeldonTrans, Bt, betac, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }  
    }

  if (!sparse || all_test)
    {
      C = D;
      MltAdd(alphac, SeldonNoTrans, A, SeldonConjTrans, Bct, betac, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }
      
      C = D;
      MltAdd(alphac, SeldonConjTrans, Act, SeldonNoTrans, B, betac, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }
      
      C = D;
      MltAdd(alphac, SeldonTrans, At, SeldonConjTrans, Bct, betac, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }
      
      
      C = D;
      MltAdd(alphac, SeldonConjTrans, Act, SeldonTrans, Bt, betac, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }  
      
      C = D;
      MltAdd(alphac, SeldonConjTrans, Act, SeldonConjTrans, Bct, betac, C);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j)) > threshold)
              || isnan(abs(C(i, j) - betac*D(i, j) - alphac*C2(i, j))))
            {
              cout << "MltAdd incorrect" << endl;
              abort();
            }  
    }
}

// checking dense symmetric matrix
template<class T, class Prop, class Storage, class Allocator,
	 class T2, class Prop2, class Storage2, class Allocator2>
void CheckSymmetricMatrix(Matrix<T, Prop, Storage, Allocator>& A,
			  Matrix<T2, Prop2, Storage2, Allocator2>& B)
{
  int m = 35, nnz = 0;
  GhostIf<false> dense;
  Matrix<T2, Prop2, Storage2, Allocator2> C, C0, AB(m, m), BA(m, m);
  
  GenerateRandomMatrix(A, m, m, nnz, dense);
  GenerateRandomMatrix(B, m, m, nnz, dense);
  GenerateRandomMatrix(C, m, m, nnz, dense);
  
  for (int i = 0; i < m; i++)
    {
      T val;
      SetComplexReal(realpart(A(i, i)), val);
      A.Set(i, i, val);
      SetComplexReal(realpart(B(i, i)), val);
      B.Set(i, i, val);
      SetComplexReal(realpart(C(i, i)), val);
      C.Set(i, i, val);
    }
  
  T alpha, beta;
  GetRandNumber(alpha);
  GetRandNumber(beta);
  
  C0 = C;
  MltAdd(SeldonLeft, alpha, A, B, beta, C);
  MltTest(A, B, AB);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(C(i, j) - beta*C0(i, j) - alpha*AB(i, j)) > threshold)
          || isnan(abs(C(i, j) - beta*C0(i, j) - alpha*AB(i, j))))
	{
          DISP(C(i, j)); DISP(beta); DISP(C0(i, j)); DISP(alpha); DISP(AB(i, j));
          DISP(abs(C(i, j) - beta*C0(i, j) - alpha*AB(i, j)));
	  cout << "MltAdd incorrect" << endl;
	  abort();
	}

  C = C0;
  MltAdd(SeldonRight, alpha, A, B, beta, C);
  MltTest(B, A, BA);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(C(i, j) - beta*C0(i, j) - alpha*BA(i, j)) > threshold)
          || isnan(abs(C(i, j) - beta*C0(i, j) - alpha*BA(i, j))))
	{
	  cout << "MltAdd incorrect" << endl;
	  abort();
	}

}
 

template<class T, class Prop, class Storage, class Allocator,
	 class T2, class Prop2, class Storage2, class Allocator2>
void CheckTriangularMatrix(Matrix<T, Prop, Storage, Allocator>& A,
			   Matrix<T2, Prop2, Storage2, Allocator2>& B, bool lower)
{
  T alpha, beta, one;
  SetComplexOne(one);
  GetRandNumber(alpha);
  GetRandNumber(beta);
  
  int m = 40, nnz = 0;
  GhostIf<false> dense;
  Matrix<T2, Prop2, Storage2, Allocator2> C, AB(m, m), BA(m, m), Ct, AtB(m, m), BAt(m, m);
  GenerateRandomMatrix(C, m, m, nnz, dense);
  GenerateRandomTriangular(A, m, m, lower);
  Transpose(C, Ct);
  
  B = C;
  MltTest(A, B, AB);
  MltTest(B, A, BA);
  B = C;
  MltTest(Ct, A, AtB);
  Transpose(AtB);

  MltTest(A, Ct, BAt);
  Transpose(BAt);

  Mlt(SeldonLeft, alpha, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*AB(i, j)) > threshold)
          || isnan(abs(B(i, j) - alpha*AB(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  beta = one / alpha;
  Solve(SeldonLeft, beta, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  B = C;
  Mlt(SeldonRight, alpha, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*BA(i, j)) > threshold)
          || isnan(abs(B(i, j) - alpha*BA(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonRight, beta, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  B = C;
  Mlt(SeldonLeft, alpha, SeldonNoTrans, SeldonNonUnit, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*AB(i, j)) > threshold)
          || isnan(abs(B(i, j) - alpha*AB(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonLeft, beta, SeldonNoTrans, SeldonNonUnit, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  B = C;
  Mlt(SeldonRight, alpha, SeldonNoTrans, SeldonNonUnit, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*BA(i, j)) > threshold)
          || isnan(abs(B(i, j) - alpha*BA(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonRight, beta, SeldonNoTrans, SeldonNonUnit, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  B = C;
  Mlt(SeldonLeft, alpha, SeldonTrans, SeldonNonUnit, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*AtB(i, j)) > threshold)
          || isnan(abs(B(i, j) - alpha*AtB(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonLeft, beta, SeldonTrans, SeldonNonUnit, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  B = C;
  Mlt(SeldonRight, alpha, SeldonTrans, SeldonNonUnit, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*BAt(i, j)) > threshold)
          || isnan(abs(B(i, j) - alpha*BAt(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonRight, beta, SeldonTrans, SeldonNonUnit, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }


  if (IsComplexMatrix(A))
    {
      Matrix<T2, Prop2, Storage2, Allocator2> Cc, AcB(m, m), BAc(m, m);
      TransposeConj(C, Cc);
      MltTest(Cc, A, AcB);
      TransposeConj(AcB);

      MltTest(A, Cc, BAc);
      TransposeConj(BAc);
      
      B = C;
      Mlt(SeldonLeft, alpha, SeldonConjTrans, SeldonNonUnit, A, B);
      
      for (int i = 0; i < m; i++)
	for (int j = 0; j < m; j++)
	  if ((abs(B(i, j) - alpha*AcB(i, j)) > threshold)
              || isnan(abs(B(i, j) - alpha*AcB(i, j))))
	    {
	      cout << "Mlt incorrect" << endl;
	      abort();
	    }
      
      Solve(SeldonLeft, beta, SeldonConjTrans, SeldonNonUnit, A, B);
      
      if (!EqualMatrix(B, C))
	{
	  cout << "Solve incorrect" << endl;
	  abort();
	}
      
      B = C;
      Mlt(SeldonRight, alpha, SeldonConjTrans, SeldonNonUnit, A, B);
      
      for (int i = 0; i < m; i++)
	for (int j = 0; j < m; j++)
	  if ((abs(B(i, j) - alpha*BAc(i, j)) > threshold) || isnan(abs(B(i, j) - alpha*BAc(i, j))))
	    {
	      cout << "Mlt incorrect" << endl;
	      abort();
	    }
      
      Solve(SeldonRight, beta, SeldonConjTrans, SeldonNonUnit, A, B);
      
      if (!EqualMatrix(B, C))
	{
	  cout << "Solve incorrect" << endl;
	  abort();
	}
    }

  // testing unit diagonal
  Vector<T> diagA(m);
  for (int i = 0; i < m; i++)
    {
      diagA(i) = A(i, i);
      A.Set(i, i, one);
    }
  
  B = C;
  MltTest(A, B, AB);
  MltTest(B, A, BA);
  B = C;
  MltTest(Ct, A, AtB);
  Transpose(AtB);

  MltTest(A, Ct, BAt);
  Transpose(BAt);

  for (int i = 0; i < m; i++)
    A.Set(i, i, diagA(i));
  
  B = C;
  Mlt(SeldonLeft, alpha, SeldonNoTrans, SeldonUnit, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*AB(i, j)) > threshold) || isnan(abs(B(i, j) - alpha*AB(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonLeft, beta, SeldonNoTrans, SeldonUnit, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  B = C;
  Mlt(SeldonRight, alpha, SeldonNoTrans, SeldonUnit, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*BA(i, j)) > threshold) || isnan(abs(B(i, j) - alpha*BA(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonRight, beta, SeldonNoTrans, SeldonUnit, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  B = C;
  Mlt(SeldonLeft, alpha, SeldonTrans, SeldonUnit, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*AtB(i, j)) > threshold) || isnan(abs(B(i, j) - alpha*AtB(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonLeft, beta, SeldonTrans, SeldonUnit, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  B = C;
  Mlt(SeldonRight, alpha, SeldonTrans, SeldonUnit, A, B);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*BAt(i, j)) > threshold) || isnan(abs(B(i, j) - alpha*BAt(i, j))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Solve(SeldonRight, beta, SeldonTrans, SeldonUnit, A, B);
  
  if (!EqualMatrix(B, C))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  if (IsComplexMatrix(A))
    {
      Matrix<T2, Prop2, Storage2, Allocator2> Cc, AcB(m, m), BAc(m, m);

      for (int i = 0; i < m; i++)
	A.Set(i, i, one);
      
      TransposeConj(C, Cc);
      MltTest(Cc, A, AcB);
      TransposeConj(AcB);

      MltTest(A, Cc, BAc);
      TransposeConj(BAc);
      
      for (int i = 0; i < m; i++)
	A.Set(i, i, diagA(i));
      
      B = C;
      Mlt(SeldonLeft, alpha, SeldonConjTrans, SeldonUnit, A, B);
      
      for (int i = 0; i < m; i++)
	for (int j = 0; j < m; j++)
	  if ((abs(B(i, j) - alpha*AcB(i, j)) > threshold) || isnan(abs(B(i, j) - alpha*AcB(i, j))))
	    {
	      cout << "Mlt incorrect" << endl;
	      abort();
	    }
      
      Solve(SeldonLeft, beta, SeldonConjTrans, SeldonUnit, A, B);
      
      if (!EqualMatrix(B, C))
	{
	  cout << "Solve incorrect" << endl;
	  abort();
	}
      
      B = C;
      Mlt(SeldonRight, alpha, SeldonConjTrans, SeldonUnit, A, B);
      
      for (int i = 0; i < m; i++)
	for (int j = 0; j < m; j++)
	  if ((abs(B(i, j) - alpha*BAc(i, j)) > threshold) || isnan(abs(B(i, j) - alpha*BAc(i, j))))
	    {
	      cout << "Mlt incorrect" << endl;
	      abort();
	    }
      
      Solve(SeldonRight, beta, SeldonConjTrans, SeldonUnit, A, B);
      
      if (!EqualMatrix(B, C))
	{
	  cout << "Solve incorrect" << endl;
	  abort();
	}
    }
}

template<class T0, class Prop0, class Storage0, class Allocator0,
         class T1, class Prop1, class Storage1, class Allocator1>
void CheckSymMatrix(Matrix<T0, Prop0, Storage0, Allocator0>& A,
                    Matrix<T1, Prop1, Storage1, Allocator1>& B)
{
  int m = 45;
  GhostIf<false> dense;
  Matrix<T0, Prop0, Storage0, Allocator0> BBt(m, m), C;
  Matrix<T1, Prop1, Storage1, Allocator1> Bt;
  GenerateRandomMatrix(A, m, m, 0, dense);
  GenerateRandomMatrix(B, m, m, 0, dense);
  Transpose(B, Bt);
  MltTest(B, Bt, BBt);
  
  T0 alpha, beta;
  GetRandNumber(alpha); GetRandNumber(beta);

  // testing product C = B B^T (symmetric matrix)
  C = A;
  MltAdd(alpha, SeldonNoTrans, B, SeldonTrans, B, beta, C);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T0 val_ref = beta*A(i, j) + alpha*BBt(i, j);
        if ((abs(val_ref - C(i, j)) > threshold) || isnan(abs(val_ref - C(i, j))))
          {
            DISP(i); DISP(j); DISP(C(i, j));
            DISP(val_ref);
            abort();
          }
      }

  C = A;
  MltAdd(alpha, SeldonTrans, Bt, SeldonTrans, B, beta, C);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T0 val_ref = beta*A(i, j) + alpha*BBt(i, j);
        if ((abs(val_ref - C(i, j)) > threshold) || isnan(abs(val_ref - C(i, j))))
          {
            DISP(i); DISP(j); DISP(C(i, j));
            DISP(val_ref);
            abort();
          }
      }

  C = A;
  MltAdd(alpha, SeldonTrans, Bt, SeldonNoTrans, Bt, beta, C);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T0 val_ref = beta*A(i, j) + alpha*BBt(i, j);
        if ((abs(val_ref - C(i, j)) > threshold) || isnan(abs(val_ref - C(i, j))))
          {
            DISP(i); DISP(j); DISP(C(i, j));
            DISP(val_ref);
            abort();
          }
      }

  C = A;
  Mlt(B, Bt, C);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T0 val_ref = BBt(i, j);
        if ((abs(val_ref - C(i, j)) > threshold) || isnan(abs(val_ref - C(i, j))))
          {
            DISP(i); DISP(j); DISP(C(i, j));
            DISP(val_ref);
            abort();
          }
      }
  
  // testing addition of symmetric matrices
  GenerateRandomMatrix(C, m, m, 0, dense);

  Matrix<T0, Prop0, Storage0, Allocator0> D(C);
  Add(alpha, A, C);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T0 val_ref = D(i, j) + alpha*A(i, j);
        if ((abs(val_ref - C(i, j)) > threshold) || isnan(abs(val_ref - C(i, j))))
          {
            DISP(i); DISP(j); DISP(C(i, j));
            DISP(val_ref);
            abort();
          }
      }
  
}


template<class T0, class Prop0, class Storage0, class Allocator0,
         class T1, class Prop1, class Storage1, class Allocator1>
void CheckHermMatrix(Matrix<T0, Prop0, Storage0, Allocator0>& A,
                     Matrix<T1, Prop1, Storage1, Allocator1>& B)
{
  int m = 15;
  GhostIf<false> dense;
  Matrix<T0, Prop0, Storage0, Allocator0> BBt(m, m), C;
  Matrix<T1, Prop1, Storage1, Allocator1> Bt;
  GenerateRandomMatrix(A, m, m, 0, dense);
  GenerateRandomMatrix(B, m, m, 0, dense);
  TransposeConj(B, Bt);
  MltTest(B, Bt, BBt);

  // testing addition of symmetric matrices
  GenerateRandomMatrix(C, m, m, 0, dense);

  Real_wp alpha(2.3);
  Matrix<T0, Prop0, Storage0, Allocator0> D(C);
  Add(alpha, A, C);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T0 val_ref = D(i, j) + alpha*A(i, j);
        if ((abs(val_ref - C(i, j)) > threshold) || isnan(abs(val_ref - C(i, j))))
          {
            DISP(i); DISP(j); DISP(C(i, j));
            DISP(val_ref);
            abort();
          }
      }
  
}

template<class T0, class Prop0, class Storage0, class Allocator0>
void CheckTriangMatrix(Matrix<T0, Prop0, Storage0, Allocator0>& A, bool lower)
{
  int m = 10;
  Matrix<T0, Prop0, Storage0, Allocator0> B, C(m, m);
  Matrix<T0, General, RowMajor> AB(m, m);
  GenerateRandomTriangular(A, m, m, lower);
  GenerateRandomTriangular(B, m, m, lower);
  
  MltTest(A, B, AB);  
  Mlt(A, B, C);
  //DISP(A); DISP(B); DISP(C); DISP(AB);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(AB(i, j) - C(i, j)) > threshold) || isnan(abs(AB(i, j) - C(i, j))))
        {
          cout << " Mlt incorrect" <<endl;
          abort();
        }
  
  T0 alpha; GetRandNumber(alpha);
  C = B;
  Add(alpha, A, B);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(B(i, j) - alpha*A(i, j) - C(i,j)) > threshold) || isnan(abs(B(i, j) - alpha*A(i, j) - C(i,j))))
        {
          cout << " Add incorrect" <<endl;
          abort();
        }
}

int main(int argc, char** argv)
{
  threshold = 1e-11;
  
  GhostIf<true> sparse;
  GhostIf<false> dense;
  
  {
    // testing function Add 
    // in configurations not treated in other unitary tests
    Matrix<Real_wp, General, ArrayRowSparse> B, C;
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    
    int m = 34, n = m, nnz = 160;
    GenerateRandomMatrix(A, m, n, nnz, sparse);
    
    GenerateRandomMatrix(B, m, n, nnz, sparse);
    
    Real_wp alpha;
    GetRandNumber(alpha);
    C = B;
    Add(alpha, A, B);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++)
        if ((abs(B(i, j) - C(i, j) - alpha*A(i, j)) > threshold) || isnan(abs(B(i, j) - C(i, j) - alpha*A(i, j))))
          {
            cout << "Add incorrect" << endl;
            abort();
          }
    
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> Bs, Cs;
    GenerateRandomMatrix(Bs, m, n, nnz+12, sparse);    
    GenerateRandomMatrix(Cs, m, n, nnz+27, sparse);    
    
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> Ac, Ac2;
    GenerateRandomMatrix(Ac, m, n, nnz-16, sparse);
    Ac2 = Ac;
    
    AddMatrix(alpha, Bs, Cs, Ac);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++)
        if ((abs(Ac(i, j) - Ac2(i, j) - alpha*Complex_wp(Bs(i, j), Cs(i, j))) > threshold)
            || isnan(abs(Ac(i, j) - Ac2(i, j) - alpha*Complex_wp(Bs(i, j), Cs(i, j)))))
          {
            cout << "Add incorrect" << endl;
            abort();
          }        

    m = 41; n = 35;
    GenerateRandomMatrix(B, m, n, nnz, sparse);
    GenerateRandomMatrix(C, m, n, nnz+23, sparse);
    
    Matrix<Complex_wp, General, ArrayRowSparse> Mc, Mc2;
    GenerateRandomMatrix(Mc, m, n, nnz-16, sparse);
    
    Mc2 = Mc;
    AddMatrix(alpha, B, C, Mc);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        if ((abs(Mc(i, j) - Mc2(i, j) - alpha*Complex_wp(B(i, j), C(i, j))) > threshold)
            || isnan(abs(Mc(i, j) - Mc2(i, j) - alpha*Complex_wp(B(i, j), C(i, j)))))
          {
            cout << "Add incorrect" << endl;
            abort();
          }    
    
    // testing function Mlt for specific configurations
    int k = 28;
    Matrix<Real_wp, General, RowMajor> Bd, Cd, Cd2(m, k);
    GenerateRandomMatrix(Bd, n, k, nnz, dense);
    Cd.Reallocate(m, k);
    Mlt(B, Bd, Cd);
    MltTest(B, Bd, Cd2);
    if (!EqualMatrix(Cd, Cd2))
      {
        cout << "Mlt incorrect" << endl;
        abort();
      }    
  }  
  
  // testing matrix-matrix products
  
  {
    Matrix<Real_wp, General, RowMajor> A;
    CheckRealMatrix(A, dense, true);
  }
  
  {
    Matrix<Real_wp, General, ColMajor> A;
    CheckRealMatrix(A, dense, true);
  }

  {
    Matrix<Complex_wp, General, RowMajor> A;
    CheckComplexMatrix(A, dense, true);
  }

  {
    Matrix<Complex_wp, General, ColMajor> A;
    CheckComplexMatrix(A, dense, true);
  }
    
  {
    Matrix<Real_wp, General, RowSparse> A, At;
  
    // product with dense matrices
    Matrix<Real_wp, General, RowMajor> B, C, C0, D, Bt;
    Real_wp alpha, beta;
    A.Clear();
    int m = 60, n = 42, k = 51, nnz = 2000;
    GenerateRandomMatrix(A, m, n, nnz, sparse);
    GenerateRandomMatrix(B, n, k, nnz, dense);
    A.WriteText("A.dat");
    B.WriteText("B.dat");
    GenerateRandomMatrix(C0, m, k, nnz, dense);
    D.Reallocate(m, k);
    GetRandNumber(alpha);
    GetRandNumber(beta);
    C = C0;
    MltAdd(alpha, A, B, beta, C);
    MltTest(A, B, D);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < k; j++)
        if ((abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j)) > threshold) || isnan(abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j))))
          {
            cout << "MltAdd incorrect" << endl;
            abort();
          }
    
    Transpose(A, At);
    Transpose(B, Bt);
    At.WriteText("At.dat");
    Bt.WriteText("Bt.dat");
    C = C0;
    MltAdd(alpha, SeldonTrans, At, SeldonNoTrans, B, beta, C);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < k; j++)
        if ((abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j)) > threshold)
            || isnan(abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j))))
          {
            cout << "MltAdd incorrect" << endl;
            abort();
          }
    
    C = C0;
    MltAdd(alpha, SeldonTrans, At, SeldonTrans, Bt, beta, C);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < k; j++)
        if ((abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j)) > threshold)
            || isnan(abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j))))
          {
            DISP(i); DISP(j); DISP(C(i, j));
            DISP(beta*C0(i, j) + alpha*D(i,j));
            cout << "MltAdd incorrect" << endl;
            abort();
          }

    C = C0;
    MltAdd(alpha, SeldonNoTrans, A, SeldonNoTrans, B, beta, C);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < k; j++)
        if ((abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j)) > threshold)
            || isnan(abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j))))
          {
            cout << "MltAdd incorrect" << endl;
            abort();
          }

    C = C0;
    MltAdd(alpha, SeldonNoTrans, A, SeldonTrans, Bt, beta, C);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < k; j++)
        if ((abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j)) > threshold)
            || isnan(abs(C(i, j) - beta*C0(i, j) - alpha*D(i,j))))
          {
            cout << "MltAdd incorrect" << endl;
            abort();
          }
    
  }

  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckRealMatrix(A, sparse, true);
  }

  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    CheckComplexMatrix(A, sparse, true);
  }

  // testing other matrix-matrix products
  {
    Matrix<Real_wp, Symmetric, RowSym> A;
    Matrix<Real_wp, General, RowMajor> B;
    CheckSymmetricMatrix(A, B);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSym> A;
    Matrix<Real_wp, General, ColMajor> B;
    CheckSymmetricMatrix(A, B);
  }
#endif

  {
    Matrix<Complex_wp, Hermitian, RowHerm> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckSymmetricMatrix(A, B);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Hermitian, ColHerm> A;
    Matrix<Complex_wp, General, ColMajor> B;
    CheckSymmetricMatrix(A, B);
  }
#endif
  
  // Multiplication and resolution with triangular systems

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColLoTriang> A;
    Matrix<Real_wp, General, ColMajor> B;
    CheckTriangularMatrix(A, B, true);
  }

  {
    Matrix<Complex_wp, General, ColLoTriang> A;
    Matrix<Complex_wp, General, ColMajor> B;
    CheckTriangularMatrix(A, B, true);
  }
#endif

  {
    Matrix<Real_wp, General, RowLoTriang> A;
    Matrix<Real_wp, General, RowMajor> B;
    CheckTriangularMatrix(A, B, true);
  }

  {
    Matrix<Complex_wp, General, RowLoTriang> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckTriangularMatrix(A, B, true);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColUpTriang> A;
    Matrix<Real_wp, General, ColMajor> B;
    CheckTriangularMatrix(A, B, false);
  }

  {
    Matrix<Complex_wp, General, ColUpTriang> A;
    Matrix<Complex_wp, General, ColMajor> B;
    CheckTriangularMatrix(A, B, false);
  }
#endif

  {
    Matrix<Real_wp, General, RowUpTriang> A;
    Matrix<Real_wp, General, RowMajor> B;
    CheckTriangularMatrix(A, B, false);
  }

  {
    Matrix<Complex_wp, General, RowUpTriang> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckTriangularMatrix(A, B, false);
  }

  // testing symmetric matrices
  {
    Matrix<Real_wp, Symmetric, RowSymPacked> A;
    Matrix<Real_wp, General, RowMajor> B;
    CheckSymMatrix(A, B);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymPacked> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckSymMatrix(A, B);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSymPacked> A;
    Matrix<Real_wp, General, RowMajor> B;
    CheckSymMatrix(A, B);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSymPacked> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckSymMatrix(A, B);
  }
#endif

  {
    Matrix<Real_wp, Symmetric, RowSym> A;
    Matrix<Real_wp, General, RowMajor> B;
    CheckSymMatrix(A, B);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSym> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckSymMatrix(A, B);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSym> A;
    Matrix<Real_wp, General, RowMajor> B;
    CheckSymMatrix(A, B);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSym> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckSymMatrix(A, B);
  }
#endif

  // testing triangular matrices
  {
    Matrix<Real_wp, General, RowLoTriang> A;
    CheckTriangMatrix(A, true);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColLoTriang> A;
    CheckTriangMatrix(A, true);
  }
#endif

  {
    Matrix<Real_wp, General, RowUpTriang> A;
    CheckTriangMatrix(A, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColUpTriang> A;
    CheckTriangMatrix(A, false);
  }
#endif

  {
    Matrix<Real_wp, General, RowLoTriangPacked> A;
    CheckTriangMatrix(A, true);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColLoTriangPacked> A;
    CheckTriangMatrix(A, true);
  }
#endif

  {
    Matrix<Real_wp, General, RowUpTriangPacked> A;
    CheckTriangMatrix(A, false);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColUpTriangPacked> A;
    CheckTriangMatrix(A, false);
  }
#endif

  // testing hermitian matrices
  {
    Matrix<Complex_wp, Hermitian, RowHerm> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckHermMatrix(A, B);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Hermitian, ColHerm> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckHermMatrix(A, B);
  }
#endif

  {
    Matrix<Complex_wp, Hermitian, RowHermPacked> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckHermMatrix(A, B);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Hermitian, ColHermPacked> A;
    Matrix<Complex_wp, General, RowMajor> B;
    CheckHermMatrix(A, B);
  }
#endif
  
  cout << "All tests passed successfully" << endl;

  return 0;
}
