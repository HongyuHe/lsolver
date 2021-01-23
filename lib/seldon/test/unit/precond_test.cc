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
                          int m, int n, int nnz)
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
      {
        DISP(i); DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void CheckSymmetricPreconditioning(Matrix<T, Prop, Storage, Allocator>& A)
{
  SorPreconditioner<T> prec;
  int rapport = int(ceil(-log10(1e-16)/16));
  Iteration<typename ClassComplexType<T>::Treal> iter(rapport*30, 0.1*threshold);
  iter.ShowFullHistory();
  
  T zero;
  SetComplexZero(zero);
  
  int n = 100, nnz = 400;
  //int n = 10, nnz = 25;
  GenerateRandomMatrix(A, n, n, nnz);
  for (int i = 0; i < n; i++)
    {
      Real_wp sum = 1.0;
      for (int j = 0; j < n; j++)
	sum += abs(A(i, j));
      
      SetComplexReal(sum, A.Get(i, i));
    }

  
  prec.InitSymmetricPreconditioning();
  prec.SetParameterRelaxation(Real_wp(0.5));
  prec.SetNumberIterations(1);
  
  Vector<T> x(n), b(n), y;
  GenerateRandomVector(x, n);
  GenerateRandomVector(b, n);
  Mlt(A, x, b);
  y = x;
  x.Fill(zero);
  
  int success = CoCg(A, x, b, prec, iter);
  
  if (success != 0)
    {
      cout << "Sor preconditioning incorrect" << endl;
      abort();
    }  
  
  if (!EqualVector(x, y, 10.0*threshold))
    {
      cout << "Cg incorrect" << endl;
      abort();
    }
  
  x.Fill(0);
  success = QmrSym(A, x, b, prec, iter);
  
  if (success != 0)
    {
      cout << "Sor preconditioning incorrect" << endl;
      abort();
    }  
  
  if (!EqualVector(x, y, 10.0*threshold))
    {
      cout << "Cg incorrect" << endl;
      abort();
    }

  
  IlutPreconditioning<T> ilut;
  ilut.SetFactorisationType(ilut.ILUT);
  if (ilut.GetFactorisationType() != ilut.ILUT)
    {
      cout << "GetFactorisationType incorrect" << endl;
      abort();
    }
  
  ilut.SetPrintLevel(3);
  if (ilut.GetPrintLevel() != 3)
    {
      cout << "GetPrintLevel incorrect" << endl;
      abort();
    }
  
  ilut.SetSymmetricAlgorithm();
  ilut.SetDroppingThreshold(0.01);
  if (ilut.GetDroppingThreshold() != 0.01)
    {
      cout << "GetDroppingThreshold incorrect" << endl;
      abort();
    }
  
  IVect perm(n); perm.Fill();
  //FindSparseOrdering(A, perm, SparseMatrixOrdering::AUTO);
  
  // testing ilut as a direct solver
  ilut.SetDroppingThreshold(0);
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);
  ilut.Solve(A, b, x);
  if (!EqualVector(x, y, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
    
  // we keep initial matrix A
  ilut.SetDroppingThreshold(0.01);
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);  
  success = BiCgcr(A, x, b, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  if (!EqualVector(x, y, 10.0*threshold))
    {
      cout << "BiCgcr incorrect" << endl;
      abort();
    }
  
  ilut.SetFactorisationType(ilut.ILU_D);
  ilut.SetDiagonalCoefficient(0.1);
  if (ilut.GetDiagonalCoefficient() != 0.1)
    {
      cout << "GetDiagonalCoefficient incorrect" << endl;
      abort();
    }

  // testing ilut as a direct solver
  ilut.SetDroppingThreshold(0);
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);
  ilut.Solve(A, b, x);
  if (!EqualVector(x, y, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  ilut.SetDroppingThreshold(0.01);
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);
  success = Symmlq(A, x, b, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, y, 10.0*threshold))
    {
      cout << "Symmlq incorrect" << endl;
      abort();
    }
  
  ilut.SetFactorisationType(ilut.ILUT_K);
  
  ilut.SetFillLevel(n);
  ilut.SetDroppingThreshold(0);
  ilut.SetAdditionalFillNumber(n);
  if (ilut.GetAdditionalFillNumber() != n)
    {
      cout << "GetAdditionalFillNumber incorrect" << endl;
      abort();
    }
  
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);
  ilut.Solve(A, b, x);
  if (!EqualVector(x, y, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }

  ilut.SetDroppingThreshold(0.01);
  ilut.SetFillLevel(2);
  ilut.SetAdditionalFillNumber(10);
  
  if (ilut.GetFillLevel() != 2)
    {
      cout << "GetFillLevel incorrect" << endl;
      abort();
    }
  
  cout << "Testing MinRes" << endl;
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);  
  success = MinRes(A, x, b, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, y, 10.0*threshold))
    {
      cout << "MinRes incorrect" << endl;
      abort();
    }
  
  if (!IsComplexMatrix(A))
    {
      ilut.SetFactorisationType(ilut.ILU_K);
      ilut.SetFillLevel(n);
      
      ilut.FactorizeMatrix(perm, A, true);
      x.Fill(zero);
      ilut.Solve(A, b, x);
      if (!EqualVector(x, y, threshold))
	{
	  cout << "Ilut preconditioning incorrect" << endl;
	  abort();
	}

      ilut.SetFillLevel(2);
      ilut.FactorizeMatrix(perm, A, true);
      x.Fill(zero);  
      success = Cg(A, x, b, ilut, iter);
  
      if (success != 0)
	{
	  cout << "Ilut preconditioning incorrect" << endl;
	  abort();
	}
      
      if (!EqualVector(x, y, 10.0*threshold))
	{
	  cout << "Cg incorrect" << endl;
	  abort();
	}
    }

  Matrix<T, Prop, Storage, Allocator> B(n, n);
  IVect col_max(n);
  col_max.Fill();
  for (int i = 0; i < n; i++)
    {
      int k = i;
      for (int j = i-1; j >= 0; j--)
	if (A(i, j) != zero)
	  k = j;
      
      for (int j = i; j < n; j++)
	if (A(i, j) != zero)
	  col_max(i) = j;
      
      for (int j = k; j <= i; j++)
	col_max(j) = max(col_max(j), i);
    }
  
  for (int i = 0; i < n; i++)
    for (int j = i; j <= col_max(i); j++)
      B.Set(i, j, A(i, j));
  
  ilut.SetFactorisationType(ilut.ILU_0);
  
  ilut.FactorizeMatrix(perm, B, true);
  x.Fill(zero);
  ilut.Solve(A, b, x);
  if (!EqualVector(x, y, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  iter.SetMaxNumberIteration(30*rapport);
  ilut.SetFactorisationType(ilut.ILU_0);
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);  
  success = BiCg(A, x, b, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  if (!EqualVector(x, y, 10.0*threshold))
    {
      cout << "BiCg incorrect" << endl;
      abort();
    }

  ilut.SetFactorisationType(ilut.MILU_0);
  ilut.FactorizeMatrix(perm, B, true);
  x.Fill(zero);
  ilut.Solve(A, b, x);
  if (!EqualVector(x, y, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  iter.SetMaxNumberIteration(100*rapport);
  ilut.SetFactorisationType(ilut.MILU_0);
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);  
  success = Cgs(A, x, b, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  if (!EqualVector(x, y, 10.0*threshold))
    {
      cout << "Cgs incorrect" << endl;
      abort();
    }


}


template<class T, class Prop, class Storage, class Allocator>
void CheckGeneralPreconditioning(Matrix<T, Prop, Storage, Allocator>& A)
{
  SorPreconditioner<T> prec;
  int rapport = int(ceil(-log10(1e-16)/16));
  Iteration<typename ClassComplexType<T>::Treal> iter(rapport*50, 0.1*threshold);
  iter.ShowFullHistory();
  
  T zero, val;
  SetComplexZero(zero);
  
  int n = 100, nnz = 800;
  Matrix<T, Prop, Storage, Allocator> C;
  GenerateRandomMatrix(A, n, n, nnz);
  for (int i = 0; i < n; i++)
    if (A.GetRowSize(i) == 0)
      {
	GetRandNumber(val);
	A.Set(i, i, val);
      }
  
  C = A;
  for (int i = 0; i < n; i++)
    {
      Real_wp sum = 1.0;
      for (int j = 0; j < n; j++)
	sum += abs(A(i, j));
      
      SetComplexReal(sum, A.Get(i, i));
    }
  
  prec.InitUnSymmetricPreconditioning();
  prec.SetParameterRelaxation(Real_wp(0.5));
  prec.SetNumberIterations(1);
  iter.SetRestart(5);
  
  Vector<T> x(n), b(n), xc(n), bc(n), y, yc;
  GenerateRandomVector(x, n);
  GenerateRandomVector(xc, n);
  Mlt(A, x, b);
  Mlt(C, xc, bc);
  y = x;
  yc = xc;
  x.Fill(zero);
  int success = BiCgStab(A, x, b, prec, iter);
  
  if (success != 0)
    {
      cout << "Sor preconditioning incorrect" << endl;
      abort();
    }  
  
  if (!EqualVector(x, y, 30.0*threshold))
    {
      cout << "BiCgStab incorrect" << endl;
      abort();
    }
  
  iter.SetMaxNumberIteration(rapport*100);
  x.Fill(zero);
  success = Cgne(A, x, b, prec, iter);
    
  if ( (success != 0) || (!EqualVector(x, y, 30.0*threshold)) )
    {
      cout << "Cgne incorrect" << endl;
      abort();
    }
  
  iter.SetMaxNumberIteration(50*rapport);
  x.Fill(zero);
  success = Cgs(A, x, b, prec, iter);
  
  if ( (success != 0) || (!EqualVector(x, y, 30.0*threshold)) )
    {
      cout << "Cgs incorrect" << endl;
      abort();
    }
  
  x.Fill(zero);
  success = BiCg(A, x, b, prec, iter);
    
  if ( (success != 0) || (!EqualVector(x, y, 30.0*threshold)) )
    {
      cout << "BiCg incorrect" << endl;
      abort();
    }

  x.Fill(zero);
  success = Qmr(A, x, b, prec, iter);
  
  if ( (success != 0) || (!EqualVector(x, y, 30.0*threshold)) )
    {
      cout << "Qmr incorrect" << endl;
      abort();
    }

  x.Fill(zero);
  success = QCgs(A, x, b, prec, iter);
  
  if ( (success != 0) || (!EqualVector(x, y, 30.0*threshold)) )
    {
      cout << "QCgs incorrect" << endl;
      abort();
    }
  
  
  iter.SetMaxNumberIteration(50*rapport);
  IlutPreconditioning<T> ilut;
  ilut.SetFactorisationType(ilut.ILUT);
  if (ilut.GetFactorisationType() != ilut.ILUT)
    {
      cout << "GetFactorisationType incorrect" << endl;
      abort();
    }

  ilut.SetPivotThreshold(0.1);
  if (ilut.GetPivotThreshold() != 0.1)
    {
      cout << "GetPivotThreshold incorrect" << endl;
      abort();
    }
  
  ilut.SetPivotBlockInteger(n);
  if (ilut.GetPivotBlockInteger() != n)
    {
      cout << "GetPivotBlockInteger incorrect" << endl;
      abort();
    }
  
  ilut.SetPrintLevel(3);
  if (ilut.GetPrintLevel() != 3)
    {
      cout << "GetPrintLevel incorrect" << endl;
      abort();
    }
  
  ilut.SetUnsymmetricAlgorithm();
  ilut.SetDroppingThreshold(0.01);
  if (ilut.GetDroppingThreshold() != 0.01)
    {
      cout << "GetDroppingThreshold incorrect" << endl;
      abort();
    }
  
  IVect perm(n); perm.Fill();
  //FindSparseOrdering(C, perm, SparseMatrixOrdering::AUTO);
  //DISP(perm);
  
  // testing ilut as a direct solver
  ilut.SetDroppingThreshold(0);
  ilut.FactorizeMatrix(perm, C, true);
  xc.Fill(zero);
  ilut.Solve(C, bc, xc);
  if (!EqualVector(xc, yc, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
    
  // we keep initial matrix A
  ilut.SetDroppingThreshold(0.01);
  ilut.FactorizeMatrix(perm, C, true);
  xc.Fill(zero);  
  success = BiCgStabl(C, xc, bc, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  if (!EqualVector(xc, yc, 30.0*threshold))
    {
      cout << "BiCgStabl incorrect" << endl;
      abort();
    }
  
  ilut.SetFactorisationType(ilut.ILU_D);
  ilut.SetDiagonalCoefficient(0.1);
  if (ilut.GetDiagonalCoefficient() != 0.1)
    {
      cout << "GetDiagonalCoefficient incorrect" << endl;
      abort();
    }

  // testing ilut as a direct solver
  ilut.SetDroppingThreshold(0);
  ilut.FactorizeMatrix(perm, C, true);
  xc.Fill(zero);
  ilut.Solve(C, bc, xc);
  if (!EqualVector(xc, yc, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  ilut.SetDroppingThreshold(0.01);
  ilut.FactorizeMatrix(perm, C, true);
  xc.Fill(zero);
  success = Gcr(C, xc, bc, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }

  if (!EqualVector(xc, yc, 30.0*threshold))
    {
      cout << "Gcr incorrect" << endl;
      abort();
    }
  
  ilut.SetFactorisationType(ilut.ILUT_K);
  
  ilut.SetFillLevel(n);
  ilut.SetDroppingThreshold(0);
  ilut.SetAdditionalFillNumber(n);
  if (ilut.GetAdditionalFillNumber() != n)
    {
      cout << "GetAdditionalFillNumber incorrect" << endl;
      abort();
    }
  
  ilut.FactorizeMatrix(perm, C, true);
  xc.Fill(zero);
  ilut.Solve(C, bc, xc);
  if (!EqualVector(xc, yc, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }

  ilut.SetDroppingThreshold(0.01);
  ilut.SetFillLevel(2);
  ilut.SetAdditionalFillNumber(n);
  
  if (ilut.GetFillLevel() != 2)
    {
      cout << "GetFillLevel incorrect" << endl;
      abort();
    }
  
  cout << "Testing Gmres" << endl;
  ilut.FactorizeMatrix(perm, C, true);
  xc.Fill(zero);  
  
  iter.SetRestart(10);
  iter.ShowFullHistory();
  iter.SaveFullHistory(string("residu.txt"));
  success = Gmres(C, xc, bc, ilut, iter);
  //exit(0);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }

  if (!EqualVector(xc, yc, 30.0*threshold))
    {
      cout << "Gmres incorrect" << endl;
      abort();
    }
  
  if (!IsComplexMatrix(A))
    {
      ilut.SetFactorisationType(ilut.ILU_K);
      ilut.SetFillLevel(n);
      
      ilut.FactorizeMatrix(perm, A, true);
      x.Fill(zero);
      ilut.Solve(A, b, x);
      if (!EqualVector(x, y, threshold))
	{
	  cout << "Ilut preconditioning incorrect" << endl;
	  abort();
	}
      
      iter.SetMaxNumberIteration(200*rapport);
      ilut.SetFillLevel(2);
      ilut.FactorizeMatrix(perm, A, true);
      x.Fill(zero);  
      success = Lsqr(A, x, b, ilut, iter);
      
      if (success != 0)
	{
	  cout << "Ilut preconditioning incorrect" << endl;
	  abort();
	}
      
      if (!EqualVector(x, y, 30.0*threshold))
	{
	  cout << "Lsqr incorrect" << endl;
	  abort();
	}

      iter.SetMaxNumberIteration(50*rapport);
    }

  Matrix<T, Prop, Storage, Allocator> B(n, n);
  IVect col_max(n), col_min(n);
  col_max.Fill();
  for (int i = 0; i < n; i++)
    {
      int k = i;
      for (int j = i-1; j >= 0; j--)
	if (A(i, j) != zero)
	  k = j;
      
      col_min(i) = k;
      for (int j = i; j < n; j++)
	if (A(i, j) != zero)
	  col_max(i) = j;
      
      for (int j = k; j <= i; j++)
	col_max(j) = max(col_max(j), i);
    }
  
  col_min.Fill(0); col_max.Fill(n-1);
  for (int i = 0; i < n; i++)
    for (int j = col_min(i); j <= col_max(i); j++)
      B.Set(i, j, A(i, j));
  
  ilut.SetFactorisationType(ilut.ILU_0);
  
  ilut.FactorizeMatrix(perm, B, true);
  x.Fill(zero);
  ilut.Solve(A, b, x);
  if (!EqualVector(x, y, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  iter.SetMaxNumberIteration(30*rapport);
  ilut.SetFactorisationType(ilut.ILU_0);
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);  
  success = TfQmr(A, x, b, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  if (!EqualVector(x, y, 30.0*threshold))
    {
      cout << "TfQmr incorrect" << endl;
      abort();
    }

  ilut.SetFactorisationType(ilut.MILU_0);
  ilut.FactorizeMatrix(perm, B, true);
  x.Fill(zero);
  ilut.Solve(A, b, x);
  if (!EqualVector(x, y, threshold))
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  iter.SetMaxNumberIteration(100*rapport);
  ilut.SetFactorisationType(ilut.MILU_0);
  ilut.FactorizeMatrix(perm, A, true);
  x.Fill(zero);  
  success = BiCg(A, x, b, ilut, iter);
  
  if (success != 0)
    {
      cout << "Ilut preconditioning incorrect" << endl;
      abort();
    }
  
  if (!EqualVector(x, y, 30.0*threshold))
    {
      cout << "BiCg incorrect" << endl;
      abort();
    }
}

int main(int argc, char** argv)
{
  threshold = 1e-12;
  
  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckSymmetricPreconditioning(A);
  }
  
  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    CheckSymmetricPreconditioning(A);
  }

  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckGeneralPreconditioning(A);
  }

  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    CheckGeneralPreconditioning(A);
  }
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
