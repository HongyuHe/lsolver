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


void GenerateRandomPermutation(int n, IVect& permut)
{
  Vector<bool> NumUsed(n);
  NumUsed.Fill(false);
  permut.Reallocate(n);
  permut.Fill(-1);
  int nb = 0;
  // premiere iteration
  for (int i = 0; i < n; i++)
    {
      int i2 = rand()%n;
      if (!NumUsed(i2))
        {
          NumUsed(i2) = true;
          permut(i) = i2;
          nb++;
        }
    }
  
  while (nb < n)
    {
      // on recupere les numeros non-selectionnes
      IVect non_selec(n-nb);
      int k = 0;
      for (int i = 0; i < n; i++)
        if (!NumUsed(i))
          non_selec(k++) = i;
      
      // iteration suivante
      for (int i = 0; i < n; i++)
        if (permut(i) == -1)
          {
            int i2 = rand()%(n-nb);
            int j = non_selec(i2);
            if (!NumUsed(j))
              {
                NumUsed(j) = true;
                permut(i) = j;
                nb++;
              }
          }
    }
}

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
void GenerateRandomVector(Vector<T, VectSparse>& x, int n, int nnz)
{
  for (int i = 0; i < nnz; i++)
    {
      int j = rand()%n;
      GetRandNumber(x.Get(j));
    }
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz, GhostIf<true>& sparse_form)
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
                              int m, int n, bool low)
{
  A.Reallocate(m, n);
  T x;
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

template<class T, class Prop, class Storage, class Allocator>
void CheckModifRowCol(Matrix<T, Prop, Storage, Allocator>& A,
                      GhostIf<true>& sparse_form)
{
  int irow = rand()%A.GetM();
  int jcol = rand()%A.GetN();
  
  Vector<T, VectSparse> row, col;
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  
  // testing GetRow and GetCol
  GetRow(A, irow, row);
  for (int i = 0; i < A.GetN(); i++)
    if ((A(irow, i) != row(i)) || isnan(A(irow, i)) || isnan(row(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }

  GetCol(A, jcol, col);
  for (int i = 0; i < A.GetM(); i++)
    if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }

  // testing SetRow and SetCol when the pattern is not modified
  A = B;
  for (int i = 0; i < row.GetM(); i++)
    GetRandNumber(row.Value(i));
  
  SetRow(row, irow, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == irow)
          {
            if ((A(irow, j) != row(j)) || isnan(A(irow, j)) || isnan(row(j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }
 
  A = B;
  for (int i = 0; i < col.GetM(); i++)
    GetRandNumber(col.Value(i));
  
  SetCol(col, jcol, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (j == jcol)
          {
            if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      } 

  // testing with a different pattern  
  GenerateRandomVector(row, A.GetN(), A.GetN()/5);
  GenerateRandomVector(col, A.GetM(), A.GetM()/5);
  
  A = B;
  SetRow(row, irow, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == irow)
          {
            if ((A(irow, j) != row(j)) || isnan(A(irow, j)) || isnan(row(j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }

  A = B;
  SetCol(col, jcol, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (j == jcol)
          {
            if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }   
}

template<class T, class Prop, class Storage, class Allocator>
void CheckUpdateMatrix(Matrix<T, Prop, Storage, Allocator>& A,
		       Matrix<T, Prop, Storage, Allocator>& B,
		       Vector<T>& row, Vector<T>& col)
{
  T alpha;
  GetRandNumber(alpha);
  DISP(alpha);  
  
  A = B;
  Rank1Update(alpha, col, row, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j) - alpha*col(i)*row(j)) > threshold)
          || isnan(A(i, j)) || isnan(B(i, j)) || isnan(alpha*col(i)*row(j)))
	{
	  cout << "Rank1Update incorrect" << endl;
	  abort();
	}
  

}

template<class T, class Prop, class Storage, class Allocator>
void CheckUpdateMatrix(Matrix<complex<T>, Prop, Storage, Allocator>& A,
		       Matrix<complex<T>, Prop, Storage, Allocator>& B,
		       Vector<complex<T> >& row, Vector<complex<T> >& col)
{
  complex<T> alpha;
  GetRandNumber(alpha);
  DISP(alpha);  
  
  A = B;
  Rank1Update(alpha, col, row, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j) - alpha*col(i)*row(j)) > threshold)
          || isnan(A(i, j)) || isnan(B(i, j)) || isnan(alpha*col(i)*row(j)))
	{
	  cout << "Rank1Update incorrect" << endl;
	  abort();
	}

  A = B;
  Rank1Update(alpha, col, SeldonUnconj, row, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j) - alpha*row(j)*col(i)) > threshold)
          || isnan(A(i, j)) || isnan(B(i, j)) || isnan(alpha*col(i)*row(j)))
	{
	  cout << "Rank1Update incorrect" << endl;
	  abort();
	}

  A = B;
  Rank1Update(alpha, col, SeldonConj, row, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j) - alpha*conjugate(row(j))*col(i)) > threshold)
          || isnan(A(i, j)) || isnan(B(i, j)) || isnan(alpha*conjugate(col(i)*row(j))))
	{
	  cout << "Rank1Update incorrect" << endl;
	  abort();
	}
  
}

template<class T, class Prop, class Storage, class Allocator>
void CheckModifRowCol(Matrix<T, Prop, Storage, Allocator>& A,
                      GhostIf<false>& sparse_form)
{
  int irow = rand()%A.GetM();
  int jcol = rand()%A.GetN();
  
  Vector<T> row, col;
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  
  GetRow(A, irow, row);
  for (int i = 0; i < A.GetN(); i++)
    if ((A(irow, i) != row(i)) || isnan(A(irow, i)) || isnan(row(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }

  GetCol(A, jcol, col);
  for (int i = 0; i < A.GetM(); i++)
    if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }
  
  GenerateRandomVector(row, A.GetN());
  GenerateRandomVector(col, A.GetM());
  
  SetRow(row, irow, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == irow)
          {
            if ((A(irow, j) != row(j)) || isnan(A(irow, j)) || isnan(row(j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }

  A = B;
  SetCol(col, jcol, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (j == jcol)
          {
            if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }  
  
  CheckUpdateMatrix(A, B, row, col);  
}


template<class T, class Prop, class Storage, class Allocator>
void CheckGeneralMatrix(Matrix<T, Prop, Storage, Allocator>& A)
{
  int m = 70, n = 36, nnz = 200;

  GhostIf<Storage::Sparse> sparse_form;
  GenerateRandomMatrix(A, m, n, nnz, sparse_form);
  
  CheckModifRowCol(A, sparse_form);

  IVect col_perm;
  GenerateRandomPermutation(n, col_perm);

  IVect row_perm;
  GenerateRandomPermutation(m, row_perm);
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  ApplyInversePermutation(A, row_perm, col_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((B(i, j) != A(row_perm(i), col_perm(j)))
          || isnan(B(i, j)) || isnan(A(row_perm(i), col_perm(j))))
        {
          DISP(i); DISP(j);
          DISP(B(i, j)); DISP(A(row_perm(i), col_perm(j)));
          cout << "ApplyInversePermutation incorrect" << endl;
          abort();
        }
  
  A = B;
  ApplyPermutation(A, row_perm, col_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(row_perm(i), col_perm(j)))
          || isnan(A(i, j)) || isnan(B(row_perm(i), col_perm(j))))
        {
          DISP(i); DISP(j);
          DISP(A(i, j)); DISP(B(row_perm(i), col_perm(j)));
          cout << "ApplyPermutation incorrect" << endl;
          abort();
        }

  VectReal_wp Drow, Dcol;
  GenerateRandomVector(Drow, m);
  GenerateRandomVector(Dcol, n);
  
  A = B;
  ScaleMatrix(A, Drow, Dcol);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Drow(i)*Dcol(j)*B(i,j)) > threshold)
          || isnan(A(i, j)) || isnan(Drow(i)*Dcol(j)*B(i,j)))
        {
          cout << "Scale Matrix incorrect" << endl;
          abort();
        }
  
  A = B;
  ScaleLeftMatrix(A, Drow);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Drow(i)*B(i,j)) > threshold)
          || isnan(A(i, j)) || isnan(Drow(i)*B(i, j)))
        {
          cout << "ScaleLeft Matrix incorrect" << endl;
          abort();
        }

  A = B;
  ScaleRightMatrix(A, Dcol);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Dcol(j)*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(Dcol(j)*B(i, j)))
        {
          cout << "ScaleRight Matrix incorrect" << endl;
          abort();
        }

  // testing Mlt and Add
  T alpha;
  GetRandNumber(alpha);
  A = B;
  Mlt(alpha, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - alpha*B(i, j)) > threshold) || isnan(A(i, j)) || isnan(alpha*B(i, j)))
        {
          cout << "Mlt incorrect" << endl;
          abort();
        }
  
  Real_wp beta;
  GetRandNumber(beta);
  A = B;
  Mlt(beta, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - beta*B(i, j)) > threshold) || isnan(A(i, j)) || isnan(beta*B(i, j)))
        {
          cout << "Mlt incorrect" << endl;
          abort();
        }
  
  B.FillRand();
  Mlt(Real_wp(1)/RAND_MAX, B);
  Matrix<T, Prop, Storage, Allocator> C;
  C = A;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i, j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  A = C;
  Add(beta, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - beta*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(beta*B(i, j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }
  
  GenerateRandomMatrix(B, m, n, nnz+23, sparse_form);
  
  A = C;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i, j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }
  
  Real_wp norm_A = MaxAbs(A);
  Real_wp val_max = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (abs(A(i, j)) > val_max)
        val_max = abs(A(i, j));
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "MaxAbs incorrect" << endl;
      abort();
    }

  val_max = 0;
  norm_A = NormInf(A);
  for (int i = 0; i < m; i++)
    {
      Real_wp sum(0);
      for (int j = 0; j < n; j++)
        sum += abs(A(i, j));
      
      if (sum > val_max)
        val_max = sum;
    }
    
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "NormInf incorrect" << endl;
      abort();
    }

  val_max = 0;
  norm_A = Norm1(A);
  for (int j = 0; j < n; j++)
    {
      Real_wp sum(0);
      for (int i = 0; i < m; i++)
        sum += abs(A(i, j));
      
      if (sum > val_max)
        val_max = sum;
    }
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "Norm1 incorrect" << endl;
      abort();
    }
  
  Transpose(A, B);
  if ( (A.GetM() != B.GetN()) || (A.GetN() != B.GetM()))
    {
      cout << "Transpose incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(j, i)) || isnan(A(i, j)) || isnan(B(j, i)))
        {
          cout << "Transpose incorrect" << endl;
          abort();
        }
  
  C = A;
  Transpose(A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((A(i, j) != B(i, j) ) || isnan(A(i, j)) || isnan(B(i, j)))
        {
          cout << "Transpose incorrect" << endl;
          abort();
        }
  
  Conjugate(A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if ((A(i, j) != conjugate(B(i, j)) ) || isnan(A(i, j)) || isnan(conjugate(B(i, j)) ))
        {
          cout << "Conjugate incorrect" << endl;
          abort();
        }
  
  TransposeConj(C);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((A(i, j) != C(i, j) ) || isnan(A(i, j)) || isnan(C(i, j)))
        {
          cout << "TransposeConj incorrect" << endl;
          abort();
        }
}

template<class T, class Prop, class Storage, class Allocator>
void CheckUpdateSym(Matrix<T, Prop, Storage, Allocator>& A,
		    Matrix<T, Prop, Storage, Allocator>& B,
		    Vector<T>& row, Vector<T>& col)
{
  T alpha;
  GetRandNumber(alpha);
  
  A = B;
  Rank1Update(alpha, row, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j) - alpha*row(i)*row(j)) > threshold)
          || isnan(A(i, j)) || isnan(B(i, j)) || isnan(alpha*row(i)*row(j)))
	{
	  cout << "Rank1Update incorrect" << endl;
	  abort();
	}


  A = B;
  Rank2Update(alpha, row, col, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j) - alpha*(row(i)*col(j) + col(i)*row(j))) > threshold)
          || isnan(A(i, j) ) || isnan(B(i, j)) || isnan(alpha*(row(i)*col(j) + col(i)*row(j))))
	{
	  cout << "Rank2Update incorrect" << endl;
	  abort();
	}
}

template<class T, class Prop, class Storage, class Allocator>
void CheckModifSym(Matrix<T, Prop, Storage, Allocator>& A,
                   GhostIf<true>& sparse_form)
{
  int irow = rand()%A.GetM();
  int jcol = rand()%A.GetN();
  
  Vector<T, VectSparse> row, col;
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  
  // testing GetRow, GetCol
  GetRow(A, irow, row);
  for (int i = 0; i < A.GetN(); i++)
    if ((A(irow, i) != row(i)) || isnan(A(irow, i)) || isnan(row(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }

  GetCol(A, jcol, col);
  for (int i = 0; i < A.GetM(); i++)
    if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }
  
  // testing SetRow, SetCol for a same pattern
  for (int i = 0; i < row.GetM(); i++)
    GetRandNumber(row.Value(i));

  A = B;
  SetRow(row, irow, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == irow)
          {
            if ((A(irow, j) != row(j)) || isnan(A(irow, j)) || isnan(row(j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else if (j == irow)
          {
            if ((A(i, irow) != row(i)) || isnan(A(i, irow)) || isnan(row(i)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }

  for (int i = 0; i < col.GetM(); i++)
    GetRandNumber(col.Value(i));
 
  A = B;
  SetCol(col, jcol, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (j == jcol)
          {
            if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else if (i == jcol)
          {
            if ((A(jcol, j) != col(j)) || isnan(A(jcol, j)) || isnan(col(j)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }  
  
  // testing with a different pattern
  GenerateRandomVector(row, A.GetN(), A.GetN()/4);
  GenerateRandomVector(col, A.GetM(), A.GetM()/4);
  
  A = B;
  SetRow(row, irow, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == irow)
          {
            if ((A(irow, j) != row(j)) || isnan(A(irow, j)) || isnan(row(j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else if (j == irow)
          {
            if ((A(i, irow) != row(i)) || isnan(A(i, irow)) || isnan(row(i)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }

  A = B;
  SetCol(col, jcol, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (j == jcol)
          {
            if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else if (i == jcol)
          {
            if ((A(jcol, j) != col(j)) || isnan(A(jcol, j)) || isnan(col(j)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }  
}

template<class T, class Prop, class Storage, class Allocator>
void CheckModifSym(Matrix<T, Prop, Storage, Allocator>& A,
                   GhostIf<false>& sparse_form)
{
  int irow = rand()%A.GetM();
  int jcol = rand()%A.GetN();
  
  Vector<T> row, col;
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  
  GetRow(A, irow, row);
  for (int i = 0; i < A.GetN(); i++)
    if ((A(irow, i) != row(i)) || isnan(A(irow, i)) || isnan(row(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }

  GetCol(A, jcol, col);
  for (int i = 0; i < A.GetM(); i++)
    if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }

  GenerateRandomVector(row, A.GetN());
  GenerateRandomVector(col, A.GetM());
  
  SetRow(row, irow, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == irow)
          {
            if ((A(irow, j) != row(j)) || isnan(A(irow, j)) || isnan(row(j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else if (j == irow)
          {
            if ((A(i, irow) != row(i)) || isnan(A(i, irow)) || isnan(row(i)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }

  A = B;
  SetCol(col, jcol, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (j == jcol)
          {
            if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else if (i == jcol)
          {
            if ((A(jcol, j) != col(j)) || isnan(A(jcol, j)) || isnan(col(j)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }  
  
  CheckUpdateSym(A, B, row, col);
}


template<class T, class Prop, class Storage, class Allocator>
void CheckSymmetricMatrix(Matrix<T, Prop, Storage, Allocator>& A)
{
  int m = 48, n = 48, nnz = 200;
  
  GhostIf<Storage::Sparse> sparse_form;
  GenerateRandomMatrix(A, m, n, nnz, sparse_form);
  
  CheckModifSym(A, sparse_form);

  IVect row_perm;
  GenerateRandomPermutation(m, row_perm);
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  ApplyInversePermutation(A, row_perm, row_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((B(i, j) != A(row_perm(i), row_perm(j)))
          || isnan(B(i, j)) || isnan(A(row_perm(i), row_perm(j))))
        {
          DISP(i); DISP(j);
          DISP(B(i, j)); DISP(A(row_perm(i), row_perm(j)));
          cout << "ApplyInversePermutation incorrect" << endl;
          abort();
        }
  
  A = B;
  ApplyPermutation(A, row_perm, row_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(row_perm(i), row_perm(j)))
          || isnan(A(i, j)) || isnan(B(row_perm(i), row_perm(j))))
        {
          DISP(i); DISP(j);
          DISP(A(i, j)); DISP(B(row_perm(i), row_perm(j)));
          cout << "ApplyPermutation incorrect" << endl;
          abort();
        }

  VectReal_wp Drow;
  GenerateRandomVector(Drow, m);
  
  A = B;
  ScaleMatrix(A, Drow, Drow);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Drow(i)*Drow(j)*B(i,j)) > threshold)
          || isnan(A(i, j)) || isnan(Drow(i)*Drow(j)*B(i,j)))
        {
          cout << "Scale Matrix incorrect" << endl;
          abort();
        }

  // testing Mlt and Add
  T alpha;
  GetRandNumber(alpha);
  A = B;
  Mlt(alpha, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(alpha*B(i, j)))
        {
          cout << "Mlt incorrect" << endl;
          abort();
        }
  
  Real_wp beta;
  GetRandNumber(beta);
  A = B;
  Mlt(beta, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - beta*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(beta*B(i,j)))
        {
          cout << "Mlt incorrect" << endl;
          abort();
        }
  
  B.FillRand();
  Mlt(Real_wp(1)/RAND_MAX, B);
  Matrix<T, Prop, Storage, Allocator> C;
  C = A;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i, j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  A = C;
  Add(beta, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - beta*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(beta*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }
  
  GenerateRandomMatrix(B, m, n, nnz+23, sparse_form);
  
  A = C;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }
  
  Real_wp norm_A = MaxAbs(A);
  Real_wp val_max = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (abs(A(i, j)) > val_max)
        val_max = abs(A(i, j));
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "MaxAbs incorrect" << endl;
      abort();
    }

  val_max = 0;
  norm_A = NormInf(A);
  for (int i = 0; i < m; i++)
    {
      Real_wp sum(0);
      for (int j = 0; j < n; j++)
        sum += abs(A(i, j));
      
      if (sum > val_max)
        val_max = sum;
    }
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "NormInf incorrect" << endl;
      abort();
    }

  val_max = 0;
  norm_A = Norm1(A);
  for (int j = 0; j < n; j++)
    {
      Real_wp sum(0);
      for (int i = 0; i < m; i++)
        sum += abs(A(i, j));
      
      if (sum > val_max)
        val_max = sum;
    }
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "Norm1 incorrect" << endl;
      abort();
    }

  Transpose(A, B);
  if ( (A.GetM() != B.GetN()) || (A.GetN() != B.GetM()))
    {
      cout << "Transpose incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(j, i)) || isnan(A(i, j)) || isnan(B(j, i)))
        {
          cout << "Transpose incorrect" << endl;
          abort();
        }
  
  C = A;
  Transpose(A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((A(i, j) != B(i, j) ) || isnan(A(i, j)) || isnan(B(i, j)))
        {
          cout << "Transpose incorrect" << endl;
          abort();
        }
  
  Conjugate(A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if ((A(i, j) != conjugate(B(i, j)) ) || isnan(A(i, j)) || isnan(B(i, j)))
        {
          cout << "Conjugate incorrect" << endl;
          abort();
        }
  
  TransposeConj(C);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((A(i, j) != C(i, j) ) || isnan(A(i, j)) || isnan(C(i, j)))
        {
          cout << "TransposeConj incorrect" << endl;
          abort();
        }
}

template<class T, class Prop, class Storage, class Allocator>
void CheckModifHerm(Matrix<T, Prop, Storage, Allocator>& A,
                    GhostIf<false>& sparse_form)
{
  int irow = rand()%A.GetM();
  int jcol = rand()%A.GetN();
  
  Vector<T> row, col;
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  
  GetRow(A, irow, row);
  for (int i = 0; i < A.GetN(); i++)
    if ((A(irow, i) != row(i)) || isnan(A(irow, i)) || isnan(row(i)))
      {
        DISP(row(i)); DISP(A(irow, i));
        cout << "GetRow incorrect" << endl;
        abort();
      }

  GetCol(A, jcol, col);
  for (int i = 0; i < A.GetM(); i++)
    if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }

  GenerateRandomVector(row, A.GetN());
  GenerateRandomVector(col, A.GetM());
  
  SetRow(row, irow, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == irow)
          {
            if ((A(irow, j) != row(j)) || isnan(A(irow, j)) || isnan(row(j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else if (j == irow)
          {
            if ((A(i, irow) != conjugate(row(i))) || isnan(A(i, irow)) || isnan(row(i)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }

  A = B;
  SetCol(col, jcol, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (j == jcol)
          {
            if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else if (i == jcol)
          {
            if ((A(jcol, j) != conjugate(col(j))) || isnan(A(jcol, j)) || isnan(col(j)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }  

  Real_wp alpha;
  GetRandNumber(alpha);
  
  A = B;
  Rank1Update(alpha, row, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j) - alpha*row(i)*conjugate(row(j))) > threshold)
          || isnan(A(i, j)) || isnan(B(i, j)) || isnan(alpha*row(i)*conjugate(row(j))))
	{
	  cout << "Rank1Update incorrect" << endl;
	  abort();
	}


  A = B;
  Rank2Update(alpha, row, col, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j) - alpha*(row(i)*conjugate(col(j)) + col(i)*conjugate(row(j)))) > threshold)
          || isnan(A(i, j)) || isnan(B(i, j)) || isnan(alpha*(row(i)*conjugate(col(j)) + col(i)*conjugate(row(j)))))
	{
	  cout << "Rank2Update incorrect" << endl;
	  abort();
	}

}


template<class T, class Prop, class Storage, class Allocator>
void CheckHermitianMatrix(Matrix<T, Prop, Storage, Allocator>& A)
{
  int m = 4, n = 4, nnz = 5;
  
  GhostIf<Storage::Sparse> sparse_form;
  GenerateRandomMatrix(A, m, n, nnz, sparse_form);
  for (int i = 0; i < m; i++)
    A.Set(i, i, T(real(A(i, i)), 0));
  
  CheckModifHerm(A, sparse_form);

  IVect row_perm;
  GenerateRandomPermutation(m, row_perm);
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  ApplyInversePermutation(A, row_perm, row_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((B(i, j) != A(row_perm(i), row_perm(j)))
          || isnan(B(i, j)) || isnan(A(row_perm(i), row_perm(j))))
        {
          DISP(i); DISP(j);
          DISP(B(i, j)); DISP(A(row_perm(i), row_perm(j)));
          cout << "ApplyInversePermutation incorrect" << endl;
          abort();
        }
  
  A = B;
  ApplyPermutation(A, row_perm, row_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(row_perm(i), row_perm(j)))
          || isnan(A(i, j)) || isnan(B(row_perm(i), row_perm(j))))
        {
          DISP(i); DISP(j);
          DISP(A(i, j)); DISP(B(row_perm(i), row_perm(j)));
          cout << "ApplyPermutation incorrect" << endl;
          abort();
        }

  VectReal_wp Drow;
  GenerateRandomVector(Drow, m);
  
  A = B;
  ScaleMatrix(A, Drow, Drow);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Drow(i)*Drow(j)*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(Drow(i)*Drow(j)*B(i, j)))
        {
          cout << "Scale Matrix incorrect" << endl;
          abort();
        }

  // testing Mlt and Add
  Real_wp alpha;
  GetRandNumber(alpha);
  A = B;
  Mlt(alpha, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(alpha*B(i, j)))
        {
          cout << "Mlt incorrect" << endl;
          abort();
        }
  
  Real_wp beta;
  GetRandNumber(beta);
  A = B;
  Mlt(beta, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - beta*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(beta*B(i, j)))
        {
          cout << "Mlt incorrect" << endl;
          abort();
        }
  
  B.FillRand();
  Mlt(Real_wp(1)/RAND_MAX, B);
  Matrix<T, Prop, Storage, Allocator> C;
  C = A;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold)
          || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  A = C;
  Add(beta, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - beta*B(i,j)) > threshold)
          || isnan(A(i, j)) || isnan(C(i, j)) || isnan(beta*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }
  
  GenerateRandomMatrix(B, m, n, nnz+23, sparse_form);
  for (int i = 0; i < m; i++)
    B.Set(i, i, T(real(B(i, i)), 0));
  
  A = C;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold)
          || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i,j)))
        {
          DISP(i); DISP(j); DISP(A(i, j)); DISP(C(i, j)); DISP(B(i, j));
          cout << "Add incorrect" << endl;
          abort();
        }
  
  Real_wp norm_A = MaxAbs(A);
  Real_wp val_max = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (abs(A(i, j)) > val_max)
        val_max = abs(A(i, j));
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      DISP(A);
      DISP(val_max); DISP(norm_A);
      cout << "MaxAbs incorrect" << endl;
      abort();
    }

  val_max = 0;
  norm_A = NormInf(A);
  for (int i = 0; i < m; i++)
    {
      Real_wp sum(0);
      for (int j = 0; j < n; j++)
        sum += abs(A(i, j));
      
      if (sum > val_max)
        val_max = sum;
    }
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "NormInf incorrect" << endl;
      abort();
    }

  val_max = 0;
  norm_A = Norm1(A);
  for (int j = 0; j < n; j++)
    {
      Real_wp sum(0);
      for (int i = 0; i < m; i++)
        sum += abs(A(i, j));
      
      if (sum > val_max)
        val_max = sum;
    }
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "Norm1 incorrect" << endl;
      abort();
    }

  Transpose(A, B);
  if ( (A.GetM() != B.GetN()) || (A.GetN() != B.GetM()))
    {
      cout << "Transpose incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(j, i)) || isnan(A(i, j)) || isnan(B(j, i)))
        {
          if ((i == j) && A(i, j) == conjugate(B(j, i)))
            {
            }
          else
            {
              cout << "Transpose incorrect" << endl;
              abort();
            }
        }
  
  C = A;
  Transpose(A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((A(i, j) != B(i, j) ) || isnan(A(i, j)) || isnan(B(i, j)))
        {
          cout << "Transpose incorrect" << endl;
          abort();
        }
  
  Conjugate(A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if ((A(i, j) != conjugate(B(i, j)) ) || isnan(A(i, j)) || isnan(B(i, j)))
        {
          cout << "Conjugate incorrect" << endl;
          abort();
        }
  
  TransposeConj(C);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((A(i, j) != C(i, j) ) || isnan(A(i, j)) || isnan(C(i, j)))
        {
          cout << "TransposeConj incorrect" << endl;
          abort();
        }

}

template<class T, class Prop, class Storage, class Allocator>
void CheckModifTriangular(Matrix<T, Prop, Storage, Allocator>& A,
                          GhostIf<false>& sparse_form, bool low)
{
  int irow = rand()%A.GetM();
  int jcol = rand()%A.GetN();
  
  Vector<T> row, col;
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  
  GetRow(A, irow, row);
  for (int i = 0; i < A.GetN(); i++)
    if ((A(irow, i) != row(i)) || isnan(A(irow, i)) || isnan(row(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }

  GetCol(A, jcol, col);
  for (int i = 0; i < A.GetM(); i++)
    if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
      {
        cout << "GetRow incorrect" << endl;
        abort();
      }
  
  GenerateRandomVector(row, A.GetN());
  GenerateRandomVector(col, A.GetM());
  
  if (!low)
    {
      for (int i = 0; i < irow; i++)
        SetComplexZero(row(i));

      for (int i = jcol+1; i < A.GetN(); i++)
        SetComplexZero(col(i));
    }
  else
    {
      for (int i = irow+1; i < A.GetN(); i++)
        SetComplexZero(row(i));

      for (int i = 0; i < jcol; i++)
        SetComplexZero(col(i));
    }
  
  SetRow(row, irow, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == irow)
          {
            if ((A(irow, j) != row(j)) || isnan(A(irow, j)) || isnan(row(j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }

  A = B;
  SetCol(col, jcol, A);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (j == jcol)
          {
            if ((A(i, jcol) != col(i)) || isnan(A(i, jcol)) || isnan(col(i)))
              {
                cout << "SetCol incorrect" << endl;
                abort();
              }
          }
        else
          {
            if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
              {
                cout << "SetRow incorrect" << endl;
                abort();
              }
          }
      }  
}

template<class T, class Prop, class Storage, class Allocator>
void CheckTriangularMatrix(Matrix<T, Prop, Storage, Allocator>& A, bool low)
{
  int m = 39, n = 39;

  GhostIf<Storage::Sparse> sparse_form;
  GenerateRandomTriangular(A, m, n, low);
  
  CheckModifTriangular(A, sparse_form, low);

  VectReal_wp Drow, Dcol;
  GenerateRandomVector(Drow, m);
  GenerateRandomVector(Dcol, n);
  
  Matrix<T, Prop, Storage, Allocator> B(A);
  
  A = B;
  ScaleMatrix(A, Drow, Dcol);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Drow(i)*Dcol(j)*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(Drow(i)*Dcol(j)*B(i,j)))
        {
          cout << "Scale Matrix incorrect" << endl;
          abort();
        }
  
  A = B;
  ScaleLeftMatrix(A, Drow);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Drow(i)*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(Drow(i)*B(i,j)))
        {
          cout << "ScaleLeft Matrix incorrect" << endl;
          abort();
        }

  A = B;
  ScaleRightMatrix(A, Dcol);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Dcol(j)*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(Dcol(j)*B(i,j)))
        {
          cout << "ScaleRight Matrix incorrect" << endl;
          abort();
        }

  // testing Mlt and Add
  T alpha;
  GetRandNumber(alpha);
  A = B;
  Mlt(alpha, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(alpha*B(i,j)))
        {
          cout << "Mlt incorrect" << endl;
          abort();
        }
  
  Real_wp beta;
  GetRandNumber(beta);
  A = B;
  Mlt(beta, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - beta*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(beta*B(i,j)))
        {
          cout << "Mlt incorrect" << endl;
          abort();
        }
  
  B.FillRand();
  Mlt(Real_wp(1)/RAND_MAX, B);
  Matrix<T, Prop, Storage, Allocator> C;
  C = A;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  A = C;
  Add(beta, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - beta*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(beta*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }
  
  GenerateRandomTriangular(B, m, n, low);
  
  A = C;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  Real_wp norm_A = MaxAbs(A);
  Real_wp val_max = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (abs(A(i, j)) > val_max)
        val_max = abs(A(i, j));
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "MaxAbs incorrect" << endl;
      abort();
    }

  val_max = 0;
  norm_A = NormInf(A);
  for (int i = 0; i < m; i++)
    {
      Real_wp sum(0);
      for (int j = 0; j < n; j++)
        sum += abs(A(i, j));
      
      if (sum > val_max)
        val_max = sum;
    }
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "NormInf incorrect" << endl;
      abort();
    }

  val_max = 0;
  norm_A = Norm1(A);
  for (int j = 0; j < n; j++)
    {
      Real_wp sum(0);
      for (int i = 0; i < m; i++)
        sum += abs(A(i, j));
      
      if (sum > val_max)
        val_max = sum;
    }
  
  if ((abs(val_max - norm_A) > threshold) || isnan(val_max) || isnan(norm_A))
    {
      cout << "Norm1 incorrect" << endl;
      abort();
    }
}

int main(int argc, char** argv)
{
  threshold = 1e-12;
  
  //srand(time(NULL));

  {
    Matrix<Real_wp, General, RowMajor> A;
    CheckGeneralMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Real_wp, General, ColMajor> A;
    CheckGeneralMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Real_wp, Symmetric, RowSymPacked> A;
    CheckSymmetricMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSymPacked> A;
    CheckSymmetricMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  } 
#endif
  
 {
    Matrix<Real_wp, Symmetric, RowSym> A;
    CheckSymmetricMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSym> A;
    CheckSymmetricMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
#endif

 {
    Matrix<Complex_wp, Hermitian, RowHermPacked> A;
    CheckHermitianMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
 {
    Matrix<Complex_wp, Hermitian, ColHermPacked> A;
    CheckHermitianMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
#endif

 {
    Matrix<Complex_wp, Hermitian, RowHerm> A;
    CheckHermitianMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
 {
    Matrix<Complex_wp, Hermitian, ColHerm> A;
    CheckHermitianMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
#endif
  
  {
    Matrix<Real_wp, General, RowSparse> A;
    CheckGeneralMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, RowSparse> A;
    CheckGeneralMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Real_wp, General, ColSparse> A;
    CheckGeneralMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, ColSparse> A;
    CheckGeneralMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckGeneralMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    CheckGeneralMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
  
#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ArrayColSparse> A;
    CheckGeneralMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, ArrayColSparse> A;
    CheckGeneralMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
#endif

  {
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    CheckSymmetricMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymSparse> A;
    CheckSymmetricMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSymSparse> A;
    CheckSymmetricMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, Symmetric, ColSymSparse> A;
    CheckSymmetricMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
#endif
  
  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckSymmetricMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    CheckSymmetricMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A;
    CheckSymmetricMatrix(A);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayColSymSparse> A;
    CheckSymmetricMatrix(A);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (!IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
#endif
  
  {
    Matrix<Real_wp, General, RowLoTriangPacked> A;
    CheckTriangularMatrix(A, true);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, RowLoTriangPacked> A;
    CheckTriangularMatrix(A, true);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Real_wp, General, RowLoTriang> A;
    CheckTriangularMatrix(A, true);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, RowLoTriang> A;
    CheckTriangularMatrix(A, true);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColLoTriangPacked> A;
    CheckTriangularMatrix(A, true);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, ColLoTriangPacked> A;
    CheckTriangularMatrix(A, true);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Real_wp, General, ColLoTriang> A;
    CheckTriangularMatrix(A, true);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, ColLoTriang> A;
    CheckTriangularMatrix(A, true);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
#endif

  {
    Matrix<Real_wp, General, RowUpTriangPacked> A;
    CheckTriangularMatrix(A, false);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, RowUpTriangPacked> A;
    CheckTriangularMatrix(A, false);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Real_wp, General, RowUpTriang> A;
    CheckTriangularMatrix(A, false);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, RowUpTriang> A;
    CheckTriangularMatrix(A, false);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ColUpTriangPacked> A;
    CheckTriangularMatrix(A, false);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, ColUpTriangPacked> A;
    CheckTriangularMatrix(A, false);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Real_wp, General, ColUpTriang> A;
    CheckTriangularMatrix(A, false);
    if (IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, ColUpTriang> A;
    CheckTriangularMatrix(A, false);
    if (!IsComplexMatrix(A))
      {
        cout << "IsComplexMatrix incorrect" << endl;
        abort();
      }

    if (IsSymmetricMatrix(A))
      {
        cout << "IsSymmetricMatrix incorrect" << endl;
        abort();
      }
  }
#endif
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
