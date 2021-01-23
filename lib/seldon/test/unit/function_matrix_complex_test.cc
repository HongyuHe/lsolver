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

template<class T, class Prop, class Storage, class Allocator,
         class Prop2, class Storage2, class Allocator2>
void CheckGeneralMatrix(Matrix<complex<T>, Prop, Storage, Allocator>& A,
                        Matrix<complex<T>, Prop2, Storage2, Allocator2>& Bcplx)
{
  int m = 70, n = 36, nnz = 200;

  GenerateRandomMatrix(A, m, n, nnz);
  
  IVect col_perm;
  GenerateRandomPermutation(n, col_perm);

  IVect row_perm;
  GenerateRandomPermutation(m, row_perm);
  
  Matrix<complex<T>, Prop, Storage, Allocator> B(A);
  ApplyInversePermutation(A, row_perm, col_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((B(i, j) != A(row_perm(i), col_perm(j))) || isnan(B(i, j)) || isnan(A(row_perm(i), col_perm(j))))
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
      if ((A(i, j) != B(row_perm(i), col_perm(j))) || isnan(A(i, j)) || isnan(B(row_perm(i), col_perm(j))))
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
      if ((abs(A(i, j) - Dcol(j)*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(Dcol(j)*B(i, j)))
        {
          cout << "ScaleRight Matrix incorrect" << endl;
          abort();
        }

  // testing Mlt and Add
  T alpha;
  GetRandNumber(alpha);
  A = B;
  Mlt(complex<T>(alpha, 0), A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (abs(A(i, j) - alpha*B(i,j)) > threshold)
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
  Matrix<complex<T>, Prop, Storage, Allocator> C;
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
  
  GenerateRandomMatrix(B, m, n, nnz+23);
  
  A = C;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i, j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  A = C;  
  GenerateRandomMatrix(Bcplx, m, n, nnz+12);
  Matrix<complex<T>, Prop2, Storage2> Bcplx2(Bcplx);
  Add(alpha, A, Bcplx);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(Bcplx(i, j) - Bcplx2(i, j) - alpha*A(i, j)) > threshold)
          || isnan(Bcplx(i, j)) || isnan(Bcplx2(i, j)) || isnan(alpha*A(i, j)))
        {
          DISP(Bcplx(i, j)); DISP(Bcplx2(i, j)); DISP(alpha*A(i, j));
          DISP(abs(Bcplx(i, j) - Bcplx2(i, j) - alpha*A(i, j)));
          cout << "Add incorrect" << endl;
          abort();
        }
  
  Matrix<Real_wp, Prop2, Storage2> Breal, Bimag;
  GenerateRandomMatrix(Breal, m, n, nnz+21);
  GenerateRandomMatrix(Bimag, m, n, nnz-39);
  AddMatrix(alpha, Breal, Bimag, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - A(i, j) - alpha*complex<T>(Breal(i,j), Bimag(i, j))) > threshold)
          || isnan(C(i, j)) || isnan(A(i, j)) || isnan(alpha*complex<T>(Breal(i,j), Bimag(i, j))))
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
  
  if ((abs(val_max - norm_A) > threshold)  || isnan(val_max) || isnan(norm_A))
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
      if ((A(i, j) != conj(B(i, j)) ) || isnan(A(i, j)) || isnan(conj(B(i, j))))
        {
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
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

template<class T, class Prop, class Storage, class Allocator,
         class Prop2, class Storage2, class Allocator2>
void CheckSymmetricMatrix(Matrix<complex<T>, Prop, Storage, Allocator>& A,
                          Matrix<complex<T>, Prop2, Storage2, Allocator2>& Bcplx)
{
  int m = 48, n = 48, nnz = 200;
  
  GenerateRandomMatrix(A, m, n, nnz);
  
  IVect row_perm;
  GenerateRandomPermutation(m, row_perm);
  
  Matrix<complex<T>, Prop, Storage, Allocator> B(A);
  ApplyInversePermutation(A, row_perm, row_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((B(i, j) != A(row_perm(i), row_perm(j))) || isnan(B(i, j)) || isnan(A(row_perm(i), row_perm(j))) )
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
      if ((A(i, j) != B(row_perm(i), row_perm(j))) || isnan(A(i, j)) || isnan(B(row_perm(i), row_perm(j))))
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
      if ((abs(A(i, j) - Drow(i)*Drow(j)*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(Drow(i)*Drow(j)*B(i,j)))
        {
          cout << "Scale Matrix incorrect" << endl;
          abort();
        }

  // testing Mlt and Add
  T alpha;
  GetRandNumber(alpha);
  A = B;
  Mlt(complex<T>(alpha, 0), A);
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
  Matrix<complex<T>, Prop, Storage, Allocator> C;
  C = A;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i,j)))
        {
          DISP(A(i, j)); DISP(C(i, j)); DISP(alpha); DISP(B(i, j));
          DISP(A(i, j) - C(i, j) - alpha*B(i,j));
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
  
  GenerateRandomMatrix(B, m, n, nnz+23);
  
  A = C;
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  A = C;  
  GenerateRandomMatrix(Bcplx, m, n, nnz+12);
  Matrix<complex<T>, Prop2, Storage2, Allocator2> Bcplx2(Bcplx);
  Add(alpha, A, Bcplx);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(Bcplx(i, j) - Bcplx2(i, j) - alpha*A(i,j)) > threshold)
          || isnan(Bcplx(i, j)) || isnan(Bcplx2(i, j)) || isnan(alpha*A(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  Matrix<T, Prop2, Storage2> Breal, Bimag;
  GenerateRandomMatrix(Breal, m, n, nnz+21);
  GenerateRandomMatrix(Bimag, m, n, nnz-39);
  AddMatrix(alpha, Breal, Bimag, C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) - A(i, j) - alpha*complex<T>(Breal(i,j), Bimag(i, j))) > threshold)
          || isnan(C(i, j)) || isnan(A(i, j)) || isnan(alpha*complex<T>(Breal(i,j), Bimag(i, j))))
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
      if ((A(i, j) != B(i, j) ) || isnan(A(i, j)) || isnan(B(i, j) ))
        {
          cout << "Transpose incorrect" << endl;
          abort();
        }
  
  Conjugate(A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if ((A(i, j) != conj(B(i, j)) ) || isnan(A(i, j)) || isnan(conj(B(i, j))))
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
void CheckGeneralMatrix(Matrix<complex<T>, Prop, Storage, Allocator>& A)
{
  int m = 70, n = 36, nnz = 200;

  GenerateRandomMatrix(A, m, n, nnz);
  
  IVect col_perm;
  GenerateRandomPermutation(n, col_perm);

  IVect row_perm;
  GenerateRandomPermutation(m, row_perm);
  
  Matrix<complex<T>, Prop, Storage, Allocator> B(A);
  /* ApplyInversePermutation(A, row_perm, col_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (B(i, j) != A(row_perm(i), col_perm(j)))
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
      if (A(i, j) != B(row_perm(i), col_perm(j)))
        {
          DISP(i); DISP(j);
          DISP(A(i, j)); DISP(B(row_perm(i), col_perm(j)));
          cout << "ApplyPermutation incorrect" << endl;
          abort();
        }
  */
  VectReal_wp Drow, Dcol;
  GenerateRandomVector(Drow, m);
  GenerateRandomVector(Dcol, n);
  
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
  Mlt(complex<T>(alpha, 0), A);
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
  
  GenerateRandomMatrix(B, m, n, nnz+23);
  
  Matrix<complex<T>, Prop, Storage, Allocator> C;
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
  
  A = C;  
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
      DISP(val_max); DISP(norm_A); DISP(val_max-norm_A);
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
      if (A(i, j) != B(i, j) )
        {
          cout << "Transpose incorrect" << endl;
          abort();
        }
  
  Conjugate(A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if ((A(i, j) != conj(B(i, j)) ) || isnan(A(i, j)) || isnan(conj(B(i, j))) )
        {
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
          cout << "Conjugate incorrect" << endl;
          abort();
        }
  
  TransposeConj(C);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((A(i, j) != C(i, j) ) || isnan(A(i, j)) || isnan(C(i, j)) )
        {
          cout << "TransposeConj incorrect" << endl;
          abort();
        }
}

template<class T, class Prop, class Storage, class Allocator>
void CheckSymmetricMatrix(Matrix<complex<T>, Prop, Storage, Allocator>& A)
{
  int m = 48, n = 48, nnz = 200;
  
  GenerateRandomMatrix(A, m, n, nnz);
  
  IVect row_perm;
  GenerateRandomPermutation(m, row_perm);
  
  Matrix<complex<T>, Prop, Storage, Allocator> B(A);
  /* ApplyInversePermutation(A, row_perm, row_perm);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (B(i, j) != A(row_perm(i), row_perm(j)))
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
      if (A(i, j) != B(row_perm(i), row_perm(j)))
        {
          DISP(i); DISP(j);
          DISP(A(i, j)); DISP(B(row_perm(i), row_perm(j)));
          cout << "ApplyPermutation incorrect" << endl;
          abort();
        }
  */
  VectReal_wp Drow;
  GenerateRandomVector(Drow, m);
  
  A = B;
  ScaleMatrix(A, Drow, Drow);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - Drow(i)*Drow(j)*B(i,j)) > threshold)
          || isnan(A(i, j) ) || isnan(Drow(i)*Drow(j)*B(i,j)))
        {
          cout << "Scale Matrix incorrect" << endl;
          abort();
        }

  // testing Mlt and Add
  T alpha;
  GetRandNumber(alpha);
  A = B;
  Mlt(complex<T>(alpha, 0), A);
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
  
  GenerateRandomMatrix(B, m, n, nnz+23);
  
  Matrix<complex<T>, Prop, Storage, Allocator> C(A);
  Add(alpha, B, A);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - C(i, j) - alpha*B(i,j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(alpha*B(i,j)))
        {
          cout << "Add incorrect" << endl;
          abort();
        }

  A = C;  
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
      if ((A(i, j) != conj(B(i, j)) ) || isnan(A(i, j)) || isnan(conj(B(i, j))))
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

int main(int argc, char** argv)
{
  threshold = 1e-12;
  
  //srand(time(NULL));

  {
    Matrix<Complex_wp, General, RowComplexSparse> A;
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
    
    CheckGeneralMatrix(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, General, ColComplexSparse> A;
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
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    Matrix<Complex_wp, General, ArrayRowSparse> B;
    CheckGeneralMatrix(A, B);
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
    Matrix<Complex_wp, General, ArrayColComplexSparse> A;
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
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
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

    CheckSymmetricMatrix(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A;
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
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> B;
    CheckSymmetricMatrix(A, B);
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
    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A;
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

  cout << "All tests passed successfully" << endl;
  
  return 0;
}
