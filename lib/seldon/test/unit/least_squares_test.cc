#include "SeldonLib.hxx"

using namespace Seldon;

typedef Vector<double> VectReal_wp;

Real_wp threshold = 1e-12;

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

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y,
		 Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i) - y(i)) > eps)
      return false;
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
			      int m, int n)
{
  T one, zero;
  SetComplexOne(one);
  SetComplexZero(zero);

  A.Reallocate(m, n);
  A.Fill(zero);
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	GetRandNumber(x);
	A.Set(i, j, x);
      }
}

template<class T, class Prop, class Storage, class Alloc>
bool CheckIdentity(const Matrix<T, Prop, Storage, Alloc>& A)
{
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == j)
          {
            if (abs(A(i, j) - 1.0) > threshold)
              return false;
          }
        else
          {
            if (abs(A(i, j)) > threshold)
              return false;
          }
      }
  
  return true;
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
      if (abs(A(i, j) - B(i, j)) > threshold)
        return false;
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void CheckGeneralMatrix(Matrix<T, Prop, Storage, Allocator>& A)
{
  T one, zero;
  SetComplexOne(one);
  SetComplexZero(zero);

  Matrix<T, Prop, Storage, Allocator> A0;
  
  ////////////////
  // First case : m > n
  //////////////////////
  
  int m = 25, n = 14;
  GenerateRandomMatrix(A, m, n);
  A0 = A;
  
  Matrix<T, General, Storage> R;
  Matrix<T, General, RowUpTriangPacked> Rt(n, n);
  Matrix<T, Prop, Storage, Allocator> Q, B(n, n);
  // we get Q, R as in a Matlab script
  GetQR(A, Q, R);
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < i; j++)
	R(i, j) = zero;
      
      for (int j = i; j < n; j++)
	Rt.Val(i, j) = R(i, j);
    }
  
  // checking orthogonality of Q
  B.Fill(zero);  
  MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, Q, zero, B);
  if (!CheckIdentity(B))
    {
      cout << "Q not orthogonal" << endl;
      abort();
    }
  
  // then checking that A = Q*R
  B.Reallocate(m, n);
  MltAdd(one, Q, R, zero, B);
  if (!EqualMatrix(B, A))
    {
      cout <<"QR factorisation incorrect" << endl;
      abort();
    }
  
  // testing GetQR/SolveQR
  A = A0;
  Vector<T> tau, b(m), x(n), y;
  GetQR(A, tau);
  
  GenerateRandomVector(b, m);
  y = b;
  SolveQR(A, tau, y);
  
  if (IsComplexMatrix(Q))
    Mlt(SeldonConjTrans, Q, b, x);
  else
    Mlt(SeldonTrans, Q, b, x);
  
  Solve(Rt, x);
  
  if (!EqualVector(x, y))
    {
      cout << "GetQR/SolveQR incorrect" << endl;
      abort();
    }
  
  // testing MltQ_FromQR
  GenerateRandomVector(x, n);
  b.Fill(zero);
  for (int i = 0; i < x.GetM(); i++)
    b(i) = x(i);
    
  MltQ_FromQR(SeldonNoTrans, A, tau, b);
  y.Reallocate(m);
  Mlt(SeldonNoTrans, Q, x, y);
  
  if (!EqualVector(b, y))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }
  
  GenerateRandomVector(x, m);
  b.Fill(zero);
  for (int i = 0; i < x.GetM(); i++)
    b(i) = x(i);
  
  if (IsComplexMatrix(Q))
    MltQ_FromQR(SeldonConjTrans, A, tau, b);
  else
    MltQ_FromQR(SeldonTrans, A, tau, b);
  
  y.Reallocate(n);
  if (IsComplexMatrix(Q))
    Mlt(SeldonConjTrans, Q, x, y);
  else
    Mlt(SeldonTrans, Q, x, y);
  
  b.Resize(n);
  if (!EqualVector(b, y))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }


  // testing MltQ_FromQR with a matrix
  Matrix<T, Prop, Storage, Allocator> C, D, C0;
  int p = 5;
  GenerateRandomMatrix(C, n, p);
  C0 = C;
  C.Resize(m, p);
  for (int i = n; i < m; i++)
    for (int j = 0; j < p; j++)
      C(i, j) = zero;
  
  MltQ_FromQR(SeldonLeft, SeldonNoTrans, A, tau, C);
  D.Reallocate(m, p);
  MltAdd(one, SeldonNoTrans, Q, SeldonNoTrans, C0, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQR incorrect" << endl;
      abort();
    }
  
  GenerateRandomMatrix(C, m, p);
  C0 = C;
  
  if (IsComplexMatrix(Q))
    MltQ_FromQR(SeldonLeft, SeldonConjTrans, A, tau, C);
  else
    MltQ_FromQR(SeldonLeft, SeldonTrans, A, tau, C);
  
  D.Reallocate(n, p);
  if (IsComplexMatrix(Q))
    MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C0, zero, D);
  else
    MltAdd(one, SeldonTrans, Q, SeldonNoTrans, C0, zero, D);
  
  C.Resize(n, p);
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }
    
  GenerateRandomMatrix(C, p, m);
  C0 = C;
  
  MltQ_FromQR(SeldonRight, SeldonNoTrans, A, tau, C);
  D.Reallocate(p, n);
  MltAdd(one, SeldonNoTrans, C0, SeldonNoTrans, Q, zero, D);
  
  C.Resize(p, n);
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }

  GenerateRandomMatrix(C, p, n);
  C0 = C;
  C.Resize(p, m);
  for (int i = n; i < m; i++)
    for (int j = 0; j < p; j++)
      C(j, i) = zero;
  
  if (IsComplexMatrix(Q))
    MltQ_FromQR(SeldonRight, SeldonConjTrans, A, tau, C);
  else
    MltQ_FromQR(SeldonRight, SeldonTrans, A, tau, C);
  
  D.Reallocate(p, m);
  if (IsComplexMatrix(Q))
    MltAdd(one, SeldonNoTrans, C0, SeldonConjTrans, Q, zero, D);
  else
    MltAdd(one, SeldonNoTrans, C0, SeldonTrans, Q, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }

  //////////////////
  // Second case : m < n
  ///////////////////////
  
  m = 15; n = 18;
  GenerateRandomMatrix(A, m, n);
  Rt.Reallocate(m, m);
  A0 = A;
  
  // extracting Q and R as in Matlab
  B.Reallocate(m, m);
  GetQR(A, Q, R);
  for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < i; j++)
	R(i, j) = zero;
      
      for (int j = i; j < m; j++)
	Rt.Val(i, j) = R(i, j);
    }
  
  // checking orthogonality of Q
  B.Fill(zero);
  MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, Q, zero, B);
  if (!CheckIdentity(B))
    {
      cout << "Q not orthogonal" << endl;
      abort();
    }
  
  // checking that Q*R = A
  B.Reallocate(m, n);
  MltAdd(one, Q, R, zero, B);
  if (!EqualMatrix(B, A))
    {
      cout <<"QR factorisation incorrect" << endl;
      abort();
    }
  
  // testing GetQR/SolveQR
  A = A0;
  b.Reallocate(m); x.Reallocate(m);
  GetQR(A, tau);
  
  GenerateRandomVector(b, m);
  y = b;
  SolveQR(A, tau, y);
  
  if (IsComplexMatrix(Q))
    Mlt(SeldonConjTrans, Q, b, x);
  else
    Mlt(SeldonConjTrans, Q, b, x);
  
  Solve(Rt, x);
  
  if (!EqualVector(x, y))
    {
      cout << "GetQR/SolveQR incorrect" << endl;
      abort();
    }
  
  // testing MltQ_FromQR
  GenerateRandomVector(x, m);
  b.Reallocate(m);
  for (int i = 0; i < m; i++)
    b(i) = x(i);
    
  MltQ_FromQR(SeldonNoTrans, A, tau, b);
  y.Reallocate(m);
  Mlt(SeldonNoTrans, Q, x, y);
  
  if (!EqualVector(b, y))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }
  
  GenerateRandomVector(x, m);
  b.Fill(zero);
  for (int i = 0; i < x.GetM(); i++)
    b(i) = x(i);
  
  if (IsComplexMatrix(Q))
    MltQ_FromQR(SeldonConjTrans, A, tau, b);
  else
    MltQ_FromQR(SeldonTrans, A, tau, b);
  
  y.Reallocate(m);
  if (IsComplexMatrix(Q))
    Mlt(SeldonConjTrans, Q, x, y);
  else
    Mlt(SeldonTrans, Q, x, y);
  
  if (!EqualVector(b, y))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }


  // testing MltQ_FromQR with a matrix
  GenerateRandomMatrix(C, m, p);
  C0 = C;
  
  MltQ_FromQR(SeldonLeft, SeldonNoTrans, A, tau, C);
  D.Reallocate(m, p);
  MltAdd(one, SeldonNoTrans, Q, SeldonNoTrans, C0, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }
  
  GenerateRandomMatrix(C, m, p);
  C0 = C;
  
  if (IsComplexMatrix(Q))
    MltQ_FromQR(SeldonLeft, SeldonConjTrans, A, tau, C);
  else
    MltQ_FromQR(SeldonLeft, SeldonTrans, A, tau, C);
  
  D.Reallocate(m, p);
  if (IsComplexMatrix(Q))
    MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C0, zero, D);
  else
    MltAdd(one, SeldonTrans, Q, SeldonNoTrans, C0, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }
    
  GenerateRandomMatrix(C, p, m);
  C0 = C;
  
  MltQ_FromQR(SeldonRight, SeldonNoTrans, A, tau, C);
  D.Reallocate(p, m);
  MltAdd(one, SeldonNoTrans, C0, SeldonNoTrans, Q, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }

  GenerateRandomMatrix(C, p, m);
  C0 = C;
  
  if (IsComplexMatrix(Q))
    MltQ_FromQR(SeldonRight, SeldonConjTrans, A, tau, C);
  else
    MltQ_FromQR(SeldonRight, SeldonTrans, A, tau, C);
  
  D.Reallocate(p, m);
  if (IsComplexMatrix(Q))
    MltAdd(one, SeldonNoTrans, C0, SeldonConjTrans, Q, zero, D);
  else
    MltAdd(one, SeldonNoTrans, C0, SeldonTrans, Q, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }
  
  //////////////////
  // LQ factorisation with m > n
  ///////////////////////
  
  m = 27; n = 21;
  GenerateRandomMatrix(A, m, n);
  Matrix<T, General, RowLoTriangPacked> Rl(n, n);
  A0 = A;
  
  // extracting L and Q as in Matlab
  Matrix<T, General, Storage> L;
  B.Reallocate(n, n); Q.Clear();
  GetLQ(A, L, Q);
  for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < j; i++)
	L(i, j) = zero;
      
      for (int i = j; i < n; i++)
	Rl.Val(i, j) = L(i, j);
    }
  
  // checking orthogonality of Q
  B.Fill(zero);
  MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, Q, zero, B);
  if (!CheckIdentity(B))
    {
      cout << "Q not orthogonal" << endl;
      abort();
    }
  
  // checking that L*Q = A
  B.Reallocate(m, n);
  MltAdd(one, L, Q, zero, B);
  if (!EqualMatrix(B, A))
    {
      cout <<"LQ factorisation incorrect" << endl;
      abort();
    }
  
  // testing GetLQ
  A = A0;
  b.Reallocate(n); x.Reallocate(n);
  GetLQ(A, tau);
  
  // testing MltQ_FromLQ
  GenerateRandomVector(x, n);
  b.Reallocate(n);
  for (int i = 0; i < n; i++)
    b(i) = x(i);
    
  MltQ_FromLQ(SeldonNoTrans, A, tau, b);
  y.Reallocate(n);
  Mlt(SeldonNoTrans, Q, x, y);
  
  if (!EqualVector(b, y))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }
  
  GenerateRandomVector(x, n);
  b.Fill(zero);
  for (int i = 0; i < x.GetM(); i++)
    b(i) = x(i);
  
  if (IsComplexMatrix(Q))
    MltQ_FromLQ(SeldonConjTrans, A, tau, b);
  else
    MltQ_FromLQ(SeldonTrans, A, tau, b);
  
  y.Reallocate(n);
  if (IsComplexMatrix(Q))
    Mlt(SeldonConjTrans, Q, x, y);
  else
    Mlt(SeldonTrans, Q, x, y);
  
  if (!EqualVector(b, y))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }

  // testing MltQ_FromLQ with a matrix
  GenerateRandomMatrix(C, n, p);
  C0 = C;
  
  MltQ_FromLQ(SeldonLeft, SeldonNoTrans, A, tau, C);
  D.Reallocate(n, p);
  MltAdd(one, SeldonNoTrans, Q, SeldonNoTrans, C0, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }
  
  GenerateRandomMatrix(C, n, p);
  C0 = C;
  
  if (IsComplexMatrix(Q))
    MltQ_FromLQ(SeldonLeft, SeldonConjTrans, A, tau, C);
  else
    MltQ_FromLQ(SeldonLeft, SeldonTrans, A, tau, C);
  
  D.Reallocate(n, p);
  if (IsComplexMatrix(Q))
    MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C0, zero, D);
  else
    MltAdd(one, SeldonTrans, Q, SeldonNoTrans, C0, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }
    
  GenerateRandomMatrix(C, p, n);
  C0 = C;
  
  MltQ_FromLQ(SeldonRight, SeldonNoTrans, A, tau, C);
  D.Reallocate(p, n);
  MltAdd(one, SeldonNoTrans, C0, SeldonNoTrans, Q, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }

  GenerateRandomMatrix(C, p, n);
  C0 = C;
  
  if (IsComplexMatrix(Q))
    MltQ_FromLQ(SeldonRight, SeldonConjTrans, A, tau, C);
  else
    MltQ_FromLQ(SeldonRight, SeldonTrans, A, tau, C);
  
  D.Reallocate(p, n);
  if (IsComplexMatrix(Q))
    MltAdd(one, SeldonNoTrans, C0, SeldonConjTrans, Q, zero, D);
  else
    MltAdd(one, SeldonNoTrans, C0, SeldonTrans, Q, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }
  
  
  //////////////
  // Last case : m < n (LQ)
  /////////////
  
  m = 12; n = 16;
  GenerateRandomMatrix(A, m, n);
  Rl.Reallocate(m, m);
  A0 = A;
  
  // extracting L and Q as in Matlab
  B.Reallocate(m, m); Q.Clear();
  GetLQ(A, L, Q);
  for (int j = 0; j < m; j++)
    {
      for (int i = 0; i < j; i++)
	L(i, j) = zero;
      
      for (int i = j; i < m; i++)
	Rl.Val(i, j) = L(i, j);
    }
  
  // checking orthogonality of Q
  B.Fill(zero);
  MltAdd(one, SeldonNoTrans, Q, SeldonConjTrans, Q, zero, B);
  if (!CheckIdentity(B))
    {
      cout << "Q not orthogonal" << endl;
      abort();
    }
  
  // checking that L*Q = A
  B.Reallocate(m, n);
  MltAdd(one, L, Q, zero, B);
  if (!EqualMatrix(B, A))
    {
      cout <<"LQ factorisation incorrect" << endl;
      abort();
    }

  // testing GetLQ/SolveLQ
  A = A0;
  b.Reallocate(n); x.Reallocate(n);
  GetLQ(A, tau);
  
  GenerateRandomVector(b, m);
  y = b;
  y.Resize(n);
  SolveLQ(A, tau, y);
  
  Solve(Rl, b);
  if (IsComplexMatrix(Q))
    Mlt(SeldonConjTrans, Q, b, x);
  else
    Mlt(SeldonTrans, Q, b, x);
    
  if (!EqualVector(x, y))
    {
      cout << "GetQR/SolveQR incorrect" << endl;
      abort();
    }

  // testing MltQ_FromLQ
  GenerateRandomVector(x, n);
  b.Reallocate(n);
  for (int i = 0; i < n; i++)
    b(i) = x(i);
    
  MltQ_FromLQ(SeldonNoTrans, A, tau, b);
  y.Reallocate(m);
  Mlt(SeldonNoTrans, Q, x, y);
  
  b.Resize(m);
  if (!EqualVector(b, y))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }
  
  GenerateRandomVector(x, m);
  b.Reallocate(n);
  b.Fill(zero);
  for (int i = 0; i < x.GetM(); i++)
    b(i) = x(i);
  
  if (IsComplexMatrix(Q))
    MltQ_FromLQ(SeldonConjTrans, A, tau, b);
  else
    MltQ_FromLQ(SeldonTrans, A, tau, b);
  
  y.Reallocate(n);
  if (IsComplexMatrix(Q))
    Mlt(SeldonConjTrans, Q, x, y);
  else
    Mlt(SeldonTrans, Q, x, y);
  
  if (!EqualVector(b, y))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }
  
  // testing MltQ_FromLQ with a matrix
  GenerateRandomMatrix(C, n, p);
  C0 = C;
  
  MltQ_FromLQ(SeldonLeft, SeldonNoTrans, A, tau, C);
  D.Reallocate(m, p);
  MltAdd(one, SeldonNoTrans, Q, SeldonNoTrans, C0, zero, D);
  
  C.Resize(m, p);
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }
  
  GenerateRandomMatrix(C, m, p);
  C0 = C;
  C.Resize(n, p);
  for (int i = m; i < n; i++)
    for (int j = 0; j < p; j++)
      C(i, j) = zero;
  
  if (IsComplexMatrix(Q))
    MltQ_FromLQ(SeldonLeft, SeldonConjTrans, A, tau, C);
  else
    MltQ_FromLQ(SeldonLeft, SeldonTrans, A, tau, C);
  
  D.Reallocate(n, p);
  if (IsComplexMatrix(Q))
    MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C0, zero, D);
  else
    MltAdd(one, SeldonTrans, Q, SeldonNoTrans, C0, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }
    
  GenerateRandomMatrix(C, p, m);
  C0 = C;
  C.Resize(p, n);
  for (int i = m; i < n; i++)
    for (int j = 0; j < p; j++)
      C(j, i) = zero;
  
  MltQ_FromLQ(SeldonRight, SeldonNoTrans, A, tau, C);
  D.Reallocate(p, n);
  MltAdd(one, SeldonNoTrans, C0, SeldonNoTrans, Q, zero, D);
  
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromLQ incorrect" << endl;
      abort();
    }

  GenerateRandomMatrix(C, p, n);
  C0 = C;
  
  if (IsComplexMatrix(Q))
    MltQ_FromLQ(SeldonRight, SeldonConjTrans, A, tau, C);
  else
    MltQ_FromLQ(SeldonRight, SeldonTrans, A, tau, C);
  
  D.Reallocate(p, m);
  if (IsComplexMatrix(Q))
    MltAdd(one, SeldonNoTrans, C0, SeldonConjTrans, Q, zero, D);
  else
    MltAdd(one, SeldonNoTrans, C0, SeldonTrans, Q, zero, D);
  
  C.Resize(p, m);
  if (!EqualMatrix(C, D))
    {
      cout << "MltQ_FromQr incorrect" << endl;
      abort();
    }
}

int main(int argc, char** argv)
{
  threshold = 1e-12;

  // srand(time(NULL));
  
  {
    Matrix<Real_wp, General, ColMajor> A;
    CheckGeneralMatrix(A);
  }

  {
    Matrix<Complex_wp, General, ColMajor> A;
    CheckGeneralMatrix(A);
  }
  
  {
    Matrix<Real_wp, General, RowMajor> A;
    CheckGeneralMatrix(A);
  }
  
  {
    Matrix<Complex_wp, General, RowMajor> A;
    CheckGeneralMatrix(A);
  }
  
  cout << "All tests passed successfully" << endl;

  return 0;
}
