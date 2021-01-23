#include "SeldonLib.hxx"

using namespace Seldon;

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

void CheckPermutation(const IVect& num2)
{
  IVect num(num2);
  Sort(num);
  for (int i = 0; i < num.GetM(); i++)
    if (num(i) != i)
      {
	cout << " Not a permutation array" << endl;
	abort();
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
    if (abs(x(i) - y(i)) > eps)
      return false;
  
  return true;
}


template<class T, class T2>
bool EqualMatrix(const Matrix<T, General, ColMajor>& x,
		 const Matrix<T2, General, ColMajor>& y,
		 Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;

  if (x.GetN() != y.GetN())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    for (int j = 0; j < x.GetN(); j++)
      if ( (abs(x(i, j) - y(i, j)) > eps) || isnan(abs(x(i, j) - y(i, j))) )
	return false;
  
  return true;
}

template<class MatrixSparse, class MatrixLU>
void CheckSolveComplexVector(const MatrixSparse& A, MatrixLU& mat_lu, Vector<Complex_wp>& x)
{
}

template<class MatrixSparse, class MatrixLU>
void CheckSolveComplexVector(const MatrixSparse& A, MatrixLU& mat_lu, Vector<Real_wp>& x)
{
  int n = A.GetM();
  Vector<Complex_wp> xc(n), bc(n), yc;
  
  GenerateRandomVector(xc, n);  
  
  Mlt(A, xc, bc);  
  yc = xc;
  
  xc = bc;
  SolveLU(mat_lu, xc);
  
  if (!EqualVector(xc, yc))
    {
      cout << "Solves by a complex vector incorrect" << endl;
      abort();
    }
}


template<class T, class Prop, class Storage, class Allocator>
void CheckDirectSolver(Matrix<T, Prop, Storage, Allocator>& A)
{
  T zero, one;
  SetComplexZero(zero);
  SetComplexOne(one);
  
  int n = 400, nnz = 5000;
  //int n = 7, nnz = 20;
  GenerateRandomMatrix(A, n, n, nnz);

  for (int i = 0; i < n; i++)
    {
      Real_wp sum = 0;
      for (int j = 0; j < n; j++)
	sum += abs(A(i, j));
      
      A.Set(i, i, 1.0 + sum);
    }
  
  Vector<T> x, y, b(n), bt(n), z(n);

  int nb_rhs = 20;
  Matrix<T, General, ColMajor> x_mat(n, nb_rhs), b_mat(n, nb_rhs), bt_mat(n, nb_rhs);
  
  for (int i = 0; i < nb_rhs; i ++)
    {
      GenerateRandomVector(x, n);
      Mlt(A, x, b);
      Mlt(SeldonTrans, A, x, bt);
      SetCol(b, i, b_mat);
      SetCol(bt, i, bt_mat);
      SetCol(x, i, x_mat);
    }

  Matrix<T, General, ColMajor> y_mat(x_mat);
  
  IVect num(n);
  GenerateRandomVector(x, n);
  GenerateRandomPermutation(n, num);

  z.Fill(zero);
  b.Fill(zero); bt.Fill(zero);
  Mlt(A, x, b);
  Mlt(SeldonTrans, A, x, bt);
  y = x;

  {
#ifdef SELDON_WITH_WSMP
    MatrixWsmp<T> mat_lu;
    mat_lu.ShowMessages();
    mat_lu.HideMessages();
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "FactorizeMatrix of Wsmp incorrect" << endl;
	abort();
      }
    
    mat_lu.HideMessages();
    
    x = b;
    mat_lu.Solve(x);
    if (!EqualVector(x, y))
      {
	cout << "Solve of Wsmp incorrect " << endl;
        DISP(x); DISP(y);
	abort();
      }    
    
    x = b;
    mat_lu.Solve(SeldonNoTrans, x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of Wsmp incorrect " << endl;
	abort();
      }

    x = bt;
    mat_lu.Solve(SeldonTrans, x);

    if (!EqualVector(x, y))
      {
        DISP(x); DISP(y);
	cout << "Solve of Wsmp incorrect " << endl;
	abort();
      }

    // testing resolution of multiple right hand-sides
    Matrix<T, General, ColMajor> Ainv(n, n);
    Ainv.SetIdentity();
    
    SolveLU(mat_lu, Ainv);
    
    Matrix<T, General, ColMajor> Adense(n, n), mat_id(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	Adense(i, j) = A(i, j);
    
    Mlt(Adense, Ainv, mat_id);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  T val_ex = 0.0;
	  if (i == j)
	    val_ex = 1.0;
	  
	  if (abs(mat_id(i, j) - val_ex) > threshold)
	    {
	      cout << "Solve of Wsmp by a matrix is incorrect" << endl;
	      abort();
	    }
	}
    
    // testing resolution with a complex vector
    CheckSolveComplexVector(A, mat_lu, x);
    
    mat_lu.Clear();
#else
    cout << "Not testing Wsmp" << endl;
#endif
  }

#ifdef SELDON_WITH_SUPERLU

  bool test_superlu = true;
#ifdef SELDON_WITH_SUPERLU_DIST
#ifdef SELDON_WITH_SUPERLU_DOUBLE
  test_superlu = !IsComplexMatrix(A);
#else
  test_superlu = IsComplexMatrix(A);
#endif
#endif

  if (test_superlu)
  {
    DISP(sizeof(T));
    MatrixSuperLU<T> mat_lu;
    mat_lu.HideMessages();
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "FactorizeMatrix of SuperLU incorrect" << endl;
	abort();
      }
    
    mat_lu.HideMessages();
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of SuperLU incorrect " << endl;
	abort();
      }    
    
    x = b;
    mat_lu.Solve(SeldonNoTrans, x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of SuperLU incorrect " << endl;
	abort();
      }

    x = bt;
    mat_lu.Solve(SeldonTrans, x);

    if (!EqualVector(x, y))
      {
	cout << "Solve of SuperLU incorrect " << endl;
	abort();
      }

    // testing resolution of multiple right hand-sides
    Matrix<T, General, ColMajor> Ainv(n, n);
    Ainv.SetIdentity();
    
    SolveLU(mat_lu, Ainv);
    
    Matrix<T, General, ColMajor> Adense(n, n), mat_id(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	Adense(i, j) = A(i, j);
    
    Mlt(Adense, Ainv, mat_id);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  T val_ex = 0.0;
	  if (i == j)
	    val_ex = 1.0;
	  
	  if (abs(mat_id(i, j) - val_ex) > threshold)
	    {
	      cout << "Solve of SuperLU by a matrix is incorrect" << endl;
	      abort();
	    }
	}
    
    // testing resolution with a complex vector
    CheckSolveComplexVector(A, mat_lu, x);
    
    mat_lu.Clear();
    
    // testing with another ordering
    mat_lu.SelectOrdering(superlu::MMD_AT_PLUS_A);
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "SelectOrdering of SuperLU incorrect" << endl;
	abort();
      }
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of SuperLU incorrect " << endl;
	abort();
      }    
    
    mat_lu.SetPermutation(num);
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "SetPermutation of SuperLU incorrect" << endl;
	abort();
      }
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of SuperLU incorrect " << endl;
	abort();
      }

#ifndef SELDON_WITH_SUPERLU_MT
#ifndef SELDON_WITH_SUPERLU_DIST 
    {
      Matrix<T, General, RowSparse> L, U;
      mat_lu.GetLU(L, U);
      IVect row_perm = mat_lu.GetRowPermutation();
      IVect col_perm = mat_lu.GetColPermutation();
      
      for (int i = 0; i < n; i++)
	{
	  int j = L.GetPtr()[i+1];
	  if ((L.GetInd()[j-1] != i) || (L.GetData()[j-1] != one))
	    {
	      cout << "L not triangular" << endl;
	      abort();
	    }
	  
	  j = U.GetPtr()[i];
	  if (U.GetInd()[j] != i)
	    {
	      cout << "U not triangular" << endl;
	      abort();
	    }
	}
      
      for (int i = 0; i < n; i++)
	x(col_perm(i)) = y(i);
      
      Mlt(U, x, z);
      Mlt(L, z, x);
      
      for (int i = 0; i < n; i++)
	z(i) = x(row_perm(i));
      
      if (!EqualVector(z, b))
	{
	  cout << "GetLU incorrect" << endl;
	  abort();
	}
    }
    
    {
      Matrix<T, General, ColSparse> L, U;
      mat_lu.GetLU(L, U);
      IVect row_perm = mat_lu.GetRowPermutation();
      IVect col_perm = mat_lu.GetColPermutation();
      //A.WriteText("A.dat");
      //L.WriteText("L.dat");
      //U.WriteText("U.dat");
      //mat_lu.GetRowPermutation().WriteText("row_perm.dat");
      //mat_lu.GetColPermutation().WriteText("col_perm.dat");
      
      for (int i = 0; i < n; i++)
	{
	  int j = L.GetPtr()[i];
	  if ((L.GetInd()[j] != i) || (L.GetData()[j] != one))
	    {
	      cout << "L not triangular" << endl;
	      abort();
	    }
	  
	  j = U.GetPtr()[i+1];
	  if (U.GetInd()[j-1] != i)
	    {
	      cout << "U not triangular" << endl;
	      abort();
	    }
	}
      
      for (int i = 0; i < n; i++)
	x(col_perm(i)) = y(i);
      
      Mlt(U, x, z);
      Mlt(L, z, x);
      
      for (int i = 0; i < n; i++)
	z(i) = x(row_perm(i));
      
      if (!EqualVector(z, b))
	{
	  cout << "GetLU incorrect" << endl;
	  abort();
	}
    }
#endif
#endif
  }
#else
    cout << "Not testing SuperLU" << endl;
#endif
  
  {
#ifdef SELDON_WITH_UMFPACK
    MatrixUmfPack<T> mat_lu;
    mat_lu.HideMessages();
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "FactorizeMatrix of UmfPack incorrect" << endl;
	abort();
      }
    
    mat_lu.HideMessages();
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of UmfPack incorrect " << endl;
	abort();
      }    
    
    x = b;
    mat_lu.Solve(SeldonNoTrans, x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of UmfPack incorrect " << endl;
	abort();
      }

    x = bt;
    mat_lu.Solve(SeldonTrans, x);

    if (!EqualVector(x, y))
      {
	cout << "Solve of UmfPack incorrect " << endl;
	abort();
      }

    // testing resolution of multiple right hand-sides
    Matrix<T, General, ColMajor> Ainv(n, n);
    Ainv.SetIdentity();
    
    SolveLU(mat_lu, Ainv);
    
    Matrix<T, General, ColMajor> Adense(n, n), mat_id(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	Adense(i, j) = A(i, j);
    
    Mlt(Adense, Ainv, mat_id);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  T val_ex = 0.0;
	  if (i == j)
	    val_ex = 1.0;
	  
	  if (abs(mat_id(i, j) - val_ex) > threshold)
	    {
	      cout << "Solve of UmfPack by a matrix is incorrect" << endl;
	      abort();
	    }
	}
    
    // testing resolution with a complex vector
    CheckSolveComplexVector(A, mat_lu, x);
    
    mat_lu.Clear();

    mat_lu.SelectOrdering(UMFPACK_ORDERING_METIS);
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "SelectOrdering of UmfPack incorrect" << endl;
	abort();
      }
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of UmfPack incorrect " << endl;
	abort();
      }    
    
    /*mat_lu.SetPermutation(num);
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "SetPermutation of UmfPack incorrect" << endl;
	abort();
      }
    */
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of UmfPack incorrect " << endl;
	abort();
      }
#else
    cout << "Not testing UmfPack" << endl;
#endif
  }

  {
#ifdef SELDON_WITH_MUMPS
    MatrixMumps<T> mat_lu;
    mat_lu.HideMessages();
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "FactorizeMatrix of Mumps incorrect" << endl;
	abort();
      }
    
    mat_lu.HideMessages();
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of Mumps incorrect " << endl;
	abort();
      }    
    
    x = b;
    mat_lu.Solve(SeldonNoTrans, x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of Mumps incorrect " << endl;
	abort();
      }

    x = bt;
    mat_lu.Solve(SeldonTrans, x);

    if (!EqualVector(x, y))
      {
	cout << "Solve of Mumps incorrect " << endl;
	abort();
      }
    
    // testing resolution of multiple right hand-sides
    Matrix<T, General, ColMajor> Ainv(n, n);
    Ainv.SetIdentity();
    
    SolveLU(mat_lu, Ainv);
    
    Matrix<T, General, ColMajor> Adense(n, n), mat_id(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	Adense(i, j) = A(i, j);
    
    Mlt(Adense, Ainv, mat_id);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  T val_ex = 0.0;
	  if (i == j)
	    val_ex = 1.0;
	  
	  if (abs(mat_id(i, j) - val_ex) > threshold)
	    {
	      cout << "Solve of Mumps by a matrix is incorrect" << endl;
	      abort();
	    }
	}
    
    // testing resolution with a complex vector
    CheckSolveComplexVector(A, mat_lu, x);
    
    mat_lu.Clear();
    
    // testing resolution with a user-defined ordering
    mat_lu.SelectOrdering(6);
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "SelectOrdering of Mumps incorrect" << endl;
	abort();
      }
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of Mumps incorrect " << endl;
	abort();
      }    
    
    mat_lu.SetPermutation(num);
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "SetPermutation of Mumps incorrect" << endl;
	abort();
      }

    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y, 100.0*threshold))
      {
	cout << "Solve of Mumps incorrect " << endl;
	abort();
      }
    
    mat_lu.SelectOrdering(5);    
    IVect perm(n);
    perm.Fill(0);
    mat_lu.FindOrdering(A,  perm, true);
    CheckPermutation(perm);

    if (n > 100)
      {
	int p = 80;
	int p2 = 0;
	IVect Index(n);
	Index.Fill(-1);
	IVect nums(p);
	nums.Fill(-1);
	while (p2 < p)
	  {
	    int i = rand()%n;
	    if (Index(i) == -1)
	      {
		Index(i) = p2;
		p2++;
	      }
	  }
	
	p2 = 0;
	for (int i = 0; i < n; i++)
	  if (Index(i) >= 0)
	    {
	      nums(p2) = i;
	      Index(i) = p2++;
	    }
	
	IVect Index2(n);
	Index2.Fill(-1);
	p2 = 0;
	for (int i = 0; i < n; i++)
	  if (Index(i) == -1)
	    Index2(i) = p2++;
	
	Matrix<T, General, ColMajor> mat_schur, A11(p, p), A21(n-p, p),
	  A12(p, n-p), A22(n-p, n-p), C(n-p, p);
	
	C.Fill(zero);
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < n; j++)
	    {
	      if (Index(i) == -1)
		{
		  if (Index(j) == -1)
		    A22(Index2(i), Index2(j)) = A(i, j);
		  else
		    A21(Index2(i), Index(j)) = A(i, j);
		}
	      else
		{
		  if (Index(j) == -1)
		    A12(Index(i), Index2(j)) = A(i, j);
		  else
		    A11(Index(i), Index(j)) = A(i, j);
		}
	    }
	
	GetInverse(A22);
	Mlt(A22, A21, C);
	MltAdd(-one, A12, C, one, A11);    
	
	mat_lu.GetSchurMatrix(A, nums, mat_schur, true);
	for (int i = 0; i < p; i++)
	  for (int j = 0; j < p; j++)
	    if (abs(mat_schur(i, j) - A11(i, j)) > 20.0*threshold)
	      {
		cout << "GetSchurMatrix incorrect" << endl;
		abort();
	      }
      }
    
#else
    cout << "Not testing Mumps" << endl;
#endif
  }


  {
#ifdef SELDON_WITH_PARDISO
    MatrixPardiso<T> mat_lu;
    mat_lu.HideMessages();
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "FactorizeMatrix of Pardiso incorrect" << endl;
	abort();
      }
    
    mat_lu.HideMessages();
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of Pardiso incorrect " << endl;
	abort();
      }    
    
    x = b;
    mat_lu.Solve(SeldonNoTrans, x);
    
    if (!EqualVector(x, y))
      {
	cout << "Solve of Pardiso incorrect " << endl;
	abort();
      }

    x = bt;
    mat_lu.Solve(SeldonTrans, x);

    if (!EqualVector(x, y))
      {
        //DISP(x); DISP(y);
	//cout << "Solve of Pardiso incorrect " << endl;
	//abort();
      }
    
    // testing resolution of multiple right hand-sides
    Matrix<T, General, ColMajor> Ainv(n, n);
    Ainv.SetIdentity();
    
    SolveLU(mat_lu, Ainv);
    
    Matrix<T, General, ColMajor> Adense(n, n), mat_id(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	Adense(i, j) = A(i, j);
    
    Mlt(Adense, Ainv, mat_id);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  T val_ex = 0.0;
	  if (i == j)
	    val_ex = 1.0;
	  
	  if (abs(mat_id(i, j) - val_ex) > threshold)
	    {
	      cout << "Solve of Pardiso by a matrix is incorrect" << endl;
	      abort();
	    }
	}
    
    // testing resolution with a complex vector
    CheckSolveComplexVector(A, mat_lu, x);
    
    mat_lu.Clear();
    
    // testing resolution with a user-defined ordering
    mat_lu.SetPermutation(num);
    mat_lu.FactorizeMatrix(A, true);
    
    if (mat_lu.GetInfoFactorization() != 0)
      {
	cout << "SetPermutation of Pardiso incorrect" << endl;
	abort();
      }

    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y, 100.0*threshold))
      {
	cout << "Solve of Pardiso incorrect " << endl;
	abort();
      }
#else
    cout << "Not testing Pardiso" << endl;
#endif
  }
  
  //WriteMatrixMarket(A, "mat_unsym.mtx");
  //y.WriteText("sol.dat");
  //b.WriteText("rhs.dat");
  
  {
#ifdef SELDON_WITH_PASTIX
    MatrixPastix<T> mat_lu;
    //if (IsSymmetricMatrix(A))
    mat_lu.DoNotRefineSolution();
    
    mat_lu.HideMessages();
    mat_lu.FactorizeMatrix(A, true);
    
    //mat_lu.HideMessages();
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y, 10.0*threshold))
      {
	cout << "Solve of Pastix incorrect " << endl;
	abort();
      }    
    
    x = b;
    mat_lu.Solve(SeldonNoTrans, x);
    
    if (!EqualVector(x, y, 10.0*threshold))
      {
	cout << "Solve of Pastix incorrect " << endl;
	abort();
      }
    
    x = bt;
    mat_lu.Solve(SeldonTrans, x);

    if (!EqualVector(x, y))
      {
	cout << "Solve of Pastix incorrect " << endl;
	abort();
      }    

    // testing resolution of multiple right hand-sides
    Matrix<T, General, ColMajor> Ainv(n, n);
    Ainv.SetIdentity();
    
    SolveLU(mat_lu, Ainv);
    
    Matrix<T, General, ColMajor> Adense(n, n), mat_id(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	Adense(i, j) = A(i, j);
    
    Mlt(Adense, Ainv, mat_id);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  T val_ex = 0.0;
	  if (i == j)
	    val_ex = 1.0;
	  
	  if (abs(mat_id(i, j) - val_ex) > threshold)
	    {
	      cout << "Solve of Pastix by a matrix is incorrect" << endl;
	      abort();
	    }
	}
    
    // testing resolution with a complex vector
    CheckSolveComplexVector(A, mat_lu, x);
    
    mat_lu.Clear();
    
    // testing with another ordering
    mat_lu.SelectOrdering(API_ORDER_SCOTCH);
    mat_lu.FactorizeMatrix(A, true);
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y, 10.0*threshold))
      {
	cout << "Solve of Pastix incorrect " << endl;
	abort();
      }    
    
    /*mat_lu.SetPermutation(num);
    mat_lu.FactorizeMatrix(A, true);
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y, 10.0*threshold))
      {
	cout << "Solve of Pastix incorrect " << endl;
	abort();
      }
    */
    mat_lu.SelectOrdering(API_ORDER_SCOTCH);    
    IVect perm(n);
    perm.Fill(0);
    mat_lu.FindOrdering(A,  perm, true);
    CheckPermutation(perm);
#else
    cout << "Not testing Pastix" << endl;
#endif
  }

  if (false)
  {
    SparseSeldonSolver<T> mat_lu;
        
    mat_lu.SetPivotThreshold(0.1);
    if (mat_lu.GetPivotThreshold() != 0.1)
      {
	cout << "GetPivotThreshold incorrect" << endl;
	abort();
      }
    
    IVect perm(n);
    perm.Fill(-1);
    FindSparseOrdering(A, perm, SparseMatrixOrdering::AMD);
    CheckPermutation(perm);

    perm.Fill(-1);
    FindSparseOrdering(A, perm, SparseMatrixOrdering::PORD);
    CheckPermutation(perm);

    perm.Fill(-1);
    FindSparseOrdering(A, perm, SparseMatrixOrdering::SCOTCH);
    CheckPermutation(perm);

    perm.Fill(-1);
    FindSparseOrdering(A, perm, SparseMatrixOrdering::METIS);
    CheckPermutation(perm);

    perm.Fill(-1);
    FindSparseOrdering(A, perm, SparseMatrixOrdering::COLAMD);
    CheckPermutation(perm);

    perm.Fill(-1);
    FindSparseOrdering(A, perm, SparseMatrixOrdering::QAMD);
    CheckPermutation(perm);

    perm.Fill(-1);
    FindSparseOrdering(A, perm, SparseMatrixOrdering::AUTO);
    CheckPermutation(perm);

    
    mat_lu.ShowMessages();
    mat_lu.FactorizeMatrix(perm, A, true);
    
    mat_lu.HideMessages();
    
    x = b;
    mat_lu.Solve(x);
    
    if (!EqualVector(x, y, threshold))
      {
	cout << "Solve of Seldon incorrect " << endl;
	abort();
      }    
    
    x = b;
    mat_lu.Solve(SeldonNoTrans, x);
    
    if (!EqualVector(x, y, threshold))
      {
	cout << "Solve of Seldon incorrect " << endl;
	abort();
      }
    
    x = bt;
    mat_lu.Solve(SeldonTrans, x);

    if (!EqualVector(x, y))
      {
	cout << "Solve of Seldon incorrect " << endl;
	abort();
      }
    
    mat_lu.Clear();
  }

  {
    SparseDirectSolver<T> mat_lu;
    
    for (int t = 0; t < 7; t++)
      if (SparseDirectSolver<T>::IsAvailableSolver(t))
	{
	  cout << "Testing Solver " << t << endl;
	  mat_lu.SelectDirectSolver(t);
	  
	  bool test_superlu = true;
#ifdef SELDON_WITH_SUPERLU_DIST
#ifdef SELDON_WITH_SUPERLU_DOUBLE
	  test_superlu = !IsComplexMatrix(A);
#else
	  test_superlu = IsComplexMatrix(A);
#endif
#endif

	  bool test_solver = true;
	  if (t == SparseDirectSolver<T>::SUPERLU)
	    test_solver = test_superlu;

          DISP(test_solver);
	  if (test_solver)
	    {
	      mat_lu.ShowMessages();
	      mat_lu.Factorize(A, true);
	      
	      int ierr = 0;
	      int info = mat_lu.GetInfoFactorization(ierr);
	      if (info != mat_lu.FACTO_OK)
		{
		  cout << "Factorize incorrect" << endl;
		  abort();
		}
	      
	      mat_lu.HideMessages();
	      
	      x = b;
	      mat_lu.Solve(x);
	      
	      if (!EqualVector(x, y, 10.0*threshold))
		{
		  cout << "Solve of SparseDirectSolver incorrect " << endl;
		  abort();
		}    
	      
	      x = b;
	      mat_lu.Solve(SeldonNoTrans, x);
	      
	      if (!EqualVector(x, y, 10.0*threshold))
		{
		  cout << "Solve of SparseDirectSolver incorrect " << endl;
		  abort();
		}
	      
	      x = bt;
	      mat_lu.Solve(SeldonTrans, x);
	      
	      if (!EqualVector(x, y))
		{
		  cout << "Solve of SparseDirectSolver incorrect " << endl;
		  abort();
		}
	      
	  // testing resolution of multiple right hand-sides
	      x_mat = b_mat;
	      mat_lu.Solve(x_mat);
	      
	      if (!EqualMatrix(x_mat, y_mat))
		{
		  cout << "Solve of SparseDirectSolver incorrect" << endl;
		  abort();
		}
	      
	      x_mat = b_mat;
	      mat_lu.Solve(SeldonNoTrans, x_mat);
	      
	      if (!EqualMatrix(x_mat, y_mat))
		{
		  cout << "Solve of SparseDirectSolver incorrect" << endl;
		  abort();
		}
	      
	      x_mat = bt_mat;
	      mat_lu.Solve(SeldonTrans, x_mat);
	      
	      if (!EqualMatrix(x_mat, y_mat))
		{
		  cout << "Solve of SparseDirectSolver incorrect" << endl;
		  abort();
		}
	      
	      mat_lu.Clear();
	    }
	}
  }

}

int main(int argc, char** argv)
{
  InitSeldon(argc, argv);

  threshold = 1e-10;
  
  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckDirectSolver(A);
  }
  
  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    CheckDirectSolver(A);
  }

  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckDirectSolver(A);
  }
  
  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    CheckDirectSolver(A);
  }
  
  cout << "All tests passed successfully" << endl;
  return FinalizeSeldon();
}

