#include "SeldonLib.hxx"

using namespace Seldon;

double threshold = 1e-10;

/************************
 * Checking eigenvalues *
 ************************/


// checking symmetric standard eigenproblems
template<class MatrixS, class Vector1, class Matrix1>
void CheckEigenvalues(const MatrixS& mat_stiff,
                      const Vector1& lambda, const Matrix1& eigen_vec)
{
  int N = mat_stiff.GetM();
  typedef typename Matrix1::entry_type T1;
  Vector<T1> X(N), Y(N);
  X.Fill(0);
  Y.Fill(0);
  for (int i = 0; i < lambda.GetM(); i++)
    {
      for (int j = 0; j < N; j++)
        X(j) = eigen_vec(j, i);
      
      Mlt(mat_stiff, X, Y);
      double err = 0;
      double normeX = sqrt(abs(DotProdConj(X, X)));
      for (int j = 0; j < N; j++)
        err += pow(abs(Y(j) - lambda(i)*X(j)), 2.0);
      
      err = sqrt(err);
      if (err > threshold*normeX)
	{
	  cout << "Error on eigenvalue " << lambda(i) << endl;
	  cout << "Error = " << err/normeX << endl;
	  abort();
	}
    }
  
  
}


// checking unsymmetric standard eigenproblems
template<class Prop, class Storage, class Matrix1>
void CheckEigenvalues(const Matrix<double, Prop, Storage>& mat_stiff,
                      const Vector<double>& lambda,
                      const Vector<double>& lambda_imag, const Matrix1& eigen_vec)
{
  int N = mat_stiff.GetM();
  Vector<double> X(N), Xi(N), Y(N), Yi(N);
  X.Fill(0); Xi.Fill(0);
  Y.Fill(0); Yi.Fill(0);
  
  int i = 0;
  while (i < lambda.GetM())
    {
      bool eigen_pair = false;
      if (i < lambda.GetM()-1)
        {
          if ( (lambda(i) == lambda(i+1)) && (lambda_imag(i) == -lambda_imag(i+1)))
            eigen_pair = true;
        }

      if ((lambda_imag(i) != 0) && (!eigen_pair))
	{
	  DISP(i); DISP(lambda.GetM());
	  cout << "Eigenpair at the end of the list" << endl;
	  break;
	}
      
      double err = 0, normeX(1);
      if (eigen_pair)
        {
          for (int j = 0; j < N; j++)
            {
              X(j) = eigen_vec(j, i);
              Xi(j) = eigen_vec(j, i+1);
            }
          
          Mlt(mat_stiff, X, Y);
          Mlt(mat_stiff, Xi, Yi);
	  normeX = sqrt(DotProd(Xi, Xi) + DotProd(X, X));
          for (int j = 0; j < N; j++)
            err += pow(abs(complex<double>(Y(j), Yi(j))
                           - complex<double>(lambda(i), lambda_imag(i))
                           *complex<double>(X(j), Xi(j)) ), 2.0);
          
          err = sqrt(err);
        }
      else
        {
          for (int j = 0; j < N; j++)
            X(j) = eigen_vec(j, i);
          
          Mlt(mat_stiff, X, Y);
	  normeX = sqrt(DotProd(X, X));
          for (int j = 0; j < N; j++)
            err += pow(Y(j) - lambda(i)*X(j), 2.0);
          
          err = sqrt(err);
        }
      
      if (err > threshold*normeX)
	{
	  cout << "Error on eigenvalue " << lambda(i) << endl;
	  cout << "Error = " << err/normeX << endl;
	  abort();
	}
      
      if (eigen_pair)
        i += 2;
      else
        i++;
    }
  
  
}

template<class MatrixS, class Matrix1>
void CheckEigenvalues(const MatrixS& mat_stiff,
                      const Vector<complex<double> >& lambda,
                      const Vector<complex<double> >& lambda_imag, const Matrix1& eigen_vec)
{
}


// checking unsymmetric generalized eigenproblems
template<class Prop0, class Storage0, class Prop1, class Storage1, class Matrix1>
void CheckEigenvalues(const Matrix<double, Prop0, Storage0>& mat_stiff,
                      const Matrix<double, Prop1, Storage1>& mat_mass,
                      const Vector<double>& lambda,
                      const Vector<double>& lambda_imag, const Matrix1& eigen_vec)
{
  int N = mat_stiff.GetM();
  Vector<double> X(N), Xi(N), Y(N), Yi(N), Mx(N), Mxi(N);
  X.Fill(0); Xi.Fill(0);
  Mx.Fill(0); Mxi.Fill(0);
  Y.Fill(0); Yi.Fill(0);
  
  int i = 0;
  while (i < lambda.GetM())
    {
      bool eigen_pair = false;
      if (i < lambda.GetM()-1)
        {
          if ( (lambda(i) == lambda(i+1)) && (lambda_imag(i) == -lambda_imag(i+1)))
            eigen_pair = true;
        }

      if ((lambda_imag(i) != 0) && (!eigen_pair))
	{
	  DISP(i); DISP(lambda.GetM());
	  cout << "Eigenpair at the end of the list" << endl;
	  break;
	}
      
      double err = 0, normeX(1);
      if (eigen_pair)
        {
          for (int j = 0; j < N; j++)
            {
              X(j) = eigen_vec(j, i);
              Xi(j) = eigen_vec(j, i+1);
            }
          
          Mlt(mat_stiff, X, Y);
          Mlt(mat_stiff, Xi, Yi);
          Mlt(mat_mass, X, Mx);
          Mlt(mat_mass, Xi, Mxi);
	  normeX = sqrt(DotProd(Xi, Xi) + DotProd(X, X));
          for (int j = 0; j < N; j++)
            err += pow(abs(complex<double>(Y(j), Yi(j))
                           - complex<double>(lambda(i), lambda_imag(i))
                           *complex<double>(Mx(j), Mxi(j)) ), 2.0);
          
          err = sqrt(err);
        }
      else
        {
          for (int j = 0; j < N; j++)
            X(j) = eigen_vec(j, i);
          
          Mlt(mat_stiff, X, Y);
          Mlt(mat_mass, X, Mx);
	  normeX = sqrt(DotProd(X, X));
          for (int j = 0; j < N; j++)
            err += pow(Y(j) - lambda(i)*Mx(j), 2.0);
          
          err = sqrt(err);
        }
      
      if (err > threshold*normeX)
        {
          cout << "Error on eigenvalue " << lambda(i) << endl;
          cout << "Error = " << err/normeX << endl;
          abort();
        }
      
      if (eigen_pair)
        i += 2;
      else
        i++;
    }
  
  
}

template<class MatrixS, class MatrixM, class Matrix1>
void CheckEigenvalues(const MatrixS& mat_stiff,
                      const MatrixM& mat_mass,
                      const Vector<complex<double> >& lambda,
                      const Vector<complex<double> >& lambda_imag, const Matrix1& eigen_vec)
{}


// checking complex generalized eigenproblems
template<class MatrixS, class MatrixM, class Vector1, class Matrix1>
void CheckEigenvalues(const MatrixS& mat_stiff, const MatrixM& mat_mass,
                      const Vector1& lambda, const Matrix1& eigen_vec)
{
  int N = mat_stiff.GetM();
  typedef typename Matrix1::entry_type T1;
  Vector<T1> X(N), Y(N), Mx(N);
  X.Fill(0);
  Y.Fill(0);
  Mx.Fill(0);
  for (int i = 0; i < lambda.GetM(); i++)
    {
      for (int j = 0; j < N; j++)
        X(j) = eigen_vec(j, i);
      
      Mlt(mat_stiff, X, Y);
      Mlt(mat_mass, X, Mx);
      double normeX = sqrt(abs(DotProdConj(X, X)));
      double err = 0;
      for (int j = 0; j < N; j++)
        err += pow(abs(Y(j) - lambda(i)*Mx(j)), 2.0);
      
      err = sqrt(err);
      if (err > threshold*normeX)
        {
          cout << "Error on eigenvalue " << lambda(i) << endl;
          cout << "Error = " << err/normeX << endl;
          abort();
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

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz)
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


// copy from a sparse matrix to a complexe dense matrix
template<class T0, class T1>
void Copy(const Matrix<T0, Symmetric, ArrayRowSymSparse>& Asp, Matrix<T1>& A)
{
  int N = Asp.GetM();
  A.Reallocate(N, N);
  A.Fill(0);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < Asp.GetRowSize(i); j++)
      {
	A(i, Asp.Index(i, j)) = Asp.Value(i, j);
	A(Asp.Index(i, j), i) = Asp.Value(i, j);
      }
}

// copy from a sparse matrix to a complexe dense matrix
template<class T0, class T1>
void Copy(const Matrix<T0, General, ArrayRowSparse>& Asp, Matrix<T1>& A)
{
  int N = Asp.GetM();
  A.Reallocate(N, N);
  A.Fill(0);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < Asp.GetRowSize(i); j++)
      A(i, Asp.Index(i, j)) = Asp.Value(i, j);
}

template<class MatrixK, class MatrixM>
void FindReferenceEigenvalues(MatrixK& K, MatrixM& M, Vector<complex<double> >& L)
{
  Matrix<complex<double> > invM, Kd, A;
  
  Copy(M, invM);
  GetInverse(invM);
  
  Copy(K, Kd);
  
  A.Reallocate(K.GetM(), K.GetM());
  Mlt(invM, Kd, A);
  
  Kd.Clear(); invM.Clear();
  GetEigenvalues(A, L);
}


template<class MatrixK, class MatrixM, class T>
void TestGeneralProblem(MatrixK& K, MatrixM& M, T&,
			bool add_pos_diag, int nb_eigenval, bool test_cplx_shift = true,
			bool regular_mode = false)
{
  Vector<T>  lambda, lambda_imag;
  Matrix<T, General, ColMajor> eigen_vec;
  
  int m = 40, n = m, nnz = 4*m;
  GenerateRandomMatrix(K, m, n, nnz);
  GenerateRandomMatrix(M, m, n, nnz);
  
  typename MatrixM::entry_type coef_diag;  
  for (int i = 0; i < M.GetM(); i++)
    {
      GetRandNumber(coef_diag);
      if (add_pos_diag)
	coef_diag += 5.0;
      
      M.AddInteraction(i, i, coef_diag);
    }
    
  K.WriteText("K.dat");
  M.WriteText("M.dat");

  SparseEigenProblem<T, MatrixK, MatrixM> var_eig;

  var_eig.SetStoppingCriterion(1e-14);
  var_eig.SetNbAskedEigenvalues(nb_eigenval);
  if (regular_mode)
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
  else
    var_eig.SetComputationalMode(var_eig.INVERT_MODE);
  
  var_eig.InitMatrix(K, M);
  
  // Large eigenvalues
  var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
  DISP(lambda); DISP(lambda_imag);
  
  if (sizeof(T) == 16)
    CheckEigenvalues(K, M, lambda, eigen_vec);
  else
    CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);

  // eigenvalues close to a real shift
  if (regular_mode)
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);

  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 0.02, var_eig.SORTED_MODULUS);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
  DISP(lambda); DISP(lambda_imag);
  
  if (sizeof(T) == 16)
    CheckEigenvalues(K, M, lambda, eigen_vec);
  else
    CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);

  // eigenvalues close to a complex shift
  if (test_cplx_shift)
    {
      complex<double> shift_cplx;
      GetRandNumber(shift_cplx); DISP(shift_cplx);
      var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES,
			      shift_cplx, var_eig.SORTED_MODULUS);
      
      GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
      DISP(lambda); DISP(lambda_imag);
      
      if (sizeof(T) == 16)
	CheckEigenvalues(K, M, lambda, eigen_vec);  
      else
	CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);  
      
      cout << endl << endl;
    }
}


template<class MatrixK, class MatrixM, class T>
void TestSymProblem(MatrixK& K, MatrixM& M, T&,
		    bool add_pos_diag, int nb_eigenval, bool add_pos_diag_stiff)
{
  Vector<T>  lambda, lambda_imag;
  Matrix<T, General, ColMajor> eigen_vec;
  
  int m = 40, n = m, nnz = 4*m;
  GenerateRandomMatrix(K, m, n, nnz);
  GenerateRandomMatrix(M, m, n, nnz);
  
  typename MatrixM::entry_type coef_diag;  
  for (int i = 0; i < M.GetM(); i++)
    {
      GetRandNumber(coef_diag);
      if (add_pos_diag)
	coef_diag += 5.0;
      
      M.AddInteraction(i, i, coef_diag);

      if (add_pos_diag_stiff)
	{
	  GetRandNumber(coef_diag);
	  coef_diag += 5.0;
	  
	  K.AddInteraction(i, i, coef_diag);
	}
    }
    
  K.WriteText("K.dat");
  M.WriteText("M.dat");

  SparseEigenProblem<T, MatrixK, MatrixM> var_eig;

  var_eig.SetStoppingCriterion(1e-14);
  var_eig.SetNbAskedEigenvalues(nb_eigenval);
  var_eig.InitMatrix(K, M);
  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 0.02, var_eig.SORTED_MODULUS);  
  
  if (add_pos_diag_stiff)
    {
      // using Buckling mode
      var_eig.SetComputationalMode(var_eig.BUCKLING_MODE);
      
      GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
      DISP(lambda); DISP(lambda_imag);
      
      CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);
    }
  
  // using Cayley mode
  var_eig.SetComputationalMode(var_eig.CAYLEY_MODE);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
  DISP(lambda); DISP(lambda_imag);
  
  CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);  
  
  cout << endl << endl;
}


template<class MatrixK, class MatrixM, class T>
void TestStandardProblem(MatrixK& K, MatrixM& M, T&,
			 bool add_pos_diag, int nb_eigenval, bool test_chol)
{
  Vector<T>  lambda, lambda_imag;
  Matrix<T, General, ColMajor> eigen_vec;
  
  int m = 40, n = m, nnz = 4*m;
  GenerateRandomMatrix(K, m, n, nnz);
  GenerateRandomMatrix(M, m, n, nnz);
  
  typename MatrixM::entry_type coef_diag;  
  for (int i = 0; i < M.GetM(); i++)
    {
      GetRandNumber(coef_diag);
      if (add_pos_diag)
	coef_diag += 5.0;
      
      M.AddInteraction(i, i, coef_diag);
    }
  
  K.WriteText("K.dat");
  M.WriteText("M.dat");
  
  SparseEigenProblem<T, MatrixK, MatrixM> var_eig;
  
  // first, we test the case without the mass matrix M
  
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(nb_eigenval);
  var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
  
  var_eig.InitMatrix(K);
  var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
  DISP(lambda); DISP(lambda_imag);
  
  if (sizeof(T) == 16)
    CheckEigenvalues(K, lambda, eigen_vec);
  else
    CheckEigenvalues(K, lambda, lambda_imag, eigen_vec);
  
  // eigenvalues close to a real shift
  var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 0.02, var_eig.SORTED_MODULUS);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
  DISP(lambda); DISP(lambda_imag);
  
  if (sizeof(T) == 16)
    CheckEigenvalues(K, lambda, eigen_vec);
  else
    CheckEigenvalues(K, lambda, lambda_imag, eigen_vec);
  
  // case with a diagonal matrix M    
  MatrixM Mdiag(n, n);
  Vector<typename MatrixM::entry_type> D;
  GenerateRandomVector(D, n);
  
  D.Write("D.dat");
  
  for (int i = 0; i < n; i++)
    Mdiag.AddInteraction(i, i, D(i));
  
  var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
  
  var_eig.InitMatrix(K, Mdiag);
  var_eig.SetDiagonalMass();
  var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
  DISP(lambda); DISP(lambda_imag);
  
  if (sizeof(T) == 16)
    CheckEigenvalues(K, Mdiag, lambda, eigen_vec);
  else
    CheckEigenvalues(K, Mdiag, lambda, lambda_imag, eigen_vec);

  // eigenvalues close to a real shift
  var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 0.02, var_eig.SORTED_MODULUS);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
  DISP(lambda); DISP(lambda_imag);
  
  if (sizeof(T) == 16)
    CheckEigenvalues(K, Mdiag, lambda, eigen_vec);
  else
    CheckEigenvalues(K, Mdiag, lambda, lambda_imag, eigen_vec);
  
  
  // case with a symmetric positive definite mass matrix
  if (test_chol)
    {
      var_eig.InitMatrix(K, M);
      var_eig.SetCholeskyFactoForMass();
      
      var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
      var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);
      
      GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
      DISP(lambda); DISP(lambda_imag);
      
      if (sizeof(T) == 16)
	CheckEigenvalues(K, M, lambda, eigen_vec);
      else
	CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);
      
      // eigenvalues close to a real shift
      var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
      var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 0.02, var_eig.SORTED_MODULUS);
      
      GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
      DISP(lambda); DISP(lambda_imag);
      
      if (sizeof(T) == 16)
	CheckEigenvalues(K, M, lambda, eigen_vec);
      else
	CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);
    }
}

int main(int argc, char** argv)
{
  InitSeldon(argc, argv);
  
  bool all_test = true;
  int nb_eigenval = 6;
  
  // testing regular mode (generalized problem)
  {
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, true, nb_eigenval, false, true);
  }

  {
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    double x;
    TestGeneralProblem(K, M, x, true, nb_eigenval, true, true);
  }

  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, true, nb_eigenval, true, true);
  }

  {
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    double x;
    TestGeneralProblem(K, M, x, true, nb_eigenval, true, true);
  }
  
  // testing buckling/cayley mode (real symmetric generalized problem)
  {
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    double x;
    TestSymProblem(K, M, x, true, nb_eigenval, true);
    TestSymProblem(K, M, x, true, nb_eigenval, false);
  }
  
  // testing standard mode
  {
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestStandardProblem(K, M, x, true, nb_eigenval, true);
  }

  {
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestStandardProblem(K, M, x, true, nb_eigenval, false);
  }

  {
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    double x;
    TestStandardProblem(K, M, x, true, nb_eigenval, true);
  }

  {
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestStandardProblem(K, M, x, true, nb_eigenval, false);
  }

  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestStandardProblem(K, M, x, true, nb_eigenval, true);
  }

  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestStandardProblem(K, M, x, true, nb_eigenval, false);
  }

  {
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    double x;
    TestStandardProblem(K, M, x, true, nb_eigenval, true);
  }

  {
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestStandardProblem(K, M, x, true, nb_eigenval, false);
  }
  
  // testing invert mode
  {
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval);
  }
  
  {
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval);
  }

  {
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<double, General, ArrayRowSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval);
  }

  {
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<complex<double>, General, ArrayRowSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval);
  }
  
  {
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    double x;
    TestGeneralProblem(K, M, x, true, nb_eigenval, false);
  }

  {
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, true, nb_eigenval);
  }

  {
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<double, General, ArrayRowSparse> M;
    double x;
    TestGeneralProblem(K, M, x, true, nb_eigenval, false);
  }

  {
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<complex<double>, General, ArrayRowSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, true, nb_eigenval);
  }

  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval);
  }

  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval);
  }

  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, General, ArrayRowSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval);
  }

  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> K;
    Matrix<complex<double>, General, ArrayRowSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval);
  }

  {
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    double x;
    TestGeneralProblem(K, M, x, true, nb_eigenval, false);
  }

  {
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, true, nb_eigenval);
  }

  {
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, General, ArrayRowSparse> M;
    double x;
    TestGeneralProblem(K, M, x, true, nb_eigenval, false);
  }

  {
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<complex<double>, General, ArrayRowSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, true, nb_eigenval);
  }
  
  if (all_test)
    cout << "All tests passed successfully" << endl;

  return FinalizeSeldon();
}
