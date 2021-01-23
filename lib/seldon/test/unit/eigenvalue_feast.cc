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
  CheckEigenvalues(mat_stiff, lambda, eigen_vec);
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

template<class MatrixS, class MatrixM, class Matrix1>
void CheckEigenvalues(const MatrixS& mat_stiff,
                      const MatrixM& mat_mass,
                      const Vector<complex<double> >& lambda,
                      const Vector<complex<double> >& lambda_imag, const Matrix1& eigen_vec)
{
  CheckEigenvalues(mat_stiff, mat_mass, lambda, eigen_vec);
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
void TestSymmetricProblem(MatrixK& K, MatrixM& M, T&,
                          bool add_pos_diag, int nb_eigenval, bool test_chol)
{
  int type_solver = TypeEigenvalueSolver::FEAST;
  Vector<T>  lambda, lambda_imag;
  Matrix<T, General, ColMajor> eigen_vec;
  
  int m = 40, n = m, nnz = 4*m;
  GenerateRandomMatrix(K, m, n, nnz);
  GenerateRandomMatrix(M, m, n, nnz);
  
  // on rajoute des coefs sur la diagonale de M
  // pour avoir une matrice definie positive
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
  var_eig.InitMatrix(K);  
  
  // eigenvalues in a given interval
  var_eig.SetIntervalSpectrum(-1.7, -1.0);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec, type_solver);
  DISP(lambda); DISP(lambda_imag);
  
  CheckEigenvalues(K, lambda, eigen_vec);

  // case with a diagonal matrix M    
  MatrixM Mdiag(n, n);
  Vector<typename MatrixM::entry_type> D;
  GenerateRandomVector(D, n);
  
  D.Write("D.dat");
  
  for (int i = 0; i < n; i++)
    Mdiag.AddInteraction(i, i, D(i));
  
  var_eig.InitMatrix(K, Mdiag);
  var_eig.SetDiagonalMass();

  // eigenvalues in a given interval
  var_eig.SetIntervalSpectrum(0.1, 1.4);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec, type_solver);
  DISP(lambda); DISP(lambda_imag);
  
  CheckEigenvalues(K, Mdiag, lambda, eigen_vec);
  
  // case with a symmetric positive definite mass matrix
  if (test_chol)
    {
      var_eig.InitMatrix(K, M);
      var_eig.SetCholeskyFactoForMass();
      
      var_eig.SetIntervalSpectrum(-0.2, -0.04);
      
      GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec, type_solver);
      DISP(lambda); DISP(lambda_imag);
      
      CheckEigenvalues(K, M, lambda, eigen_vec);
    }
}


template<class MatrixK, class MatrixM, class T>
void TestGeneralProblem(MatrixK& K, MatrixM& M, T&,
                        bool add_pos_diag, int nb_eigenval, bool test_chol)
{
  int type_solver = TypeEigenvalueSolver::FEAST;
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
  var_eig.SetPrintLevel(1);
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(nb_eigenval);  
  var_eig.InitMatrix(K);  
  
  var_eig.SetCircleSpectrum(complex<double>(1.2, 0.4), 1.0);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec, type_solver);
  DISP(lambda); DISP(lambda_imag);
  
  CheckEigenvalues(K, lambda, lambda_imag, eigen_vec);
  
  // then we test with the mass matrix
  var_eig.InitMatrix(K, M);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec, type_solver);
  DISP(lambda); DISP(lambda_imag);
  
  CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);
  
  return;
  

}

class HermitianEigenProblem : public SparseEigenProblem<complex<double>, Matrix<complex<double>, General, ArrayRowSparse>, Matrix<complex<double>, General, ArrayRowSparse> >
{
public:
  bool IsHermitianProblem() const { return true; }

};


template<class MatrixK, class MatrixM>
void TestHermitianProblem(MatrixK& K, MatrixM& M, int nb_eigenval)
{
  int type_solver = TypeEigenvalueSolver::FEAST;
  Vector<complex<double> >  lambda, lambda_imag;
  Matrix<complex<double>, General, ColMajor> eigen_vec;
  
  int m = 40, n = m, nnz = 4*m;
  GenerateRandomMatrix(K, m, n, nnz);
  GenerateRandomMatrix(M, m, n, nnz);
  
  MatrixK K2 = K; MatrixM M2 = M;
  
  // on rend les matrices hermitiennes
  for (int i = 0; i < M2.GetM(); i++)
    for (int j = 0; j < M2.GetRowSize(i); j++)
      M.AddInteraction(M2.Index(i, j), i, conjugate(M2.Value(i, j)));

  for (int i = 0; i < K2.GetM(); i++)
    for (int j = 0; j < K2.GetRowSize(i); j++)
      K.AddInteraction(K2.Index(i, j), i, conjugate(K2.Value(i, j)));

  // on rajoute des coefs sur la diagonale de M
  // pour avoir une matrice definie positive
  double coef_diag;
  for (int i = 0; i < M.GetM(); i++)
    {
      GetRandNumber(coef_diag);
      coef_diag += 3.0;
      
      M.AddInteraction(i, i, coef_diag);
    }
  
  K.WriteText("K.dat");
  M.WriteText("M.dat");
  
  HermitianEigenProblem var_eig;
  
  // first, we test the case without the mass matrix M
  
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(nb_eigenval);  
  var_eig.InitMatrix(K);  
  var_eig.SetPrintLevel(1);
  
  // eigenvalues in a given interval
  var_eig.SetIntervalSpectrum(-2.5, -1.0);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec, type_solver);
  DISP(lambda); DISP(lambda_imag);
  
  CheckEigenvalues(K, lambda, eigen_vec);
  
  // generalized problem
  var_eig.InitMatrix(K, M);  
  var_eig.SetIntervalSpectrum(-0.8, -0.2);
  
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec, type_solver);
  DISP(lambda); DISP(lambda_imag);
  
  CheckEigenvalues(K, M, lambda, eigen_vec);
}



int main(int argc, char** argv)
{
  InitSeldon(argc, argv);
  
  bool all_test = true;

  srand(0);
  
  // testing non-symmetric solver
  {
    int nb_eigenval = 20;
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<double, General, ArrayRowSparse> M;
    double x;
    TestGeneralProblem(K, M, x, false, nb_eigenval, false);
  }
  
  srand(0);

  {
    int nb_eigenval = 20;
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<complex<double>, General, ArrayRowSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval, false);
  }
  
  srand(0);
  {
    int nb_eigenval = 20;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> K;
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> M;
    complex<double> x;
    TestGeneralProblem(K, M, x, false, nb_eigenval, false);
  }
 
  srand(0);
  
  // testing symmetric solver
  {
    int nb_eigenval = 6;
    Matrix<double, Symmetric, ArrayRowSymSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    double x;
    TestSymmetricProblem(K, M, x, true, nb_eigenval, true);
  }

  srand(0);

  {
    int nb_eigenval = 12;
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<complex<double>, General, ArrayRowSparse> M;
    TestHermitianProblem(K, M, nb_eigenval);
  }
    
  if (all_test)
    cout << "All tests passed successfully" << endl;

  return FinalizeSeldon();
}
