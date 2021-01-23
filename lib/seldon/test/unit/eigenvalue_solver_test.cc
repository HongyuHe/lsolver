#include "SeldonLib.hxx"

using namespace Seldon;

double threshold = 1e-6;

/******************************
 * Matrix defined by the user *
 ******************************/


// basic class defining a matrix
template<class T>
class Matrix_Laplacian1D : public VirtualMatrix<T>
{
protected :
  double L;
  
public :
  double dx;
  
  void Init(int n_, double L_)
  {
    this->m_ = n_;
    this->n_ = n_;
    L = L_;
    dx = L/(n_+1);
  }

  // matrix vector product Y = A*X
  void MltVector(const Vector<T>& X, Vector<T>& Y) const
  {
    int n = this->n_;
    Y(0) = 2.0*X(0) - X(1);
    Y(n-1) = 2.0*X(n-1) - X(n-2);
    for (int i = 1; i < n-1; i++)
      Y(i) = 2.0*X(i) - X(i-1) - X(i+1);
    
    Mlt(1.0/(dx*dx), Y);
  }
  
  bool IsSymmetric() const { return true; }
  
};


// returns true if the matrix is complex
template<class T>
bool IsComplexMatrix(const Matrix_Laplacian1D<T>& A)
{
  return false;
}


// returns true if the matrix is symmetric
template<class T>
bool IsSymmetricMatrix(const Matrix_Laplacian1D<T>& A)
{
  return true;
}


/************************
 * Checking eigenvalues *
 ************************/


// checking symmetric standard eigenproblems
template<class Matrix1>
void CheckEigenvalues(const VirtualMatrix<double>& mat_stiff,
                      const Vector<double>& lambda, const Matrix1& eigen_vec)
{
  int N = mat_stiff.GetM();
  Vector<double> X(N), Y(N);
  X.Fill(0);
  Y.Fill(0);
  for (int i = 0; i < lambda.GetM(); i++)
    {
      for (int j = 0; j < N; j++)
        X(j) = eigen_vec(j, i);
      
      mat_stiff.MltVector(X, Y);
      double err = 0;
      double normeX = sqrt(DotProd(X, X));
      for (int j = 0; j < N; j++)
        err += pow(Y(j) - lambda(i)*X(j), 2.0);
      
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


// checking complex standard eigenproblems
template<class Prop, class Storage, class Matrix1>
void CheckEigenvalues(const Matrix<complex<double>, Prop, Storage>& mat_stiff,
                      const Vector<complex<double> >& lambda, const Matrix1& eigen_vec)
{
  int N = mat_stiff.GetM();
  Vector<complex<double> > X(N), Y(N);
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
          cout << "Error = " << err << endl;
          abort();
        }
    }
  
  
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


// checking symmetric generalized eigenproblems
template<class Prop0, class Storage0,
         class Prop1, class Storage1, class Matrix1>
void CheckEigenvalues(const Matrix<double, Prop0, Storage0>& mat_stiff,
                      const Matrix<double, Prop1, Storage1>& mat_mass,
                      const Vector<double>& lambda, const Matrix1& eigen_vec)
{
  int N = mat_stiff.GetM();
  Vector<double> X(N), Y(N), Mx(N);
  X.Fill(0);
  Y.Fill(0);
  Mx.Fill(0);
  for (int i = 0; i < lambda.GetM(); i++)
    {
      for (int j = 0; j < N; j++)
        X(j) = eigen_vec(j, i);
      
      Mlt(mat_stiff, X, Y);
      Mlt(mat_mass, X, Mx);
      double err = 0;
      double normeX = sqrt(abs(DotProdConj(X, X)));
      for (int j = 0; j < N; j++)
        err += pow(Y(j) - lambda(i)*Mx(j), 2.0);
      
      err = sqrt(err);
      if (err > threshold*normeX)
        {
          cout << "Error on eigenvalue " << lambda(i) << endl;
          cout << "Error = " << err/normeX << endl;
          abort();
        }
    }  
}


// checking complex generalized eigenproblems
template<class Prop0, class Storage0,
         class Prop1, class Storage1, class Matrix1>
void CheckEigenvalues(const Matrix<complex<double>, Prop0, Storage0>& mat_stiff,
                      const Matrix<double, Prop1, Storage1>& mat_mass,
                      const Vector<complex<double> >& lambda, const Matrix1& eigen_vec)
{
  int N = mat_stiff.GetM();
  Vector<complex<double> > X(N), Y(N), Mx(N);
  X.Fill(0);
  Y.Fill(0);
  Mx.Fill(0);
  for (int i = 0; i < lambda.GetM(); i++)
    {
      for (int j = 0; j < N; j++)
        X(j) = eigen_vec(j, i);
      
      Mlt(mat_stiff, X, Y);
      Mlt(mat_mass, X, Mx);
      double err = 0;
      double normeX = sqrt(abs(DotProdConj(X, X)));
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

int main(int argc, char** argv)
{
  InitSeldon(argc, argv);
  
  bool all_test = true;
  int nb_eigenval = 10;
  
  {
    // testing dense real symmetric eigenvalue problems
    Matrix<double, Symmetric, RowSymPacked> mat_stiff, K;
    Matrix<double, Symmetric, RowSymPacked> mat_mass, M;
    
    mat_stiff.ReadText("test/matrix/MatStiffPyramidHcurl.dat");
    mat_mass.ReadText("test/matrix/MatMassPyramidHcurl.dat");
    K = mat_stiff; M = mat_mass;
    
    // setting eigenvalue problem
    DenseEigenProblem<double, double, Symmetric, RowSymPacked> var_eig;
    var_eig.SetStoppingCriterion(1e-14);
    var_eig.SetNbAskedEigenvalues(nb_eigenval);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    
    // finding large eigenvalues of K
    var_eig.InitMatrix(mat_stiff);
    Vector<double> lambda, lambda_imag;
    Matrix<double, General, ColMajor> eigen_vec;
    
    cout << "Testing computation of large eigenvalues of symmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);
    
    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues close to 1e-3
    var_eig.InitMatrix(mat_stiff);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 1e-3, var_eig.SORTED_REAL);
    
    cout << "Testing computation of clustered eigenvalues of symmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda);
        
    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding large eigenvalues of M^-1 K
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_REAL);
    
    cout << "Testing computation of large eigenvalues of generalized symmetric eigenproblem..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues of M^-1 K close to 5
    mat_mass = M;
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.INVERT_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 5.0, var_eig.SORTED_REAL);

    cout << "Testing computation of clustered eigenvalues of generalized symmetric eigenproblem..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues of M^-1 K close to 100
    mat_mass = M;
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 100.0, var_eig.SORTED_REAL);

    cout << "Testing computation of clustered eigenvalues of generalized symmetric eigenproblem (shift-invert mode)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues of M^-1 K close to 100
    mat_mass = M;
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.BUCKLING_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 100.0, var_eig.SORTED_REAL);

    cout << "Testing computation of clustered eigenvalues of generalized symmetric eigenproblem (Buckling mode)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues of M^-1 K close to 5
    mat_mass = M;
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.CAYLEY_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 5.0, var_eig.SORTED_REAL);

    cout << "Testing computation of clustered eigenvalues of generalized symmetric eigenproblem (Cayley mode)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // use of Cholesky factorisation to reduce
    // generalized eigenvalue problem to standard one
    mat_mass = M;
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetCholeskyFactoForMass(true);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_REAL);
    
    cout << "Testing computation of large eigenvalues of generalized symmetric eigenproblem (Cholesky facto)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    mat_mass = M;
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetCholeskyFactoForMass(true);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 50, var_eig.SORTED_REAL);

    cout << "Testing computation of clustered eigenvalues of generalized symmetric eigenproblem (Cholesky facto)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    

    var_eig.SetCholeskyFactoForMass(false);        
    // standard problem too if mass matrix is diagonal
    mat_stiff.ReadText("test/matrix/MatStiffHexaHcurl.dat");
    mat_mass.ReadText("test/matrix/MatMassHexaHcurl.dat");
    K = mat_stiff; M = mat_mass;
    
    mat_mass = M;
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetDiagonalMass(true);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_REAL);
    cout << "Testing computation of large eigenvalues of generalized symmetric eigenproblem (diagonal mass)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    mat_mass = M;
    mat_stiff = K;
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetDiagonalMass(true);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 50, var_eig.SORTED_REAL);
    
    cout << "Testing computation of clustered eigenvalues of generalized symmetric eigenproblem (diagonal mass)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    DISP(lambda);
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
  }
  
  {
    cout << endl << endl << "Unsymmetric dense matrices " << endl;
    // testing dense real unsymmetric eigenvalue problem
    Matrix<double, General, RowMajor> mat_stiff;
    Matrix<double, Symmetric, RowSymPacked> mat_mass;
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    
    K.ReadText("test/matrix/MatStiffDenseUnsym.dat");
    M.ReadText("test/matrix/MatMassDenseUnsym.dat");
    Copy(K, mat_stiff);
    Copy(M, mat_mass);
    
    // setting eigenvalue problem
    DenseEigenProblem<double, double, General, RowMajor> var_eig;
    var_eig.SetStoppingCriterion(1e-12);
    var_eig.SetNbAskedEigenvalues(nb_eigenval);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    
    // finding large eigenvalues of K
    var_eig.InitMatrix(mat_stiff);
    Vector<double> lambda, lambda_imag;
    Matrix<double, General, ColMajor> eigen_vec;
    
    cout << "Testing computation of large eigenvalues of unsymmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, lambda, lambda_imag, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues closest to a complex shift
    Copy(K, mat_stiff);
    Copy(M, mat_mass);    
    var_eig.InitMatrix(mat_stiff);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, complex<double>(0, 0.104), var_eig.SORTED_MODULUS);
    
    cout << "Testing computation of clustered eigenvalues of unsymmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, lambda, lambda_imag, eigen_vec);
    cout << endl << endl;
    
    // finding large eigenvalues of M^-1 K
    Copy(K, mat_stiff);
    Copy(M, mat_mass);
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);

    cout << "Testing computation of large eigenvalues of generalized unsymmetric eigenproblem..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues of M^-1 K close to a complex shift
    Copy(K, mat_stiff);
    Copy(M, mat_mass);    
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, complex<double>(0, 0.0121), var_eig.SORTED_MODULUS);

    cout << "Testing computation of clustered eigenvalues of generalized unsymmetric eigenproblem..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);
    cout << endl << endl;
  }
  
  {
    cout << endl << endl << "Unsymmetric complex matrices " << endl;
    // testing dense complex unsymmetric eigenvalue problem
    Matrix<complex<double>, General, RowMajor> mat_stiff, K, Mc;
    Matrix<double, Symmetric, RowSymPacked> mat_mass, M;
    
    K.Read("test/matrix/MatStiffDenseComplex.dat");
    Mc.Read("test/matrix/MatMassDenseComplex.dat");
    Copy(K, mat_stiff);
    M.Reallocate(Mc.GetM(), Mc.GetM());
    for (int i = 0; i < Mc.GetM(); i++)
      for (int j = i; j < Mc.GetM(); j++)
        M(i, j) = abs(Mc(i, j));
    
    // setting eigenvalue problem
    DenseEigenProblem<complex<double>, complex<double>, General, RowMajor> var_eig;
    var_eig.SetStoppingCriterion(1e-12);
    var_eig.SetNbAskedEigenvalues(nb_eigenval);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    
    // finding large eigenvalues of K
    var_eig.InitMatrix(mat_stiff);
    Vector<complex<double> >  lambda, lambda_imag;
    Matrix<complex<double>, General, ColMajor> eigen_vec;

    cout << "Testing computation of large eigenvalues of complex matrix ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda);
    
    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues closest to a complex shift
    Copy(K, mat_stiff);
    Copy(M, mat_mass);    
    var_eig.InitMatrix(mat_stiff);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES,
                            complex<double>(0.01, 1e-4), var_eig.SORTED_MODULUS);

    cout << "Testing computation of clustered eigenvalues of complex matrix ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda);
    
    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding large eigenvalues of M^-1 K
    Copy(K, mat_stiff);
    Copy(M, mat_mass);    
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);

    cout << "Testing computation of large eigenvalues of generalized complex eigenproblem ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;

    // finding eigenvalues of M^-1 K close to a complex shift
    Copy(K, mat_stiff);
    Copy(M, mat_mass);    
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, complex<double>(1.0,1e-3), var_eig.SORTED_MODULUS);

    cout << "Testing computation of clustered eigenvalues of generalized complex eigenproblem ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // of course, you can also use shift-invert mode to reduce the problem to a standard one
    Copy(K, mat_stiff);
    Copy(M, mat_mass);    
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.INVERT_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, complex<double>(1.0,1e-3), var_eig.SORTED_MODULUS);

    cout << "Testing computation of clustered eigenvalues of generalized complex eigenproblem ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
  }
  
  {
    cout << endl << endl << "Symmetric sparse matrices " << endl;
    // sparse symmetric eigenvalue problems
    Matrix<double, Symmetric, ArrayRowSymSparse> K, M;
    Vector<double> lambda, lambda_imag;
    Matrix<double, General, ColMajor> eigen_vec;
    
    M.ReadText("test/matrix/SparseMassSymmetric.dat");
    K.ReadText("test/matrix/SparseStiffnessSymmetric.dat");
    
    SparseEigenProblem<double, Matrix<double, Symmetric, ArrayRowSymSparse> > var_eig;
    var_eig.SetStoppingCriterion(1e-12);
    var_eig.SetNbAskedEigenvalues(nb_eigenval);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    
    // finding large eigenvalues of K
    var_eig.InitMatrix(K);
    
    cout << "Testing computation of large eigenvalues of symmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);

    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;
    
    // then eigenvalues of K close to a shift
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 1.0);
    
    cout << "Testing computation of clustered eigenvalues of symmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);

    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;
    
    // large eigenvalues of M^-1 K
    var_eig.InitMatrix(K, M);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0.0);
    
    cout << "Testing computation of large eigenvalues of generalized symmetric "
         << "eigenproblem ..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);

    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;

    // clustered eigenvalues of M^-1 K
    var_eig.InitMatrix(K, M);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 100.0);
    
    cout << "Testing computation of clustered eigenvalues of generalized symmetric "
         << "eigenproblem ..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);

    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;

    var_eig.SetComputationalMode(var_eig.BUCKLING_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 100.0);
    
    cout << "Testing computation of clustered eigenvalues of generalized symmetric "
         << "eigenproblem (buckling) ..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);

    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;

    var_eig.SetComputationalMode(var_eig.CAYLEY_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 100.0);
    
    cout << "Testing computation of clustered eigenvalues of generalized symmetric "
         << "eigenproblem (cayley) ..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);

    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // we can exploit that mass matrix is diagonal (when it is the case)
    var_eig.SetDiagonalMass();
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 100.0);
    
    cout << "Testing computation of clustered eigenvalues of generalized symmetric "
         << "eigenproblem (diagonal)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);

    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // or use Cholesky factorisation if mass matrix is not diagonal
    K.ReadText("test/matrix/SparseStiffnessTriSymmetric.dat");
    M.ReadText("test/matrix/SparseMassTriSymmetric.dat");
    var_eig.InitMatrix(K, M);
    var_eig.SetDiagonalMass(false);
    var_eig.SetCholeskyFactoForMass(true);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 100.0);
    
    cout << "Testing computation of clustered eigenvalues of generalized symmetric "
         << "eigenproblem (cholesky)..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);

    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
  }
  
  {
    cout << endl << endl << "Unsymmetric sparse matrices " << endl;
    // sparse unsymmetric eigenvalue problem
    Matrix<double, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    
    K.ReadText("test/matrix/MatStiffDenseUnsym.dat");
    M.ReadText("test/matrix/MatMassDenseUnsym.dat");
    
    // setting eigenvalue problem
    SparseEigenProblem<double, Matrix<double, General, ArrayRowSparse> > var_eig;
    var_eig.SetStoppingCriterion(1e-12);
    var_eig.SetNbAskedEigenvalues(nb_eigenval);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    
    // finding large eigenvalues of K
    var_eig.InitMatrix(K);
    Vector<double> lambda, lambda_imag;
    Matrix<double, General, ColMajor> eigen_vec;
    
    cout << "Testing computation of large eigenvalues of unsymmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, lambda, lambda_imag, eigen_vec);
    cout << endl << endl;

    // finding eigenvalues closest to a complex shift
    var_eig.InitMatrix(K);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES,
                            complex<double>(0, 0.104), var_eig.SORTED_MODULUS);
    
    cout << "Testing computation of clustered eigenvalues of unsymmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, lambda, lambda_imag, eigen_vec);
    cout << endl << endl;
    
    // finding large eigenvalues of M^-1 K
    var_eig.InitMatrix(K, M);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);

    cout << "Testing computation of large eigenvalues of generalized unsymmetric eigenproblem..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues of M^-1 K close to a complex shift
    var_eig.InitMatrix(K, M);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES,
                            complex<double>(0, 0.0121), var_eig.SORTED_MODULUS);
    
    cout << "Testing computation of clustered eigenvalues of generalized unsymmetric eigenproblem..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, lambda_imag, eigen_vec);
    cout << endl << endl;
  }
    
  {
    cout << endl << endl << "Unsymmetric complex matrices " << endl;
    // testing dense complex unsymmetric eigenvalue problem
    Matrix<complex<double>, General, ArrayRowSparse> K;
    Matrix<double, Symmetric, ArrayRowSymSparse> M;
    
    K.ReadText("test/matrix/MatStiffDenseUnsym.dat");
    M.ReadText("test/matrix/MatMassDenseUnsym.dat");
    //K.ReadText("test/matrix/SparseComplexMass.dat");
    //M.ReadText("test/matrix/SparseComplexStiffness.dat");
    
    // setting eigenvalue problem
    SparseEigenProblem<complex<double>,
      Matrix<complex<double>, General, ArrayRowSparse> > var_eig;
    
    var_eig.SetStoppingCriterion(1e-12);
    var_eig.SetNbAskedEigenvalues(nb_eigenval);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    
    // finding large eigenvalues of K
    var_eig.InitMatrix(K);
    Vector<complex<double> >  lambda, lambda_imag;
    Matrix<complex<double>, General, ColMajor> eigen_vec;

    cout << "Testing computation of large eigenvalues of complex matrix ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda);
    
    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding eigenvalues closest to a complex shift
    var_eig.InitMatrix(K);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES,
                            complex<double>(0, 0.104), var_eig.SORTED_MODULUS);

    cout << "Testing computation of clustered eigenvalues of complex matrix ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda);
    
    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;
    
    // finding large eigenvalues of M^-1 K
    var_eig.InitMatrix(K, M);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);

    cout << "Testing computation of large eigenvalues of generalized complex eigenproblem ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;

    // finding eigenvalues of M^-1 K close to a complex shift
    var_eig.InitMatrix(K, M);
    var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, complex<double>(0, 0.0121), var_eig.SORTED_MODULUS);

    cout << "Testing computation of clustered eigenvalues of generalized complex eigenproblem ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
    
    // of course, you can also use shift-invert mode to reduce the problem to a standard one
    var_eig.InitMatrix(K, M);
    var_eig.SetComputationalMode(var_eig.INVERT_MODE);
    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, complex<double>(0, 0.0121), var_eig.SORTED_MODULUS);

    cout << "Testing computation of clustered eigenvalues of generalized complex eigenproblem ... " << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    DISP(lambda); DISP(lambda_imag);
    
    CheckEigenvalues(K, M, lambda, eigen_vec);
    cout << endl << endl;
  }

  {
    // testing matrix-free class (defined by the user)
    Matrix_Laplacian1D<double> K;
    K.Init(200, 2.0);
    
    // setting eigenvalue problem
    VirtualEigenProblem<double, double, double> var_eig;
    var_eig.SetStoppingCriterion(1e-12);
    var_eig.SetNbAskedEigenvalues(nb_eigenval);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.SMALL_EIGENVALUES, 0, var_eig.SORTED_MODULUS);

    // finding large eigenvalues of K
    var_eig.InitMatrix(K);
    Vector<double> lambda, lambda_imag;
    Matrix<double, General, ColMajor> eigen_vec;
    
    cout << "Testing computation of large eigenvalues of symmetric matrix..." << endl;
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);    
    DISP(lambda);
    
    CheckEigenvalues(K, lambda, eigen_vec);
    cout << endl << endl;

  }
  
  if (all_test)
    cout << "All tests passed successfully" << endl;

  return FinalizeSeldon();
}
