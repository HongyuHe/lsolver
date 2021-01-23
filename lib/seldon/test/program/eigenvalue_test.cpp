// C library for time function and for randomization
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <stdint.h>

#include "SeldonLib.hxx"

using namespace Seldon;

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
    x = complex<T>(0, rand())/T(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/T(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/T(RAND_MAX);
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

int main(int argc, char** argv)
{
#ifdef SELDON_WITH_MPI
  MPI_Init(&argc, &argv);
#endif
  
  cout.precision(15);
  
  {
    // simple example with a symmetric matrix
    Vector<double> LambdaRealRef, LambdaImagRef;
    Matrix<double, Symmetric, ArrayRowSymSparse> A;
    Matrix<double, Symmetric, RowSymPacked> Adense;
    
    // generation of the sparse symmetric matrix
    int n = 40, nnz = 4*n;  
    GenerateRandomMatrix(A, n, n, nnz);
    
    A.WriteText("A.dat");
    
    Copy(A, Adense);
    
    // the reference eigenvalues and eigenvectors are computed
    Matrix<double> EigenvecRef;
    GetEigenvaluesEigenvectors(Adense, LambdaRealRef, EigenvecRef);
    
    // setting eigenvalue problem
    int nb_eigenval = 5;
    SparseEigenProblem<double, Matrix<double, Symmetric, ArrayRowSymSparse> > var_eig;
    var_eig.SetStoppingCriterion(1e-12);
    var_eig.SetNbAskedEigenvalues(nb_eigenval);
    
    // by default large eigenvalues are researched
    // if you want to search small eigenvalues or eigenvalues close to a given
    // value, you should call var_eig.SetTypeSpectrum(...)
    
#ifdef SELDON_WITH_ARPACK
    // Arpack parameters
    // Regular mode uses only matrix vector products to obtain large eigenvalues
    // to obtain eigenvalues close to a given value, SHIFTED_MODE should be used
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    
    // if you want to change the maximal number of iterations
    // var_eig.SetNbMaximumIterations(20000);
#endif
    
#ifdef SELDON_WITH_ANASAZI
    // Anasazi parameters
    // Regular mode uses only matrix vector products to obtain large eigenvalues
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    
    // BKS (Block Krylov Schur) solver should be used for non-symmetric system
    //var_eig.SetEigensolverType(var_eig.SOLVER_BKS);

    // for symmetric system, you can use SOLVER_LOBPCG or SOLVER_BD
    var_eig.SetEigensolverType(var_eig.SOLVER_LOBPCG);
    
    // used for BD or BKS
    var_eig.SetNbBlocks(4);

    // size of a block
    var_eig.SetNbArnoldiVectors(7);

    // you can change the maximal number of iterations and restarts
    //var_eig.SetNbMaximumRestarts(200);
    //var_eig.SetNbMaximumIterations(20000);
#endif
    
#ifdef SELDON_WITH_FEAST
    var_eig.SetPrintLevel(1);
    // Feast parameters
    // the eigenvalues are searched in the interval [a, b]
    // you must provide a and b
    // as a result large eigenvalues are not reachable with feast
    var_eig.SetIntervalSpectrum(0.1, 0.3);
#endif
    
    // eigenproblem is A x = lambda x
    var_eig.InitMatrix(A);

    // if you have a generalized eigenproblem, you can write
    // var_eig.InitMatrix(A, B);
    // in that case the eigenproblem is A x = lambda B x

    // eigenvalues and eigenvectors
    Vector<double> lambda, lambda_imag;
    Matrix<double, General, ColMajor> eigen_vec;
    
    // main function 
    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);

    // eigenvalues are written
    DISP(lambda); DISP(lambda_imag);
    DISP(LambdaRealRef); DISP(LambdaImagRef);
    
    // checking that eigenvalues are correct
    for (int i = 0; i < lambda.GetM(); i++)
      {
	int j = -1;
	for (int k = 0; k < A.GetM(); k++)
	  if (abs(lambda(i) - LambdaRealRef(k)) <= 1e-10)
	    j = k;
	
	Vector<double> x(A.GetM()), xref(A.GetM());
	for (int k = 0; k < A.GetM(); k++)
	  {
	    x(k) = eigen_vec(k, i);
	    xref(k) = EigenvecRef(k, j);
	  }
	
	double norm_diff = 0.0;
	if (x(0)*xref(0) < 0)
	  for (int k = 0; k < A.GetM(); k++)
	    norm_diff = max(norm_diff, abs(x(k) + xref(k)));
	else
	  for (int k = 0; k < A.GetM(); k++)
	    norm_diff = max(norm_diff, abs(x(k) - xref(k)));
	
	DISP(i); DISP(norm_diff);
      }   
  }
  
#ifdef SELDON_WITH_MPI
  MPI_Finalize();
#endif

  return 0;
}

