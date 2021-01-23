
#ifndef SELDON_FILE_ARPACK_HXX
#define SELDON_FILE_ARPACK_HXX

// Headers from ARPACK.

// Modifications (by Lin Wu):
// Replacements:
//    integer       --> ARPACK_INTEGER
//    real          --> ARPACK_REAL
//    doublereal    --> ARPACK_DOUBLEREAL
//    complex       --> ARPACK_COMPLEX
//    doublecomplex --> ARPACK_DOUBLECOMPLEX
//    logical       --> ARPACK_LOGICAL
 

/*! Reverse communication interface for the Implicitly Restarted Arnoldi
  Iteration. */

// Arpack interface.
#undef ARPACK_INTEGER
#define ARPACK_INTEGER int
#undef ARPACK_REAL
#define ARPACK_REAL float
#undef ARPACK_DOUBLEREAL
#define ARPACK_DOUBLEREAL double
#undef ARPACK_COMPLEX
#define ARPACK_COMPLEX void
#undef ARPACK_DOUBLECOMPLEX
#define ARPACK_DOUBLECOMPLEX void
#undef ARPACK_LOGICAL
#define ARPACK_LOGICAL int

extern "C"
{
#include "carpack.h"
}

namespace Seldon
{

#ifdef SELDON_WITH_MPI
  void saupd(ARPACK_INTEGER comm, ARPACK_INTEGER& ido, char bmat,
	     ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_DOUBLEREAL& tol, ARPACK_DOUBLEREAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
	     ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void saupd(ARPACK_INTEGER comm, ARPACK_INTEGER& ido, char bmat,
	     ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_REAL& tol, ARPACK_REAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
	     ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void seupd(ARPACK_INTEGER comm,
	     ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
	     ARPACK_DOUBLEREAL* d, ARPACK_DOUBLEREAL* Z,
	     ARPACK_INTEGER ldz, ARPACK_DOUBLEREAL sigma, char bmat,
	     ARPACK_INTEGER n, char* which, ARPACK_INTEGER nev,
	     ARPACK_DOUBLEREAL tol, ARPACK_DOUBLEREAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
	     ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void seupd(ARPACK_INTEGER comm, ARPACK_LOGICAL rvec,
	     char HowMny, ARPACK_LOGICAL *select,
	     ARPACK_REAL* d, ARPACK_REAL* Z, ARPACK_INTEGER ldz,
	     ARPACK_REAL sigma, char bmat, ARPACK_INTEGER n,
	     char* which, ARPACK_INTEGER nev, ARPACK_REAL tol,
	     ARPACK_REAL* resid, ARPACK_INTEGER ncv, ARPACK_REAL* V,
	     ARPACK_INTEGER ldv, ARPACK_INTEGER* iparam,
	     ARPACK_INTEGER* ipntr, ARPACK_REAL* workd, ARPACK_REAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void naupd(ARPACK_INTEGER comm,
	     ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_DOUBLEREAL& tol,
	     ARPACK_DOUBLEREAL* resid, ARPACK_INTEGER ncv,
	     ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv, ARPACK_INTEGER* iparam,
	     ARPACK_INTEGER* ipntr, ARPACK_DOUBLEREAL* workd,
	     ARPACK_DOUBLEREAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void naupd(ARPACK_INTEGER comm,
	     ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_REAL& tol, ARPACK_REAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
	     ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void neupd(ARPACK_INTEGER comm,
	     ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
	     ARPACK_DOUBLEREAL* dr, ARPACK_DOUBLEREAL* di,
	     ARPACK_DOUBLEREAL* Z, ARPACK_INTEGER ldz, ARPACK_DOUBLEREAL sigmar,
	     ARPACK_DOUBLEREAL sigmai, ARPACK_DOUBLEREAL* workv, char bmat,
	     ARPACK_INTEGER n, char* which, ARPACK_INTEGER nev,
	     ARPACK_DOUBLEREAL tol, ARPACK_DOUBLEREAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
	     ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void neupd(ARPACK_INTEGER comm,
	     ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
	     ARPACK_REAL* dr, ARPACK_REAL* di, ARPACK_REAL* Z,
	     ARPACK_INTEGER ldz, ARPACK_REAL sigmar, ARPACK_REAL sigmai,
	     ARPACK_REAL* workv, char bmat, ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_REAL tol, ARPACK_REAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
	     ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
#else
  void saupd(ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_DOUBLEREAL& tol, ARPACK_DOUBLEREAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
	     ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void saupd(ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_REAL& tol, ARPACK_REAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
	     ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void seupd(ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
	     ARPACK_DOUBLEREAL* d, ARPACK_DOUBLEREAL* Z,
	     ARPACK_INTEGER ldz, ARPACK_DOUBLEREAL sigma, char bmat,
	     ARPACK_INTEGER n, char* which, ARPACK_INTEGER nev,
	     ARPACK_DOUBLEREAL tol, ARPACK_DOUBLEREAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
	     ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void seupd(ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
	     ARPACK_REAL* d, ARPACK_REAL* Z, ARPACK_INTEGER ldz,
	     ARPACK_REAL sigma, char bmat, ARPACK_INTEGER n,
	     char* which, ARPACK_INTEGER nev, ARPACK_REAL tol,
	     ARPACK_REAL* resid, ARPACK_INTEGER ncv, ARPACK_REAL* V,
	     ARPACK_INTEGER ldv, ARPACK_INTEGER* iparam,
	     ARPACK_INTEGER* ipntr, ARPACK_REAL* workd, ARPACK_REAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void naupd(ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_DOUBLEREAL& tol,
	     ARPACK_DOUBLEREAL* resid, ARPACK_INTEGER ncv,
	     ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv, ARPACK_INTEGER* iparam,
	     ARPACK_INTEGER* ipntr, ARPACK_DOUBLEREAL* workd,
	     ARPACK_DOUBLEREAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void naupd(ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_REAL& tol, ARPACK_REAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
	     ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void neupd(ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
	     ARPACK_DOUBLEREAL* dr, ARPACK_DOUBLEREAL* di,
	     ARPACK_DOUBLEREAL* Z, ARPACK_INTEGER ldz, ARPACK_DOUBLEREAL sigmar,
	     ARPACK_DOUBLEREAL sigmai, ARPACK_DOUBLEREAL* workv, char bmat,
	     ARPACK_INTEGER n, char* which, ARPACK_INTEGER nev,
	     ARPACK_DOUBLEREAL tol, ARPACK_DOUBLEREAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
	     ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
	     ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
  
  void neupd(ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
	     ARPACK_REAL* dr, ARPACK_REAL* di, ARPACK_REAL* Z,
	     ARPACK_INTEGER ldz, ARPACK_REAL sigmar, ARPACK_REAL sigmai,
	     ARPACK_REAL* workv, char bmat, ARPACK_INTEGER n, char* which,
	     ARPACK_INTEGER nev, ARPACK_REAL tol, ARPACK_REAL* resid,
	     ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
	     ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
	     ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info);
#endif
  
#ifdef SELDON_WITH_VIRTUAL
  template<class T>
  class EigenProblem_Base;
  
  template<class T>
  void FindEigenvaluesArpack(EigenProblem_Base<T>& var,
                             Vector<T>& eigen_values,
                             Vector<T>& lambda_imag,
                             Matrix<T, General, ColMajor>& eigen_vectors);
#else  
  // main function to find eigenvalues and eigenvectors with Arpack
  template<class EigenProblem, class T, class Allocator1,
	   class Allocator2, class Allocator3>
  void FindEigenvaluesArpack(EigenProblem& var,
                             Vector<T, VectFull, Allocator1>& eigen_values,
                             Vector<T, VectFull, Allocator2>& lambda_imag,
                             Matrix<T, General, ColMajor, Allocator3>& eigen_vectors);
#endif
  
  
  
  template<class T, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Alloc5>
  void CallArpack(int comm, int& ido, char& bmat, int& n, string& which, int& nev,
		  T& tol, Vector<T, VectFull, Allocator1>& resid,
		  int& ncv, Matrix<T, General, ColMajor, Allocator2>& v,
		  int& ldv, Vector<int>& iparam, Vector<int>& ipntr, int sym,
		  Vector<T, VectFull, Allocator3>& workd,
		  Vector<T, VectFull, Allocator4>& workl,
		  int& lworkl, Vector<T, VectFull, Alloc5>& rwork, int& info);

  template<class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Alloc5>
  void CallArpack(int comm, int& ido, char& bmat, int& n,
		  string& which, int& nev, double& tol,
		  Vector<complex<double>, VectFull, Allocator1>& resid, int& ncv,
		  Matrix<complex<double>, General, ColMajor, Allocator2>& v,
		  int& ldv, Vector<int>& iparam, Vector<int>& ipntr, int sym,
		  Vector<complex<double>, VectFull, Allocator3>& workd,
		  Vector<complex<double>, VectFull, Allocator4>& workl,
		  int& lworkl, Vector<double, VectFull, Alloc5>& rwork, int& info);

  template<class T, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Allocator5,
	   class Allocator6, class Allocator7, class Alloc8>
  void CallArpack(int comm, int& rvec, char& howmny, Vector<int>& selec,
		  Vector<T, VectFull, Allocator1>& lambda,
		  Vector<T, VectFull, Allocator2>& lambda_i,
		  Matrix<T, General, ColMajor, Allocator3>& eigen_vec,
		  int& ldz, T& shiftr, T& shifti, char& bmat, int& n,
		  string& which, int& nev, T& tol,
		  Vector<T, VectFull, Allocator4>& resid,
		  int& ncv, Matrix<T, General, ColMajor, Allocator5>&v,
		  int& ldv, Vector<int>& iparam, Vector<int>& ipntr,
		  int sym, Vector<T, VectFull, Allocator6>& workd,
		  Vector<T, VectFull, Allocator7>& workl,
		  int& lworkl, Vector<T, VectFull, Alloc8>& rwork, int& info);

  template<class Allocator1, class Allocator2, class Allocator3,
	   class Allocator4, class Allocator5, class Allocator6,
	   class Allocator7, class Alloc8>
  void CallArpack(int comm, int& rvec, char& howmny, Vector<int>& selec,
		  Vector<complex<double>, VectFull, Allocator1>& lambda,
		  Vector<complex<double>, VectFull, Allocator2>& lambda_i,
		  Matrix<complex<double>, General, ColMajor, Allocator3>& eigen_vectors,
		  int& ldz, complex<double>& shiftr, complex<double>& shifti,
		  char& bmat, int& n, string& which, int& nev, double& tol,
		  Vector<complex<double>, VectFull, Allocator4>& resid,
		  int& ncv, Matrix<complex<double>, General, ColMajor, Allocator5>& v,
		  int& ldv, Vector<int>& iparam, Vector<int>& ipntr, int sym,
		  Vector<complex<double>, VectFull, Allocator6>& workd,
		  Vector<complex<double>, VectFull, Allocator7>& workl,
		  int& lworkl, Vector<double, VectFull, Alloc8>& rwork, int& info);
  
}


#endif
