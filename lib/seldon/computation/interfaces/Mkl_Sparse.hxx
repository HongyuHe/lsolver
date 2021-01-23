#ifndef SELDON_FILE_MKL_SPARSE_HXX

extern "C"
{
  float  cblas_sdoti(const int N, const float *X, const int *indx,
		     const float *Y);
  
  double cblas_ddoti(const int N, const double *X, const int *indx,
		     const double *Y);
  
  void cblas_cdotui_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);
  
  void cblas_cdotci_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);

  void cblas_zdotui_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);

  void cblas_zdotci_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);
  
  void cblas_saxpby(const int N, const float alpha, const float *X,
                    const int incX, const float beta, const float *Y,
                    const int incY);

  void cblas_daxpby(const int N, const double alpha, const double *X,
                    const int incX, const double beta, const double *Y,
                    const int incY);
  
  void cblas_caxpby(const int N, const void *alpha, const void *X,
                    const int incX, const void *beta, const void *Y,
                    const int incY);

  void cblas_zaxpby(const int N, const void *alpha, const void *X,
                    const int incX, const void *beta, const void *Y,
                    const int incY);
  
  void cblas_saxpyi(const int N, const float alpha, const float *X,
		    const int *indx, float *Y);

  void cblas_daxpyi(const int N, const double alpha, const double *X,
		    const int *indx, double *Y);
  
  void cblas_caxpyi(const int N, const void *alpha, const void *X,
		    const int *indx, void *Y);
  
  void cblas_zaxpyi(const int N, const void *alpha, const void *X,
		    const int *indx, void *Y);  
  
  void cblas_sgthr(const int N, const float* y, float* x, const int* indx);

  void cblas_dgthr(const int N, const double* y, double* x, const int* indx);

  void cblas_cgthr(const int N, const void* y, void* x, const int* indx);

  void cblas_zgthr(const int N, const void* y, void* x, const int* indx);

  void cblas_sgthrz(const int N, float* y, float* x, const int* indx);

  void cblas_dgthrz(const int N, double* y, double* x, const int* indx);

  void cblas_cgthrz(const int N, void* y, void* x, const int* indx);

  void cblas_zgthrz(const int N, void* y, void* x, const int* indx);
  
  void cblas_sroti(const int N, float* x, const int* indx, float* y,
		   const float c, const float s);

  void cblas_droti(const int N, double* x, const int* indx, double* y,
		   const double c, const double s);
  
  void cblas_ssctr(const int N, float* x, const int* indx, float* y);
  
  void cblas_dsctr(const int N, double* x, const int* indx, double* y);
  
  void cblas_csctr(const int N, void* x, const int* indx, void* y);
  
  void cblas_zsctr(const int N, void* x, const int* indx, void* y);
  
  // Blas Level 2
  
  void mkl_cspblas_scsrgemv(char *transa, int *m, float *a, int *ia,
                            int *ja, float *x, float *y);

  void mkl_cspblas_dcsrgemv(char *transa, int *m, double *a, int *ia,
                            int *ja, double *x, double *y);

  void mkl_cspblas_ccsrgemv(char *transa, int *m, void *a, int *ia,
                            int *ja, void *x, void *y);

  void mkl_cspblas_zcsrgemv(char *transa, int *m, void *a, int *ia,
                            int *ja, void *x, void *y);

  void mkl_cspblas_scsrsymv(char *uplo, int *m, float *a, int *ia,
                            int *ja, float *x, float *y);

  void mkl_cspblas_dcsrsymv(char *uplo, int *m, double *a, int *ia,
                            int *ja, double *x, double *y);

  void mkl_cspblas_ccsrsymv(char *uplo, int *m, void *a, int *ia,
                            int *ja, void *x, void *y);

  void mkl_cspblas_zcsrsymv(char *uplo, int *m, void *a, int *ia,
                            int *ja, void *x, void *y);

  void mkl_cspblas_scsrtrsv(char *uplo, char *transa, char *diag, int *m,
                            float *a, int *ia, int *ja, float *x, float *y);

  void mkl_cspblas_dcsrtrsv(char *uplo, char *transa, char *diag, int *m,
                            double *a, int *ia, int *ja, double *x, double *y);

  void mkl_cspblas_ccsrtrsv(char *uplo, char *transa, char *diag, int *m,
                            void *a, int *ia, int *ja, void *x, void *y);

  void mkl_cspblas_zcsrtrsv(char *uplo, char *transa, char *diag, int *m,
                            void *a, int *ia, int *ja, void *x, void *y);
  
  void mkl_scsrmv(char *transa, int *m, int *n, const float *alpha,
                  char *matdescra, float *a, int *inda, int *ptrb, int *ptre,
                  float *x, const float *beta, float *y);
  
  void mkl_dcsrmv(char *transa, int *m, int *n, const double *alpha,
                  char *matdescra, double *a, int *inda, int *ptrb, int *ptre,
                  double *x, const double *beta, double *y);

  void mkl_ccsrmv(char *transa, int *m, int *n, const void *alpha,
                  char *matdescra, void *a, int *inda, int *ptrb, int *ptre,
                  void *x, const void *beta, void *y);
  
  void mkl_zcsrmv(char *transa, int *m, int *n, const void *alpha,
                  char *matdescra, void *a, int *inda, int *ptrb, int *ptre,
                  void *x, const void *beta, void *y);
  
  // Blas Level 3
  
  void mkl_scsrmm(char *transa, int *m, int *n, int *k, const float *alpha,
                  char *matdescra, float *a, int *inda, int *ptrb, int *ptre,
                  float *b, int *ldb, const float *beta, float *c, int *ldc);

  void mkl_dcsrmm(char *transa, int *m, int *n, int *k, const double *alpha,
                  char *matdescra, double *a, int *inda, int *ptrb, int *ptre,
                  double *b, int *ldb, const double *beta, double *c, int *ldc);

  void mkl_ccsrmm(char *transa, int *m, int *n, int *k, const void *alpha,
                  char *matdescra, void *a, int *inda, int *ptrb, int *ptre,
                  void *b, int *ldb, const void *beta, void *c, int *ldc);

  void mkl_zcsrmm(char *transa, int *m, int *n, int *k, const void *alpha,
                  char *matdescra, void *a, int *inda, int *ptrb, int *ptre,
                  void *b, int *ldb, const void *beta, void *c, int *ldc);

  void mkl_scsrsm(char *transa, int *m, int *n, const float *alpha,
                  char *matdescra, float *a, int *inda, int *ptrb, int *ptre,
                  float *b, int *ldb, float *c, int *ldc);

  void mkl_dcsrsm(char *transa, int *m, int *n, const double *alpha,
                  char *matdescra, double *a, int *inda, int *ptrb, int *ptre,
                  double *b, int *ldb, double *c, int *ldc);

  void mkl_ccsrsm(char *transa, int *m, int *n, const void *alpha,
                  char *matdescra, void *a, int *inda, int *ptrb, int *ptre,
                  void *b, int *ldb, void *c, int *ldc);

  void mkl_zcsrsm(char *transa, int *m, int *n, const void *alpha,
                  char *matdescra, void *a, int *inda, int *ptrb, int *ptre,
                  void *b, int *ldb, void *c, int *ldc);
  
  void mkl_scsradd(char *transa, int *job, int *sort, int *m, int *n,
                   float *a, int *ja, int *ia, const float *beta,
                   float *b, int *jb, int *ib, float *c,
                   int *jc, int *ic, int *nzmax, int *ierr);

  void mkl_dcsradd(char *transa, int *job, int *sort, int *m, int *n,
                   double *a, int *ja, int *ia, const double *beta,
                   double *b, int *jb, int *ib, double *c,
                   int *jc, int *ic, int *nzmax, int *ierr);

  void mkl_ccsradd(char *transa, int *job, int *sort, int *m, int *n,
                   void *a, int *ja, int *ia, const void *beta,
                   void *b, int *jb, int *ib, void *c,
                   int *jc, int *ic, int *nzmax, int *ierr);

  void mkl_zcsradd(char *transa, int *job, int *sort, int *m, int *n,
                   void *a, int *ja, int *ia, const void *beta,
                   void *b, int *jb, int *ib, void *c,
                   int *jc, int *ic, int *nzmax, int *ierr);
  
  void mkl_scsrmultcsr(char *transa, int *job, int *sort, int *m, int *n, int *k,
                       float *a, int *ja, int *ia, float *b, int *jb, int *ib,
                       float *c, int *jc, int *ic, int *nz, int *ierr);

  void mkl_dcsrmultcsr(char *transa, int *job, int *sort, int *m, int *n, int *k,
                       double *a, int *ja, int *ia, double *b, int *jb, int *ib,
                       double *c, int *jc, int *ic, int *nz, int *ierr);

  void mkl_ccsrmultcsr(char *transa, int *job, int *sort, int *m, int *n, int *k,
                       void *a, int *ja, int *ia, void *b, int *jb, int *ib,
                       void *c, int *jc, int *ic, int *nz, int *ierr);

  void mkl_zcsrmultcsr(char *transa, int *job, int *sort, int *m, int *n, int *k,
                       void *a, int *ja, int *ia, void *b, int *jb, int *ib,
                       void *c, int *jc, int *ic, int *nz, int *ierr);
  
}

namespace Seldon
{

  template<class Alloc1, class Alloc2>
  void Add(const float& alpha, const Vector<float, VectFull, Alloc1>& X,
           const float& beta, Vector<float, VectFull, Alloc2>& Y);
 
  template<class Alloc1, class Alloc2>
  void Add(const double& alpha, const Vector<double, VectFull, Alloc1>& X,
           const double& beta, Vector<double, VectFull, Alloc2>& Y);

  template<class Alloc1, class Alloc2>
  void Add(const complex<float>& alpha,
           const Vector<complex<float>, VectFull, Alloc1>& X,
           const complex<float>& beta,
           Vector<complex<float>, VectFull, Alloc2>& Y);
  
  template<class Alloc1, class Alloc2>
  void Add(const complex<double>& alpha,
           const Vector<complex<double>, VectFull, Alloc1>& X,
           const complex<double>& beta,
           Vector<complex<double>, VectFull, Alloc2>& Y);
  
  /***********************
   * Sparse Blas Level 1 *
   ***********************/
  
  
  template<class Alloc1, class Alloc2>
  void Add(const float& alpha, const Vector<float, VectSparse, Alloc1>& X,
	   Vector<float, VectFull, Alloc2>& Y);

  template<class Alloc1, class Alloc2>
  void Add(const double& alpha, const Vector<double, VectSparse, Alloc1>& X,
	   Vector<double, VectFull, Alloc2>& Y);
  
  template<class Alloc1, class Alloc2>
  void Add(const complex<float>& alpha,
           const Vector<complex<float>, VectSparse, Alloc1>& X,
           Vector<complex<float>, VectFull, Alloc2>& Y);

  template<class Alloc1, class Alloc2>
  void Add(const complex<double>& alpha,
           const Vector<complex<double>, VectSparse, Alloc1>& X,
           Vector<complex<double>, VectFull, Alloc2>& Y);

  template<class Alloc1, class Alloc2>
  float DotProdVector(const Vector<float, VectSparse, Alloc1>& x,
                      const Vector<float, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  double DotProdVector(const Vector<double, VectSparse, Alloc1>& x,
                       const Vector<double, VectFull, Alloc2>& y);
  
  template<class Alloc1, class Alloc2>
  complex<float> DotProdVector(const Vector<complex<float>, VectSparse, Alloc1>& x,
                               const Vector<complex<float>, VectFull, Alloc2>& y);
  
  template<class Alloc1, class Alloc2>
  complex<double> DotProdVector(const Vector<complex<double>, VectSparse, Alloc1>& x,
                                const Vector<complex<double>, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  float DotProdConjVector(const Vector<float, VectSparse, Alloc1>& x,
                          const Vector<float, VectFull, Alloc2>& y);
  
  template<class Alloc1, class Alloc2>
  double DotProdConjVector(const Vector<double, VectSparse, Alloc1>& x,
                           const Vector<double, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  complex<float> DotProdConjVector(const Vector<complex<float>, VectSparse, Alloc1>& x,
                                   const Vector<complex<float>, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  complex<double> DotProdConjVector(const Vector<complex<double>, VectSparse, Alloc1>& x,
                                    const Vector<complex<double>, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<float, VectFull, Alloc1>& y,
			 Vector<float, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<double, VectFull, Alloc1>& y,
			 Vector<double, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<complex<float>, VectFull, Alloc1>& y,
			 Vector<complex<float>, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<complex<double>, VectFull, Alloc1>& y,
			 Vector<complex<double>, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<float, VectFull, Alloc1>& y,
			     Vector<float, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<double, VectFull, Alloc1>& y,
			     Vector<double, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<complex<float>, VectFull, Alloc1>& y,
			     Vector<complex<float>, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<complex<double>, VectFull, Alloc1>& y,
			     Vector<complex<double>, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void ApplyRot(Vector<float, VectSparse, Alloc1>& x,
		Vector<float, VectFull, Alloc2>& y,
		const float& c, const float& s);

  template<class Alloc1, class Alloc2>
  void ApplyRot(Vector<double, VectSparse, Alloc1>& x,
		Vector<double, VectFull, Alloc2>& y,
		const double& c, const double& s);

  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<float, VectSparse, Alloc1>& x,
			  Vector<float, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<double, VectSparse, Alloc1>& x,
			  Vector<double, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<complex<float>, VectSparse, Alloc1>& x,
			  Vector<complex<float>, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<complex<double>, VectSparse, Alloc1>& x,
			  Vector<complex<double>, VectFull, Alloc2>& y);
  
  
  /***********************
   * Sparse Blas Level 2 *
   ***********************/

  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<float, General, RowSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<double, General, RowSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const float& alpha,
                 const Matrix<float, General, RowSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const double& alpha,
                 const Matrix<double, General, RowSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const complex<float>& alpha,
                 const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const complex<double>& alpha,
                 const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const float& alpha,
                 const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const double& alpha,
                 const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const complex<float>& alpha,
                 const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const complex<double>& alpha,
                 const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<float, General, RowSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<double, General, RowSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const float& alpha,
                    const Matrix<float, General, RowSparse, Alloc1>& A,
                    const Vector<float, VectFull, Alloc2>& X,
                    const float& beta,
                    Vector<float, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const double& alpha,
                    const Matrix<double, General, RowSparse, Alloc1>& A,
                    const Vector<double, VectFull, Alloc2>& X,
                    const double& beta,
                    Vector<double, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<float>& alpha,
                    const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                    const Vector<complex<float>, VectFull, Alloc2>& X,
                    const complex<float>& beta,
                    Vector<complex<float>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<double>& alpha,
                    const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                    const Vector<complex<double>, VectFull, Alloc2>& X,
                    const complex<double>& beta,
                    Vector<complex<double>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const float& alpha,
                    const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<float, VectFull, Alloc2>& X,
                    const float& beta,
                    Vector<float, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const double& alpha,
                    const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<double, VectFull, Alloc2>& X,
                    const double& beta,
                    Vector<double, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<float>& alpha,
                    const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<complex<float>, VectFull, Alloc2>& X,
                    const complex<float>& beta,
                    Vector<complex<float>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<double>& alpha,
                    const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<complex<double>, VectFull, Alloc2>& X,
                    const complex<double>& beta,
                    Vector<complex<double>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const float& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<float, General, RowSparse, Alloc1>& A,
                    const Vector<float, VectFull, Alloc2>& X,
                    const float& beta,
                    Vector<float, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const double& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<double, General, RowSparse, Alloc1>& A,
                    const Vector<double, VectFull, Alloc2>& X,
                    const double& beta,
                    Vector<double, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<float>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                    const Vector<complex<float>, VectFull, Alloc2>& X,
                    const complex<float>& beta,
                    Vector<complex<float>, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<double>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                    const Vector<complex<double>, VectFull, Alloc2>& X,
                    const complex<double>& beta,
                    Vector<complex<double>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const float& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<float, VectFull, Alloc2>& X,
                    const float& beta,
                    Vector<float, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const double& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<double, VectFull, Alloc2>& X,
                    const double& beta,
                    Vector<double, VectFull, Alloc3>& Y);

  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<float>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<complex<float>, VectFull, Alloc2>& X,
                    const complex<float>& beta,
                    Vector<complex<float>, VectFull, Alloc3>& Y);
  
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<double>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<complex<double>, VectFull, Alloc2>& X,
                    const complex<double>& beta,
                    Vector<complex<double>, VectFull, Alloc3>& Y);
  
  template <class Prop0, class Allocator0, class Allocator1>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<float, Prop0, RowSparse, Allocator0>& A,
             const Vector<float, VectFull, Allocator1>& X,
             Vector<float, VectFull, Allocator1>& Y);
  
  template <class Prop0, class Allocator0, class Allocator1>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<double, Prop0, RowSparse, Allocator0>& A,
             const Vector<double, VectFull, Allocator1>& X,
             Vector<double, VectFull, Allocator1>& Y);

  template <class Prop0, class Allocator0, class Allocator1>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<float>, Prop0, RowSparse, Allocator0>& A,
             const Vector<complex<float>, VectFull, Allocator1>& X,
             Vector<complex<float>, VectFull, Allocator1>& Y);

  template <class Prop0, class Allocator0, class Allocator1>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<double>, Prop0, RowSparse, Allocator0>& A,
             const Vector<complex<double>, VectFull, Allocator1>& X,
             Vector<complex<double>, VectFull, Allocator1>& Y);
  
  template <class Prop0, class Allocator0, class Allocator1>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<float, Prop0, RowSymSparse, Allocator0>& A,
             const Vector<float, VectFull, Allocator1>& X,
             Vector<float, VectFull, Allocator1>& Y);
  
  
  template <class Prop0, class Allocator0, class Allocator1>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<double, Prop0, RowSymSparse, Allocator0>& A,
             const Vector<double, VectFull, Allocator1>& X,
             Vector<double, VectFull, Allocator1>& Y);
  
  template <class Prop0, class Allocator0, class Allocator1>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<float>, Prop0, RowSymSparse, Allocator0>& A,
             const Vector<complex<float>, VectFull, Allocator1>& X,
             Vector<complex<float>, VectFull, Allocator1>& Y);
  
  template <class Prop0, class Allocator0, class Allocator1>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<double>, Prop0, RowSymSparse, Allocator0>& A,
             const Vector<complex<double>, VectFull, Allocator1>& X,
             Vector<complex<double>, VectFull, Allocator1>& Y);
  
  
  /***********************
   * Sparse Blas Level 3 *
   ***********************/

  
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const float& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<float, General, RowSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<float, General, RowMajor, Alloc1>& B,
              const float& beta,
                    Matrix<float, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const double& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<double, General, RowSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<double, General, RowMajor, Alloc1>& B,
                    const double& beta,
                    Matrix<double, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const complex<float>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<float>, General, RowSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<complex<float>, General, RowMajor, Alloc1>& B,
                    const complex<float>& beta,
                    Matrix<complex<float>, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const complex<double>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<double>, General, RowSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<complex<double>, General, RowMajor, Alloc1>& B,
                    const complex<double>& beta,
                    Matrix<complex<double>, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const float& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<float, Symmetric, RowSymSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<float, General, RowMajor, Alloc1>& B,
                    const float& beta,
                    Matrix<float, General, RowMajor, Alloc2>& C);
  
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const double& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<double, Symmetric, RowSymSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<double, General, RowMajor, Alloc1>& B,
                    const double& beta,
                    Matrix<double, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const complex<float>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<complex<float>, General, RowMajor, Alloc1>& B,
                    const complex<float>& beta,
                    Matrix<complex<float>, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const complex<double>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<complex<double>, General, RowMajor, Alloc1>& B,
                    const complex<double>& beta,
                    Matrix<complex<double>, General, RowMajor, Alloc2>& C);
  
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<float, General, RowSparse, Alloc0>& A,
             const Matrix<float, General, RowMajor, Alloc1>& B,
             Matrix<float, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<double, General, RowSparse, Alloc0>& A,
             const Matrix<double, General, RowMajor, Alloc1>& B,
             Matrix<double, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<float>, General, RowSparse, Alloc0>& A,
             const Matrix<complex<float>, General, RowMajor, Alloc1>& B,
             Matrix<complex<float>, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<double>, General, RowSparse, Alloc0>& A,
             const Matrix<complex<double>, General, RowMajor, Alloc1>& B,
             Matrix<complex<double>, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<float, Symmetric, RowSymSparse, Alloc0>& A,
             const Matrix<float, General, RowMajor, Alloc1>& B,
             Matrix<float, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<double, Symmetric, RowSymSparse, Alloc0>& A,
             const Matrix<double, General, RowMajor, Alloc1>& B,
             Matrix<double, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc0>& A,
             const Matrix<complex<float>, General, RowMajor, Alloc1>& B,
             Matrix<complex<float>, General, RowMajor, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc0>& A,
             const Matrix<complex<double>, General, RowMajor, Alloc1>& B,
             Matrix<complex<double>, General, RowMajor, Alloc2>& C);
  
  template<class Alloc0, class Alloc1>
  void AddMatrix(const float& alpha,
                 const Matrix<float, General, RowSparse, Alloc0>& A,
                 Matrix<float, General, RowSparse, Alloc1>& B);
  
  template<class Alloc0, class Alloc1>
  void AddMatrix(const double& alpha,
                 const Matrix<double, General, RowSparse, Alloc0>& A,
                 Matrix<double, General, RowSparse, Alloc1>& B);
  
  template<class Alloc0, class Alloc1>
  void AddMatrix(const complex<float>& alpha,
                 const Matrix<complex<float>, General, RowSparse, Alloc0>& A,
                 Matrix<complex<float>, General, RowSparse, Alloc1>& B);

  template<class Alloc0, class Alloc1>
  void AddMatrix(const complex<double>& alpha,
                 const Matrix<complex<double>, General, RowSparse, Alloc0>& A,
                 Matrix<complex<double>, General, RowSparse, Alloc1>& B);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltMatrix(const Matrix<float, General, RowSparse, Alloc0>& A,
                 const Matrix<float, General, RowSparse, Alloc1>& B,
                 Matrix<float, General, RowSparse, Alloc2>& C);
  
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltMatrix(const Matrix<double, General, RowSparse, Alloc0>& A,
                 const Matrix<double, General, RowSparse, Alloc1>& B,
                 Matrix<double, General, RowSparse, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltMatrix(const Matrix<complex<float>, General, RowSparse, Alloc0>& A,
                 const Matrix<complex<float>, General, RowSparse, Alloc1>& B,
                 Matrix<complex<float>, General, RowSparse, Alloc2>& C);

  template<class Alloc0, class Alloc1, class Alloc2>
  void MltMatrix(const Matrix<complex<double>, General, RowSparse, Alloc0>& A,
           const Matrix<complex<double>, General, RowSparse, Alloc1>& B,
           Matrix<complex<double>, General, RowSparse, Alloc2>& C);
  
} // end namespace

#define SELDON_FILE_MKL_SPARSE_HXX
#endif
