#ifndef SELDON_FILE_MKL_SPARSE_CXX

#include "Mkl_Sparse.hxx"

namespace Seldon
{

  //! computes Y = beta Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const float& alpha, const Vector<float, VectFull, Alloc1>& X,
           const float& beta, Vector<float, VectFull, Alloc2>& Y)
  {
    cblas_saxpby(X.GetM(), alpha, X.GetData(), 1, beta, Y.GetData(), 1);
  }


  //! computes Y = beta Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const double& alpha, const Vector<double, VectFull, Alloc1>& X,
           const double& beta, Vector<double, VectFull, Alloc2>& Y)
  {
    cblas_daxpby(X.GetM(), alpha, X.GetData(), 1, beta, Y.GetData(), 1);
  }


  //! computes Y = beta Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const complex<float>& alpha,
           const Vector<complex<float>, VectFull, Alloc1>& X,
           const complex<float>& beta,
           Vector<complex<float>, VectFull, Alloc2>& Y)
  {
    cblas_caxpby(X.GetM(), &alpha, X.GetDataVoid(), 1,
                 &beta, Y.GetDataVoid(), 1);
  }


  //! computes Y = beta Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const complex<double>& alpha,
           const Vector<complex<double>, VectFull, Alloc1>& X,
           const complex<double>& beta,
           Vector<complex<double>, VectFull, Alloc2>& Y)
  {
    cblas_zaxpby(X.GetM(), &alpha, X.GetDataVoid(), 1,
                 &beta, Y.GetDataVoid(), 1);
  }

  
  /***********************
   * Sparse Blas Level 1 *
   ***********************/
  
  
  //! computes Y = Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const float& alpha, const Vector<float, VectSparse, Alloc1>& X,
	   Vector<float, VectFull, Alloc2>& Y)
  {
    cblas_saxpyi(X.GetM(), alpha, X.GetData(), X.GetIndex(), Y.GetData());
  }

  
  //! computes Y = Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const double& alpha, const Vector<double, VectSparse, Alloc1>& X,
	   Vector<double, VectFull, Alloc2>& Y)
  {
    cblas_daxpyi(X.GetM(), alpha, X.GetData(), X.GetIndex(), Y.GetData());
  }


  //! computes Y = Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const complex<float>& alpha,
           const Vector<complex<float>, VectSparse, Alloc1>& X,
	   Vector<complex<float>, VectFull, Alloc2>& Y)
  {
    cblas_caxpyi(X.GetM(), &alpha, X.GetDataVoid(), X.GetIndex(), Y.GetDataVoid());
  }


  //! computes Y = Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const complex<double>& alpha,
           const Vector<complex<double>, VectSparse, Alloc1>& X,
	   Vector<complex<double>, VectFull, Alloc2>& Y)
  {
    cblas_zaxpyi(X.GetM(), &alpha, X.GetDataVoid(), X.GetIndex(), Y.GetDataVoid());
  }

  
  //! returns x^T y
  template<class Alloc1, class Alloc2>
  float DotProdVector(const Vector<float, VectSparse, Alloc1>& x,
                      const Vector<float, VectFull, Alloc2>& y)
  {
    return cblas_sdoti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  double DotProdVector(const Vector<double, VectSparse, Alloc1>& x,
                       const Vector<double, VectFull, Alloc2>& y)
  {
    return cblas_ddoti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  complex<float> DotProdVector(const Vector<complex<float>, VectSparse, Alloc1>& x,
                               const Vector<complex<float>, VectFull, Alloc2>& y)
  {
    complex<float> scal;
    cblas_cdotui_sub(x.GetM(), x.GetDataVoid(), x.GetIndex(),
		     y.GetDataVoid(), reinterpret_cast<void*>(&scal));
    
    return scal;
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  complex<double> DotProdVector(const Vector<complex<double>, VectSparse, Alloc1>& x,
                                const Vector<complex<double>, VectFull, Alloc2>& y)
  {
    complex<double> scal;
    cblas_zdotui_sub(x.GetM(), x.GetDataVoid(), x.GetIndex(),
		     y.GetDataVoid(), reinterpret_cast<void*>(&scal));
    
    return scal;
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  float DotProdConjVector(const Vector<float, VectSparse, Alloc1>& x,
                          const Vector<float, VectFull, Alloc2>& y)
  {
    return cblas_sdoti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  double DotProdConjVector(const Vector<double, VectSparse, Alloc1>& x,
                           const Vector<double, VectFull, Alloc2>& y)
  {
    return cblas_ddoti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! returns x^H y
  template<class Alloc1, class Alloc2>
  complex<float> DotProdConjVector(const Vector<complex<float>, VectSparse, Alloc1>& x,
                                   const Vector<complex<float>, VectFull, Alloc2>& y)
  {
    complex<float> scal;
    cblas_cdotci_sub(x.GetM(), x.GetDataVoid(), x.GetIndex(),
		     y.GetDataVoid(), reinterpret_cast<void*>(&scal));
    
    return scal;
  }


  //! returns x^H y
  template<class Alloc1, class Alloc2>
  complex<double> DotProdConjVector(const Vector<complex<double>, VectSparse, Alloc1>& x,
                                    const Vector<complex<double>, VectFull, Alloc2>& y)
  {
    complex<double> scal;
    cblas_zdotci_sub(x.GetM(), x.GetDataVoid(), x.GetIndex(),
		     y.GetDataVoid(), reinterpret_cast<void*>(&scal));
    
    return scal;
  }

  
  //! x = y(index)
  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<float, VectFull, Alloc1>& y,
			 Vector<float, VectSparse, Alloc2>& x)
  {
    cblas_sgthr(x.GetM(), y.GetData(), x.GetData(), x.GetIndex());
  }


  //! x = y(index)
  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<double, VectFull, Alloc1>& y,
			 Vector<double, VectSparse, Alloc2>& x)
  {
    cblas_dgthr(x.GetM(), y.GetData(), x.GetData(), x.GetIndex());
  }


  //! x = y(index)
  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<complex<float>, VectFull, Alloc1>& y,
			 Vector<complex<float>, VectSparse, Alloc2>& x)
  {
    cblas_cgthr(x.GetM(), y.GetDataVoid(), x.GetDataVoid(), x.GetIndex());
  }


  //! x = y(index)
  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<complex<double>, VectFull, Alloc1>& y,
			 Vector<complex<double>, VectSparse, Alloc2>& x)
  {
    cblas_zgthr(x.GetM(), y.GetDataVoid(), x.GetDataVoid(), x.GetIndex());
  }


  //! x = y(index) and y(index) = 0
  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<float, VectFull, Alloc1>& y,
			 Vector<float, VectSparse, Alloc2>& x)
  {
    cblas_sgthrz(x.GetM(), y.GetData(), x.GetData(), x.GetIndex());
  }


  //! x = y(index) and y(index) = 0
  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<double, VectFull, Alloc1>& y,
			     Vector<double, VectSparse, Alloc2>& x)
  {
    cblas_dgthrz(x.GetM(), y.GetData(), x.GetData(), x.GetIndex());
  }


  //! x = y(index) and y(index) = 0
  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<complex<float>, VectFull, Alloc1>& y,
			     Vector<complex<float>, VectSparse, Alloc2>& x)
  {
    cblas_cgthrz(x.GetM(), y.GetDataVoid(), x.GetDataVoid(), x.GetIndex());
  }


  //! x = y(index) and y(index) = 0
  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<complex<double>, VectFull, Alloc1>& y,
			     Vector<complex<double>, VectSparse, Alloc2>& x)
  {
    cblas_zgthrz(x.GetM(), y.GetDataVoid(), x.GetDataVoid(), x.GetIndex());
  }

  
  //! applies rotation to vectors x and y
  template<class Alloc1, class Alloc2>
  void ApplyRot(Vector<float, VectSparse, Alloc1>& x,
		Vector<float, VectFull, Alloc2>& y,
		const float& c, const float& s)
  {
    cblas_sroti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData(), c, s);
  }


  //! applies rotation to vectors x and y
  template<class Alloc1, class Alloc2>
  void ApplyRot(Vector<double, VectSparse, Alloc1>& x,
		Vector<double, VectFull, Alloc2>& y,
		const double& c, const double& s)
  {
    cblas_droti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData(), c, s);
  }
  
  
  //! y(index) = x
  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<float, VectSparse, Alloc1>& x,
			  Vector<float, VectFull, Alloc2>& y)
  {
    cblas_ssctr(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! y(index) = x
  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<double, VectSparse, Alloc1>& x,
			  Vector<double, VectFull, Alloc2>& y)
  {
    cblas_dsctr(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! y(index) = x
  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<complex<float>, VectSparse, Alloc1>& x,
			  Vector<complex<float>, VectFull, Alloc2>& y)
  {
    cblas_csctr(x.GetM(), x.GetDataVoid(), x.GetIndex(), y.GetDataVoid());
  }


  //! y(index) = x
  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<complex<double>, VectSparse, Alloc1>& x,
			  Vector<complex<double>, VectFull, Alloc2>& y)
  {
    cblas_zsctr(x.GetM(), x.GetDataVoid(), x.GetIndex(), y.GetDataVoid());
  }


  /***********************
   * Sparse Blas Level 2 *
   ***********************/
  
  
  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<float, General, RowSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA('N');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_scsrgemv(&transA, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                         X.GetData(), Y.GetData());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<double, General, RowSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA('N');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_dcsrgemv(&transA, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                         X.GetData(), Y.GetData());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA('N');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_ccsrgemv(&transA, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                         X.GetDataVoid(), Y.GetDataVoid());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA('N');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_zcsrgemv(&transA, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                         X.GetDataVoid(), Y.GetDataVoid());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_scsrsymv(&uplo, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                         X.GetData(), Y.GetData());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif
    
    char uplo('U');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_dcsrsymv(&uplo, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                         X.GetData(), Y.GetData());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_ccsrsymv(&uplo, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                         X.GetDataVoid(), Y.GetDataVoid());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_zcsrsymv(&uplo, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                         X.GetDataVoid(), Y.GetDataVoid());
  }
  

  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const float& alpha,
                 const Matrix<float, General, RowSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA('N');
    int ma = A.GetM();
    float zero(0), one(1);
    if (alpha == zero)
      Y.Zero();
    else
      {
        Vector<int> ptr; A.FillPtrInt(ptr);
        mkl_cspblas_scsrgemv(&transA, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                             X.GetData(), Y.GetData());
        
        if (alpha != one)
          Mlt(alpha, Y);
      }
  }


  //! Y = alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const double& alpha,
                 const Matrix<double, General, RowSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA('N');
    int ma = A.GetM();
    double zero(0), one(1);
    if (alpha == zero)
      Y.Zero();
    else
      {
        Vector<int> ptr; A.FillPtrInt(ptr);
        mkl_cspblas_dcsrgemv(&transA, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                             X.GetData(), Y.GetData());
        
        if (alpha != one)
          Mlt(alpha, Y);
      }        
  }


  //! Y = alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const complex<float>& alpha,
                 const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA('N');
    int ma = A.GetM();
    complex<float> zero(0, 0), one(1, 0);
    if (alpha == zero)
      Y.Zero();
    else
      {
        Vector<int> ptr; A.FillPtrInt(ptr);
        mkl_cspblas_ccsrgemv(&transA, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                             X.GetDataVoid(), Y.GetDataVoid());
        
        if (alpha != one)
          Mlt(alpha, Y);
      }        
  }


  //! Y = alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const complex<double>& alpha,
                 const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA('N');
    int ma = A.GetM();
    complex<double> zero(0, 0), one(1, 0);
    if (alpha == zero)
      Y.Zero();
    else
      {
        Vector<int> ptr; A.FillPtrInt(ptr);
        mkl_cspblas_zcsrgemv(&transA, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                             X.GetDataVoid(), Y.GetDataVoid());
        
        if (alpha != one)
          Mlt(alpha, Y);
      }        
  }


  //! Y = alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const float& alpha,
                 const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    float zero(0), one(1);
    if (alpha == zero)
      Y.Zero();
    else
      {
        Vector<int> ptr; A.FillPtrInt(ptr);
        mkl_cspblas_scsrsymv(&uplo, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                             X.GetData(), Y.GetData());
        
        if (alpha != one)
          Mlt(alpha, Y);
      }        
  }


  //! Y = alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const double& alpha,
                 const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    double zero(0), one(1);
    if (alpha == zero)
      Y.Zero();
    else
      {
        Vector<int> ptr; A.FillPtrInt(ptr);
        mkl_cspblas_dcsrsymv(&uplo, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                             X.GetData(), Y.GetData());
        
        if (alpha != one)
          Mlt(alpha, Y);
      }        
  }


  //! Y = alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const complex<float>& alpha,
                 const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    complex<float> zero(0, 0), one(1, 0);
    if (alpha == zero)
      Y.Zero();
    else
      {
        Vector<int> ptr; A.FillPtrInt(ptr);
        mkl_cspblas_ccsrsymv(&uplo, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                             X.GetDataVoid(), Y.GetDataVoid());

        if (alpha != one)
          Mlt(alpha, Y);
      }
  }


  //! Y = alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const complex<double>& alpha,
                 const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    complex<double> zero(0, 0), one(1, 0);
    if (alpha == zero)
      Y.Zero();
    else
      {
        Vector<int> ptr; A.FillPtrInt(ptr);
        mkl_cspblas_zcsrsymv(&uplo, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                             X.GetDataVoid(), Y.GetDataVoid());
        
        if (alpha != one)
          Mlt(alpha, Y);
      }
  }


  //! Y = A X or A^T X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<float, General, RowSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(trans, A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA(trans.Char());
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    if ( (ma == A.GetN()) || (trans.NoTrans()))
      mkl_cspblas_scsrgemv(&transA, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
			   X.GetData(), Y.GetData());
    else
      {
        Y.Zero();
        int na = A.GetN(); float one(1), zero(0);
        char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
        mkl_scsrmv(&transA, &ma, &na, &one, matdescra, A.GetData(), 
                   A.GetInd(), ptr.GetData(), ptr.GetData()+1,
                   X.GetData(), &zero, Y.GetData());
      }
  }


  //! Y = A X or A^T X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<double, General, RowSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(trans, A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA(trans.Char());
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    if ( (ma == A.GetN()) || (trans.NoTrans()))
      mkl_cspblas_dcsrgemv(&transA, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
			   X.GetData(), Y.GetData());
    else
      {
        Y.Zero();
        int na = A.GetN(); double one(1), zero(0);
        char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
        mkl_dcsrmv(&transA, &ma, &na, &one, matdescra, A.GetData(), 
                   A.GetInd(), ptr.GetData(), ptr.GetData()+1,
                   X.GetData(), &zero, Y.GetData());        
      }
  }


  //! Y = A X or A^T X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(trans, A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA(trans.Char());
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    if ( (ma == A.GetN()) || (trans.NoTrans()))
      mkl_cspblas_ccsrgemv(&transA, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
			   X.GetDataVoid(), Y.GetDataVoid());
    else
      {
        Y.Zero();
        int na = A.GetN(); complex<float> one(1, 0), zero(0, 0);
        char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
        mkl_ccsrmv(&transA, &ma, &na, &one, matdescra, A.GetDataVoid(), 
                   A.GetInd(), ptr.GetData(), ptr.GetData()+1,
                   X.GetDataVoid(), &zero, Y.GetDataVoid());        
      }
  }


  //! Y = A X or A^T X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(trans, A, X, Y, "MltVector(A, X, Y)");
#endif

    char transA(trans.Char());
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    if ( (ma == A.GetN()) || (trans.NoTrans()))
      mkl_cspblas_zcsrgemv(&transA, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
			   X.GetDataVoid(), Y.GetDataVoid());
    else
      {
        Y.Zero();
        int na = A.GetN(); complex<double> one(1, 0), zero(0, 0);
        char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
        mkl_zcsrmv(&transA, &ma, &na, &one, matdescra, A.GetDataVoid(), 
                   A.GetInd(), ptr.GetData(), ptr.GetData()+1,
                   X.GetDataVoid(), &zero, Y.GetDataVoid());        
      }
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<float, VectFull, Alloc2>& X,
                 Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_scsrsymv(&uplo, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                         X.GetData(), Y.GetData());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<double, VectFull, Alloc2>& X,
                 Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    char uplo('U');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_dcsrsymv(&uplo, &ma, A.GetData(), ptr.GetData(), A.GetInd(),
                         X.GetData(), Y.GetData());
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<float>, VectFull, Alloc2>& X,
                 Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    Vector<complex<float>, VectFull, Alloc2>& Xv
      = const_cast<Vector<complex<float>, VectFull, Alloc2>& >(X);
    
    if (trans.ConjTrans())
      Conjugate(Xv);

    char uplo('U');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_ccsrsymv(&uplo, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                         X.GetDataVoid(), Y.GetDataVoid());

    if (trans.ConjTrans())
      {
	Conjugate(Xv);
	Conjugate(Y);
      }
  }


  //! Y = A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltVector(const SeldonTranspose& trans,
                 const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                 const Vector<complex<double>, VectFull, Alloc2>& X,
                 Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltVector(A, X, Y)");
#endif

    Vector<complex<double>, VectFull, Alloc2>& Xv
      = const_cast<Vector<complex<double>, VectFull, Alloc2>& >(X);

    if (trans.ConjTrans())
      Conjugate(Xv);

    char uplo('U');
    int ma = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_zcsrsymv(&uplo, &ma, A.GetDataVoid(), ptr.GetData(), A.GetInd(),
                         X.GetDataVoid(), Y.GetDataVoid());

    if (trans.ConjTrans())
      {
	Conjugate(Xv);
	Conjugate(Y);
      }
  }
  

  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const float& alpha,
                    const Matrix<float, General, RowSparse, Alloc1>& A,
                    const Vector<float, VectFull, Alloc2>& X,
                    const float& beta,
                    Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAddVector(alpha, A, X, beta, Y)");
#endif

    char transa('N');
    int m = A.GetM(), n = A.GetN();
    Vector<int> ptr; A.FillPtrInt(ptr);
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    mkl_scsrmv(&transa, &m, &n, &alpha, matdescra, A.GetData(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetData(), &beta, Y.GetData());
  }
  
  
  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const double& alpha,
                    const Matrix<double, General, RowSparse, Alloc1>& A,
                    const Vector<double, VectFull, Alloc2>& X,
                    const double& beta,
                    Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAddVector(alpha, A, X, beta, Y)");
#endif

    char transa('N');
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_dcsrmv(&transa, &m, &n, &alpha, matdescra, A.GetData(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetData(), &beta, Y.GetData());
  }


  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<float>& alpha,
                    const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                    const Vector<complex<float>, VectFull, Alloc2>& X,
                    const complex<float>& beta,
                    Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAddVector(alpha, A, X, beta, Y)");
#endif

    char transa('N');
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_ccsrmv(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetDataVoid(), &beta, Y.GetDataVoid());
  }
  
  
  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<double>& alpha,
                    const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                    const Vector<complex<double>, VectFull, Alloc2>& X,
                    const complex<double>& beta,
                    Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAddVector(alpha, A, X, beta, Y)");
#endif

    char transa('N');
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_zcsrmv(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetDataVoid(), &beta, Y.GetDataVoid());
  }


  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const float& alpha,
                    const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<float, VectFull, Alloc2>& X,
                    const float& beta,
                    Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAddVector(alpha, A, X, beta, Y)");
#endif

    char transa('N');
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_scsrmv(&transa, &m, &n, &alpha, matdescra, A.GetData(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetData(), &beta, Y.GetData());
  }
  
  
  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const double& alpha,
                    const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<double, VectFull, Alloc2>& X,
                    const double& beta,
                    Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAddVector(alpha, A, X, beta, Y)");
#endif

    char transa('N');
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_dcsrmv(&transa, &m, &n, &alpha, matdescra, A.GetData(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetData(), &beta, Y.GetData());

  }


  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<float>& alpha,
                    const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<complex<float>, VectFull, Alloc2>& X,
                    const complex<float>& beta,
                    Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAddVector(alpha, A, X, beta, Y)");
#endif

    char transa('N');
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_ccsrmv(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetDataVoid(), &beta, Y.GetDataVoid());
  }
  
  
  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<double>& alpha,
                    const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<complex<double>, VectFull, Alloc2>& X,
                    const complex<double>& beta,
                    Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAddVector(alpha, A, X, beta, Y)");
#endif

    char transa('N');
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_zcsrmv(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetDataVoid(), &beta, Y.GetDataVoid());
  }


  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const float& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<float, General, RowSparse, Alloc1>& A,
                    const Vector<float, VectFull, Alloc2>& X,
                    const float& beta,
                    Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAddVector(alpha, trans, A, X, beta, Y)");
#endif

    char transa(TransA.Char());
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_scsrmv(&transa, &m, &n, &alpha, matdescra, A.GetData(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetData(), &beta, Y.GetData());
  }
  
  
  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const double& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<double, General, RowSparse, Alloc1>& A,
                    const Vector<double, VectFull, Alloc2>& X,
                    const double& beta,
                    Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAddVector(alpha, trans, A, X, beta, Y)");
#endif

    char transa(TransA.Char());
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_dcsrmv(&transa, &m, &n, &alpha, matdescra, A.GetData(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetData(), &beta, Y.GetData());
  }


  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<float>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<float>, General, RowSparse, Alloc1>& A,
                    const Vector<complex<float>, VectFull, Alloc2>& X,
                    const complex<float>& beta,
                    Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAddVector(alpha, trans, A, X, beta, Y)");
#endif

    char transa(TransA.Char());
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_ccsrmv(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetDataVoid(), &beta, Y.GetDataVoid());
  }
  
  
  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<double>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<double>, General, RowSparse, Alloc1>& A,
                    const Vector<complex<double>, VectFull, Alloc2>& X,
                    const complex<double>& beta,
                    Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAddVector(alpha, trans, A, X, beta, Y)");
#endif

    char transa(TransA.Char());
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_zcsrmv(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetDataVoid(), &beta, Y.GetDataVoid());
  }


  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const float& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<float, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<float, VectFull, Alloc2>& X,
                    const float& beta,
                    Vector<float, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAddVector(alpha, trans, A, X, beta, Y)");
#endif

    char transa(TransA.Char());
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_scsrmv(&transa, &m, &n, &alpha, matdescra, A.GetData(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetData(), &beta, Y.GetData());
  }
  
  
  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const double& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<double, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<double, VectFull, Alloc2>& X,
                    const double& beta,
                    Vector<double, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAddVector(alpha, trans, A, X, beta, Y)");
#endif

    char transa(TransA.Char());
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_dcsrmv(&transa, &m, &n, &alpha, matdescra, A.GetData(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetData(), &beta, Y.GetData());

  }


  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<float>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<complex<float>, VectFull, Alloc2>& X,
                    const complex<float>& beta,
                    Vector<complex<float>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAddVector(alpha, trans, A, X, beta, Y)");
#endif

    char transa(TransA.Char());
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_ccsrmv(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetDataVoid(), &beta, Y.GetDataVoid());
  }
  
  
  //! Y = beta Y + alpha A X
  template<class Alloc1, class Alloc2, class Alloc3>
  void MltAddVector(const complex<double>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc1>& A,
                    const Vector<complex<double>, VectFull, Alloc2>& X,
                    const complex<double>& beta,
                    Vector<complex<double>, VectFull, Alloc3>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAddVector(alpha, trans, A, X, beta, Y)");
#endif

    char transa(TransA.Char());
    int m = A.GetM(), n = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_zcsrmv(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(), 
               A.GetInd(), ptr.GetData(), ptr.GetData()+1,
               X.GetDataVoid(), &beta, Y.GetDataVoid());
  }

  
  //! solves by triangular upper part or lower part of A
  template <class Prop0, class Allocator0, class Allocator1>
  void
  Solve(const SeldonUplo& Uplo,
        const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<float, Prop0, RowSparse, Allocator0>& A,
	const Vector<float, VectFull, Allocator1>& X,
        Vector<float, VectFull, Allocator1>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "Solve(uplo, trans, diag, A, X)");
#endif
    
    char uplo(Uplo.Char()), transa(TransA.Char()), diaga(DiagA.Char());
    int m = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_scsrtrsv(&uplo, &transa, &diaga, &m, A.GetData(),
                         ptr.GetData(), A.GetInd(), X.GetData(), Y.GetData());
  }

  
  //! solves by triangular upper part or lower part of A
  template <class Prop0, class Allocator0, class Allocator1>
  void
  Solve(const SeldonUplo& Uplo,
        const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<double, Prop0, RowSparse, Allocator0>& A,
	const Vector<double, VectFull, Allocator1>& X,
        Vector<double, VectFull, Allocator1>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "Solve(uplo, trans, diag, A, X)");
#endif
    
    char uplo(Uplo.Char()), transa(TransA.Char()), diaga(DiagA.Char());
    int m = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_dcsrtrsv(&uplo, &transa, &diaga, &m, A.GetData(),
                         ptr.GetData(), A.GetInd(), X.GetData(), Y.GetData());
  }


  //! solves by triangular upper part or lower part of A
  template <class Prop0, class Allocator0, class Allocator1>
  void
  Solve(const SeldonUplo& Uplo,
        const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowSparse, Allocator0>& A,
	const Vector<complex<float>, VectFull, Allocator1>& X,
        Vector<complex<float>, VectFull, Allocator1>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "Solve(uplo, trans, diag, A, X)");
#endif
    
    char uplo(Uplo.Char()), transa(TransA.Char()), diaga(DiagA.Char());
    int m = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_ccsrtrsv(&uplo, &transa, &diaga, &m, A.GetDataVoid(),
                         ptr.GetData(), A.GetInd(), X.GetDataVoid(), Y.GetDataVoid());
  }

  
  //! solves by triangular upper part or lower part of A
  template <class Prop0, class Allocator0, class Allocator1>
  void
  Solve(const SeldonUplo& Uplo,
        const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0, RowSparse, Allocator0>& A,
	const Vector<complex<double>, VectFull, Allocator1>& X,
        Vector<complex<double>, VectFull, Allocator1>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "Solve(uplo, trans, diag, A, X)");
#endif
    
    char uplo(Uplo.Char()), transa(TransA.Char()), diaga(DiagA.Char());
    int m = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_zcsrtrsv(&uplo, &transa, &diaga, &m, A.GetDataVoid(),
                         ptr.GetData(), A.GetInd(), X.GetDataVoid(), Y.GetDataVoid());
  }


  //! solves by triangular upper part of A
  template <class Prop0, class Allocator0, class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<float, Prop0, RowSymSparse, Allocator0>& A,
	const Vector<float, VectFull, Allocator1>& X,
        Vector<float, VectFull, Allocator1>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "Solve(trans, diag, A, X)");
#endif
    
    char uplo('U'), transa(TransA.Char()), diaga(DiagA.Char());
    int m = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_scsrtrsv(&uplo, &transa, &diaga, &m, A.GetData(),
                         ptr.GetData(), A.GetInd(), X.GetData(), Y.GetData());
  }


  //! solves by triangular upper part of A
  template <class Prop0, class Allocator0, class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<double, Prop0, RowSymSparse, Allocator0>& A,
	const Vector<double, VectFull, Allocator1>& X,
        Vector<double, VectFull, Allocator1>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "Solve(trans, diag, A, X)");
#endif
    
    char uplo('U'), transa(TransA.Char()), diaga(DiagA.Char());
    int m = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_dcsrtrsv(&uplo, &transa, &diaga, &m, A.GetData(),
                         ptr.GetData(), A.GetInd(), X.GetData(), Y.GetData());
  }


  //! solves by triangular upper part of A
  template <class Prop0, class Allocator0, class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowSymSparse, Allocator0>& A,
	const Vector<complex<float>, VectFull, Allocator1>& X,
        Vector<complex<float>, VectFull, Allocator1>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "Solve(trans, diag, A, X)");
#endif
    
    char uplo('U'), transa(TransA.Char()), diaga(DiagA.Char());
    int m = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_ccsrtrsv(&uplo, &transa, &diaga, &m, A.GetDataVoid(),
                         ptr.GetData(), A.GetInd(), X.GetDataVoid(), Y.GetDataVoid());
  }

  
  //! solves by triangular upper part of A
  template <class Prop0, class Allocator0, class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0, RowSymSparse, Allocator0>& A,
	const Vector<complex<double>, VectFull, Allocator1>& X,
        Vector<complex<double>, VectFull, Allocator1>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "Solve(trans, diag, A, X)");
#endif
    
    char uplo('U'), transa(TransA.Char()), diaga(DiagA.Char());
    int m = A.GetM();
    Vector<int> ptr; A.FillPtrInt(ptr);
    mkl_cspblas_zcsrtrsv(&uplo, &transa, &diaga, &m, A.GetDataVoid(),
                         ptr.GetData(), A.GetInd(), X.GetDataVoid(), Y.GetDataVoid());
  }


  /***********************
   * Sparse Blas Level 3 *
   ***********************/

  
  //! C = alpha A B + beta C
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const float& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<float, General, RowSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<float, General, RowMajor, Alloc1>& B,
                    const float& beta,
                    Matrix<float, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C, "MltAddVector(alpha, transA, A, transB, B, beta, C)");
#endif
    
    if (!B.NoTrans())
      {
	cout << "Only product with non-tranpose B is available" << endl;
	abort();
      }

    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN(), k = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_scsrmm(&transa, &m, &n, &k, &alpha, matdescra, A.GetData(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetData(),
               &n, &beta, C.GetData(), &n);
  }


  //! C = alpha A B + beta C
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const double& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<double, General, RowSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<double, General, RowMajor, Alloc1>& B,
                    const double& beta,
                    Matrix<double, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAddVector(alpha, transA, A, transB, B, beta, C)");
#endif

    if (!TransB.NoTrans())
      {
	cout << "Only product with non-tranpose B is available" << endl;
	abort();
      }
    
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN(), k = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_dcsrmm(&transa, &m, &n, &k, &alpha, matdescra, A.GetData(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetData(),
               &n, &beta, C.GetData(), &n);
  }


  //! C = alpha A B + beta C
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const complex<float>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<float>, General, RowSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<complex<float>, General, RowMajor, Alloc1>& B,
                    const complex<float>& beta,
                    Matrix<complex<float>, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAddVector(alpha, transA, A, transB, B, beta, C)");
#endif

    if (!TransB.NoTrans())
      {
	cout << "Only product with non-tranpose B is available" << endl;
	abort();
      }

    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN(), k = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_ccsrmm(&transa, &m, &n, &k, &alpha, matdescra, A.GetDataVoid(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetDataVoid(),
               &n, &beta, C.GetDataVoid(), &n);
  }


  //! C = alpha A B + beta C
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const complex<double>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<double>, General, RowSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<complex<double>, General, RowMajor, Alloc1>& B,
                    const complex<double>& beta,
                    Matrix<complex<double>, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAddVector(alpha, transA, A, transB, B, beta, C)");
#endif

    if (!TransB.NoTrans())
      {
	cout << "Only product with non-tranpose B is available" << endl;
	abort();
      }

    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN(), k = A.GetN();
    char matdescra[6] = {'G', '0', '0', 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_zcsrmm(&transa, &m, &n, &k, &alpha, matdescra, A.GetDataVoid(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetDataVoid(),
               &n, &beta, C.GetDataVoid(), &n);
  }


  //! C = alpha A B + beta C
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const float& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<float, Symmetric, RowSymSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<float, General, RowMajor, Alloc1>& B,
                    const float& beta,
                    Matrix<float, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAddVector(alpha, transA, A, transB, B, beta, C)");
#endif

    if (!TransB.NoTrans())
      {
	cout << "Only product with non-tranpose B is available" << endl;
	abort();
      }

    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN(), k = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_scsrmm(&transa, &m, &n, &k, &alpha, matdescra, A.GetData(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetData(),
               &n, &beta, C.GetData(), &n);
  }


  //! C = alpha A B + beta C
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const double& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<double, Symmetric, RowSymSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<double, General, RowMajor, Alloc1>& B,
                    const double& beta,
                    Matrix<double, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAddVector(alpha, transA, A, transB, B, beta, C)");
#endif
    
    if (!TransB.NoTrans())
      {
	cout << "Only product with non-tranpose B is available" << endl;
	abort();
      }

    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN(), k = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_dcsrmm(&transa, &m, &n, &k, &alpha, matdescra, A.GetData(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetData(),
               &n, &beta, C.GetData(), &n);
  }


  //! C = alpha A B + beta C
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const complex<float>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<complex<float>, General, RowMajor, Alloc1>& B,
                    const complex<float>& beta,
                    Matrix<complex<float>, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAddVector(alpha, transA, A, transB, B, beta, C)");
#endif

    if (!TransB.NoTrans())
      {
	cout << "Only product with non-tranpose B is available" << endl;
	abort();
      }

    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN(), k = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_ccsrmm(&transa, &m, &n, &k, &alpha, matdescra, A.GetDataVoid(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetDataVoid(),
               &n, &beta, C.GetDataVoid(), &n);
  }


  //! C = alpha A B + beta C
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltAddMatrix(const complex<double>& alpha,
                    const SeldonTranspose& TransA,
                    const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc0>& A,
                    const SeldonTranspose& TransB,
                    const Matrix<complex<double>, General, RowMajor, Alloc1>& B,
                    const complex<double>& beta,
                    Matrix<complex<double>, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAddVector(alpha, transA, A, transB, B, beta, C)");
#endif

    if (!TransB.NoTrans())
      {
	cout << "Only product with non-tranpose B is available" << endl;
	abort();
      }

    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN(), k = A.GetN();
    char matdescra[6] = {'S', 'U', 'N', 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_zcsrmm(&transa, &m, &n, &k, &alpha, matdescra, A.GetDataVoid(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetDataVoid(),
               &n, &beta, C.GetDataVoid(), &n);
  }

  
  //! C = inv(T) B where T is a triangular part of A
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<float, General, RowSparse, Alloc0>& A,
             const Matrix<float, General, RowMajor, Alloc1>& B,
             Matrix<float, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "Solve(A, B, C)");
#endif

    float alpha(1);
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN();
    char matdescra[6] = {'T', Uplo.Char(), DiagA.Char(), 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_scsrsm(&transa, &m, &n, &alpha, matdescra, A.GetData(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetData(),
               &n, C.GetData(), &n);
  }


  //! C = inv(T) B where T is a triangular part of A
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<double, General, RowSparse, Alloc0>& A,
             const Matrix<double, General, RowMajor, Alloc1>& B,
             Matrix<double, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "Solve(A, B, C)");
#endif

    double alpha(1);
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN();
    char matdescra[6] = {'T', Uplo.Char(), DiagA.Char(), 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_dcsrsm(&transa, &m, &n, &alpha, matdescra, A.GetData(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetData(),
               &n, C.GetData(), &n);
  }


  //! C = inv(T) B where T is a triangular part of A
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<float>, General, RowSparse, Alloc0>& A,
             const Matrix<complex<float>, General, RowMajor, Alloc1>& B,
             Matrix<complex<float>, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "Solve(A, B, C)");
#endif

    complex<float> alpha(1, 0);
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN();
    char matdescra[6] = {'T', Uplo.Char(), DiagA.Char(), 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_ccsrsm(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetDataVoid(),
               &n, C.GetDataVoid(), &n);
  }


  //! C = inv(T) B where T is a triangular part of A
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonUplo& Uplo,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<double>, General, RowSparse, Alloc0>& A,
             const Matrix<complex<double>, General, RowMajor, Alloc1>& B,
             Matrix<complex<double>, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "Solve(A, B, C)");
#endif

    complex<double> alpha(1, 0);
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN();
    char matdescra[6] = {'T', Uplo.Char(), DiagA.Char(), 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_zcsrsm(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetDataVoid(),
               &n, C.GetDataVoid(), &n);
  }


  //! C = inv(T) B where T is a triangular part of A
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<float, Symmetric, RowSymSparse, Alloc0>& A,
             const Matrix<float, General, RowMajor, Alloc1>& B,
             Matrix<float, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "Solve(A, B, C)");
#endif

    float alpha(1);
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN();
    char matdescra[6] = {'T', 'U', DiagA.Char(), 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_scsrsm(&transa, &m, &n, &alpha, matdescra, A.GetData(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetData(),
               &n, C.GetData(), &n);
  }


  //! C = inv(T) B where T is a triangular part of A
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<double, Symmetric, RowSymSparse, Alloc0>& A,
             const Matrix<double, General, RowMajor, Alloc1>& B,
             Matrix<double, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "Solve(A, B, C)");
#endif

    double alpha(1);
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN();
    char matdescra[6] = {'T', 'U', DiagA.Char(), 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_dcsrsm(&transa, &m, &n, &alpha, matdescra, A.GetData(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetData(),
               &n, C.GetData(), &n);
  }


  //! C = inv(T) B where T is a triangular part of A
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<float>, Symmetric, RowSymSparse, Alloc0>& A,
             const Matrix<complex<float>, General, RowMajor, Alloc1>& B,
             Matrix<complex<float>, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "Solve(A, B, C)");
#endif

    complex<float> alpha(1, 0);
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN();
    char matdescra[6] = {'T', 'U', DiagA.Char(), 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_ccsrsm(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetDataVoid(),
               &n, C.GetDataVoid(), &n);
  }


  //! C = inv(T) B where T is a triangular part of A
  template<class Alloc0, class Alloc1, class Alloc2>
  void Solve(const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<complex<double>, Symmetric, RowSymSparse, Alloc0>& A,
             const Matrix<complex<double>, General, RowMajor, Alloc1>& B,
             Matrix<complex<double>, General, RowMajor, Alloc2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "Solve(A, B, C)");
#endif

    complex<double> alpha(1, 0);
    char transa(TransA.Char());
    int m = A.GetM(), n = C.GetN();
    char matdescra[6] = {'T', 'U', DiagA.Char(), 'C', '0', '0'};
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    mkl_zcsrsm(&transa, &m, &n, &alpha, matdescra, A.GetDataVoid(),
               A.GetInd(), ptrA.GetData(), ptrA.GetData()+1, B.GetDataVoid(),
               &n, C.GetDataVoid(), &n);
  }


  //! B = B + alpha A
  template<class Alloc0, class Alloc1>
  void AddMatrix(const float& alpha,
                 const Matrix<float, General, RowSparse, Alloc0>& A,
                 Matrix<float, General, RowSparse, Alloc1>& B)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (A.GetM() != B.GetM() || A.GetN() != B.GetN())
      throw WrongDim("Add(alpha, const Matrix<RowSparse>& A, "
                     "Matrix<RowSparse>& B)",
                     "Unable to add a " + to_str(A.GetM()) + " x "
                     + to_str(A.GetN()) + " matrix with a "
                     + to_str(B.GetM()) + " x " + to_str(B.GetN())
                     + " matrix.");
#endif

    // switching to 1-based indexing 
    int m = B.GetM(), n = B.GetN();
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    Vector<int> ptrB; B.FillPtrInt(ptrB);
    int* indB = B.GetInd(); int* indA = A.GetInd();
    for (int i = 0; i <= m; i++)
      {
        ptrA(i)++;
        ptrB(i)++;
      }
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]++;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]++;

    char transa('N');
    Vector<int> PtrC(m+1);
    
    // mkl_?csradd is called a first time to know the number of non-zero entries
    int request = 1, sort = 0, nz = 0, ierr = 0; 
    mkl_scsradd(&transa, &request, &sort, &m, &n, B.GetData(),
                B.GetInd(), ptrB.GetData(), &alpha, A.GetData(), A.GetInd(),
                ptrA.GetData(), NULL, NULL, PtrC.GetData(), &nz, &ierr);
    
    nz = PtrC(m)-1;
    Vector<float, VectFull, Alloc1> DataC(nz);
    Vector<int> IndC(nz);
    
    // computes C = A + alpha B
    request = 2;
    mkl_scsradd(&transa, &request, &sort, &m, &n, B.GetData(),
                B.GetInd(), ptrB.GetData(), &alpha, A.GetData(), A.GetInd(),
                ptrA.GetData(), DataC.GetData(), IndC.GetData(),
                PtrC.GetData(), &nz, &ierr);
    
    // switching back to 0-based indexing
    Vector<long> ptrC(m+1);
    for (int i = 0; i <= m; i++)
      ptrC(i) = PtrC(i)-1;
    
    for (int i = 0; i < nz; i++)
      IndC(i)--;

    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]--;
    
    // then B is replaced by C
    B.SetData(m, n, DataC, ptrC, IndC);
  }

  
  //! B = B + alpha A  
  template<class Alloc0, class Alloc1>
  void AddMatrix(const double& alpha,
                 const Matrix<double, General, RowSparse, Alloc0>& A,
                 Matrix<double, General, RowSparse, Alloc1>& B)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (A.GetM() != B.GetM() || A.GetN() != B.GetN())
      throw WrongDim("Add(alpha, const Matrix<RowSparse>& A, "
                     "Matrix<RowSparse>& B)",
                     "Unable to add a " + to_str(A.GetM()) + " x "
                     + to_str(A.GetN()) + " matrix with a "
                     + to_str(B.GetM()) + " x " + to_str(B.GetN())
                     + " matrix.");
#endif

    // switching to 1-based indexing 
    int m = B.GetM(), n = B.GetN();
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    Vector<int> ptrB; B.FillPtrInt(ptrB);
    int* indB = B.GetInd(); int* indA = A.GetInd();
    for (int i = 0; i <= m; i++)
      {
        ptrA(i)++;
        ptrB(i)++;
      }
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]++;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]++;

    char transa('N');
    Vector<int> PtrC(m+1);
    
    // mkl_?csradd is called a first time to know the number of non-zero entries
    int request = 1, sort = 0, nz = 0, ierr = 0; 
    mkl_dcsradd(&transa, &request, &sort, &m, &n, B.GetData(),
                B.GetInd(), ptrB.GetData(), &alpha, A.GetData(), A.GetInd(),
                ptrA.GetData(), NULL, NULL, PtrC.GetData(), &nz, &ierr);
    
    nz = PtrC(m)-1;
    Vector<double, VectFull, Alloc1> DataC(nz);
    Vector<int> IndC(nz);
    
    // computes C = A + alpha B
    request = 2;
    mkl_dcsradd(&transa, &request, &sort, &m, &n, B.GetData(),
                B.GetInd(), ptrB.GetData(), &alpha, A.GetData(), A.GetInd(),
                ptrA.GetData(), DataC.GetData(), IndC.GetData(),
                PtrC.GetData(), &nz, &ierr);
    
    // switching back to 0-based indexing
    Vector<long> ptrC(m+1);
    for (int i = 0; i <= m; i++)
      ptrC(i) = PtrC(i)-1;
    
    for (int i = 0; i < nz; i++)
      IndC(i)--;

    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]--;
    
    // then B is replaced by C
    B.SetData(m, n, DataC, ptrC, IndC);
  }


  //! B = B + alpha A  
  template<class Alloc0, class Alloc1>
  void AddMatrix(const complex<float>& alpha,
                 const Matrix<complex<float>, General, RowSparse, Alloc0>& A,
                 Matrix<complex<float>, General, RowSparse, Alloc1>& B)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (A.GetM() != B.GetM() || A.GetN() != B.GetN())
      throw WrongDim("Add(alpha, const Matrix<RowSparse>& A, "
                     "Matrix<RowSparse>& B)",
                     "Unable to add a " + to_str(A.GetM()) + " x "
                     + to_str(A.GetN()) + " matrix with a "
                     + to_str(B.GetM()) + " x " + to_str(B.GetN())
                     + " matrix.");
#endif

    // switching to 1-based indexing 
    int m = B.GetM(), n = B.GetN();
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    Vector<int> ptrB; B.FillPtrInt(ptrB);
    int* indB = B.GetInd(); int* indA = A.GetInd();
    for (int i = 0; i <= m; i++)
      {
        ptrA(i)++;
        ptrB(i)++;
      }
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]++;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]++;

    char transa('N');
    Vector<int> PtrC(m+1);
    
    // mkl_?csradd is called a first time to know the number of non-zero entries
    int request = 1, sort = 0, nz = 0, ierr = 0; 
    mkl_ccsradd(&transa, &request, &sort, &m, &n, B.GetDataVoid(),
                B.GetInd(), ptrB.GetData(), &alpha, A.GetDataVoid(), A.GetInd(),
                ptrA.GetData(), NULL, NULL, PtrC.GetData(), &nz, &ierr);
    
    nz = PtrC(m)-1;
    Vector<complex<float>, VectFull, Alloc1> DataC(nz);
    Vector<int> IndC(nz);
    
    // computes C = A + alpha B
    request = 2;
    mkl_ccsradd(&transa, &request, &sort, &m, &n, B.GetDataVoid(),
                B.GetInd(), ptrB.GetData(), &alpha, A.GetDataVoid(), A.GetInd(),
                ptrA.GetData(), DataC.GetDataVoid(), IndC.GetData(),
                PtrC.GetData(), &nz, &ierr);
    
    // switching back to 0-based indexing
    Vector<long> ptrC(m+1);
    for (int i = 0; i <= m; i++)
      ptrC(i) = PtrC(i)-1;
    
    for (int i = 0; i < nz; i++)
      IndC(i)--;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]--;
    
    // then B is replaced by C
    B.SetData(m, n, DataC, ptrC, IndC);
  }
  
  
  //! B = B + alpha A
  template<class Alloc0, class Alloc1>
  void AddMatrix(const complex<double>& alpha,
                 const Matrix<complex<double>, General, RowSparse, Alloc0>& A,
                 Matrix<complex<double>, General, RowSparse, Alloc1>& B)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (A.GetM() != B.GetM() || A.GetN() != B.GetN())
      throw WrongDim("Add(alpha, const Matrix<RowSparse>& A, "
                     "Matrix<RowSparse>& B)",
                     "Unable to add a " + to_str(A.GetM()) + " x "
                     + to_str(A.GetN()) + " matrix with a "
                     + to_str(B.GetM()) + " x " + to_str(B.GetN())
                     + " matrix.");
#endif

    // switching to 1-based indexing 
    int m = B.GetM(), n = B.GetN();
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    Vector<int> ptrB; B.FillPtrInt(ptrB);
    int* indB = B.GetInd(); int* indA = A.GetInd();
    for (int i = 0; i <= m; i++)
      {
        ptrA(i)++;
        ptrB(i)++;
      }
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]++;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]++;

    char transa('N');
    Vector<int> PtrC(m+1);
    
    // mkl_?csradd is called a first time to know the number of non-zero entries
    int request = 1, sort = 0, nz = 0, ierr = 0; 
    mkl_zcsradd(&transa, &request, &sort, &m, &n, B.GetDataVoid(),
                B.GetInd(), ptrB.GetData(), &alpha, A.GetDataVoid(), A.GetInd(),
                ptrA.GetData(), NULL, NULL, PtrC.GetData(), &nz, &ierr);
    
    nz = PtrC(m)-1;
    Vector<complex<double>, VectFull, Alloc1> DataC(nz);
    Vector<int> IndC(nz);
    
    // computes C = A + alpha B
    request = 2;
    mkl_zcsradd(&transa, &request, &sort, &m, &n, B.GetDataVoid(),
                B.GetInd(), ptrB.GetData(), &alpha, A.GetDataVoid(), A.GetInd(),
                ptrA.GetData(), DataC.GetDataVoid(), IndC.GetData(),
                PtrC.GetData(), &nz, &ierr);
    
    // switching back to 0-based indexing
    Vector<long> ptrC(m+1);
    for (int i = 0; i <= m; i++)
      ptrC(i) = PtrC(i)-1;
    
    for (int i = 0; i < nz; i++)
      IndC(i)--;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]--;
    
    // then B is replaced by C
    B.SetData(m, n, DataC, ptrC, IndC);
  }

  
  //! C = A B
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltMatrix(const Matrix<float, General, RowSparse, Alloc0>& A,
                 const Matrix<float, General, RowSparse, Alloc1>& B,
                 Matrix<float, General, RowSparse, Alloc2>& C)
  {
    // switching to 1-based indexing 
    int m = A.GetM(), n = A.GetN(), k = B.GetN();
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    Vector<int> ptrB; B.FillPtrInt(ptrB);
    int* indB = B.GetInd(); int* indA = A.GetInd();
    for (int i = 0; i <= m; i++)
      ptrA(i)++;
    
    for (int i = 0; i <= n; i++)
      ptrB(i)++;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]++;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]++;
    
    // mkl_?csrmultcsr is called a first time to know the number
    // of non-zero entries of the matrix C
    char transa('N');
    Vector<int> PtrC(m+1);
    int request = 1, sort = 0, nz = 0, info;
    mkl_scsrmultcsr(&transa, &request, &sort, &m, &n, &k, A.GetData(),
                    A.GetInd(), ptrA.GetData(),
                    B.GetData(), B.GetInd(), ptrB.GetData(),
                    NULL, NULL, PtrC.GetData(), &nz, &info);
    
    // computes product A B
    nz = PtrC(m)-1;
    Vector<float, VectFull, Alloc2> DataC(nz);
    Vector<int> IndC(nz);
    
    request = 2;
    mkl_scsrmultcsr(&transa, &request, &sort, &m, &n, &k, A.GetData(),
                    A.GetInd(), ptrA.GetData(),
                    B.GetData(), B.GetInd(), ptrB.GetData(),
                    DataC.GetData(), IndC.GetData(), PtrC.GetData(), &nz, &info);

    // switching back to 0-based indexing
    Vector<long> ptrC(m+1);
    for (int i = 0; i <= m; i++)
      ptrC(i) = PtrC(i)-1;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]--;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]--;
    
    for (int i = 0; i < nz; i++)
      IndC(i)--;
    
    // then C is filled with values containing A B
    C.SetData(m, k, DataC, ptrC, IndC);    
  }


  //! C = A B
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltMatrix(const Matrix<double, General, RowSparse, Alloc0>& A,
                 const Matrix<double, General, RowSparse, Alloc1>& B,
                 Matrix<double, General, RowSparse, Alloc2>& C)
  {
    // switching to 1-based indexing 
    int m = A.GetM(), n = A.GetN(), k = B.GetN();
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    Vector<int> ptrB; B.FillPtrInt(ptrB);
    int* indB = B.GetInd(); int* indA = A.GetInd();
    for (int i = 0; i <= m; i++)
      ptrA(i)++;
    
    for (int i = 0; i <= n; i++)
      ptrB(i)++;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]++;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]++;
    
    // mkl_?csrmultcsr is called a first time to know the number
    // of non-zero entries of the matrix C
    char transa('N');
    Vector<int> PtrC(m+1);
    int request = 1, sort = 0, nz = 0, info;
    mkl_dcsrmultcsr(&transa, &request, &sort, &m, &n, &k, A.GetData(),
                    A.GetInd(), ptrA.GetData(),
                    B.GetData(), B.GetInd(), ptrB.GetData(),
                    NULL, NULL, PtrC.GetData(), &nz, &info);
    
    // computes product A B
    nz = PtrC(m)-1;
    Vector<double, VectFull, Alloc2> DataC(nz);
    Vector<int> IndC(nz);
    
    request = 2;
    mkl_dcsrmultcsr(&transa, &request, &sort, &m, &n, &k, A.GetData(),
                    A.GetInd(), ptrA.GetData(),
                    B.GetData(), B.GetInd(), ptrB.GetData(),
                    DataC.GetData(), IndC.GetData(), PtrC.GetData(), &nz, &info);

    // switching back to 0-based indexing
    Vector<long> ptrC(m+1);
    for (int i = 0; i <= m; i++)
      ptrC(i) = PtrC(i)-1;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]--;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]--;
    
    for (int i = 0; i < nz; i++)
      IndC(i)--;
    
    // then C is filled with values containing A B
    C.SetData(m, k, DataC, ptrC, IndC);    
  }


  //! C = A B
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltMatrix(const Matrix<complex<float>, General, RowSparse, Alloc0>& A,
                 const Matrix<complex<float>, General, RowSparse, Alloc1>& B,
                 Matrix<complex<float>, General, RowSparse, Alloc2>& C)
  {
    // switching to 1-based indexing 
    int m = A.GetM(), n = A.GetN(), k = B.GetN();
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    Vector<int> ptrB; B.FillPtrInt(ptrB);
    int* indB = B.GetInd(); int* indA = A.GetInd();
    for (int i = 0; i <= m; i++)
      ptrA(i)++;
    
    for (int i = 0; i <= n; i++)
      ptrB(i)++;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]++;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]++;
    
    // mkl_?csrmultcsr is called a first time to know the number
    // of non-zero entries of the matrix C
    char transa('N');
    Vector<int> PtrC(m+1);
    int request = 1, sort = 0, nz = 0, info;
    mkl_ccsrmultcsr(&transa, &request, &sort, &m, &n, &k, A.GetDataVoid(),
                    A.GetInd(), ptrA.GetData(), B.GetDataVoid(), B.GetInd(),
                    ptrB.GetData(), NULL, NULL, PtrC.GetData(), &nz, &info);
    
    // computes product A B
    nz = PtrC(m)-1;
    Vector<complex<float>, VectFull, Alloc2> DataC(nz);
    Vector<int> IndC(nz);
    
    request = 2;
    mkl_ccsrmultcsr(&transa, &request, &sort, &m, &n, &k, A.GetDataVoid(),
                    A.GetInd(), ptrA.GetData(), B.GetDataVoid(), B.GetInd(), ptrB.GetData(),
                    DataC.GetDataVoid(), IndC.GetData(), PtrC.GetData(), &nz, &info);

    // switching back to 0-based indexing
    Vector<long> ptrC(m+1);
    for (int i = 0; i <= m; i++)
      ptrC(i) = PtrC(i)-1;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]--;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]--;
    
    for (int i = 0; i < nz; i++)
      IndC(i)--;
    
    // then C is filled with values containing A B
    C.SetData(m, k, DataC, ptrC, IndC);    
  }


  //! C = A B
  template<class Alloc0, class Alloc1, class Alloc2>
  void MltMatrix(const Matrix<complex<double>, General, RowSparse, Alloc0>& A,
                 const Matrix<complex<double>, General, RowSparse, Alloc1>& B,
                 Matrix<complex<double>, General, RowSparse, Alloc2>& C)
  {
    // switching to 1-based indexing 
    int m = A.GetM(), n = A.GetN(), k = B.GetN();
    Vector<int> ptrA; A.FillPtrInt(ptrA);
    Vector<int> ptrB; B.FillPtrInt(ptrB);
    int* indB = B.GetInd(); int* indA = A.GetInd();
    for (int i = 0; i <= m; i++)
      ptrA(i)++;
    
    for (int i = 0; i <= n; i++)
      ptrB(i)++;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]++;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]++;
    
    // mkl_?csrmultcsr is called a first time to know the number
    // of non-zero entries of the matrix C
    char transa('N');
    Vector<int> PtrC(m+1);
    int request = 1, sort = 0, nz = 0, info;
    mkl_zcsrmultcsr(&transa, &request, &sort, &m, &n, &k, A.GetDataVoid(),
                    A.GetInd(), ptrA.GetData(), B.GetDataVoid(), B.GetInd(),
                    ptrB.GetData(), NULL, NULL, PtrC.GetData(), &nz, &info);
    
    // computes product A B
    nz = PtrC(m)-1;
    Vector<complex<double>, VectFull, Alloc2> DataC(nz);
    Vector<int> IndC(nz);
    
    request = 2;
    mkl_zcsrmultcsr(&transa, &request, &sort, &m, &n, &k, A.GetDataVoid(),
                    A.GetInd(), ptrA.GetData(), B.GetDataVoid(), B.GetInd(), ptrB.GetData(),
                    DataC.GetDataVoid(), IndC.GetData(), PtrC.GetData(), &nz, &info);

    // switching back to 0-based indexing
    Vector<long> ptrC(m+1);
    for (int i = 0; i <= m; i++)
      ptrC(i) = PtrC(i)-1;
    
    for (int i = 0; i < A.GetIndSize(); i++)
      indA[i]--;

    for (int i = 0; i < B.GetIndSize(); i++)
      indB[i]--;
    
    for (int i = 0; i < nz; i++)
      IndC(i)--;

    // then C is filled with values containing A B
    C.SetData(m, k, DataC, ptrC, IndC);    
  }
  
} // end namespace

#define SELDON_FILE_MKL_SPARSE_CXX
#endif
