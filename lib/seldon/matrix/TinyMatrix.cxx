#ifndef MONTJOIE_FILE_TINY_MATRIX_CXX

#include "TinyMatrix.hxx"

namespace Seldon
{

    /****************
   * 2x2 matrices *
   ****************/
  
  
  //! B = A^-1
  template<class T>
  void GetInverse(const TinyMatrix<T, General, 2, 2> & A,
		  TinyMatrix<T, General, 2, 2> & B)
  {
    T one; SetComplexOne(one);
    T d = Det(A);
    T invDet = one/d;
    B(0,0) = A(1,1)*invDet;
    B(1,0) = -A(1,0)*invDet;
    B(0,1) = -A(0,1)*invDet;
    B(1,1) = A(0,0)*invDet;
  }
  
  
  //! replaces B by its inverse
  template<class T>
  void GetInverse(TinyMatrix<T, General, 2, 2> & B)
  {
    T one; SetComplexOne(one);
    T d = Det(B);
    T invDet = one/d;
    d = B(0,0);
    B(0,0) = B(1,1)*invDet;
    B(1,0) *= -invDet;
    B(0,1) *= -invDet;
    B(1,1) = d*invDet;
  }
  
  
  //! B = A^-1
  template<class T>
  void GetInverse(const TinyMatrix<T, Symmetric, 2, 2> & A,
		  TinyMatrix<T, Symmetric, 2, 2> & B)
  {
    T one; SetComplexOne(one);
    T d = Det(A);
    T invDet = one/d;
    B(0,0) = A(1,1)*invDet;
    B(0,1) = -A(0,1)*invDet;
    B(1,1) = A(0,0)*invDet;
  }
  
  
  //! replaces B by its inverse
  template<class T>
  void GetInverse(TinyMatrix<T, Symmetric, 2, 2> & B)
  {
    T one; SetComplexOne(one);
    T d = Det(B);
    T invDet = one/d;
    d = B(0,0);
    B(0,0) = B(1,1)*invDet;
    B(0,1) *= -invDet;
    B(1,1) = d*invDet;
  }
  

  //! computes eigenvalues of A (real matrix)
  template<class T>
  void GetEigenvalues(TinyMatrix<T, General, 2, 2>& A,
		      TinyVector<T, 2> & LambdaR, TinyVector<T, 2>& LambdaI)
  {
    T trA = A(0, 0) + A(1, 1);
    T delta = trA*trA - 4.0*(A(0,0)*A(1,1) - A(0,1)*A(1,0));
    if (delta < 0)
      {
	delta = sqrt(-delta);
	LambdaR(0) = 0.5*trA; LambdaR(1) = 0.5*trA;
	LambdaI(0) = 0.5*delta; LambdaI(1) = -0.5*delta;
      }
    else
      {
	delta = sqrt(delta);
	LambdaR(0) = 0.5*(trA - delta);
	LambdaR(1) = 0.5*(trA + delta);
	LambdaI(0) = 0; LambdaI(1) = 0;
      }
  }


  //! computes eigenvalues of A (complex matrix)
  template<class T>
  void GetEigenvalues(TinyMatrix<complex<T>, General, 2, 2>& A,
		      TinyVector<complex<T>, 2> & Lambda)
  {
    complex<T> trA = A(0, 0) + A(1, 1);
    complex<T> delta = trA - 4.0*(A(0,0)*A(1,1) - A(0,1)*A(1,0));
    delta = sqrt(delta);
    Lambda(0) = 0.5*(trA - delta);
    Lambda(1) = 0.5*(trA + delta);
  }
  
  
  //! replaces A by its square root
  template<class T>
  void GetSquareRoot(TinyMatrix<T, Symmetric, 2, 2>& A)
  {
    T a = A(0,0); T b = A(1,1); T c = A(1,0);
    // discrimant of equation x^2 - (a+b) x + ab - c^2
    T delta = sqrt((a-b)*(a-b) + 4.0*c*c);
    // the two eigenvalues
    T l1 = T(0.5)*(a+b-delta);  T l2 = T(0.5)*(a+b+delta); 
    // then square root is equal to
    // 1/(sqrt(l1)+sqrt(l2)) [ a+sqrt(l1*l2) , c ; c, b+sqrt(l1*l2)]
    T root_l1 = sqrt(l1); T root_l2 = sqrt(l2);
    // l1 and l2 are respectively inverse of the sum of eigenvalues'
    // square root, and product of eigenvalues' square root
    l1 = 1.0/(root_l1 + root_l2); l2 = root_l1*root_l2;
    A(0,0) = (a+l2)*l1; A(1,0) = c*l1; A(1,1) = (b+l2)*l1;
  }
  
  
  /****************
   * 3x3 matrices *
   ****************/
  
  
  //! returns determinant of matrix A
  template<class T, class Prop>
  T Det(const TinyMatrix<T, Prop, 3, 3> & A)
  {
    T d = A(0,0)*(A(1,1)*A(2,2) - A(2,1)*A(1,2))
      - A(0,1)*(A(1,0)*A(2,2) - A(2,0)*A(1,2))
      + A(0,2)*(A(1,0)*A(2,1) - A(2,0)*A(1,1));
    
    return d;
  }

  
  //! B = A^-1
  template<class T>
  void GetInverse(const TinyMatrix<T, General, 3, 3> & A,
		  TinyMatrix<T, General, 3, 3> & B)
  {
    T one; SetComplexOne(one);
    T d = Det(A);
    T invDet = one/d;
    B(0,0) = (A(1,1)*A(2,2)-A(2,1)*A(1,2))*invDet;
    B(0,1) = -(A(0,1)*A(2,2)-A(2,1)*A(0,2))*invDet;
    B(0,2) = (A(0,1)*A(1,2)-A(1,1)*A(0,2))*invDet;
    B(1,0) = -(A(1,0)*A(2,2)-A(2,0)*A(1,2))*invDet;
    B(1,1) = (A(0,0)*A(2,2)-A(2,0)*A(0,2))*invDet;
    B(1,2) = -(A(0,0)*A(1,2)-A(1,0)*A(0,2))*invDet;
    B(2,0) = (A(1,0)*A(2,1)-A(2,0)*A(1,1))*invDet;
    B(2,1) = -(A(0,0)*A(2,1)-A(2,0)*A(0,1))*invDet;
    B(2,2) = (A(0,0)*A(1,1)-A(1,0)*A(0,1))*invDet;
  }
  
  
  //! replaces B by its inverse
  template<class T>
  void GetInverse(TinyMatrix<T, General, 3, 3> & B)
  {
    TinyMatrix<T, General, 3, 3> A = B;
    GetInverse(A, B);
  }
  
  
  //! B = A^-1
  template<class T>
  void GetInverse(const TinyMatrix<T, Symmetric, 3, 3> & A,
		  TinyMatrix<T, Symmetric, 3, 3> & B)
  {
    T one; SetComplexOne(one);
    T d = Det(A);
    T invDet = one/d;
    B(0,0) = (A(1,1)*A(2,2)-A(2,1)*A(1,2))*invDet;
    B(0,1) = -(A(0,1)*A(2,2)-A(2,1)*A(0,2))*invDet;
    B(0,2) = (A(0,1)*A(1,2)-A(1,1)*A(0,2))*invDet;
    B(1,1) = (A(0,0)*A(2,2)-A(2,0)*A(0,2))*invDet;
    B(1,2) = -(A(0,0)*A(1,2)-A(1,0)*A(0,2))*invDet;
    B(2,2) = (A(0,0)*A(1,1)-A(1,0)*A(0,1))*invDet;
  }
  
  
  //! replaces B by its inverse
  template<class T>
  void GetInverse(TinyMatrix<T, Symmetric, 3, 3> & B)
  {
    TinyMatrix<T, Symmetric, 3, 3> A = B;
    GetInverse(A, B);
  }
  
  
  //! computation of I - n n^T where n is the normale
  //! this is the projector to tangential plane
  template<class T0, class T1>
  void GetTangentialProjector(const TinyVector<T0, 3>& n,
			      TinyMatrix<T1, Symmetric, 3, 3>& P)
  {
    P(0, 0) = n(1)*n(1) + n(2)*n(2);
    P(0, 1) = -n(0)*n(1);
    P(0, 2) = -n(0)*n(2);
    P(1, 1) = n(0)*n(0) + n(2)*n(2);
    P(1, 2) = -n(1)*n(2);
    P(2, 2) = n(0)*n(0) + n(1)*n(1);
  }

  
  //! computation of n n^T where n is the normale
  //! this is the projector to normale
  template<class T0, class T1>
  void GetNormalProjector(const TinyVector<T0, 3>& n,
			  TinyMatrix<T1, Symmetric, 3, 3>& P)
  {
    P(0, 0) = n(0)*n(0);
    P(0, 1) = n(0)*n(1);
    P(0, 2) = n(0)*n(2);
    P(1, 1) = n(1)*n(1);
    P(1, 2) = n(1)*n(2);
    P(2, 2) = n(2)*n(2);
  }


  //! computation of n n^T where n is the normale
  //! this is the projector to normale
  template<class T0, class T1>
  void GetNormalProjector(const TinyVector<T0, 3>& n,
			  TinyMatrix<T1, General, 3, 3>& P)
  {
    P(0, 0) = n(0)*n(0);
    P(0, 1) = n(0)*n(1);
    P(0, 2) = n(0)*n(2);
    P(1, 1) = n(1)*n(1);
    P(1, 2) = n(1)*n(2);
    P(2, 2) = n(2)*n(2);
    P(1, 0) = P(0, 1);
    P(2, 0) = P(0, 2);
    P(2, 1) = P(1, 2);
  }


  /**************************
   * Functions for any size *
   **************************/
  
  
  //! returns norm of a column of matrix A
  template<class T, class Property, int m, int n>
  typename ClassComplexType<T>::Treal
  Norm2_Column(const TinyMatrix<T, Property, m, n>& A,
	       int first_row, int index_col)
  {
    typename ClassComplexType<T>::Treal sum(0);
    for (int i = first_row; i < m; i++)
      sum += absSquare(A(i, index_col));

    sum = sqrt(sum);
    return sum;
  }


  //! Returns 2-norm of a column of matrix A
  template<class T, class Property, class Storage, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm2_Column(const Matrix<T, Property, Storage, Allocator>& A,
	       int first_row, int index_col)
  {
    typename ClassComplexType<T>::Treal sum(0);
    for (int i = first_row; i < A.GetM(); i++)
      sum += absSquare(A(i, index_col));
    
    sum = sqrt(sum);
    return sum;
  }

    
  //! A is replaced by its inverse
  /*!
    LU factorisation is performed, then inverse is formed
    All the loops are unrolled, so that this function is very efficient
    for small values of m, more efficient than Lapack when m < 10.
    Partial pivoting is used so that this algorithm always works if A
    is invertible.
   */
  template<class T, int m>
  void GetInverse(TinyMatrix<T, General, m, m>& A)
  {
    // storing permutation for Gauss pivot
    TinyVector<int, m> pivot;
    TinyMatrixLoop<m-1>::PivotGauss(A, pivot);
    
    // dividing last row by diagonal coefficient
    T one; SetComplexOne(one);
    T coef = one/A(m-1, m-1);
    TinyMatrixDoubleLoop<m-1, m>::MltRow(A, coef);
    A(m-1, m-1) = coef;
    
    // inverting by upper matrix U
    TinyMatrixLoop<m-1>::SolveUpper(A);
    
    // permuting columns with pivot
    TinyMatrixLoop<m-1>::PermuteColumn(A, pivot);
  }
  
  
  //! B = A^-1
  template<class T, int m>
  void GetInverse(const TinyMatrix<T, General, m, m>& A,
		  TinyMatrix<T, General, m, m>& B)
  {
    B = A;
    GetInverse(B);
  }
  
  
  //! replaces A by its inverse
  /*!
    For symmetric matrices, we use the unsymmetric algorithm.
    It should be more efficient to write a specified algorithm
    for symmetric matrices.
   */
  template<class T, int m>
  void GetInverse(TinyMatrix<T, Symmetric, m, m> & A)
  {
    TinyMatrix<T, General, m, m> B;
    B = A;
    GetInverse(B);
    A = B;
  }
  
  
  //! B = A^-1
  template<class T, int m>
  void GetInverse(const TinyMatrix<T, Symmetric, m, m>& A,
		  TinyMatrix<T, Symmetric, m, m>& B)
  {
    B = A;
    GetInverse(B);
  }

  
  /***************************
   * GetCholesky for any size *
   ***************************/
  
  
  //! replaces A by its Cholesky factorisation
  template<class T, int m>
  void GetCholesky(TinyMatrix<T, Symmetric, m, m>& A)
  {
    TinyMatrixLoop<m>::GetCholesky(A);
  }

  
  //! solves L x = b, assuming that GetCholesky has been called
  template<class T, class T2, int m>
  void SolveCholesky(const class_SeldonNoTrans& trans,
                     const TinyMatrix<T, Symmetric, m, m>& A,
                     TinyVector<T2, m>& x)
  {
    TinyMatrixLoop<m>::SolveCholesky(trans, A, x);
  }


  //! solves L^T x = b, assuming that GetCholesky has been called
  template<class T, class T2, int m>
  void SolveCholesky(const class_SeldonTrans& trans,
                     const TinyMatrix<T, Symmetric, m, m>& A,
                     TinyVector<T2, m>& x)
  {
    TinyMatrixLoop<m>::SolveCholesky(trans, A, x);
  }


  //! replaces x by L x, assuming that GetCholesky has been called
  template<class T, class T2, int m>
  void MltCholesky(const class_SeldonNoTrans& trans,
		   const TinyMatrix<T, Symmetric, m, m>& A,
		   TinyVector<T2, m>& x)
  {
    TinyMatrixLoop<m>::MltCholesky(trans, A, x);
  }


  //! replaces x by L^T x, assuming that GetCholesky has been called
  template<class T, class T2, int m>
  void MltCholesky(const class_SeldonTrans& trans,
		   const TinyMatrix<T, Symmetric, m, m>& A,
		   TinyVector<T2, m>& x)
  {
    TinyMatrixLoop<m>::MltCholesky(trans, A, x);
  }
  
  
  /********************************************************
   * Eigenvalues, Square root, absolute value of matrices *
   ********************************************************/
  
  
  //! computes eigenvalues and eigenvectors of matrix A
  template<int m, class T>
  void GetEigenvaluesEigenvectors(TinyMatrix<T, Symmetric, m, m>& A,
				  TinyVector<T, m>& w,
				  TinyMatrix<T, General, m, m>& z)
  {
    Matrix<T, Symmetric, RowSymPacked> A2(m, m);
    Matrix<T, General, RowMajor> z2(m, m);
    Vector<T> w2(m);
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
	A2(i, j) = A(i, j);
    
    GetEigenvaluesEigenvectors(A2, w2, z2); 

    for (int i = 0; i < m; i++)
      {
	w(i) = w2(i);
	
	for (int j = 0; j < m; j++)
	  z(i, j) = z2(i, j);
      }
  }


  //! computes eigenvalues of matrix A
  template<int m, class T>
  void GetEigenvalues(TinyMatrix<T, Symmetric, m, m>& A,
		      TinyVector<T, m>& w)
  {
    Matrix<T, Symmetric, RowSymPacked> A2(m, m);
    Vector<T> w2(m);
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
	A2(i, j) = A(i, j);
    
    GetEigenvalues(A2, w2); 
    
    for (int i = 0; i < m; i++)
      w(i) = w2(i);
  }
  
  
  //! replaces A by its square root
  template<class T, int m>
  void GetSquareRoot(TinyMatrix<T, Symmetric, m, m>& A)
  {
    TinyVector<T, m> w;
    TinyMatrix<T, General, m, m> z, inv_z;
    GetEigenvaluesEigenvectors(A, w, z);
        
    // eigenvectors orthonormals, inverse of z is its transpose
    Transpose(z, inv_z);
    // we compute P D^{1/2} P^{-1}  
    for (int i = 0; i < m; i++)
      {
	w(i) = sqrt(w(i));
	for (int j = 0; j < m; j++)
	  inv_z(i, j) *= w(i);
      }
    
    Mlt(z, inv_z, A);
  }
  
} // end namespace

#define MONTJOIE_FILE_TINY_MATRIX_CXX
#endif
