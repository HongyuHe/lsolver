#ifndef SELDON_FILE_FUNCTIONS_BASE_HXX

// in this file, we put functions such as Add, Mlt, MltAdd
// with a reduced number of templates in order
// to forbid undesired instantiations of generic functions

namespace Seldon
{
  
  // Inline functions

  // Mlt with a scalar
  template<class T, class Allocator>
  void Mlt(const T& alpha, Array3D<T, Allocator>& A);

  template<class T, class Allocator>
  void Mlt(const T& alpha, Array3D<complex<T>, Allocator>& A);
  
#ifdef SELDON_WITH_REDUCED_TEMPLATE
  template<class T, class Storage, class Allocator>
  void Mlt(const T& alpha, Vector<T, Storage, Allocator>& A);
#else
  template<class T, class T0, class Storage, class Allocator>
  void Mlt(const T& alpha, Vector<T0, Storage, Allocator>& A);
#endif
  
  template<class T, class Storage, class Allocator>
  void Mlt(const T& alpha, Vector<complex<T>, Storage, Allocator>& A);
    
  // Add with two vectors
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  void Add(const T& alpha, const Vector<T, Storage1, Allocator1>& X,
	   Vector<T, Storage2, Allocator2>& Y);

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  void Add(const T& alpha, const Vector<complex<T>, Storage1, Allocator1>& X,
	   Vector<complex<T>, Storage2, Allocator2>& Y);

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  void Add(const T& alpha, const Vector<T, Storage1, Allocator1>& X,
	   const T& beta, Vector<T, Storage2, Allocator2>& Y);
  
  // Copy of two vectors
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  void Copy(const Vector<T, Storage1, Allocator1>& X,
	    Vector<T, Storage2, Allocator2>& Y);

  // Scalar product between vectors
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  T DotProd(const Vector<T, Storage1, Allocator1>& X,
	    const Vector<T, Storage2, Allocator2>& Y);
  
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  T DotProdConj(const Vector<T, Storage1, Allocator1>& X,
		const Vector<T, Storage2, Allocator2>& Y);

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  complex<T> DotProdConj(const Vector<complex<T>, Storage1, Allocator1>& X,
			 const Vector<complex<T>, Storage2, Allocator2>& Y);
    
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  complex<T> DotProdConj(const Vector<complex<T>, Storage1, Allocator1>& X,
			 const Vector<T, Storage2, Allocator2>& Y);

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  complex<T> DotProdConj(const Vector<T, Storage1, Allocator1>& X,
			 const Vector<complex<T>, Storage2, Allocator2>& Y);

  template<class T, class E>
  typename ClassComplexType<T>::Treal Norm1(const VectorExpression<T, E>& X);
  
  template<class T, class E>
  typename ClassComplexType<T>::Treal Norm2(const VectorExpression<T, E>& X);

  // matrix-vector product Mlt
#ifdef SELDON_WITH_REDUCED_TEMPLATE
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Mlt(const Matrix<T, Prop0, Storage0, Allocator0>& M,
	   const Vector<T, Storage1, Allocator1>& X,
	   Vector<T, Storage2, Allocator2>& Y);
#else
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	   const Vector<T1, Storage1, Allocator1>& X,
	   Vector<T2, Storage2, Allocator2>& Y);
#endif
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Mlt(const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
	   const Vector<T, Storage1, Allocator1>& X,
	   Vector<T, Storage2, Allocator2>& Y);
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Mlt(const Matrix<T, Prop0, Storage0, Allocator0>& M,
	   const Vector<complex<T>, Storage1, Allocator1>& X,
	   Vector<complex<T>, Storage2, Allocator2>& Y);
  
#ifdef SELDON_WITH_REDUCED_TEMPLATE
  template <class T, class Prop1, class Storage1, class Allocator1,
	    class Storage2, class Allocator2,
	    class Storage3, class Allocator3>
  void Mlt(const T& alpha,
	   const Matrix<T, Prop1, Storage1, Allocator1>& M,
	   const Vector<T, Storage2, Allocator2>& X,
	   Vector<T, Storage3, Allocator3>& Y);
#else
  template <class T0,
            class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T0& alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y);
#endif
  
  template <class T, class Prop1, class Storage1, class Allocator1,
	    class Storage2, class Allocator2,
	    class Storage3, class Allocator3>
  void Mlt(const complex<T>& alpha,
	   const Matrix<T, Prop1, Storage1, Allocator1>& M,
	   const Vector<complex<T>, Storage2, Allocator2>& X,
	   Vector<complex<T>, Storage3, Allocator3>& Y);
  
  template <class T, class Prop1, class Storage1, class Allocator1,
	    class Storage2, class Allocator2,
	    class Storage3, class Allocator3>
  void Mlt(const T& alpha,
	   const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	   const Vector<T, Storage2, Allocator2>& X,
	   Vector<T, Storage3, Allocator3>& Y);

  // matrix vector product MltAdd
#ifdef SELDON_WITH_REDUCED_TEMPLATE
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<T, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<T, Storage4, Allocator4>& Y);
#else
  template<class T,
	   class T1, class Prop1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
           class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T& alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3& beta,
	      Vector<T4, Storage4, Allocator4>& Y);
#endif

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha,
	      const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y);
  
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y);

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha,
	      const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	      const Vector<T, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<T, Storage4, Allocator4>& Y);
  
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const complex<T>& alpha,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const complex<T>& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y);

  // SOR relaxation
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const Matrix<T, Prop0, Storage0, Allocator0>& M,
	   Vector<T, Storage2, Allocator2>& Y,
	   const Vector<T, Storage1, Allocator1>& X,
	   const typename ClassComplexType<T>::Treal& omega,
	   int iter, int type_ssor = 2);
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const SeldonTranspose&,
	   const Matrix<T, Prop0, Storage0, Allocator0>& M,
	   Vector<T, Storage2, Allocator2>& Y,
	   const Vector<T, Storage1, Allocator1>& X,
	   const typename ClassComplexType<T>::Treal& omega, int iter, int type_ssor = 3);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const Matrix<T, Prop0, Storage0, Allocator0>& M,
	   Vector<complex<T>, Storage2, Allocator2>& Y,
	   const Vector<complex<T>, Storage1, Allocator1>& X,
	   const typename ClassComplexType<T>::Treal& omega,
	   int iter, int type_ssor = 2);
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const SeldonTranspose&,
	   const Matrix<T, Prop0, Storage0, Allocator0>& M,
	   Vector<complex<T>, Storage2, Allocator2>& Y,
	   const Vector<complex<T>, Storage1, Allocator1>& X,
	   const typename ClassComplexType<T>::Treal& omega, int iter, int type_ssor = 3);
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
	   Vector<T, Storage2, Allocator2>& Y,
	   const Vector<T, Storage1, Allocator1>& X,
	   const typename ClassComplexType<T>::Treal& omega,
	   int iter, int type_ssor = 2);
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const SeldonTranspose&,
	   const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
	   Vector<T, Storage2, Allocator2>& Y,
	   const Vector<T, Storage1, Allocator1>& X,
	   const typename ClassComplexType<T>::Treal& omega, int iter, int type_ssor = 3);

  // triangular solution
  template <class T, class Prop0, class Allocator0,
            class Allocator1>
  void
  Solve(const SeldonUplo& Uplo,
	const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<T, Prop0, RowSparse, Allocator0>& A,
	const Vector<T, VectFull, Allocator1>& X,
	Vector<T, VectFull, Allocator1>& Y);
  
  // SolveLU without pivot
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  void SolveLU(const Matrix<T, Prop0, Storage0, Allocator0>& M,
	       Vector<T, Storage1, Allocator1>& Y);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  void SolveLU(const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
	       Vector<T, Storage1, Allocator1>& Y);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  void SolveLU(const Matrix<T, Prop0, Storage0, Allocator0>& M,
	       Vector<complex<T>, Storage1, Allocator1>& Y);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  void SolveLU(const SeldonTranspose& transA,
	       const Matrix<T, Prop0, Storage0, Allocator0>& M,
	       Vector<T, Storage1, Allocator1>& Y);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  void SolveLU(const SeldonTranspose& transA,
	       const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
	       Vector<T, Storage1, Allocator1>& Y);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  void SolveLU(const SeldonTranspose& transA,
	       const Matrix<T, Prop0, Storage0, Allocator0>& M,
	       Vector<complex<T>, Storage1, Allocator1>& Y);
  
  // SolveLU with pivot
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<T0>& x);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void SolveLU(const Matrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<T0>& x);

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<T0, General, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<T0>& x);

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<complex<T0>, General, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<T0>& x);

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<T0, Symmetric, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<T0>& x);

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<complex<T0>, Symmetric, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<T0>& x);
  
  // product of a scalar with a matrix
  template<class T, class Prop, class Storage, class Allocator>
  void Mlt(const T& alpha, Matrix<T, Prop, Storage, Allocator>& A);
  
  template<class T, class Prop, class Storage, class Allocator>
  void Mlt(const T& alpha, Matrix<complex<T>, Prop, Storage, Allocator>& A);
  
  // matrix matrix product Mlt
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2,
	    class Prop3, class Storage3, class Allocator3>
  void Mlt(const T& alpha,
	   const Matrix<T, Prop1, Storage1, Allocator1>& A,
	   const Matrix<T, Prop2, Storage2, Allocator2>& B,
	   Matrix<T, Prop3, Storage3, Allocator3>& C);
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Mlt(const Matrix<T, Prop0, Storage0, Allocator0>& A,
	   const Matrix<T, Prop1, Storage1, Allocator1>& B,
	   Matrix<T, Prop2, Storage2, Allocator2>& C);
  
  // matrix matrix product MltAdd
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2,
	    class Prop3, class Storage3, class Allocator3>
  void MltAdd(const T& alpha,
	      const Matrix<T, Prop1, Storage1, Allocator1>& A,
	      const Matrix<T, Prop2, Storage2, Allocator2>& B,
	      const T& beta,
	      Matrix<T, Prop3, Storage3, Allocator3>& C);
  
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2,
	    class Prop3, class Storage3, class Allocator3>
  void MltAdd(const T& alpha, const SeldonTranspose& transA,
	      const Matrix<T, Prop1, Storage1, Allocator1>& A,
	      const SeldonTranspose& transB,
	      const Matrix<T, Prop2, Storage2, Allocator2>& B,
	      const T& beta,
	      Matrix<T, Prop3, Storage3, Allocator3>& C);

  // Adds two matrices
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Add(const T& alpha,
	   const Matrix<T, Prop1, Storage1, Allocator1>& A,
	   Matrix<T, Prop2, Storage2, Allocator2>& B);

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Add(const complex<T>& alpha,
	   const Matrix<T, Prop1, Storage1, Allocator1>& A,
	   Matrix<complex<T>, Prop2, Storage2, Allocator2>& B);

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Add(const T& alpha,
	   const Matrix<complex<T>, Prop1, Storage1, Allocator1>& A,
	   Matrix<complex<T>, Prop2, Storage2, Allocator2>& B);
  
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Add(const T& alpha,
	   const Matrix<complex<T>, Prop1, Storage1, Allocator1>& A,
	   Matrix<T, Prop2, Storage2, Allocator2>& B);

  // Copies two matrices
  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Prop2, class Storage2, class Allocator2>
  void Copy(const Matrix<T, Prop1, Storage1, Allocator1>& A,
	    Matrix<T, Prop2, Storage2, Allocator2>& B);

  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Prop2, class Storage2, class Allocator2>
  void Copy(const Matrix<complex<T>, Prop1, Storage1, Allocator1>& A,
	    Matrix<T, Prop2, Storage2, Allocator2>& B);

  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Prop2, class Storage2, class Allocator2>
  void Copy(const Matrix<T, Prop1, Storage1, Allocator1>& A,
	    Matrix<complex<T>, Prop2, Storage2, Allocator2>& B);
  
  template<class T, class Prop, class Storage, class Allocator>
  bool IsSymmetricMatrix(const Matrix<T, Prop, Storage, Allocator>& A);

  template<class T, class Storage, class Allocator>
  bool IsSymmetricMatrix(const Matrix<T, Symmetric, Storage, Allocator>& A);
    
  template<class T, class Prop, class Storage, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, Storage, Allocator>& A);

  template<class T, class Prop, class Storage, class Allocator>
  bool IsComplexMatrix(const Matrix<complex<T>, Prop, Storage, Allocator>& A);

  // non-inline
  
  // matrix vector product Mlt and MltAdd
  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  void Mlt(const SeldonTranspose& Trans,
	   const Matrix<T, Prop1, Storage1, Allocator1>& M,
	   const Vector<T, Storage2, Allocator2>& X,
	   Vector<T, Storage3, Allocator3>& Y);

  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  void Mlt(const SeldonTranspose& Trans,
	   const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	   const Vector<T, Storage2, Allocator2>& X,
	   Vector<T, Storage3, Allocator3>& Y);

  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  void Mlt(const SeldonTranspose& Trans,
	   const Matrix<T, Prop1, Storage1, Allocator1>& M,
	   const Vector<complex<T>, Storage2, Allocator2>& X,
	   Vector<complex<T>, Storage3, Allocator3>& Y);
  
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha, const SeldonTranspose& Trans,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<T, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<T, Storage4, Allocator4>& Y);

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha, const SeldonTranspose& Trans,
	      const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y);
  
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha, const SeldonTranspose& Trans,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y);

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha, const SeldonTranspose& Trans,
	      const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	      const Vector<T, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<T, Storage4, Allocator4>& Y);

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const complex<T>& alpha, const SeldonTranspose& Trans,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const complex<T>& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y);

  // SolveLU with pivot
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x);
  
  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<T0, General, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x);

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<T0, Symmetric, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x);
  
}

#define SELDON_FILE_FUNCTIONS_BASE_HXX
#endif
