#ifndef SELDON_FILE_FUNCTIONS_BASE_INLINE_CXX

// in this file, we put functions such as Add, Mlt, MltAdd
// with a reduced number of templates in order
// to forbid undesired instantiations of generic functions

/*
  Functions defined in this file:
  
  M X -> Y
  Mlt(M, X, Y)

  alpha M X -> Y
  Mlt(alpha, M, X, Y)

  M X -> Y
  M^T X -> Y
  Mlt(Trans, M, X, Y)

*/

namespace Seldon
{
  
  /************************
   * Functions in Array3D *
   ************************/
  
  template<class T, class Allocator>
  inline void Mlt(const T& alpha, Array3D<T, Allocator>& A)
  {
    MltScalar(alpha, A);
  }

  template<class T, class Allocator>
  inline void Mlt(const T& alpha, Array3D<complex<T>, Allocator>& A)
  {
    MltScalar(alpha, A);
  }
  
  /*********************************
   * Functions in Functions_Vector *
   *********************************/

#ifdef SELDON_WITH_REDUCED_TEMPLATE
  template<class T, class Storage, class Allocator>
  inline void Mlt(const T& alpha, Vector<T, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }
#else
  template<class T, class T0, class Storage, class Allocator>
  inline void Mlt(const T& alpha, Vector<T0, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }
#endif

  template<class T, class Storage, class Allocator>
  inline void Mlt(const T& alpha, Vector<complex<T>, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }
  
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline void Add(const T& alpha, const Vector<T, Storage1, Allocator1>& X,
		  Vector<T, Storage2, Allocator2>& Y)
  {
    AddVector(alpha, X, Y);
  }

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline void Add(const T& alpha, const Vector<complex<T>, Storage1, Allocator1>& X,
		  Vector<complex<T>, Storage2, Allocator2>& Y)
  {
    AddVector(alpha, X, Y);
  }

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline void Add(const T& alpha, const Vector<T, Storage1, Allocator1>& X,
		  const T& beta, Vector<T, Storage2, Allocator2>& Y)
  {
    AddVector(alpha, X, beta, Y);
  }
  
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline void Copy(const Vector<T, Storage1, Allocator1>& X,
		   Vector<T, Storage2, Allocator2>& Y)
  {
    CopyVector(X, Y);
  }

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline T DotProd(const Vector<T, Storage1, Allocator1>& X,
		   const Vector<T, Storage2, Allocator2>& Y)
  {
    return DotProdVector(X, Y);
  }

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline T DotProdConj(const Vector<T, Storage1, Allocator1>& X,
		       const Vector<T, Storage2, Allocator2>& Y)
  {
    return DotProdVector(X, Y);
  }

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline complex<T> DotProdConj(const Vector<complex<T>, Storage1, Allocator1>& X,
				const Vector<complex<T>, Storage2, Allocator2>& Y)
  {
    return DotProdConjVector(X, Y);
  }

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline complex<T> DotProdConj(const Vector<complex<T>, Storage1, Allocator1>& X,
				const Vector<T, Storage2, Allocator2>& Y)
  {
    return DotProdConjVector(X, Y);
  }

  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2>
  inline complex<T> DotProdConj(const Vector<T, Storage1, Allocator1>& X,
				const Vector<complex<T>, Storage2, Allocator2>& Y)
  {
    return DotProdVector(X, Y);
  }

  //! returns 1-norm of an expression with vectors
  template<class T, class E>
  inline typename ClassComplexType<T>::Treal
  Norm1(const VectorExpression<T, E>& X)
  {
    typename ClassComplexType<T>::Treal value(0);
    
    for (long i = 0; i < X.GetSize(); i++)
      value += ComplexAbs(X(i));
    
    return value;
  }

  //! returns 2-norm of an expression with vectors
  template<class T, class E>
  typename ClassComplexType<T>::Treal
  Norm2(const VectorExpression<T, E>& X)
  {
    typename ClassComplexType<T>::Treal value(0);
    
    for (long i = 0; i < X.GetSize(); i++)
      value += absSquare(X(i));
    
    return sqrt(value);
  }


  /**********************************
   * Functions in Functions_MatVect *
   **********************************/
  
#ifdef SELDON_WITH_REDUCED_TEMPLATE
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline void Mlt(const Matrix<T, Prop0, Storage0, Allocator0>& M,
		  const Vector<T, Storage1, Allocator1>& X,
		  Vector<T, Storage2, Allocator2>& Y)
  {
    MltVector(M, X, Y);
  }
#else
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  inline void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		  const Vector<T1, Storage1, Allocator1>& X,
		  Vector<T2, Storage2, Allocator2>& Y)
  {
    MltVector(M, X, Y);
  }
#endif

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline void Mlt(const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
		  const Vector<T, Storage1, Allocator1>& X,
		  Vector<T, Storage2, Allocator2>& Y)
  {
    throw WrongArgument("Mlt", "Incompatible matrix-vector product");
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline void Mlt(const Matrix<T, Prop0, Storage0, Allocator0>& M,
		  const Vector<complex<T>, Storage1, Allocator1>& X,
		  Vector<complex<T>, Storage2, Allocator2>& Y)
  {
    MltVector(M, X, Y);
  }


  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  void Mlt(const SeldonTranspose& Trans,
	   const Matrix<T, Prop1, Storage1, Allocator1>& M,
	   const Vector<T, Storage2, Allocator2>& X,
	   Vector<T, Storage3, Allocator3>& Y)
  {
    MltVector(Trans, M, X, Y);
  }


  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  inline void Mlt(const SeldonTranspose& Trans,
		  const Matrix<T, Prop1, Storage1, Allocator1>& M,
		  const Vector<complex<T>, Storage2, Allocator2>& X,
		  Vector<complex<T>, Storage3, Allocator3>& Y)
  {
    MltVector(Trans, M, X, Y);
  }
  
  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  inline void Mlt(const SeldonTranspose& Trans,
		  const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
		  const Vector<T, Storage2, Allocator2>& X,
		  Vector<T, Storage3, Allocator3>& Y)
  {
    throw WrongArgument("Mlt", "Incompatible matrix-vector product");
  }

#ifdef SELDON_WITH_REDUCED_TEMPLATE
  template <class T, class Prop1, class Storage1, class Allocator1,
	    class Storage2, class Allocator2,
	    class Storage3, class Allocator3>
  inline void Mlt(const T& alpha,
		  const Matrix<T, Prop1, Storage1, Allocator1>& M,
		  const Vector<T, Storage2, Allocator2>& X,
		  Vector<T, Storage3, Allocator3>& Y)
  {
    MltVector(M, X, Y);
    Mlt(alpha, Y);
  }
#else
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  inline void Mlt(const T0& alpha,
		  const Matrix<T1, Prop1, Storage1, Allocator1>& M,
		  const Vector<T2, Storage2, Allocator2>& X,
		  Vector<T3, Storage3, Allocator3>& Y)
  {
    MltVector(M, X, Y);
    Mlt(alpha, Y);
  }
#endif

  template <class T, class Prop1, class Storage1, class Allocator1,
	    class Storage2, class Allocator2,
	    class Storage3, class Allocator3>
  inline void Mlt(const complex<T>& alpha,
		  const Matrix<T, Prop1, Storage1, Allocator1>& M,
		  const Vector<complex<T>, Storage2, Allocator2>& X,
		  Vector<complex<T>, Storage3, Allocator3>& Y)
  {
    MltVector(M, X, Y);
    Mlt(alpha, Y);
  }
    
  template <class T, class Prop1, class Storage1, class Allocator1,
	    class Storage2, class Allocator2,
	    class Storage3, class Allocator3>
  inline void Mlt(const T& alpha,
		  const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
		  const Vector<T, Storage2, Allocator2>& X,
		  Vector<T, Storage3, Allocator3>& Y)
  {
    throw WrongArgument("Mlt", "Incompatible matrix-vector product");
  }
        
#ifdef SELDON_WITH_REDUCED_TEMPLATE
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha,
		     const Matrix<T, Prop1, Storage1, Allocator1>& M,
		     const Vector<T, Storage2, Allocator2>& X,
		     const T& beta,
		     Vector<T, Storage4, Allocator4>& Y)
  {
    MltAddVector(alpha, M, X, beta, Y);
  }
#else
  template<class T,
	   class T1, class Prop1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3,
	   class T4, class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha,
		     const Matrix<T1, Prop1, Storage1, Allocator1>& M,
		     const Vector<T2, Storage2, Allocator2>& X,
		     const T3& beta,
		     Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAddVector(alpha, M, X, beta, Y);
  }
#endif

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha,
		     const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
		     const Vector<complex<T>, Storage2, Allocator2>& X,
		     const T& beta,
		     Vector<complex<T>, Storage4, Allocator4>& Y)
  {
    MltAddVector(alpha, M, X, beta, Y);
  }
  
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha,
		     const Matrix<T, Prop1, Storage1, Allocator1>& M,
		     const Vector<complex<T>, Storage2, Allocator2>& X,
		     const T& beta,
		     Vector<complex<T>, Storage4, Allocator4>& Y)
  {
    MltAddVector(alpha, M, X, beta, Y);
  }

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha,
		     const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
		     const Vector<T, Storage2, Allocator2>& X,
		     const T& beta,
		     Vector<T, Storage4, Allocator4>& Y)
  {
    throw WrongArgument("MltAdd", "Incompatible matrix-vector product");
  }

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const complex<T>& alpha,
		     const Matrix<T, Prop1, Storage1, Allocator1>& M,
		     const Vector<complex<T>, Storage2, Allocator2>& X,
		     const complex<T>& beta,
		     Vector<complex<T>, Storage4, Allocator4>& Y)
  {
    MltAddVector(alpha, M, X, beta, Y);
  }


  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha, const SeldonTranspose& Trans,
		     const Matrix<T, Prop1, Storage1, Allocator1>& M,
		     const Vector<T, Storage2, Allocator2>& X,
		     const T& beta,
		     Vector<T, Storage4, Allocator4>& Y)
  {
    MltAddVector(alpha, Trans, M, X, beta, Y);
  }


  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha, const SeldonTranspose& Trans,
		     const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
		     const Vector<complex<T>, Storage2, Allocator2>& X,
		     const T& beta,
		     Vector<complex<T>, Storage4, Allocator4>& Y)
  {
    MltAddVector(alpha, Trans, M, X, beta, Y);
  }
  
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha, const SeldonTranspose& Trans,
		     const Matrix<T, Prop1, Storage1, Allocator1>& M,
		     const Vector<complex<T>, Storage2, Allocator2>& X,
		     const T& beta,
		     Vector<complex<T>, Storage4, Allocator4>& Y)
  {
    MltAddVector(alpha, Trans, M, X, beta, Y);
  }

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const T& alpha, const SeldonTranspose& Trans,
		     const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
		     const Vector<T, Storage2, Allocator2>& X,
		     const T& beta,
		     Vector<T, Storage4, Allocator4>& Y)
  {
    throw WrongArgument("MltAdd", "Incompatible matrix-vector product");
  }

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  inline void MltAdd(const complex<T>& alpha, const SeldonTranspose& Trans,
		     const Matrix<T, Prop1, Storage1, Allocator1>& M,
		     const Vector<complex<T>, Storage2, Allocator2>& X,
		     const complex<T>& beta,
		     Vector<complex<T>, Storage4, Allocator4>& Y)
  {

    MltAddVector(alpha, Trans, M, X, beta, Y);
  }
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const Matrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const typename ClassComplexType<T>::Treal& omega,
		  int iter, int type_ssor)
  {
    SorVector(M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const SeldonTranspose& trans,
		  const Matrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const typename ClassComplexType<T>::Treal& omega, int iter, int type_ssor)
  {
    SorVector(trans, M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const Matrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<complex<T>, Storage2, Allocator2>& Y,
		  const Vector<complex<T>, Storage1, Allocator1>& X,
		  const typename ClassComplexType<T>::Treal& omega,
		  int iter, int type_ssor)
  {
    SorVector(M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const SeldonTranspose& trans,
		  const Matrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<complex<T>, Storage2, Allocator2>& Y,
		  const Vector<complex<T>, Storage1, Allocator1>& X,
		  const typename ClassComplexType<T>::Treal& omega, int iter, int type_ssor)
  {
    SorVector(trans, M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const typename ClassComplexType<T>::Treal& omega,
		  int iter, int type_ssor)
  {
    throw WrongArgument("SOR", "incompatible types");
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const SeldonTranspose&,
		  const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const typename ClassComplexType<T>::Treal& omega, int iter, int type_ssor)
  {
    throw WrongArgument("SOR", "incompatible types");
  }


  template <class T, class Prop0, class Allocator0,
            class Allocator1>
  void
  inline Solve(const SeldonUplo& Uplo,
	       const SeldonTranspose& TransA,
	       const SeldonDiag& DiagA,
	       const Matrix<T, Prop0, RowSparse, Allocator0>& A,
	       const Vector<T, VectFull, Allocator1>& X,
	       Vector<T, VectFull, Allocator1>& Y)
  {
    SolveTriangular(Uplo, TransA, DiagA, A, X, Y);
  }

  // SolveLU without pivot
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  inline void SolveLU(const Matrix<T, Prop0, Storage0, Allocator0>& M,
		      Vector<T, Storage1, Allocator1>& Y)
  {
    SolveLuVector(M, Y);
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  inline void SolveLU(const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
		      Vector<T, Storage1, Allocator1>& Y)
  {
    throw WrongArgument("SolveLU", "incompatible types");
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  inline void SolveLU(const Matrix<T, Prop0, Storage0, Allocator0>& M,
		      Vector<complex<T>, Storage1, Allocator1>& Y)
  {
    SolveLuVector(M, Y);
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  inline void SolveLU(const SeldonTranspose& transA,
		      const Matrix<T, Prop0, Storage0, Allocator0>& M,
		      Vector<T, Storage1, Allocator1>& Y)
  {
    SolveLuVector(transA, M, Y);
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  inline void SolveLU(const SeldonTranspose& transA,
		      const Matrix<complex<T>, Prop0, Storage0, Allocator0>& M,
		      Vector<T, Storage1, Allocator1>& Y)
  {
    throw WrongArgument("SolveLU", "incompatible types");
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1>
  inline void SolveLU(const SeldonTranspose& transA,
		      const Matrix<T, Prop0, Storage0, Allocator0>& M,
		      Vector<complex<T>, Storage1, Allocator1>& Y)
  {
    SolveLuVector(transA, M, Y);
  }

  // SolveLU with pivot

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		      const Vector<int>& pivot, Vector<T0>& x)
  {
    SolveLuVector(A, pivot, x);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void SolveLU(const Matrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		      const Vector<int>& pivot, Vector<T0>& x)
  {
    throw WrongArgument("SolveLU", "incompatible types");
  }

  template<class T0, class Storage0, class Allocator0>
  inline void SolveLU(const SeldonTranspose& trans,
		      const Matrix<T0, General, Storage0, Allocator0>& A,
		      const Vector<int>& pivot, Vector<T0>& x)
  {
    SolveLuVector(trans, A, pivot, x);
  }

  template<class T0, class Storage0, class Allocator0>
  inline void SolveLU(const SeldonTranspose& trans,
		      const Matrix<complex<T0>, General, Storage0, Allocator0>& A,
		      const Vector<int>& pivot, Vector<T0>& x)
  {
    throw WrongArgument("SolveLU", "incompatible types");
  }


  template<class T0, class Storage0, class Allocator0>
  inline void SolveLU(const SeldonTranspose& trans,
		      const Matrix<T0, Symmetric, Storage0, Allocator0>& A,
		      const Vector<int>& pivot, Vector<T0>& x)
  {
    SolveLuVector(A, pivot, x);
  }

  template<class T0, class Storage0, class Allocator0>
  inline void SolveLU(const SeldonTranspose& trans,
		      const Matrix<complex<T0>, Symmetric, Storage0, Allocator0>& A,
		      const Vector<int>& pivot, Vector<T0>& x)
  {
    throw WrongArgument("SolveLU", "incompatible types");
  }

  /*************************************
   * Functions in Functions_Matrix.cxx *
   *************************************/

  template<class T, class Prop, class Storage, class Allocator>
  inline void Mlt(const T& alpha, Matrix<T, Prop, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }

  template<class T, class Prop, class Storage, class Allocator>
  inline void Mlt(const T& alpha, Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2,
	    class Prop3, class Storage3, class Allocator3>
  inline void Mlt(const T& alpha,
		  const Matrix<T, Prop1, Storage1, Allocator1>& A,
		  const Matrix<T, Prop2, Storage2, Allocator2>& B,
		  Matrix<T, Prop3, Storage3, Allocator3>& C)
  {
    T zero;
    SetComplexZero(zero);
    C.Zero();
    MltAddMatrix(alpha, A, B, zero, C);
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Mlt(const Matrix<T, Prop0, Storage0, Allocator0>& A,
		  const Matrix<T, Prop1, Storage1, Allocator1>& B,
		  Matrix<T, Prop2, Storage2, Allocator2>& C)
  {
    T one, zero;
    SetComplexZero(zero);
    SetComplexOne(one);
    C.Zero();
    MltAddMatrix(one, A, B, zero, C);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2,
	    class Prop3, class Storage3, class Allocator3>
  inline void MltAdd(const T& alpha,
		     const Matrix<T, Prop1, Storage1, Allocator1>& A,
		     const Matrix<T, Prop2, Storage2, Allocator2>& B,
		     const T& beta,
		     Matrix<T, Prop3, Storage3, Allocator3>& C)
  {
    MltAddMatrix(alpha, A, B, beta, C);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2,
	    class Prop3, class Storage3, class Allocator3>
  inline void MltAdd(const T& alpha, const SeldonTranspose& transA,
		     const Matrix<T, Prop1, Storage1, Allocator1>& A,
		     const SeldonTranspose& transB,
		     const Matrix<T, Prop2, Storage2, Allocator2>& B,
		     const T& beta,
		     Matrix<T, Prop3, Storage3, Allocator3>& C)
  {
    MltAddMatrix(alpha, transA, A, transB, B, beta, C);
  }


  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const Matrix<T, Prop1, Storage1, Allocator1>& A,
		  Matrix<T, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const complex<T>& alpha,
		  const Matrix<T, Prop1, Storage1, Allocator1>& A,
		  Matrix<complex<T>, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }


  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const Matrix<complex<T>, Prop1, Storage1, Allocator1>& A,
		  Matrix<complex<T>, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const Matrix<complex<T>, Prop1, Storage1, Allocator1>& A,
		  Matrix<T, Prop2, Storage2, Allocator2>& B)
  {
    throw WrongArgument("Add", "incompatible types");    
  }


  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Prop2, class Storage2, class Allocator2>
  inline void Copy(const Matrix<T, Prop1, Storage1, Allocator1>& A,
		   Matrix<T, Prop2, Storage2, Allocator2>& B)
  {
    CopyMatrix(A, B);
  }

  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Prop2, class Storage2, class Allocator2>
  inline void Copy(const Matrix<complex<T>, Prop1, Storage1, Allocator1>& A,
		   Matrix<T, Prop2, Storage2, Allocator2>& B)
  {
    throw WrongArgument("Copy", "incompatible types");    
  }

  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Prop2, class Storage2, class Allocator2>
  inline void Copy(const Matrix<T, Prop1, Storage1, Allocator1>& A,
		   Matrix<complex<T>, Prop2, Storage2, Allocator2>& B)
  {
    CopyMatrix(A, B);
  }


  //! returns true if the matrix is symmetric
  template<class T, class Prop, class Storage, class Allocator>
  inline bool IsSymmetricMatrix(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    return false;
  }


  //! returns true if the matrix is symmetric
  template<class T, class Storage, class Allocator>
  inline bool IsSymmetricMatrix(const Matrix<T, Symmetric, Storage, Allocator>& A)
  {
    return true;
  }
  
  
  //! returns true if the matrix is complex
  template<class T, class Prop, class Storage, class Allocator>
  inline bool IsComplexMatrix(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    return false;
  }


  //! returns true if the matrix is complex
  template<class T, class Prop, class Storage, class Allocator>
  inline bool IsComplexMatrix(const Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    return true;
  }
  
}

#define SELDON_FILE_FUNCTIONS_BASE_INLINE_CXX
#endif
