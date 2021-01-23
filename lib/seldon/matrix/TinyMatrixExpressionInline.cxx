#ifndef SELDON_FILE_TINY_MATRIX_EXPRESSION_INLINE_CXX

namespace Seldon
{
  
  //! returns the element i of expression
  template<class T, int m, int n, class E>
  inline const T TinyMatrixExpression<T, m, n, E>::operator()(int i, int j) const
  {
    return static_cast<const E&>(*this)(i, j);
  }
    
  
  //! Constructor u-v with two expressions u and v
  template<class T, int m, int n, class E1, class E2> template<class T1, class T2>
  inline TinyMatrixDifference<T, m, n, E1, E2>::
  TinyMatrixDifference(const TinyMatrixExpression<T1, m, n, E1>& u,
		       const TinyMatrixExpression<T2, m, n, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the difference
  template<class T, int m, int n, class E1, class E2> 
  inline const T TinyMatrixDifference<T, m, n, E1, E2>::operator()(int i, int j) const
  {
    return u_(i, j) - v_(i, j);
  }
  
  
  //! Constructor u+v with two expressions u and v
  template<class T, int m, int n, class E1, class E2> template<class T1, class T2>
  inline TinyMatrixSum<T, m, n, E1, E2>
  ::TinyMatrixSum(const TinyMatrixExpression<T1, m, n, E1>& u,
		  const TinyMatrixExpression<T2, m, n, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the sum  
  template<class T, int m, int n, class E1, class E2>
  inline const T TinyMatrixSum<T, m, n, E1, E2>::operator()(int i, int j) const
  {
    return u_(i, j) + v_(i, j);
  }
  

  //! Constructor u*v with two expressions u and v
  template<class T, int m, int n, class E1, class E2> template<class T1, class T2>
  inline TinyMatrixProduct<T, m, n, E1, E2>
  ::TinyMatrixProduct(const TinyMatrixExpression<T1, m, n, E1>& u,
		      const TinyMatrixExpression<T2, m, n, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the element-wise product 
  template<class T, int m, int n, class E1, class E2>
  inline const T TinyMatrixProduct<T, m, n, E1, E2>::operator()(int i, int j) const
  {
    return u_(i, j) * v_(i, j);
  }


  //! Constructor u / v with two expressions u and v
  template<class T, int m, int n, class E1, class E2> template<class T1, class T2>
  inline TinyMatrixDivision<T, m, n, E1, E2>
  ::TinyMatrixDivision(const TinyMatrixExpression<T1, m, n, E1>& u,
		       const TinyMatrixExpression<T2, m, n, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the element-wise division
  template<class T, int m, int n, class E1, class E2>
  inline const T TinyMatrixDivision<T, m, n, E1, E2>::operator()(int i, int j) const
  {
    return u_(i, j) / v_(i, j);
  }


  //! Constructor alpha * u with a scalar alpha and an expression u
  template<class T, int m, int n, class T0, class E> template<class T1>
  inline TinyMatrixScaled<T, m, n, T0, E>
  ::TinyMatrixScaled(const T0& alpha,
		     const TinyMatrixExpression<T1, m, n, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha*u
  template<class T, int m, int n, class T0, class E>
  inline const T TinyMatrixScaled<T, m, n, T0, E>::operator()(int i, int j) const
  {
    return alpha_*u_(i, j);
  }


  //! Constructor alpha / u with a scalar alpha and an expression u
  template<class T, int m, int n, class T0, class E> template<class T1>
  inline TinyMatrixScalDiv<T, m, n, T0, E>::TinyMatrixScalDiv(const T0& alpha,
							      const TinyMatrixExpression<T1, m, n, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha / u
  template<class T, int m, int n, class T0, class E>
  inline const T TinyMatrixScalDiv<T, m, n, T0, E>::operator()(int i, int j) const
  {
    return alpha_ / u_(i, j);
  }


  //! Constructor alpha + u with a scalar alpha and an expression u
  template<class T, int m, int n, class T0, class E> template<class T1>
  inline TinyMatrixScalSum<T, m, n, T0, E>
  ::TinyMatrixScalSum(const T0& alpha,
		      const TinyMatrixExpression<T1, m, n, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha + u
  template<class T, int m, int n, class T0, class E>
  inline const T TinyMatrixScalSum<T, m, n, T0, E>::operator()(int i, int j) const
  {
    return alpha_ + u_(i, j);
  }


  //! Constructor alpha - u with a scalar alpha and an expression u
  template<class T, int m, int n, class T0, class E> template<class T1>
  inline TinyMatrixScalDiff<T, m, n, T0, E>
  ::TinyMatrixScalDiff(const T0& alpha,
		       const TinyMatrixExpression<T1, m, n, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha - u
  template<class T, int m, int n, class T0, class E>
  inline const T TinyMatrixScalDiff<T, m, n, T0, E>::operator()(int i, int j) const
  {
    return alpha_ - u_(i, j);
  }


  //! Constructor -u with an expression u
  template<class T, int m, int n, class E>
  inline TinyMatrixOpposite<T, m, n, E>::TinyMatrixOpposite(const TinyMatrixExpression<T, m, n, E>& u)
    : u_(u)
  {
  }
  
    
  //! returns the i-th element of -u
  template<class T, int m, int n, class E>
  inline const T TinyMatrixOpposite<T, m, n, E>::operator()(int i, int j) const
  {
    return -u_(i, j);
  }


  //! Constructor transpose(u) with an expression u
  template<class T, int m, int n, class E>
  inline TinyMatrixTranspose<T, m, n, E>::TinyMatrixTranspose(const TinyMatrixExpression<T, n, m, E>& u)
    : u_(u)
  {
  }
  
    
  //! returns the i-th element of -u
  template<class T, int m, int n, class E>
  inline const T TinyMatrixTranspose<T, m, n, E>::operator()(int i, int j) const
  {
    return u_(j, i);
  }
  
  
  /*************
   * Operators *
   *************/
  
  
  //! returns u+v
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixSum<T, m, n, E1, E2> 
  operator +(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v)
  {
    return TinyMatrixSum<T, m, n, E1, E2>(u, v);
  }
  

  //! returns u+v
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixSum<complex<T>, m, n, E1, E2> 
  operator +(const TinyMatrixExpression<complex<T>, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v)
  {
    return TinyMatrixSum<complex<T>, m, n, E1, E2>(u, v);
  }
  

  //! returns u+v
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixSum<complex<T>, m, n, E1, E2> 
  operator +(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<complex<T>, m, n, E2>& v)
  {
    return TinyMatrixSum<complex<T>, m, n, E1, E2>(u, v);
  }
 

  //! returns u-v
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixDifference<T, m, n, E1, E2>
  operator -(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v)
  {
    return TinyMatrixDifference<T, m, n, E1, E2>(u, v);
  }


  //! returns u-v
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixDifference<complex<T>, m, n, E1, E2>
  operator -(const TinyMatrixExpression<complex<T>, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v)
  {
    return TinyMatrixDifference<complex<T>, m, n, E1, E2>(u, v);
  }


  //! returns u-v
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixDifference<complex<T>, m, n, E1, E2>
  operator -(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<complex<T>, m, n, E2>& v)
  {
    return TinyMatrixDifference<complex<T>, m, n, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixProduct<T, m, n, E1, E2> 
  operator *(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v)
  {
    return TinyMatrixProduct<T, m, n, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixProduct<complex<T>, m, n, E1, E2> 
  operator *(const TinyMatrixExpression<complex<T>, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v)
  {
    return TinyMatrixProduct<complex<T>, m, n, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixProduct<complex<T>, m, n, E1, E2> 
  operator *(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<complex<T>, m, n, E2>& v)
  {
    return TinyMatrixProduct<complex<T>, m, n, E1, E2>(u, v);
  }
  

  //! returns u/v element-wise
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixDivision<T, m, n, E1, E2>
  operator /(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v)
  {
    return TinyMatrixDivision<T, m, n, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixDivision<complex<T>, m, n, E1, E2>
  operator /(const TinyMatrixExpression<complex<T>, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v)
  {
    return TinyMatrixDivision<complex<T>, m, n, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixDivision<complex<T>, m, n, E1, E2>
  operator /(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<complex<T>, m, n, E2>& v)
  {
    return TinyMatrixDivision<complex<T>, m, n, E1, E2>(u, v);
  }
  
  
  //! returns alpha*u
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<T, m, n, T, E> operator *(const T& alpha,
						       const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixScaled<T, m, n, T, E>(alpha, u);
  }


  //! returns alpha*u
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<complex<T>, m, n, complex<T>, E>
  operator *(const complex<T>& alpha,
	     const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixScaled<complex<T>, m, n, complex<T>, E>(alpha, u);
  }

  
  //! returns alpha*u
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<complex<T>, m, n, T, E>
  operator *(const T& alpha,
	     const TinyMatrixExpression<complex<T>, m, n, E>& u)
  {
    return TinyMatrixScaled<complex<T>, m, n, T, E>(alpha, u);
  }


  //! returns u*alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<T, m, n, T, E> operator *(const TinyMatrixExpression<T, m, n, E>& u,
						       const T& alpha)
  {
    return TinyMatrixScaled<T, m, n, T, E>(alpha, u);
  }


  //! returns u*alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<complex<T>, m, n, complex<T>, E>
  operator *(const TinyMatrixExpression<T, m, n, E>& u,
	     const complex<T>& alpha)
  {
    return TinyMatrixScaled<complex<T>, m, n, complex<T>, E>(alpha, u);
  }

  
  //! returns u*alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<complex<T>, m, n, T, E>
  operator *(const TinyMatrixExpression<complex<T>, m, n, E>& u,
	     const T& alpha)
  {
    return TinyMatrixScaled<complex<T>, m, n, T, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalDiv<T, m, n, T, E> operator /(const T& alpha,
							const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixScalDiv<T, m, n, T, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalDiv<complex<T>, m, n, complex<T>, E>
  operator /(const complex<T>& alpha,
	     const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixScalDiv<complex<T>, m, n, complex<T>, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalDiv<complex<T>, m, n, T, E>
  operator /(const T& alpha,
	     const TinyMatrixExpression<complex<T>, m, n, E>& u)
  {
    return TinyMatrixScalDiv<complex<T>, m, n, T, E>(alpha, u);
  }

  
  //! returns u / alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<T, m, n, T, E> operator /(const TinyMatrixExpression<T, m, n, E>& u,
						       const T& alpha)
  {
    T one; SetComplexOne(one);
    return TinyMatrixScaled<T, m, n, T, E>(one/alpha, u);
  }


  //! returns u / alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<complex<T>, m, n, complex<T>, E>
  operator /(const TinyMatrixExpression<T, m, n, E>& u,
	     const complex<T>& alpha)
  {
    T one; SetComplexOne(one);
    return TinyMatrixScaled<complex<T>, m, n, complex<T>, E>(one/alpha, u);
  }


  //! returns u / alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScaled<complex<T>, m, n, T, E> 
  operator /(const TinyMatrixExpression<complex<T>, m, n, E>& u,
	     const T& alpha)
  {
    T one; SetComplexOne(one);
    return TinyMatrixScaled<complex<T>, m, n, T, E>(one/alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<T, m, n, T, E> operator +(const T& alpha,
							const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixScalSum<T, m, n, T, E>(alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>
  operator +(const complex<T>& alpha,
	     const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>(alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<complex<T>, m, n, T, E>
  operator +(const T& alpha,
	     const TinyMatrixExpression<complex<T>, m, n, E>& u)
  {
    return TinyMatrixScalSum<complex<T>, m, n, T, E>(alpha, u);
  }

  
  //! returns u + alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<T, m, n, T, E> operator +(const TinyMatrixExpression<T, m, n, E>& u,
							const T& alpha)
  {
    return TinyMatrixScalSum<T, m, n, T, E>(alpha, u);
  }


  //! returns u + alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>
  operator +(const TinyMatrixExpression<T, m, n, E>& u,
	     const complex<T>& alpha)
  {
    return TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>(alpha, u);
  }


  //! returns u + alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<complex<T>, m, n, T, E>
  operator +(const TinyMatrixExpression<complex<T>, m, n, E>& u,
	     const T& alpha)
  {
    return TinyMatrixScalSum<complex<T>, m, n, T, E>(alpha, u);
  }

    
  //! returns u - alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<T, m, n, T, E> operator -(const TinyMatrixExpression<T, m, n, E>& u,
							const T& alpha)
  {
    return TinyMatrixScalSum<T, m, n, T, E>(-alpha, u);
  }


  //! returns u - alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>
  operator -(const TinyMatrixExpression<T, m, n, E>& u,
	     const complex<T>& alpha)
  {
    return TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>(-alpha, u);
  }


  //! returns u - alpha
  template<class T, int m, int n, class E>
  inline const TinyMatrixScalSum<complex<T>, m, n, T, E>
  operator -(const TinyMatrixExpression<complex<T>, m, n, E>& u,
	     const T& alpha)
  {
    return TinyMatrixScalSum<complex<T>, m, n, T, E>(-alpha, u);
  }

    
  //! returns alpha - u
  template<class T, int m, int n, class E> 
  inline const TinyMatrixScalDiff<T, m, n, T, E> operator -(const T& alpha,
						      const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixScalDiff<T, m, n, T, E>(alpha, u);
  }


  //! returns alpha - u
  template<class T, int m, int n, class E> 
  inline const TinyMatrixScalDiff<complex<T>, m, n, complex<T>, E>
  operator -(const complex<T>& alpha,
	     const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixScalDiff<complex<T>, m, n, complex<T>, E>(alpha, u);
  }


  //! returns alpha - u
  template<class T, int m, int n, class E> 
  inline const TinyMatrixScalDiff<complex<T>, m, n, T, E>
  operator -(const T& alpha,
	     const TinyMatrixExpression<complex<T>, m, n, E>& u)
  {
    return TinyMatrixScalDiff<complex<T>, m, n, T, E>(alpha, u);
  }
  
  
  //! returns -u
  template<class T, int m, int n, class E>
  inline const TinyMatrixOpposite<T, m, n, E> operator-(const TinyMatrixExpression<T, m, n, E>& u)
  {
    return TinyMatrixOpposite<T, m, n, E>(u);
  }


  //! returns transpose(u)
  template<class T, int m, int n, class E>
  inline const TinyMatrixTranspose<T, m, n, E> transpose(const TinyMatrixExpression<T, n, m, E>& u)
  {
    return TinyMatrixTranspose<T, m, n, E>(u);
  }

}

#define SELDON_FILE_TINY_MATRIX_EXPRESSION_INLINE_CXX
#endif
