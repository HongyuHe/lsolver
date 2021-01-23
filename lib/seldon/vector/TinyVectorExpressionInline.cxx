#ifndef SELDON_FILE_TINY_VECTOR_EXPRESSION_INLINE_CXX

namespace Seldon
{
  
  //! returns the element i of expression
  template<class T, int m, class E>
  inline const T TinyVectorExpression<T, m, E>::operator()(int i) const
  {
    return static_cast<const E&>(*this)(i);
  }
    
  
  //! Constructor u-v with two expressions u and v
  template<class T, int m, class E1, class E2> template<class T1, class T2>
  inline TinyVectorDifference<T, m, E1, E2>::
  TinyVectorDifference(const TinyVectorExpression<T1, m, E1>& u,
		       const TinyVectorExpression<T2, m, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the difference
  template<class T, int m, class E1, class E2> 
  inline const T TinyVectorDifference<T, m, E1, E2>::operator()(int i) const
  {
    return u_(i) - v_(i);
  }
  
  
  //! Constructor u+v with two expressions u and v
  template<class T, int m, class E1, class E2> template<class T1, class T2>
  inline TinyVectorSum<T, m, E1, E2>
  ::TinyVectorSum(const TinyVectorExpression<T1, m, E1>& u,
		  const TinyVectorExpression<T2, m, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the sum  
  template<class T, int m, class E1, class E2>
  inline const T TinyVectorSum<T, m, E1, E2>::operator()(int i) const
  {
    return u_(i) + v_(i);
  }
  

  //! Constructor u*v with two expressions u and v
  template<class T, int m, class E1, class E2> template<class T1, class T2>
  inline TinyVectorProduct<T, m, E1, E2>
  ::TinyVectorProduct(const TinyVectorExpression<T1, m, E1>& u,
		      const TinyVectorExpression<T2, m, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the element-wise product 
  template<class T, int m, class E1, class E2>
  inline const T TinyVectorProduct<T, m, E1, E2>::operator()(int i) const
  {
    return u_(i) * v_(i);
  }


  //! Constructor u / v with two expressions u and v
  template<class T, int m, class E1, class E2> template<class T1, class T2>
  inline TinyVectorDivision<T, m, E1, E2>
  ::TinyVectorDivision(const TinyVectorExpression<T1, m, E1>& u,
		       const TinyVectorExpression<T2, m, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the element-wise division
  template<class T, int m, class E1, class E2>
  inline const T TinyVectorDivision<T, m, E1, E2>::operator()(int i) const
  {
    return u_(i) / v_(i);
  }


  //! Constructor alpha * u with a scalar alpha and an expression u
  template<class T, int m, class T0, class E> template<class T1>
  inline TinyVectorScaled<T, m, T0, E>::TinyVectorScaled(const T0& alpha,
							 const TinyVectorExpression<T1, m, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha*u
  template<class T, int m, class T0, class E>
  inline const T TinyVectorScaled<T, m, T0, E>::operator()(int i) const
  {
    return alpha_*u_(i);
  }


  //! Constructor alpha / u with a scalar alpha and an expression u
  template<class T, int m, class T0, class E> template<class T1>
  inline TinyVectorScalDiv<T, m, T0, E>::TinyVectorScalDiv(const T0& alpha,
							   const TinyVectorExpression<T1, m, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha / u
  template<class T, int m, class T0, class E>
  inline const T TinyVectorScalDiv<T, m, T0, E>::operator()(int i) const
  {
    return alpha_ / u_(i);
  }


  //! Constructor alpha + u with a scalar alpha and an expression u
  template<class T, int m, class T0, class E> template<class T1>
  inline TinyVectorScalSum<T, m, T0, E>::TinyVectorScalSum(const T0& alpha,
							   const TinyVectorExpression<T1, m, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha + u
  template<class T, int m, class T0, class E>
  inline const T TinyVectorScalSum<T, m, T0, E>::operator()(int i) const
  {
    return alpha_ + u_(i);
  }


  //! Constructor alpha - u with a scalar alpha and an expression u
  template<class T, int m, class T0, class E> template<class T1>
  inline TinyVectorScalDiff<T, m, T0, E>::TinyVectorScalDiff(const T0& alpha,
							     const TinyVectorExpression<T1, m, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha - u
  template<class T, int m, class T0, class E>
  inline const T TinyVectorScalDiff<T, m, T0, E>::operator()(int i) const
  {
    return alpha_ - u_(i);
  }


  //! Constructor -u with an expression u
  template<class T, int m, class E>
  inline TinyVectorOpposite<T, m, E>::TinyVectorOpposite(const TinyVectorExpression<T, m, E>& u)
    : u_(u)
  {
  }
  
    
  //! returns the i-th element of -u
  template<class T, int m, class E>
  inline const T TinyVectorOpposite<T, m, E>::operator()(int i) const
  {
    return -u_(i);
  }
  
  
  /*************
   * Operators *
   *************/
  
  
  //! returns u+v
  template<class T, int m, class E1, class E2>
  inline const TinyVectorSum<T, m, E1, E2> 
  operator +(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorSum<T, m, E1, E2>(u, v);
  }
  

  //! returns u+v
  template<class T, int m, class E1, class E2>
  inline const TinyVectorSum<complex<T>, m, E1, E2> 
  operator +(const TinyVectorExpression<complex<T>, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorSum<complex<T>, m, E1, E2>(u, v);
  }
  

  //! returns u+v
  template<class T, int m, class E1, class E2>
  inline const TinyVectorSum<complex<T>, m, E1, E2> 
  operator +(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<complex<T>, m, E2>& v)
  {
    return TinyVectorSum<complex<T>, m, E1, E2>(u, v);
  }
 

  //! returns u-v
  template<class T, int m, class E1, class E2>
  inline const TinyVectorDifference<T, m, E1, E2>
  operator -(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorDifference<T, m, E1, E2>(u, v);
  }


  //! returns u-v
  template<class T, int m, class E1, class E2>
  inline const TinyVectorDifference<complex<T>, m, E1, E2>
  operator -(const TinyVectorExpression<complex<T>, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorDifference<complex<T>, m, E1, E2>(u, v);
  }


  //! returns u-v
  template<class T, int m, class E1, class E2>
  inline const TinyVectorDifference<complex<T>, m, E1, E2>
  operator -(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<complex<T>, m, E2>& v)
  {
    return TinyVectorDifference<complex<T>, m, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, class E1, class E2>
  inline const TinyVectorProduct<T, m, E1, E2> 
  operator *(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorProduct<T, m, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, class E1, class E2>
  inline const TinyVectorProduct<complex<T>, m, E1, E2> 
  operator *(const TinyVectorExpression<complex<T>, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorProduct<complex<T>, m, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, class E1, class E2>
  inline const TinyVectorProduct<complex<T>, m, E1, E2> 
  operator *(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<complex<T>, m, E2>& v)
  {
    return TinyVectorProduct<complex<T>, m, E1, E2>(u, v);
  }
  

  //! returns u/v element-wise
  template<class T, int m, class E1, class E2>
  inline const TinyVectorDivision<T, m, E1, E2>
  operator /(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorDivision<T, m, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, int m, class E1, class E2>
  inline const TinyVectorDivision<complex<T>, m, E1, E2>
  operator /(const TinyVectorExpression<complex<T>, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorDivision<complex<T>, m, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, int m, class E1, class E2>
  inline const TinyVectorDivision<complex<T>, m, E1, E2>
  operator /(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<complex<T>, m, E2>& v)
  {
    return TinyVectorDivision<complex<T>, m, E1, E2>(u, v);
  }
  
  
  //! returns alpha*u
  template<class T, int m, class E>
  inline const TinyVectorScaled<T, m, T, E> operator *(const T& alpha,
						       const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorScaled<T, m, T, E>(alpha, u);
  }


  //! returns alpha*u
  template<class T, int m, class E>
  inline const TinyVectorScaled<complex<T>, m, complex<T>, E>
  operator *(const complex<T>& alpha,
	     const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorScaled<complex<T>, m, complex<T>, E>(alpha, u);
  }

  
  //! returns alpha*u
  template<class T, int m, class E>
  inline const TinyVectorScaled<complex<T>, m, T, E>
  operator *(const T& alpha,
	     const TinyVectorExpression<complex<T>, m, E>& u)
  {
    return TinyVectorScaled<complex<T>, m, T, E>(alpha, u);
  }


  //! returns u*alpha
  template<class T, int m, class E>
  inline const TinyVectorScaled<T, m, T, E> operator *(const TinyVectorExpression<T, m, E>& u,
						       const T& alpha)
  {
    return TinyVectorScaled<T, m, T, E>(alpha, u);
  }


  //! returns u*alpha
  template<class T, int m, class E>
  inline const TinyVectorScaled<complex<T>, m, complex<T>, E>
  operator *(const TinyVectorExpression<T, m, E>& u,
	     const complex<T>& alpha)
  {
    return TinyVectorScaled<complex<T>, m, complex<T>, E>(alpha, u);
  }

  
  //! returns u*alpha
  template<class T, int m, class E>
  inline const TinyVectorScaled<complex<T>, m, T, E>
  operator *(const TinyVectorExpression<complex<T>, m, E>& u,
	     const T& alpha)
  {
    return TinyVectorScaled<complex<T>, m, T, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, class E>
  inline const TinyVectorScalDiv<T, m, T, E> operator /(const T& alpha,
							const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorScalDiv<T, m, T, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, class E>
  inline const TinyVectorScalDiv<complex<T>, m, complex<T>, E>
  operator /(const complex<T>& alpha,
	     const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorScalDiv<complex<T>, m, complex<T>, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, class E>
  inline const TinyVectorScalDiv<complex<T>, m, T, E>
  operator /(const T& alpha,
	     const TinyVectorExpression<complex<T>, m, E>& u)
  {
    return TinyVectorScalDiv<complex<T>, m, T, E>(alpha, u);
  }

  
  //! returns u / alpha
  template<class T, int m, class E>
  inline const TinyVectorScaled<T, m, T, E> operator /(const TinyVectorExpression<T, m, E>& u,
						       const T& alpha)
  {
    T one; SetComplexOne(one);
    return TinyVectorScaled<T, m, T, E>(one/alpha, u);
  }


  //! returns u / alpha
  template<class T, int m, class E>
  inline const TinyVectorScaled<complex<T>, m, complex<T>, E>
  operator /(const TinyVectorExpression<T, m, E>& u,
	     const complex<T>& alpha)
  {
    T one; SetComplexOne(one);
    return TinyVectorScaled<complex<T>, m, complex<T>, E>(one/alpha, u);
  }


  //! returns u / alpha
  template<class T, int m, class E>
  inline const TinyVectorScaled<complex<T>, m, T, E> 
  operator /(const TinyVectorExpression<complex<T>, m, E>& u,
	     const T& alpha)
  {
    T one; SetComplexOne(one);
    return TinyVectorScaled<complex<T>, m, T, E>(one/alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, class E>
  inline const TinyVectorScalSum<T, m, T, E> operator +(const T& alpha,
							const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorScalSum<T, m, T, E>(alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, class E>
  inline const TinyVectorScalSum<complex<T>, m, complex<T>, E>
  operator +(const complex<T>& alpha,
	     const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorScalSum<complex<T>, m, complex<T>, E>(alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, class E>
  inline const TinyVectorScalSum<complex<T>, m, T, E>
  operator +(const T& alpha,
	     const TinyVectorExpression<complex<T>, m, E>& u)
  {
    return TinyVectorScalSum<complex<T>, m, T, E>(alpha, u);
  }

  
  //! returns u + alpha
  template<class T, int m, class E>
  inline const TinyVectorScalSum<T, m, T, E> operator +(const TinyVectorExpression<T, m, E>& u,
							const T& alpha)
  {
    return TinyVectorScalSum<T, m, T, E>(alpha, u);
  }


  //! returns u + alpha
  template<class T, int m, class E>
  inline const TinyVectorScalSum<complex<T>, m, complex<T>, E>
  operator +(const TinyVectorExpression<T, m, E>& u,
	     const complex<T>& alpha)
  {
    return TinyVectorScalSum<complex<T>, m, complex<T>, E>(alpha, u);
  }


  //! returns u + alpha
  template<class T, int m, class E>
  inline const TinyVectorScalSum<complex<T>, m, T, E>
  operator +(const TinyVectorExpression<complex<T>, m, E>& u,
	     const T& alpha)
  {
    return TinyVectorScalSum<complex<T>, m, T, E>(alpha, u);
  }

    
  //! returns u - alpha
  template<class T, int m, class E>
  inline const TinyVectorScalSum<T, m, T, E> operator -(const TinyVectorExpression<T, m, E>& u,
							const T& alpha)
  {
    return TinyVectorScalSum<T, m, T, E>(-alpha, u);
  }


  //! returns u - alpha
  template<class T, int m, class E>
  inline const TinyVectorScalSum<complex<T>, m, complex<T>, E>
  operator -(const TinyVectorExpression<T, m, E>& u,
	     const complex<T>& alpha)
  {
    return TinyVectorScalSum<complex<T>, m, complex<T>, E>(-alpha, u);
  }


  //! returns u - alpha
  template<class T, int m, class E>
  inline const TinyVectorScalSum<complex<T>, m, T, E>
  operator -(const TinyVectorExpression<complex<T>, m, E>& u,
	     const T& alpha)
  {
    return TinyVectorScalSum<complex<T>, m, T, E>(-alpha, u);
  }

    
  //! returns alpha - u
  template<class T, int m, class E> 
  inline const TinyVectorScalDiff<T, m, T, E> operator -(const T& alpha,
						      const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorScalDiff<T, m, T, E>(alpha, u);
  }


  //! returns alpha - u
  template<class T, int m, class E> 
  inline const TinyVectorScalDiff<complex<T>, m, complex<T>, E>
  operator -(const complex<T>& alpha,
	     const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorScalDiff<complex<T>, m, complex<T>, E>(alpha, u);
  }


  //! returns alpha - u
  template<class T, int m, class E> 
  inline const TinyVectorScalDiff<complex<T>, m, T, E>
  operator -(const T& alpha,
	     const TinyVectorExpression<complex<T>, m, E>& u)
  {
    return TinyVectorScalDiff<complex<T>, m, T, E>(alpha, u);
  }
  
  
  //! returns -u
  template<class T, int m, class E>
  inline const TinyVectorOpposite<T, m, E> operator-(const TinyVectorExpression<T, m, E>& u)
  {
    return TinyVectorOpposite<T, m, E>(u);
  }
  
}

#define SELDON_FILE_TINY_VECTOR_EXPRESSION_INLINE_CXX
#endif
