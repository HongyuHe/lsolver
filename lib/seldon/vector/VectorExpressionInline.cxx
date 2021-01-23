#ifndef SELDON_FILE_VECTOR_EXPRESSION_INLINE_CXX

namespace Seldon
{
  
  //! returns the size of the associated vector
  template<class T, class E>
  inline long VectorExpression<T, E>::GetSize() const
  {
    return static_cast<const E&>(*this).GetSize();
  }

  
  //! returns the element i of expression
  template<class T, class E>
  inline const T VectorExpression<T, E>::operator()(long i) const
  {
    return static_cast<const E&>(*this)(i);
  }
    
  
  //! Constructor u-v with two expressions u and v
  template<class T, class E1, class E2> template<class T1, class T2>
  inline VectorDifference<T, E1, E2>::
  VectorDifference(const VectorExpression<T1, E1>& u,
		   const VectorExpression<T2, E2>& v)
    : u_(u), v_(v)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (u_.GetSize() != v_.GetSize())
      throw WrongDim("VectorDifference",
		     string("Cannot subtract u and v because the sizes ")
		     +to_str(u.GetSize()) + " and " + to_str(v.GetSize())
		     + " are different");
#endif
  }
  
  
  //! returns the size of the associated vectors
  template<class T, class E1, class E2>
  inline long VectorDifference<T, E1, E2>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of the difference
  template<class T, class E1, class E2>
  inline const T VectorDifference<T, E1, E2>::operator()(long i) const
  {
    return u_(i) - v_(i);
  }
  
  
  //! Constructor u+v with two expressions u and v
  template<class T, class E1, class E2> template<class T1, class T2>
  inline VectorSum<T, E1, E2>::VectorSum(const VectorExpression<T1, E1>& u,
					 const VectorExpression<T2, E2>& v)
    : u_(u), v_(v)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (long(u_.GetSize()) != long(v_.GetSize()))
      throw WrongDim("VectorDifference",
                     string("Cannot Add u and v because the sizes ")
                     +to_str(u.GetSize()) + " and " + to_str(v.GetSize())
		     + " are different");
#endif
  }
  
  
  //! returns the size of the associated vectors
  template<class T, class E1, class E2>
  inline long VectorSum<T, E1, E2>::GetSize() const
  {
    return u_.GetSize();
  }
  

  //! returns the i-th element of the sum  
  template<class T, class E1, class E2>
  inline const T VectorSum<T, E1, E2>::operator()(long i) const
  {
    return u_(i) + v_(i);
  }
  

  //! Constructor u*v with two expressions u and v
  template<class T, class E1, class E2> template<class T1, class T2>
  inline VectorProduct<T, E1, E2>
  ::VectorProduct(const VectorExpression<T1, E1>& u,
		  const VectorExpression<T2, E2>& v)
    : u_(u), v_(v)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (u_.GetSize() != v_.GetSize())
      throw WrongDim("VectorDifference",
                     string("Cannot multiply u and v because the sizes ")
                     +to_str(u.GetSize()) + " and " + to_str(v.GetSize())
		     + " are different");
#endif
  }
  
  
  //! returns the size of the associated vectors
  template<class T, class E1, class E2>
  inline long VectorProduct<T, E1, E2>::GetSize() const
  {
    return u_.GetSize();
  }
  

  //! returns the i-th element of the element-wise product 
  template<class T, class E1, class E2>
  inline const T VectorProduct<T, E1, E2>::operator()(long i) const
  {
    return u_(i) * v_(i);
  }


  //! Constructor u / v with two expressions u and v
  template<class T, class E1, class E2> template<class T1, class T2>
  inline VectorDivision<T, E1, E2>
  ::VectorDivision(const VectorExpression<T1, E1>& u,
		   const VectorExpression<T2, E2>& v)
    : u_(u), v_(v)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (u_.GetSize() != v_.GetSize())
      throw WrongDim("VectorDifference",
                     string("Cannot divide u and v because the sizes ")
                     +to_str(u.GetSize()) + " and " + to_str(v.GetSize())
		     + " are different");
#endif
  }
  
  
  //! returns the size of the associated vectors
  template<class T, class E1, class E2>
  inline long VectorDivision<T, E1, E2>::GetSize() const
  {
    return u_.GetSize();
  }
  

  //! returns the i-th element of the element-wise division
  template<class T, class E1, class E2>
  inline const T VectorDivision<T, E1, E2>::operator()(long i) const
  {
    return u_(i) / v_(i);
  }


  //! Constructor alpha * u with a scalar alpha and an expression u
  template<class T, class T0, class E> template<class T1>
  inline VectorScaled<T, T0, E>::VectorScaled(const T0& alpha,
					      const VectorExpression<T1, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the size of the associated vector
  template<class T, class T0, class E>
  inline long VectorScaled<T, T0, E>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of alpha*u
  template<class T, class T0, class E>
  inline const T VectorScaled<T, T0, E>::operator()(long i) const
  {
    return alpha_*u_(i);
  }


  //! Constructor alpha / u with a scalar alpha and an expression u
  template<class T, class T0, class E> template<class T1>
  inline VectorScalDiv<T, T0, E>::VectorScalDiv(const T0& alpha,
						const VectorExpression<T1, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the size of the associated vector
  template<class T, class T0, class E>
  inline long VectorScalDiv<T, T0, E>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of alpha / u
  template<class T, class T0, class E>
  inline const T VectorScalDiv<T, T0, E>::operator()(long i) const
  {
    return alpha_ / u_(i);
  }


  //! Constructor alpha + u with a scalar alpha and an expression u
  template<class T, class T0, class E> template<class T1>
  inline VectorScalSum<T, T0, E>::VectorScalSum(const T0& alpha,
						const VectorExpression<T1, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the size of the associated vector
  template<class T, class T0, class E>
  inline long VectorScalSum<T, T0, E>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of alpha + u
  template<class T, class T0, class E>
  inline const T VectorScalSum<T, T0, E>::operator()(long i) const
  {
    return alpha_ + u_(i);
  }


  //! Constructor alpha - u with a scalar alpha and an expression u
  template<class T, class T0, class E> template<class T1>
  inline VectorScalDiff<T, T0, E>::VectorScalDiff(const T0& alpha,
						  const VectorExpression<T1, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the size of the associated vector
  template<class T, class T0, class E>
  inline long VectorScalDiff<T, T0, E>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of alpha - u
  template<class T, class T0, class E>
  inline const T VectorScalDiff<T, T0, E>::operator()(long i) const
  {
    return alpha_ - u_(i);
  }


  //! Constructor -u with an expression u
  template<class T, class E>
  inline VectorOpposite<T, E>::VectorOpposite(const VectorExpression<T, E>& u)
    : u_(u)
  {
  }
  
    
  //! returns the size of the associated vector
  template<class T, class E>
  inline long VectorOpposite<T, E>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of -u
  template<class T, class E>
  inline const T VectorOpposite<T, E>::operator()(long i) const
  {
    return -u_(i);
  }
  
  
  /*************
   * Operators *
   *************/
  
  
  //! returns u+v
  template<class T, class E1, class E2>
  inline const VectorSum<T, E1, E2> 
  operator +(const VectorExpression<T, E1>& u,
	     const VectorExpression<T, E2>& v)
  {
    return VectorSum<T, E1, E2>(u, v);
  }


  //! returns u+v
  template<class T, class E1, class E2>
  inline const VectorSum<complex<T>, E1, E2> 
  operator +(const VectorExpression<complex<T>, E1>& u,
	     const VectorExpression<T, E2>& v)
  {
    return VectorSum<complex<T>, E1, E2>(u, v);
  }


  //! returns u+v
  template<class T, class E1, class E2>
  inline const VectorSum<complex<T>, E1, E2> 
  operator +(const VectorExpression<T, E1>& u,
	     const VectorExpression<complex<T>, E2>& v)
  {
    return VectorSum<complex<T>, E1, E2>(u, v);
  }
  

  //! returns u-v
  template<class T, class E1, class E2>
  inline const VectorDifference<T, E1, E2>
  operator -(const VectorExpression<T, E1>& u,
	     const VectorExpression<T, E2>& v)
  {
    return VectorDifference<T, E1, E2>(u, v);
  }


  //! returns u-v
  template<class T, class E1, class E2>
  inline const VectorDifference<complex<T>, E1, E2>
  operator -(const VectorExpression<complex<T>, E1>& u,
	     const VectorExpression<T, E2>& v)
  {
    return VectorDifference<complex<T>, E1, E2>(u, v);
  }


  //! returns u-v
  template<class T, class E1, class E2>
  inline const VectorDifference<complex<T>, E1, E2>
  operator -(const VectorExpression<T, E1>& u,
	     const VectorExpression<complex<T>, E2>& v)
  {
    return VectorDifference<complex<T>, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, class E1, class E2>
  inline const VectorProduct<T, E1, E2> 
  operator *(const VectorExpression<T, E1>& u,
	     const VectorExpression<T, E2>& v)
  {
    return VectorProduct<T, E1, E2>(u, v);
  }
  

  //! returns u*v element-wise
  template<class T, class E1, class E2>
  inline const VectorProduct<complex<T>, E1, E2> 
  operator *(const VectorExpression<complex<T>, E1>& u,
	     const VectorExpression<T, E2>& v)
  {
    return VectorProduct<complex<T>, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, class E1, class E2>
  inline const VectorProduct<complex<T>, E1, E2> 
  operator *(const VectorExpression<T, E1>& u,
	     const VectorExpression<complex<T>, E2>& v)
  {
    return VectorProduct<complex<T>, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, class E1, class E2>
  inline const VectorDivision<T, E1, E2>
  operator /(const VectorExpression<T, E1>& u,
	     const VectorExpression<T, E2>& v)
  {
    return VectorDivision<T, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, class E1, class E2>
  inline const VectorDivision<complex<T>, E1, E2>
  operator /(const VectorExpression<complex<T>, E1>& u,
	     const VectorExpression<T, E2>& v)
  {
    return VectorDivision<complex<T>, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, class E1, class E2>
  inline const VectorDivision<complex<T>, E1, E2>
  operator /(const VectorExpression<T, E1>& u,
	     const VectorExpression<complex<T>, E2>& v)
  {
    return VectorDivision<complex<T>, E1, E2>(u, v);
  }
  
  
  //! returns alpha*u
  template<class T, class E>
  inline const VectorScaled<T, T, E>
  operator *(const T& alpha,
	     const VectorExpression<T, E>& u)
  {
    return VectorScaled<T, T, E>(alpha, u);
  }


  //! returns alpha*u
  template<class T, class E>
  inline const VectorScaled<complex<T>, complex<T>, E>
  operator *(const complex<T>& alpha,
	     const VectorExpression<T, E>& u)
  {
    return VectorScaled<complex<T>, complex<T>, E>(alpha, u);
  }


  //! returns alpha*u
  template<class T, class E>
  inline const VectorScaled<complex<T>, T, E>
  operator *(const T& alpha,
	     const VectorExpression<complex<T>, E>& u)
  {
    return VectorScaled<complex<T>, T, E>(alpha, u);
  }

  
  //! returns u*alpha
  template<class T, class E>
  inline const VectorScaled<T, T, E> operator *(const VectorExpression<T, E>& u,
						const T& alpha)
  {
    return VectorScaled<T, T, E>(alpha, u);
  }


  //! returns u*alpha
  template<class T, class E>
  inline const VectorScaled<complex<T>, complex<T>, E>
  operator *(const VectorExpression<T, E>& u,
	     const complex<T>& alpha)
  {
    return VectorScaled<complex<T>, complex<T>, E>(alpha, u);
  }


  //! returns u*alpha
  template<class T, class E>
  inline const VectorScaled<complex<T>, T, E>
  operator *(const VectorExpression<complex<T>, E>& u,
	     const T& alpha)
  {
    return VectorScaled<complex<T>, T, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, class E>
  inline const VectorScalDiv<T, T, E> operator /(const T& alpha,
						 const VectorExpression<T, E>& u)
  {
    return VectorScalDiv<T, T, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, class E>
  inline const VectorScalDiv<complex<T>, complex<T>, E>
  operator /(const complex<T>& alpha,
	     const VectorExpression<T, E>& u)
  {
    return VectorScalDiv<complex<T>, complex<T>, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, class E>
  inline const VectorScalDiv<complex<T>, T, E>
  operator /(const T& alpha,
	     const VectorExpression<complex<T>, E>& u)
  {
    return VectorScalDiv<complex<T>, T, E>(alpha, u);
  }

  
  //! returns u / alpha
  template<class T, class E>
  inline const VectorScaled<T, T, E> operator /(const VectorExpression<T, E>& u,
						const T& alpha)
  {
    T one; SetComplexOne(one);
    return VectorScaled<T, T, E>(one/alpha, u);
  }


  //! returns u / alpha
  template<class T, class E>
  inline const VectorScaled<complex<T>, complex<T>, E>
  operator /(const VectorExpression<T, E>& u,
	     const complex<T>& alpha)
  {
    T one; SetComplexOne(one);
    return VectorScaled<complex<T>, complex<T>, E>(one/alpha, u);
  }


  //! returns u / alpha
  template<class T, class E>
  inline const VectorScaled<complex<T>, T, E>
  operator /(const VectorExpression<complex<T>, E>& u,
	     const T& alpha)
  {
    T one; SetComplexOne(one);
    return VectorScaled<complex<T>, T, E>(one/alpha, u);
  }


  //! returns alpha + u
  template<class T, class E>
  inline const VectorScalSum<T, T, E> operator +(const T& alpha,
						 const VectorExpression<T, E>& u)
  {
    return VectorScalSum<T, T, E>(alpha, u);
  }


  //! returns alpha + u
  template<class T, class E>
  inline const VectorScalSum<complex<T>, complex<T>, E>
  operator +(const complex<T>& alpha,
	     const VectorExpression<T, E>& u)
  {
    return VectorScalSum<complex<T>, complex<T>, E>(alpha, u);
  }


  //! returns alpha + u
  template<class T, class E>
  inline const VectorScalSum<complex<T>, T, E>
  operator +(const T& alpha,
	     const VectorExpression<complex<T>, E>& u)
  {
    return VectorScalSum<complex<T>, T, E>(alpha, u);
  }

  
  //! returns u + alpha
  template<class T, class E>
  inline const VectorScalSum<T, T, E> operator +(const VectorExpression<T, E>& u,
						 const T& alpha)
  {
    return VectorScalSum<T, T, E>(alpha, u);
  }


  //! returns u + alpha
  template<class T, class E>
  inline const VectorScalSum<complex<T>, complex<T>, E>
  operator +(const VectorExpression<T, E>& u,
	     const complex<T>& alpha)
  {
    return VectorScalSum<complex<T>, complex<T>, E>(alpha, u);
  }


  //! returns u + alpha
  template<class T, class E>
  inline const VectorScalSum<complex<T>, T, E>
  operator +(const VectorExpression<complex<T>, E>& u,
	     const T& alpha)
  {
    return VectorScalSum<complex<T>, T, E>(alpha, u);
  }

    
  //! returns u - alpha 
  template<class T, class E>
  inline const VectorScalSum<T, T, E> operator -(const VectorExpression<T, E>& u,
					      const T& alpha)
  {
    return VectorScalSum<T, T, E>(-alpha, u);
  }

    
  //! returns u - alpha 
  template<class T, class E>
  inline const VectorScalSum<complex<T>, complex<T>, E>
  operator -(const VectorExpression<T, E>& u,
	     const complex<T>& alpha)
  {
    return VectorScalSum<complex<T>, complex<T>, E>(-alpha, u);
  }

    
  //! returns u - alpha 
  template<class T, class E>
  inline const VectorScalSum<complex<T>, T, E>
  operator -(const VectorExpression<complex<T>, E>& u,
	     const T& alpha)
  {
    return VectorScalSum<complex<T>, T, E>(-alpha, u);
  }

    
  //! returns alpha - u
  template<class T, class E>
  inline const VectorScalDiff<T, T, E> operator -(const T& alpha,
						  const VectorExpression<T, E>& u)
  {
    return VectorScalDiff<T, T, E>(alpha, u);
  }


  //! returns alpha - u
  template<class T, class E>
  inline const VectorScalDiff<complex<T>, complex<T>, E>
  operator -(const complex<T>& alpha,
	     const VectorExpression<T, E>& u)
  {
    return VectorScalDiff<complex<T>, complex<T>, E>(alpha, u);
  }


  //! returns alpha - u
  template<class T, class E>
  inline const VectorScalDiff<complex<T>, T, E> 
  operator -(const T& alpha,
	     const VectorExpression<complex<T>, E>& u)
  {
    return VectorScalDiff<complex<T>, T, E>(alpha, u);
  }
  
  
  //! returns -u
  template<class T, class E>
  inline const VectorOpposite<T, E> operator-(const VectorExpression<T, E>& u)
  {
    return VectorOpposite<T, E>(u);
  }
  
}

#define SELDON_FILE_VECTOR_EXPRESSION_INLINE_CXX
#endif
