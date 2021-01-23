#ifndef SELDON_FILE_TINY_VECTOR_EXPRESSION_HXX

namespace Seldon
{
  
  //! Expression between vectors
  template<class T, int m, class E>
  class TinyVectorExpression
  {
  public:
    const T operator()(int) const;
    
    inline operator E&() { return static_cast<E&>(*this); }
    inline operator E const&() const { return static_cast<const E&>(*this); }
    
  };


  //! Difference between two expressions
  template<class T, int m, class E1, class E2>
  class TinyVectorDifference : public TinyVectorExpression<T, m, TinyVectorDifference<T, m, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyVectorDifference(const TinyVectorExpression<T1, m, E1>&, const TinyVectorExpression<T2, m, E2>&);
    
    const T operator()(int) const;
    
  };


  //! Sum between two expressions
  template<class T, int m, class E1, class E2>
  class TinyVectorSum : public TinyVectorExpression<T, m, TinyVectorSum<T, m, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyVectorSum(const TinyVectorExpression<T1, m, E1>&, const TinyVectorExpression<T2, m, E2>&);
    
    const T operator()(int) const;
    
  };


  //! Product between two expressions (elementwise)
  template<class T, int m, class E1, class E2>
  class TinyVectorProduct : public TinyVectorExpression<T, m, TinyVectorProduct<T, m, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyVectorProduct(const TinyVectorExpression<T1, m, E1>&, const TinyVectorExpression<T2, m, E2>&);
    
    const T operator()(int) const;
    
  };


  //! Division between two expressions (elementwise)
  template<class T, int m, class E1, class E2>
  class TinyVectorDivision : public TinyVectorExpression<T, m, TinyVectorDivision<T, m, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyVectorDivision(const TinyVectorExpression<T1, m, E1>&, const TinyVectorExpression<T2, m, E2>&);
    
    const T operator()(int) const;
    
  };


  //! Multiplication between a scalar and a vector
  template<class T, int m, class T0, class E>
  class TinyVectorScaled : public TinyVectorExpression<T, m, TinyVectorScaled<T, m, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyVectorScaled(const T0&, const TinyVectorExpression<T1, m, E>&);
    
    const T operator()(int) const;
    
  };


  //! Division between a scalar and a vector
  template<class T, int m, class T0, class E>
  class TinyVectorScalDiv : public TinyVectorExpression<T, m, TinyVectorScalDiv<T, m, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyVectorScalDiv(const T0&, const TinyVectorExpression<T1, m, E>&);
    
    const T operator()(int) const;
    
  };


  //! Addition between a scalar and a vector
  template<class T, int m, class T0, class E>
  class TinyVectorScalSum : public TinyVectorExpression<T, m, TinyVectorScalSum<T, m, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyVectorScalSum(const T0&, const TinyVectorExpression<T1, m, E>&);
    
    const T operator()(int) const;
    
  };


  //! Subtraction between a scalar and a vector
  template<class T, int m, class T0, class E>
  class TinyVectorScalDiff : public TinyVectorExpression<T, m, TinyVectorScalDiff<T, m, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyVectorScalDiff(const T0&, const TinyVectorExpression<T1, m, E>&);
    
    const T operator()(int) const;
    
  };


  //! Opposite of a vector
  template<class T, int m, class E>
  class TinyVectorOpposite : public TinyVectorExpression<T, m, TinyVectorOpposite<T, m, E> >
  {
    const E& u_;
    
  public:
    TinyVectorOpposite(const TinyVectorExpression<T, m, E>&);
    
    const T operator()(int) const;
    
  };
  
  
  /*************
   * Operators *
   *************/  
  
  
  template<class T, int m, class E1, class E2>
  const TinyVectorSum<T, m, E1, E2> 
  operator +(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorSum<complex<T>, m, E1, E2> 
  operator +(const TinyVectorExpression<complex<T>, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorSum<complex<T>, m, E1, E2> 
  operator +(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<complex<T>, m, E2>& v);

  template<class T, int m, class E1, class E2>
  inline const TinyVectorDifference<T, m, E1, E2>
  operator -(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  inline const TinyVectorDifference<complex<T>, m, E1, E2>
  operator -(const TinyVectorExpression<complex<T>, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorDifference<complex<T>, m, E1, E2>
  operator -(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<complex<T>, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorProduct<T, m, E1, E2> 
  operator *(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorProduct<complex<T>, m, E1, E2> 
  operator *(const TinyVectorExpression<complex<T>, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorProduct<complex<T>, m, E1, E2> 
  operator *(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<complex<T>, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorDivision<T, m, E1, E2>
  operator /(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorDivision<complex<T>, m, E1, E2>
  operator /(const TinyVectorExpression<complex<T>, m, E1>& u,
	     const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  const TinyVectorDivision<complex<T>, m, E1, E2>
  operator /(const TinyVectorExpression<T, m, E1>& u,
	     const TinyVectorExpression<complex<T>, m, E2>& v);

  template<class T, int m, class E>
  const TinyVectorScaled<T, m, T, E> operator *(const T& alpha,
						const TinyVectorExpression<T, m, E>& u);

  template<class T, int m, class E>
  const TinyVectorScaled<complex<T>, m, complex<T>, E>
  operator *(const complex<T>& alpha,
	     const TinyVectorExpression<T, m, E>& u);

  template<class T, int m, class E>
  const TinyVectorScaled<complex<T>, m, T, E>
  operator *(const T& alpha,
	     const TinyVectorExpression<complex<T>, m, E>& u);

  template<class T, int m, class E>
  const TinyVectorScaled<T, m, T, E> operator *(const TinyVectorExpression<T, m, E>& u,
						const T& alpha);

  template<class T, int m, class E>
  const TinyVectorScaled<complex<T>, m, complex<T>, E>
  operator *(const TinyVectorExpression<T, m, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, class E>
  const TinyVectorScaled<complex<T>, m, T, E>
  operator *(const TinyVectorExpression<complex<T>, m, E>& u,
	     const T& alpha);

  template<class T, int m, class E>
  const TinyVectorScalDiv<T, m, T, E> operator /(const T& alpha,
						 const TinyVectorExpression<T, m, E>& u);

  template<class T, int m, class E>
  const TinyVectorScalDiv<complex<T>, m, complex<T>, E>
  operator /(const complex<T>& alpha,
	     const TinyVectorExpression<T, m, E>& u);

  template<class T, int m, class E>
  const TinyVectorScalDiv<complex<T>, m, T, E>
  operator /(const T& alpha,
	     const TinyVectorExpression<complex<T>, m, E>& u);
  
  template<class T, int m, class E>
  const TinyVectorScaled<T, m, T, E> operator /(const TinyVectorExpression<T, m, E>& u,
						const T& alpha);

  template<class T, int m, class E>
  const TinyVectorScaled<complex<T>, m, complex<T>, E>
  operator /(const TinyVectorExpression<T, m, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, class E>
  const TinyVectorScaled<complex<T>, m, T, E> 
  operator /(const TinyVectorExpression<complex<T>, m, E>& u,
	     const T& alpha);

  template<class T, int m, class E>
  const TinyVectorScalSum<T, m, T, E> operator +(const T& alpha,
						 const TinyVectorExpression<T, m, E>& u);

  template<class T, int m, class E>
  const TinyVectorScalSum<complex<T>, m, complex<T>, E>
  operator +(const complex<T>& alpha,
	     const TinyVectorExpression<T, m, E>& u);

  template<class T, int m, class E>
  const TinyVectorScalSum<complex<T>, m, T, E>
  operator +(const T& alpha,
	     const TinyVectorExpression<complex<T>, m, E>& u);

  template<class T, int m, class E>
  const TinyVectorScalSum<T, m, T, E> operator +(const TinyVectorExpression<T, m, E>& u,
						 const T& alpha);

  template<class T, int m, class E>
  const TinyVectorScalSum<complex<T>, m, complex<T>, E>
  operator +(const TinyVectorExpression<T, m, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, class E>
  const TinyVectorScalSum<complex<T>, m, T, E>
  operator +(const TinyVectorExpression<complex<T>, m, E>& u,
	     const T& alpha);

  template<class T, int m, class E>
  const TinyVectorScalSum<T, m, T, E> operator -(const TinyVectorExpression<T, m, E>& u,
						 const T& alpha);

  template<class T, int m, class E>
  const TinyVectorScalSum<complex<T>, m, complex<T>, E>
  operator -(const TinyVectorExpression<T, m, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, class E>
  const TinyVectorScalSum<complex<T>, m, T, E>
  operator -(const TinyVectorExpression<complex<T>, m, E>& u,
	     const T& alpha);

  template<class T, int m, class E> 
  const TinyVectorScalDiff<T, m, T, E> operator -(const T& alpha,
						  const TinyVectorExpression<T, m, E>& u);

  template<class T, int m, class E> 
  const TinyVectorScalDiff<complex<T>, m, complex<T>, E>
  operator -(const complex<T>& alpha,
	     const TinyVectorExpression<T, m, E>& u);

  template<class T, int m, class E> 
  const TinyVectorScalDiff<complex<T>, m, T, E>
  operator -(const T& alpha,
	     const TinyVectorExpression<complex<T>, m, E>& u);

  template<class T, int m, class E>
  inline const TinyVectorOpposite<T, m, E> operator-(const TinyVectorExpression<T, m, E>& u);
  
}

#define SELDON_FILE_TINY_VECTOR_EXPRESSION_HXX
#endif
