#ifndef SELDON_FILE_TINY_MATRIX_EXPRESSION_HXX

namespace Seldon
{
  
  //! Expression between vectors
  template<class T, int m, int n, class E>
  class TinyMatrixExpression
  {
  public:
    const T operator()(int, int) const;
    
    inline operator E&() { return static_cast<E&>(*this); }
    inline operator E const&() const { return static_cast<const E&>(*this); }
    
  };


  //! Difference between two expressions
  template<class T, int m, int n, class E1, class E2>
  class TinyMatrixDifference : public TinyMatrixExpression<T, m, n, TinyMatrixDifference<T, m, n, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyMatrixDifference(const TinyMatrixExpression<T1, m, n, E1>&, const TinyMatrixExpression<T2, m, n, E2>&);
    
    const T operator()(int, int) const;
    
  };


  //! Sum between two expressions
  template<class T, int m, int n, class E1, class E2>
  class TinyMatrixSum : public TinyMatrixExpression<T, m, n, TinyMatrixSum<T, m, n, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyMatrixSum(const TinyMatrixExpression<T1, m, n, E1>&, const TinyMatrixExpression<T2, m, n, E2>&);
    
    const T operator()(int, int) const;
    
  };


  //! Product between two expressions (elementwise)
  template<class T, int m, int n, class E1, class E2>
  class TinyMatrixProduct : public TinyMatrixExpression<T, m, n, TinyMatrixProduct<T, m, n, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyMatrixProduct(const TinyMatrixExpression<T1, m, n, E1>&, const TinyMatrixExpression<T2, m, n, E2>&);
    
    const T operator()(int, int) const;
    
  };


  //! Division between two expressions (elementwise)
  template<class T, int m, int n, class E1, class E2>
  class TinyMatrixDivision : public TinyMatrixExpression<T, m, n, TinyMatrixDivision<T, m, n, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyMatrixDivision(const TinyMatrixExpression<T1, m, n, E1>&, const TinyMatrixExpression<T2, m, n, E2>&);
    
    const T operator()(int, int) const;
    
  };


  //! Multiplication between a scalar and a matrix
  template<class T, int m, int n, class T0, class E>
  class TinyMatrixScaled : public TinyMatrixExpression<T, m, n, TinyMatrixScaled<T, m, n, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyMatrixScaled(const T0&, const TinyMatrixExpression<T1, m, n, E>&);
    
    const T operator()(int, int) const;
    
  };


  //! Division between a scalar and a matrix
  template<class T, int m, int n, class T0, class E>
  class TinyMatrixScalDiv : public TinyMatrixExpression<T, m, n, TinyMatrixScalDiv<T, m, n, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyMatrixScalDiv(const T0&, const TinyMatrixExpression<T1, m, n, E>&);
    
    const T operator()(int, int) const;
    
  };


  //! Addition between a scalar and a matrix
  template<class T, int m, int n, class T0, class E>
  class TinyMatrixScalSum : public TinyMatrixExpression<T, m, n, TinyMatrixScalSum<T, m, n, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyMatrixScalSum(const T0&, const TinyMatrixExpression<T1, m, n, E>&);
    
    const T operator()(int, int) const;
    
  };


  //! Subtraction between a scalar and a matrix
  template<class T, int m, int n, class T0, class E>
  class TinyMatrixScalDiff : public TinyMatrixExpression<T, m, n, TinyMatrixScalDiff<T, m, n, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyMatrixScalDiff(const T0&, const TinyMatrixExpression<T1, m, n, E>&);
    
    const T operator()(int, int) const;
    
  };


  //! Opposite of a matrix
  template<class T, int m, int n, class E>
  class TinyMatrixOpposite : public TinyMatrixExpression<T, m, n, TinyMatrixOpposite<T, m, n, E> >
  {
    const E& u_;
    
  public:
    TinyMatrixOpposite(const TinyMatrixExpression<T, m, n, E>&);
    
    const T operator()(int, int) const;
    
  };


  //! Transpose of a matrix
  template<class T, int m, int n, class E>
  class TinyMatrixTranspose : public TinyMatrixExpression<T, m, n, TinyMatrixTranspose<T, m, n, E> >
  {
    const E& u_;
    
  public:
    TinyMatrixTranspose(const TinyMatrixExpression<T, n, m, E>&);
    
    const T operator()(int, int) const;
    
  };

  
  /*************
   * Operators *
   *************/  
  
  
  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixSum<T, m, n, E1, E2> 
  operator +(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixSum<complex<T>, m, n, E1, E2> 
  operator +(const TinyMatrixExpression<complex<T>, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixSum<complex<T>, m, n, E1, E2> 
  operator +(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<complex<T>, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixDifference<T, m, n, E1, E2>
  operator -(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  inline const TinyMatrixDifference<complex<T>, m, n, E1, E2>
  operator -(const TinyMatrixExpression<complex<T>, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixDifference<complex<T>, m, n, E1, E2>
  operator -(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<complex<T>, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixProduct<T, m, n, E1, E2> 
  operator *(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixProduct<complex<T>, m, n, E1, E2> 
  operator *(const TinyMatrixExpression<complex<T>, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixProduct<complex<T>, m, n, E1, E2> 
  operator *(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<complex<T>, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixDivision<T, m, n, E1, E2>
  operator /(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixDivision<complex<T>, m, n, E1, E2>
  operator /(const TinyMatrixExpression<complex<T>, m, n, E1>& u,
	     const TinyMatrixExpression<T, m, n, E2>& v);

  template<class T, int m, int n, class E1, class E2>
  const TinyMatrixDivision<complex<T>, m, n, E1, E2>
  operator /(const TinyMatrixExpression<T, m, n, E1>& u,
	     const TinyMatrixExpression<complex<T>, m, n, E2>& v);

  template<class T, int m, int n, class E>
  const TinyMatrixScaled<T, m, n, T, E> operator *(const T& alpha,
						const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E>
  const TinyMatrixScaled<complex<T>, m, n, complex<T>, E>
  operator *(const complex<T>& alpha,
	     const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E>
  const TinyMatrixScaled<complex<T>, m, n, T, E>
  operator *(const T& alpha,
	     const TinyMatrixExpression<complex<T>, m, n, E>& u);

  template<class T, int m, int n, class E>
  const TinyMatrixScaled<T, m, n, T, E> operator *(const TinyMatrixExpression<T, m, n, E>& u,
						const T& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScaled<complex<T>, m, n, complex<T>, E>
  operator *(const TinyMatrixExpression<T, m, n, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScaled<complex<T>, m, n, T, E>
  operator *(const TinyMatrixExpression<complex<T>, m, n, E>& u,
	     const T& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScalDiv<T, m, n, T, E> operator /(const T& alpha,
						 const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E>
  const TinyMatrixScalDiv<complex<T>, m, n, complex<T>, E>
  operator /(const complex<T>& alpha,
	     const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E>
  const TinyMatrixScalDiv<complex<T>, m, n, T, E>
  operator /(const T& alpha,
	     const TinyMatrixExpression<complex<T>, m, n, E>& u);
  
  template<class T, int m, int n, class E>
  const TinyMatrixScaled<T, m, n, T, E> operator /(const TinyMatrixExpression<T, m, n, E>& u,
						const T& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScaled<complex<T>, m, n, complex<T>, E>
  operator /(const TinyMatrixExpression<T, m, n, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScaled<complex<T>, m, n, T, E> 
  operator /(const TinyMatrixExpression<complex<T>, m, n, E>& u,
	     const T& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<T, m, n, T, E> operator +(const T& alpha,
						 const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>
  operator +(const complex<T>& alpha,
	     const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<complex<T>, m, n, T, E>
  operator +(const T& alpha,
	     const TinyMatrixExpression<complex<T>, m, n, E>& u);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<T, m, n, T, E> operator +(const TinyMatrixExpression<T, m, n, E>& u,
						 const T& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>
  operator +(const TinyMatrixExpression<T, m, n, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<complex<T>, m, n, T, E>
  operator +(const TinyMatrixExpression<complex<T>, m, n, E>& u,
	     const T& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<T, m, n, T, E> operator -(const TinyMatrixExpression<T, m, n, E>& u,
						 const T& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<complex<T>, m, n, complex<T>, E>
  operator -(const TinyMatrixExpression<T, m, n, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, int n, class E>
  const TinyMatrixScalSum<complex<T>, m, n, T, E>
  operator -(const TinyMatrixExpression<complex<T>, m, n, E>& u,
	     const T& alpha);

  template<class T, int m, int n, class E> 
  const TinyMatrixScalDiff<T, m, n, T, E> operator -(const T& alpha,
						  const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E> 
  const TinyMatrixScalDiff<complex<T>, m, n, complex<T>, E>
  operator -(const complex<T>& alpha,
	     const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E> 
  const TinyMatrixScalDiff<complex<T>, m, n, T, E>
  operator -(const T& alpha,
	     const TinyMatrixExpression<complex<T>, m, n, E>& u);

  template<class T, int m, int n, class E>
  inline const TinyMatrixOpposite<T, m, n, E> operator-(const TinyMatrixExpression<T, m, n, E>& u);

  template<class T, int m, int n, class E>
  inline const TinyMatrixTranspose<T, m, n, E> transpose(const TinyMatrixExpression<T, n, m, E>& u);
  
}

#define SELDON_FILE_TINY_MATRIX_EXPRESSION_HXX
#endif
