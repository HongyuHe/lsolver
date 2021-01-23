#ifndef SELDON_FILE_VECTOR_EXPRESSION_HXX

namespace Seldon
{
  
  //! Expression between vectors
  template<class T, class E>
  class VectorExpression
  {
  public:
    long GetSize() const;
    const T operator()(long) const;
    
    inline operator E&() { return static_cast<E&>(*this); }
    inline operator E const&() const { return static_cast<const E&>(*this); }
    
  };


  //! Difference between two expressions
  template<class T, class E1, class E2>
  class VectorDifference : public VectorExpression<T, VectorDifference<T, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    VectorDifference(const VectorExpression<T1, E1>&, const VectorExpression<T2, E2>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };


  //! Sum between two expressions
  template<class T, class E1, class E2>
  class VectorSum : public VectorExpression<T, VectorSum<T, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    VectorSum(const VectorExpression<T1, E1>&, const VectorExpression<T2, E2>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };


  //! Product between two expressions (elementwise)
  template<class T, class E1, class E2>
  class VectorProduct : public VectorExpression<T, VectorProduct<T, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    VectorProduct(const VectorExpression<T1, E1>&, const VectorExpression<T2, E2>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };


  //! Division between two expressions (elementwise)
  template<class T, class E1, class E2>
  class VectorDivision : public VectorExpression<T, VectorDivision<T, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    VectorDivision(const VectorExpression<T1, E1>&, const VectorExpression<T2, E2>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };


  //! Multiplication between a scalar and a vector
  template<class T, class T0, class E>
  class VectorScaled : public VectorExpression<T, VectorScaled<T, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    VectorScaled(const T0&, const VectorExpression<T1, E>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };


  //! Division between a scalar and a vector
  template<class T, class T0, class E>
  class VectorScalDiv : public VectorExpression<T, VectorScalDiv<T, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    VectorScalDiv(const T0&, const VectorExpression<T1, E>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };


  //! Addition between a scalar and a vector
  template<class T, class T0, class E>
  class VectorScalSum : public VectorExpression<T, VectorScalSum<T, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    VectorScalSum(const T0&, const VectorExpression<T1, E>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };


  //! Subtraction between a scalar and a vector
  template<class T, class T0, class E>
  class VectorScalDiff : public VectorExpression<T, VectorScalDiff<T, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    VectorScalDiff(const T0&, const VectorExpression<T1, E>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };


  //! Opposite of a vector
  template<class T, class E>
  class VectorOpposite : public VectorExpression<T, VectorOpposite<T, E> >
  {
    const E& u_;
    
  public:
    VectorOpposite(const VectorExpression<T, E>&);
    
    long GetSize() const;
    const T operator()(long) const;
    
  };
  
  
  /*************
   * Operators *
   *************/
  
  
  template<class T, class E1, class E2>
  const VectorSum<T, E1, E2> operator+(const VectorExpression<T, E1>&, const VectorExpression<T, E2>&);

  template<class T, class E1, class E2>
  const VectorSum<complex<T>, E1, E2> operator+(const VectorExpression<complex<T>, E1>&,
						const VectorExpression<T, E2>&);

  template<class T, class E1, class E2>
  const VectorSum<complex<T>, E1, E2> operator+(const VectorExpression<T, E1>&,
						const VectorExpression<complex<T>, E2>&);


  template<class T, class E1, class E2>
  const VectorDifference<T, E1, E2> operator-(const VectorExpression<T, E1>&, const VectorExpression<T, E2>&);

  template<class T, class E1, class E2>
  const VectorDifference<complex<T>, E1, E2> operator-(const VectorExpression<complex<T>, E1>&,
						       const VectorExpression<T, E2>&);

  template<class T, class E1, class E2>
  const VectorDifference<complex<T>, E1, E2> operator-(const VectorExpression<T, E1>&,
						       const VectorExpression<complex<T>, E2>&);


  template<class T, class E1, class E2>
  const VectorProduct<T, E1, E2> operator*(const VectorExpression<T, E1>&, const VectorExpression<T, E2>&);

  template<class T, class E1, class E2>
  const VectorProduct<complex<T>, E1, E2> operator*(const VectorExpression<complex<T>, E1>&,
						    const VectorExpression<T, E2>&);

  template<class T, class E1, class E2>
  const VectorProduct<complex<T>, E1, E2> operator*(const VectorExpression<T, E1>&,
						    const VectorExpression<complex<T>, E2>&);


  template<class T, class E1, class E2>
  const VectorDivision<T, E1, E2> operator/(const VectorExpression<T, E1>&, const VectorExpression<T, E2>&);

  template<class T, class E1, class E2>
  const VectorDivision<complex<T>, E1, E2> operator/(const VectorExpression<complex<T>, E1>&,
						     const VectorExpression<T, E2>&);

  template<class T, class E1, class E2>
  const VectorDivision<complex<T>, E1, E2> operator/(const VectorExpression<T, E1>&,
						     const VectorExpression<complex<T>, E2>&);


  template<class T, class E>
  const VectorScaled<T, T, E> operator*(const T&, const VectorExpression<T, E>&);

  template<class T, class E>
  const VectorScaled<complex<T>, complex<T>, E> operator*(const complex<T>&,
							  const VectorExpression<T, E>&);

  template<class T, class E>
  const VectorScaled<complex<T>, T, E> operator*(const T&,
						 const VectorExpression<complex<T>, E>&);
  
  template<class T, class E>
  const VectorScaled<T, T, E> operator*(const VectorExpression<T, E>&, const T&);

  template<class T, class E>
  const VectorScaled<complex<T>, complex<T>, E> operator*(const VectorExpression<T, E>&,
							  const complex<T>&);

  template<class T, class E>
  const VectorScaled<complex<T>, T, E> operator*(const VectorExpression<complex<T>, E>&, const T&);


  template<class T, class E>
  const VectorScalDiv<T, T, E> operator/(const T&, const VectorExpression<T, E>&);

  template<class T, class E>
  const VectorScalDiv<complex<T>, complex<T>, E> operator/(const complex<T>&,
							   const VectorExpression<T, E>&);

  template<class T, class E>
  const VectorScalDiv<complex<T>, T, E> operator/(const T&,
						 const VectorExpression<complex<T>, E>&);
  
  template<class T, class E>
  const VectorScaled<T, T, E> operator/(const VectorExpression<T, E>&, const T&);
  
  template<class T, class E>
  const VectorScaled<complex<T>, complex<T>, E> operator/(const VectorExpression<T, E>&,
							  const complex<T>&);

  template<class T, class E>
  const VectorScaled<complex<T>, T, E> operator/(const VectorExpression<complex<T>, E>&, const T&);


  template<class T, class E>
  const VectorScalSum<T, T, E> operator+(const T&, const VectorExpression<T, E>&);

  template<class T, class E>
  const VectorScalSum<complex<T>, complex<T>, E> operator+(const complex<T>&,
							  const VectorExpression<T, E>&);

  template<class T, class E>
  const VectorScalSum<complex<T>, T, E> operator+(const T&,
						 const VectorExpression<complex<T>, E>&);
  
  template<class T, class E>
  const VectorScalSum<T, T, E> operator+(const VectorExpression<T, E>&, const T&);

  template<class T, class E>
  const VectorScalSum<complex<T>, complex<T>, E> operator+(const VectorExpression<T, E>&,
							  const complex<T>&);

  template<class T, class E>
  const VectorScalSum<complex<T>, T, E> operator+(const VectorExpression<complex<T>, E>&, const T&);


  template<class T, class E>
  const VectorScalDiff<T, T, E> operator-(const T&, const VectorExpression<T, E>&);

  template<class T, class E>
  const VectorScalDiff<complex<T>, complex<T>, E> operator-(const complex<T>&,
							    const VectorExpression<T, E>&);

  template<class T, class E>
  const VectorScalDiff<complex<T>, T, E> operator-(const T&,
						   const VectorExpression<complex<T>, E>&);
  
  template<class T, class E>
  const VectorScalSum<T, T, E> operator-(const VectorExpression<T, E>&, const T&);

  template<class T, class E>
  const VectorScalSum<complex<T>, complex<T>, E> operator-(const VectorExpression<T, E>&,
							  const complex<T>&);

  template<class T, class E>
  const VectorScalSum<complex<T>, T, E> operator-(const VectorExpression<complex<T>, E>&, const T&);


  template<class T, class E>
  const VectorOpposite<T, E> operator-(const VectorExpression<T, E>&);
  
}

#define SELDON_FILE_VECTOR_EXPRESSION_HXX
#endif
