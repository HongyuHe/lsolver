#ifndef SELDON_FILE_TINY_VECTOR_HXX

#include "TinyVectorExpression.hxx"

namespace Seldon
{
  
  /**********************
   * TinyVector classes *
   **********************/

 
  template<int n>
  class TinyVectorLoop;

  template<class T, int m>
  class TinyVector;
  
  //! Class storing a tiny vector whose size is known at compilation time
  template <class T, int m_>
  class TinyVector_Base : public TinyVectorExpression<T, m_, TinyVector<T, m_> >
  {
  public :
    typedef T value_type;
    
    template<int n>
    friend class TinyVectorLoop;
    
  protected: 
    // pointer data_;
    T data_[m_];
    
  public :
    TinyVector_Base();
    
    explicit TinyVector_Base(const T& a);
    explicit TinyVector_Base(const T& a, const T& b);
    TinyVector_Base(const T& a, const T& b, const T& c);
    TinyVector_Base(const T& a, const T& b, const T& c, const T& d);
    TinyVector_Base(const T& a, const T& b, const T& c, const T& d, const T& e);
    
    void Init(const T & a);
    void Init(const T & a, const T& b);
    void Init(const T & a, const T& b, const T& c);
    void Init(const T & a, const T& b, const T& c, const T& d);
    void Init(const T & a, const T& b, const T& c, const T& d, const T& e);
    
    // Basic Functions
    int GetM() const;
    int GetLength() const;
    int GetSize() const;
    size_t GetMemorySize() const;
    
    TinyVector<T, m_>& GetLeafClass();
    const TinyVector<T, m_>& GetLeafClass() const;
    
    // Element access and affection
    const T& operator()(int i) const;
    T& operator()(int i);

    // Convenient Functions
    void Print() const; // we display the coordinates
    void Print(ostream&) const; // we display the coordinates
    void Zero() ;  // we put zero on coordinates
    void Fill();
    
    template<class T0>
    void Fill(const T0& x);
    void FillRand();
    bool IsZero() const;
    
    // Norms
    typename ClassComplexType<T>::Treal DistanceSquare(const TinyVector<T, m_> & p) const;
    typename ClassComplexType<T>::Treal Distance(const TinyVector<T, m_> & p) const;
    
  };

  
  //! vector with real/complex components
  template <class T, int m_>
  class TinyVector : public TinyVector_Base<T, m_>
  {
  public :
    static typename ClassComplexType<T>::Treal threshold;
    
    TinyVector();
    
    explicit TinyVector(const T& a);
    explicit TinyVector(const T& a, const T& b);
    TinyVector(const T& a, const T& b, const T& c);
    TinyVector(const T& a, const T& b, const T& c, const T& d);
    TinyVector(const T& a, const T& b, const T& c, const T& d, const T& e);

    template<class E>
    TinyVector(const TinyVectorExpression<T, m_, E>& u);
    
    // operation this = u
    template<class T1, class E>
    TinyVector<T, m_>& operator =(const TinyVectorExpression<T1, m_, E> & u);

    // this += u
    template<class T1, class E>
    TinyVector<T, m_>& operator +=(const TinyVectorExpression<T1, m_, E> & u);

    // this -= u
    template<class T1, class E>
    TinyVector<T, m_>& operator -=(const TinyVectorExpression<T1, m_, E> & u);
    
    // fills vector with x
    TinyVector<T, m_>& operator =(const T& x);

    // operation this = a*this where a is scalar
    template<class T1>
    TinyVector<T, m_>& operator *=( const T1& a );
    
    // operation this = 1/a*this where a is scalar
    template<class T1>
    TinyVector<T, m_>& operator /=( const T1& a );
            
  };


  //! vector with integer components
  template<int m_>
  class TinyVector<int, m_> : public TinyVector_Base<int, m_>
  {
  public :
    TinyVector();
    
    explicit TinyVector(const int& a);
    explicit TinyVector(const int& a, const int& b);
    TinyVector(const int& a, const int& b, const int& c);
    TinyVector(const int& a, const int& b, const int& c, const int& d);
    TinyVector(const int& a, const int& b, const int& c, const int& d, const int& e);

  };


  /*************
   * Operators *
   *************/


  // Operators to compare expressions
  template<int m>
  bool operator==(const TinyVector<int, m>& u, const TinyVector<int, m>& v);

  template<int m>
  bool operator!=(const TinyVector<int, m> & u, const TinyVector<int, m>& v);

  template<int m>
  bool operator<=(const TinyVector<int, m> & u, const TinyVector<int, m>& v);

  template<int m>
  bool operator<(const TinyVector<int, m> & u, const TinyVector<int, m>& v);

  template<int m>
  bool operator>=(const TinyVector<int, m> & u, const TinyVector<int, m>& v);

  template<int m>
  bool operator>(const TinyVector<int, m> & u, const TinyVector<int, m>& v);

  template<class T, int m, class E1, class E2>
  bool operator==(const TinyVectorExpression<T, m, E1>& u,
		  const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  bool operator!=(const TinyVectorExpression<T, m, E1> & u,
		  const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  bool operator<=(const TinyVectorExpression<T, m, E1> & u,
		  const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  bool operator<(const TinyVectorExpression<T, m, E1> & u,
		 const TinyVectorExpression<T, m, E2>& v);
    
  template<class T, int m, class E1, class E2>
  bool operator>=(const TinyVectorExpression<int, m, E1> & u,
		  const TinyVectorExpression<int, m, E2>& v);

  template<class T, int m, class E1, class E2>
  bool operator>(const TinyVectorExpression<T, m, E1> & u,
		 const TinyVectorExpression<T, m, E2>& v);

  template <class T, int m>
  ostream& operator <<(ostream& out,const TinyVector<T,m>& V);

  template <class T, int m>
  istream& operator >>(istream& out, TinyVector<T,m>& V);
  
  // threshold for comparing vectors
#ifdef MONTJOIE_WITH_MPFR
  template<class T, int m>
  typename ClassComplexType<T>::Treal TinyVector<T, m>::threshold(0, MONTJOIE_DEFAULT_PRECISION);
#else
  template<class T, int m>
  typename ClassComplexType<T>::Treal TinyVector<T, m>::threshold(0);
#endif


  /**************************
   * TinyVectorLoop classes *
   **************************/


  //! class used for unrolling loops for operators/functions of TinyVector
  template<int n>
  class TinyVectorLoop
  {
  public :
    template<int m, class T0>
    static void Zero(TinyVector<T0, m>& x);
    
    template<int m, class T0>
    static void FillGreater(TinyVector<T0, m>& x, const T0& a);

    template<int m, class T0>
    static void InitValue(TinyVector<T0, m>& x, const T0& a);

    template<int m, class T0, class T1>
    static void Fill(TinyVector<T0, m>& x, const T1& a);
    
    template<int m, class T0>
    static void Fill(TinyVector<T0, m>& x);
    
    template<int m, class T0>
    static void FillRand(TinyVector<T0, m>& x);
    
    template<int m, class T0, class E0>
    static bool IsZero(const TinyVectorExpression<T0, m, E0>& x);
    
    template<int m, class T1, class E, class T0>
    static void Copy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y);

    template<int m, class T1, class E, class T0>
    static void AddCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y);

    template<int m, class T1, class E, class T0>
    static void DiffCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y);

    template<int m, class T0, class T1>
    static void 
    CopyVector(const TinyVector<Vector<T0>, m>& uvec, int j, TinyVector<T1, m>& u);
    
    template<int m, class T0, class T1>
    static void 
    CopyVector(const TinyVector<T1, m>& u, int j, TinyVector<Vector<T0>, m>& uvec);
    
    template<int m, class T0, class T1>
    static void CopyVector(const Vector<T0>& Vh, int j, TinyVector<T1, m>& tmp);

    template<int m, class T>
    static void CopyVector(const Vector<complex<T> >& Vh, int j, TinyVector<T, m>& vr, TinyVector<T, m>& vi);

    template<int m, class T0, class T1>
    static void CopyVector(const TinyVector<T1, m>& tmp, int j, Vector<T0>& Vh);

    template<int m, class T0, class T1>
    static void AddVector(const Vector<T0>& Vh, int j, TinyVector<T1, m>& tmp);

    template<int m, class T0, class T1>
    static void AddVector(const TinyVector<T1, m>& tmp, int j, Vector<T0>& Vh);

    template<int m, class T0, class T1>
    static void ExtractVector(const Vector<T0>& Vh, int j, TinyVector<T1, m>& tmp);

    template<int m, class T0, class T1>
    static void ExtractVector(const TinyVector<T1, m>& tmp, int j, Vector<T0>& Vh);

    template<int m, class T0, class T1>
    static void ExtractVector(const Vector<Vector<T0> >& u,
			      int j, int offset, TinyVector<T1, m>& v);

    template<int m, class T0, class T1>
    static void CopyVector(const TinyVector<T1, m>& v,
			   int j, int offset, Vector<Vector<T0> >& u);
    
    template<int p, int m, class T0, class T1>
    static void ExtractVector(const Vector<Vector<TinyVector<T0, p> > >& u,
			      int j, int k, int offset, TinyVector<T1, m>& v);
    
    template<int m, class T0, class T1>
    static void ExtractVector(const Vector<T0>& Vh, const IVect& row_num,
                              int j, TinyVector<T1, m>& tmp);

    template<int m, class T0, class T1>
    static void ExtractVector(const TinyVector<T1, m>& tmp,
                              const IVect& row_num, int j, Vector<T0>& Vh);

    template<class T0, class T1>
    static void ExtractVector(const Vector<T0>& x_src, int offset_source,
			      Vector<T1>& x_dest, int offset_dest);
    
    template<int m, class T0, class T1>
    static void AddVector(const T0& alpha, const TinyVector<T1, m>& tmp,
                          const IVect& row_num, int j, Vector<T0>& Vh);

    template<int m, class T0, class T1>
    static void ExtractVectorConj(const Vector<T0>& Vh, const IVect& row_num,
				  int j, TinyVector<T1, m>& tmp);

    template<int m, class T0, class T1>
    static void ExtractVectorConj(const TinyVector<T1, m>& tmp,
                                  const IVect& row_num, int j, Vector<T0>& Vh);

    template<int m, class T0, class T1>
    static void AddVectorConj(const T0& alpha, const TinyVector<T1, m>& tmp,
                              const IVect& row_num, int j, Vector<T0>& Vh);
    
    template<int m, class T0, class T1>
    static void MltVector(const T1& alpha, int j, TinyVector<Vector<T0>, m>& uvec);

    template<int m, class T0, class T1, class E0, class E1>
    static bool Equal(const TinyVectorExpression<T0, m, E0>& x,
		      const TinyVectorExpression<T1, m, E1>& y);
    
    template<int m, class T0, class T1, class E0, class E1>
    static bool LessOrEqual(const TinyVectorExpression<T0, m, E0>& x,
			    const TinyVectorExpression<T1, m, E1>& y);
    
    template<int m, class T0, class T1, class E0, class E1>
    static bool Less(const TinyVectorExpression<T0, m, E0>& x,
		     const TinyVectorExpression<T1, m, E1>& y);
    
    template<int m, class T0, class T1>
    static bool EqualInt(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y);
    
    template<int m, class T0, class T1>
    static bool LessOrEqualInt(const TinyVector<T0, m>& x,
			       const TinyVector<T1, m>& y);
    
    template<int m, class T0, class T1>
    static bool LessInt(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y);
    
    template<int m, class T0, class T1>
    static void Mlt(const T1& alpha, TinyVector<T0, m>& x);
    
    template<int m, class T0, class T1, class T2, class E0, class E1>
    static void DotProd(const TinyVectorExpression<T0, m, E0>& u,
			const TinyVectorExpression<T1, m, E1>& v, T2& scal);
    
    template<int m, class T0, class T1, class T2, class E0, class E1>
    static void DotProdConj(const TinyVectorExpression<T0, m, E0>& u,
			    const TinyVectorExpression<T1, m, E1>& v, T2& scal);

    template<int m, class T0, class T1, class E0>
    static void Abs(const TinyVectorExpression<T0, m, E0>& u, T1& scal);

    template<int m, class T0, class T1, class E0>
    static void NormInf(const TinyVectorExpression<T0, m, E0>& u, T1& scal);

    template<int m, class T0, class T1, class E0>
    static void AbsSquare(const TinyVectorExpression<T0, m, E0>& u, T1& scal);
    
    template<int m, class T>
    static void Print(const TinyVector<T, m>& V, ostream& out);

    template<int m, class T>
    static ostream& WriteText(const TinyVector<T, m>& V, ostream& out);

    template<int m, class T>
    static istream& ReadText(TinyVector<T, m>& V, istream& in);

    template<class T, int m>
    static void UpdateMinimum(const TinyVector<T, m>& u, TinyVector<T, m>& v);

    template<class T, int m>
    static void UpdateMaximum(const TinyVector<T, m>& u, TinyVector<T, m>& v);
    
  };
  
  //! Class used to terminate the loop
  template<>
  class TinyVectorLoop<0>
  {
  public :
    template<int m, class T0>
    static inline void Zero(TinyVector<T0, m>& x) {}

    template<int m, class T0, class E0>
    static inline bool IsZero(const TinyVectorExpression<T0, m, E0>& x) { return true; }
    
    template<int m, class T0>
    static inline void FillGreater(TinyVector<T0, m>& x, const T0& a){}

    template<int m, class T0>
    static inline void InitValue(TinyVector<T0, m>& x, const T0& a){}

    template<int m, class T0, class T1>
    static inline void Fill(TinyVector<T0, m>& x, const T1& a){}
    
    template<int m, class T0>
    static inline void Fill(TinyVector<T0, m>& x){}
    
    template<int m, class T0>
    static inline void FillRand(TinyVector<T0, m>& x){}
    
    template<int m, class T1, class E, class T0>
    static inline void Copy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y) {} 

    template<int m, class T1, class E, class T0>
    static inline void AddCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y) {} 

    template<int m, class T1, class E, class T0>
    static inline void DiffCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y) {} 
    
    template<int m, class T0, class T1>
    static inline void CopyVector(const TinyVector<Vector<T0>, m>& uvec,
                                  int j, TinyVector<T1, m>& u) { }
    
    template<int m, class T0, class T1>
    static inline void CopyVector(const TinyVector<T1, m>& u, int j,
                                  TinyVector<Vector<T0>, m>& uvec) {}
    
    template<int m, class T0, class T1>
    static inline void CopyVector(const Vector<T0>& Vh, int j, TinyVector<T1, m>& tmp) {}

    template<int m, class T>
    static inline void CopyVector(const Vector<complex<T> >& Vh, int j,
				  TinyVector<T, m>& vr, TinyVector<T, m>& vi){}
    
    template<int m, class T0, class T1>
    static inline void CopyVector(const TinyVector<T1, m>& tmp, int j, Vector<T0>& Vh) {}

    template<int m, class T0, class T1>
    static inline void AddVector(const Vector<T0>& Vh, int j, TinyVector<T1, m>& tmp) {}

    template<int m, class T0, class T1>
    static inline void AddVector(const TinyVector<T1, m>& tmp, int j, Vector<T0>& Vh) {}

    template<int m, class T0, class T1>
    static inline void ExtractVector(const Vector<T0>& Vh, int j, TinyVector<T1, m>& tmp) {}

    template<int m, class T0, class T1>
    static inline void ExtractVector(const TinyVector<T1, m>& tmp, int j, Vector<T0>& Vh) {}

    template<int m, class T0, class T1>
    static inline void ExtractVector(const Vector<T0>& Vh, const IVect& row_num,
                              int j, TinyVector<T1, m>& tmp) {}

    template<int m, class T0, class T1>
    static inline void ExtractVector(const Vector<Vector<T0> >& u,
				     int j, int offset, TinyVector<T1, m>& v) {}
    
    template<int m, class T0, class T1>
    static inline void CopyVector(const TinyVector<T1, m>& v,
				  int j, int offset, Vector<Vector<T0> >& u) {}
    
    template<int p, int m, class T0, class T1>
    static inline void ExtractVector(const Vector<Vector<TinyVector<T0, p> > >& u,
				     int j, int k, int offset, TinyVector<T1, m>& v) {}

    template<int m, class T0, class T1>
    static inline void ExtractVector(const TinyVector<T1, m>& tmp,
                                     const IVect& row_num, int j, Vector<T0>& Vh) {}

    template<class T0, class T1>
    static inline void ExtractVector(const Vector<T0>& x_src, int offset_source,
				     Vector<T1>& x_dest, int offset_dest) {}
    
    template<int m, class T0, class T1>
    static inline void AddVector(const T0& alpha, const TinyVector<T1, m>& tmp,
                                 const IVect& row_num, int j, Vector<T0>& Vh) {}

    template<int m, class T0, class T1>
    static inline void ExtractVectorConj(const Vector<T0>& Vh, const IVect& row_num,
                                         int j, TinyVector<T1, m>& tmp) {}
    
    template<int m, class T0, class T1>
    static inline void ExtractVectorConj(const TinyVector<T1, m>& tmp,
                                         const IVect& row_num, int j, Vector<T0>& Vh) {}
    
    template<int m, class T0, class T1>
    static inline void AddVectorConj(const T0& alpha, const TinyVector<T1, m>& tmp,
                                     const IVect& row_num, int j, Vector<T0>& Vh) {}
    
    template<int m, class T0, class T1>
    static inline void MltVector(const T1& alpha, int j, TinyVector<Vector<T0>, m>& uvec) {}
    
    template<int m, class T0, class T1, class E0, class E1>
    static inline bool Equal(const TinyVectorExpression<T0, m, E0>& x,
			     const TinyVectorExpression<T1, m, E1>& y)
    {
      return true;
    }

    template<int m, class T0, class T1, class E0, class E1>
    static inline bool LessOrEqual(const TinyVectorExpression<T0, m, E0>& x,
				   const TinyVectorExpression<T1, m, E1>& y)
    {
      return true;
    }

    template<int m, class T0, class T1, class E0, class E1>
    static inline bool Less(const TinyVectorExpression<T0, m, E0>& x,
			    const TinyVectorExpression<T1, m, E1>& y)
    {
      return false;
    }

    template<int m, class T0, class T1>
    static inline bool EqualInt(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y)
    {
      return true;
    }

    template<int m, class T0, class T1>
    static inline bool LessOrEqualInt(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y)
    {
      return true;
    }

    template<int m, class T0, class T1>
    static inline bool LessInt(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y)
    {
      return false;
    }

    template<int m, class T0, class T1>
    static inline void Mlt(const T1& alpha, TinyVector<T0, m>& x) {}
    
    template<int m, class T0, class T1, class T2, class E0, class E1>
    static inline void DotProd(const TinyVectorExpression<T0, m, E0>& u,
			       const TinyVectorExpression<T1, m, E1>& v, T2& scal) {}
    
    template<int m, class T0, class T1, class E0, class E1, class T2>
    static inline void DotProdConj(const TinyVectorExpression<T0, m, E0>& u,
				   const TinyVectorExpression<T1, m, E1>& v, T2& scal) {}

    template<int m, class T0, class T1, class E0>
    static inline void Abs(const TinyVectorExpression<T0, m, E0>& u, T1& scal) {}
    
    template<int m, class T0, class T1, class E0>
    static inline void NormInf(const TinyVectorExpression<T0, m, E0>& u, T1& scal) {}

    template<int m, class T0, class T1, class E0>
    static inline void AbsSquare(const TinyVectorExpression<T0, m, E0>& u, T1& scal) {}
    
    template<int m, class T>
    static inline void Print(const TinyVector<T, m>& V, ostream& out) { }

    template<int m, class T>
    static inline ostream& WriteText(const TinyVector<T, m>& V, ostream& out) { return out; }

    template<int m, class T>
    static inline istream& ReadText(const TinyVector<T, m>& V, istream& in) { return in; }

    template<class T, int m>
    static inline void UpdateMinimum(const TinyVector<T, m>& u, TinyVector<T, m>& v) {}

    template<class T, int m>
    static inline void UpdateMaximum(const TinyVector<T, m>& u, TinyVector<T, m>& v) {}
    
  };

    
  /***********************
   * Blas-like functions *
   ***********************/
  
  // returns (u,v)
  template<class T, int m, class E1, class E2>
  T DotProd(const TinyVectorExpression<T, m, E1> & u,
	    const TinyVectorExpression<T, m, E2> & v);
    
  // returns (u,v) where v is complex and u real
  template<class T, int m, class E1, class E2>
  complex<T> DotProd(const TinyVectorExpression<T, m, E1> & u,
		     const TinyVectorExpression<complex<T>, m, E2> & v);
  
  template<class T, int m, class E1, class E2>
  complex<T> DotProd(const TinyVectorExpression<complex<T>, m, E1> & u,
		     const TinyVectorExpression<T, m, E2> & v);

  // returns (conj(u), v) where u  and v are complex
  template<class T, int m, class E1, class E2>
  complex<T> DotProdConj(const TinyVectorExpression<complex<T>, m, E1> & u,
			 const TinyVectorExpression<complex<T>, m, E2> & v);
  
  // returns || u ||^2 
  template<class T, int m, class E>
  typename ClassComplexType<T>::Treal AbsSquare(const TinyVectorExpression<T, m, E> & u);

  // returns || u ||_1
  template<class T, int m, class E>
  typename ClassComplexType<T>::Treal Norm1(const TinyVectorExpression<T, m, E> & u);

  // returns || u ||_inf
  template<class T, int m, class E>
  typename ClassComplexType<T>::Treal NormInf(const TinyVectorExpression<T, m, E> & u);
  
  template <class T, int m>
  T Distance(const TinyVector<T, m>& u, const TinyVector<T, m> & v);

  // returns || u ||_2
  template<class T, int m, class E>
  typename ClassComplexType<T>::Treal Norm2(const TinyVectorExpression<T, m, E> & u);

  template<class T0, class T1, class T2, int m>
  void Add(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y, TinyVector<T2, m>& z);
  
  template<class T0, class T1, class T2, int m>
  void Add(const T0& alpha, const TinyVector<T1, m>& y, TinyVector<T2, m>& x);
  
  template<class T0, class T, int m>
  void Mlt(const T0& alpha, TinyVector<T, m>& p);


  /***********************
   * Extractions functions *
   ***********************/

  
  template<int p, int m, class T0, class T1>
  void ExtractVector(const Vector<Vector<TinyVector<T0, p> > >& u,
		     int j, int k, int offset, TinyVector<T1, m>& v);

  template<int m, class T0, class T1>
  void ExtractVector(const Vector<Vector<T0> >& u,
		     int j, int offset, TinyVector<T1, m>& v);

  template<int m, class T0, class T1>
  void CopyVector(const TinyVector<T1, m>& v,
		  int j, int offset, Vector<Vector<T0> >& u);

  template<int m, class T0, class T1>
  void ExtractVector(const Vector<T0>& Vh, const IVect& row_num,
		     int nb, TinyVector<T1, m>& tmp);

  template<int m, class T0, class T1>
  void ExtractVector(const TinyVector<T1, m>& tmp,
		     const IVect& row_num, int nb, Vector<T0>& Vh);

  template<int m, class T0, class T1>
  void AddVector(const T0& alpha, const TinyVector<T1, m>& tmp,
		 const IVect& row_num, int nb, Vector<T0>& Vh);

  template<int m, class T0, class T1>
  void ExtractVectorConj(const Vector<T0>& Vh, const IVect& row_num,
			 int nb, TinyVector<T1, m>& tmp);

  template<int m, class T0, class T1>
  void ExtractVectorConj(const TinyVector<T1, m>& tmp,
			 const IVect& row_num, int nb, Vector<T0>& Vh);

  template<int m, class T0, class T1>
  void AddVectorConj(const T0& alpha, const TinyVector<T1, m>& tmp,
		     const IVect& row_num, int nb, Vector<T0>& Vh);


  template<int m, class T0, class T1>
  void CopyVector(const TinyVector<Vector<T0>, m>& uvec,
		  int j, TinyVector<T1, m>& u);

  template<int m, class T0, class T1>
  void CopyVector(const TinyVector<T1, m>& u,
		  int j, TinyVector<Vector<T0>, m>& uvec);
    
  template<int m, class T0, class T1>
  void CopyVector(const TinyVector<T1, m>& u, int j, Vector<T0>& uvec);

  template<int m, class T0, class T1>
  void CopyVector(const Vector<T0>& uvec, int j, TinyVector<T1, m>& u);

  template<int m, class T>
  void CopyVector(const Vector<complex<T> >& Vh, int j, TinyVector<T, m>& vr, TinyVector<T, m>& vi);

  template<int m, class T0, class T1>
  void AddVector(const TinyVector<T1, m>& u, int offset, Vector<T0>& uvec);

  template<int m, class T0, class T1>
  void AddVector(const Vector<T0>& uvec, int offset, TinyVector<T1, m>& u);

  template<int m, class T0, class T1>
  void ExtractVector(const TinyVector<T1, m>& u, int offset, Vector<T0>& uvec);

  template<int m, class T0, class T1>
  void ExtractVector(const Vector<T0>& uvec, int offset, TinyVector<T1, m>& u);

  
  /*******************
   * Other functions *
   *******************/
  

  template<class T0, int m>
  void SymmetrizePointPlane(const TinyVector<T0, m>& u, const TinyVector<T0, m>& normale,
			    const T0& d, TinyVector<T0, m>& v);

  template<class T, int m>
  void UpdateMinimum(const TinyVector<T, m>& u, TinyVector<T, m>& v);

  template<class T, int m>
  void UpdateMaximum(const TinyVector<T, m>& u, TinyVector<T, m>& v);

  template<int m, class T0, class T1>
  void MltVector(const T1& alpha, int j, TinyVector<Vector<T0>, m>& uvec);
  
  template<class T>
  int IntersectionEdges(const TinyVector<T,2>& pointA, const TinyVector<T,2>& pointB,
                        const TinyVector<T,2>& pt1, const TinyVector<T,2>& pt2,
                        TinyVector<T,2>& res, const T& threshold);
  
  template<class T>
  int IntersectionDroites(const TinyVector<T,2>& pointA, const TinyVector<T,2>& pointB,
                          const TinyVector<T,2>& pt1, const TinyVector<T,2>& pt2,
                          TinyVector<T,2>& res, const T& threshold);
  
  template<class T, int p>
  void Mlt(const T& alpha, Vector<TinyVector<T, p> >& X);
  
  template<class T, class Prop, class Storage, int p>
  void Mlt(const SeldonTranspose&, const Matrix<T, Prop, Storage>& A,
	   const Vector<TinyVector<T, p> >& x, Vector<TinyVector<T, p> >& y);
  
  template<class T, class Prop, class Storage, int p>
  void Mlt(const Matrix<T, Prop, Storage>& A,
	   const Vector<TinyVector<T, p> >& x, Vector<TinyVector<T, p> >& y);


  /************************************
   * Functions for 2-D or 3-D vectors *
   ************************************/

  // res = u^v where u,v and res are vectors
  template<class T1, class T2, class T3>
  void TimesProd(const TinyVector<T1, 3> & u, const TinyVector<T2, 3> & v, TinyVector<T3, 3>& res);

  template<class T>
  void ForceZeroVdotN(const TinyVector<T, 2>& normale, TinyVector<T, 2>& v0);

  template<class T>
  void ForceZeroVdotN(const TinyVector<T, 3>& normale, TinyVector<T, 3>& v0);
  
  template<class T>
  bool PointInsideBoundingBox(const TinyVector<T, 2>&, const TinyVector<TinyVector<T, 2>, 2>&);

  template<class T>
  bool PointInsideBoundingBox(const TinyVector<T, 3>&, const TinyVector<TinyVector<T, 3>, 2>&);

  //! generates four points respecting the symmetry of the reference pyramid with one point
  template<class T0, class T1, class T2, class T3>
  void GenerateSymPts(const TinyVector<T0,3> & u, TinyVector<T1,3> & v,
                      TinyVector<T2,3> & w, TinyVector<T3,3> & t);
  
  // res = (u^v) . w = det(u,v,w) where u,v and w are vectors
  template<class T>
  T Determinant(const TinyVector<T,3> & u, const TinyVector<T,3> & v, const TinyVector<T,3>& w);

  //! constructing tangent vectors u1 and u2, from normale to a plane
  template<class T>
  void GetVectorPlane(const TinyVector<T,3>& normale,
                      TinyVector<T, 3>& u1, TinyVector<T, 3>& u2);

  void FillZero(double& X);
  void FillZero(complex<double>& X);
  void FillZero(int& X);
  void FillZero(float& X);
  void FillZero(complex<float>& X);
  void FillZero(bool& X);
  
  template<class T, class Storage, class Allocator>
  void FillZero(Vector<T, Storage, Allocator>& X);

  template<class T, class Prop, class Storage, class Allocator>
  void FillZero(Matrix<T, Prop, Storage, Allocator>& X);
  
  template<class T, class Prop, class Allocator>
  void FillZero(Matrix<T, Prop, RowSparse, Allocator>& X);

  template<class T, int m>
  void FillZero(TinyVector<T, m>& X);
  
  template<class T>
  void FillZero(T& x);

  template<class T>
  void FillZero(complex<T>& x);

  template<class T, int p>
  void SetComplexZero(TinyVector<T, p>& x);
  
  template<class T>
  bool IsComplexZero(const T& x);

  template<class T>
  bool IsComplexZero(const complex<T>& x);

} // end namespace

#define SELDON_FILE_TINY_VECTOR_HXX
#endif


