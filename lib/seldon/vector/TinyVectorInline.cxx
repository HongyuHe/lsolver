#ifndef SELDON_FILE_TINY_VECTOR_INLINE_CXX

#include "TinyVectorExpressionInline.cxx"

namespace Seldon
{
  
  /*******************
   * TinyVector_Base *
   *******************/
  
  
  //! all components are set to 0
  template <class T, int m>
  inline TinyVector_Base<T, m>::TinyVector_Base()
  {
    this->Zero();
  }
  
  
  //! all components are initialized to a
  template <class T, int m>
  inline TinyVector_Base<T, m>::TinyVector_Base(const T& a)
  {
    this->Fill(a);
  }
  
  
  //! first two components are initialized to a and b, others to 0
  template <class T, int m>
  inline TinyVector_Base<T, m>::TinyVector_Base(const T& a, const T& b)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::InitValue(this->GetLeafClass(), b);
    TinyVectorLoop<2*(m>2)>::FillGreater(this->GetLeafClass(), T(0));
  }
  
  
  //! first three components are initialized to a, b and c, others to 0
  template <class T, int m>
  inline TinyVector_Base<T, m>::TinyVector_Base(const T& a, const T& b ,const T& c)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::InitValue(this->GetLeafClass(), b);
    TinyVectorLoop<2*(m>2)>::InitValue(this->GetLeafClass(), c);
    TinyVectorLoop<3*(m>3)>::FillGreater(this->GetLeafClass(), T(0));
  }

  
  //! first four components are initialized to a, b, c, and d, others to 0
  template <class T, int m>
  inline TinyVector_Base<T, m>::TinyVector_Base(const T& a, const T& b ,const T& c, const T& d)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::InitValue(this->GetLeafClass(), b);
    TinyVectorLoop<2*(m>2)>::InitValue(this->GetLeafClass(), c);
    TinyVectorLoop<3*(m>3)>::InitValue(this->GetLeafClass(), d);
    TinyVectorLoop<4*(m>4)>::FillGreater(this->GetLeafClass(), T(0));
  }


  //! first four components are initialized to a, b, c, and d, others to 0
  template <class T, int m>
  inline TinyVector_Base<T, m>::TinyVector_Base(const T& a, const T& b ,const T& c, const T& d, const T& e)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::InitValue(this->GetLeafClass(), b);
    TinyVectorLoop<2*(m>2)>::InitValue(this->GetLeafClass(), c);
    TinyVectorLoop<3*(m>3)>::InitValue(this->GetLeafClass(), d);
    TinyVectorLoop<4*(m>4)>::InitValue(this->GetLeafClass(), e);
    TinyVectorLoop<5*(m>5)>::FillGreater(this->GetLeafClass(), T(0));
  }


  //! sets first component to a
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Init(const T& a) 
  {
#ifdef SELDON_CHECK_BOUNDS
    if (m < 1)
      throw WrongIndex("TinyVector::TinyVector",
		       string("size should be greater than 1,")+
		       "but is equal to" + to_str(m) + ".");
#endif

    this->data_[0] = a;
  }

  
  //! sets two first components of the vector
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Init(const T& a, const T& b) 
  {
#ifdef SELDON_CHECK_BOUNDS
    if (m < 2)
      throw WrongIndex("TinyVector::TinyVector",
		       string("size should be greater than 2,")+
		       "but is equal to" + to_str(m) + ".");
#endif

    this->data_[0] = a;
    this->data_[1] = b;
  }
  
  
  //! sets three first components of the vector
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Init(const T& a, const T& b, const T& c) 
  {
#ifdef SELDON_CHECK_BOUNDS
    if (m < 3)
      throw WrongIndex("TinyVector::TinyVector",
		       string("size should be greater than 3,")+
		       "but is equal to" + to_str(m) + ".");
#endif

    this->data_[0] = a;
    this->data_[1] = b;
    this->data_[2] = c;
  }
  
  
  //! sets four first components of the vector
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Init(const T& a, const T& b, const T& c, const T& d) 
  {
#ifdef SELDON_CHECK_BOUNDS
    if (m < 4)
      throw WrongIndex("TinyVector::TinyVector",
                       string("size should be greater than 4,")
		       +"but is equal to" + to_str(m) + ".");
#endif

    this->data_[0] = a;
    this->data_[1] = b;
    this->data_[2] = c;
    this->data_[3] = d;
  }
  

  //! sets five first components of the vector
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Init(const T& a, const T& b,
					  const T& c, const T& d, const T& e) 
  {
#ifdef SELDON_CHECK_BOUNDS
    if (m < 5)
      throw WrongIndex("TinyVector::TinyVector",
		       string("size should be greater than 5,")+
		       "but is equal to" + to_str(m) + ".");
#endif

    this->data_[0] = a;
    this->data_[1] = b;
    this->data_[2] = c;
    this->data_[3] = d;
    this->data_[4] = e;
  }
  

  //! Returns the number of stored elements.
  template<class T, int m_>
  inline int TinyVector_Base<T, m_>::GetM() const
  {
    return m_;
  }
  
  
  //! Returns the number of stored elements.
  template<class T, int m_>
  inline int TinyVector_Base<T, m_>::GetLength() const
  {
    return m_;
  }
  
  
  //! Returns the number of stored elements.
  template<class T, int m_>
  inline int TinyVector_Base<T, m_>::GetSize() const
  {
    return m_;
  }

  
  //! returns the size used by the object in bytes
  template<class T, int m_>
  inline size_t TinyVector_Base<T, m_>::GetMemorySize() const
  {
    return m_*sizeof(T);
  }
  
  
  //! returns leaf class
  template<class T, int m_>
  inline TinyVector<T, m_>& TinyVector_Base<T, m_>::GetLeafClass()
  {
    return static_cast<TinyVector<T, m_>& >(*this);
  }


  //! returns leaf class
  template<class T, int m_>
  inline const TinyVector<T, m_>& TinyVector_Base<T, m_>::GetLeafClass() const
  {
    return static_cast<const TinyVector<T, m_>& >(*this);
  }
  
    
  //! returns x(i)
  template <class T, int m_>
  inline const T& TinyVector_Base<T, m_>::operator()(int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, m_, "TinVector");
#endif
    
    return this->data_[i];
  }
  
  
  //! returns x(i)
  template <class T, int m_>
  inline T& TinyVector_Base<T, m_>::operator()(int i)
  {
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, m_, "TinVector");
#endif
    
    return this->data_[i];
  }


  //! displays the vector without brackets
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Print() const
  {
    Print(cout);
  }
  
  
  //! displays the vector without brackets
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Print(ostream& out) const
  {
    TinyVectorLoop<m>::Print(this->GetLeafClass(), out);
  }

    
  //! Sets all components to 0
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Zero()
  {
    TinyVectorLoop<m>::Zero(this->GetLeafClass());
  } 
  
  
  //! Sets the vector to [0, 1, 2, 3, ..., n-1]
  template <class T, int m>
  inline void TinyVector_Base<T, m>::Fill()
  {
    TinyVectorLoop<m>::Fill(this->GetLeafClass());
  }
  
  
  //! Sets all components to a given value
  template<class T, int m> template<class T0>
  inline void TinyVector_Base<T, m>::Fill(const T0& x)
  {
    TinyVectorLoop<m>::Fill(this->GetLeafClass(), x);
  }
  
  
  //! sets randomly the components with function rand()
  template <class T, int m>
  inline void TinyVector_Base<T, m>::FillRand()
  {
    TinyVectorLoop<m>::FillRand(this->GetLeafClass());
  }


  //! returns true if all components are 0
  template <class T, int m>
  inline bool TinyVector_Base<T, m>::IsZero() const
  {
    return TinyVectorLoop<m>::IsZero(*this);
  } 

      
  //! returns |u-v|^2 (usual euclidian norm)
  template <class T, int m>
  inline typename ClassComplexType<T>::Treal TinyVector_Base<T, m>::DistanceSquare(const TinyVector<T, m> & p) const
  {
    return AbsSquare(*this - p);
  }
  
  
  //! returns |u-v| (usual euclidian norm)
  template <class T, int m>
  inline typename ClassComplexType<T>::Treal TinyVector_Base<T, m>::Distance(const TinyVector<T, m> & p) const
  {
    return sqrt(DistanceSquare(p));
  }
  

  /*******************
   * TinyVector<T> *
   *******************/
  
  
  //! Sets all components to 0
  template <class T, int m_>
  inline TinyVector<T, m_>::TinyVector() : TinyVector_Base<T, m_>()
  {
  }

    
  //! Sets the first component to a
  template <class T, int m_>
  inline TinyVector<T, m_>::TinyVector(const T& a) : TinyVector_Base<T, m_>(a)
  {
  }
  
  
  //! Sets the two first components
  template <class T, int m_>
  inline TinyVector<T, m_>::TinyVector(const T& a, const T& b)  : TinyVector_Base<T, m_>(a, b) 
  {
  }
  

  //! Sets the three first components
  template <class T, int m_>
  inline TinyVector<T, m_>::TinyVector(const T& a, const T& b, const T& c) : TinyVector_Base<T, m_>(a, b, c)
  {
  }
  

  //! Sets the four first components
  template <class T, int m_>
  inline TinyVector<T, m_>::TinyVector(const T& a, const T& b, const T& c, const T& d)
    :  TinyVector_Base<T, m_>(a, b, c, d)
  {
  }
  

  //! Sets the five first components
  template <class T, int m_>
  inline TinyVector<T, m_>::TinyVector(const T& a, const T& b, const T& c, const T& d, const T& e)
    :  TinyVector_Base<T, m_>(a, b, c, d, e)
  {
  }


  //! Constructor taking an expression with vectors
  template<class T, int m_> template<class E>
  inline TinyVector<T, m_>::TinyVector(const TinyVectorExpression<T, m_, E> & u)
  {
    TinyVectorLoop<m_>::Copy(u, *this);
  }


  //! Operator *this = expression with vectors
  template<class T, int m_> template<class T1, class E>
  inline TinyVector<T, m_>& TinyVector<T, m_>::operator =(const TinyVectorExpression<T1, m_, E> & u)
  {
    TinyVectorLoop<m_>::Copy(u, *this);
    
    return *this;
  }


  //! Operator *this += expression with vectors
  template<class T, int m_> template<class T1, class E>
  inline TinyVector<T, m_>& TinyVector<T, m_>::operator +=(const TinyVectorExpression<T1, m_, E> & u)
  {
    TinyVectorLoop<m_>::AddCopy(u, *this);
    
    return *this;
  }


  //! Operator *this -= expression with vectors
  template<class T, int m_> template<class T1, class E>
  inline TinyVector<T, m_>& TinyVector<T, m_>::operator -=(const TinyVectorExpression<T1, m_, E> & u)
  {
    TinyVectorLoop<m_>::DiffCopy(u, *this);
    
    return *this;
  }

  
  //! Sets all elements of the vector to a given value x
  template<class T, int m_>
  TinyVector<T, m_>& TinyVector<T, m_>::operator =(const T& x)
  {
    this->Fill(x);
    return *this;
  }

  
  //! operation this = this * a where a is scalar
  template <class T, int m> template<class T1>
  inline TinyVector<T, m> & TinyVector<T, m>::operator *=( const T1& a )
  {
    TinyVectorLoop<m>::Mlt(a, *this);
    return *this;
  }
  

  //! operation this = this / a where a is scalar
  template <class T, int m> template<class T1>
  inline TinyVector<T, m> & TinyVector<T, m>::operator /=( const T1& a )
  {
    T one; SetComplexOne(one);
    T a_; SetComplexReal(a, a_);
    TinyVectorLoop<m>::Mlt(one/a_, *this);
    return *this;
  }
 

  /*******************
   * TinyVector<int> *
   *******************/

  
  //! All components are set to 0
  template<int m_>
  inline TinyVector<int, m_>::TinyVector() : TinyVector_Base<int, m_>()
  {
  }

  
  //! Sets the first component of the vector
  template<int m_>
  inline TinyVector<int, m_>::TinyVector(const int& a) : TinyVector_Base<int, m_>(a)
  {
  }

  
  //! Sets the two first components of the vector
  template<int m_>
  inline TinyVector<int, m_>::TinyVector(const int& a, const int& b)  : TinyVector_Base<int, m_>(a, b)
  {
  }

  
  //! Sets the three first components of the vector
  template<int m_>
  inline TinyVector<int, m_>::TinyVector(const int& a, const int& b, const int& c)
    : TinyVector_Base<int, m_>(a, b, c) 
  {
  }

  
  //! Sets the fourth first components of the vector
  template<int m_>
  inline TinyVector<int, m_>::TinyVector(const int& a, const int& b, const int& c, const int& d)
    :  TinyVector_Base<int, m_>(a, b, c, d)
  {
  }

  
  //! Sets the first five components of the vector
  template<int m_>
  inline TinyVector<int, m_>::TinyVector(const int& a, const int& b, const int& c, const int& d, const int& e)
    :  TinyVector_Base<int, m_>(a, b, c, d, e)
  {
  }


  /*************
   * Operators *
   *************/


  //! returns true if u == v, false otherwise
  template<int m>
  inline bool operator==(const TinyVector<int, m>& u, const TinyVector<int, m>& v)
  {
    return TinyVectorLoop<m>::EqualInt(u, v);
  }
  
  
  //! returns false if u == v, true otherwise
  template<int m>
  inline bool operator!=(const TinyVector<int, m> & u, const TinyVector<int, m>& v)
  {
    return !TinyVectorLoop<m>::EqualInt(u, v);
  }
  
    
  //! returns true if u <= v, false otherwise
  template<int m>
  inline bool operator<=(const TinyVector<int, m> & u, const TinyVector<int, m>& v)
  {
    return TinyVectorLoop<m>::LessOrEqualInt(u, v);
  }

  //! returns *this < u
  template<int m>
  inline bool operator<(const TinyVector<int, m> & u, const TinyVector<int, m>& v)
  {
    return TinyVectorLoop<m>::LessInt(u, v);
  }
  
    
  //! returns true if u >= v, false otherwise
  template<int m>
  inline bool operator>=(const TinyVector<int, m> & u, const TinyVector<int, m>& v)
  {
    return !TinyVectorLoop<m>::LessInt(u, v);
  }

    
  //! returns true if u > v, false otherwise
  template<int m>
  inline bool operator>(const TinyVector<int, m> & u, const TinyVector<int, m>& v)
  {
    return !TinyVectorLoop<m>::LessOrEqualInt(u, v);
  }


  //! returns true if u == v, false otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator==(const TinyVectorExpression<T, m, E1>& u,
			 const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorLoop<m>::Equal(u, v);
  }
  
  
  //! returns false if u == v, true otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator!=(const TinyVectorExpression<T, m, E1> & u,
			 const TinyVectorExpression<T, m, E2>& v)
  {
    return !TinyVectorLoop<m>::Equal(u, v);
  }
  
    
  //! returns true if u <= v, false otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator<=(const TinyVectorExpression<T, m, E1> & u,
			 const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorLoop<m>::LessOrEqual(u, v);
  }

  //! returns *this < u
  template<class T, int m, class E1, class E2>
  inline bool operator<(const TinyVectorExpression<T, m, E1> & u,
			const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorLoop<m>::Less(u, v);
  }
  
    
  //! returns true if u >= v, false otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator>=(const TinyVectorExpression<T, m, E1> & u,
			 const TinyVectorExpression<T, m, E2>& v)
  {
    return !TinyVectorLoop<m>::Less(u, v);
  }

    
  //! returns true if u > v, false otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator>(const TinyVectorExpression<T, m, E1> & u,
			const TinyVectorExpression<T, m, E2>& v)
  {
    return !TinyVectorLoop<m>::LessOrEqual(u, v);
  }


  //! displays V
  template <class T, int m>
  inline ostream& operator <<(ostream& out, const TinyVector<T,m>& V)
  {
    out<<"(";
    return TinyVectorLoop<m>::WriteText(V, out);
  }


  //! reads V
  template <class T, int m>
  inline istream& operator >>(istream& in, TinyVector<T,m>& V)
  {
    return TinyVectorLoop<m>::ReadText(V, in);
  }
    

  /******************
   * TinyVectorLoop *
   ******************/


  //! Sets all components to 0
  template<int n> template<int m, class T0>
  inline void TinyVectorLoop<n>::Zero(TinyVector<T0, m>& x)
  {
    FillZero(x(n-1));
    TinyVectorLoop<n-1>::Zero(x);
  }

    
  //! Fills with a in ascendant order (instead of descendant order for Fill)
  template <int n> template<int m, class T0>
  inline void TinyVectorLoop<n>::FillGreater(TinyVector<T0, m>& x, const T0& a)
  {
    x.data_[n] = a;
    TinyVectorLoop<(n+1)*(m>n+1)>::FillGreater(x, a);
  }

  
  //! initializes a single value
  template <int n> template<int m, class T0>
  inline void TinyVectorLoop<n>::InitValue(TinyVector<T0, m>& x, const T0& a)
  {
    x.data_[n] = a;
  }
  

  //! Sets all components to a given value
  template <int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::Fill(TinyVector<T0, m>& x, const T1& a)
  {
    SetComplexReal(a, x.data_[n-1]);
    TinyVectorLoop<n-1>::Fill(x, a);
  }
  
  
  //! Sets the vector to [0, 1, 2, 3, ..., n-1]
  template <int n> template<int m, class T0>
  inline void TinyVectorLoop<n>::Fill(TinyVector<T0, m>& x)
  {
    SetComplexReal(n-1, x.data_[n-1]);
    TinyVectorLoop<n-1>::Fill(x);
  }

  
  //! sets randomly the components with function rand()
  template <int n> template<int m, class T0>
  inline void TinyVectorLoop<n>::FillRand(TinyVector<T0, m>& x)
  {
    SetComplexReal(rand(), x.data_[n-1]);
    TinyVectorLoop<n-1>::FillRand(x);
  }

  
  //! returns true if all components are 0
  template<int n> template<int m, class T0, class E0>
  inline bool TinyVectorLoop<n>::IsZero(const TinyVectorExpression<T0, m, E0>& x)
  {
    if (abs(x(n-1)) > TinyVector<T0, m>::threshold)
      return false;
    
    return TinyVectorLoop<n-1>::IsZero(x);
  }


  //! y = x
  template<int n> template<int m, class T1, class E, class T0>
  inline void TinyVectorLoop<n>::Copy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y)
  {
    y(n-1) = x(n-1);
    TinyVectorLoop<n-1>::Copy(x, y);
  }
  
  
  //! y += x
  template<int n> template<int m, class T1, class E, class T0>
  inline void TinyVectorLoop<n>::AddCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y)
  {
    y(n-1) += x(n-1);
    TinyVectorLoop<n-1>::AddCopy(x, y); 
  }
  

  //! y -= x  
  template<int n> template<int m, class T1, class E, class T0>
  inline void TinyVectorLoop<n>::DiffCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y)
  {
    y(n-1) -= x(n-1);
    TinyVectorLoop<n-1>::DiffCopy(x, y); 
  }
  

  //! u(:) = uvec(:)(j)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  CopyVector(const TinyVector<Vector<T0>, m>& uvec, int j, TinyVector<T1, m>& u)
  {
    u(n-1) = uvec(n-1)(j);
    TinyVectorLoop<n-1>::CopyVector(uvec, j, u);
  }

  
  //! uvec(:)(j) = u(:)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  CopyVector(const TinyVector<T1, m>& u, int j, TinyVector<Vector<T0>, m>& uvec)
  {
    uvec(n-1)(j) = u(n-1);
    TinyVectorLoop<n-1>::CopyVector(u, j, uvec);
  }
  
  
  //! tmp(:) = Vh(m*j + :)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  CopyVector(const Vector<T0>& Vh, int j, TinyVector<T1, m>& tmp)
  {
    tmp(n-1) = Vh(m*j + n-1);
    TinyVectorLoop<n-1>::CopyVector(Vh, j, tmp);
  }


  //! vr(:) = real(Vh(m*j + :)) and vi(:) = imag(Vh(m*j + :))
  template<int n> template<int m, class T>
  inline void TinyVectorLoop<n>::
  CopyVector(const Vector<complex<T> >& Vh, int j, TinyVector<T, m>& vr, TinyVector<T, m>& vi)
  {
    vr(n-1) = realpart(Vh(m*j + n-1));
    vi(n-1) = imagpart(Vh(m*j + n-1));
    TinyVectorLoop<n-1>::CopyVector(Vh, j, vr, vi);
  }
  
  
  //! Vh(m*j + :) = tmp(:)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  CopyVector(const TinyVector<T1, m>& tmp, int j, Vector<T0>& Vh)
  {
    Vh(m*j + n-1) = tmp(n-1);
    TinyVectorLoop<n-1>::CopyVector(tmp, j, Vh);
  }
  

  //! tmp(:) += Vh(offset + :)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  AddVector(const Vector<T0>& Vh, int offset, TinyVector<T1, m>& tmp)
  {
    tmp(n-1) += Vh(offset + n-1);
    TinyVectorLoop<n-1>::AddVector(Vh, offset, tmp);
  }
  
  
  //! Vh(offset + :) += tmp(:)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  AddVector(const TinyVector<T1, m>& tmp, int offset, Vector<T0>& Vh)
  {
    Vh(offset + n-1) += tmp(n-1);
    TinyVectorLoop<n-1>::AddVector(tmp, offset, Vh);
  }

  
  //! tmp(:) = Vh(offset + :)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  ExtractVector(const Vector<T0>& Vh, int offset, TinyVector<T1, m>& tmp)
  {
    tmp(n-1) = Vh(offset + n-1);
    TinyVectorLoop<n-1>::ExtractVector(Vh, offset, tmp);
  }
  
  
  //! Vh(offset + :) = tmp(:)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  ExtractVector(const TinyVector<T1, m>& tmp, int offset, Vector<T0>& Vh)
  {
    Vh(offset + n-1) = tmp(n-1);
    TinyVectorLoop<n-1>::ExtractVector(tmp, offset, Vh);
  }


  //! v(:) = u(offset+:)(j)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::ExtractVector(const Vector<Vector<T0> >& u,
		  int j, int offset, TinyVector<T1, m>& v)
  {
    v(n-1) = u(offset+n-1)(j);
    TinyVectorLoop<n-1>::ExtractVector(u, j, offset, v);
  }


  //! u(offset+:)(j) = v(:)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::CopyVector(const TinyVector<T1, m>& v,
	       int j, int offset, Vector<Vector<T0> >& u)
  {
    u(offset+n-1)(j) = v(n-1);
    TinyVectorLoop<n-1>::CopyVector(v, j, offset, u);
  }


  //! v(:) = u(offset+:)(j)(k)
  template<int n> template<int p, int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::ExtractVector(const Vector<Vector<TinyVector<T0, p> > >& u,
		  int j, int k, int offset, TinyVector<T1, m>& v)
  {
    v(n-1) = u(offset+n-1)(j)(k);
    TinyVectorLoop<n-1>::ExtractVector(u, j, k, offset, v);
  }
  
  
  //! tmp = Vh(row_num(nb+:))
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::ExtractVector(const Vector<T0>& Vh, const IVect& row_num,
                  int nb, TinyVector<T1, m>& tmp)
  {
    tmp(n-1) = Vh(row_num(nb + n-1));
    TinyVectorLoop<n-1>::ExtractVector(Vh, row_num, nb, tmp);
  }
  
  
  //! Vh(row_num(nb+:)) = tmp
  template<int n>  template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::ExtractVector(const TinyVector<T1, m>& tmp,
                  const IVect& row_num, int nb, Vector<T0>& Vh)
  {
    Vh(row_num(nb + n-1)) = tmp(n-1);
    TinyVectorLoop<n-1>::ExtractVector(tmp, row_num, nb, Vh);
  }
  

  //! x_dest(offset_dest+:) = x_src(offset_source+:)
  template<int n> template<class T0, class T1>
  inline void TinyVectorLoop<n>
  ::ExtractVector(const Vector<T0>& x_src, int offset_source,
		  Vector<T1>& x_dest, int offset_dest)
  {
    x_dest(offset_dest + n-1) = x_src(offset_source + n-1);
    TinyVectorLoop<n-1>::ExtractVector(x_src, offset_source, x_dest, offset_dest);
  }
  
  
  //! Vh(row_num(nb+:)) += alpha * tmp
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::AddVector(const T0& alpha, const TinyVector<T1, m>& tmp,
              const IVect& row_num, int nb, Vector<T0>& Vh)
  {
    Vh(row_num(nb + n-1)) += alpha*tmp(n-1);
    TinyVectorLoop<n-1>::AddVector(alpha, tmp, row_num, nb, Vh);
  }


  //! tmp = conj(Vh(row_num(nb+:)))
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::ExtractVectorConj(const Vector<T0>& Vh, const IVect& row_num,
		      int nb, TinyVector<T1, m>& tmp)
  {
    tmp(n-1) = conjugate(Vh(row_num(nb + n-1)));
    TinyVectorLoop<n-1>::ExtractVectorConj(Vh, row_num, nb, tmp);
  }

  
  //! Vh(row_num(nb+:)) = conj(tmp)
  template<int n>  template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::ExtractVectorConj(const TinyVector<T1, m>& tmp,
                      const IVect& row_num, int nb, Vector<T0>& Vh)
  {
    Vh(row_num(nb + n-1)) = conjugate(tmp(n-1));
    TinyVectorLoop<n-1>::ExtractVectorConj(tmp, row_num, nb, Vh);
  }
  
  
  //! Vh(row_num(nb+:)) += alpha * conj(tmp)
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>
  ::AddVectorConj(const T0& alpha, const TinyVector<T1, m>& tmp,
                  const IVect& row_num, int nb, Vector<T0>& Vh)
  {
    Vh(row_num(nb + n-1)) += alpha*conjugate(tmp(n-1));
    TinyVectorLoop<n-1>::AddVectorConj(alpha, tmp, row_num, nb, Vh);
  }

       
  //! u(:)(j) *= alpha
  template<int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::
  MltVector(const T1& alpha, int j, TinyVector<Vector<T0>, m>& uvec)
  {
    uvec(n-1)(j) *= alpha;
    TinyVectorLoop<n-1>::MltVector(alpha, j, uvec);
  }
  
  
  //! returns x == y
  template <int n> template<int m, class T0, class T1, class E0, class E1>
  inline bool TinyVectorLoop<n>::
  Equal(const TinyVectorExpression<T0, m, E0>& x,
	const TinyVectorExpression<T1, m, E1>& y)
  {
    if (abs(x(m-n)-y(m-n)) > TinyVector<T0, m>::threshold)
      return false;
    
    return TinyVectorLoop<n-1>::Equal(x, y);
  }

  
  //! returns x <= y
  template <int n> template<int m, class T0, class T1, class E0, class E1>
  inline bool TinyVectorLoop<n>::LessOrEqual(const TinyVectorExpression<T0, m, E0>& x,
                                             const TinyVectorExpression<T1, m, E1>& y)
  {
    if (x(m-n) < (y(m-n)-TinyVector<T0, m>::threshold) )
      return true;
    else if ( x(m-n) > (y(m-n) + TinyVector<T0, m>::threshold))
      return false;
    
    return TinyVectorLoop<n-1>::LessOrEqual(x, y);
  }
  

  //! returns x < y
  template <int n> template<int m, class T0, class T1, class E0, class E1>
  inline bool TinyVectorLoop<n>::Less(const TinyVectorExpression<T0, m, E0>& x,
				      const TinyVectorExpression<T1, m, E1>& y)
  {
    if (x(m-n) < (y(m-n)-TinyVector<T0, m>::threshold) )
      return true;
    else if (x(m-n) > (y(m-n) + TinyVector<T0, m>::threshold))
      return false;

    return TinyVectorLoop<n-1>::Less(x, y);
  }
  

  //! returns x == y
  template <int n> template<int m, class T0, class T1>
  inline bool TinyVectorLoop<n>::
  EqualInt(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y)
  {
    if (x.data_[m-n] != y.data_[m-n])
      return false;
    
    return TinyVectorLoop<n-1>::EqualInt(x, y);
  }

  
  //! returns x <= y
  template <int n> template<int m, class T0, class T1>
  inline bool TinyVectorLoop<n>::LessOrEqualInt(const TinyVector<T0, m>& x,
						const TinyVector<T1, m>& y)
  {
    if (x.data_[m-n] < y.data_[m-n] )
      return true;
    else if ( x.data_[m-n] > y.data_[m-n])
      return false;
    
    return TinyVectorLoop<n-1>::LessOrEqualInt(x, y);
  }
  

  //! returns x < y
  template <int n> template<int m, class T0, class T1>
  inline bool TinyVectorLoop<n>::LessInt(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y)
  {
    if (x.data_[m-n] < y.data_[m-n])
      return true;
    else if (x.data_[m-n] > y.data_[m-n])
      return false;

    return TinyVectorLoop<n-1>::LessInt(x, y);
  }
  

  //! Multiplies x by alpha
  template <int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::Mlt(const T1& alpha, TinyVector<T0, m>& x)
  {
    x(n-1) *= alpha;
    TinyVectorLoop<n-1>::Mlt(alpha, x);
  }
  

  //! scal = u.v
  template <int n> template<int m, class T0, class T1, class T2, class E0, class E1>
  inline void TinyVectorLoop<n>::
  DotProd(const TinyVectorExpression<T0, m, E0>& u,
	  const TinyVectorExpression<T1, m, E1>& v, T2& scal)
  {
    scal += u(n)*v(n);  
    TinyVectorLoop<n-1>::DotProd(u, v, scal);
  }
  
  
  //! returns (conj(u), v)
  template <int n> template<int m, class T0, class T1, class T2, class E0, class E1>
  inline void TinyVectorLoop<n>::
  DotProdConj(const TinyVectorExpression<T0, m, E0>& u,
	      const TinyVectorExpression<T1, m, E1>& v, T2& scal)
  {
    scal += conjugate(u(n))*v(n);  
    TinyVectorLoop<n-1>::DotProdConj(u, v, scal);
  }
  

  //! computes scal = \sum |u_i|^2
  template <int n> template<int m, class T0, class T1, class E0>
  inline void TinyVectorLoop<n>::Abs(const TinyVectorExpression<T0, m, E0>& u, T1& scal)
  {
    scal += abs(u(n));
    TinyVectorLoop<n-1>::Abs(u, scal);
  }


  //! computes scal = max |u_i|^2
  template <int n> template<int m, class T0, class T1, class E0>
  inline void TinyVectorLoop<n>::NormInf(const TinyVectorExpression<T0, m, E0>& u, T1& scal)
  {
    scal = max(scal, abs(u(n)));
    TinyVectorLoop<n-1>::NormInf(u, scal);
  }
  
  
  //! computes scal = \sum |u_i|^2
  template <int n> template<int m, class T0, class T1, class E0>
  inline void TinyVectorLoop<n>::AbsSquare(const TinyVectorExpression<T0, m, E0>& u, T1& scal)
  {
    scal += absSquare(u(n));
    TinyVectorLoop<n-1>::AbsSquare(u, scal);
  }


  //! displays V without brackets
  template<int n> template<int m, class T>
  inline void TinyVectorLoop<n>::Print(const TinyVector<T, m>& V, ostream& out)
  {
    if (n==1)
      out << V(m-n);  
    else
      out<<V(m-n)<<" ";  
    
    return TinyVectorLoop<n-1>::Print(V, out);
  }
  
  
  //! displays V
  template<int n> template<int m, class T>
  inline ostream& TinyVectorLoop<n>::WriteText(const TinyVector<T, m>& V, ostream& out)
  {
    if (n==1)
      out<<V(m-n)<<')';  
    else
      out<<V(m-n)<<", ";  
    
    return TinyVectorLoop<n-1>::WriteText(V, out);
  }
  

  //! reads V
  template<int n> template<int m, class T>
  inline istream& TinyVectorLoop<n>::ReadText(TinyVector<T, m>& V, istream& in)
  {
    in>>V(m-n);  
    return TinyVectorLoop<n-1>::ReadText(V, in);
  }
    
  
  //! sets v = min(v, u)
  template<int n> template<class T, int m>
  inline void TinyVectorLoop<n>::UpdateMinimum(const TinyVector<T, m>& u, TinyVector<T, m>& v)
  {
    v.data_[n-1] = min(v.data_[n-1], u.data_[n-1]);
    TinyVectorLoop<n-1>::UpdateMinimum(u, v);
  }

  
  //! sets v = max(v, u)
  template<int n> template<class T, int m>
  inline void TinyVectorLoop<n>::UpdateMaximum(const TinyVector<T, m>& u, TinyVector<T, m>& v)
  {
    v.data_[n-1] = max(v.data_[n-1], u.data_[n-1]);
    TinyVectorLoop<n-1>::UpdateMaximum(u, v);
  }


  /******************
   * Blas Functions *
   ******************/
  

  //! returns scalar product (u,v)
  template<class T, int m, class E1, class E2>
  inline T DotProd(const TinyVectorExpression<T, m, E1> & u,
		   const TinyVectorExpression<T, m, E2> & v)
  {
    T scal = u(0)*v(0);
    TinyVectorLoop<m-1>::DotProd(u, v, scal);
    
    return scal;
  }

  
  //! returns scalar product (u,v)
  template<class T, int m, class E1, class E2>
  inline complex<T> DotProd(const TinyVectorExpression<T, m, E1> & u,
			    const TinyVectorExpression<complex<T>, m, E2> & v)
  {
    complex<T> scal = u(0)*v(0);
    TinyVectorLoop<m-1>::DotProd(u, v, scal);
    
    return scal;
  }
  
  
  //! returns (u,v) where u is complex and v real
  template<class T, int m, class E1, class E2>
  inline complex<T> DotProd(const TinyVectorExpression<complex<T>, m, E1> & u,
			    const TinyVectorExpression<T, m, E2> & v)
  {
    complex<T> scal = u(0)*v(0);
    TinyVectorLoop<m-1>::DotProd(u, v, scal);
    
    return scal;
  }
  
  
  //! returns (conj(u), v)
  template<class T, int m, class E1, class E2>
  inline complex<T> DotProdConj(const TinyVectorExpression<complex<T>, m, E1> & u,
				const TinyVectorExpression<complex<T>, m, E2> & v)
  {
    complex<T> scal = conjugate(u(0))*v(0);
    TinyVectorLoop<m-1>::DotProdConj(u, v, scal);
    
    return scal;
  }

  
  //! returns || u ||^2
  template<class T, int m, class E>
  inline typename ClassComplexType<T>::Treal AbsSquare(const TinyVectorExpression<T, m, E> & u)
  {
    typename ClassComplexType<T>::Treal scal = absSquare(u(0));
    TinyVectorLoop<m-1>::AbsSquare(u, scal);
    
    return scal;
  }


  //! returns || u ||_1
  template<class T, int m, class E>
  inline typename ClassComplexType<T>::Treal Norm1(const TinyVectorExpression<T, m, E> & u)
  {
    typename ClassComplexType<T>::Treal scal = abs(u(0));
    TinyVectorLoop<m-1>::Abs(u, scal);
    
    return scal;
  }


  //! returns || u ||_\infty
  template<class T, int m, class E>
  inline typename ClassComplexType<T>::Treal NormInf(const TinyVectorExpression<T, m, E> & u)
  {
    typename ClassComplexType<T>::Treal scal = abs(u(0));
    TinyVectorLoop<m-1>::NormInf(u, scal);
    
    return scal;
  }

  
  //! returns || u-v ||
  template <class T, int m>
  inline T Distance(const TinyVector<T, m>& u, const TinyVector<T, m> & v)
  {
    return u.Distance(v);
  }

  
  //! returns || v ||_2
  template<class T, int m, class E>
  inline typename ClassComplexType<T>::Treal Norm2(const TinyVectorExpression<T, m, E> & p)
  {
    return sqrt(AbsSquare(p));
  }

  
  //! z = x + y
  template<class T0, class T1, class T2, int m>
  void Add(const TinyVector<T0, m>& x, const TinyVector<T1, m>& y, TinyVector<T2, m>& z)
  {
    z = x + y;
  }

  
  //! x = x + alpha*y
  template<class T0, class T1, class T2, int m>
  inline void Add(const T0& alpha, const TinyVector<T1, m>& y, TinyVector<T2, m>& x)
  {
    x += alpha*y;
  }


  //! x = alpha * x
  template<class T0, class T, int m>
  inline void Mlt(const T0& alpha, TinyVector<T, m>& p)
  {
    p *= alpha;
  }


  /************************
   * Extraction Functions *
   ************************/
 
    
  //! v(:) = u(offset+:)(j)(k)
  template<int p, int m, class T0, class T1>
  inline void ExtractVector(const Vector<Vector<TinyVector<T0, p> > >& u,
			    int j, int k, int offset, TinyVector<T1, m>& v)
  {
    TinyVectorLoop<m>::ExtractVector(u, j, k, offset, v);
  }


  //! v(:) = u(offset+:)(j)
  template<int m, class T0, class T1>
  inline void ExtractVector(const Vector<Vector<T0> >& u,
			    int j, int offset, TinyVector<T1, m>& v)
  {
    TinyVectorLoop<m>::ExtractVector(u, j, offset, v);
  }


  //! u(offset+:)(j) = v(:)
  template<int m, class T0, class T1>
  inline void CopyVector(const TinyVector<T1, m>& v,
			 int j, int offset, Vector<Vector<T0> >& u)
  {
    TinyVectorLoop<m>::CopyVector(v, j, offset, u);
  }
  
  
  //! tmp = Vh(row_num(nb+:))  
  template<int m, class T0, class T1>
  inline void ExtractVector(const Vector<T0>& Vh, const IVect& row_num,
                            int nb, TinyVector<T1, m>& tmp)
  {
    TinyVectorLoop<m>::ExtractVector(Vh, row_num, nb, tmp);
  }
  

  //! Vh(row_num(nb+:)) = tmp
  template<int m, class T0, class T1>
  inline void ExtractVector(const TinyVector<T1, m>& tmp,
                            const IVect& row_num, int nb, Vector<T0>& Vh)
  {
    TinyVectorLoop<m>::ExtractVector(tmp, row_num, nb, Vh);
  }


  //! Vh(row_num(nb+:)) += alpha * tmp
  template<int m, class T0, class T1>
  inline void AddVector(const T0& alpha, const TinyVector<T1, m>& tmp,
                        const IVect& row_num, int nb, Vector<T0>& Vh)
  {
    TinyVectorLoop<m>::AddVector(alpha, tmp, row_num, nb, Vh);
  }


  //! tmp = conj(Vh(row_num(nb+:)))
  template<int m, class T0, class T1>
  inline void ExtractVectorConj(const Vector<T0>& Vh, const IVect& row_num,
                                int nb, TinyVector<T1, m>& tmp)
  {
    TinyVectorLoop<m>::ExtractVectorConj(Vh, row_num, nb, tmp);
  }

  
  //! Vh(row_num(nb+:)) = conj(tmp)
  template<int m, class T0, class T1>
  inline void ExtractVectorConj(const TinyVector<T1, m>& tmp,
                                const IVect& row_num, int nb, Vector<T0>& Vh)
  {
    TinyVectorLoop<m>::ExtractVectorConj(tmp, row_num, nb, Vh);
  }


  //! Vh(row_num(nb+:)) += alpha * conj(tmp)
  template<int m, class T0, class T1>
  inline void AddVectorConj(const T0& alpha, const TinyVector<T1, m>& tmp,
                            const IVect& row_num, int nb, Vector<T0>& Vh)
  {
    TinyVectorLoop<m>::AddVectorConj(alpha, tmp, row_num, nb, Vh);
  }
  

  //! u(:) = uvec(:)(j)
  template<int m, class T0, class T1>
  inline void CopyVector(const TinyVector<Vector<T0>, m>& uvec,
                         int j, TinyVector<T1, m>& u)
  {
    TinyVectorLoop<m>::CopyVector(uvec, j, u);
  }
  
  
  //! uvec(:)(j) = u(:)
  template<int m, class T0, class T1>
  inline void CopyVector(const TinyVector<T1, m>& u,
                         int j, TinyVector<Vector<T0>, m>& uvec)
  {
    TinyVectorLoop<m>::CopyVector(u, j, uvec);
  }

  
  //! uvec(m*j + :) = u(:)
  template<int m, class T0, class T1>
  inline void CopyVector(const TinyVector<T1, m>& u, int j, Vector<T0>& uvec)
  {
    TinyVectorLoop<m>::ExtractVector(u, m*j, uvec);
  }
  
  
  //! u(:) = uvec(m*j + :)
  template<int m, class T0, class T1>
  inline void CopyVector(const Vector<T0>& uvec, int j, TinyVector<T1, m>& u)
  {
    TinyVectorLoop<m>::ExtractVector(uvec, m*j, u);
  }
  

  //! vr(:) = real(Vh(m*j + :)) and vi(:) = imag(Vh(m*j + :))
  template<int m, class T>
  inline void CopyVector(const Vector<complex<T> >& Vh, int j,
			 TinyVector<T, m>& vr, TinyVector<T, m>& vi)
  {
    TinyVectorLoop<m>::CopyVector(Vh, j, vr, vi);
  }
  
  
  //! uvec(offset + :) += u(:)
  template<int m, class T0, class T1>
  inline void AddVector(const TinyVector<T1, m>& u, int offset, Vector<T0>& uvec)
  {
    TinyVectorLoop<m>::AddVector(u, offset, uvec);
  }
  
  
  //! u(:) += uvec(offset + :)
  template<int m, class T0, class T1>
  inline void AddVector(const Vector<T0>& uvec, int offset, TinyVector<T1, m>& u)
  {
    TinyVectorLoop<m>::AddVector(uvec, offset, u);
  }

  
  //! uvec(offset + :) = u(:)
  template<int m, class T0, class T1>
  inline void ExtractVector(const TinyVector<T1, m>& u, int offset, Vector<T0>& uvec)
  {
    TinyVectorLoop<m>::ExtractVector(u, offset, uvec);
  }
  
  
  //! u(:) = uvec(offset + :)
  template<int m, class T0, class T1>
  inline void ExtractVector(const Vector<T0>& uvec, int offset, TinyVector<T1, m>& u)
  {
    TinyVectorLoop<m>::ExtractVector(uvec, offset, u);
  }


  /*******************
   * Other Functions *
   *******************/
 

  //! Applies the symmetry with respect to the hyperplane a x + b y + c z + .. + d = 0
  //! v = S(u) where S is this symmetry, normale contains (a, b, c, ...)
  template<class T0, int m>
  inline void SymmetrizePointPlane(const TinyVector<T0, m>& u, const TinyVector<T0, m>& normale,
				   const T0& d, TinyVector<T0, m>& v)
  {
    T0 scal = 2*(-d - DotProd(u, normale));
    v = u;
    Add(scal, normale, v);
  }
  
  
  //! sets v = min(v, u) elementwise
  template<class T, int m>
  inline void UpdateMinimum(const TinyVector<T, m>& u, TinyVector<T, m>& v)
  {
    TinyVectorLoop<m>::UpdateMinimum(u, v);
  }


  //! sets v = max(v, u) elementwise
  template<class T, int m>
  inline void UpdateMaximum(const TinyVector<T, m>& u, TinyVector<T, m>& v)
  {
    TinyVectorLoop<m>::UpdateMaximum(u, v);
  }

  
  //! uvec(:)(j) *= alpha
  template<int m, class T0, class T1>
  inline void MltVector(const T1& alpha, int j, TinyVector<Vector<T0>, m>& uvec)
  {
    TinyVectorLoop<m>::MltVector(alpha, j, uvec);
  }


  template<class T, int p>
  inline void Mlt(const T& alpha, Vector<TinyVector<T, p> >& X)
  {
    Vector<T> Xvec;
    Xvec.SetData(X.GetM()*p, reinterpret_cast<T*>(X.GetData()));
    Mlt(alpha, Xvec);
    Xvec.Nullify();
  }

  
  template<class T, class Prop, class Storage, int p>
  inline void Mlt(const SeldonTranspose& trans, const Matrix<T, Prop, Storage>& A,
		  const Vector<TinyVector<T, p> >& x, Vector<TinyVector<T, p> >& y)
  {
    if (trans.NoTrans())
      MltVector(SeldonNoTrans, A, x, y);
    else
      MltVector(SeldonTrans, A, x, y);
  }
  

  template<class T, class Prop, class Storage, int p>
  inline void Mlt(const Matrix<T, Prop, Storage>& A,
		  const Vector<TinyVector<T, p> >& x, Vector<TinyVector<T, p> >& y)
  {
    MltVector(A, x, y);
  }

  
  /*****************************
   * Functions for 3-D vectors *
   *****************************/
  
  
  //! cross product res = u^v where u,v and res are vectors
  template<class T1, class T2, class T3>
  inline void TimesProd(const TinyVector<T1, 3> & u,
                        const TinyVector<T2, 3> & v, TinyVector<T3, 3>& res )
  {
    res(0) = u(1)*v(2) - u(2)*v(1);
    res(1) = u(2)*v(0) - u(0)*v(2);
    res(2) = u(0)*v(1) - u(1)*v(0);
  }

  
  //! cross product for 2-D vector
  template<class T1, class T2>
  inline T1 TimesProd(const TinyVector<T1, 2> & u, const TinyVector<T2, 2> & v)
  {
    return u(0)*v(1) - u(1)*v(0);
  }
  
  
  //! modifies v0 such that v0 \cdot normale = 0
  template<class T>
  inline void ForceZeroVdotN(const TinyVector<T, 2>& normale, TinyVector<T, 2>& v0)
  {
    if (abs(normale(0)) > abs(normale(1)))
      v0.Init(-normale(1)*v0(1)/normale(0), v0(1));
    else
      v0.Init(v0(0), -normale(0)*v0(0)/normale(1));
  }


  //! modifies v0 such that v0 \cdot normale = 0
  template<class T>
  inline void ForceZeroVdotN(const TinyVector<T, 3>& normale, TinyVector<T, 3>& v0)
  {
    if (abs(normale(0)) > max(abs(normale(1)), abs(normale(2))))
      v0.Init(-(normale(2)*v0(2)+normale(1)*v0(1))/normale(0), v0(1), v0(2));
    else if (abs(normale(1)) > max(abs(normale(0)), abs(normale(2))))
      v0.Init(v0(0), -(normale(0)*v0(0)+normale(2)*v0(2))/normale(1), v0(2));
    else
      v0.Init(v0(0), v0(1), -(normale(0)*v0(0)+normale(1)*v0(1))/normale(2));
  }


  //! returns true if x is inside the bounding box s
  template<class T>
  inline bool PointInsideBoundingBox(const TinyVector<T, 2>& x, const TinyVector<TinyVector<T, 2>, 2>& s)
  {
    //DISP(x); DISP(s); DISP((TinyVector<T, 2>::threshold));
    if ((x(0) >= s(0)(0)-TinyVector<T, 2>::threshold) && (x(0) <= s(1)(0)+TinyVector<T, 2>::threshold)
	&& (x(1) >= s(0)(1)-TinyVector<T, 2>::threshold) && (x(1) <= s(1)(1)+TinyVector<T, 2>::threshold))
      return true;
    
    return false;
  }
  

  //! returns true if x is inside the bounding box s
  template<class T>
  inline bool PointInsideBoundingBox(const TinyVector<T, 3>& x, const TinyVector<TinyVector<T, 3>, 2>& s)
  {
    if ((x(0) >= s(0)(0)-TinyVector<T, 3>::threshold) && (x(0) <= s(1)(0)+TinyVector<T, 3>::threshold)
	&& (x(1) >= s(0)(1)-TinyVector<T, 3>::threshold) && (x(1) <= s(1)(1)+TinyVector<T, 3>::threshold)
	&& (x(2) >= s(0)(2)-TinyVector<T, 3>::threshold) && (x(2) <= s(1)(2)+TinyVector<T, 3>::threshold) )
      return true;
    
    return false;
  }

  
  //! generates four points respecting the symmetry of the reference pyramid with one point
  template<class T0, class T1, class T2, class T3>
  inline void GenerateSymPts(const TinyVector<T0, 3> & u, TinyVector<T1, 3> & v,
                             TinyVector<T2, 3> & w, TinyVector<T3, 3> & t)
  {
    v(0) = -u(1); v(1) = u(0); v(2) = u(2);
    w(0) = -v(1); w(1) = v(0); w(2) = v(2);
    t(0) = -w(1); t(1) = w(0); t(2) = w(2);
  }
  
  
  //! res = (u^v) . w = det(u,v,w) where u,v and w are vectors
  template<class T>
  inline T Determinant(const TinyVector<T, 3> & u,
                       const TinyVector<T, 3> & v, const TinyVector<T,3>& w )
  {
    T res = u(0)*(v(1)*w(2)-v(2)*w(1)) - u(1)*(v(0)*w(2)-v(2)*w(0)) + u(2)*(v(0)*w(1)-v(1)*w(0));
    return res;
  }
  
  
  /************
   * FillZero *
   ************/
  
  //! Sets x to 0
  inline void FillZero(double& X)
  {
    X = 0;
  }

  //! Sets x to 0
  inline void FillZero(complex<double>& X)
  {
    X = 0;
  }
  
  //! Sets x to 0
  inline void FillZero(int& X)
  {
    X = 0;
  }

  //! Sets x to 0
  inline void FillZero(float& X)
  {
    X = 0;
  }

  //! Sets x to 0
  inline void FillZero(complex<float>& X)
  {
    X = 0;
  }
  
  //! Sets x to 0
  inline void FillZero(bool& X)
  {
    X = false;
  }
  
  //! fills all the elements in vector X with 0
  template<class T, class Storage, class Allocator>
  inline void FillZero(Vector<T, Storage, Allocator>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      FillZero(X(i));
  }

  //! fills all the elements in vector X with 0
  template<class T, class Prop, class Storage, class Allocator>
  inline void FillZero(Matrix<T, Prop, Storage, Allocator>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      for (int j = 0; j < X.GetN(); j++)
	FillZero(X(i, j));
  }

  //! Sets x to 0
  template<class T, class Prop, class Allocator>
  inline void FillZero(Matrix<T, Prop, RowSparse, Allocator>& X)
  {
  }

  //! Sets x to 0
  template<class T, class Prop, class Allocator>
  inline void FillZero(Matrix<T, Prop, ArrayRowSparse, Allocator>& X)
  {
  }
  
  //! Sets x to 0
  template<class T, int m>
  inline void FillZero(TinyVector<T, m>& X)
  {
    X.Zero();
  }

  //! Sets x to 0
  template<class T>
  inline void FillZero(T& x)
  {
    x = T(0);
  }

  //! Sets x to 0
  template<class T>
  inline void FillZero(complex<T>& x)
  {
    x = complex<T>(0, 0);
  }

  template<class T, int p>
  inline void SetComplexZero(TinyVector<T, p>& x)
  {
    T zero; SetComplexZero(zero);
    x.Fill(zero);
  }
  
  template<class T>
  inline bool IsComplexZero(const T& x)
  {
    if (x == T(0))
      return true;
    
    return false;
  }

  template<class T>
  inline bool IsComplexZero(const complex<T>& x)
  {
    if (x == complex<T>(0, 0))
      return true;
    
    return false;
  }

} // end namespace

#define SELDON_FILE_TINY_VECTOR_INLINE_CXX
#endif
