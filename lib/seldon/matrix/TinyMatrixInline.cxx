#ifndef MONTJOIE_FILE_TINY_MATRIX_INLINE_CXX

#include "TinyMatrix.hxx"

#include "TinyMatrixExpressionInline.cxx"

namespace Seldon
{
  
  /***********************
   * TinyMatrix<General> *
   ***********************/
  
  
  //! Default constructor : all elements of the matrix are set to 0
  template<class T, int m, int n>
  inline TinyMatrix<T, General, m, n>::TinyMatrix()
  {
    this->Zero();
  }
  
  
  //! Default constructor : all elements of the matrix are set to 0
  template<class T, int m, int n>
  inline TinyMatrix<T, General, m, n>::TinyMatrix(int i)
  {
    this->Zero();
  }


  //! Default constructor : all elements of the matrix are set to 0
  template<class T, int m, int n>
  inline TinyMatrix<T, General, m, n>::TinyMatrix(int i, int j)
  {
    this->Zero();
  }


  //! Constructor taking columns as arguments
  template<class T, int m, int n> template<class T0>
  inline TinyMatrix<T, General, m, n>::
  TinyMatrix(const TinyVector<TinyVector<T0, m>, n>& A)
  {
    TinyMatrixLoop<n>::Init(A, *this);
  }
  

  //! Constructor taking an expression
  template<class T, int m, int n> template<class E>
  inline TinyMatrix<T, General, m, n>
  ::TinyMatrix(const TinyMatrixExpression<T, m, n, E>& A)
  {
    TinyMatrixLoop<m*n>::Copy(A, *this);
  }
  
  
  //! returns the number of rows
  template<class T, int m, int n>
  inline int TinyMatrix<T, General, m, n>::GetM()
  {
    return m;
  }
  
  
  //! returns the number of columns
  template<class T, int m, int n>
  inline int TinyMatrix<T, General, m, n>::GetN()
  {
    return n;
  }
  
  
  //! returns the number of stored elements
  template<class T, int m, int n>
  inline int TinyMatrix<T, General, m, n>::GetSize()
  {
    return m*n;
  }
  

  //! returns the pointer to the C-array
  template<class T, int m, int n>
  inline T* TinyMatrix<T, General, m, n>::GetData()
  {
    return data_;
  }  


  //! *this = *this*a
  template<class T, int m, int n> template<class T1>
  inline TinyMatrix<T, General, m, n> & TinyMatrix<T, General, m, n>::
  operator *=(const T1 & a )
  {
    TinyMatrixLoop<m*n>::MltScal(a, *this);    
    return *this;
  }
  
  
  //! *this = *this + A
  template<class T, int m, int n> template<class T1, class E>
  inline TinyMatrix<T, General, m, n>& TinyMatrix<T, General, m, n>::
  operator +=(const TinyMatrixExpression<T1, m, n, E> & B) 
  {
    TinyMatrixLoop<m*n>::AddCopy(B, *this);
    return *this;
  }
  
  
  //! *this = *this - B
  template<class T, int m, int n> template<class T1, class E>
  inline TinyMatrix<T, General, m, n>& TinyMatrix<T, General, m, n>::
  operator -=(const TinyMatrixExpression<T1, m, n, E> & B) 
  {
    TinyMatrixLoop<m*n>::DiffCopy(B, *this);
    return *this;
  }


  //! *this = x (all values are set to x)
  template<class T, int m, int n>
  inline TinyMatrix<T, General, m, n> & TinyMatrix<T, General, m, n>::operator =(const T& x)
  {
    this->Fill(x);
    return *this;
  }

  
  //! *this = A
  template<class T, int m, int n> template<class T0, class E>
  inline TinyMatrix<T, General, m, n> & TinyMatrix<T, General, m, n>::
  operator =(const TinyMatrixExpression<T0, m, n, E>& A)
  {
    TinyMatrixLoop<m*n>::Copy(A, *this);
    return *this;
  }
  
  
  //! returns A(a, b)
  template<class T, int m, int n>
  inline T& TinyMatrix<T, General, m, n>::operator()(int a, int b)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(a, b, m, n, "TinyMatrix");
#endif
    
    return this->data_[a*n+b];
  }
  
  
  //! returns A(a, b)
  template<class T, int m, int n>
  inline const T& TinyMatrix<T, General, m, n>::operator()(int a, int b) const
  {
    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(a, b, m, n, "TinyMatrix");
#endif
    
    return this->data_[a*n+b];
  }


  //! sets all elements of the matrix to 0
  template<class T, int m, int n>
  inline void TinyMatrix<T, General, m, n>::Zero()
  {
    TinyMatrixLoop<m*n>::Zero(*this);
  }
  

  //! sets matrix to the identity matrix
  template<class T, int m, int n>
  inline void TinyMatrix<T, General, m, n>::SetIdentity()
  {
    
#ifdef SELDON_CHECK_BOUNDS
    if (m != n)
      throw WrongDim("TinyMatrix::SetIdentity()",
		     "The number of rows should be equal to the number of columns");
#endif
    
    TinyMatrixLoop<m*n>::SetIdentity(*this);
  }
  

  //! sets matrix to the identity matrix multiplied by a coefficient alpha
  template<class T, int m, int n> template<class T0>
  inline void TinyMatrix<T, General, m, n>::SetDiagonal(const T0& alpha)
  {
    
#ifdef SELDON_CHECK_BOUNDS
    if (m != n)
      throw WrongDim("TinyMatrix::SetDiagonal(alpha)",
		     "The number of rows should be equal to the number of columns");
#endif
    
    TinyMatrixLoop<m*n>::SetDiagonal(*this, alpha);
  }


  //! sets matrix to [0, 1, 2; 3, 4, 5 ...]
  template<class T, int m, int n>
  inline void TinyMatrix<T, General, m, n>::Fill()
  {
    TinyMatrixLoop<m*n>::Fill(*this);
  }
  
  
  //! sets randomly all the elements of A
  template<class T, int m, int n>
  inline void TinyMatrix<T, General, m, n>::FillRand()
  {
    TinyMatrixLoop<m*n>::FillRand(*this);
  }


  //! sets all the elements of A to a given value
  template<class T, int m, int n> template<class T0>
  inline void TinyMatrix<T, General, m, n>::Fill(const T0& a)
  {
    TinyMatrixLoop<m*n>::Fill(*this, a);
  }
  

  //! returns true if the matrix is null
  template<class T, int m, int n>
  inline bool TinyMatrix<T, General, m, n>::IsZero() const
  {
    return TinyMatrixLoop<m*n>::IsZero(*this);
  }
  

  //! writing matrix in binary format
  template<class T, int m, int n>
  inline void TinyMatrix<T, General, m, n>::Write(const string& file_name) const
  {
    ofstream out(file_name.data());
    Write(out);
    out.close();
  }
  
  
  //! writing matrix in binary format
  template<class T, int m, int n>
  inline void TinyMatrix<T, General, m, n>::Write(ostream& out) const
  {
    int itmp = m;
    out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
    itmp = n;
    out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
    
    return TinyMatrixLoop<m>::Write(out, *this);    
  }


  /************************
   * TinyMatrixTripleLoop *
   ************************/

  
  //!  A(j, :) += val*A(i, :)
  template<int i, int j, int k> template<class T, int m>
  inline void TinyMatrixTripleLoop<i, j, k>::
  AddRow(TinyMatrix<T, General, m, m>& A, const T& val)
  {
    A(j, k-1) += val*A(i, k-1);
    TinyMatrixTripleLoop<i, j, k-1>::AddRow(A, val);
  }
  

  //! Loop with m between 0 and j-1
  template<int j, int k, int m> template<class T, int p>
  inline void TinyMatrixTripleLoop<j, k, m>
  ::ModifyUpperCholesky(TinyMatrix<T, Symmetric, p, p>& A, T& invVal, T& vloc)
  {
    vloc -= A(k, m-1)*A(m-1, j);
    TinyMatrixTripleLoop<j, k, m-1>::ModifyUpperCholesky(A, invVal, vloc);    
  }
  
  
  /************************
   * TinyMatrixDoubleLoop *
   ************************/

  
  //! Conversion from vector of vectors to matrix
  template<int p, int q> template<int m, int n, class T0, class T1>
  inline void TinyMatrixDoubleLoop<p, q>::
  Init(const TinyVector<T0, m>& x, TinyMatrix<T1, General, m, n>& A)
  {
    A(q-1, p-1) = x(q-1);
    TinyMatrixDoubleLoop<p, q-1>::Init(x, A);
  }
    
  
  //! writes matrix A in ascii format
  template<int p, int q> template<int m, int n, class T, class E>
  inline void TinyMatrixDoubleLoop<p, q>::
  WriteText(ostream& out, const TinyMatrixExpression<T, m, n, E>& A)
  {
    out << A(m-p, n-q) << "  ";
    TinyMatrixDoubleLoop<p, q-1>::WriteText(out, A);
  }
  
  
  //! writes matrix A in binary format
  template<int p, int q> template<int m, int n, class T, class Prop>
  inline void TinyMatrixDoubleLoop<p, q>::
  Write(ostream& out, const TinyMatrix<T, Prop, m, n>& A)
  {    
    out.write(reinterpret_cast<const char*>(&A(m-p, n-q)), sizeof(T));
    
    TinyMatrixDoubleLoop<p, q-1>::Write(out, A);
  }
  

  template<int p, int q> template<int m, int n, class T, class E1, class E2>
  inline bool TinyMatrixDoubleLoop<p, q>
  ::IsEqual(const TinyMatrixExpression<T, m, n, E1>& A,
	    const TinyMatrixExpression<T, m, n, E2>& B)
  {
    if (abs(A(p-1, q-1) - B(p-1, q-1)) > TinyVector<T, m>::threshold)
      return false;

    return TinyMatrixDoubleLoop<p, q-1>::IsEqual(A, B);
  }
  

  //! y = A*x
  template<int i, int j> template<int m, int n, class T0, class E0,
				  class T1, class E1, class T2>
  inline void TinyMatrixDoubleLoop<i, j>::
  Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
      const TinyVectorExpression<T1, n, E1>& x, T2& val)
  {
    val += A(i-1, j)*x(j);
    TinyMatrixDoubleLoop<i, j-1>::Mlt(A, x, val);
  }
  

  //! C = A*B
  template<int p, int q>
  template<int m, int n, int k, class T0, class E0,
	   class T1, class E1, class T2, class Prop2>
  inline void TinyMatrixDoubleLoop<p, q>::
  Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
      const TinyMatrixExpression<T1, n, k, E1>& B,
      TinyMatrix<T2, Prop2, m, k>& C)
  {
    C(TinyMatrixNode<Prop2, m, k, p-1 >::i,
      TinyMatrixNode<Prop2, m, k, p-1 >::j)
      += A(TinyMatrixNode<Prop2, m, k, p-1 >::i, q)
	   *B(q, TinyMatrixNode<Prop2, m, k, p-1 >::j);
    
    TinyMatrixDoubleLoop<p, q-1>::Mlt(A, B, C);
  }
  

  //! A = A + alpha * x y^T
  template<int p, int q> 
  template<int m, int n, class T0, class T1, class E1, class T2, class E2, class T3>
  inline void TinyMatrixDoubleLoop<p, q>::
  Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
              const TinyVectorExpression<T2, n, E2>& y, TinyMatrix<T3, General, m, n>& A)
  {
    A(p-1, q-1) += alpha*x(p-1)*y(q-1);
    TinyMatrixDoubleLoop<p, q-1>::Rank1Update(alpha, x, y, A);
  }


  //! A = x y^T
  template<int p, int q> 
  template<int m, int n, class T1, class E1, class T2, class E2, class T3>
  inline void TinyMatrixDoubleLoop<p, q>::
  Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
              const TinyVectorExpression<T2, n, E2>& y, TinyMatrix<T3, General, m, n>& A)
  {
    A(p-1, q-1) = x(p-1)*y(q-1);
    TinyMatrixDoubleLoop<p, q-1>::Rank1Matrix(x, y, A);
  }
  

  //! A = A + alpha * x x^T
  template<int p, int q> 
  template<int m, class T0, class T1, class E1, class T3>
  inline void TinyMatrixDoubleLoop<p, q>::
  Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
              TinyMatrix<T3, Symmetric, m, m>& A)
  {
    A(p-1, q-1) += alpha*x(p-1)*x(q-1);
    TinyMatrixDoubleLoop<p, q-1>::Rank1Update(alpha, x, A);
  }


  //! A = x x^T
  template<int p, int q> 
  template<int m, class T1, class E1, class T3>
  inline void TinyMatrixDoubleLoop<p, q>::
  Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
              TinyMatrix<T3, Symmetric, m, m>& A)
  {
    A(p-1, q-1) = x(p-1)*x(q-1);
    TinyMatrixDoubleLoop<p, q-1>::Rank1Matrix(x, A);
  }
  

  //! returns index jmax for which |A(i, j)| is maximal
  template<int i, int j> template<class T, int m>
  inline void TinyMatrixDoubleLoop<i, j>::GetMaximumColumn(TinyMatrix<T, General, m, m>& A,
                                                           int& jmax, T& val)
  {
    if (abs(A(j, i)) > abs(val))
      {
	jmax = j;
	val = A(j, i);
      }
    
    TinyMatrixDoubleLoop<i, (j-1)*(j-1 > i)>::GetMaximumColumn(A, jmax, val);
  }
  
  
  //! swapping rows i and i2 (for pivoting)
  template<int i, int j> template<class T, int m>
  inline void TinyMatrixDoubleLoop<i, j>::
  SwapRow(TinyMatrix<T, General, m, m>& A, int i2, T& val)
  {
    val = A(i, j-1);
    A(i, j-1) = A(i2, j-1);
    A(i2, j-1) = val;
    TinyMatrixDoubleLoop<i, j-1>::SwapRow(A, i2, val);
  }

  
  //! swapping columns i and i2 (for pivoting)
  template<int i, int j> template<class T, int m>
  inline void TinyMatrixDoubleLoop<i, j>::
  SwapColumn(TinyMatrix<T, General, m, m>& A, int i2, T& val)
  {
    val = A(j-1, i);
    A(j-1, i) = A(j-1, i2);
    A(j-1, i2) = val;
    TinyMatrixDoubleLoop<i, j-1>::SwapColumn(A, i2, val);
  }

  
  //! A(i, :) = A(i, :)*coef
  template<int i, int j> template<class T, int m>
  inline void TinyMatrixDoubleLoop<i, j>::
  MltRow(TinyMatrix<T, General, m, m>& A, const T& coef)
  {
    A(i, j-1) *= coef;
    TinyMatrixDoubleLoop<i, j-1>::MltRow(A, coef);
  }
  

  //! performs Lj = Lj - a_ji/aii Li to eliminate element a_ij
  template<int i, int j> template<class T, int m>
  inline void TinyMatrixDoubleLoop<i, j>::
  PerformElimination(TinyMatrix<T, General, m, m>& A, const T& coef, T& val)
  {
    // val =  -a_ji / a_ii
    val = -A(j, i)*coef;
    
    // replacing Lj by Lj + val * Li
    TinyMatrixTripleLoop<i, j, m>::AddRow(A, val);
    
    // storing val in a_ij
    A(j, i) = val;
    
    TinyMatrixDoubleLoop<i, (j-1)*(j-1>i) >::PerformElimination(A, coef, val);
  }
  

  //! solving by upper matrix
  template<int i, int j1> template<class T, int m>
  inline void TinyMatrixDoubleLoop<i, j1>::
  PerformSolve(TinyMatrix<T, General, m, m>& A, T& val)
  {
    // j = m - 1 - j1 + i + 1
    // val =  -a_ij
    val = -A(i, m - 1 - j1 + i + 1);
    
    // replacing Li by Li + val * Lj
    TinyMatrixTripleLoop<m - 1 - j1 + i + 1, i, m>::AddRow(A, val);
    
    A(i, m - 1 - j1 + i + 1) = val*A(m - 1 - j1 + i + 1, m - 1 - j1 + i + 1);
    TinyMatrixDoubleLoop<i, (j1-1)*(j1-1>i) >::PerformSolve(A, val);
  }
  

    //! Loop val -= A(k, j)*A(k, j) with k between 0 and j-1
  template<int j, int k> template<class T, int m>
  inline void TinyMatrixDoubleLoop<j, k>
  ::GetDiagonalCholesky(TinyMatrix<T, Symmetric, m, m>& A, T& val)
  {
    val -= A(k-1, j)*A(k-1, j);
    TinyMatrixDoubleLoop<j, k-1>::GetDiagonalCholesky(A, val);
  }


  //! Loop with k between j+1 and n-1
  template<int j, int k> template<class T, int m>
  inline void TinyMatrixDoubleLoop<j, k>
  ::ModifyUpperCholesky(TinyMatrix<T, Symmetric, m, m>& A, T& invVal, T& vloc)
  {
    vloc = A(j, m-k);
    TinyMatrixTripleLoop<j, m-k, j>::ModifyUpperCholesky(A, invVal, vloc);
    
    A(j, m-k) = vloc*invVal;
    TinyMatrixDoubleLoop<j, k-1>::ModifyUpperCholesky(A, invVal, vloc);
  }
  

  //! double loop for Cholesky resolution and NoTranspose
  template<int i, int k> template<class T, class T2, int m>
  inline void TinyMatrixDoubleLoop<i, k>
  ::SolveCholesky(const class_SeldonNoTrans& trans,
                  const TinyMatrix<T, Symmetric, m, m>& A,
                  TinyVector<T2, m>& x)
  {
    x(i+k) -= A(i, i+k)*x(i);
    TinyMatrixDoubleLoop<i, k-1>::SolveCholesky(trans, A, x);
  }
  

  //! double loop for Cholesky resolution and Transpose
  template<int i, int k> template<class T, class T2, int m>
  inline void TinyMatrixDoubleLoop<i, k>
  ::SolveCholesky(const class_SeldonTrans& trans,
                  const TinyMatrix<T, Symmetric, m, m>& A,
                  TinyVector<T2, m>& x, T2& val)
  {
    val -= A(i, i+k)*x(i+k);
    TinyMatrixDoubleLoop<i, k-1>::SolveCholesky(trans, A, x, val);
  }
  

  //! double loop for Cholesky multiplication and NoTranspose
  template<int i, int k> template<class T, class T2, int m>
  inline void TinyMatrixDoubleLoop<i, k>
  ::MltCholesky(const class_SeldonNoTrans& trans,
                const TinyMatrix<T, Symmetric, m, m>& A,
                TinyVector<T2, m>& x)
  {
    x(i+k) += A(i, i+k)*x(i);
    TinyMatrixDoubleLoop<i, k-1>::MltCholesky(trans, A, x);
  }
  

  //! double loop for Cholesky multiplication and Transpose
  template<int i, int k> template<class T, class T2, int m>
  inline void TinyMatrixDoubleLoop<i, k>
  ::MltCholesky(const class_SeldonTrans& trans,
                const TinyMatrix<T, Symmetric, m, m>& A,
                TinyVector<T2, m>& x, T2& val)
  {
    val += A(i, i+k)*x(i+k);
    TinyMatrixDoubleLoop<i, k-1>::MltCholesky(trans, A, x, val);
  }
  
  
  /******************
   * TinyMatrixLoop *
   ******************/
  

  //! sets all elements of the matrix to 0
  template<int p> template<int m, int n, class T, class Prop>
  inline void TinyMatrixLoop<p>::Zero(TinyMatrix<T, Prop, m, n>& A)
  {
    FillZero(A.data_[p-1]);
    TinyMatrixLoop<p-1>::Zero(A);
  }
  
  
  //! Conversion from vector of vectors to matrix
  template<int p> template<int m, int n, class T0, class T1>
  inline void TinyMatrixLoop<p>::
  Init(const TinyVector<TinyVector<T0, m>, n>& x, TinyMatrix<T1, General, m, n>& A)
  {
    TinyMatrixDoubleLoop<p, m>::Init(x(p-1), A);
    TinyMatrixLoop<p-1>::Init(x, A);
  }

  
  //! A = alpha*A
  template<int k> template<int m, int n, class T0,
			   class Prop, class T1>
  inline void TinyMatrixLoop<k>::MltScal(const T0& alpha,
					 TinyMatrix<T1, Prop, m, n>& A)
  {
    A.data_[k-1] *= alpha;
    TinyMatrixLoop<k-1>::MltScal(alpha, A);
  }
  

  //! y = x
  template<int p> template<int m, int n, class T1, class E, class T0, class Prop>
  void TinyMatrixLoop<p>::Copy(const TinyMatrixExpression<T1, m, n, E>& x,
			       TinyMatrix<T0, Prop, m, n>& y)
  {
    y(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j)
      = x(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j);

    TinyMatrixLoop<p-1>::Copy(x, y);
  }


  //! y += x
  template<int p> template<int m, int n, class T1, class E, class T0, class Prop>
  void TinyMatrixLoop<p>::AddCopy(const TinyMatrixExpression<T1, m, n, E>& x,
				  TinyMatrix<T0, Prop, m, n>& y)
  {
    y(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j)
      += x(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j);
    
    TinyMatrixLoop<p-1>::AddCopy(x, y);
  }


  //! y -= x
  template<int p> template<int m, int n, class T1, class E, class T0, class Prop>
  void TinyMatrixLoop<p>::DiffCopy(const TinyMatrixExpression<T1, m, n, E>& x,
				   TinyMatrix<T0, Prop, m, n>& y)
  {
    y(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j)
      -= x(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j);
    
    TinyMatrixLoop<p-1>::DiffCopy(x, y);
  }
    

  //! sets matrix to the identity matrix
  template<int p> template<int m, int n, class T, class Prop>
  inline void TinyMatrixLoop<p>::SetIdentity(TinyMatrix<T, Prop, m, n>& A)
  {
    SetComplexReal((TinyMatrixNode<Prop, m, n, p-1>::i == TinyMatrixNode<Prop, m, n, p-1>::j),
                   A(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j));
    
    return TinyMatrixLoop<p-1>::SetIdentity(A);
  }
  

  //! sets matrix to the identity matrix multiplied by a coefficient alpha
  template<int p> template<int m, int n, class T, class Prop, class T0>
  inline void TinyMatrixLoop<p>::SetDiagonal(TinyMatrix<T, Prop, m, n>& A, const T0& alpha)
  {
    if (TinyMatrixNode<Prop, m, n, p-1>::i == TinyMatrixNode<Prop, m, n, p-1>::j)
      SetComplexReal(alpha, A(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j));
    else
      SetComplexZero(A(TinyMatrixNode<Prop, m, n, p-1>::i, TinyMatrixNode<Prop, m, n, p-1>::j));
    
    return TinyMatrixLoop<p-1>::SetDiagonal(A, alpha);
  }


  //! sets matrix to [0, 1, 2; 3, 4, 5 ...]
  template<int p> template<int m, int n, class T, class Prop>
  inline void TinyMatrixLoop<p>::Fill(TinyMatrix<T, Prop, m, n>& A)
  {
    SetComplexReal(p-1, A.data_[p-1]);
    TinyMatrixLoop<p-1>::Fill(A);
  }

  
  //! sets randomly all the elements of A
  template<int p> template<int m, int n, class T, class Prop>
  inline void TinyMatrixLoop<p>::FillRand(TinyMatrix<T, Prop, m, n>& A)
  {
    SetComplexReal(rand(), A.data_[p-1]);
    TinyMatrixLoop<p-1>::FillRand(A);
  }


    //! sets all the elements of A to a given value
  template<int p> template<int m, int n, class T, class Prop, class T0>
  inline void TinyMatrixLoop<p>::Fill(TinyMatrix<T, Prop, m, n>& A, const T0& alpha)
  {
    SetComplexReal(alpha, A.data_[p-1]);
    TinyMatrixLoop<p-1>::Fill(A, alpha);
  }


  //! returns true if the matrix is null
  template<int p> template<int m, int n, class Prop, class T>
  inline bool TinyMatrixLoop<p>::IsZero(const TinyMatrix<T, Prop, m, n>& A)
  {
    if (!IsComplexZero(A.data_[p-1]))
      return false;
    
    return TinyMatrixLoop<p-1>::IsZero(A);
  }
  
  
  //! writes matrix A in ascii format
  template<int p> template<int m, int n, class T, class E>
  inline void TinyMatrixLoop<p>::
  WriteText(ostream& out, const TinyMatrixExpression<T, m, n, E>& A)
  {
    TinyMatrixDoubleLoop<p, n>::WriteText(out, A);
    out << '\n';
    TinyMatrixLoop<p-1>::WriteText(out, A);
  }

  
  //! writes matrix A in binary format
  template<int p> template<int m, int n, class T, class Prop>
  inline void TinyMatrixLoop<p>::Write(ostream& out, const TinyMatrix<T, Prop, m, n>& A)
  {
    TinyMatrixDoubleLoop<p, n>::Write(out, A);
    TinyMatrixLoop<p-1>::Write(out, A);
  }


  //! returns true if A == B
  template<int p> template<int m, int n, class T, class E1, class E2>
  inline bool TinyMatrixLoop<p>
  ::IsEqual(const TinyMatrixExpression<T, m, n, E1>& A,
	    const TinyMatrixExpression<T, m, n, E2>& B)
  {
    if (!TinyMatrixDoubleLoop<p, n>::IsEqual(A, B))
      return false;
    
    return TinyMatrixLoop<p-1>::IsEqual(A, B);
  }


  //! y = A*x
  template<int i> template<int m, int n, class T0, class E0,
			   class T1, class E1, class T2>
  inline void TinyMatrixLoop<i>::
  Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
      const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y)
  {
    y(i-1) = A(i-1, 0)*x(0);
    TinyMatrixDoubleLoop<i, n-1>::Mlt(A, x, y(i-1));
    TinyMatrixLoop<i-1>::Mlt(A, x, y);
  }
  

  //! y = y + A*x
  template<int i> template<int m, int n, class T0, class E0,
			   class T1, class E1, class T2>
  inline void TinyMatrixLoop<i>::
  MltAdd(const TinyMatrixExpression<T0, m, n, E0>& A,
         const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y)
  {
    y(i-1) += A(i-1, 0)*x(0);
    TinyMatrixDoubleLoop<i, n-1>::Mlt(A, x, y(i-1));
    TinyMatrixLoop<i-1>::MltAdd(A, x, y);
  }
  

  //! y = y + alpha*A*x
  template<int i> template<int m, int n, class T0, class E0,
			   class T1, class E1, class T2, class T3>
  inline void TinyMatrixLoop<i>::
  MltAdd(const T3& alpha, const TinyMatrixExpression<T0, m, n, E0>& A,
	 const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y)
  {
    T2 val = A(i-1, 0)*x(0);
    TinyMatrixDoubleLoop<i, n-1>::Mlt(A, x, val);
    y(i-1) += alpha*val;
    TinyMatrixLoop<i-1>::MltAdd(alpha, A, x, y);
  }

  
  //! C = A*B
  template<int p> template<int m, int n, int k, class T0, class E0,
			   class T1, class E1, class T2, class Prop2>
  inline void TinyMatrixLoop<p>::
  Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
      const TinyMatrixExpression<T1, n, k, E1>& B,
      TinyMatrix<T2, Prop2, m, k>& C)
  {
    C(TinyMatrixNode<Prop2, m, k, p-1 >::i, TinyMatrixNode<Prop2, m, k, p-1 >::j)
      = A(TinyMatrixNode<Prop2, m, k, p-1 >::i, 0)*B(0, TinyMatrixNode<Prop2, m, k, p-1 >::j);
    
    TinyMatrixDoubleLoop<p, n-1>::Mlt(A, B, C);
    TinyMatrixLoop<p-1>::Mlt(A, B, C);
  }
  

  //! A = A + alpha * x y^T
  template<int p> template<int m, int n, class T0, class T1, class E1,
			   class T2, class E2, class T3>
  inline void TinyMatrixLoop<p>::
  Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
              const TinyVectorExpression<T2, n, E2>& y, TinyMatrix<T3, General, m, n>& A)
  {
    TinyMatrixDoubleLoop<p, n>::Rank1Update(alpha, x, y, A);
    TinyMatrixLoop<p-1>::Rank1Update(alpha, x, y, A);
  }


  //! A = x y^T
  template<int p> template<int m, int n, class T1, class E1,
			   class T2, class E2, class T3>
  inline void TinyMatrixLoop<p>::
  Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
              const TinyVectorExpression<T2, n, E2>& y, TinyMatrix<T3, General, m, n>& A)
  {
    TinyMatrixDoubleLoop<p, n>::Rank1Matrix(x, y, A);
    TinyMatrixLoop<p-1>::Rank1Matrix(x, y, A);
  }
  

    //! A = A + alpha * x x^T
  template<int p> template<int m, class T0, class T1, class E1, class T3>
  inline void TinyMatrixLoop<p>::
  Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
              TinyMatrix<T3, Symmetric, m, m>& A)
  {
    TinyMatrixDoubleLoop<p, p>::Rank1Update(alpha, x, A);
    TinyMatrixLoop<p-1>::Rank1Update(alpha, x, A);
  }


  //! A = x x^T
  template<int p> template<int m, class T1, class E1, class T3>
  inline void TinyMatrixLoop<p>::
  Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
              TinyMatrix<T3, Symmetric, m, m>& A)
  {
    TinyMatrixDoubleLoop<p, p>::Rank1Matrix(x, A);
    TinyMatrixLoop<p-1>::Rank1Matrix(x, A);
  }

  
  //! x = A(:, k)
  template<int p> template<int m, int n, class T1, class E1>
  inline void TinyMatrixLoop<p>::
  GetCol(const TinyMatrixExpression<T1, m, n, E1>& A, int k, TinyVector<T1, m>& x)
  {
    x(p-1) = A(p-1, k);
    TinyMatrixLoop<p-1>::GetCol(A, k ,x);
  }
  
  
  //! x = A(k, :)
  template<int p> template<int m, int n, class T0, class E0, class T1>
  inline void TinyMatrixLoop<p>::
  GetRow(const TinyMatrixExpression<T0, m, n, E0>& A, int k, TinyVector<T1, n>& x)
  {
    x(p-1) = A(k, p-1);
    TinyMatrixLoop<p-1>::GetRow(A, k ,x);
  }
  
  
  //! A(:, k) = x
  template<int p> template<int m, int n, class T1, class E1, class Prop>
  inline void TinyMatrixLoop<p>::
  SetCol(const TinyVectorExpression<T1, m, E1>& x, int k, TinyMatrix<T1, Prop, m, n>& A)
  {
    A(p-1, k) = x(p-1);
    TinyMatrixLoop<p-1>::SetCol(x, k, A);
  }
  
  
  //! A(k, :) = x
  template<int p> template<int m, int n, class T1, class E1, class Prop>
  inline void TinyMatrixLoop<p>::
  SetRow(const TinyVectorExpression<T1, n, E1>& x, int k, TinyMatrix<T1, Prop, m, n>& A)
  {
    A(k,p-1) = x(p-1);
    TinyMatrixLoop<p-1>::SetRow(x, k, A);
  }
  

    //! computes the maximal element of matrix A
  template<int p> template<int m, int n, class T, class Prop, class T0>
  inline void TinyMatrixLoop<p>::GetMaxAbs(const TinyMatrix<T, Prop, m, n>& A, T0& amax)
  {
    amax = max(amax, abs(A.data_[p-1]));
    TinyMatrixLoop<p-1>::GetMaxAbs(A, amax);
  }
  

  //! step i1 of Gauss elimination
  template<int i1> template<class T, int m>
  inline void TinyMatrixLoop<i1>::
  PivotGauss(TinyMatrix<T, General, m, m>& A, TinyVector<int, m>& pivot)
  {
    // finding the maximum coefficient in the column i
    // i = m -1 - i1
    int jmax = m -1 - i1; T val = A(m -1 - i1, m -1 - i1), coef;
    TinyMatrixDoubleLoop<m -1 - i1, m-1>::GetMaximumColumn(A, jmax, val);
    pivot(m -1 - i1) = jmax;
    
    // now swapping rows i and jmax
    if (m -1 - i1 != jmax)
      TinyMatrixDoubleLoop<m -1 - i1, m>::SwapRow(A, jmax, val);
    
    // loop over each row i+1 -> m, to make elimination (we solve by L as well)
    T one; SetComplexOne(one);
    coef = one/A(m -1 - i1, m -1 - i1);
    TinyMatrixDoubleLoop<m -1 - i1, m-1>::PerformElimination(A, coef, val);
    
    // multiplying row i with coef, and storing coef in a_ii
    TinyMatrixDoubleLoop<m - 1 - i1, m>::MltRow(A, coef);
    A(m -1 - i1, m -1 - i1) = coef;
    
    TinyMatrixLoop<i1-1>::PivotGauss(A, pivot);
  }


  //! solving by upper matrix
  template<int i1> template<class T, int m>
  inline void TinyMatrixLoop<i1>::SolveUpper(TinyMatrix<T, General, m, m>& A)
  {
    // i = m -1 - i1  
    T val;
    TinyMatrixDoubleLoop<m - 1 - i1, m-1>::PerformSolve(A, val);
    TinyMatrixLoop<i1-1>::SolveUpper(A);
  }
  

  //! swapping columns for Gauss pivoting
  template<int i> template<class T, int m>
  inline void TinyMatrixLoop<i>::
  PermuteColumn(TinyMatrix<T, General, m, m>& A, const TinyVector<int, m>& pivot)
  {
    T val;
    if (pivot(i-1) != i-1)
      TinyMatrixDoubleLoop<i-1, m>::SwapColumn(A, pivot(i-1), val);
    
    TinyMatrixLoop<i-1>::PermuteColumn(A, pivot);
  }

  
  //! main loop for Cholesky factorisation
  template<int j> template<class T, int m>
  inline void TinyMatrixLoop<j>::GetCholesky(TinyMatrix<T, Symmetric, m, m>& A)
  {
    T val = A(m-j, m-j);
    TinyMatrixDoubleLoop<m-j, m-j>::GetDiagonalCholesky(A, val);

    val = sqrt(val);
    A(m-j, m-j) = val;
    T one; SetComplexOne(one);
    T invVal = one / val, vloc;
    
    TinyMatrixDoubleLoop<m-j, j-1>::ModifyUpperCholesky(A, invVal, vloc);
    
    TinyMatrixLoop<j-1>::GetCholesky(A);
  }


  //! simple loop for Cholesky resolution and NoTranspose
  template<int i> template<class T, class T2, int m>
  inline void TinyMatrixLoop<i>::SolveCholesky(const class_SeldonNoTrans& trans,
                                               const TinyMatrix<T, Symmetric, m, m>& A,
                                               TinyVector<T2, m>& x)
  {
    x(m-i) /= A(m-i, m-i);
    TinyMatrixDoubleLoop<m-i, i-1>::SolveCholesky(trans, A, x);
    TinyMatrixLoop<i-1>::SolveCholesky(trans, A, x);
  }
  
  
  //! simple loop for Cholesky resolution and Transpose
  template<int i> template<class T, class T2, int m>
  inline void TinyMatrixLoop<i>::SolveCholesky(const class_SeldonTrans& trans,
                                               const TinyMatrix<T, Symmetric, m, m>& A,
                                               TinyVector<T2, m>& x)
  {
    T2 val = x(i-1);
    TinyMatrixDoubleLoop<i-1, m-i>::SolveCholesky(trans, A, x, val);
    
    x(i-1) = val / A(i-1, i-1);
    TinyMatrixLoop<i-1>::SolveCholesky(trans, A, x);
  }


  //! simple loop for Cholesky multiplication and No Transpose
  template<int i> template<class T, class T2, int m>
  inline void TinyMatrixLoop<i>::MltCholesky(const class_SeldonNoTrans& trans,
                                             const TinyMatrix<T, Symmetric, m, m>& A,
                                             TinyVector<T2, m>& x)
  {
    TinyMatrixDoubleLoop<i-1, m-i>::MltCholesky(trans, A, x);
    
    x(i-1) *= A(i-1, i-1);
    TinyMatrixLoop<i-1>::MltCholesky(trans, A, x);
  }


  //! simple loop for Cholesky multiplication and Transpose
  template<int i> template<class T, class T2, int m>
  inline void TinyMatrixLoop<i>::MltCholesky(const class_SeldonTrans& trans,
                                             const TinyMatrix<T, Symmetric, m, m>& A,
                                             TinyVector<T2, m>& x)
  {
    T val = x(m-i)*A(m-i, m-i);
    TinyMatrixDoubleLoop<m-i, i-1>::MltCholesky(trans, A, x, val);
    x(m-i) = val;
    TinyMatrixLoop<i-1>::MltCholesky(trans, A, x);
  }
  
  
  /*************************
   * TinyMatrix<Symmetric> *
   *************************/
  
  
  //! Default constructor : all elements of the matrix are set to 0
  template<class T, int m>
  inline TinyMatrix<T, Symmetric, m, m>::TinyMatrix()
  {
    this->Zero();
  }
  
  
  //! Default constructor : all elements of the matrix are set to 0
  template<class T, int m>
  inline TinyMatrix<T, Symmetric, m, m>::TinyMatrix(int i)
  {
    this->Zero();
  }


  //! Default constructor : all elements of the matrix are set to 0
  template<class T, int m>
  inline TinyMatrix<T, Symmetric, m, m>::TinyMatrix(int i, int j)
  {
    this->Zero();
  }
  

  //! Constructor taking an expression
  template<class T, int m> template<class E>
  inline TinyMatrix<T, Symmetric, m, m>
  ::TinyMatrix(const TinyMatrixExpression<T, m, m, E>& A)
  {
    TinyMatrixLoop<m*(m+1)/2>::Copy(A, *this);
  }
  
  
  //! returns the number of rows
  template<class T, int m>
  inline int TinyMatrix<T, Symmetric, m, m>::GetM()
  {
    return m;
  }
  
  
  //! returns the number of columns
  template<class T, int m>
  inline int TinyMatrix<T, Symmetric, m, m>::GetN()
  {
    return m;
  }
  
  
  //! returns the number of stored elements
  template<class T, int m>
  inline int TinyMatrix<T, Symmetric, m, m>::GetSize()
  {
    return m*(m+1)/2;
  }

  
  //! returns the pointer to the C-array
  template<class T, int m>
  inline T* TinyMatrix<T, Symmetric, m, m>::GetData()
  {
    return data_;
  }  
  

  //! *this = *this*a
  template<class T, int m> template<class T1>
  inline TinyMatrix<T, Symmetric, m, m> & TinyMatrix<T, Symmetric, m, m>::
  operator *=(const T1& a )
  {
    TinyMatrixLoop<m*(m+1)/2>::MltScal(a, *this);
    
    return *this;
  }
  
  
  //! *this = *this + A
  template<class T, int m> template<class T1, class E>
  inline TinyMatrix<T, Symmetric, m, m> & TinyMatrix<T, Symmetric, m, m>::
  operator +=(const TinyMatrixExpression<T1, m, m, E>& B )
  {
    TinyMatrixLoop<m*(m+1)/2>::AddCopy(B, *this);
    return *this;
  }
  
  
  //! *this = *this -B
  template<class T, int m> template<class T1, class E>
  inline TinyMatrix<T, Symmetric, m, m> & TinyMatrix<T, Symmetric, m, m>::
  operator -=(const TinyMatrixExpression<T1, m, m, E>& B )
  {
    TinyMatrixLoop<m*(m+1)/2>::DiffCopy(B, *this);
    return *this;
  }
  
  
  //! *this = x (sets all values to x)
  template<class T, int m>
  inline TinyMatrix<T, Symmetric, m, m> & 
  TinyMatrix<T, Symmetric, m, m>::operator =(const T& x)
  {
    this->Fill(x);
    return *this;
  }
  
  
  //! *this = x
  template<class T, int m> template<class T0, class E>
  inline TinyMatrix<T, Symmetric, m, m> & TinyMatrix<T, Symmetric, m, m>::
  operator =(const TinyMatrixExpression<T0, m, m, E>& x)
  {
    TinyMatrixLoop<m*(m+1)/2>::Copy(x, *this);
    return *this;
  }
  
  
  //! returns A(i, j)
  template<class T, int m>
  inline T& TinyMatrix<T, Symmetric, m, m>::operator()(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, m, m, "TinyMatrix");
#endif
    
    return this->data_[j > i ? i*m - (i*(i+1))/2 + j: j*m - (j*(j+1))/2 + i];
  }
  
  
  //! returns A(i, j)
  template<class T, int m>
  inline const T& TinyMatrix<T, Symmetric, m, m>::operator()(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, m, m, "TinyMatrix");
#endif
    
    return this->data_[j > i ? i*m - (i*(i+1))/2 + j: j*m - (j*(j+1))/2 + i ];
  }


    //! sets all elements of the matrix to 0
  template<class T, int m>
  inline void TinyMatrix<T, Symmetric, m, m>::Zero()
  {
    TinyMatrixLoop<m*(m+1)/2>::Zero(*this);
  }
  
  
  //! sets the matrix to the identity matrix
  template<class T, int m>
  inline void TinyMatrix<T, Symmetric, m, m>::SetIdentity()
  {
    return TinyMatrixLoop<m*(m+1)/2>::SetIdentity(*this);
  }


  //! sets the matrix to the identity matrix multiplied by a coefficient alpha
  template<class T, int m> template<class T0>
  inline void TinyMatrix<T, Symmetric, m, m>::SetDiagonal(const T0& alpha)
  {
    return TinyMatrixLoop<m*(m+1)/2>::SetDiagonal(*this, alpha);
  }
  
  
  //! sets matrix to [0, 1, 2; 1, 3, 4 ...]
  template<class T, int m>
  inline void TinyMatrix<T, Symmetric, m, m>::Fill()
  {
    return TinyMatrixLoop<m*(m+1)/2>::Fill(*this);
  }
  
  
  //! sets randomly all the elements of A
  template<class T, int m>
  inline void TinyMatrix<T, Symmetric, m, m>::FillRand()
  {
    return TinyMatrixLoop<m*(m+1)/2>::FillRand(*this);
  }
  
  
  //! sets matrix to [0, 1, 2; 1, 3, 4 ...]
  template<class T, int m> template<class T0>
  inline void TinyMatrix<T, Symmetric, m, m>::Fill(const T0& a)
  {
    return TinyMatrixLoop<m*(m+1)/2>::Fill(*this, a);
  }
  

  //! returns true if the matrix is null
  template<class T, int m>
  inline bool TinyMatrix<T, Symmetric, m, m>::IsZero() const
  {
    return TinyMatrixLoop<m*(m+1)/2>::IsZero(*this);
  }
  

  //! writing matrix in ascii format
  template<class T, int m>
  inline void TinyMatrix<T, Symmetric, m, m>::WriteText(const string& file_name) const
  {
    ofstream out(file_name.data());
    out.precision(cout.precision());
    out << *this;
    out.close();
  }


  //! writing matrix in binary format
  template<class T, int m>
  inline void TinyMatrix<T, Symmetric, m, m>::Write(const string& file_name) const
  {
    ofstream out(file_name.data());
    Write(out);
    out.close();
  }
  
  
  //! writing matrix in binary format
  template<class T, int m>
  inline void TinyMatrix<T, Symmetric, m, m>::Write(ostream& out) const
  {
    int itmp = m;
    out.write(reinterpret_cast<char*>(&itmp), sizeof(int));    
    out.write(reinterpret_cast<char*>(&itmp), sizeof(int));    
    
    return TinyMatrixLoop<m>::Write(out, *this);    
  }


  /*************
   * Operators *
   *************/

  
  //! returns true if *this == u
  template<class T, int m, int n, class E1, class E2>
  inline bool operator==(const TinyMatrixExpression<T, m, n, E1> & u,
			 const TinyMatrixExpression<T, m, n, E2> & v)
  {
    return TinyMatrixLoop<m>::IsEqual(u, v);
  }
  
  
  //! returns true if *this != u
  template<class T, int m, int n, class E1, class E2>
  inline bool operator!=(const TinyMatrixExpression<T, m, n, E1> & u,
			 const TinyMatrixExpression<T, m, n, E2> & v)
  {
    return !TinyMatrixLoop<m>::IsEqual(u, v);
  }
  

  //! displays matrix
  template <class T, int m, int n, class E>
  inline ostream& operator <<(ostream& out, const TinyMatrixExpression<T, m, n, E> & A)
  {
    TinyMatrixLoop<m>::WriteText(out, A);
    return out;
  }
  
  
  /***************************
   * Matrix-vector functions *
   ***************************/

  
  //! returns A*x
  template<class T, int m, int n, class E, class T1, class E1> inline TinyVector<T1, m>
  dot(const TinyMatrixExpression<T, m, n, E>& A,
      const TinyVectorExpression<T1, n, E1>& x)
  {
    TinyVector<T1, m> b;
    TinyMatrixLoop<m>::Mlt(A, x, b);
    
    return b;
  }
  
  
  //! returns A*B
  template<class T, int m, int n, class E, class T1, int k, class E1>
  inline TinyMatrix<T, General, m, k>
  dot(const TinyMatrixExpression<T, m, n, E>& A, 
      const TinyMatrixExpression<T1, n, k, E1> & B)
  {
    TinyMatrix<T, General, m, k> C;
    TinyMatrixLoop<m*k>::Mlt(A, B, C);
    
    return C;
  }

  
  //! y = A*x
  template<class T0, class E0, class T1, class E1, class T2, int m, int n>
  inline void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
                  const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y)
  {
    TinyMatrixLoop<m>::Mlt(A, x, y);
  }
  
  
  //! y = y + A*x
  template<class T1, class E1, class T2, class E2, class T3, int m, int n>
  inline void MltAdd(const TinyMatrixExpression<T1, m, n, E1>& A,
                     const TinyVectorExpression<T2, n, E2>& x, TinyVector<T3, m>& y)
  {
    TinyMatrixLoop<m>::MltAdd(A, x, y);
  }

  
  //! y = y + alpha*A*x
  template<class T0, class T1, class E1, class T2, class E2, class T3, int m, int n>
  inline void MltAdd(const T0& alpha, const TinyMatrixExpression<T1, m, n, E1>& A,
                     const TinyVectorExpression<T2, n, E2>& x, TinyVector<T3, m>& y)
  {
    TinyMatrixLoop<m>::MltAdd(alpha, A, x, y);
  }
   
  
  //! y = A^t * x
  template<class T0, class E0, class T1, class E1, class T2, int m, int n>
  inline void MltTrans(const TinyMatrixExpression<T0, m, n, E0>& A,
                       const TinyVectorExpression<T1, m, E1>& x, TinyVector<T2, n>& y)
  {
    TinyMatrixLoop<n>::Mlt(transpose(A), x, y);
  }


  //! y = A^t * x
  template<class T0, class E0, class T1, class E1, class T2, int m, int n>
  inline void Mlt(const class_SeldonTrans&,
		  const TinyMatrixExpression<T0, m, n, E0>& A,
		  const TinyVectorExpression<T1, m, E1>& x, TinyVector<T2, n>& y)
  {
    TinyMatrixLoop<n>::Mlt(transpose(A), x, y);
  }
  
  
  //! A = A + alpha * x y^T
  template<int m, int n, class T0, class T1, class E1, class T2, class E2, class T3>
  inline void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                          const TinyVectorExpression<T2, n, E2>& y,
			  TinyMatrix<T3, General, m, n>& A)
  {
    TinyMatrixLoop<m>::Rank1Update(alpha, x, y, A);
  }


  //! A = x y^T
  template<int m, int n, class T1, class E1, class T2, class E2, class T3>
  inline void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                          const TinyVectorExpression<T2, n, E2>& y,
			  TinyMatrix<T3, General, m, n>& A)
  {
    TinyMatrixLoop<m>::Rank1Matrix(x, y, A);
  }
  

  //! A = A + alpha * x y^T
  template<int m, class T0, class T1, class E1, class T3>
  inline void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                          TinyMatrix<T3, Symmetric, m, m>& A)
  {
    TinyMatrixLoop<m>::Rank1Update(alpha, x, A);
  }


  //! A = x x^T
  template<int m, class T1, class E1, class T3>
  inline void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                          TinyMatrix<T3, Symmetric, m, m>& A)
  {
    TinyMatrixLoop<m>::Rank1Matrix(x, A);
  }

  
  //! x = A(:, k)
  template<int m, int n, class T1, class E1>
  inline void GetCol(const TinyMatrixExpression<T1, m, n, E1>& A, int k, TinyVector<T1, m>& x)
  {
    TinyMatrixLoop<m>::GetCol(A, k, x);
  }
  
  
  //! x = A(k, :)
  template<int m, int n, class T0, class E0, class T1>
  inline void GetRow(const TinyMatrixExpression<T0, m, n, E0>& A, int k, TinyVector<T1, n>& x)
  {
    TinyMatrixLoop<n>::GetRow(A, k, x);
  }
  
  
  //! A(:, k) = x
  template<int m, int n, class T1, class E1, class Prop>
  inline void SetCol(const TinyVectorExpression<T1, m, E1>& x, int k, TinyMatrix<T1, Prop, m, n>& A)
  {
    TinyMatrixLoop<m>::SetCol(x, k, A);
  }
  
  
  //! A(k, :) = x
  template<int m, int n, class T1, class E1, class Prop>
  inline void SetRow(const TinyVectorExpression<T1, n, E1>& x, int k, TinyMatrix<T1, Prop, m, n>& A)
  {
    TinyMatrixLoop<n>::SetRow(x, k, A);
  }
  
  
  /****************************
   * Matrix-Matrix operations *
   ****************************/
  
  
  //! B = A
  template<class T, int m, int n, class E, class Prop>
  inline void Copy(const TinyMatrixExpression<T, m, n, E>& A, TinyMatrix<T, Prop, m, n>& B)
  {
    B = A;
  }
  
  
  //! C = A + B
  template<class T0, class E0, class T1, class E1,
	   class T2, class Prop2, int m, int n>
  inline void Add(const TinyMatrixExpression<T0, m, n, E0>& A,
		  const TinyMatrixExpression<T1, m, n, E1>& B,
		  TinyMatrix<T2, Prop2, m, n>& C)
  {
    C = A + B;
  }
  
  
  //! B = B + alpha*A
  template<class T1, class E1, class T2,
	   class Prop2, int m, int n>
  inline void Add(const T1& alpha, const TinyMatrixExpression<T1, m, n, E1>& A,
		  TinyMatrix<T2, Prop2, m, n>& B)
  {
    B += alpha*A;
  }
  
  
  //! B = B+A
  template<class T0, class E0, class T1, class Prop1, int m, int n>
  inline void Add(const TinyMatrixExpression<T0, m, n, E0>& A, TinyMatrix<T1, Prop1, m, n>& B)
  {
    B += A;
  }
  
  
  //! A = A*alpha
  template<class T0, class T1, class Prop, int m, int n>
  inline void Mlt(const T0& alpha, TinyMatrix<T1, Prop, m, n>& A)
  {
    A *= alpha;
  }
  
    
  //! C = A*B
  template<class T0, class E0, class T1, class E1,
	   class T2, class Prop2, int m, int n, int k>
  inline void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
		  const TinyMatrixExpression<T1, n, k, E1>& B,
		  TinyMatrix<T2, Prop2, m, k>& C)
  {
    TinyMatrixLoop<TinyMatrix<T2, Prop2, m, k>::size_ >::Mlt(A, B, C);
  }
  
  
  //! C = A*B^T
  template<class T0, class E0, class T1, class E1,
	   class T2, class Prop2, int m, int n, int k>
  inline void MltTrans(const TinyMatrixExpression<T0, m, n, E0>& A,
		       const TinyMatrixExpression<T1, k, n, E1>& B,
		       TinyMatrix<T2, Prop2, m, k>& C)
  {
    TinyMatrixLoop<TinyMatrix<T2, Prop2, m, k>::size_ >::Mlt(A, transpose(B), C);
  }


  //! C = beta C + alpha A^T B
  template<class T3, class T0, class E0, class T1, class E1,
	   class T2, class Prop2, class T4, int m, int n, int k>
  inline void MltAdd(const T3& alpha, const class_SeldonTrans&,
		     const TinyMatrixExpression<T0, n, m, E0>& A,
		     const class_SeldonNoTrans&, const TinyMatrixExpression<T1, n, k, E1>& B,
		     const T4& beta, TinyMatrix<T2, Prop2, m, k>& C)
  {
    TinyMatrix<T2, Prop2, m, k> D;
    Mlt(transpose(A), B, D);
    C = beta*C + alpha*D;
  }
  
  
  //! B = A^T
  template<class T, int m, int n, class E>
  inline void Transpose(const TinyMatrixExpression<T, m, n, E> & A,
                        TinyMatrix<T, General, n, m> & B)
  {
    B = transpose(A);
  }
  
  
  //! replacing B by its transposed
  template<class T, int m>
  inline void Transpose(TinyMatrix<T, General, m, m>& B)
  {
    TinyMatrix<T, General, m, m> A(B);
    B = transpose(A);
  }

  
  template<class T, class Prop, int m, int n> inline
  typename ClassComplexType<T>::Treal MaxAbs(const TinyMatrix<T, Prop, m, n>& A)
  {
    typename ClassComplexType<T>::Treal amax(0);
    TinyMatrixLoop<TinyMatrix<T, Prop, m, n>::size_>::GetMaxAbs(A, amax);
    return amax;
  }
  
  
  /****************
   * 1x1 matrices *
   ****************/
  
  
  //! returns determinant of matrix A
  template<class T, class Prop>
  inline T Det(const TinyMatrix<T, Prop, 1, 1> & A)
  {
    return A(0, 0);
  }
  
  
  //! B = A^-1
  template<class T>
  inline void GetInverse(const TinyMatrix<T, General, 1, 1> & A,
			 TinyMatrix<T, General, 1, 1> & B)
  {
    T one; SetComplexOne(one);
    B(0,0) = one/A(0, 0);
  }
  
  
  //! B = A^-1
  template<class T>
  inline void GetInverse(const TinyMatrix<T, Symmetric, 1, 1> & A,
			 TinyMatrix<T, Symmetric, 1, 1> & B)
  {
    T one; SetComplexOne(one);
    B(0,0) = one/A(0, 0);
  }

  
  //! replaces B by its inverse
  template<class T>
  inline void GetInverse(TinyMatrix<T, General, 1, 1> & B)
  {
    T one; SetComplexOne(one);
    B(0,0) = one/B(0, 0);
  }

  
  //! replaces B by its inverse
  template<class T>
  inline void GetInverse(TinyMatrix<T, Symmetric, 1, 1> & B)
  {
    T one; SetComplexOne(one);
    B(0,0) = one/B(0, 0);
  }


  //! computation of n n^T where n is the normale
  //! this is the projector to normale
  template<class T0, class T1>
  inline void GetNormalProjector(const TinyVector<T0, 1>& n, TinyMatrix<T1, Symmetric, 1, 1>& P)    
  {
    P(0,0) = n(0)*n(0);
  }

  
  //! computation of n n^T where n is the normale
  //! this is the projector to normale
  template<class T0, class T1>
  inline void GetNormalProjector(const T0& n, TinyMatrix<T1, Symmetric, 1, 1>& P)    
  {
    P(0,0) = n*n;
  }

  
  /****************
   * 2x2 matrices *
   ****************/
  
  
  //! returns determinant of matrix A
  template<class T, class Prop>
  inline T Det(const TinyMatrix<T, Prop, 2, 2> & A)
  {
    return A(0,0)*A(1,1) - A(0,1)*A(1,0);
  }

  
  //! computation of I - n n^T where n is the normale
  //! this is the projector to tangential plane
  template<class T0, class T1>
  inline void GetTangentialProjector(const TinyVector<T0, 2>& n,
                                     TinyMatrix<T1, Symmetric, 2, 2>& P)
  {
    P(0, 0) = n(1)*n(1);
    P(0, 1) = -n(0)*n(1);
    P(1, 1) = n(0)*n(0);
  }
  
  
  //! computation of n n^T where n is the normale
  //! this is the projector to normale
  template<class T0, class T1>
  inline void GetNormalProjector(const TinyVector<T0, 2>& n, TinyMatrix<T1, Symmetric, 2, 2>& P)
  {
    P(0, 0) = n(0)*n(0);
    P(0, 1) = n(0)*n(1);
    P(1, 1) = n(1)*n(1);
  }


  //! computation of n n^T where n is the normale
  //! this is the projector to normale
  template<class T0, class T1>
  inline void GetNormalProjector(const TinyVector<T0, 2>& n, TinyMatrix<T1, General, 2, 2>& P)
  {
    P(0, 0) = n(0)*n(0);
    P(0, 1) = n(0)*n(1);
    P(1, 1) = n(1)*n(1);
    P(1, 0) = P(0, 1);
  }
  
  
  /**************************
   * Functions for any size *
   **************************/
  
  
  //! tangential projector P = I - n n^T for any size
  template<class T0, class T1, int m>
  inline void GetTangentialProjector(const TinyVector<T0, m>& n,
                                     TinyMatrix<T1, Symmetric, m, m>& P)
  {
    T0 one; SetComplexOne(one);
    P.SetIdentity();
    Rank1Update(-one, n, P);
  }
  
  
  //! normal projector P = n n^T for any size
  template<class T0, class T1, int m>
  inline void GetNormalProjector(const TinyVector<T0, m>& n,
                                 TinyMatrix<T1, Symmetric, m, m>& P)
  {
    Rank1Matrix(n, P);
  }

    
  template<class T, class Prop, int m, int n>
  inline void FillZero(TinyMatrix<T, Prop, m, n>& X)
  {
    X.Zero();
  }
  
  
  //! res = res + A . B, where the scalar product is performed between columns of A and B
  template<class T, int m, int n>
  inline void DotProdCol(const TinyMatrix<T, General, m, n> &A,
                         const TinyMatrix<T, General, m, n>& B, TinyVector<T, n> &res)
  {    
    TinyVector<T, m> vectA, vectB;
    for (int p = 0; p < n ; p ++)
      {
        GetCol(A, p, vectA);
        GetCol(B, p, vectB);
        res(p) += DotProd(vectA, vectB);
      }
  }
  
} // end namespace

#define MONTJOIE_FILE_TINY_MATRIX_INLINE_CXX
#endif
