// Copyright (C) 2014 INRIA
// Author(s): Marc Duruflé
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.

#ifndef SELDON_FILE_BAND_MATRIX_INLINE_CXX

#include "BandMatrix.hxx"

namespace Seldon
{

  /***************
   * Matrix_Band *
   ***************/
  
  
  //! returns the number of rows
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Band<T, Prop, Storage, Allocator>::GetM() const
  {
    return this->m_;
  }


  //! returns the number of extra-diagonals in lower part of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Band<T, Prop, Storage, Allocator>::GetKL() const
  {
    return kl_;
  }


  //! returns the number of extra-diagonals in upper part of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Band<T, Prop, Storage, Allocator>::GetKU() const
  {
    return ku_;
  }

  
  //! returns the number of rows
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Band<T, Prop, Storage, Allocator>::GetN() const
  {
    return this->n_;
  }
  

  //! returns the number of elements stored in the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_Band<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return data_.GetDataSize();
  }
  
  
  //! returns the memory used by the object in bytes
  template <class T, class Prop, class Storage, class Allocator>
  inline size_t Matrix_Band<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    return sizeof(*this) + data_.GetMemorySize() - sizeof(data_);
  }
  
  
  //! fills non-zero entries with 0
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>::Zero()
  {
    data_.Zero();
  }
  
  
  //! present for compatibility
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>::HideMessages()
  {
  }
  
  
  //! changes the size of the matrix, previous entries are lost
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::Reallocate(int m, int n)
  {
    Reallocate(m, n, 0, 0);
  }

  
  //! returns a pointer to the array containing values
  template <class T, class Prop, class Storage, class Allocator>
  inline T* Matrix_Band<T, Prop, Storage, Allocator>::GetData() const
  {
    return data_.GetData();
  }
  
  
  //! multiplication by a scalar
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix<T, Prop, Storage, Allocator>& 
  Matrix_Band<T, Prop, Storage, Allocator>::operator *=(const T& alpha)
  {
    data_ *= alpha;
    return static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);
  }
  
  
  //! returns a reference to A(i, j) 
  template <class T, class Prop, class Storage, class Allocator>
  inline T& Matrix_Band<T, Prop, Storage, Allocator>::Val(int i, int j)
  {
    return Get(i, j);
  }
  
  
  //! returns a reference to A(i, j) 
  template <class T, class Prop, class Storage, class Allocator>
  inline const T& Matrix_Band<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {
    return Get(i, j);
  }
  
  
  //! sets A(i, j)   
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }


  //! sets A(i, j)   
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::SetEntry(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }

    
#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Tcplx>& x, const Vector<Tcplx>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    Seldon::
      MltAdd(alpha,
	     static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	     x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    Seldon::
      MltAdd(alpha,
	     static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	     x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    Seldon::
      MltAdd(alpha, trans,
	     static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	     x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    Seldon::
      MltAdd(alpha, trans,
	     static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	     x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Band<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_Band<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return false;
  }
#endif


  /****************
   * Matrix_Arrow *
   ****************/
  

  //! default constructor
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Arrow<T, Prop, Storage, Allocator>::Matrix_Arrow()
    : Matrix_Band<T, Prop, Storage, Allocator>()
  {
  }
  
  
  //! returns the number of rows
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Arrow<T, Prop, Storage, Allocator>::GetM() const
  {
    return this->m_ + last_rows_.GetM();
  }
  
  
  //! returns the number of columns
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Arrow<T, Prop, Storage, Allocator>::GetN() const
  {
    return this->n_ + last_columns_.GetN();
  }

  
  //! returns the number of dense rows placed at the end of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Arrow<T, Prop, Storage, Allocator>::GetNbLastRow() const
  {
    return last_rows_.GetM();
  }
  
  
  //! returns the number of dense columns placed at the end of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Arrow<T, Prop, Storage, Allocator>::GetNbLastCol() const
  {
    return last_columns_.GetN();
  }
  

  //! returns the number of non-zero entries
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_Arrow<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return this->data_.GetDataSize() + this->last_rows_.GetDataSize() 
      + this->last_columns_.GetDataSize() + this->last_block_.GetDataSize();
  }


  //! returns the memory used by the object in bytes
  template <class T, class Prop, class Storage, class Allocator>
  inline size_t Matrix_Arrow<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    return Matrix_Band<T, Prop, Storage, Allocator>::GetMemorySize()
      + last_rows_.GetMemorySize() + last_columns_.GetMemorySize()
      + last_block_.GetMemorySize();
  }
  
  
  //! present for compatibility
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Arrow<T, Prop, Storage, Allocator>::HideMessages()
  {
  }
  
  
  //! changes the size of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Arrow<T, Prop, Storage, Allocator>::Reallocate(int m, int n)
  {
    Reallocate(m, n, 0, 0);
  }
  
  
  //! returns a reference to A(i, j)   
  template <class T, class Prop, class Storage, class Allocator>
  inline T& Matrix_Arrow<T, Prop, Storage, Allocator>::Val(int i, int j)
  {
    return Get(i, j);
  }
  
  
  //! returns a reference to A(i, j) 
  template <class T, class Prop, class Storage, class Allocator>
  inline const T& Matrix_Arrow<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {
    return Get(i, j);
  }
    

  //! sets A(i, j)     
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Arrow<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }


  //! sets A(i, j)     
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Arrow<T, Prop, Storage, Allocator>
  ::SetEntry(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }

  
  /**************
   * Functions *
   *************/


  //! B = B + alpha*A
  template<class T0, class T1, class Prop1, class Allocator1,
	   class T2, class Prop2, class Allocator2>
  inline void AddMatrix(const T0& alpha,
			const Matrix<T1, Prop1, BandedCol, Allocator1>& A,
			Matrix<T2, Prop2, BandedCol, Allocator2>& B)           
  {
    B.Add_(alpha, A);
  }

  
  //! conversion from ArrayRowSparse to band matrix
  template<class T, class Allocator>
  inline void Copy(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
		   Matrix<T, General, BandedCol, Allocator>& B)
  {
    B.Copy(A);
  }
  
  
  //! resolution of A x = b, once GetLU has been called
  template<class T, class Allocator>
  inline void SolveLU(const Matrix<T, General, BandedCol, Allocator>& mat_lu,
		      Vector<T>& x)
  {
    mat_lu.Solve(x);
  }


  //! resolution of A x = b, once GetLU has been called
  template<class T, class Allocator>
  inline void SolveLU(const Matrix<T, General, BandedCol, Allocator>& mat_lu,
		      Vector<complex<T> >& x)
  {
    mat_lu.Solve(x);
  }
  

  //! B = B + alpha*A
  template<class T0, class T1, class Allocator>
  inline void AddVector(const T0& alpha,
			const Matrix<T1, General, BandedCol, Allocator>& A,
			Matrix<T1, General, BandedCol, Allocator>& B)           
  {
    B.Add_(alpha, A);
  }
  
  
  //! y = beta*y + alpha*A*x
  template<class T0, class T1, class T, class T2, class Allocator>
  inline void MltAddVector(const T0& alpha,
			   const Matrix<T, General, BandedCol, Allocator>& A,
			   const Vector<T2>& x, const T1& beta, Vector<T2>& y)
  {
    Mlt(beta, y);
    A.MltAdd(alpha, SeldonNoTrans, x, y);
  }


  //! y = beta*y + alpha*A*x
  template<class T0, class T1, class T, class T2, class Allocator>
  inline void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
			   const Matrix<T, General, BandedCol, Allocator>& A,
			   const Vector<T2>& x, const T1& beta, Vector<T2>& y)
  {
    Mlt(beta, y);
    A.MltAdd(alpha, trans, x, y);
  }


  //! y = A*x
  template<class T, class Allocator, class T1>
  inline void MltVector(const Matrix<T, General, BandedCol, Allocator>& A,
			const Vector<T1>& x, Vector<T1>& y)
  {
    T1 zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    y.Fill(zero);
    A.MltAdd(one, SeldonNoTrans, x, y);
  }


  //! y = A*x
  template<class T, class Allocator, class T1>
  inline void MltVector(const SeldonTranspose& trans,
			const Matrix<T, General, BandedCol, Allocator>& A,
			const Vector<T1>& x, Vector<T1>& y)
  {
    T1 zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    y.Fill(zero);
    A.MltAdd(one, trans, x, y);
  }

  
  //! A = alpha*A
  template<class T0, class T1, class Allocator>
  inline void MltScalar(const T0& alpha,
			Matrix<T1, General, BandedCol, Allocator>& A)
  {
    A *= alpha;
  }
  
  
  //! displays matrix
  template<class T, class Allocator>
  inline ostream& operator<<(ostream& out,
			     const Matrix<T, General, BandedCol, Allocator>& A)
  {
    A.WriteText(out);
    return out;
  }


  //! LU factorisation
  template<class T, class Allocator>
  inline void GetLU(Matrix<T, General, ArrowCol, Allocator>& A,
		    Matrix<T, General, ArrowCol, Allocator>& mat_lu,
		    bool keep_matrix)
  {
    mat_lu = A;
    if (!keep_matrix)
      A.Clear();
    
    mat_lu.Factorize();
  }
  

  //! LU factorisation
  template<class T, class Allocator>
  inline void GetLU(Matrix<T, General, ArrowCol, Allocator>& A)
  {
    A.Factorize();
  }  
  
  
  //! resolution of A x = b, once GetLU has been called
  template<class T, class Allocator>
  inline void SolveLU(const Matrix<T, General, ArrowCol, Allocator>& mat_lu,
		      Vector<T>& x)
  {
    mat_lu.Solve(x);
  }


  //! resolution of A x = b, once GetLU has been called
  template<class T, class Allocator>
  inline void SolveLU(const Matrix<T, General, ArrowCol, Allocator>& mat_lu,
		      Vector<complex<T> >& x)
  {
    mat_lu.Solve(x);
  }
  

  //! B = B + alpha*A
  template<class T0, class T1, class Prop1, class Allocator1,
	   class T2, class Prop2, class Allocator2>
  inline void AddMatrix(const T0& alpha,
			const Matrix<T1, Prop1, ArrowCol, Allocator1>& A,
			Matrix<T2, Prop2, ArrowCol, Allocator2>& B)           
  {
    B.Add_(alpha, A);
  }
  
  
  //! y = beta*y + alpha*A*x
  template<class T0, class T1, class T, class T2, class Allocator>
  inline void MltAddVector(const T0& alpha,
			   const Matrix<T, General, ArrowCol, Allocator>& A,
			   const Vector<T2>& x, const T1& beta, Vector<T2>& y)
  {
    Mlt(beta, y);
    A.MltAdd(alpha, SeldonNoTrans, x, y);
  }


  //! y = beta*y + alpha*A*x
  template<class T0, class T1, class T, class T2, class Allocator>
  inline void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
			   const Matrix<T, General, ArrowCol, Allocator>& A,
			   const Vector<T2>& x, const T1& beta, Vector<T2>& y)
  {
    Mlt(beta, y);
    A.MltAdd(alpha, trans, x, y);
  }


  //! y = A*x
  template<class T, class Allocator, class T1>
  inline void MltVector(const Matrix<T, General, ArrowCol, Allocator>& A,
			const Vector<T1>& x, Vector<T1>& y)
  {
    T1 zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    y.Fill(zero);
    A.MltAdd(one, SeldonNoTrans, x, y);
  }


  //! y = A*x
  template<class T, class Allocator, class T1>
  inline void MltVector(const SeldonTranspose& trans,
			const Matrix<T, General, ArrowCol, Allocator>& A,
			const Vector<T1>& x, Vector<T1>& y)
  {
    T1 zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    y.Fill(zero);
    A.MltAdd(one, trans, x, y);
  }


  //! A = alpha*A
  template<class T0, class T1, class Allocator>
  inline void MltScalar(const T0& alpha,
			Matrix<T1, General, ArrowCol, Allocator>& A)
  {
    A *= alpha;
  }
  
  
  //! displays matrix
  template<class T, class Allocator>
  inline ostream& operator<<(ostream& out,
			     const Matrix<T, General, ArrowCol, Allocator>& A)
  {
    A.WriteText(out);
    return out;
  }
  
}

#define SELDON_FILE_BAND_MATRIX_INLINE_CXX
#endif
