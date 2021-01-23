// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_INLINE_CXX

#include "Matrix_ArrayComplexSparse.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_ArrayComplexSparse()
    : VirtualMatrix<T>(), val_real_(), val_imag_()
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_ArrayComplexSparse(int i, int j): VirtualMatrix<T>(i, j),
    val_real_(Storage::GetFirst(i, j)), val_imag_(Storage::GetFirst(i, j))
  {
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ~Matrix_ArrayComplexSparse()
  {
    Clear();
  }


  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix is
    empty (0 by 0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->m_ = 0;
    this->n_ = 0;
    val_real_.Clear();
    val_imag_.Clear();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory (real+imaginary part).
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetNonZeros() const
  {
    return (GetRealNonZeros()+GetImagNonZeros());
  }


  //! Returns the number of elements stored in memory (real part).
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealDataSize() const
  {
    return GetRealNonZeros();
  }


  //! Returns the number of elements stored in memory (imaginary part).
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagDataSize() const
  {
    return GetImagNonZeros();
  }


  //! Returns the number of elements stored in memory (real+imaginary part).
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetDataSize() const
  {
    return (GetRealNonZeros()+GetImagNonZeros());
  }


  //! Returns column indices of non-zero entries in row (real part).
  /*!
    \param[in] i row number.
    \return The array of column indices of non-zero entries
    of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealInd(int i) const
  {
    return val_real_(i).GetIndex();
  }


  //! Returns values of non-zero entries of a row (real part).
  /*!
    \param[in] i row number.
    \return The array of values of non-zero entries of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type*
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::GetRealData(int i)
    const
  {
    return val_real_(i).GetData();
  }


  //! Returns column indices of non-zero entries in row (imaginary part).
  /*!
    \param[in] i row number.
    \return the array of column indices of non-zero entries
    of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagInd(int i) const
  {
    return val_imag_(i).GetIndex();
  }


  //! Returns values of non-zero entries of a row (imaginary part).
  /*!
    \param[in] i row number.
    \return The array of values of non-zero entries of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type*
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::GetImagData(int i)
    const
  {
    return val_imag_(i).GetData();
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename 
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::entry_type
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::operator()
    (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    return entry_type(this->val_real_(Storage::GetFirst(i, j))
		      (Storage::GetSecond(i, j)),
		      this->val_imag_(Storage::GetFirst(i, j))
		      (Storage::GetSecond(i, j)) );
  }


  //! Unavailable access method.
  /*! This method is declared for consistency with other classes, but it is
    not defined because no reference can possibly be returned.
    \param[in] i row index.
    \param[in] j column index.
    \return Raises an exception.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::entry_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Val(int i, int j)
  {
    throw Undefined("Matrix_ArrayComplexSparse::Val(int i, int j)");
  }
  
  
  //! Unavailable access method.
  /*! This method is declared for consistency with other classes, but it is
    not defined because no reference can possibly be returned.
    \param[in] i row index.
    \param[in] j column index.
    \return Raises an exception.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::entry_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Val(int i, int j) const
  {
    throw Undefined("Matrix_ArrayComplexSparse::Val(int i, int j)");
  }

  
  //! Unavailable access method.
  /*! This method is declared for consistency with other classes, but it is
    not defined because no reference can possibly be returned.
    \param[in] i row index.
    \param[in] j column index.
    \return Raises an exception.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::entry_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Get(int i, int j)
  {
    throw Undefined("Matrix_ArrayComplexSparse::Get(int i, int j)");
  }
  
  
  //! Unavailable access method.
  /*! This method is declared for consistency with other classes, but it is
    not defined because no reference can possibly be returned.
    \param[in] i row index.
    \param[in] j column index.
    \return Raises an exception.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::entry_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Get(int i, int j) const
  {
    throw Undefined("Matrix_ArrayComplexSparse::Get(int i, int j)");
  }

  
    //! Returns acces to real part of A(i, j)
  /*! 
    \param[in] i row index.
    \param[in] j column index.
    \return real part of A(i, j)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::ValReal(int i, int j)
  {
    return val_real_(Storage::GetFirst(i, j)).Val(Storage::GetSecond(i, j));
  }

  
  //! Returns acces to real part of A(i, j)
  /*! 
    \param[in] i row index.
    \param[in] j column index.
    \return real part of A(i, j)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::ValReal(int i, int j) const
  {
    return val_real_(Storage::GetFirst(i, j)).Val(Storage::GetSecond(i, j));
  }
  
  
  //! Returns acces to imaginary part of A(i, j)
  /*! 
    \param[in] i row index.
    \param[in] j column index.
    \return imaginary part of A(i, j)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::ValImag(int i, int j)
  {
    return val_imag_(Storage::GetFirst(i, j)).Val(Storage::GetSecond(i, j));
  }

  
  //! Returns acces to imaginary part of A(i, j)
  /*! 
    \param[in] i row index.
    \param[in] j column index.
    \return imaginary part of A(i, j)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::ValImag(int i, int j) const
  {
    return val_imag_(Storage::GetFirst(i, j)).Val(Storage::GetSecond(i, j));
  }
  
  
  //! Returns acces to real part of A(i, j)
  /*! 
    \param[in] i row index.
    \param[in] j column index.
    \return real part of A(i, j)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::GetReal(int i, int j)
  {
    return val_real_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j));
  }

  
  //! Returns acces to real part of A(i, j)
  /*! 
    \param[in] i row index.
    \param[in] j column index.
    \return real part of A(i, j)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::GetReal(int i, int j) const
  {
    return val_real_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j));
  }
  
  
  //! Returns acces to imaginary part of A(i, j)
  /*! 
    \param[in] i row index.
    \param[in] j column index.
    \return imaginary part of A(i, j)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::GetImag(int i, int j)
  {
    return val_imag_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j));
  }

  
  //! Returns acces to imaginary part of A(i, j)
  /*! 
    \param[in] i row index.
    \param[in] j column index.
    \return imaginary part of A(i, j)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::GetImag(int i, int j) const
  {
    return val_imag_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j));
  }
  
  
  //! Returns j-th non-zero value of row i (real part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ValueReal(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_real_(i).GetM(), "Matrix_ArrayComplexSparse");
#endif

    return val_real_(i).Value(j);
  }


  //! Returns j-th non-zero value of row i (real part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ValueReal(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_real_(i).GetM(), "Matrix_ArrayComplexSparse");
#endif

    return val_real_(i).Value(j);
  }


  //! Returns column number of j-th non-zero value of row i (real part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return Column number of j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  IndexReal(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_real_(i).GetM(), "Matrix_ArrayComplexSparse");
#endif

    return val_real_(i).Index(j);
  }


  //! Returns column number of j-th non-zero value of row i (real part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return column number of j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int& Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  IndexReal(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_real_(i).GetM(), "Matrix_ArrayComplexSparse");
#endif

    return val_real_(i).Index(j);
  }


  //! Returns j-th non-zero value of row i (imaginary part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ValueImag(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_imag_(i).GetM(), "Matrix_ArrayComplexSparse");
#endif

    return val_imag_(i).Value(j);
  }


  //! Returns j-th non-zero value of row i (imaginary part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ValueImag (int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_imag_(i).GetM(), "Matrix_ArrayComplexSparse");
#endif

    return val_imag_(i).Value(j);
  }


  //! Returns column number of j-th non-zero value of row i (imaginary part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return Column number of j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  IndexImag(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_imag_(i).GetM(), "Matrix_ArrayComplexSparse");
#endif

    return val_imag_(i).Index(j);
  }


  //! Returns column number of j-th non-zero value of row i (imaginary part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return column number of j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int& Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  IndexImag(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_imag_(i).GetM(), "Matrix_ArrayComplexSparse");
#endif

    return val_imag_(i).Index(j);
  }


  //! Redefines a row/column of the matrix
  /*!
    \param[in] i row/col number
    \param[in] n number of non-zero entries in the row
    \param[in] val values
    \param[in] ind column numbers
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  SetRealData(int i, int n, value_type* val, int* ind)
  {
    val_real_(i).SetData(n, val, ind);
  }


  //! Redefines a row/column of the matrix
  /*!
    \param[in] i row/col number
    \param[in] n number of non-zero entries in the row
    \param[in] val values
    \param[in] ind column numbers
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  SetImagData(int i, int n, value_type* val, int* ind)
  {
    val_imag_(i).SetData(n, val, ind);
  }


  //!  Clears a row without releasing memory.
  /*!
    On exit, the row is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::NullifyReal(int i)
  {
    val_real_(i).Nullify();
  }


  //!  Clears a row without releasing memory.
  /*!
    On exit, the row is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::NullifyImag(int i)
  {
    val_imag_(i).Nullify();
  }


  //! Redefines the real part of the matrix.
  /*!
    \param[in] m new number of rows.
    \param[in] n new number of columns.
    \param[in] val array of sparse rows/columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  SetRealData(int m, int n, Vector<value_type, VectSparse, Allocator>* val)
  {
    this->m_ = m;
    this->n_ = n;
    val_real_.SetData(Storage::GetFirst(m, n), val);
  }


  //! Redefines the imaginary part of the matrix.
  /*!
    \param[in] m new number of rows.
    \param[in] n new number of columns.
    \param[in] val array of sparse rows/columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  SetImagData(int m, int n, Vector<value_type, VectSparse, Allocator>* val)
  {
    this->m_ = m;
    this->n_ = n;
    val_imag_.SetData(Storage::GetFirst(m, n), val);
  }


  //!  Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::NullifyReal()
  {
    this->m_ = 0;
    this->n_ = 0;
    val_real_.Nullify();
  }


  //!  Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::NullifyImag()
  {
    this->m_ = 0;
    this->n_ = 0;
    val_imag_.Nullify();
  }


#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Tcplx>& x, const Vector<Tcplx>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return IsSymmetricMatrix(static_cast<const Matrix<T, Prop, Storage,
			     Allocator>& >(*this));
  }
#endif

  
  ////////////////////////////////////
  // MATRIX<ARRAY_COLCOMPLEXSPARSE> //
  ////////////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColComplexSparse, Allocator>::Matrix()  :
    Matrix_ArrayComplexSparse<T, Prop, ArrayColComplexSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::Matrix(int i, int j) :
    Matrix_ArrayComplexSparse<T, Prop, ArrayColComplexSparse, Allocator>(i, j)
  {
  }


  //! Clears column i.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::ClearRealColumn(int i)
  {
    this->val_real_(i).Clear();
  }


  //! Clears column i.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::ClearImagColumn(int i)
  {
    this->val_imag_(i).Clear();
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Alloc> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Alloc>
  ::ReallocateRealColumn(int i, int j)
  {
    this->val_real_(i).Reallocate(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Alloc> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Alloc>
  ::ReallocateImagColumn(int i, int j)
  {
    this->val_imag_(i).Reallocate(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::ResizeRealColumn(int i, int j)
  {
    this->val_real_(i).Resize(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::ResizeImagColumn(int i, int j)
  {
    this->val_imag_(i).Resize(j);
  }


  //! Swaps two columns.
  /*!
    \param[in] i first column number.
    \param[in] j second column number.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::SwapRealColumn(int i, int j)
  {
    Swap(this->val_real_(i), this->val_real_(j));
  }


  //! Swaps two columns.
  /*!
    \param[in] i first column number.
    \param[in] j second column number.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::SwapImagColumn(int i, int j)
  {
    Swap(this->val_imag_(i), this->val_imag_(j));
  }


  //! Sets row numbers of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \param[in] new_index new row numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColComplexSparse, Allocator>::
  ReplaceRealIndexColumn(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_real_(i).GetM(); j++)
      this->val_real_(i).Index(j) = new_index(j);
  }


  //! Sets row numbers of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \param[in] new_index new row numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColComplexSparse, Allocator>::
  ReplaceImagIndexColumn(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_imag_(i).GetM(); j++)
      this->val_imag_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \return The number of non-zero entries of the column i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColComplexSparse, Allocator>::
  GetRealColumnSize(int i) const
  {
    return this->val_real_(i).GetSize();
  }


  //! Returns the number of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \return The number of non-zero entries of the column i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColComplexSparse, Allocator>::
  GetImagColumnSize(int i) const
  {
    return this->val_imag_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::PrintRealColumn(int i) const
  {
    this->val_real_(i).Print();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::PrintImagColumn(int i) const
  {
    this->val_imag_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::AssembleRealColumn(int i)
  {
    this->val_real_(i).Assemble();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>
  ::AssembleImagColumn(int i)
  {
    this->val_imag_(i).Assemble();
  }


  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColComplexSparse, Allocator>::
  AddInteraction(int i, int j, const entry_type& val)
  {
    if (real(val) != value_type(0))
      this->val_real_(j).AddInteraction(i, real(val));

    if (imag(val) != value_type(0))
      this->val_imag_(j).AddInteraction(i, imag(val));
  }

  
  ////////////////////////////////////
  // MATRIX<ARRAY_ROWCOMPLEXSPARSE> //
  ////////////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::Matrix()  :
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::Matrix(int i, int j):
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>(i, j)
  {
  }


  //! Clears a row
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::ClearRealRow(int i)
  {
    this->val_real_(i).Clear();
  }


  //! Clears a row
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::ClearImagRow(int i)
  {
    this->val_imag_(i).Clear();
  }


  //! Clears a row
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::ClearRow(int i)
  {
    ClearRealRow(i);
    ClearImagRow(i);
  }


  //! Changes the size of a row.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries of the row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  ReallocateRealRow(int i, int j)
  {
    this->val_real_(i).Reallocate(j);
  }


  //! Changes the size of a row.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries of the row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  ReallocateImagRow(int i, int j)
  {
    this->val_imag_(i).Reallocate(j);
  }


  //! Changes the size of a row.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries of the row.
    \note Data is kept.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::ResizeRealRow(int i, int j)
  {
    this->val_real_(i).Resize(j);
  }


  //! Changes the size of a row.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries of the row.
    \note Data is kept.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::ResizeImagRow(int i, int j)
  {
    this->val_imag_(i).Resize(j);
  }


  //! Swaps two rows
  /*!
    \param[in] i first row number.
    \param[in] j second row number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::SwapRealRow(int i,int j)
  {
    Swap(this->val_real_(i), this->val_real_(j));
  }


  //! Swaps two rows
  /*!
    \param[in] i first row number.
    \param[in] j second row number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::SwapImagRow(int i, int j)
  {
    Swap(this->val_imag_(i), this->val_imag_(j));
  }


  //! Sets column numbers of non-zero entries of a row.
  /*!
    \param[in] i column number.
    \param[in] new_index new column numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  ReplaceRealIndexRow(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_real_(i).GetM(); j++)
      this->val_real_(i).Index(j) = new_index(j);
  }


  //! Sets column numbers of non-zero entries of a row.
  /*!
    \param[in] i column number.
    \param[in] new_index new column numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  ReplaceImagIndexRow(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_imag_(i).GetM(); j++)
      this->val_imag_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \return The number of non-zero entries of the row i.
  */
  template <class T, class Prop, class Allocator> inline
  int Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::GetRealRowSize(int i) const
  {
    return this->val_real_(i).GetSize();
  }


  //! Returns the number of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \return The number of non-zero entries of the row i.
  */
  template <class T, class Prop, class Allocator> inline
  int Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::GetImagRowSize(int i) const
  {
    return this->val_imag_(i).GetSize();
  }


  //! Displays non-zero values of a row.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::PrintRealRow(int i) const
  {
    this->val_real_(i).Print();
  }


  //! Displays non-zero values of a row.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::PrintImagRow(int i) const
  {
    this->val_imag_(i).Print();
  }


  //! Assembles a row.
  /*!
    \param[in] i row number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::AssembleRealRow(int i)
  {
    this->val_real_(i).Assemble();
  }


  //! Assembles a row.
  /*!
    \param[in] i row number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
  ::AssembleImagRow(int i)
  {
    this->val_imag_(i).Assemble();
  }


  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  AddInteraction(int i, int j, const entry_type& val)
  {
    if (real(val) != value_type(0))
      this->val_real_(i).AddInteraction(j, real(val));

    if (imag(val) != value_type(0))
      this->val_imag_(i).AddInteraction(j, imag(val));
  }

  
  ///////////////////////////////////////
  // MATRIX<ARRAY_COLSYMCOMPLEXSPARSE> //
  ///////////////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::Matrix()  :
    Matrix_ArrayComplexSparse<T, Prop, ArrayColSymComplexSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArrayComplexSparse<T, Prop, ArrayColSymComplexSparse, Allocator>(i, j)
  {
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename 
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::entry_type
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::operator() (int i, int j)
    const
  {
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i <= j)
      return entry_type(this->val_real_(j)(i), this->val_imag_(j)(i));

    return entry_type(this->val_real_(i)(j), this->val_imag_(i)(j));
  }

  
  //! Returns access to real part of element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return real part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Allocator::value_type&
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::ValReal(int i, int j)
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i > j)
      return this->val_real_(i).Val(j);
    
    return this->val_real_(j).Val(i);
  }
  
  
  //! Returns access to real part element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return real part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename Allocator::value_type&
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>
  ::ValReal(int i, int j) const
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i > j)
      return this->val_real_(j).Val(i);
    
    return this->val_real_(j).Val(i);
  }

  
  //! Returns access to real part element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return real part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Allocator::value_type&
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::GetReal(int i, int j)
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i <= j)
      return this->val_real_(j).Get(i);
    
    return this->val_real_(i).Get(j);
  }
  
  
  //! Returns access to real part of element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return real part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename Allocator::value_type&
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>
  ::GetReal(int i, int j) const
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i <= j)
      return this->val_real_(j).Get(i);
    
    return this->val_real_(i).Get(j);
  }
  
  
  //! Returns access to imaginary part of element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return imaginary part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Allocator::value_type&
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::ValImag(int i, int j)
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i > j)
      return this->val_imag_(i).Val(j);
    
    return this->val_imag_(j).Val(i);
  }
  
  
  //! Returns access to imaginary part of element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return imaginary part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename Allocator::value_type&
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>
  ::ValImag(int i, int j) const
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    if (i > j)
      return this->val_imag_(i).Val(j);
    
    return this->val_imag_(j).Val(i);
  }

  
  //! Returns access to imaginary part of element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return imaginary part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Allocator::value_type&
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::GetImag(int i, int j)
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i <= j)
      return this->val_imag_(j).Get(i);
    
    return this->val_imag_(i).Get(j);
  }
  
  
  //! Returns access to imaginary part element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return imaginary part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename Allocator::value_type&
  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>
  ::GetImag(int i, int j) const
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i <= j)
      return this->val_imag_(j).Val(i);
    
    return this->val_imag_(i).Val(j);
  }
  
  
  //! Clears a column.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::ClearRealColumn(int i)
  {
    this->val_real_(i).Clear();
  }


  //! Clears a column.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::ClearImagColumn(int i)
  {
    this->val_imag_(i).Clear();
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  ReallocateRealColumn(int i, int j)
  {
    this->val_real_(i).Reallocate(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  ReallocateImagColumn(int i, int j)
  {
    this->val_imag_(i).Reallocate(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  ResizeRealColumn(int i, int j)
  {
    this->val_real_(i).Resize(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  ResizeImagColumn(int i, int j)
  {
    this->val_imag_(i).Resize(j);
  }


  //! Swaps two columns.
  /*!
    \param[in] i first column number.
    \param[in] j second column number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  SwapRealColumn(int i, int j)
  {
    Swap(this->val_real_(i), this->val_real_(j));
  }


  //! Swaps two columns.
  /*!
    \param[in] i first column number.
    \param[in] j second column number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  SwapImagColumn(int i, int j)
  {
    Swap(this->val_imag_(i), this->val_imag_(j));
  }


  //! Sets row numbers of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \param[in] new_index new row numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  ReplaceRealIndexColumn(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_real_(i).GetM(); j++)
      this->val_real_(i).Index(j) = new_index(j);
  }


  //! Sets row numbers of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \param[in] new_index new row numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  ReplaceImagIndexColumn(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_imag_(i).GetM(); j++)
      this->val_imag_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \return The number of non-zero entries of the column i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  GetRealColumnSize(int i) const
  {
    return this->val_real_(i).GetSize();
  }


  //! Returns the number of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \return The number of non-zero entries of the column i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  GetImagColumnSize(int i) const
  {
    return this->val_imag_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  PrintRealColumn(int i) const
  {
    this->val_real_(i).Print();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  PrintImagColumn(int i) const
  {
    this->val_imag_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  AssembleRealColumn(int i)
  {
    this->val_real_(i).Assemble();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  AssembleImagColumn(int i)
  {
    this->val_imag_(i).Assemble();
  }


  ///////////////////////////////////////
  // MATRIX<ARRAY_ROWSYMCOMPLEXSPARSE> //
  ///////////////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::Matrix() :
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowSymComplexSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::Matrix(int i, int j):
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowSymComplexSparse, Allocator>(i, j)
  {
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename 
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::entry_type
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::operator() (int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i <= j)
      return entry_type(this->val_real_(i)(j), this->val_imag_(i)(j));

    return entry_type(this->val_real_(j)(i), this->val_imag_(j)(i));
  }

  
  //! Returns access to real part of element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return real part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Allocator::value_type&
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::ValReal(int i, int j)
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i > j)
      return this->val_real_(j).Val(i);
    
    return this->val_real_(i).Val(j);
  }
  
  
  //! Returns access to real part element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return real part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename Allocator::value_type&
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::ValReal(int i, int j) const
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    if (i > j)
      return this->val_real_(j).Val(i);
    
    return this->val_real_(i).Val(j);
  }

  
  //! Returns access to real part element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return real part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Allocator::value_type&
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::GetReal(int i, int j)
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i <= j)
      return this->val_real_(i).Get(j);
    
    return this->val_real_(j).Get(i);
  }
  
  
  //! Returns access to real part of element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return real part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename Allocator::value_type&
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::GetReal(int i, int j) const
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i <= j)
      return this->val_real_(i).Get(j);
    
    return this->val_real_(j).Get(i);
  }
  
  
  //! Returns access to imaginary part of element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return imaginary part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Allocator::value_type&
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::ValImag(int i, int j)
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    if (i > j)
      return this->val_imag_(j).Val(i);
    
    return this->val_imag_(i).Val(j);
  }
  
  
  //! Returns access to imaginary part of element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return imaginary part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename Allocator::value_type&
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::ValImag(int i, int j) const
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i > j)
      return this->val_imag_(j).Val(i);
    
    return this->val_imag_(i).Val(j);
  }

  
  //! Returns access to imaginary part of element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return imaginary part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Allocator::value_type&
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::GetImag(int i, int j)
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i <= j)
      return this->val_imag_(i).Get(j);
    
    return this->val_imag_(j).Get(i);
  }
  
  
  //! Returns access to imaginary part element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return imaginary part of element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const typename Allocator::value_type&
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::GetImag(int i, int j) const
  {    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i <= j)
      return this->val_imag_(i).Val(j);
    
    return this->val_imag_(j).Val(i);
  }
  
  
  //! Clears a row.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::ClearRealRow(int i)
  {
    this->val_real_(i).Clear();
  }


  //! Clears a row.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::ClearImagRow(int i)
  {
    this->val_imag_(i).Clear();
  }


  //! Clears a row
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::ClearRow(int i)
  {
    ClearRealRow(i);
    ClearImagRow(i);
  }
  

  //! Reallocates row i.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries in the row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  ReallocateRealRow(int i,int j)
  {
    this->val_real_(i).Reallocate(j);
  }


  //! Reallocates row i.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries in the row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  ReallocateImagRow(int i,int j)
  {
    this->val_imag_(i).Reallocate(j);
  }


  //! Reallocates row i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the row.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  ResizeRealRow(int i,int j)
  {
    this->val_real_(i).Resize(j);
  }


  //! Reallocates row i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the row.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  ResizeImagRow(int i,int j)
  {
    this->val_imag_(i).Resize(j);
  }


  //! Swaps two rows.
  /*!
    \param[in] i first row number.
    \param[in] j second row number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  SwapRealRow(int i,int j)
  {
    Swap(this->val_real_(i), this->val_real_(j));
  }


  //! Swaps two rows.
  /*!
    \param[in] i first row number.
    \param[in] j second row number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  SwapImagRow(int i,int j)
  {
    Swap(this->val_imag_(i), this->val_imag_(j));
  }


  //! Sets column numbers of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \param[in] new_index new column numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  ReplaceRealIndexRow(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_real_(i).GetM(); j++)
      this->val_real_(i).Index(j) = new_index(j);
  }


  //! Sets column numbers of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \param[in] new_index new column numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  ReplaceImagIndexRow(int i,IVect& new_index)
  {
    for (int j = 0; j < this->val_imag_(i).GetM(); j++)
      this->val_imag_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \return The number of non-zero entries of the row i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::GetRealRowSize(int i)
    const
  {
    return this->val_real_(i).GetSize();
  }


  //! Returns the number of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \return The number of non-zero entries of the row i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::GetImagRowSize(int i) const
  {
    return this->val_imag_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::PrintRealRow(int i) const
  {
    this->val_real_(i).Print();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::PrintImagRow(int i) const
  {
    this->val_imag_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::AssembleRealRow(int i)
  {
    this->val_real_(i).Assemble();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::AssembleImagRow(int i)
  {
    this->val_imag_(i).Assemble();
  }
  
} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_INLINE_CXX
#endif
