// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2003-2011 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_SYMSPARSE_INLINE_CXX

#include "Matrix_SymSparse.hxx"

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
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::Matrix_SymSparse():
    Matrix_Base<T, Allocator>()
  {
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymSparse(int i, int j): Matrix_Base<T, Allocator>()
  {
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;

    Reallocate(i, j);
  }


  //! Constructor.
  /*! Builds a sparse matrix of size i x j , with nz non-zero (stored)
    elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized. Indices of non-zero entries
    are not initialized either.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Matrix_SymSparse(int i, int j, long nz):
    Matrix_Base<T, Allocator>()
  {
    this->nz_ = 0;
    ind_ = NULL;
    ptr_ = NULL;

    Reallocate(i, j, nz);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::~Matrix_SymSparse()
  {
    this->Clear();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory (the number
    of non-zero elements).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymSparse<T, Prop, Storage, Allocator>::GetNonZeros() const
  {
    return nz_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory (the number
    of non-zero elements).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymSparse<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return nz_;
  }


  //! Returns (row or column) start indices.
  /*!
    Returns the array ('ptr_') of start indices.
    \return The array of start indices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long* Matrix_SymSparse<T, Prop, Storage, Allocator>::GetPtr() const
  {
    return ptr_;
  }


  //! Returns (row or column) indices of non-zero entries.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries. This array defines non-zero entries
    indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_SymSparse<T, Prop, Storage, Allocator>::GetInd() const
  {
    return ind_;
  }


  //! Returns the length of the array of start indices.
  /*!
    \return The length of the array of start indices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_SymSparse<T, Prop, Storage, Allocator>::GetPtrSize() const
  {
    return (this->m_ + 1);
  }


  //! Returns the length of the array of (column or row) indices.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored). This array defines non-zero
    entries indices (that are stored) if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymSparse<T, Prop, Storage, Allocator>::GetIndSize() const
  {
    return nz_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access method.
  /*! Returns reference to element (\a i, \a j) 
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Matrix_SymSparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_SymSparse<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    return Val(i, j);
  }
  
  
  //! Add a value to a non-zero entry.
  /*! This function adds \a val to the element (\a i, \a j), provided that
    this element is already a non-zero entry. Otherwise 
    a non-zero entry is inserted equal to \a val.
    \param[in] i row index.
    \param[in] j column index.
    \param[in] val value to be added to the element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      Get(i, j) += val;
  }


  //! Adds values to several non-zero entries on a given row
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::AddInteractionRow(int i, int nb, const Vector<int>& col,
		      const Vector<T>& val, bool sorted)
  {
    throw Undefined("AddInteractionRow", "Not implemented");
  }


  //! Sets an element (i, j) to a value
  /*! This function sets \a val to the element (\a i, \a j)
    \param[in] i row index.
    \param[in] j column index.
    \param[in] val A(i, j) = val
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }


  //! Sets an element (i, j) to a value
  /*! This function sets \a val to the element (\a i, \a j)
    \param[in] i row index.
    \param[in] j column index.
    \param[in] val A(i, j) = val
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::SetEntry(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>&
  Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_SymSparse<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Tcplx>& x, const Vector<Tcplx>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_SymSparse<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return true;
  }
#endif

  
  //////////////////////////
  // MATRIX<COLSYMSPARSE> //
  //////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymSparse, Allocator>::Matrix():
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>(i, j, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with nz non-zero (stored) elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymSparse, Allocator>::Matrix(int i, int j, long nz):
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>(i, j, nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row start indices.
    \param ind column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix<T, Prop, ColSymSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& values,
	 Vector<long, Storage1, Allocator1>& ptr,
	 Vector<int, Storage2, Allocator2>& ind):
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>(i, j, values, ptr, ind)
  {
  }



  //////////////////////////
  // MATRIX<ROWSYMSPARSE> //
  //////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymSparse, Allocator>::Matrix():
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>(i, j, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with nz non-zero (stored) elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymSparse, Allocator>::Matrix(int i, int j, long nz):
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>(i, j, nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr column start indices.
    \param ind row indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix<T, Prop, RowSymSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& values,
	 Vector<long, Storage1, Allocator1>& ptr,
	 Vector<int, Storage2, Allocator2>& ind):
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>(i, j, values, ptr, ind)
  {
  }
    
  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMSPARSE_INLINE_CXX
#endif
