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


#ifndef SELDON_FILE_MATRIX_ARRAY_SPARSE_INLINE_CXX

#include "Matrix_ArraySparse.hxx"

namespace Seldon
{

  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocator>::Matrix_ArraySparse()
    : VirtualMatrix<T>(), val_()
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Matrix_ArraySparse(int i, int j) :
    VirtualMatrix<T>(i, j), val_(Storage::GetFirst(i, j))
  {
  }


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocat>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocat>::~Matrix_ArraySparse()
  {
    this->m_ = 0;
    this->n_ = 0;
  }


  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix is
    empty (0 by 0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->val_.Clear();
    this->m_ = 0;
    this->n_ = 0;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetDataSize()
    const
  {
    return GetNonZeros();
  }


  //! Returns (row or column) indices of non-zero entries in row
  /*!
    \param[in] i row (or column) number.
    \return The array of column (or row) indices of non-zero entries
    of row (or column) i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetIndex(int i)
    const
  {
    return val_(i).GetIndex();
  }


  //! Returns values of non-zero entries of a row/column.
  /*!
    \param[in] i row (or column) number.
    \return The array of values of non-zero entries of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T*
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetData(int i) const
  {
    return val_(i).GetData();
  }


  //! Returns values of non-zero entries.
  /*!
    \return Array of sparse rows
    There is a different array for each row/column.
  */
  template <class T, class Prop, class Storage, class Allocat>
  inline Vector<T, VectSparse, Allocat>*
  Matrix_ArraySparse<T, Prop, Storage, Allocat>::GetData() const
  {
    return val_.GetData();
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/
  
  
#ifdef SELDON_WITH_MODIFIABLE_PARENTHESIS_OPERATOR
  template <class T, class Prop, class Storage, class Allocator>
  inline T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::operator() (int i, int j)  
  {
    return Get(i, j);
  }
#endif
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const T
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::operator() (int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_ArraySparse");
#endif

    return this->val_(Storage::GetFirst(i, j))(Storage::GetSecond(i, j));
  }


  //! Access to element (i, j).
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_ArraySparse");
#endif

    return this->val_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j));
  }

  
  //! Access to element (i, j).
  /*!
    Returns the value of element (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_ArraySparse");
#endif

    return this->val_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j));
  }
  
  
  //! Access method.
  /*! Returns the value of element (\a i, \a j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_ArraySparse");
#endif

    return
      this->val_(Storage::GetFirst(i, j)).Val(Storage::GetSecond(i, j));
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_ArraySparse");
#endif

    return
      this->val_(Storage::GetFirst(i, j)).Val(Storage::GetSecond(i, j));
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& x)
  {
    this->Get(i, j) = x;
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::SetEntry(int i, int j, const T& x)
  {
    this->Get(i, j) = x;
  }


  //! Returns j-th non-zero value of row/column i.
  /*!
    \param[in] i row/column number.
    \param[in] j local number.
    \return j-th non-zero entry of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const T& Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Value (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_(i).GetM(), "Matrix_ArraySparse");
#endif

    return val_(i).Value(j);
  }


  //! Returns j-th non-zero value of row/column i.
  /*!
    \param[in] i row/column number.
    \param[in] j local number.
    \return j-th non-zero entry of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Value (int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_(i).GetM(), "Matrix_ArraySparse");
#endif

    return val_(i).Value(j);
  }


  //! Returns column/row number of j-th non-zero value of row/column i.
  /*!
    \param[in] i row/column number.
    \param[in] j local number.
    \return Column/row number of j-th non-zero value of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator> inline
  int Matrix_ArraySparse<T, Prop, Storage, Allocator>::Index(int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_(i).GetM(), "Matrix_ArraySparse");
#endif

    return val_(i).Index(j);
  }


  //! Returns column/row number of j-th non-zero value of row/column i.
  /*!
    \param[in] i row/column number.
    \param[in] j local number.
    \return Column/row number of j-th non-zero value of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator> inline
  int& Matrix_ArraySparse<T, Prop, Storage, Allocator>::Index(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, Storage::GetFirst(this->m_, this->n_),
		this->val_(i).GetM(), "Matrix_ArraySparse");
#endif

    return val_(i).Index(j);
  }


  //! Redefines a row/column of the matrix
  /*!
    \param[in] i row/col number
    \param[in] n number of non-zero entries in the row
    \param[in] val values
    \param[in] ind column numbers
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  SetData(int i, int n, T* val, int* ind)
  {
    val_(i).SetData(n, val, ind);
  }


  //!  Clears a row without releasing memory.
  /*!
    On exit, the row is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Nullify(int i)
  {
    val_(i).Nullify();
  }


  //! Redefines the matrix.
  /*!
    \param[in] m new number of rows.
    \param[in] n new number of columns.
    \param[in] val array of sparse rows/columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  SetData(int m, int n, Vector<T, VectSparse, Allocator>* val)
  {
    this->m_ = m;
    this->n_ = n;
    val_.SetData(Storage::GetFirst(m, n), val);
  }


  //!  Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->n_ = 0;
    val_.Nullify();
  }


#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Tcplx>& x, const Vector<Tcplx>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_ArraySparse<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return IsSymmetricMatrix(static_cast<const Matrix<T, Prop, Storage,
			     Allocator>& >(*this));
  }
#endif


  /////////////////////////////
  // MATRIX<ARRAY_COLSPARSE> //
  /////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSparse, Allocator>::Matrix():
    Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>(i, j)
  {
  }


  //! Clears column i.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::ClearColumn(int i)
  {
    this->val_(i).Clear();
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Alloc> inline
  void Matrix<T, Prop, ArrayColSparse, Alloc>::ReallocateColumn(int i,int j)
  {
    this->val_(i).Reallocate(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::ResizeColumn(int i,int j)
  {
    this->val_(i).Resize(j);
  }


  //! Swaps two columns.
  /*!
    \param[in] i first column number.
    \param[in] j second column number.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::SwapColumn(int i,int j)
  {
    Swap(this->val_(i), this->val_(j));
  }


  //! Sets row numbers of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \param[in] new_index new row numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  ReplaceIndexColumn(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_(i).GetM(); j++)
      this->val_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \return The number of non-zero entries of the column i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColSparse, Allocator>::
  GetColumnSize(int i) const
  {
    return this->val_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::PrintColumn(int i) const
  {
    this->val_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::AssembleColumn(int i)
  {
    this->val_(i).Assemble();
  }


  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    this->val_(j).AddInteraction(i, val);
  }


  /////////////////////////////
  // MATRIX<ARRAY_ROWSPARSE> //
  /////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSparse, Allocator>::Matrix():
    Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>(i, j)
  {
  }


  //! Clears a row
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::ClearRow(int i)
  {
    this->val_(i).Clear();
  }


  //! Changes the size of a row.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries of the row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  ReallocateRow(int i, int j)
  {
    this->val_(i).Reallocate(j);
  }


  //! Changes the size of a row.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries of the row.
    \note Data is kept.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::ResizeRow(int i, int j)
  {
    this->val_(i).Resize(j);
  }


  //! Swaps two rows
  /*!
    \param[in] i first row number.
    \param[in] j second row number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::SwapRow(int i,int j)
  {
    Swap(this->val_(i), this->val_(j));
  }


  //! Sets column numbers of non-zero entries of a row.
  /*!
    \param[in] i column number.
    \param[in] new_index new column numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  ReplaceIndexRow(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_(i).GetM(); j++)
      this->val_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \return The number of non-zero entries of the row i.
  */
  template <class T, class Prop, class Allocator> inline
  int Matrix<T, Prop, ArrayRowSparse, Allocator>::GetRowSize(int i) const
  {
    return this->val_(i).GetSize();
  }


  //! Displays non-zero values of a row.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::PrintRow(int i) const
  {
    this->val_(i).Print();
  }


  //! Assembles a row.
  /*!
    \param[in] i row number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::AssembleRow(int i)
  {
    this->val_(i).Assemble();
  }

  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    this->val_(i).AddInteraction(j, val);
  }

  
  ////////////////////////////////
  // MATRIX<ARRAY_COLSYMSPARSE> //
  ////////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSymSparse, Allocator>::Matrix():
    Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>(i, j)
  {
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


#ifdef SELDON_WITH_MODIFIABLE_PARENTHESIS_OPERATOR
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::operator() (int i, int j)  
  {
    return Get(i, j);
  }
#endif
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::operator() (int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i < j)
      return this->val_(j)(i);

    return this->val_(i)(j);
  }


  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i < j)
      return this->val_(j).Get(i);
    
    return this->val_(i).Get(j);
  }

  
  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i < j)
      return this->val_(j).Get(i);
    
    return this->val_(i).Get(j);
  }

  
  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i > j)
      return this->val_(i).Val(j);
    
    return this->val_(j).Val(i);
  }

  
  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i > j)
      return this->val_(i).Val(j);
    
    return this->val_(j).Val(i);
  }
  
  
  //! Sets element (i, j) of the matrix
  /*!
    \param i row index.
    \param j column index.
    \param x A(i, j) = x
  */  
  template <class T, class Prop, class Allocator>
  inline void
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Set(int i, int j, const T& x)
  {
    if (i < j)
      this->val_(j).Get(i) = x;
    else
      this->val_(i).Get(j) = x;
  }


  //! Sets element (i, j) of the matrix
  /*!
    \param i row index.
    \param j column index.
    \param x A(i, j) = x
  */  
  template <class T, class Prop, class Allocator>
  inline void
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::SetEntry(int i, int j, const T& x)
  {
    if (i < j)
      this->val_(j).Get(i) = x;
    else
      this->val_(i).Get(j) = x;
  }
  
  
  //! Clears a column.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSymSparse, Allocator>::ClearColumn(int i)
  {
    this->val_(i).Clear();
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ReallocateColumn(int i, int j)
  {
    this->val_(i).Reallocate(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ResizeColumn(int i, int j)
  {
    this->val_(i).Resize(j);
  }


  //! Swaps two columns.
  /*!
    \param[in] i first column number.
    \param[in] j second column number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  SwapColumn(int i, int j)
  {
    Swap(this->val_(i), this->val_(j));
  }


  //! Sets row numbers of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \param[in] new_index new row numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ReplaceIndexColumn(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_(i).GetM(); j++)
      this->val_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \return The number of non-zero entries of the column i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  GetColumnSize(int i) const
  {
    return this->val_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  PrintColumn(int i) const
  {
    this->val_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AssembleColumn(int i)
  {
    this->val_(i).Assemble();
  }


  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      this->val_(j).AddInteraction(i, val);
  }

  
  ////////////////////////////////
  // MATRIX<ARRAY_ROWSYMSPARSE> //
  ////////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Matrix():
    Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>(i, j)
  {
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


#ifdef SELDON_WITH_MODIFIABLE_PARENTHESIS_OPERATOR
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::operator() (int i, int j)  
  {
    return Get(i, j);
  }
#endif

  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::operator() (int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i < j)
      return this->val_(i)(j);

    return this->val_(j)(i);
  }


  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i < j)
      return this->val_(i).Get(j);
    
    return this->val_(j).Get(i);
  }

  
  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif

    if (i < j)
      return this->val_(i).Get(j);
    
    return this->val_(j).Get(i);
  }

  
  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i > j)
      return this->val_(j).Val(i);
    
    return this->val_(i).Val(j);
  }

  
  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix");
#endif
    
    if (i > j)
      return this->val_(j).Val(i);
    
    return this->val_(i).Val(j);
  }
  
  
  //! Sets element (i, j) of the matrix
  /*!
    \param i row index.
    \param j column index.
    \param x A(i, j) = x
  */  
  template <class T, class Prop, class Allocator>
  inline void
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Set(int i, int j, const T& x)
  {
    if (i < j)
      this->val_(i).Get(j) = x;
    else
      this->val_(j).Get(i) = x;
  }


  //! Sets element (i, j) of the matrix
  /*!
    \param i row index.
    \param j column index.
    \param x A(i, j) = x
  */  
  template <class T, class Prop, class Allocator>
  inline void
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::SetEntry(int i, int j, const T& x)
  {
    if (i < j)
      this->val_(i).Get(j) = x;
    else
      this->val_(j).Get(i) = x;
  }
  
  
  //! Clears a row.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::ClearRow(int i)
  {
    this->val_(i).Clear();
  }


  //! Reallocates row i.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries in the row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ReallocateRow(int i,int j)
  {
    this->val_(i).Reallocate(j);
  }


  //! Reallocates row i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the row.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ResizeRow(int i,int j)
  {
    this->val_(i).Resize(j);
  }


  //! Swaps two rows.
  /*!
    \param[in] i first row number.
    \param[in] j second row number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  SwapRow(int i,int j)
  {
    Swap(this->val_(i), this->val_(j));
  }


  //! Sets column numbers of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \param[in] new_index new column numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ReplaceIndexRow(int i,IVect& new_index)
  {
    for (int j = 0; j < this->val_(i).GetM(); j++)
      this->val_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \return The number of non-zero entries of the row i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayRowSymSparse, Allocator>::GetRowSize(int i)
    const
  {
    return this->val_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::PrintRow(int i)
    const
  {
    this->val_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>
  ::AssembleRow(int i)
  {
    this->val_(i).Assemble();
  }


  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      this->val_(i).AddInteraction(j, val);
  }


  template <class T, class Prop, class Allocator>
  inline ostream& operator <<(ostream& out,
			      const Matrix<T, Prop, ArrayRowSparse, Allocator>& A)
  {
    A.WriteText(out);

    return out;
  }


  template <class T, class Prop, class Allocator>
  inline ostream& operator <<(ostream& out,
			      const Matrix<T, Prop, ArrayColSparse, Allocator>& A)
  {
    A.WriteText(out);

    return out;
  }


  template <class T, class Prop, class Allocator>
  inline ostream& operator <<(ostream& out,
			      const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    A.WriteText(out);

    return out;
  }


  template <class T, class Prop, class Allocator>
  inline ostream& operator <<(ostream& out,
			      const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A)
  {
    A.WriteText(out);

    return out;
  }
  
} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_SPARSE_INLINE_CXX
#endif
