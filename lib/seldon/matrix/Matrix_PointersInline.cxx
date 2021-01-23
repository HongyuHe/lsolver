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


#ifndef SELDON_FILE_MATRIX_POINTERS_INLINE_CXX

#include "Matrix_Pointers.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Pointers<T, Prop, Storage, Allocator>::Matrix_Pointers():
    Matrix_Base<T, Allocator>()
  {
    me_ = NULL;
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Pointers<T, Prop, Storage, Allocator>::~Matrix_Pointers()
  {
    this->Clear();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of rows multiplied by the number of columns
    because the matrix is full.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_Pointers<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return long(this->m_) * long(this->n_);
  }


  //! Returns size of A in bytes used to store the matrix.
  template<class T, class Prop, class Storage, class Allocator>
  inline size_t Matrix_Pointers<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + size_t(this->GetDataSize())*sizeof(T);
    taille += size_t(Storage::GetFirst(this->m_, this->n_))*sizeof(pointer);
    return taille;
  }
  

  //! Returns the pointer 'me_'.
  /*! Returns the pointer 'me_' that defines an array pointing to the first
    row or column elements, so that 'me_[1]' points to the first element of
    the second row or column.
    \return The pointer 'me_'.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename
  Matrix_Pointers<T, Prop, Storage, Allocator>::pointer*
  Matrix_Pointers<T, Prop, Storage, Allocator>::GetMe() const
  {
    return me_;
  }


  //! Returns the leading dimension.
  /*!
    \return The leading dimension.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Pointers<T, Prop, Storage, Allocator>::GetLD() const
  {
    return Storage::GetSecond(this->m_, this->n_);
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Returns a pointer to a data element.
  /*!
    \param i index along dimension #1.
    \param j index along dimension #2.
    \return A pointer to the data element.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::pointer
  Matrix_Pointers<T, Prop, Storage, Allocator>::GetDataPointer(int i, int j)
    const
  {
    long lgth = Storage::GetSecond(this->m_, this->n_);
    return this->data_ + long(Storage::GetFirst(i, j)) * lgth
      + long(Storage::GetSecond(i, j));
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::operator() (int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_Pointers");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>
  ::const_reference Matrix_Pointers<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_Pointers");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_Pointers");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }

  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
    return Val(i, j);
  }
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_Pointers");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    return Val(i, j);
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::operator[] (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->GetDataSize(), "Matrix_Pointers");
#endif

    return this->data_[i];
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::operator[] (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->GetDataSize(), "Matrix_Pointers");
#endif

    return this->data_[i];
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& val)
  {
    this->Val(i, j) = val;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Pointers<T, Prop, Storage, Allocator>&
  Matrix_Pointers<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_Pointers<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Copy(const Matrix_Pointers<T, Prop, Storage, Allocator>& A)
  {
    this->Reallocate(A.GetM(), A.GetN());

    Allocator::memorycpy(this->data_, A.GetData(), this->GetDataSize());
  }


#ifdef SELDON_WITH_VIRTUAL
  //! returns true if the matrix is symmetric
  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_Pointers<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return false;
  }


  //! adds val to A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const T& val)
  {
    // not implemented since T can be a complex structure
    abort();
  }


  //! Adds several values on row i
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::AddInteractionRow(int i, int n, const Vector<int>& col,
		      const Vector<T>& val, bool sorted)
  {
    // not implemented since T can be a complex structure
    abort();
  }


  //! Adds several values on column i
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::AddInteractionColumn(int i, int n, const Vector<int>& row,
			 const Vector<T>& val, bool sorted)
  {
    // not implemented since T can be a complex structure
    abort();
  }


  //! sets to zero the row i
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::ClearRow(int i)
  {
    // not implemented since T can be a complex structure
    abort();
  }

  
  //! Computes y = beta + alpha A x.
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    // not implemented since T can be a complex structure
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = beta + alpha A x.
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    // not implemented since T can be a complex structure
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = beta + alpha A x.
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose&,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    // not implemented since T can be a complex structure
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = beta + alpha A x.
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose&,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    // not implemented since T can be a complex structure
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = A x
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    // not implemented since T can be a complex structure
    throw Undefined("MltVector", "Not implemented");
  }


  //! Computes y = A x
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    // not implemented since T can be a complex structure
    throw Undefined("MltVector", "Not implemented");
  }
  

  //! Computes y = A^T x
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::MltVector(const SeldonTranspose&,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    // not implemented since T can be a complex structure
    throw Undefined("MltVector with transpose", "Not implemented");
  }


  //! Computes y = A^T x
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::MltVector(const SeldonTranspose&,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    // not implemented since T can be a complex structure
    throw Undefined("MltVector with transpose", "Not implemented");
  }

#endif


  //////////////////////
  // MATRIX<COLMAJOR> //
  //////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajor, Allocator>::Matrix():
    Matrix_Pointers<T, Prop, ColMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i by j full column-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajor, Allocator>::Matrix(int i, int j):
    Matrix_Pointers<T, Prop, ColMajor, Allocator>(i, j)
  {
  }


  //! Copy constructor.
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajor, Allocator>
  ::Matrix(const Matrix<T, Prop, ColMajor, Allocator>& A):
    Matrix_Pointers<T, Prop, ColMajor, Allocator>(A)
  {
  }


  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColMajor, Allocator>&
  Matrix<T, Prop, ColMajor, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }

  
  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajor, Allocator>&
  Matrix<T, Prop, ColMajor, Allocator>
  ::operator= (const Matrix<T, Prop, ColMajor, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Multiplies the matrix by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Prop, class Allocator> template<class T0>
  inline Matrix<T, Prop, ColMajor, Allocator>&
  Matrix<T, Prop, ColMajor, Allocator>::operator*= (const T0& alpha)
  {
    for (long i = 0; i < long(this->m_) * long(this->n_); i++)
      this->data_[i] *= alpha;

    return *this;
  }


  //////////////////////
  // MATRIX<ROWMAJOR> //
  //////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajor, Allocator>::Matrix():
    Matrix_Pointers<T, Prop, RowMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i by j full row-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajor, Allocator>::Matrix(int i, int j):
    Matrix_Pointers<T, Prop, RowMajor, Allocator>(i, j)
  {
  }


  //! Copy constructor.
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajor, Allocator>
  ::Matrix(const Matrix<T, Prop, RowMajor, Allocator>& A):
    Matrix_Pointers<T, Prop, RowMajor, Allocator>(A)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowMajor, Allocator>&
  Matrix<T, Prop, RowMajor, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }
  
  
  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajor, Allocator>&
  Matrix<T, Prop, RowMajor, Allocator>
  ::operator= (const Matrix<T, Prop, RowMajor, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Multiplies the matrix by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Prop, class Allocator> template<class T0>
  inline Matrix<T, Prop, RowMajor, Allocator>&
  Matrix<T, Prop, RowMajor, Allocator>::operator*= (const T0& alpha)
  {
    for (long i = 0; i < long(this->m_)*long(this->n_); i++)
      this->data_[i] *= alpha;

    return *this;
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_POINTERS_INLINE_CXX
#endif
