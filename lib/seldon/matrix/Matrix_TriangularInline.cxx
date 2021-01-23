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


#ifndef SELDON_FILE_MATRIX_TRIANGULAR_INLINE_CXX

#include "Matrix_Triangular.hxx"

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
  inline Matrix_Triangular<T, Prop, Storage, Allocator>::Matrix_Triangular():
    Matrix_Base<T, Allocator>()
  {
    me_ = NULL;
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Triangular<T, Prop, Storage, Allocator>::~Matrix_Triangular()
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
  inline long Matrix_Triangular<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return long(this->m_) * long(this->n_);
  }


  //! Returns size of A in bytes used to store the matrix.
  template <class T, class Prop, class Storage, class Allocator>
  inline size_t Matrix_Triangular<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + size_t(GetDataSize())*sizeof(T);
    taille += sizeof(pointer)*this->m_;
    return taille;
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
  Matrix_Triangular<T, Prop, Storage, Allocator>::value_type
  Matrix_Triangular<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_Triangular");
#endif
    
    T zero;
    SetComplexZero(zero);
    
    if (Storage::UpLo())
      {
	if (i > j)
	  return zero;
	else
	  return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
      }
    else
      {
	if (i < j)
	  return zero;
	else
	  return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
      }
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
    \note It returns the value stored i, memory. Even if the value is in the
    part of the matrix that should be filled with zeros, the value could be
    non-zero.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Triangular<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Triangular<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBoundsTriang(i, j, this->m_, this->n_,
		      Storage::UpLo(), "Matrix_Triangular");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
    \note It returns the value stored i, memory. Even if the value is in the
    part of the matrix that should be filled with zeros, the value could be
    non-zero.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Triangular<T, Prop, Storage, Allocator>::reference
  Matrix_Triangular<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBoundsTriang(i, j, this->m_, this->n_,
		      Storage::UpLo(), "Matrix_Triangular");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Returns the element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Triangular<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Triangular<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    return this->Val(i, j);
  }


  //! Returns the element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Triangular<T, Prop, Storage, Allocator>::reference
  Matrix_Triangular<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
    return this->Val(i, j);
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Triangular<T, Prop, Storage, Allocator>::reference
  Matrix_Triangular<T, Prop, Storage, Allocator>::operator[] (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->GetDataSize(), "Matrix_Triangular");
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
  inline typename Matrix_Triangular<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Triangular<T, Prop, Storage, Allocator>::operator[] (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->GetDataSize(), "Matrix_Triangular");
#endif

    return this->data_[i];
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Triangular<T, Prop, Storage, Allocator>&
  Matrix_Triangular<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_Triangular<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& x)
  {
    this->Val(i, j) = x;
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>
  ::Copy(const Matrix_Triangular<T, Prop, Storage, Allocator>& A)
  {
    this->Reallocate(A.GetM(), A.GetN());

    Allocator::memorycpy(this->data_, A.GetData(), this->GetDataSize());
  }


#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Triangular<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_Triangular<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return false;
  }
#endif


  /////////////////////////
  // MATRIX<COLUPTRIANG> //
  /////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColUpTriang, Allocator>::Matrix():
    Matrix_Triangular<T, Prop, ColUpTriang, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j full column-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColUpTriang, Allocator>::Matrix(int i, int j):
    Matrix_Triangular<T, Prop, ColUpTriang, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills the matrix with a given value.
  /*!
    On exit, the matrix is filled with 'x'.
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColUpTriang, Allocator>&
  Matrix<T, Prop, ColUpTriang, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, ColUpTriang, Allocator>&
  Matrix<T, Prop, ColUpTriang, Allocator>::operator=
  (const Matrix<T, Prop, ColUpTriang, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColUpTriang, Allocator>&
  Matrix<T, Prop, ColUpTriang, Allocator>::operator*= (const T0& x)
  {
    for (long i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }



  /////////////////////////
  // MATRIX<COLLOTRIANG> //
  /////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColLoTriang, Allocator>::Matrix():
    Matrix_Triangular<T, Prop, ColLoTriang, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j full column-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColLoTriang, Allocator>::Matrix(int i, int j):
    Matrix_Triangular<T, Prop, ColLoTriang, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills the matrix with a given value.
  /*!
    On exit, the matrix is filled with 'x'.
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColLoTriang, Allocator>&
  Matrix<T, Prop, ColLoTriang, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, ColLoTriang, Allocator>&
  Matrix<T, Prop, ColLoTriang, Allocator>::operator=
  (const Matrix<T, Prop, ColLoTriang, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColLoTriang, Allocator>&
  Matrix<T, Prop, ColLoTriang, Allocator>::operator*= (const T0& x)
  {
    for (long i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }



  /////////////////////////
  // MATRIX<ROWUPTRIANG> //
  /////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowUpTriang, Allocator>::Matrix():
    Matrix_Triangular<T, Prop, RowUpTriang, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j full row-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowUpTriang, Allocator>::Matrix(int i, int j):
    Matrix_Triangular<T, Prop, RowUpTriang, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills the matrix with a given value.
  /*!
    On exit, the matrix is filled with 'x'.
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowUpTriang, Allocator>&
  Matrix<T, Prop, RowUpTriang, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, RowUpTriang, Allocator>&
  Matrix<T, Prop, RowUpTriang, Allocator>::operator=
  (const Matrix<T, Prop, RowUpTriang, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowUpTriang, Allocator>&
  Matrix<T, Prop, RowUpTriang, Allocator>::operator*= (const T0& x)
  {
    for (long i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }



  /////////////////////////
  // MATRIX<ROWLOTRIANG> //
  /////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowLoTriang, Allocator>::Matrix():
    Matrix_Triangular<T, Prop, RowLoTriang, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j full row-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowLoTriang, Allocator>::Matrix(int i, int j):
    Matrix_Triangular<T, Prop, RowLoTriang, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills the matrix with a given value.
  /*!
    On exit, the matrix is filled with 'x'.
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowLoTriang, Allocator>&
  Matrix<T, Prop, RowLoTriang, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, RowLoTriang, Allocator>&
  Matrix<T, Prop, RowLoTriang, Allocator>::operator=
  (const Matrix<T, Prop, RowLoTriang, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowLoTriang, Allocator>&
  Matrix<T, Prop, RowLoTriang, Allocator>::operator*= (const T0& x)
  {
    for (long i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }

  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_TRIANGULAR_INLINE_CXX
#endif
