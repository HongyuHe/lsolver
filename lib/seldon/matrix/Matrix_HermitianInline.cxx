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


#ifndef SELDON_FILE_MATRIX_HERMITIAN_INLINE_CXX

#include "Matrix_Hermitian.hxx"

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
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>::Matrix_Hermitian():
    Matrix_Base<T, Allocator>()
  {
    me_ = NULL;
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>::~Matrix_Hermitian()
  {
    Clear();
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
  inline long Matrix_Hermitian<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return long(this->m_) * long(this->n_);
  }


  //! Returns size of A in bytes used to store the matrix.
  template <class T, class Prop, class Storage, class Allocator>
  inline size_t Matrix_Hermitian<T, Prop, Storage, Allocator>::GetMemorySize() const
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
  Matrix_Hermitian<T, Prop, Storage, Allocator>::value_type
  Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_Hermitian");
#endif

    if (i > j)
      return conj(me_[Storage::GetSecond(i, j)][Storage::GetFirst(i, j)]);
    else
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
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBoundsSym(i, j, this->m_, this->n_, "Matrix_Hermitian");
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
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBoundsSym(i, j, this->m_, this->n_, "Matrix_Hermitian");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }

  
  //! Returns access to an element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBoundsSym(i, j, this->m_, this->n_, "Matrix_Hermitian");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Returns the element (i, j).
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBoundsSym(i, j, this->m_, this->n_, "Matrix_Hermitian");
#endif
    
    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }
  
  
  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::operator[] (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->GetDataSize(), "Matrix_Hermitian");
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
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::operator[] (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->GetDataSize(), "Matrix_Hermitian");
#endif

    return this->data_[i];
  }

  
  //! Sets an element of the matrix
  /*!
    \param i row index
    \param j column index
    \param x sets a(i, j) = x
   */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& x)
  {
    if (i > j)
      this->Val(j, i) = conj(x);
    else
      this->Val(i, j) = x;
  }
  

  //! Duplicates a matrix (assignement operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>&
  Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::operator=(const Matrix_Hermitian<T, Prop, Storage, Allocator>& A)
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
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Copy(const Matrix_Hermitian<T, Prop, Storage, Allocator>& A)
  {
    this->Reallocate(A.GetM(), A.GetN());

    Allocator::memorycpy(this->data_, A.GetData(), this->GetDataSize());
  }


#ifdef SELDON_WITH_VIRTUAL
  //! adds val to A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      this->Val(i, j) += val;
  }


  //! Adds several values on row i
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::AddInteractionRow(int i, int n, const Vector<int>& col,
		      const Vector<T>& val, bool sorted)
  {
    for (int k = 0; k < n; k++)
      if (i <= col(k))
	this->Val(i, col(k)) += val(k);
  }


  //! Adds several values on column i
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::AddInteractionColumn(int i, int n, const Vector<int>& row,
			 const Vector<T>& val, bool sorted)
  {
    for (int k = 0; k < n; k++)
      if (row(k) <= i)
	this->Val(row(k), i) += val(k);
  }


  //! sets to zero the row i
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::ClearRow(int i)
  {
    for (int k = 0; k < this->n_; k++)
      SetComplexZero(this->Val(i, k));
  }


  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Tcplx>& x, const Vector<Tcplx>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_Hermitian<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return false;
  }
#endif

  
  ////////////////////
  // MATRIX<COLHERM> //
  ////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColHerm, Allocator>::Matrix():
    Matrix_Hermitian<T, Prop, ColHerm, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j full column-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColHerm, Allocator>::Matrix(int i, int j):
    Matrix_Hermitian<T, Prop, ColHerm, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills a matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator= (const Matrix<T, Prop,
                                                  ColHerm,
                                                  Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
    \warning The imaginary part of x should be null to keep an Hermitian
    matrix.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator*= (const T0& x)
  {
    long size = this->GetDataSize();
    for (long i = 0; i < size;i++)
      this->data_[i] *= x;

    return *this;
  }



  ////////////////////
  // MATRIX<ROWHERM> //
  ////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowHerm, Allocator>::Matrix():
    Matrix_Hermitian<T, Prop, RowHerm, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i by j full row-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowHerm, Allocator>::Matrix(int i, int j):
    Matrix_Hermitian<T, Prop, RowHerm, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills a matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator= (const Matrix<T, Prop,
                                                       RowHerm,
                                                       Allocator>& A)
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
  inline Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator*= (const T0& x)
  {
    long size = this->GetDataSize();
    for (long i = 0; i < size; i++)
      this->data_[i] *= x;

    return *this;
  }

} // namespace Seldon.

#define SELDON_FILE_MATRIX_HERMITIAN_INLINE_CXX
#endif
