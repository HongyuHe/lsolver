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


#ifndef SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_INLINE_CXX

#include "Matrix_SymComplexSparse.hxx"

namespace Seldon
{
  
  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::~Matrix_SymComplexSparse()
  {
    this->Clear();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the cumulated number of non-zero entries of both the real and
    the imaginary part.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetNonZeros() const
  {
    return real_nz_ + imag_nz_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the cumulated number of non-zero entries of both the real and
    the imaginary part.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetDataSize() const
  {
    return real_nz_ + imag_nz_;
  }


  //! Returns (row or column) start indices for the real part.
  /*!
    Returns the array ('ptr_') of start indices for the real part.
    \return The array of start indices for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealPtr() const
  {
    return real_ptr_;
  }


  //! Returns (row or column) start indices for the imaginary part.
  /*!
    Returns the array ('ptr_') of start indices for the imaginary part.
    \return The array of start indices for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagPtr() const
  {
    return imag_ptr_;
  }


  //! Returns (row or column) indices of non-zero entries for the real part.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries for the real part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealInd() const
  {
    return real_ind_;
  }


  //! Returns (row or column) indices of non-zero entries for
  //! the imaginary part.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries for the imaginary part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagInd() const
  {
    return imag_ind_;
  }


  //! Returns the length of the array of start indices for the real part.
  /*!
    \return The length of the array of start indices for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealPtrSize() const
  {
    return (this->m_ + 1);
  }


  //! Returns the length of the array of start indices for the imaginary part.
  /*!
    \return The length of the array of start indices for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagPtrSize() const
  {
    return (this->m_ + 1);
  }


  //! Returns the length of the array of (column or row) indices for
  //! the real part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored) for the real part. This array
    defines non-zero entries indices if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices for
    the real part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealIndSize() const
  {
    return real_nz_;
  }


  //! Returns the length of the array of (column or row) indices
  //! for the imaginary part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored) for the imaginary part. This array
    defines non-zero entries indices if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices
    for the imaginary part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagIndSize() const
  {
    return imag_nz_;
  }

  
  //! Returns the length of the array of (column or row) indices for
  //! the real part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored) for the real part. This array
    defines non-zero entries indices if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices for
    the real part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealDataSize() const
  {
    return real_nz_;
  }


  //! Returns the length of the array of (column or row) indices
  //! for the imaginary part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored) for the imaginary part. This array
    defines non-zero entries indices if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices
    for the imaginary part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline long Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagDataSize() const
  {
    return imag_nz_;
  }

  
  //! Returns the array of values of the real part.
  /*!
    \return The array 'real_data_' of values of the real part..
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type*
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::GetRealData() const
  {
    return real_data_;
  }


  //! Returns the array of values of the imaginary part.
  /*!
    \return The array 'imag_data_' of values of the imaginary part..
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type*
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::GetImagData() const
  {
    return imag_data_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access method
  /*! Returns the real part of element (\a i, \a j)
    if it can be returned as a reference. If the non-zero entry
    does not exit, it is created
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::GetReal(int i, int j) const
  {
    return ValReal(i, j);
  }
  
  
  //! Access method
  /*! Returns the imaginary part of element (\a i, \a j)
    if it can be returned as a reference. If the non-zero entry
    does not exit, it is created
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::GetImag(int i, int j) const
  {
    return ValImag(i, j);
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
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const entry_type& val)
  {
    if (i <= j)
      {
        if (real(val) != value_type(0))
          GetReal(i, j) += real(val);
        
        if (imag(val) != value_type(0))
          GetImag(i, j) += imag(val);
      }
  }

  
  //! Adds values to several non-zero entries on a given row
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::AddInteractionRow(int i, int nb, const Vector<int>& col,
		      const Vector<entry_type>& val, bool sorted)
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
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const entry_type& val)
  {
    GetReal(i, j) = real(val);
    GetImag(i, j) = imag(val);
  }
  

  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymComplexSparse<T, Prop, Storage, Allocator>&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_SymComplexSparse<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Tcplx>& x, const Vector<Tcplx>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    cout << "Not defined for real numbers" << endl;
    abort();
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_SymComplexSparse<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return true;
  }
#endif

  
  /////////////////////////////////
  // MATRIX<COLSYMCOMPLEXSPARSE> //
  /////////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymComplexSparse, Allocator>::Matrix() :
    Matrix_SymComplexSparse<T, Prop, ColSymComplexSparse, Allocator>()
  {
  }


  //! Builds a i by j matrix.
  /*!
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymComplexSparse, Allocator>
  ::Matrix(int i, int j):
    Matrix_SymComplexSparse<T, Prop, ColSymComplexSparse, Allocator>(i, j,
                                                                     0, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with real_nz and imag_nz non-zero (stored)
    elements for the real part and the imaginary part respectively.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements that are stored
    for the real part.
    \param imag_nz number of non-zero elements that are stored
    for the imaginary part.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator> inline
  Matrix<T, Prop, ColSymComplexSparse, Allocator>::Matrix(int i, int j,
							  long real_nz,
							  long imag_nz):
    Matrix_SymComplexSparse<T, Prop,
    ColSymComplexSparse, Allocator>(i, j,
				    real_nz, imag_nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix<T, Prop, ColSymComplexSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<value_type, Storage0, Allocator0>& real_values,
	 Vector<long, Storage1, Allocator1>& real_ptr,
	 Vector<int, Storage2, Allocator2>& real_ind,
	 Vector<value_type, Storage0, Allocator0>& imag_values,
	 Vector<long, Storage1, Allocator1>& imag_ptr,
	 Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_SymComplexSparse<T, Prop,
    ColSymComplexSparse, Allocator>(i, j,
				    real_values,
				    real_ptr,
				    real_ind,
				    imag_values,
				    imag_ptr,
				    imag_ind)
  {
  }



  /////////////////////////////////
  // MATRIX<ROWSYMCOMPLEXSPARSE> //
  /////////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymComplexSparse, Allocator>::Matrix() :
    Matrix_SymComplexSparse<T, Prop, RowSymComplexSparse, Allocator>()
  {
  }


  //! Builds a i by j matrix.
  /*!
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymComplexSparse, Allocator>
  ::Matrix(int i, int j):
    Matrix_SymComplexSparse<T, Prop, RowSymComplexSparse, Allocator>(i, j,
                                                                     0, 0)
  {
  }


  /*! Builds a i by j matrix with real_nz and imag_nz non-zero (stored)
    elements for the real part and the imaginary part respectively.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements that are stored
    for the real part.
    \param imag_nz number of non-zero elements that are store
    for the imaginary part.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymComplexSparse, Allocator>
  ::Matrix(int i, int j, long real_nz, long imag_nz):
    Matrix_SymComplexSparse<T, Prop, RowSymComplexSparse, Allocator>(i, j,
								     real_nz,
								     imag_nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix<T, Prop, RowSymComplexSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<value_type, Storage0, Allocator0>& real_values,
	 Vector<long, Storage1, Allocator1>& real_ptr,
	 Vector<int, Storage2, Allocator2>& real_ind,
	 Vector<value_type, Storage0, Allocator0>& imag_values,
	 Vector<long, Storage1, Allocator1>& imag_ptr,
	 Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_SymComplexSparse<T, Prop,
    RowSymComplexSparse, Allocator>(i, j,
				    real_values,
				    real_ptr,
				    real_ind,
				    imag_values,
				    imag_ptr,
				    imag_ind)
  {
  }
  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_INLINE_CXX
#endif
