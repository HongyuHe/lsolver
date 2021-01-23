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


#ifndef SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_CXX

#include "Matrix_SymComplexSparse.hxx"

namespace Seldon
{

  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymComplexSparse(): Matrix_Base<T, Allocator>()
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymComplexSparse(int i, int j): Matrix_Base<T, Allocator>()
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;

    Reallocate(i, i);
  }


  //! Constructor.
  /*! Builds a sparse matrix of size i by j , with real_nz
    non-zero (stored) elements in the real part of the matrix and imag_nz
    non-zero elements in the imaginary part of the matrix.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements that are stored
    for the real part.
    \param imag_nz number of non-zero elements that are stored
    for the imaginary part.
    \note Matrix values are not initialized. Indices of non-zero entries
    are not initialized either.
  */
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymComplexSparse(int i, int j, long real_nz, long imag_nz):
    Matrix_Base<T, Allocator>()
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
    
    Reallocate(i, i, real_nz, imag_nz);
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
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_SymComplexSparse(int i, int j,
			  Vector<value_type, Storage0, Allocator0>& real_values,
			  Vector<long, Storage1, Allocator1>& real_ptr,
			  Vector<int, Storage2, Allocator2>& real_ind,
			  Vector<value_type, Storage0, Allocator0>& imag_values,
			  Vector<long, Storage1, Allocator1>& imag_ptr,
			  Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_Base<T, Allocator>(i, j)
  {
    real_nz_ = real_values.GetM();
    imag_nz_ = imag_values.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.

    if (real_ind.GetM() != real_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(real_nz_)
		       + " values (real part) but "
		       + to_str(real_ind.GetM())
		       + " row or column indices.");
      }

    if (imag_ind.GetM() != imag_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(imag_nz_)
		       + " values (imaginary part) but "
		       + to_str(imag_ind.GetM())
		       + " row or column indices.");
      }

    if (real_ptr.GetM() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (real part)")
		       + " contains " + to_str(real_ptr.GetM()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (imag_ptr.GetM() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (imaginary part)")
		       + " contains " + to_str(imag_ptr.GetM()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if ( (static_cast<long int>(2 * real_nz_ - 2)
	  / static_cast<long int>(i + 1)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(2 * imag_nz_ - 2)
	  / static_cast<long int>(i + 1)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(real_values.GetM())
		       + " values for the real part and "
		       + to_str(real_values.GetM()) + " values for"
		       + string(" the imaginary part) than elements in the")
		       + " matrix (" + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

    this->real_ptr_ = real_ptr.GetData();
    this->imag_ptr_ = imag_ptr.GetData();
    this->real_ind_ = real_ind.GetData();
    this->imag_ind_ = imag_ind.GetData();
    this->real_data_ = real_values.GetData();
    this->imag_data_ = imag_values.GetData();

    real_ptr.Nullify();
    imag_ptr.Nullify();
    real_ind.Nullify();
    imag_ind.Nullify();
    real_values.Nullify();
    imag_values.Nullify();
  }


  //! Copy constructor
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymComplexSparse(const Matrix_SymComplexSparse<T, Prop,
			    Storage, Allocator>& A)
    : Matrix_Base<T, Allocator>()
  {
    this->m_ = 0;
    this->n_ = 0;
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;

    this->Copy(A);
  }

  
  /*********************
   * MEMORY MANAGEMENT *
   *********************/


//! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix
    is empty (0x0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::Clear()
  {
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (real_ptr_ != NULL)
	  {
	    AllocatorLong::deallocate(real_ptr_, this->m_+1);
	    real_ptr_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	real_ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (imag_ptr_ != NULL)
	  {
	    AllocatorLong::deallocate(imag_ptr_, this->m_+1);
	    imag_ptr_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	imag_ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (real_ind_ != NULL)
	  {
	    AllocatorInt::deallocate(real_ind_, this->real_nz_);
	    real_ind_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	real_ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (imag_ind_ != NULL)
	  {
	    AllocatorInt::deallocate(imag_ind_, this->imag_nz_);
	    imag_ind_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	imag_ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->real_data_ != NULL)
	  {
	    Allocator::deallocate(this->real_data_, real_nz_);
	    this->real_data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->real_nz_ = 0;
	this->real_data_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->imag_data_ != NULL)
	  {
	    Allocator::deallocate(this->imag_data_, imag_nz_);
	    this->imag_data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->imag_nz_ = 0;
	this->imag_data_ = NULL;
      }
#endif

    this->m_ = 0;
    this->n_ = 0;
    this->real_nz_ = 0;
    this->imag_nz_ = 0;
  }


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by
    'real_values' (values of the real part), 'real_ptr'
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
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j,
	  Vector<value_type, Storage0, Allocator0>& real_values,
	  Vector<long, Storage1, Allocator1>& real_ptr,
	  Vector<int, Storage2, Allocator2>& real_ind,
	  Vector<value_type, Storage0, Allocator0>& imag_values,
	  Vector<long, Storage1, Allocator1>& imag_ptr,
	  Vector<int, Storage2, Allocator2>& imag_ind)
  {
    this->Clear();
    this->m_ = i;
    this->n_ = i;
    real_nz_ = real_values.GetM();
    imag_nz_ = imag_values.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.

    if (real_ind.GetM() != real_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(real_nz_)
		       + " values (real part) but "
		       + to_str(real_ind.GetM())
		       + " row or column indices.");
      }

    if (imag_ind.GetM() != imag_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(imag_nz_)
		       + " values (imaginary part) but "
		       + to_str(imag_ind.GetM())
		       + " row or column indices.");
      }

    if (real_ptr.GetM() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (real part)")
		       + " contains " + to_str(real_ptr.GetM() - 1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (imag_ptr.GetM() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (imaginary part)")
		       + " contains " + to_str(imag_ptr.GetM()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if ( (static_cast<long int>(2 * real_nz_ - 2) / static_cast<long int>(i)
	  >= static_cast<long int>(i + 1)) ||
	 (static_cast<long int>(2 * imag_nz_ - 2) / static_cast<long int>(i)
	  >= static_cast<long int>(i + 1)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(real_values.GetM())
		       + " values for the real part and "
		       + to_str(real_values.GetM()) + string(" values")
		       + string(" for the imaginary part) than elements in")
		       + " the matrix (" + to_str(i)
		       + " by " + to_str(i) + ").");
      }
#endif

    this->real_ptr_ = real_ptr.GetData();
    this->imag_ptr_ = imag_ptr.GetData();
    this->real_ind_ = real_ind.GetData();
    this->imag_ind_ = imag_ind.GetData();
    this->real_data_ = real_values.GetData();
    this->imag_data_ = imag_values.GetData();

    real_ptr.Nullify();
    imag_ptr.Nullify();
    real_ind.Nullify();
    imag_ind.Nullify();
    real_values.Nullify();
    imag_values.Nullify();
  }


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by arrays
    'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part).
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero entries (real part).
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_nz number of non-zero entries (imaginary part).
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning On exit, arrays 'real_values', 'real_ptr', 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are managed by the matrix.
    For example, it means that the destructor will release those arrays;
    therefore, the user mustn't release those arrays.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j, long real_nz,
	  typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
	  ::pointer real_values,
	  long* real_ptr, int* real_ind, long imag_nz,
	  typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
	  ::pointer imag_values,
	  long* imag_ptr, int* imag_ind)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = i;

    this->real_nz_ = real_nz;
    this->imag_nz_ = imag_nz;

    real_data_ = real_values;
    imag_data_ = imag_values;
    real_ind_ = real_ind;
    imag_ind_ = imag_ind;
    real_ptr_ = real_ptr;
    imag_ptr_ = imag_ptr;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::Nullify()
  {
    this->data_ = NULL;
    this->m_ = 0;
    this->n_ = 0;
    real_nz_ = 0;
    real_ptr_ = NULL;
    real_ind_ = NULL;
    imag_nz_ = 0;
    imag_ptr_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
  }

  
  //! Initialization of an empty matrix i x j
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j)
  {
    // previous entries are removed
    Clear();
    
    this->m_ = i;
    this->n_ = i;
    
    // we try to allocate real_ptr_ and imag_ptr_
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	real_ptr_
	  = reinterpret_cast<long*>( AllocatorLong::allocate(i+1) );

	imag_ptr_
	  = reinterpret_cast<long*>( AllocatorLong::allocate(i+1) );
        
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
        real_nz_ = 0;
        real_ptr_ = NULL;
        real_ind_ = NULL;
        imag_nz_ = 0;
        imag_ptr_ = NULL;
        imag_ind_ = NULL;
        real_data_ = NULL;
        imag_data_ = NULL;
      }
    if ((real_ptr_ == NULL) || (imag_ptr_ == NULL))
      {
	this->m_ = 0;
	this->n_ = 0;
        real_nz_ = 0;
        real_ptr_ = NULL;
        real_ind_ = NULL;
        imag_nz_ = 0;
        imag_ptr_ = NULL;
        imag_ind_ = NULL;
        real_data_ = NULL;
        imag_data_ = NULL;
      }
    if (((real_ptr_ == NULL) || (imag_ptr_ == NULL)) && i != 0 && j != 0)
      throw NoMemory("Matrix_SymComplexSparse::Reallocate(int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i+1) )
		     + " bytes to store " + to_str(i+1)
		     + " row or column start indices, for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif
    
    // then filing real_ptr_ with 0
    for (int k = 0; k <= i; k++)
      {
        real_ptr_[k] = 0;
        imag_ptr_[k] = 0;
      }
  }
  
  
  //! Initialization of a matrix i x j with n non-zero entries
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j, long real_nz, long imag_nz)
  {
    // previous entries are removed
    Clear();
    
    this->m_ = i;
    this->n_ = j;
    this->real_nz_ = real_nz;
    this->imag_nz_ = imag_nz;
    
#ifdef SELDON_CHECK_DIMENSIONS
    if ( (static_cast<long int>(2 * real_nz_ - 2) / static_cast<long int>(i+1)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(2 * imag_nz_ - 2) / static_cast<long int>(i+1)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + "Matrix_SymComplexSparse(int, int, int, int)",
		       string("There are more values to be stored (")
		       + to_str(real_nz) + " values for the real part and "
		       + to_str(imag_nz) + string(" values for the imaginary")
		       + " part) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	real_ptr_ =
	  reinterpret_cast<long*>( AllocatorLong::allocate(i + 1) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	imag_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i + 1)) + " bytes to store "
		     + to_str(i + 1) + string(" row or column")
		     + " start indices (for the real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	imag_ptr_ =
	  reinterpret_cast<long*>( AllocatorLong::allocate(i + 1) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	real_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i + 1) ) + " bytes to store "
		     + to_str(i + 1) + string(" row or column")
		     + " start indices (for the imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	real_ind_ = reinterpret_cast<int*>( AllocatorInt::allocate(real_nz_) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz) + " bytes to store "
		     + to_str(real_nz)
		     + " row or column indices (real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	imag_ind_ =
	  reinterpret_cast<int*>( AllocatorInt::allocate(imag_nz_) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(imag_ind_, imag_nz);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(imag_ind_, imag_nz);
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ind_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz) + " bytes to store "
		     + to_str(imag_nz)
		     + " row or column indices (imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->real_data_ = Allocator::allocate(real_nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(real_ind_, real_nz);
	AllocatorInt::deallocate(imag_ind_, imag_nz);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(real_ind_, real_nz);
	AllocatorInt::deallocate(imag_ind_, imag_nz);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	imag_data_ = NULL;
      }
    if (real_data_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz) + " bytes to store "
		     + to_str(real_nz) + " values (real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->imag_data_ = Allocator::allocate(imag_nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(real_ind_, real_nz);
	AllocatorInt::deallocate(imag_ind_, imag_nz);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	Allocator::deallocate(this->real_data_, real_nz_);
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(real_ind_, real_nz);
	AllocatorInt::deallocate(imag_ind_, imag_nz);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	Allocator::deallocate(this->real_data_, real_nz_);
	real_data_ = NULL;
      }
    if (imag_data_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz) + " bytes to store "
		     + to_str(imag_nz) + " values (imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

    // then filing real_ptr_ with 0
    for (int k = 0; k <= Storage::GetFirst(i, j); k++)
      {
        real_ptr_[k] = 0;
        imag_ptr_[k] = 0;
      }
  }

  
  //! Changing the number of rows and columns
  /*!
    \param i number of rows
    \param j number of columns    
   */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Resize(int i, int j)
  {
    if (i < this->m_)
      Resize(i, i, real_ptr_[i], imag_ptr_[i]);
    else
      Resize(i, i, real_nz_, imag_nz_);
  }
  
  
  //! Changing the number of rows, columns and non-zero entries
  /*!
    \param i number of rows
    \param j number of columns
    Previous entries are kept during the operation
   */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Resize(int i, int j, long real_nz, long imag_nz)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (real_nz < 0 || imag_nz < 0)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + "Resize(int, int, int, int)",
		       "Invalid number of non-zero elements: "
                       + to_str(real_nz) + " in the real part and "
                       + to_str(imag_nz) + " in the imaginary part.");
      }
    if ((real_nz > 0
         && (j == 0
             || static_cast<long int>(real_nz-1) / static_cast<long int>(j)
             >= static_cast<long int>(i)))
        ||
        (imag_nz > 0
         && (j == 0
             || static_cast<long int>(imag_nz-1) / static_cast<long int>(j)
             >= static_cast<long int>(i))))
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + "Resize(int, int, int, int)",
		       string("There are more values (") + to_str(real_nz)
		       + " values for the real part and " + to_str(imag_nz)
		       + string(" values for the imaginary part) than")
		       + " elements in the matrix (" + to_str(i) + " by "
		       + to_str(j) + ").");
      }
#endif

    if (this->m_ != i)
      {
#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif
            
            real_ptr_
              = reinterpret_cast<long*>( AllocatorLong::
					reallocate(real_ptr_, (i+1)) );
	    
#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            this->m_ = 0;
            this->n_ = 0;
            real_nz_ = 0;
            imag_nz_ = 0;
            real_ptr_ = NULL;
            imag_ptr_ = NULL;
            real_ind_ = NULL;
            imag_ind_ = NULL;
            this->real_data_ = NULL;
            this->imag_data_ = NULL;
          }
        if (real_ptr_ == NULL)
          {
            this->m_ = 0;
            this->n_ = 0;
            real_nz_ = 0;
            imag_nz_ = 0;
            imag_ptr_ = 0;
            real_ind_ = NULL;
            imag_ind_ = NULL;
            this->real_data_ = NULL;
            this->imag_data_ = NULL;
          }
        if (real_ptr_ == NULL && i != 0 && j != 0)
          throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Resize(int, int, int, int)",
                         string("Unable to allocate ")
                         + to_str(sizeof(int) * (i+1))
                         + " bytes to store " + to_str(i+1)
                         + string(" row or column start indices (for the real")
                         + " part), for a "
                         + to_str(i) + " by " + to_str(j) + " matrix.");
#endif
        
#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif
            
            imag_ptr_ = 
	      reinterpret_cast<long*>( AllocatorLong::
				      reallocate(imag_ptr_, (i+1)) );
            
#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            this->m_ = 0;
            this->n_ = 0;
            real_nz_ = 0;
            imag_nz_ = 0;
            AllocatorLong::deallocate(real_ptr_, i+1);
            real_ptr_ = NULL;
            imag_ptr_ = NULL;
            real_ind_ = NULL;
            imag_ind_ = NULL;
            this->real_data_ = NULL;
            this->imag_data_ = NULL;
          }
        if (imag_ptr_ == NULL)
          {
            this->m_ = 0;
            this->n_ = 0;
            real_nz_ = 0;
            imag_nz_ = 0;
            AllocatorLong::deallocate(real_ptr_, i+1);
            real_ptr_ = 0;
            real_ind_ = NULL;
            imag_ind_ = NULL;
            this->real_data_ = NULL;
            this->imag_data_ = NULL;
          }
        if (imag_ptr_ == NULL && i != 0 && j != 0)
          throw NoMemory(string("Matrix_SymComplexSparse::")
                         + "Resize(int, int, int, int)",
                         string("Unable to allocate ")
                         + to_str(sizeof(int) * (i+1))
                         + " bytes to store " + to_str(i+1)
                         + string(" row or column start indices (for the")
                         + string(" imaginary part), for a ")
                         + to_str(i) + " by " + to_str(j) + " matrix.");
#endif
      }
    
    if (real_nz != real_nz_)
      {
#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif
            
            real_ind_ =
              reinterpret_cast<int*>( AllocatorInt::
				      reallocate(real_ind_, real_nz) );

#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            this->m_ = 0;
            this->n_ = 0;
            real_nz_ = 0;
            imag_nz_ = 0;
            AllocatorLong::deallocate(real_ptr_, i+1);
            AllocatorLong::deallocate(imag_ptr_, i+1);
            real_ptr_ = NULL;
            imag_ptr_ = NULL;
            real_ind_ = NULL;
            imag_ind_ = NULL;
            this->real_data_ = NULL;
            this->imag_data_ = NULL;
          }
        if (real_ind_ == NULL)
          {
            this->m_ = 0;
            this->n_ = 0;
            real_nz_ = 0;
            imag_nz_ = 0;
            AllocatorLong::deallocate(real_ptr_, i+1);
            AllocatorLong::deallocate(imag_ptr_, i+1);
            real_ptr_ = NULL;
            imag_ptr_ = NULL;
            this->real_data_ = NULL;
            this->imag_data_ = NULL;
          }
        if (real_ind_ == NULL && i != 0 && j != 0)
          throw NoMemory(string("Matrix_SymComplexSparse::")
                         + "Resize(int, int, int, int)",
                         string("Unable to allocate ")
                         + to_str(sizeof(int) * real_nz)
                         + " bytes to store " + to_str(real_nz)
                         + " row or column indices (real part), for a "
                         + to_str(i) + " by " + to_str(j) + " matrix.");
#endif
      }
    
    if (imag_nz != imag_nz_)
      {
#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif

            imag_ind_ =
	      reinterpret_cast<int*>( AllocatorInt::
				      reallocate(imag_ind_, imag_nz) );
            
#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            this->m_ = 0;
            this->n_ = 0;
            real_nz_ = 0;
            imag_nz_ = 0;
	    AllocatorLong::deallocate(real_ptr_, i+1);
            AllocatorLong::deallocate(imag_ptr_, i+1);
            real_ptr_ = NULL;
            imag_ptr_ = NULL;
            AllocatorInt::deallocate(imag_ind_, imag_nz);
            real_ind_ = NULL;
            imag_ind_ = NULL;
            this->real_data_ = NULL;
            this->imag_data_ = NULL;
          }
        if (imag_ind_ == NULL)
          {
            this->m_ = 0;
            this->n_ = 0;
            real_nz_ = 0;
            imag_nz_ = 0;
            AllocatorLong::deallocate(real_ptr_, i+1);
            AllocatorLong::deallocate(imag_ptr_, i+1);
            real_ptr_ = NULL;
            imag_ptr_ = NULL;
            AllocatorInt::deallocate(imag_ind_, imag_nz);
            imag_ind_ = NULL;
            this->real_data_ = NULL;
            this->imag_data_ = NULL;
          }
        if (imag_ind_ == NULL && i != 0 && j != 0)
          throw NoMemory(string("Matrix_SymComplexSparse::")
                         + "Resize(int, int, int, int)",
                         string("Unable to allocate ")
                         + to_str(sizeof(int) * imag_nz)
                         + " bytes to store " + to_str(imag_nz)
                         + " row or column indices (imaginary part), for a "
                         + to_str(i) + " by " + to_str(j) + " matrix.");
#endif
      }
    
    if (real_nz != real_nz_)
      {
        Vector<value_type, VectFull, Allocator> val;
        val.SetData(real_nz_, real_data_);
        val.Resize(real_nz);
        
        real_data_ = val.GetData();
        val.Nullify();
      }

    if (imag_nz != imag_nz_)
      {
        Vector<value_type, VectFull, Allocator> val;
        val.SetData(imag_nz_, imag_data_);
        val.Resize(imag_nz);
        
        imag_data_ = val.GetData();
        val.Nullify();
      }

    // then filing last values of ptr_ with nz_
    for (int k = this->m_; k <= i; k++)
      {
        real_ptr_[k] = real_nz_;
        imag_ptr_[k] = imag_nz_;
      }
    
    this->m_ = i;
    this->n_ = i;
    real_nz_ = real_nz;
    imag_nz_ = imag_nz;
  }
  
  
  //! Copies a matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  Copy(const Matrix_SymComplexSparse<T, Prop, Storage, Allocator>& A)
  {
    this->Clear();
    int i = A.m_;
    int j = A.n_;
    real_nz_ = A.real_nz_;
    imag_nz_ = A.imag_nz_;
    this->m_ = i;
    this->n_ = j;
    if ((i == 0)||(j == 0))
      {
	this->m_ = 0;
	this->n_ = 0;
	this->real_nz_ = 0;
	this->imag_nz_ = 0;
	return;
      }

#ifdef SELDON_CHECK_DIMENSIONS
    if ( (static_cast<long int>(2 * real_nz_ - 2) / static_cast<long int>(i+1)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(2 * imag_nz_ - 2) / static_cast<long int>(i+1)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + "Matrix_SymComplexSparse(int, int, int, int)",
		       string("There are more values to be stored (")
		       + to_str(real_nz_) + " values for the real part and "
		       + to_str(imag_nz_) + string(" values for the imaginary")
		       + " part) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	real_ptr_ = 
	  reinterpret_cast<long*>( AllocatorLong::allocate(i + 1) );
	
	AllocatorLong::
	  memorycpy(this->real_ptr_, A.real_ptr_, (i + 1));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	imag_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i + 1)) + " bytes to store "
		     + to_str(i + 1) + string(" row or column")
		     + " start indices (for the real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	imag_ptr_ =
	  reinterpret_cast<long*>( AllocatorLong::allocate(i + 1) );
	
	AllocatorLong::memorycpy(this->imag_ptr_, A.imag_ptr_, (i + 1));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	real_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i + 1) ) + " bytes to store "
		     + to_str(i + 1) + string(" row or column")
		     + " start indices (for the imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	real_ind_ =
	  reinterpret_cast<int*>( AllocatorInt::allocate(real_nz_) );
	
	AllocatorInt::memorycpy(this->real_ind_, A.real_ind_, real_nz_);
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz_) + " bytes to store "
		     + to_str(real_nz_)
		     + " row or column indices (real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	imag_ind_ =
	  reinterpret_cast<int*>( AllocatorInt::allocate(imag_nz_) );
	
	AllocatorInt::memorycpy(this->imag_ind_, A.imag_ind_, imag_nz_);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(imag_ind_, A.imag_nz_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(imag_ind_, A.imag_nz_);
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ind_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz_) + " bytes to store "
		     + to_str(imag_nz_)
		     + " row or column indices (imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->real_data_ = Allocator::allocate(real_nz_, this);
	Allocator::memorycpy(this->real_data_, A.real_data_, real_nz_);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(real_ind_, A.real_nz_);
	AllocatorInt::deallocate(imag_ind_, A.imag_nz_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(real_ind_, A.real_nz_);
	AllocatorInt::deallocate(imag_ind_, A.imag_nz_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	imag_data_ = NULL;
      }
    if (real_data_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz_) + " bytes to store "
		     + to_str(real_nz_) + " values (real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->imag_data_ = Allocator::allocate(imag_nz_, this);
	Allocator::memorycpy(this->imag_data_, A.imag_data_, imag_nz_);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(real_ind_, A.real_nz_);
	AllocatorInt::deallocate(imag_ind_, A.imag_nz_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	Allocator::deallocate(this->real_data_, real_nz_);
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(real_ptr_, i+1);
	AllocatorLong::deallocate(imag_ptr_, i+1);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	AllocatorInt::deallocate(real_ind_, A.real_nz_);
	AllocatorInt::deallocate(imag_ind_, A.imag_nz_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	Allocator::deallocate(this->real_data_, real_nz_);
	real_data_ = NULL;
      }
    if (imag_data_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz_) + " bytes to store "
		     + to_str(imag_nz_) + " values (imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

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
  const
  complex<typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
          ::value_type>
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymComplexSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymComplexSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long real_k, imag_k; int l;
    long real_a, real_b;
    long imag_a, imag_b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    real_a = real_ptr_[Storage::GetFirst(i, j)];
    real_b = real_ptr_[Storage::GetFirst(i, j) + 1];

    imag_a = imag_ptr_[Storage::GetFirst(i, j)];
    imag_b = imag_ptr_[Storage::GetFirst(i, j) + 1];

    if (real_a != real_b)
      {
	l = Storage::GetSecond(i, j);
	for (real_k = real_a;
	     (real_k < real_b - 1) && (real_ind_[real_k] < l);
	     real_k++);
	if (imag_a != imag_b)
	  {
	    for (imag_k = imag_a;
		 (imag_k < imag_b - 1) && (imag_ind_[imag_k] < l);
		 imag_k++);
	    if (real_ind_[real_k] == l)
	      {
		if (imag_ind_[imag_k] == l)
		  return entry_type(real_data_[real_k], imag_data_[imag_k]);
		else
		  return entry_type(real_data_[real_k], value_type(0));
	      }
	    else
	      if (imag_ind_[imag_k] == l)
		return entry_type(value_type(0), imag_data_[imag_k]);
	      else
		return entry_type(value_type(0), value_type(0));
	  }
	else
	  {
	    if (real_ind_[real_k] == l)
	      return entry_type(real_data_[real_k], value_type(0));
	    else
	      return entry_type(value_type(0), value_type(0));
	  }
      }
    else
      {
	if (imag_a != imag_b)
	  {
	    l = Storage::GetSecond(i, j);
	    for (imag_k = imag_a;
		 (imag_k < imag_b - 1) && (imag_ind_[imag_k] < l);
		 imag_k++);
	    if (imag_ind_[imag_k] == l)
	      return entry_type(value_type(0), imag_data_[imag_k]);
	    else
	      return entry_type(value_type(0), value_type(0));
	  }
	else
	  return entry_type(value_type(0), value_type(0));
      }
  }

  
  //! Access method
  /*! Returns the real value of element (\a i, \a j)
    if it can be returned as a reference.
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the matrix).
  */
  template <class T, class Prop, class Storage, class Allocator>
  typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::ValReal(int i, int j)
  {
   
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymComplexSparse::ValReal(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymComplexSparse::ValReal(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k; int l;
    long a, b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = real_ptr_[Storage::GetFirst(i, j)];
    b = real_ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      throw WrongArgument("Matrix_SymComplexSparse::ValReal(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");

    l = Storage::GetSecond(i, j);

    for (k = a; (k < b-1) && (real_ind_[k] < l); k++);

    if (real_ind_[k] == l)
      return this->real_data_[k];
    else
      throw WrongArgument("Matrix_SymComplexSparse::ValReal(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry."); 
  }

  
  //! Access method
  /*! Returns the real value of element (\a i, \a j)
    if it can be returned as a reference.
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the matrix).
  */
  template <class T, class Prop, class Storage, class Allocator>
  const typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::ValReal(int i, int j) const
  {
   
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymComplexSparse::ValReal(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymComplexSparse::ValReal(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k; int l;
    long a, b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = real_ptr_[Storage::GetFirst(i, j)];
    b = real_ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      throw WrongArgument("Matrix_SymComplexSparse::ValReal(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");

    l = Storage::GetSecond(i, j);

    for (k = a; (k < b-1) && (real_ind_[k] < l); k++);

    if (real_ind_[k] == l)
      return this->real_data_[k];
    else
      throw WrongArgument("Matrix_SymComplexSparse::ValReal(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry."); 
  }

  
  //! Access method
  /*! Returns the imaginary part of element (\a i, \a j)
    if it can be returned as a reference.
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the matrix).
  */
  template <class T, class Prop, class Storage, class Allocator>
  typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::ValImag(int i, int j)
  {
   
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymComplexSparse::ValImag(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymComplexSparse::ValImag(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k; int l;
    long a, b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }
    
    a = imag_ptr_[Storage::GetFirst(i, j)];
    b = imag_ptr_[Storage::GetFirst(i, j) + 1];
    
    if (a == b)
      throw WrongArgument("Matrix_SymComplexSparse::ValImag(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");

    l = Storage::GetSecond(i, j);

    for (k = a; (k < b-1) && (imag_ind_[k] < l); k++);

    if (imag_ind_[k] == l)
      return this->imag_data_[k];
    else
      throw WrongArgument("Matrix_SymComplexSparse::ValImag(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry."); 
  }

  
  //! Access method
  /*! Returns the imaginary part of element (\a i, \a j)
    if it can be returned as a reference.
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the matrix).
  */
  template <class T, class Prop, class Storage, class Allocator>
  const typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::ValImag(int i, int j) const
  {
   
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymComplexSparse::ValImag(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymComplexSparse::ValImag(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k; int l;
    long a, b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = imag_ptr_[Storage::GetFirst(i, j)];
    b = imag_ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      throw WrongArgument("Matrix_SymComplexSparse::ValImag(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");

    l = Storage::GetSecond(i, j);

    for (k = a; (k < b-1) && (imag_ind_[k] < l); k++);

    if (imag_ind_[k] == l)
      return this->imag_data_[k];
    else
      throw WrongArgument("Matrix_SymComplexSparse::ValImag(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry."); 
  }

  
  //! Access method
  /*! Returns the real part of element (\a i, \a j)
    if it can be returned as a reference. If the non-zero entry
    does not exit, it is created
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::GetReal(int i, int j)
  {
   
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymComplexSparse::GetReal(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymComplexSparse::GetReal(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k; int l;
    long a, b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = real_ptr_[Storage::GetFirst(i, j)];
    b = real_ptr_[Storage::GetFirst(i, j) + 1];

    if (a < b)
      {
        l = Storage::GetSecond(i, j);
        
        for (k = a; (k < b) && (real_ind_[k] < l); k++);

        if ( (k < b) && (real_ind_[k] == l))
          return this->real_data_[k];
      }
    else
      k = a;
    
    // adding a non-zero entry
    Resize(this->m_, this->n_, real_nz_+1, imag_nz_);
    
    for (int m = Storage::GetFirst(i, j)+1;
         m <= Storage::GetFirst(this->m_, this->n_); m++)
      real_ptr_[m]++;
    
    for (long m = real_nz_-1; m >= k+1; m--)
      {
        real_ind_[m] = real_ind_[m-1];
        this->real_data_[m] = this->real_data_[m-1];
      }
    
    real_ind_[k] = Storage::GetSecond(i, j);

    // value of new non-zero entry is set to 0
    SetComplexZero(this->real_data_[k]);
    
    return this->real_data_[k];
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
  typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::GetImag(int i, int j)
  {
   
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymComplexSparse::GetImag(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymComplexSparse::GetImag(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k; int l;
    long a, b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = imag_ptr_[Storage::GetFirst(i, j)];
    b = imag_ptr_[Storage::GetFirst(i, j) + 1];

    if (a < b)
      {
        l = Storage::GetSecond(i, j);
        
        for (k = a; (k < b) && (imag_ind_[k] < l); k++);

        if ( (k < b) && (imag_ind_[k] == l))
          return this->imag_data_[k];
      }
    else
      k = a;
    
    // adding a non-zero entry
    Resize(this->m_, this->n_, real_nz_, imag_nz_+1);
    
    for (int m = Storage::GetFirst(i, j)+1;
         m <= Storage::GetFirst(this->m_, this->n_); m++)
      imag_ptr_[m]++;
    
    for (long m = imag_nz_-1; m >= k+1; m--)
      {
        imag_ind_[m] = imag_ind_[m-1];
        this->imag_data_[m] = this->imag_data_[m-1];
      }
    
    imag_ind_[k] = Storage::GetSecond(i, j);

    // value of new non-zero entry is set to 0
    SetComplexZero(this->imag_data_[k]);
    
    return this->imag_data_[k];
  }
  
  
  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! returns size of A in bytes used to store the matrix
  template<class T, class Prop, class Storage, class Allocator>
  size_t Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + 2*this->GetRealPtrSize()*sizeof(long);
    int coef = sizeof(value_type) + sizeof(int); // for each non-zero entry
    taille += coef*size_t(this->real_nz_ + this->imag_nz_);
    return taille;
  }
  
  
  //! Resets all non-zero entries to 0-value.
  /*! The sparsity pattern remains unchanged. */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::Zero()
  {
    Allocator::memoryset(this->real_data_, char(0),
			 this->real_nz_ * sizeof(value_type));
    
    Allocator::memoryset(this->imag_data_, char(0),
			 this->imag_nz_ * sizeof(value_type));
  }


  //! Sets the matrix to identity.
  /*! This method fills the diagonal of the matrix with ones.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::SetIdentity()
  {
    int m = this->m_;
    int n = this->n_;
    int nz = min(m, n);

    if (nz == 0)
      return;

    Clear();

    Vector<value_type, VectFull, Allocator> real_values(nz), imag_values;
    Vector<long> real_ptr(m + 1);
    Vector<int> real_ind(nz);
    Vector<long> imag_ptr(real_ptr);
    Vector<int> imag_ind;
    
    real_values.Fill(value_type(1));
    real_ind.Fill();
    imag_ind.Zero();
    real_ptr.Fill();
    imag_ptr.Zero();
    
    SetData(m, n, real_values, real_ptr, real_ind,
            imag_values, imag_ptr, imag_ind);
  }


  //! Fills the non-zero entries with 0, 1, 2, ...
  /*! On exit, the non-zero entries are 0, 1, 2, 3, ... The order of the
    numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::Fill()
  {
    for (long i = 0; i < this->real_nz_; i++)
      this->real_data_[i] = i;

    for (long i = 0; i < this->imag_nz_; i++)
      this->imag_data_[i] = value_type(0);
  }


  //! Fills the non-zero entries with a given value.
  /*!
    \param x the value to set the non-zero entries to.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Fill(const entry_type& x)
  {
    for (long i = 0; i < this->real_nz_; i++)
      this->real_data_[i] = real(x);

    for (long i = 0; i < this->imag_nz_; i++)
      this->imag_data_[i] = imag(x);
  }


  //! Fills the non-zero entries randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::FillRand()
  {
#ifndef SELDON_WITHOUT_REINIT_RANDOM
    srand(time(NULL));
#endif
    for (long i = 0; i < this->real_nz_; i++)
      SetComplexReal(rand(), this->real_data_[i]);

    for (long i = 0; i < this->imag_nz_; i++)
      SetComplexReal(rand(), this->imag_data_[i]);
  }

  
  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }
  
  
  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in binary format.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_SymComplexSparse::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();
  }
  

  //! Writes the matrix to an output stream.
  /*!
    Stores the matrix in an output stream in binary format.
    \param FileStream output stream.
  */  
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_SymComplexSparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif
    
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<long*>(&this->real_nz_)),
		     sizeof(long));
    FileStream.write(reinterpret_cast<char*>(const_cast<long*>(&this->imag_nz_)),
		     sizeof(long));
    
    FileStream.write(reinterpret_cast<char*>(this->real_ptr_),
		     sizeof(long)*(Storage::GetFirst(this->m_, this->n_)+1));
    FileStream.write(reinterpret_cast<char*>(this->real_ind_),
		     sizeof(int)*this->real_nz_);
    FileStream.write(reinterpret_cast<char*>(this->real_data_),
		     sizeof(value_type)*this->real_nz_);

    FileStream.write(reinterpret_cast<char*>(this->imag_ptr_),
		     sizeof(long)*(Storage::GetFirst(this->m_, this->n_)+1));
    FileStream.write(reinterpret_cast<char*>(this->imag_ind_),
		     sizeof(int)*this->imag_nz_);
    FileStream.write(reinterpret_cast<char*>(this->imag_data_),
		     sizeof(value_type)*this->imag_nz_);
  }
  
  
  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in ascii format.
    The entries are written in coordinate format (row column value)
    1-index convention is used
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  WriteText(string FileName, bool cplx) const
  {
    ofstream FileStream; FileStream.precision(14);
    FileStream.open(FileName.c_str());

    // changing precision
    FileStream.precision(cout.precision());
    
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_SymComplexSparse::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream, cplx);

    FileStream.close();
  }

  
  //! Writes the matrix to an output stream.
  /*!
    Stores the matrix in a file in ascii format.
    The entries are written in coordinate format (row column value)
    1-index convention is used
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  WriteText(ostream& FileStream, bool cplx) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_SymComplexSparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    // conversion in coordinate format (1-index convention)
    const Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this);

    entry_type zero; int index = 1;
    WriteCoordinateMatrix(leaf_class, FileStream, zero, index, cplx);
  }
  
  
  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_SymComplexSparse::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    \param FileStream input stream
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_SymComplexSparse::Read(istream& FileStream)",
		    "Stream is not ready.");
#endif
    
    int m, n, real_nz, imag_nz;
    FileStream.read(reinterpret_cast<char*>(&m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&n), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&real_nz), sizeof(long));
    FileStream.read(reinterpret_cast<char*>(&imag_nz), sizeof(long));
    
    Reallocate(m, n, real_nz, imag_nz);

    FileStream.read(reinterpret_cast<char*>(real_ptr_),
                    sizeof(long)*(Storage::GetFirst(m, n)+1));
    FileStream.read(reinterpret_cast<char*>(real_ind_), sizeof(int)*real_nz);
    FileStream.read(reinterpret_cast<char*>(this->real_data_),
		    sizeof(value_type)*real_nz);

    FileStream.read(reinterpret_cast<char*>(imag_ptr_),
                    sizeof(long)*(Storage::GetFirst(m, n)+1));
    FileStream.read(reinterpret_cast<char*>(imag_ind_), sizeof(int)*imag_nz);
    FileStream.read(reinterpret_cast<char*>(this->imag_data_),
		    sizeof(value_type)*imag_nz);
    
#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_SymComplexSparse::Read(istream& FileStream)",
                    string("Input operation failed.")
		    + string(" The input file may have been removed")
		    + " or may not contain enough data.");
#endif

  }

  
  //! Reads the matrix from a file.
  /*!
    Reads the matrix from a file in text format.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::ReadText(string FileName, bool cplx)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_SymComplexSparse::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream, cplx);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix from a stream in text format.
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::ReadText(istream& FileStream, bool cplx)
  {
    Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);
    
    entry_type zero; int index = 1;
    ReadCoordinateMatrix(leaf_class, FileStream, zero, index, -1, cplx);
  }

  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_CXX
#endif
