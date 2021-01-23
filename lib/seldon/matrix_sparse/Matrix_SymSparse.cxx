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


#ifndef SELDON_FILE_MATRIX_SYMSPARSE_CXX

#include "Matrix_SymSparse.hxx"

namespace Seldon
{

  
  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to i so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Matrix_SymSparse(int i, int j,
		   Vector<T, Storage0, Allocator0>& values,
		   Vector<long, Storage1, Allocator1>& ptr,
		   Vector<int, Storage2, Allocator2>& ind):
    Matrix_Base<T, Allocator>(i, j)
  {
    nz_ = values.GetLength();

#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.

    if (ind.GetLength() != nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(nz_) + " values but "
		       + to_str(ind.GetLength()) + " row or column indices.");
      }

    if (ptr.GetLength() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices contains ")
		       + to_str(ptr.GetLength()-1) + string(" row or column ")
		       + string("start  indices (plus the number of non-zero")
		       + " entries) but there are " + to_str(i)
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (static_cast<long int>(2 * nz_ - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(values.GetLength())
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

    this->ptr_ = ptr.GetData();
    this->ind_ = ind.GetData();
    this->data_ = values.GetData();

    ptr.Nullify();
    ind.Nullify();
    values.Nullify();
  }


  //! Copy constructor
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Matrix_SymSparse(const Matrix_SymSparse<T, Prop, Storage, Allocator>& A)
    : Matrix_Base<T, Allocator>()
  {
    this->m_ = 0;
    this->n_ = 0;
    this->nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
    this->Copy(A);
  }


  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit,
    the matrix is empty (0x0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Clear()
  {
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (ptr_ != NULL)
	  {
	    AllocatorLong::deallocate(ptr_, this->m_+1);
	    ptr_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (ind_ != NULL)
	  {
	    AllocatorInt::deallocate(ind_, this->nz_);
	    ind_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    Allocator::deallocate(this->data_, nz_);
	    this->data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->nz_ = 0;
	this->data_ = NULL;
      }
#endif

    this->m_ = 0;
    this->n_ = 0;
    this->nz_ = 0;
  }


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by
    'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j,
	  Vector<T, Storage0, Allocator0>& values,
	  Vector<long, Storage1, Allocator1>& ptr,
	  Vector<int, Storage2, Allocator2>& ind)
  {
    this->Clear();
    this->m_ = i;
    this->n_ = i;
    this->nz_ = values.GetLength();

#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.

    if (ind.GetM() != nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Reallocate(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(nz_) + " values but "
		       + to_str(ind.GetLength()) + " row or column indices.");
      }

    if (ptr.GetM()-1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Reallocate(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices contains ")
		       + to_str(ptr.GetLength()-1) + string(" row or column")
		       + string(" start indices (plus the number of")
		       + string(" non-zero entries) but there are ")
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (static_cast<long int>(2 * nz_ - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Reallocate(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(values.GetLength())
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

    this->ptr_ = ptr.GetData();
    this->ind_ = ind.GetData();
    this->data_ = values.GetData();

    ptr.Nullify();
    ind.Nullify();
    values.Nullify();
  }


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by arrays
    'values' (values), 'ptr' (pointers) and 'ind' (indices).
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero entries that are stored.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning On exit, arrays 'values', 'ptr' and 'ind' are managed by the
    matrix.
    For example, it means that the destructor will released those arrays;
    therefore, the user mustn't release those arrays.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::SetData(int i, int j, long nz,
	    typename Matrix_SymSparse<T, Prop, Storage, Allocator>
	    ::pointer values, long* ptr, int* ind)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = i;

    this->nz_ = nz;

    this->data_ = values;
    ind_ = ind;
    ptr_ = ptr;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Nullify()
  {
    this->data_ = NULL;
    this->m_ = 0;
    this->n_ = 0;
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
  }


  //! Initialization of an empty sparse matrix with i rows and j columns
  /*!
    \param i number of rows
    \param j number of columns
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Reallocate(int i, int j)
  {
    // clearing previous entries
    Clear();

    this->m_ = i;
    this->n_ = i;

    // we try to allocate ptr_
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ptr_ = reinterpret_cast<long*>( AllocatorLong::
				       allocate(i+1, this) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_SymSparse::Reallocate(int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i+1) )
		     + " bytes to store " + to_str(i+1)
		     + " row or column start indices, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

    // then filing ptr_ with 0
    for (int k = 0; k <= i; k++)
      ptr_[k] = 0;
  }


  //! Initialization of a sparse matrix with i rows and j columns
  /*!
    \param i number of rows
    \param j number of columns
    \param nz number of non-zero entries
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j, long nz)
  {
    // clearing previous entries
    Clear();

    this->nz_ = nz;
    this->m_ = i;
    this->n_ = i;

#ifdef SELDON_CHECK_DIMENSIONS
    if (static_cast<long int>(2 * nz_ - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		       string("There are more values (") + to_str(nz)
		       + string(" values) than elements in the upper")
		       + " part of the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ptr_ = reinterpret_cast<long*>( AllocatorLong::
				       allocate(i + 1, this) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i+1) ) + " bytes to store "
		     + to_str(i+1) + " row or column start indices, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ind_ = reinterpret_cast<int*>( AllocatorInt::
				       allocate(nz_, this) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	AllocatorLong::deallocate(ptr_, i+1);
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	AllocatorLong::deallocate(ptr_, i+1);
	ptr_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz)
		     + " row or column indices, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = Allocator::allocate(nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(ptr_, i+1);
	ptr_ = NULL;
	AllocatorInt::deallocate(ind_, nz);
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(ptr_, i+1);
	ptr_ = NULL;
	AllocatorInt::deallocate(ind_, nz);
	ind_ = NULL;
      }
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz) + " values, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

    // then filing ptr_ with 0
    for (int k = 0; k <= i; k++)
      ptr_[k] = 0;
  }


  //! Changing the number of rows and columns
  /*!
    \param i number of rows
    \param j number of columns
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Resize(int i, int j)
  {
    if (i < this->m_)
      Resize(i, i, ptr_[i]);
    else
      Resize(i, i, nz_);
  }


  //! Changing the number of rows, columns and non-zero entries
  /*!
    \param i number of rows
    \param j number of columns
    Previous entries are kept during the operation
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Resize(int i, int j, long nz)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (static_cast<long int>(2 * nz - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		       string("There are more values (") + to_str(nz)
		       + string(" values) than elements in the upper")
		       + " part of the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

    if (nz != nz_)
      {
        // trying to resize ind_ and data_
#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif

            ind_ = reinterpret_cast<int*>( AllocatorInt::reallocate(ind_, nz) );
	    
#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            this->m_ = 0;
            this->n_ = 0;
            nz_ = 0;
            AllocatorLong::deallocate(ptr_, i+1);
            ptr_ = NULL;
            ind_ = NULL;
            this->data_ = NULL;
          }
        if (ind_ == NULL)
          {
            this->m_ = 0;
            this->n_ = 0;
            nz_ = 0;
            AllocatorLong::deallocate(ptr_, i+1);
            ptr_ = NULL;
            this->data_ = NULL;
          }
        if (ind_ == NULL && i != 0 && j != 0)
          throw NoMemory("Matrix_SymSparse::Resize(int, int, int)",
                         string("Unable to allocate ") + to_str(sizeof(int) * nz)
                         + " bytes to store " + to_str(nz)
                         + " row or column indices, for a "
                         + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

        Vector<T, VectFull, Allocator> val;
        val.SetData(nz_, this->data_);
        val.Resize(nz);

        this->data_ = val.GetData();
        nz_ = nz;
        val.Nullify();
      }


    if (this->m_ != i)
      {
#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif
            // trying to resize ptr_
            ptr_ = reinterpret_cast<long*>( AllocatorLong::reallocate(ptr_, i+1) );

#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            this->m_ = 0;
            this->n_ = 0;
            nz_ = 0;
            ptr_ = NULL;
            ind_ = NULL;
            this->data_ = NULL;
          }
        if (ptr_ == NULL)
          {
            this->m_ = 0;
            this->n_ = 0;
            nz_ = 0;
            ind_ = NULL;
            this->data_ = NULL;
          }
        if (ptr_ == NULL && i != 0 && j != 0)
          throw NoMemory("Matrix_SymSparse::Resize(int, int)",
                         string("Unable to allocate ")
                         + to_str(sizeof(int) * (i+1) )
                         + " bytes to store " + to_str(i+1)
                         + " row or column start indices, for a "
                         + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

        // then filing last values of ptr_ with nz_
        for (int k = this->m_; k <= i; k++)
          ptr_[k] = this->nz_;
      }

    this->m_ = i;
    this->n_ = i;
  }


  //! Copies a matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Copy(const Matrix_SymSparse<T, Prop, Storage, Allocator>& A)
  {
    this->Clear();
    long nz = A.GetNonZeros();
    int i = A.GetM();
    int j = A.GetN();
    this->nz_ = nz;
    this->m_ = i;
    this->n_ = j;
    if ((i == 0)||(j == 0))
      {
	this->m_ = 0;
	this->n_ = 0;
	this->nz_ = 0;
	return;
      }

#ifdef SELDON_CHECK_DIMENSIONS
    if (static_cast<long int>(2 * nz_ - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		       string("There are more values (") + to_str(nz)
		       + string(" values) than elements in the upper")
		       + " part of the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ptr_ = reinterpret_cast<long*>( AllocatorLong::allocate(i + 1, this) );
	AllocatorLong::memorycpy(this->ptr_, A.ptr_, (i + 1));
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i+1) ) + " bytes to store "
		     + to_str(i+1) + " row or column start indices, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ind_ = reinterpret_cast<int*>( AllocatorInt::allocate(nz_, this) );
	AllocatorInt::memorycpy(this->ind_, A.ind_, nz_);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	AllocatorLong::deallocate(ptr_, i+1);
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	AllocatorLong::deallocate(ptr_, i+1);
	ptr_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz)
		     + " row or column indices, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = Allocator::allocate(nz_, this);
	Allocator::memorycpy(this->data_, A.data_, nz_);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(ptr_, i+1);
	ptr_ = NULL;
	AllocatorInt::deallocate(ind_, nz);
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	AllocatorLong::deallocate(ptr_, i+1);
	ptr_ = NULL;
	AllocatorInt::deallocate(ind_, nz);
	ind_ = NULL;
      }
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz) + " values, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

  }
  

  //! returns size of matrix in bytes
  template<class T, class Prop, class Storage, class Allocator>
  size_t Matrix_SymSparse<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + this->GetPtrSize()*sizeof(long);
    int coef = sizeof(T) + sizeof(int); // for each non-zero entry
    taille += coef*size_t(this->nz_);
    return taille;
  }


  //! Fills vector ptr with integers instead of longs
  template<class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::FillPtrInt(Vector<int>& Ptr) const
  {
    Ptr.Reallocate(this->m_+1);
    for (int i = 0; i <= this->m_; i++)
      Ptr(i) = ptr_[i];
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
  const typename
  Matrix_SymSparse<T, Prop, Storage, Allocator>::value_type
  Matrix_SymSparse<T, Prop, Storage, Allocator>::operator() (int i,
							     int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k, l;
    long a, b;
    T zero;
    SetComplexZero(zero);

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      return zero;

    l = Storage::GetSecond(i, j);

    for (k = a; (k<b-1) && (ind_[k]<l); k++);

    if (ind_[k] == l)
      return this->data_[k];
    else
      return zero;
  }


  //! Access method.
  /*! Returns the value of element (\a i, \a j) if it can be returned as a
    reference.
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the matrix).
  */
  template <class T, class Prop, class Storage, class Allocator>
  typename Matrix_SymSparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_SymSparse<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymSparse::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymSparse::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k, l;
    long a, b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      throw WrongArgument("Matrix_SymSparse::Val(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");

    l = Storage::GetSecond(i, j);

    for (k = a; (k<b-1) && (ind_[k]<l); k++);

    if (ind_[k] == l)
      return this->data_[k];
    else
      throw WrongArgument("Matrix_SymSparse::Val(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");
  }


  //! Access method.
  /*! Returns the value of element (\a i, \a j) if it can be returned as a
    reference.
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the matrix).
  */
  template <class T, class Prop, class Storage, class Allocator>
  const typename Matrix_SymSparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_SymSparse<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymSparse::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymSparse::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k, l;
    long a, b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      throw WrongArgument("Matrix_SymSparse::Val(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");

    l = Storage::GetSecond(i, j);

    for (k = a; (k<b-1) && (ind_[k]<l); k++);

    if (ind_[k] == l)
      return this->data_[k];
    else
      throw WrongArgument("Matrix_SymSparse::Val(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");
  }

  
  //! Access method.
  /*! Returns reference to element (\a i, \a j) 
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    If the element does not belong to sparsity pattern of the matrix,
    the matrix is resized.
  */
  template <class T, class Prop, class Storage, class Allocator>
  typename Matrix_SymSparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_SymSparse<T, Prop, Storage, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymSparse::Get(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymSparse::Get(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    long k, l;
    long a, b;
    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a < b)
      {
        l = Storage::GetSecond(i, j);
        
        for (k = a; (k < b) && (ind_[k] < l); k++);

        if ( (k < b) && (ind_[k] == l))
          return this->data_[k];
      }
    else
      k = a;
    
    // adding a non-zero entry
    Resize(this->m_, this->n_, nz_+1);
    
    for (int m = Storage::GetFirst(i, j)+1; m <= this->m_; m++)
      ptr_[m]++;
    
    for (long m = nz_-1; m >= k+1; m--)
      {
        ind_[m] = ind_[m-1];
        this->data_[m] = this->data_[m-1];
      }
    
    ind_[k] = Storage::GetSecond(i, j);
    
    // value of new non-zero entry is set to 0
    SetComplexZero(this->data_[k]);
    
    return this->data_[k];
  }
  
  
  /************************
   * CONVENIENT FUNCTIONS *
   ************************/

  
  //! Resets all non-zero entries to 0-value.
  /*! The sparsity pattern remains unchanged. */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Zero()
  {
    Allocator::memoryset(this->data_, char(0),
			 this->nz_ * sizeof(value_type));
  }

  
  //! Sets the matrix to identity.
  /*! This method fills the diagonal of the matrix with ones.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::SetIdentity()
  {
    int m = this->m_;
    int nz = this->m_;
    T one;
    SetComplexOne(one);
    
    if (nz == 0)
      return;
    
    Clear();

    Vector<T, VectFull, Allocator> values(nz);
    Vector<long> ptr(m + 1);
    Vector<int> ind(nz);
    
    values.Fill(one);
    ind.Fill();
    ptr.Fill();
    
    SetData(m, m, values, ptr, ind);
  }


  //! Fills the non-zero entries with 0, 1, 2, ...
  /*! On exit, the non-zero entries are 0, 1, 2, 3, ... The order of the
    numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Fill()
  {
    for (long i = 0; i < this->GetDataSize(); i++)
      SetComplexReal(i, this->data_[i]);
  }


  //! Fills the non-zero entries with a given value.
  /*!
    \param x the value to set the non-zero entries to.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    T x_;
    SetComplexReal(x, x_);
    for (long i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = x_;
  }


  //! Fills the non-zero entries randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::FillRand()
  {
#ifndef SELDON_WITHOUT_REINIT_RANDOM
    srand(time(NULL));
#endif
    for (long i = 0; i < this->GetDataSize(); i++)
      SetComplexReal(rand(), this->data_[i]);
  }

  
  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Print() const
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
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_SymSparse::Write(string FileName)",
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
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_SymSparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif
    
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<long*>(&this->nz_)),
		     sizeof(long));
    
    FileStream.write(reinterpret_cast<char*>(this->ptr_),
		     sizeof(long)*(this->m_+1));
    FileStream.write(reinterpret_cast<char*>(this->ind_),
		     sizeof(int)*this->nz_);
    FileStream.write(reinterpret_cast<char*>(this->data_),
		     sizeof(T)*this->nz_);
  }
  
  
  //! Writes the matrix in a file.
  /*! Stores the matrix in a file in ascii format. The entries are written in
    coordinate format (row index, column index, value).  Row and column
    indexes start at 1.
    \param FileName output file name.
    \param cplx if true the real part and imaginary part are written
          in two separate columns, otherwise the complex values
          are written (a,b)
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::WriteText(string FileName, bool cplx) const
  {
    ofstream FileStream; FileStream.precision(14);
    FileStream.open(FileName.c_str());

    // changing precision
    FileStream.precision(cout.precision());
    
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_SymSparse::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream, cplx);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*! Stores the matrix in a file in ascii format. The entries are written in
    coordinate format (row index, column index, value).  Row and column
    indexes start at 1.
    \param FileStream output file name.
    \param cplx if true the real part and imaginary part are given
    in two separate columns, otherwise the complex values
    are written (a,b)
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream, bool cplx) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_SymSparse::WriteText(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    // Conversion to coordinate format (1-index convention).
    const Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this);
    
    T zero; int index = 1;
    WriteCoordinateMatrix(leaf_class, FileStream, zero, index, cplx);
  }

  
  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_SymSparse::Read(string FileName)",
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
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_SymSparse::Read(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    int m, n; long nz;
    FileStream.read(reinterpret_cast<char*>(&m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&n), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&nz), sizeof(long));

    Reallocate(m, m, nz);

    FileStream.read(reinterpret_cast<char*>(ptr_),
                    sizeof(long)*(m+1));
    FileStream.read(reinterpret_cast<char*>(ind_), sizeof(int)*nz);
    FileStream.read(reinterpret_cast<char*>(this->data_), sizeof(T)*nz);

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_SymSparse::Read(istream& FileStream)",
                    string("Input operation failed.")
		    + string(" The input file may have been removed")
		    + " or may not contain enough data.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads the matrix from a file in text format.
    \param FileName input file name.
    \param cplx if true the real part and imaginary part are given
          in two separate columns, otherwise the complex values
          are written (a,b)
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::ReadText(string FileName, bool cplx)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_SymSparse::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream, cplx);
    
    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix from a stream in text format.
    \param FileStream input stream.
    \param cplx if true the real part and imaginary part are given
          in two separate columns, otherwise the complex values
          are written (a,b)
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::
  ReadText(istream& FileStream, bool cplx)
  {
    Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);
    
    T zero; int index = 1;
    ReadCoordinateMatrix(leaf_class, FileStream, zero, index, -1, cplx);
  }
  
  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMSPARSE_CXX
#endif
