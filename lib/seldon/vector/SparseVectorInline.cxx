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


#ifndef SELDON_FILE_SPARSE_VECTOR_INLINE_CXX

#include "SparseVector.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>::Vector():
    Vector<T, VectFull, Allocator>()
  {
    index_ = NULL;
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param i length of the vector.
  */
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>::Vector(size_t i):
    Vector<T, VectFull, Allocator>(i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->index_ = AllocatorInt::allocate(i, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->index_ = NULL;
	this->data_ = NULL;
      }

    if (this->index_ == NULL)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }

    if (this->data_ == NULL && i != 0)
      throw NoMemory("Vector<VectSparse>::Vector(int)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(i * sizeof(T)) + " bytes ("
		     + to_str(i) + " elements).");
#endif

  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param V vector to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>::
  Vector(const Vector<T, VectSparse, Allocator>& V) :
    Vector<T, VectFull, Allocator>()
  {
    this->index_ = NULL;
    Copy(V);
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, previous non-zero entries may be
    lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectSparse, Allocator>::Reallocate(size_t i)
  {
    ReallocateVector(i);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>::~Vector()
  {
    Clear();
  }
  

  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    \param i index.
    \return The value of the non-zero element #\a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::reference
  Vector<T, VectSparse, Allocator>::Value(int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "Vector<VectSparse>");
#endif

    return this->data_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The value of the non-zero element #\a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::const_reference
  Vector<T, VectSparse, Allocator>::Value(int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "Vector<VectSparse>");
#endif

    return this->data_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The index of the non-zero element #\a i.
  */
  template <class T, class Allocator>
  inline int& Vector<T, VectSparse, Allocator>::Index(int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "Vector<VectSparse>");
#endif

    return this->index_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The row number of the non-zero element #\a i.
  */
  template <class T, class Allocator>
  inline int Vector<T, VectSparse, Allocator>::Index(int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "Vector<VectSparse>");
#endif

    return this->index_[i];
  }


  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: \a X is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>& Vector<T, VectSparse, Allocator>
  ::operator= (const Vector<T, VectSparse, Allocator>& X)
  {
    this->Copy(X);

    return *this;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  /*! \brief Returns a pointer to the array containing the indices of the
    non-zero entries. */
  /*!
    \return A pointer to the array of the indices of the non-zero entries.
  */
  template <class T, class Allocator>
  inline int* Vector<T, VectSparse, Allocator>::GetIndex() const
  {
    return this->index_;
  }


  //! Returns the memory used by the object in bytes.
  /*!
    In this method, the type T is assumed to be "static"
    such that sizeof(T) provides the correct size
  */
  template <class T, class Allocator>
  inline size_t Vector<T, VectSparse, Allocator>::GetMemorySize() const
  {
    return sizeof(*this) + (sizeof(T) + sizeof(int))*size_t(this->m_);
  }
  
} // namespace Seldon.

#define SELDON_FILE_SPARSE_VECTOR_INLINE_CXX
#endif
