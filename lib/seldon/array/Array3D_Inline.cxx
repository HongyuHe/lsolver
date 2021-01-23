// Copyright (C) 2001-2009 Vivien Mallet
// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_ARRAY3D_INLINE_CXX

#include "Array3D.hxx"

namespace Seldon
{
  
  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Array3D<T, Allocator>::~Array3D()
  {
    Clear();
  }


  /*****************
   * BASIC METHODS *
   *****************/


  //! Returns the length in dimension #1.
  /*!
    \return The length in dimension #1.
  */
  template <class T, class Allocator>
  inline int Array3D<T, Allocator>::GetLength1() const
  {
    return length1_;
  }


  //! Returns the length in dimension #2.
  /*!
    \return The length in dimension #2.
  */
  template <class T, class Allocator>
  inline int Array3D<T, Allocator>::GetLength2() const
  {
    return length2_;
  }


  //! Returns the length in dimension #3.
  /*!
    \return The length in dimension #3.
  */
  template <class T, class Allocator>
  inline int Array3D<T, Allocator>::GetLength3() const
  {
    return length3_;
  }


  //! Returns the number of elements in the 3D array.
  /*!
    Returns the number of elements stored by the 3D array, i.e.
    the product of the lengths in the three dimensions.
    \return The number of elements in the 3D array.
  */
  template <class T, class Allocator>
  inline long Array3D<T, Allocator>::GetSize() const
  {
    return length1_ * length23_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory by
    the array, i.e. the product of lengths in the three
    dimensions.
    \return The number of elements stored in the array.
  */
  template <class T, class Allocator>
  inline long Array3D<T, Allocator>::GetDataSize() const
  {
    return length1_ * length23_;
  }


  //! Returns a pointer to the data array.
  /*!
    Returns a pointer to data, i.e. the data array 'data_' which stores the
    values.
    \return A pointer to the data array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::pointer Array3D<T, Allocator>
  ::GetData() const
  {
    return data_;
  }


  //! Returns a pointer to an element of data array.
  /*!
    Returns a pointer to an element of data array.
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return A pointer to the data array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::pointer Array3D<T, Allocator>
  ::GetDataPointer(int i, int j, int k) const
  {
    return data_ + i * length23_ + j * length3_ + k;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::reference
  Array3D<T, Allocator>::operator() (int i, int j, int k)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, k, length1_, length2_, length3_, "Array3D");
#endif

    return data_[i * length23_ + j * length3_ + k];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::const_reference
  Array3D<T, Allocator>::operator() (int i, int j, int k) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, k, length1_, length2_, length3_, "Array3D");
#endif

    return data_[i*length23_ + j*length3_ + k];
  }

  //! Duplicates a 3D array (assignment operator).
  /*!
    \param A 3D array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Array3D<T, Allocator>& Array3D<T, Allocator>::operator=
  (const Array3D<T, Allocator>& A)
  {
    Copy(A);

    return *this;
  }

  //! Duplicates a 3D array.
  /*!
    \param A 3D array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Array3D<T, Allocator>::Copy(const Array3D<T, Allocator>& A)
  {
    Reallocate(A.GetLength1(), A.GetLength2(), A.GetLength3());

    Allocator::memorycpy(data_, A.GetData(), GetDataSize());
  }
  
} // namespace Seldon.

#define SELDON_FILE_ARRAY3D_INLINE_CXX
#endif
