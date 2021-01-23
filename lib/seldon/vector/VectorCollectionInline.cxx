// Copyright (C) 2010, INRIA
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef SELDON_FILE_VECTOR_VECTORCOLLECTION_INLINE_CXX


#include "VectorCollection.hxx"


namespace Seldon
{


  //////////////////////
  // VECTORCOLLECTION //
  //////////////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Default constructor.
  /*!
    Nothing is allocated. The vector length is set to zero.
  */
  template <class T, class Allocator>
  inline Vector<T, Collection, Allocator>::Vector():
    Vector_Base<T, Allocator>(), label_map_(), label_vector_()
  {
    Nvector_ = 0;
  }


  //! Main constructor.
  /*! Builds a vector collection of a given size.
    \param[in] i length of the vector.
  */
  template <class T, class Allocator>
  inline Vector<T, Collection, Allocator>::Vector(int i):
    Vector_Base<T, Allocator>(0), length_(i), length_sum_(i), vector_(i),
    label_map_(), label_vector_()
  {
    Nvector_ = i;
    for (int k = 0; k < i; k++)
      length_(k) = length_sum_(k) = 0;
  }


  //! Copy constructor.
  /*! Builds a copy of a vector collection.
    \param[in] V vector collection to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, Collection, Allocator>::
  Vector(const Vector<T, Collection, Allocator>& V):
    Vector_Base<T, Allocator>(V), Nvector_(0)
  {
    Copy(V);
  }


  /*****************
   * BASIC METHODS *
   *****************/


  //! Returns the total number of elements.
  /*!
    \return The total length of the vector.
  */
  template <class T, class Allocator >
  inline int Vector<T, Collection, Allocator>::GetM() const
  {
    return this->m_;
  }


  //! Returns the total number of elements.
  /*!
    \return The total length of the vector.
  */
  template <class T, class Allocator >
  inline int Vector<T, Collection, Allocator>::GetLength() const
  {
    return this->m_;
  }


  //! Returns the number of aggregated vectors.
  /*!
    \return The total number of aggregated vectors.
  */
  template <class T, class Allocator >
  inline int Vector<T, Collection, Allocator>::GetNvector() const
  {
    return Nvector_;
  }


  //! Returns the length vector of the underlying vectors.
  /*!
    \return The lengths of the underlying vectors.
  */
  template <class T, class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<T, Collection, Allocator>::GetVectorLength() const
  {
    return length_;
  }


  //! Returns the cumulative sum of the lengths of the underlying vectors.
  /*!
    \return The cumulative sum of the lengths of the underlying vectors.
  */
  template <class T, class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<T, Collection, Allocator>::GetLengthSum() const
  {
    return length_sum_;
  }


  //! Returns the list of vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::collection_reference
  Vector<T, Collection, Allocator>::GetVector()
  {
    return vector_;
  }


  //! Returns the list of vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::const_collection_reference
  Vector<T, Collection, Allocator>::GetVector() const
  {
    return vector_;
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::vector_reference
  Vector<T, Collection, Allocator>::GetVector(int i)
  {
    return vector_(i);
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, class Allocator >
  inline typename
  Vector<T, Collection, Allocator>::const_vector_reference
  Vector<T, Collection, Allocator>::GetVector(int i) const
  {
    return vector_(i);
  }


  //! Multiplies a vector collection by a scalar.
  /*!
    \param[in] alpha scalar.
  */
  template <class T, class Allocator>
  template<class T0>
  inline Vector<T, Collection, Allocator>&
  Vector<T, Collection, Allocator>::operator*= (const T0& alpha)
  {
    for (int i = 0; i < this->Nvector_; i++)
      this->vector_(i) *= alpha;

    return *this;
  }
  
  
} // namespace Seldon.


#define SELDON_FILE_VECTOR_VECTORCOLLECTION_INLINE_CXX
#endif
