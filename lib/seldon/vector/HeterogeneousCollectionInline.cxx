// Copyright (C) 2001-2009 INRIA
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


#ifndef SELDON_FILE_VECTOR_HETEROGENEOUSCOLLECTION_INLINE_CXX


#include "HeterogeneousCollection.hxx"


namespace Seldon
{


  ////////////////////////////////////
  // VECTOR HETEROGENEOUSCOLLECTION //
  ////////////////////////////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Default constructor.
  /*!
    Nothing is allocated. The number of vectors is set to zero.
  */
  template <class T, template <class U> class Allocator >
  inline Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::Vector(): Vector_Base<T, Allocator<T> >(), label_map_(),
	      label_vector_()
  {
    Nvector_ = 0;
  }


  //! Copy constructor.
  /*! Builds a copy of a vector collection.
    \param[in] V vector collection to be copied.
  */
  template <class T, template <class U> class Allocator >
  inline Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::Vector(const
	   Vector<FloatDouble, DenseSparseCollection, Allocator<T> >& V):
    Vector_Base<T, Allocator<T> >(V), label_map_(), label_vector_()
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
  template <class T, template <class U> class Allocator >
  inline int Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetM()
    const
  {
    return this->m_;
  }


  //! Returns the total number of elements.
  /*!
    \return The total length of the vector.
  */
  template <class T, template <class U> class Allocator >
  inline int Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetLength() const
  {
    return this->m_;
  }


  //! Returns the number of aggregated vectors.
  /*!
    \return The total number of aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline int Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetNvector() const
  {
    return Nvector_;
  }


  //! Returns the length vector of the underlying vectors.
  /*!
    \return The lengths of the underlying vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetVectorLength() const
  {
    return length_;
  }


  //! Returns the cumulative sum of the lengths of the underlying vectors.
  /*!
    \return The cumulative sum of the lengths of the underlying vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetLengthSum() const
  {
    return length_sum_;
  }


  //! Returns the collection indexes of the underlying vectors.
  /*!
    \return The collection index vector of the underlying vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetCollectionIndex() const
  {
    return collection_;
  }


  //! Returns the index of the underlying vectors in the inner collection.
  /*!
    \return The index vector of the underlying vectors in the inner
    collection.
  */
  template <class T, template <class U> class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetSubvectorIndex() const
  {
    return length_sum_;
  }


  //! Returns the list of float dense vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::float_dense_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetFloatDense()
  {
    return float_dense_c_;
  }


  //! Returns the list of float dense vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const
  typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::float_dense_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetFloatDense() const
  {
    return float_dense_c_;
  }


  //! Returns the list of float sparse vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline typename
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::float_sparse_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetFloatSparse()
  {
    return float_sparse_c_;
  }


  //! Returns the list of float sparse vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const
  typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::float_sparse_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetFloatSparse() const
  {
    return float_sparse_c_;
  }


  //! Returns the list of double dense vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::double_dense_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetDoubleDense()
  {
    return double_dense_c_;
  }


  //! Returns the list of double dense vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline
  const typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::double_dense_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetDoubleDense() const
  {
    return double_dense_c_;
  }


  //! Returns the list of double sparse vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::double_sparse_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetDoubleSparse()
  {
    return double_sparse_c_;
  }


  //! Returns the list of double sparse vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline
  const typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::double_sparse_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetDoubleSparse() const
  {
    return double_sparse_c_;
  }


  //! Duplicates a vector collection (assignment operator).
  /*!
    \param[in] X vector collection to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, template <class U> class Allocator >
  inline Vector<FloatDouble, DenseSparseCollection, Allocator<T> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::operator=
  (const Vector<FloatDouble, DenseSparseCollection, Allocator<T> >& X)
  {
    this->Copy(X);
    return *this;
  }
  
} // namespace Seldon.


#define SELDON_FILE_VECTOR_HETEROGENEOUSCOLLECTION_INLINE_CXX
#endif
