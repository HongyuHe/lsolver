// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu
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


#ifndef SELDON_FILE_HETEROGENEOUS_MATRIX_COLLECTION_INLINE_CXX

#include "HeterogeneousMatrixCollection.hxx"

namespace Seldon
{


  ///////////////////////////////////
  // HETEROGENEOUSMATRIXCOLLECTION //
  ///////////////////////////////////


  //! Copy constructor.
  /*!
    \param[in] A matrix collection to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::HeterogeneousMatrixCollection(const HeterogeneousMatrixCollection<Prop0,
				  Storage0, Prop1, Storage1, Allocator>& A):
    Matrix_Base<double, Allocator<double> >()
  {
    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::~HeterogeneousMatrixCollection()
  {
    this->Clear();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of rows.
  /*!
    \return the total number of rows. It is the sum of the number of rows in
    the underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetM() const
  {
    return this->m_;
  }


  //! Returns the number of rows.
  /*!
    \return the total number of rows. It is the sum of the number of rows
    in the underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetMmatrix() const
  {
    return Mmatrix_;
  }


  //! Returns the number of rows in an underlying matrix.
  /*!
    \param[in] i row index of the underlying matrix.
    \return The number of rows in the underlying matrices with row index \a i.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetM(int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::GetM()",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
		     + to_str(i) + ".");
#endif

    return Mlocal_(i);
  }


  //! Returns the number of columns.
  /*!
    \return the total number of columns. It is the sum of the number of
    columns in the underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetN() const
  {
    return this->n_;
  }


  //! Returns the number of columns.
  /*!
    \return the total number of columns. It is the sum of the number of
    columns in the underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetNmatrix() const
  {
    return Nmatrix_;
  }


  //! Returns the number of columns in an underlying matrix.
  /*!
    \param[in] j column index of the underlying matrix.
    \return The number of columns in the underlying matrices with column index
    \a j.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetN(int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::GetN()",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    return Nlocal_(j);
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetSize() const
  {
    return this->m_ * this->n_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDataSize() const
  {
    return nz_;
  }


  //! Returns size of A in bytes used to store the matrix.
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  inline size_t 
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetMemorySize() const
  {
    // to implement
    return 0;
  }



  //! Returns the type of a given underlying matrix.
  /*!
    Type 0 refers to a float dense matrice.
    Type 1 refers to a float sparse matrice.
    Type 2 refers to a double dense matrice.
    Type 3 refers to a double sparse matrice.
    \param[in] i row of the given underlying matrix.
    \param[in] j column of the given underlying matrix.
    \return The type of the underlying matrix.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator>
  inline int
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetType(int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("HeterogeneousMatrixCollection::GetType()",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("HeterogeneousMatrixCollection::GetType()",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif
    return collection_(i, j);
  }



  //! Returns the collection of float dense underlying matrices.
  /*!
    \return the collection of float dense underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::float_dense_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetFloatDense()
  {
    return float_dense_c_;
  }


  //! Returns the collection of float dense underlying matrices.
  /*!
    \return the collection of float dense underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline const typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::float_dense_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetFloatDense() const
  {
    return float_dense_c_;
  }


  //! Returns the collection of float sparse underlying matrices.
  /*!
    \return the collection of float sparse underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::float_sparse_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetFloatSparse()
  {
    return float_sparse_c_;
  }


  //! Returns the collection of float sparse underlying matrices.
  /*!
    \return the collection of float sparse underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline const typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::float_sparse_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetFloatSparse() const
  {
    return float_sparse_c_;
  }


  //! Returns the collection of double dense underlying matrices.
  /*!
    \return the collection of double dense underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::double_dense_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDoubleDense()
  {
    return double_dense_c_;
  }


  //! Returns the collection of double dense underlying matrices.
  /*!
    \return the collection of double dense underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline const typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::double_dense_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDoubleDense() const
  {
    return double_dense_c_;
  }


  //! Returns the collection of double sparse underlying matrices.
  /*!
    \return the collection of double sparse underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::double_sparse_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDoubleSparse()
  {
    return double_sparse_c_;
  }


  //! Returns the collection of double sparse underlying matrices.
  /*!
    \return the collection of double sparse underlying matrices.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> inline const typename
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::double_sparse_c&
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::GetDoubleSparse() const
  {
    return double_sparse_c_;
  }


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <template <class U> class Allocator>
  inline
  Matrix<FloatDouble, General, DenseSparseCollection, Allocator<double> >
  ::Matrix():
    HeterogeneousMatrixCollection<General, RowMajor, General,
				  RowSparse, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <template <class U> class Allocator>
  inline
  Matrix<FloatDouble, General, DenseSparseCollection, Allocator<double> >
  ::Matrix(int i, int j):
    HeterogeneousMatrixCollection<General, RowMajor, General,
				  RowSparse, Allocator>(i, j)
  {
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_HETEROGENEOUS_COLLECTION_INLINE_CXX
#endif
