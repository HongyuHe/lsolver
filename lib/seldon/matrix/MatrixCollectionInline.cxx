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


#ifndef SELDON_FILE_MATRIX_COLLECTION_INLINE_CXX

#include "MatrixCollection.hxx"


namespace Seldon
{


  //////////////////////
  // MATRIXCOLLECTION //
  //////////////////////


  //! Copy constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline MatrixCollection<T, Prop, Storage, Allocator>
  ::MatrixCollection(const MatrixCollection<T, Prop, Storage, Allocator>& A)
    : Matrix_Base<T, Allocator>()
  {
    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline MatrixCollection<T, Prop, Storage, Allocator>::~MatrixCollection()
  {
    Clear();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of rows.
  /*!
    \return the total number of rows. It is the sum of the number of rows in
    the underlying matrices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetM() const
  {
    return this->m_;
  }


  //! Returns the number of rows.
  /*!
    \return the total number of rows. It is the sum of the number of rows
    in the underlying matrices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetMmatrix() const
  {
    return Mmatrix_;
  }


  //! Returns the number of rows in an underlying matrix.
  /*!
    \param[in] i row index of the underlying matrix.
    \return The number of rows in the underlying matrices with row index \a i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetM(int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::GetM()",
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
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetN() const
  {
    return this->n_;
  }


  //! Returns the number of columns.
  /*!
    \return the total number of columns. It is the sum of the number of
    columns in the underlying matrices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetNmatrix() const
  {
    return Nmatrix_;
  }


  //! Returns the number of columns in an underlying matrix.
  /*!
    \param[in] j column index of the underlying matrix.
    \return The number of columns in the underlying matrices with column index
    \a j.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetN(int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::GetN()",
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
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetSize() const
  {
    return this->m_ * this->n_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return nz_;
  }

  
  //! Returns size of A in bytes used to store the matrix.
  template <class T, class Prop, class Storage, class Allocator>
  inline size_t MatrixCollection<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    size_t taille = size_t(nz_)*sizeof(T);
    return taille;
  }
  

  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access to an underlying matrix.
  /*!
    Returns the underlying matrix (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return The matrix collection (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline
  typename MatrixCollection<T, Prop, Storage, Allocator>::matrix_reference
  MatrixCollection<T, Prop, Storage, Allocator>::GetMatrix(int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::GetMatrix(int, int)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::GetMatrix(int, int)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    return matrix_(i, j);
  }


  //! Access to an underlying matrix.
  /*!
    Returns the underlying matrix (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return The matrix collection (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename MatrixCollection<T, Prop, Storage, Allocator>
  ::const_matrix_reference
  MatrixCollection<T, Prop, Storage, Allocator>::GetMatrix(int i,
							   int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::GetMatrix(int, int)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::GetMatrix(int, int)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    return matrix_(i, j);
  }

  
  ////////////////////////
  // COLMAJORCOLLECTION //
  ////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajorCollection, Allocator>::Matrix():
    MatrixCollection<T, Prop, ColMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajorCollection, Allocator>::Matrix(int i, int j):
    MatrixCollection<T, Prop, ColMajor, Allocator>(i, j)
  {
  }


  ////////////////////////
  // ROWMAJORCOLLECTION //
  ////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajorCollection, Allocator>::Matrix():
    MatrixCollection<T, Prop, RowMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajorCollection, Allocator>::Matrix(int i, int j):
    MatrixCollection<T, Prop, RowMajor, Allocator>(i, j)
  {
  }


  ////////////////////////////
  // COLSYMPACKEDCOLLECTION //
  ////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymPackedCollection, Allocator>::Matrix():
    MatrixCollection<T, Prop, ColSymPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymPackedCollection, Allocator>
  ::Matrix(int i, int j):
    MatrixCollection<T, Prop, ColSymPacked, Allocator>(i, j)
  {
  }


  ////////////////////////////
  // ROWSYMPACKEDCOLLECTION //
  ////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymPackedCollection, Allocator>::Matrix():
    MatrixCollection<T, Prop, RowSymPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymPackedCollection, Allocator>
  ::Matrix(int i, int j):
    MatrixCollection<T, Prop, RowSymPacked, Allocator>(i, j)
  {
  }


} // namespace Seldon.


#define SELDON_FILE_MATRIX_COLLECTION_INLINE_CXX
#endif
