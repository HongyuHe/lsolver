// Copyright (C) 2012 INRIA
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


#ifndef SELDON_FILE_MATRIX_PETSCMATRIX_INLINE_CXX

#include "PetscMatrix.hxx"


namespace Seldon
{


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline PetscMatrix<T, Prop, Storage, Allocator>::PetscMatrix():
    Matrix_Base<T, Allocator>()
  {
    mpi_communicator_ = MPI_COMM_WORLD;
    MatCreate(PETSC_COMM_WORLD, &petsc_matrix_);
    petsc_matrix_deallocated_ = false;
  }


  //! Main constructor.
  /*! Builds a i x j full matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline PetscMatrix<T, Prop, Storage, Allocator>::PetscMatrix(int i, int j):
    Matrix_Base<T, Allocator>(i, j)
  {
    int ierr;
    mpi_communicator_ = MPI_COMM_WORLD;
    MatCreate(mpi_communicator_, &petsc_matrix_);
    petsc_matrix_deallocated_ = false;
  }


  //! Copy constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline PetscMatrix<T, Prop, Storage, Allocator>
  ::PetscMatrix(Mat& A): Matrix_Base<T, Allocator>()
  {
    Copy(A);
  }


  //! Sets the MPI communicator.
  /*!
    \param[in] mpi_communicator the mpi communicator to be set.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>
  ::SetCommunicator(MPI_Comm mpi_communicator)
  {
    mpi_communicator_ = mpi_communicator;
  }


  //! Returns the MPI communicator of the current PETSc matrix.
  /*!
   \return The MPI communicator of the current PETSc matrix.
   */
  template <class T, class Prop, class Storage, class Allocator>
  inline MPI_Comm PetscMatrix<T, Prop, Storage, Allocator>
  ::GetCommunicator() const
  {
    return mpi_communicator_;
  }


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline PetscMatrix<T, Prop, Storage, Allocator>
  ::~PetscMatrix()
  {
    Clear();
  }


  //! Clears the matrix.
  /*!
    Destructs the matrix.
    \warning On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>
  ::Clear()
  {
    if (petsc_matrix_deallocated_)
      return;
    int ierr;
    ierr = MatDestroy(&petsc_matrix_);
    CHKERRABORT(mpi_communicator_, ierr);
    petsc_matrix_deallocated_ = true;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
    \warning Memory is not released except for me_.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>
  ::Nullify()
  {
    throw Undefined("PetscMatrix<T, Prop, Storage, Allocator>"
                    "::Nullify()");
  }


  //! Returns a reference on the inner petsc matrix.
  /*!
    \return A reference on the inner petsc matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Mat& PetscMatrix<T, Prop, Storage, Allocator>
  ::GetPetscMatrix()
  {
    return petsc_matrix_;
  }


  //! Returns a const reference on the inner petsc matrix.
  /*!
    \return A const reference on the inner petsc matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const Mat& PetscMatrix<T, Prop, Storage, Allocator>
  ::GetPetscMatrix() const
  {
    return petsc_matrix_;
  }


  //! Returns size of A in bytes used to store the matrix.
  template <class T, class Prop, class Storage, class Allocator>
  inline size_t PetscMatrix<T, Prop, Storage, Allocator>
  ::GetMemorySize() const
  {
    return 0;
  }
  

  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning The previous entries are kept, extra-entries may not be
    initialized (depending of the allocator).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>
  ::Resize(int i, int j)
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::Resize(int i, int j)");
  }


  //! Changes the size of the matrix and sets its data array
  //! (low level method).
  /*!
    The matrix is first cleared (memory is freed). The matrix is then resized
    to a i x j matrix, and the data array of the matrix is set to 'data'.
    'data' elements are not duplicated: the new data array of the matrix is
    the 'data' array. It is useful to create a matrix from pre-existing data.
    \param i new number of rows.
    \param j new number of columns.
    \param data new array storing elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The matrix
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>
  ::SetData(int i, int j, pointer data)
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::SetData(int i, int j, pointer data)");
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename PetscMatrix<T, Prop, Storage, Allocator>::const_reference
  PetscMatrix<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {
    throw Undefined("PetscMatrix::Val(int i, int j) const");
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename PetscMatrix<T, Prop, Storage, Allocator>::reference
  PetscMatrix<T, Prop, Storage, Allocator>::Val(int i, int j)
  {
    throw Undefined("PetscMatrix::Val(int i, int j)");
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename PetscMatrix<T, Prop, Storage, Allocator>::reference
  PetscMatrix<T, Prop, Storage, Allocator>::operator[] (int i)
  {
    throw Undefined("PetscMatrix::operator[] (int i)");
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename PetscMatrix<T, Prop, Storage, Allocator>::const_reference
  PetscMatrix<T, Prop, Storage, Allocator>::operator[] (int i) const
  {
    throw Undefined("PetscMatrix::operator[] (int i) const");
  }


  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>::Set(int i, int j, T value)
  {
    throw Undefined("PetscMatrix::Set(int i, int j, T value)");
  }


  //! Inserts or adds values into certain locations of a matrix.
  /*! \warning These values may be cached, so 'Flush' must be called after all
    calls to SetBuffer() have been completed.
    \param[in] i row index where to insert the value.
    \param[in] i column index where to insert the value.
    \param[in] value the value to insert.
    \param[in] insert_mode either INSERT_VALUES or ADD_VALUES, where
    ADD_VALUES adds the value to the entry, and INSERT_VALUES replaces
    existing entry with new value.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>
  ::SetBuffer(int i, int j, T value, InsertMode insert_mode = INSERT_VALUES)
  {
    int ierr;
    ierr = MatSetValues(petsc_matrix_, 1, &i, 1,
                        &j, &value, insert_mode);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Assembles the PETSc matrix.
  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>::Flush() const
  {
    int ierr;
    ierr = MatAssemblyBegin(petsc_matrix_, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(mpi_communicator_, ierr);
    ierr = MatAssemblyEnd(petsc_matrix_, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Returns the range of row indices owned by this processor.
  /*! The matrix is laid out with the first \f$n_1\f$ rows on the first
    processor, next \f$n_2\f$ rows on the second, etc. If the current
    processor is \f$k\f$, this method returns \f$n_k\f$ in \a i and
    \f$n_{k+1}\f$ in \a j. If \a i is set to PETSC_NULL on entry, it is not
    modified by this function. Same is true for \a j.
    \param[in,out] i the index of the first local row.
    \param[in,out] j the index of the last local row, plus 1.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void PetscMatrix<T, Prop, Storage, Allocator>
  ::GetProcessorRowRange(int& i, int& j) const
  {
    int ierr;
    ierr = MatGetOwnershipRange(petsc_matrix_, &i, &j);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  ///////////////////////////
  // Matrix<PETScSeqDense> //
  ///////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScSeqDense, Allocator>::Matrix():
    PetscMatrix<T, Prop, RowMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScSeqDense, Allocator>::Matrix(int i, int j):
    PetscMatrix<T, Prop, RowMajor, Allocator>(i, j)
  {
    MatSetType(this->petsc_matrix_, MATSEQDENSE);
    MatSetSizes(this->petsc_matrix_, PETSC_DECIDE, PETSC_DECIDE, i, j);
  }


  //! Copy constructor.
  /*!
    \param[in] A PETCs matrix to be copied.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScSeqDense, Allocator>::Matrix(Mat& A):
  PetscMatrix<T, Prop, RowMajor, Allocator>(A)
  {
  }


  //! Copy constructor.
  /*!
    \param[in] A matrix to be copied.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScSeqDense, Allocator>
  ::Matrix(const Matrix<T, Prop, PETScSeqDense, Allocator>& A)
  {
    this->petsc_matrix_deallocated_ = true;
    Copy(A);
  }


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x i matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, PETScSeqDense, Allocator>
  ::Reallocate(int i, int j)
  {
    this->Clear();
    int ierr;
    MatCreate(MPI_COMM_SELF, &this->petsc_matrix_);
    MatSetSizes(this->petsc_matrix_, i, j, i, j);
    MatSetType(this->petsc_matrix_, MATSEQDENSE);
    this->petsc_matrix_deallocated_ = false;
    this->Flush();
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, PETScSeqDense, Allocator>::value_type
  Matrix<T, Prop, PETScSeqDense, Allocator>::operator() (int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("PetscMatrix<PETScSeqDense>::operator()",
                     string("Index should be in [0, ")
                     + to_str(this->m_-1)
                     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("PetscMatrix<PETScSeqDense>::operator()",
                     string("Index should be in [0, ") + to_str(this->n_-1)
                     + "], but is equal to " + to_str(j) + ".");
#endif
    PetscInt idxm[1] = {i};
    PetscInt idxn[1] = {j};
    PetscScalar v[1];
    MatGetValues(this->petsc_matrix_, 1, idxm, 1, idxn, v);
    return v[0];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, PETScSeqDense, Allocator>::value_type
  Matrix<T, Prop, PETScSeqDense, Allocator>::operator() (int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("PetscMatrix<PETScSeqDense>::operator()",
                     string("Index should be in [0, ")
                     + to_str(this->m_-1)
                     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("PetscMatrix<PETScSeqDense>::operator()",
                     string("Index should be in [0, ")
                     + to_str(this->n_-1)
                     + "], but is equal to " + to_str(j) + ".");
#endif
    PetscInt idxm[1] = {i};
    PetscInt idxn[1] = {j};
    PetscScalar v[1];
    MatGetValues(this->petsc_matrix_, 1, idxm, 1, idxn, v);
    return v[0];
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, PETScSeqDense, Allocator>
  ::Copy(const Matrix<T, Prop, PETScSeqDense, Allocator>& A)
  {
    this->mpi_communicator_ = A.mpi_communicator_;
    PetscMatrix<T, Prop, RowMajor, Allocator>::Copy(A.GetPetscMatrix());
    this->petsc_matrix_deallocated_ = false;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScSeqDense, Allocator>&
  Matrix<T, Prop, PETScSeqDense, Allocator>
  ::operator= (const Matrix<T, Prop, PETScSeqDense, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  ///////////////////////////
  // Matrix<PETScMPIDense> //
  ///////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIDense, Allocator>::Matrix():
    PetscMatrix<T, Prop, RowMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIDense, Allocator>::Matrix(int i, int j):
    PetscMatrix<T, Prop, RowMajor, Allocator>(i, j)
  {
    MatSetType(this->petsc_matrix_, MATMPIDENSE);
    MatSetSizes(this->petsc_matrix_, PETSC_DECIDE, PETSC_DECIDE, i, j);
  }


  //! Copy constructor.
  /*!
    \param[in] A PETCs matrix to be copied.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIDense, Allocator>::Matrix(Mat& A):
  PetscMatrix<T, Prop, RowMajor, Allocator>(A)
  {
  }


  //! Copy constructor.
  /*!
    \param[in] A matrix to be copied.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIDense, Allocator>
  ::Matrix(const Matrix<T, Prop, PETScMPIDense, Allocator>& A)
  {
    this->petsc_matrix_deallocated_ = true;
    Copy(A);
  }


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x i matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, PETScMPIDense, Allocator>
  ::Reallocate(int i, int j, int local_i, int local_j)
  {
    this->Clear();
    int ierr;
    MatCreate(this->mpi_communicator_, &this->petsc_matrix_);
    MatSetType(this->petsc_matrix_, MATMPIDENSE);
    MatSetSizes(this->petsc_matrix_, local_i, local_j, i, j);
    this->m_ = i;
    this->n_ = j;
    this->petsc_matrix_deallocated_ = false;
    this->Flush();
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, PETScMPIDense, Allocator>::value_type
  Matrix<T, Prop, PETScMPIDense, Allocator>::operator() (int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    int start, end;
    this->GetProcessorRowRange(start, end);
    if (i < start || i >= end)
      throw WrongRow("PetscMatrix<PETScMPIDense>::operator()",
                     string("Index should be in [")
                     + to_str(start)
                     + ", " + to_str(end - 1)
                     + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("PetscMatrix<PETScMPIDense>::operator()",
                     string("Index should be in [0, ") + to_str(this->n_-1)
                     + "], but is equal to " + to_str(j) + ".");
#endif
    PetscInt idxm[1] = {i};
    PetscInt idxn[1] = {j};
    PetscScalar v[1];
    MatGetValues(this->petsc_matrix_, 1, idxm, 1, idxn, v);
    return v[0];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, PETScMPIDense, Allocator>::value_type
  Matrix<T, Prop, PETScMPIDense, Allocator>::operator() (int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    int start, end;
    this->GetProcessorRowRange(start, end);
    if (i < start || i >= end)
      throw WrongRow("PetscMatrix<PETScMPIDense>::operator()",
                     string("Index should be in [")
                     + to_str(start)
                     + ", " + to_str(end - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("PetscMatrix<PETScMPIDense>::operator()",
                     string("Index should be in [0, ")
                     + to_str(this->n_-1)
                     + "], but is equal to " + to_str(j) + ".");
#endif
    PetscInt idxm[1] = {i};
    PetscInt idxn[1] = {j};
    PetscScalar v[1];
    MatGetValues(this->petsc_matrix_, 1, idxm, 1, idxn, v);
    return v[0];
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, PETScMPIDense, Allocator>
  ::Copy(const Matrix<T, Prop, PETScMPIDense, Allocator>& A)
  {
    this->mpi_communicator_ = A.mpi_communicator_;
    PetscMatrix<T, Prop, RowMajor, Allocator>::Copy(A.GetPetscMatrix());
    this->petsc_matrix_deallocated_ = false;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIDense, Allocator>&
  Matrix<T, Prop, PETScMPIDense, Allocator>
  ::operator= (const Matrix<T, Prop, PETScMPIDense, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  /////////////////////////
  // Matrix<PETScMPIAIJ> //
  /////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIAIJ, Allocator>::Matrix():
    PetscMatrix<T, Prop, RowMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIAIJ, Allocator>::Matrix(int i, int j):
    PetscMatrix<T, Prop, RowMajor, Allocator>(i, j)
  {
    MatSetType(this->petsc_matrix_, MATMPIAIJ);
    MatSetSizes(this->petsc_matrix_, PETSC_DECIDE, PETSC_DECIDE, i, j);
  }


  //! Copy constructor.
  /*!
    \param[in] A PETCs matrix to be copied.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIAIJ, Allocator>::Matrix(Mat& A):
  PetscMatrix<T, Prop, RowMajor, Allocator>(A)
  {
  }


  //! Copy constructor.
  /*!
    \param[in] A matrix to be copied.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIAIJ, Allocator>
  ::Matrix(const Matrix<T, Prop, PETScMPIAIJ, Allocator>& A)
  {
    this->petsc_matrix_deallocated_ = true;
    Copy(A);
  }


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x i matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, PETScMPIAIJ, Allocator>
  ::Reallocate(int i, int j, int local_i, int local_j)
  {
    this->Clear();
    int ierr;
    MatCreate(this->mpi_communicator_, &this->petsc_matrix_);
    MatSetType(this->petsc_matrix_, MATMPIAIJ);
    MatSetSizes(this->petsc_matrix_, local_i, local_j, i, j);
    this->m_ = i;
    this->n_ = j;
    this->petsc_matrix_deallocated_ = false;
    this->Flush();
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, PETScMPIAIJ, Allocator>::value_type
  Matrix<T, Prop, PETScMPIAIJ, Allocator>::operator() (int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    int start, end;
    this->GetProcessorRowRange(start, end);
    if (i < start || i >= end)
      throw WrongRow("PetscMatrix<PETScMPIAIJ>::operator()",
                     string("Index should be in [")
                     + to_str(start)
                     + ", " + to_str(end - 1)
                     + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("PetscMatrix<PETScMPIAIJ>::operator()",
                     string("Index should be in [0, ") + to_str(this->n_-1)
                     + "], but is equal to " + to_str(j) + ".");
#endif
    PetscInt idxm[1] = {i};
    PetscInt idxn[1] = {j};
    PetscScalar v[1];
    MatGetValues(this->petsc_matrix_, 1, idxm, 1, idxn, v);
    return v[0];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, PETScMPIAIJ, Allocator>::value_type
  Matrix<T, Prop, PETScMPIAIJ, Allocator>::operator() (int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    int start, end;
    this->GetProcessorRowRange(start, end);
    if (i < start || i >= end)
      throw WrongRow("PetscMatrix<PETScMPIAIJ>::operator()",
                     string("Index should be in [")
                     + to_str(start)
                     + ", " + to_str(end - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("PetscMatrix<PETScMPIAIJ>::operator()",
                     string("Index should be in [0, ")
                     + to_str(this->n_-1)
                     + "], but is equal to " + to_str(j) + ".");
#endif
    PetscInt idxm[1] = {i};
    PetscInt idxn[1] = {j};
    PetscScalar v[1];
    MatGetValues(this->petsc_matrix_, 1, idxm, 1, idxn, v);
    return v[0];
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, PETScMPIAIJ, Allocator>
  ::Copy(const Matrix<T, Prop, PETScMPIAIJ, Allocator>& A)
  {
    this->mpi_communicator_ = A.mpi_communicator_;
    PetscMatrix<T, Prop, RowMajor, Allocator>::Copy(A.GetPetscMatrix());
    this->petsc_matrix_deallocated_ = false;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, PETScMPIAIJ, Allocator>&
  Matrix<T, Prop, PETScMPIAIJ, Allocator>
  ::operator= (const Matrix<T, Prop, PETScMPIAIJ, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }
  
}


#define SELDON_FILE_MATRIX_PETSCMATRIX_INLINE_CXX
#endif
