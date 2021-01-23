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


#ifndef SELDON_FILE_MATRIX_PETSCMATRIX_CXX

#include "PetscMatrix.hxx"


namespace Seldon
{
  
  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>
  ::Copy(const Mat& A)
  {
    Clear();
    int ierr;
    ierr = MatDuplicate(A, MAT_COPY_VALUES, &petsc_matrix_);
    CHKERRABORT(mpi_communicator_, ierr);
    MatGetSize(A, &this->m_, &this->n_);
    mpi_communicator_ = MPI_COMM_WORLD;
    petsc_matrix_deallocated_ = false;
  }


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the matrix stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>::Zero()
  {
    MatScale(petsc_matrix_, T(0));
  }


  //! Sets the matrix to the identity.
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>::SetIdentity()
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::SetIdentity()");
  }


  //! Fills the matrix with 0, 1, 2, ...
  /*!
    On exit, the matrix is filled with 0, 1, 2, 3, ... The order of
    those numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>::Fill()
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::Fill()");
  }


  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  void PetscMatrix<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::Fill(const T0& x)");
  }


  //! Fills a matrix randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>::FillRand()
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::FillRand()");
  }


  //! Displays a sub-matrix on the standard output.
  /*!
    The sub-matrix is defined by its upper-left corner (a, b)
    and its bottom-right corner (m, n). So, elements with indices
    in [a, m] x [b, n] are displayed on the standard output,
    in text format. Each row is displayed on a single line and
    elements of a row are delimited by tabulations.
    \param a row index of the upper-left corner.
    \param b column index of the upper-left corner.
    \param m row index of the bottom-right corner.
    \param n column index of the bottom-right corner.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>
  ::Print(int a, int b, int m, int n) const
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::Print(int a, int b, int m, int n) const");
  }


  //! Displays a square sub-matrix on the standard output.
  /*!
    The sub-matrix is defined by its bottom-right corner (l, l).
    So, elements with indices in [0, 0] x [l, l] are displayed
    on the standard output, in text format. Each row is displayed
    on a single line and elements of a row are delimited
    by tabulations.
    \param l dimension of the square matrix to be displayed.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>::Print(int l) const
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::Print(int l) const");
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("PetscMatrix::Write(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::Write(ostream& FileStream, bool with_size) const");
  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in text format.
    Only matrix elements are written (not dimensions).
    Each row is written on a single line and elements of
    a row are delimited by tabulations.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("PetscMatrix::WriteText(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in text format.
    Only matrix elements are written (not dimensions).
    Each row is written on a single line and elements of
    a row are delimited by tabulations.
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::WriteText(ostream& FileStream) const");
  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("PetscMatrix::Read(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream)
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::Read(istream& FileStream) const");
  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in text format in a file.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>::ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::ReadText(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in text format from an input stream.
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void PetscMatrix<T, Prop, Storage, Allocator>
  ::ReadText(istream& FileStream)
  {
    throw Undefined("void PetscMatrix<T, Prop, Storage, Allocator>"
                    "::ReadText(istream& FileStream)");
  }


  ///////////////////////////
  // Matrix<PETScSeqDense> //
  ///////////////////////////


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, PETScSeqDense, Allocator>::Print() const
  {
    int ierr;
    ierr = MatView(this->petsc_matrix_, PETSC_VIEWER_STDOUT_SELF);
    CHKERRABORT(this->mpi_communicator_, ierr);
  }


  ///////////////////////////
  // Matrix<PETScMPIDense> //
  ///////////////////////////

  
  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, PETScMPIDense, Allocator>::Print() const
  {
    int ierr;
    ierr = MatView(this->petsc_matrix_, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRABORT(this->mpi_communicator_, ierr);
  }


  /////////////////////////
  // Matrix<PETScMPIAIJ> //
  /////////////////////////


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  template <class T0,  class Allocator0>
  void Matrix<T, Prop, PETScMPIAIJ, Allocator>
  ::Copy(const Matrix<T0, General, RowSparse, Allocator0>& A)
  {
    this->Clear();

    int ierr;
    ierr = MatCreate(this->mpi_communicator_, &this->petsc_matrix_);
    CHKERRABORT(this->mpi_communicator_, ierr);

    int ma = A.GetM();
    int na = A.GetN();
    int nnz = A.GetDataSize();
    double *value = A.GetData();
    int *column = A.GetInd();
    int *ptr = A.GetPtr();

    this->m_ = ma;
    this->n_ = na;

    ierr = MatSetType(this->petsc_matrix_, MATMPIAIJ);
    CHKERRABORT(this->mpi_communicator_, ierr);
    ierr = MatSetSizes(this->petsc_matrix_, PETSC_DECIDE, PETSC_DECIDE,
                       this->m_, this->n_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    int loc_start, loc_end;
    ierr = MatGetOwnershipRange(this->petsc_matrix_,
                                &loc_start, &loc_end);
    CHKERRABORT(this->mpi_communicator_, ierr);
    for (int i = loc_start; i < loc_end; i++)
      for (int j = ptr[i]; j < ptr[i + 1]; j++)
        ierr = MatSetValues(this->petsc_matrix_, 1, &i, 1, &column[j],
                            &value[j], INSERT_VALUES);
    ierr = MatAssemblyBegin(this->petsc_matrix_,MAT_FINAL_ASSEMBLY);
    CHKERRABORT(this->mpi_communicator_, ierr);
    ierr = MatAssemblyEnd(this->petsc_matrix_,MAT_FINAL_ASSEMBLY);
    CHKERRABORT(this->mpi_communicator_, ierr);
  }


  //! Returns the number of local rows of the inner PETSc matrix.
  /*!
    \return The number of local rows of the inner PETSc matrix.
   */
  template <class T, class Prop, class Allocator>
  int Matrix<T, Prop, PETScMPIAIJ, Allocator>::GetLocalM() const
  {
    int ierr;
    int m;
    ierr = MatGetLocalSize(this->petsc_matrix_, &m, PETSC_NULL);
    CHKERRABORT(this->mpi_communicator_, ierr);
    return m;
  }


  //! Gets the number of local columns of the inner PETSc matrix.
  /*!
    \return The number of local columns of the inner PETSc matrix.
   */
  template <class T, class Prop, class Allocator>
  int Matrix<T, Prop, PETScMPIAIJ, Allocator>::GetLocalN() const
  {
    int ierr;
    int n;
    ierr = MatGetLocalSize(this->petsc_matrix_, PETSC_NULL, &n);
    CHKERRABORT(this->mpi_communicator_, ierr);
    return n;
  }


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, PETScMPIAIJ, Allocator>::Print() const
  {
    int ierr;
    ierr = MatView(this->petsc_matrix_, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRABORT(this->mpi_communicator_, ierr);
  }


}


#define SELDON_FILE_MATRIX_PETSCMATRIX_CXX
#endif
