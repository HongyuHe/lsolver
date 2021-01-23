// Copyright (C) 2011-2012, INRIA
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


#ifndef SELDON_FILE_VECTOR_PETSCVECTOR_CXX


#include "PetscVector.hxx"


namespace Seldon
{

  //! Duplicates a vector.
  /*!
    \param petsc_vector vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  void PETScVector<T, Allocator>
  ::Copy(const Vec& petsc_vector)
  {
    Clear();

    int ierr;
    ierr = VecGetSize(petsc_vector, &this->m_);
    CHKERRABORT(mpi_communicator_, ierr);
    PetscObjectGetComm(reinterpret_cast<PetscObject>(petsc_vector),
                       &mpi_communicator_);
    CHKERRABORT(mpi_communicator_, ierr);
    ierr = VecDuplicate(petsc_vector, &petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
    ierr = VecCopy(petsc_vector, petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
    petsc_vector_deallocated_ = false;
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the vector stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Allocator>
  void PETScVector<T, Allocator>::Zero()
  {
    int ierr;
    ierr = VecSet(petsc_vector_, 0.);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Fills the vector with 0, 1, 2, ...
  template <class T, class Allocator>
  void PETScVector<T, Allocator>::Fill()
  {
    int ierr;
    Vector<int> index(this->m_);
    Vector<double> data(this->m_);
    index.Fill();
    data.Fill();
    ierr = VecSetValues(petsc_vector_, this->m_, index.GetData(),
                        data.GetData(), INSERT_VALUES);
    CHKERRABORT(mpi_communicator_, ierr);
    Flush();
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  void PETScVector<T, Allocator>::Fill(const T0& x)
  {
    int ierr;
    ierr = VecSet(petsc_vector_, double(x));
    CHKERRABORT(mpi_communicator_, ierr);
    Flush();
  }


  //! Fills the vector randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Allocator>
  void PETScVector<T, Allocator>::FillRand()
  {
    int ierr;
    Vector<int> index(this->m_);
    Vector<double> data(this->m_);
    index.Fill();
    data.FillRand();
    ierr = VecSetValues(petsc_vector_, this->m_, index.GetData(),
                        data.GetData(), INSERT_VALUES);
    CHKERRABORT(mpi_communicator_, ierr);
    Flush();
  }


  /*********
   * NORMS *
   *********/


  //! Returns the infinite norm.
  /*!
    \return The infinite norm.
  */
  template <class T, class Allocator>
  typename PETScVector<T, Allocator>::value_type
  PETScVector<T, Allocator>::GetNormInf() const
  {
    int ierr;
    value_type res;
    ierr = VecNorm(petsc_vector_, NORM_INFINITY, &res);
    CHKERRABORT(mpi_communicator_, ierr);
    return res;
  }


  //! Returns the index of the highest absolute value.
  /*!
    \return The index of the element that has the highest absolute value.
  */
  template <class T, class Allocator>
  int PETScVector<T, Allocator>::GetNormInfIndex() const
  {
    throw Undefined("PETScVector<T, Allocator>::GetNormInfIndex()");
  }


  /////////////////////////////
  // SEQUENTIAL PETSC VECTOR //
  /////////////////////////////


  /**************************
   * OUTPUT/INPUT FUNCTIONS *
   **************************/


  //! Writes the vector in a file.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::Write(string FileName, bool with_size) const");
  }


  //! Writes the vector in a file stream.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::Write(string FileName, bool with_size) const");
  }


  //! Writes the vector in a file.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<T, PETScSeq, Allocator>"
                    "::WriteText(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif
    this->WriteText(FileStream);
    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>::WriteText(ostream& FileStream) const
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::WriteText(ostream& FileStream) const");
  }


  //! Sets the vector from a file.
  /*!
    Sets the vector according to a binary file that stores the length of the
    vector (integer) and all elements.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not
    available in the file. In this case, the current size N of the vector is
    unchanged, and N elements are read in the file.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>
  ::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<T, PETScSeq, Allocator>::Read(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif
    this->Read(FileStream, with_size);
    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets the vector according to a binary file stream that stores the length
    of the vector (integer) and all elements.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not
    available in the stream. In this case, the current size N of the vector is
    unchanged, and N elements are read in the stream.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>
  ::Read(istream& FileStream, bool with_size)
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::Read(istream& FileStream, bool with_size)");
  }


  //! Sets the vector from a file.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>::ReadText(string FileName)
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::ReadText(string FileName)");
  }


  //! Sets the vector from a file stream.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>::ReadText(istream& FileStream)
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::ReadText(istream& FileStream)");
  }


  //! operator<< overloaded for vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
                        const Vector<T, PETScSeq, Allocator>& V)
  {
    out << '\t';
    for (int i = 0; i < V.GetLength() - 1; i++)
      out << V(i) << '\t';
    if (V.GetLength() != 0)
      out << V(V.GetLength() - 1);
    return out;
  }


  ///////////////////////////
  // PARALLEL PETSC VECTOR //
  ///////////////////////////


  /**************************
   * OUTPUT/INPUT FUNCTIONS *
   **************************/


  //! Writes the vector in a file.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    throw Undefined("Vector<T, PETScPar, Allocator>"
                    "::Write(string FileName, bool with_size) const");
  }


  //! Writes the vector in a file stream.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {
    int local_n;
    VecGetLocalSize(this->petsc_vector_, &local_n);

    Vector<int> index(local_n);
    index.Fill();
    int i_start, i_end;
    this->GetProcessorRange(i_start, i_end);
    for (int i = 0; i < local_n; i++)
      index(i) += i_start;
    Vector<T> data(local_n);
    VecGetValues(this->petsc_vector_, local_n, index.GetData(),
                 data.GetData());
    data.Write(FileStream, with_size);
  }


  //! Writes the vector in a file.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<T, PETScPar, Allocator>"
                    "::WriteText(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif
    this->WriteText(FileStream);
    FileStream.close();

  }


  //! Writes the vector in a file stream.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>::WriteText(ostream& FileStream) const
  {
    int local_n;
    VecGetLocalSize(this->petsc_vector_, &local_n);

    Vector<int> index(local_n);
    index.Fill();
    int i_start, i_end;
    this->GetProcessorRange(i_start, i_end);
    for (int i = 0; i < local_n; i++)
      index(i) += i_start;
    Vector<T> data(local_n);
    VecGetValues(this->petsc_vector_, local_n, index.GetData(),
                 data.GetData());
    data.WriteText(FileStream);
  }


  //! Sets the vector from a file.
  /*!
    Sets the vector according to a binary file that stores the length of the
    vector (integer) and all elements.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not
    available in the file. In this case, the current size N of the vector is
    unchanged, and N elements are read in the file.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>
  ::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<T, PETScPar, Allocator>::Read(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif
    this->Read(FileStream, with_size);
    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets the vector according to a binary file stream that stores the length
    of the vector (integer) and all elements.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not
    available in the stream. In this case, the current size N of the vector is
    unchanged, and N elements are read in the stream.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>
  ::Read(istream& FileStream, bool with_size)
  {
    throw Undefined("PETScVector<T, PETScPar, Allocator>"
                    "::Read(istream& FileStream, bool with_size)");
  }


  //! Sets the vector from a file.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>::ReadText(string FileName)
  {
    throw Undefined("PETScVector<T, PETScPar, Allocator>"
                    "::ReadText(string FileName)");
  }


  //! Sets the vector from a file stream.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>::ReadText(istream& FileStream)
  {
    throw Undefined("PETScVector<T, PETScPar, Allocator>"
                    "::ReadText(istream& FileStream)");
  }


  //! operator<< overloaded for vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
                        const Vector<T, PETScPar, Allocator>& V)
  {
    out << '\t';
    for (int i = 0; i < V.GetLength() - 1; i++)
      out << V(i) << '\t';
    if (V.GetLength() != 0)
      out << V(V.GetLength() - 1);
    return out;
  }


} // namespace Seldon.


#define SELDON_FILE_PETSCVECTOR_CXX
#endif
