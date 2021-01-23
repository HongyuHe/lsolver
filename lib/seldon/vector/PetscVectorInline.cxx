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


#ifndef SELDON_FILE_VECTOR_PETSCVECTOR_INLINE_CXX


#include "PetscVector.hxx"


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
  inline PETScVector<T, Allocator>::PETScVector():
    Vector_Base<T, Allocator>()
  {
    this->m_ = 0;
    petsc_vector_deallocated_ = false;
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param[in] i length of the vector.
    \param[in] mpi_communicator MPI communicator to use.
  */
  template <class T, class Allocator>
  inline PETScVector<T, Allocator>::PETScVector(int i, MPI_Comm mpi_communicator)
    :Vector_Base<T, Allocator>(i)
  {
    this->m_ = i;
    petsc_vector_deallocated_ = false;
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param[in] petsc_vector vector to be copied.
  */
  template <class T, class Allocator>
  inline PETScVector<T, Allocator>::
  PETScVector(Vec& petsc_vector)
  {
    petsc_vector_deallocated_ = true;
    Copy(petsc_vector);
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param[in] V vector to be copied.
  */
  template <class T, class Allocator>
  inline PETScVector<T, Allocator>::
  PETScVector(const PETScVector<T, Allocator>& V)
  {
    Copy(V);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline PETScVector<T, Allocator>::~PETScVector()
  {
    Clear();
  }


  //! Returns a reference on the inner petsc vector.
  /*!
    \return a reference on the inner petsc vector.
  */
  template <class T, class Allocator>
  inline Vec& PETScVector<T, Allocator>::GetPetscVector()
  {
    return petsc_vector_;
  }


  //! Returns a const reference on the inner petsc vector.
  /*!
    \return a const reference on the inner petsc vector.
  */
  template <class T, class Allocator>
  inline const Vec& PETScVector<T, Allocator>::GetPetscVector() const
  {
    return petsc_vector_;
  }


  //! Sets the MPI communicator.
  /*!
    \param[in] mpi_communicator the mpi communicator to be set.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>::SetCommunicator(MPI_Comm mpi_communicator)
  {
    mpi_communicator_ = mpi_communicator;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Clears the vector.
  /*!
    Destructs the vector.
    \warning On exit, the vector is an empty vector.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>::Clear()
  {
    if (petsc_vector_deallocated_)
      return;
    int ierr;
    ierr = VecDestroy(&petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
    petsc_vector_deallocated_ = true;
  }


  //! Changes the length of the vector, and keeps previous values.
  /*!
    Reallocates the vector to size i. Previous values are kept.
    \param[in] n new length of the vector.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>::Resize(int n)
  {
    throw Undefined("PETScVector<T, Allocator>::Resize(int n)");
  }


  /*! \brief Changes the length of the vector and sets its data array (low
    level method). */
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param[in] i new length of the vector.
    \param[in] data the new data array. 'data' contains the new elements of
    the vector and must therefore contain 'i' elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The vector
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::SetData(int i, typename PETScVector<T, Allocator>::pointer data)
  {
    throw Undefined("PETScVector<T, Allocator>::SetData(int i, "
                    "typename PETScVector<T, Allocator>::pointer data)");
  }


  //! Clears the vector without releasing memory.
  /*!
    On exit, the vector is empty and the memory has not been released.
    It is useful for low level manipulations on a Vector instance.
    \warning Memory is not released.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>::Nullify()
  {
    throw Undefined("PETScVector<T, Allocator>::Nullify()");
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator>
  inline typename PETScVector<T, Allocator>::value_type
  PETScVector<T, Allocator>::operator() (int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<PETSc>::operator()",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif
    int ierr;
    value_type ret[1];
    int index[1];
    index[0] = i;
    ierr = VecGetValues(petsc_vector_, 1, index, ret);
    CHKERRABORT(mpi_communicator_, ierr);
    return ret[0];
  }


  //! Inserts or adds values into certain locations of a vector.
  /*! \warning These values may be cached, so 'Flush' must be called after
    all calls to SetBuffer() have been completed.
    \param[in] i index where to insert the value.
    \param[in] value the value to insert.
    \param[in] insert_mode either INSERT_VALUES or ADD_VALUES, where
    ADD_VALUES adds the value to the entry, and INSERT_VALUES replaces
    existing entry with new value.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::SetBuffer(int i, T value, InsertMode insert_mode)
  {
    int ierr;
    int ix[1] = {i};
    T data[1] = {value};
    ierr = VecSetValues(petsc_vector_, 1, ix,
                        data, INSERT_VALUES);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Assembles the PETSc vector.
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::Flush()
  {
    int ierr;
    ierr = VecAssemblyBegin(petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
    ierr = VecAssemblyEnd(petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Returns the range of indices owned by this processor.
  /*! The vectors are laid out with the first \f$n_1\f$ elements on the first
    processor, next \f$n_2\f$ elements on the second, etc. If the current
    processor is \f$k\f$, this method returns \f$n_k\f$ in \a i and
    \f$n_{k+1}\f$ in \a j. If \a i is set to PETSC_NULL on entry, it is not
    modified by this function. Same is true for \a j.
    \param[in,out] i the index of the first local element.
    \param[in,out] j the index of the last local element, plus 1.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::GetProcessorRange(int& i, int& j) const
  {
    int ierr;
    ierr = VecGetOwnershipRange(petsc_vector_, &i, &j);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::Copy(const PETScVector<T, Allocator>& X)
  {
    Copy(X.GetPetscVector());
  }


  //! Appends an element to the vector.
  /*!
    \param x element to be appended.
    \warning This method will only work if the allocator preserves the
    elements while reallocating.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>::Append(const T& x)
  {
    throw Undefined("PETScVector<T, Allocator>::Append(const T& x)");
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Allocator>
  inline int PETScVector<T, Allocator>::GetDataSize() const
  {
    return this->m_;
  }


  //! Returns the number of elements stored.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Allocator>
  inline int PETScVector<T, Allocator>::GetLocalM() const
  {
    int size;
    VecGetLocalSize(petsc_vector_, &size);
    return size;
  }


  /////////////////////////////
  // SEQUENTIAL PETSC VECTOR //
  /////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScSeq, Allocator>::Vector():
    PETScVector<T, Allocator>()
  {
    this->mpi_communicator_ = MPI_COMM_WORLD;
    this->m_ = 0;
    int ierr;
    ierr = VecCreateSeq(PETSC_COMM_SELF, 0, &this->petsc_vector_);
    this->petsc_vector_deallocated_ = false;
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param[in] i length of the vector.
    \param[in] mpi_communicator MPI communicator to use.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScSeq, Allocator>::Vector(int i, MPI_Comm mpi_communicator)
    :PETScVector<T, Allocator>(i)
  {
    int ierr;
    this->mpi_communicator_ = mpi_communicator;
    this->m_ = i;
    ierr = VecCreateSeq(PETSC_COMM_SELF, i, &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    ierr = VecSet(this->petsc_vector_, 0.);
    CHKERRABORT(this->mpi_communicator_, ierr);
    this->Flush();
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param|in] petsc_vector vector to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScSeq, Allocator>::
  Vector(Vec& petsc_vector): PETScVector<T, Allocator>(petsc_vector)
  {
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param V vector to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScSeq, Allocator>::
  Vector(const Vector<T, PETScSeq, Allocator>& V)
  {
    Copy(V);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Vector<T, PETScSeq, Allocator>::~Vector()
  {
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, PETScSeq, Allocator>
  ::Copy(const Vector<T, PETScSeq, Allocator>& X)
  {
    Copy(X.GetPetscVector());
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, PETScSeq, Allocator>
  ::Copy(const Vec& X)
  {
    PETScVector<T, Allocator>::Copy(X);
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, initial elements of the vector may
    be lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, PETScSeq, Allocator>
  ::Reallocate(int i)
  {
    this->Clear();
    int ierr;
    this->m_ = i;
    ierr = VecCreateSeq(PETSC_COMM_SELF, i, &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    Fill(T(0));
    this->Flush();
    this->petsc_vector_deallocated_ = false;
  }


  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScSeq, Allocator>& Vector<T, PETScSeq, Allocator>
  ::operator= (const Vector<T, PETScSeq, Allocator>& X)
  {
    this->Copy(X);
    return *this;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  inline Vector<T, PETScSeq, Allocator>&
  Vector<T, PETScSeq, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);
    return *this;
  }


  //! Multiplies a vector by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Allocator>
  template<class T0>
  inline Vector<T, PETScSeq, Allocator>& Vector<T, PETScSeq, Allocator>
  ::operator*= (const T0& alpha)
  {
    int ierr;
    ierr = VecScale(this->petsc_vector_, alpha);
    CHKERRABORT(this->mpi_communicator_, ierr);
    return *this;
  }


  //! Displays the vector.
  template <class T, class Allocator>
  inline void Vector<T, PETScSeq, Allocator>::Print() const
  {
    int ierr;
    ierr = VecView(this->petsc_vector_, PETSC_VIEWER_STDOUT_SELF);
  }


  ///////////////////////////
  // PARALLEL PETSC VECTOR //
  ///////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScPar, Allocator>::Vector():
    PETScVector<T, Allocator>()
  {
    this->mpi_communicator_ = MPI_COMM_WORLD;
    int ierr;
    ierr = VecCreateMPI(this->mpi_communicator_,
                        PETSC_DECIDE, 0, &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param[in] i length of the vector.
    \param[in] mpi_communicator MPI communicator to use.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScPar, Allocator>::Vector(int i, MPI_Comm mpi_communicator):
    PETScVector<T, Allocator>(i)
  {
    int ierr;
    this->mpi_communicator_ = mpi_communicator;
    ierr = VecCreateMPI(this->mpi_communicator_, PETSC_DECIDE, i,
                        &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    ierr = VecSet(this->petsc_vector_, 0.);
    CHKERRABORT(this->mpi_communicator_, ierr);
    this->Flush();
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param[in] i length of the vector.
    \param[in] Nlocal size of local vector.
    \param[in] mpi_communicator MPI communicator to use.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScPar, Allocator>::Vector(int i, int Nlocal,
						MPI_Comm mpi_communicator):
    PETScVector<T, Allocator>(i)
  {
    int ierr;
    this->mpi_communicator_ = mpi_communicator;
    ierr = VecCreateMPI(this->mpi_communicator_, Nlocal,
                        i, &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    ierr = VecSet(this->petsc_vector_, 0.);
    CHKERRABORT(this->mpi_communicator_, ierr);
    this->Flush();
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param[in] petsc_vector vector to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScPar, Allocator>::
  Vector(Vec& petsc_vector): PETScVector<T, Allocator>(petsc_vector)
  {
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param[in] V vector to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScPar, Allocator>::
  Vector(const Vector<T, PETScPar, Allocator>& V)
  {
    Copy(V);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Vector<T, PETScPar, Allocator>::~Vector()
  {
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, PETScPar, Allocator>
  ::Copy(const Vector<T, PETScPar, Allocator>& X)
  {
    Copy(X.GetPetscVector());
    this->mpi_communicator_ = X.mpi_communicator_;
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, PETScPar, Allocator>
  ::Copy(const Vec& X)
  {
    PETScVector<T, Allocator>::Copy(X);
  }


  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScPar, Allocator>& Vector<T, PETScPar, Allocator>
  ::operator= (const Vector<T, PETScPar, Allocator>& X)
  {
    this->Copy(X);
    return *this;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  Vector<T, PETScPar, Allocator>&
  Vector<T, PETScPar, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);
    return *this;
  }


  //! Multiplies a vector by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Allocator>
  template<class T0>
  inline Vector<T, PETScPar, Allocator>& Vector<T, PETScPar, Allocator>
  ::operator*= (const T0& alpha)
  {
    int ierr;
    ierr = VecScale(this->petsc_vector_, alpha);
    CHKERRABORT(this->mpi_communicator_, ierr);
    return *this;
  }


  //! Displays the vector.
  template <class T, class Allocator>
  inline void Vector<T, PETScPar, Allocator>::Print() const
  {
    int ierr;
    ierr = VecView(this->petsc_vector_, PETSC_VIEWER_STDOUT_WORLD);
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, initial elements of the vector may
    be lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, PETScPar, Allocator>
  ::Reallocate(int i, int local_size)
  {
    this->Clear();
    int ierr;
    this->m_ = i;
    ierr = VecCreateMPI(this->mpi_communicator_, local_size, i,
                        &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    Fill(T(0));
    this->Flush();
    this->petsc_vector_deallocated_ = false;
  }

  
} // namespace Seldon.


#define SELDON_FILE_PETSCVECTOR_INLINE_CXX
#endif
