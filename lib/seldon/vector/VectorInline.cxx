// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2003-2015 Marc Duruflé
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


#ifndef SELDON_FILE_VECTOR_INLINE_CXX

#include "Vector.hxx"
#include "VectorExpressionInline.cxx"

namespace Seldon
{


  ///////////////////
  // C_VECTOR_BASE //
  ///////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Nothing is allocated. Vector length is set to zero.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::Vector_Base()
  {
    m_ = 0;
    data_ = NULL;
  }


  //! Main constructor.
  /*!
    \param i length.
    \warning Nothing is allocated.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::Vector_Base(size_t i)
  {
    m_ = i;
    data_ = NULL;
  }


  //! Copy constructor.
  /*!
    \param A base vector to be copied.
    \warning Only the length is copied.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::
  Vector_Base(const Vector_Base<T, Allocator>& A)
  {
    m_ = A.GetM();
    data_ = NULL;
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::~Vector_Base()
  {
    this->Clear();
  }


  //! Releases memory used by the vector.
  template <class T, class Allocator>
  inline void Vector_Base<T, Allocator>::Clear()
  {
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (data_ != NULL)
	  {
	    Allocator::deallocate(data_, m_);
	    m_ = 0;
	    data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	m_ = 0;
	data_ = NULL;
      }
#endif

  }

  
  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements.
  /*!
    \return The length of the vector.
  */
  template <class T, class Allocator>
  inline long Vector_Base<T, Allocator>::GetM() const
  {
    return m_;
  }


  //! Returns the number of elements.
  /*!
    \return The length of the vector.
  */
  template <class T, class Allocator>
  inline size_t Vector_Base<T, Allocator>::GetLength() const
  {
    return m_;
  }


  //! Returns the number of elements stored.
  /*!
    \return The length of the vector stored.
  */
  template <class T, class Allocator>
  inline size_t Vector_Base<T, Allocator>::GetSize() const
  {
    return m_;
  }


  //! Returns the memory used by the object in bytes.
  /*!
    In this method, the type T is assumed to be "static"
    such that sizeof(T) provides the correct size
  */
  template <class T, class Allocator>
  inline size_t Vector_Base<T, Allocator>::GetMemorySize() const
  {
    return sizeof(*this) + sizeof(T)*size_t(GetSize());
  }

  
  //! Returns a pointer to data_ (stored data).
  /*!
    \return A pointer to the data_, i.e. the data array.
  */
  template <class T, class Allocator>
  inline typename Vector_Base<T, Allocator>::pointer
  Vector_Base<T, Allocator>::GetData() const
  {
    return data_;
  }


  //! Returns a const pointer to data_ (stored data).
  /*!
    \return A const pointer to the data_, i.e. the data array.
  */
  template <class T, class Allocator>
  inline typename Vector_Base<T, Allocator>::const_pointer
  Vector_Base<T, Allocator>::GetDataConst() const
  {
    return reinterpret_cast<typename Vector_Base<T,
      Allocator>::const_pointer>(data_);
  }


  //! Returns a pointer of type "void*" to the data array (data_).
  /*!
    \return A pointer of type "void*" to the data array.
  */
  template <class T, class Allocator>
  inline void* Vector_Base<T, Allocator>::GetDataVoid() const
  {
    return reinterpret_cast<void*>(data_);
  }


  //! Returns a pointer of type "const void*" to the data array (data_).
  /*!
    \return A pointer of type "const void*" to the data array.
  */
  template <class T, class Allocator>
  inline const void* Vector_Base<T, Allocator>::GetDataConstVoid() const
  {
    return reinterpret_cast<const void*>(data_);
  }


  ///////////////////////
  // VECTOR<VECT_FULL> //
  ///////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>::Vector():
    Vector_Base<T, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param i length of the vector.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>::Vector(size_t i):
    Vector_Base<T, Allocator>(i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = Allocator::allocate(i, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      this->m_ = 0;
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Vector<VectFull>::Vector(int)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(i*sizeof(T)) + " bytes ("
		     + to_str(i) + " elements).");
#endif

  }


  //! Builds a vector using pre-existing data.
  /*!
    \param i length of the vector.
    \param data the data array. \a data contains the elements of the vector
    and must therefore contain \a i elements.
    \warning \a data has to be used carefully outside the object. Unless you
    use 'Nullify', \a data will be freed by the destructor, which means that
    \a data must have been allocated carefully. The vector allocator should be
    compatible.
    \note This constructor should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>
  ::Vector(size_t i, typename Vector<T, VectFull, Allocator>::pointer data):
    Vector_Base<T, Allocator>()
  {
    SetData(i, data);
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param V vector to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>::
  Vector(const Vector<T, VectFull, Allocator>& V):
    Vector_Base<T, Allocator>(V)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = Allocator::allocate(V.GetM(), this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      this->m_ = 0;
    if (this->data_ == NULL && V.GetM() != 0)
      throw NoMemory("Vector<VectFull>::Vector(Vector<VectFull>&)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(V.GetM()*sizeof(T)) + " bytes ("
		     + to_str(V.GetM()) + " elements).");
#endif

    Allocator::memorycpy(this->data_, V.GetData(), V.GetM());

  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>::~Vector()
  {
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
  inline void Vector<T, VectFull, Allocator>::Clear()
  {
    Vector_Base<T, Allocator>::Clear();
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, initial elements of the vector may
    be lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::Reallocate(size_t i)
  {
    ReallocateVector(i);
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, initial elements of the vector may
    be lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::ReallocateVector(size_t i)
  {
    // function implemented in the aim that explicit specialization
    // of Reallocate can call ReallocateVector
    if (i != this->m_)
      {

	this->m_ = i;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(Allocator::
					reallocate(this->data_, i, this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->data_ = NULL;
	    return;
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    return;
	  }
#endif

      }
  }


  /*! \brief Changes the length of the vector and sets its data array (low
    level method). */
  /*!
    Deallocates memory allocated for the current data array, changes the
    vector size and sets the new data array. It is useful to create a vector
    from pre-existing data.
    \param i new length of the vector.
    \param data the new data array. 'data' contains the new elements of the
    vector and must therefore contain 'i' elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The vector
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>
  ::SetData(size_t i, typename Vector<T, VectFull, Allocator>::pointer data)
  {
    this->Clear();

    this->m_ = i;

    this->data_ = data;
  }


  //! Lets the current vector point to the data of another vector.
  /*! Deallocates memory allocated for the current data array. Then sets the
    length and the data of the current vector to that of \a V. On exit, the
    current vector shares it data array with \a V.
    \param V vector whose data should be shared with current instance.
    \warning On exit, \a V and the current instance share the same data array,
    and they are likely to free it. As a consequence, before the vectors are
    destructed, it is necessary to call 'Nullify' on either \a V or the
    current instance. In addition, if the current instance is to deallocate
    the data array, its allocator should be compatible with the allocator that
    allocated the data array (probably the allocator of \a V).
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  template <class Allocator0>
  inline void Vector<T, VectFull, Allocator>
  ::SetData(const Vector<T, VectFull, Allocator0>& V)
  {
    SetData(V.GetLength(), V.GetData());
  }


  //! Clears the vector without releasing memory.
  /*!
    On exit, the vector is empty and the memory has not been released.
    It is useful for low level manipulations on a Vector instance.
    \warning Memory is not released.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->data_ = NULL;
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
  inline typename Vector<T, VectFull, Allocator>::reference
  Vector<T, VectFull, Allocator>::operator() (long i)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "Vector<VectFull>");
#endif

    return this->data_[i];
  }


  //! Access to element \a i.
  /*!
    \param i index.
    \return The value of the vector at \a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectFull, Allocator>::reference
  Vector<T, VectFull, Allocator>::Get(long i)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "Vector<VectFull>");
#endif

    return this->data_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectFull, Allocator>::const_reference
  Vector<T, VectFull, Allocator>::operator() (long i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "Vector<VectFull>");
#endif

    return this->data_[i];
  }


  //! Access to element \a i.
  /*!
    \param i index.
    \return The value of the vector at \a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectFull, Allocator>::const_reference
  Vector<T, VectFull, Allocator>::Get(long i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "Vector<VectFull>");
#endif

    return this->data_[i];
  }


  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>& Vector<T, VectFull, Allocator>
  ::operator= (const Vector<T, VectFull, Allocator>& X)
  {
    this->Copy(X);

    return *this;
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>
  ::Copy(const Vector<T, VectFull, Allocator>& X)
  {
    this->Reallocate(X.GetLength());

    Allocator::memorycpy(this->data_, X.GetData(), this->m_);
  }


  //! Duplicates a vector.
  /*!
    \return A copy of the vector.
    \note Memory is duplicated: the returned vector is therefore independent
    from the current instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>
  Vector<T, VectFull, Allocator>::Copy() const
  {
    return Vector<T, VectFull, Allocator>(*this);
  }


  //! Multiplies a vector by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Allocator> template<class T0>
  inline Vector<T, VectFull, Allocator>& Vector<T, VectFull, Allocator>
  ::operator*= (const T0& alpha)
  {
    for (size_t i = 0; i < this->m_; i++)
      this->data_[i] *= alpha;

    return *this;
  }


  //! Operator =, *this = X where X is an expression
  template<class T, class Allocator> template<class E>
  inline Vector<T, VectFull, Allocator>& Vector<T, VectFull, Allocator>
  ::operator =(const VectorExpression<T, E>& X)
  {
    this->Reallocate(X.GetSize());
    for (size_t i = 0; i < this->m_; i++)
      this->data_[i] = X(i);
    
    return *this;
  }


  //! Operator +=, *this = *this + X where X is an expression
  template<class T, class Allocator> template<class E>
  inline Vector<T, VectFull, Allocator>& Vector<T, VectFull, Allocator>
  ::operator +=(const VectorExpression<T, E>& X)
  {
    for (size_t i = 0; i < this->m_; i++)
      this->data_[i] += X(i);
    
    return *this;
  }


  //! Operator -=, *this = *this - X where X is an expression
  template<class T, class Allocator> template<class E>
  inline Vector<T, VectFull, Allocator>& Vector<T, VectFull, Allocator>
  ::operator -=(const VectorExpression<T, E>& X)
  {
    for (size_t i = 0; i < this->m_; i++)
      this->data_[i] -= X(i);
    
    return *this;
  }


  //! Appends an element to the vector.
  /*!
    \param x element to be appended.
    \warning This method will only work if the allocator preserves the
    elements while reallocating.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::Append(const T& x)
  {
    size_t i = this->GetLength();
    this->Reallocate(i + 1);
    this->data_[i] = x;
  }


  //! Appends an element at the end of the vector.
  /*!
    \param x element to be appended.
  */
  template <class T, class Allocator> template<class T0>
  inline void Vector<T, VectFull, Allocator>::PushBack(const T0& x)
  {
    Resize(this->m_+1);
    this->data_[this->m_-1] = x;
  }


  //! Appends a vector X at the end of the vector.
  /*!
    \param X vector to be appended.
  */
  template <class T, class Allocator> template<class Allocator0>
  inline void Vector<T, VectFull, Allocator>
  ::PushBack(const Vector<T, VectFull, Allocator0>& X)
  {
    size_t Nold = this->m_;
    Resize(this->m_ + X.GetSize());
    for (size_t i = 0; i < X.GetSize(); i++)
      this->data_[Nold+i] = X(i);
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Allocator>
  inline size_t Vector<T, VectFull, Allocator>::GetDataSize()
  {
    return this->m_;
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
  inline void Vector<T, VectFull, Allocator>::Zero()
  {
    Allocator::memoryset(this->data_, char(0),
			 this->GetDataSize() * sizeof(value_type));
  }


  //! Fills the vector with 0, 1, 2, ...
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::Fill()
  {
    for (size_t i = 0; i < this->m_; i++)
      SetComplexReal(i, this->data_[i]);
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  inline void Vector<T, VectFull, Allocator>::Fill(const T0& x)
  {
    T x_;
    SetComplexReal(x, x_);
    for (size_t i = 0; i < this->m_; i++)
      this->data_[i] = x_;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>&
  Vector<T, VectFull, Allocator>::operator= (const T& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Fills the vector randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::FillRand()
  {
#ifndef SELDON_WITHOUT_REINIT_RANDOM
    srand(time(NULL));
#endif
    for (size_t i = 0; i < this->m_; i++)
      SetComplexReal(rand(), this->data_[i]);
  }


  //! Displays the vector.
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::Print() const
  {
    for (size_t i = 0; i < this->GetLength(); i++)
      cout << (*this)(i) << "\t";
    cout << endl;
  }
  
  
} // namespace Seldon.

#define SELDON_FILE_VECTOR_INLINE_CXX
#endif
