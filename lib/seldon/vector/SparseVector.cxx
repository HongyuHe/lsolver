// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_SPARSE_VECTOR_CXX

#include "SparseVector.hxx"

namespace Seldon
{

  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Clears the vector.
  /*!
    Destructs the vector.
    \warning On exit, the vector is an empty vector.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Clear()
  {
    // 'data_' is released.
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	if (this->data_ != NULL)
	  {
	    Allocator::deallocate(this->data_, this->m_);
	    this->data_ = NULL;
	  }

	if (index_ != NULL)
	  {
	    AllocatorInt::deallocate(index_, this->m_);
	    index_ = NULL;
	  }

	this->m_ = 0;

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
	index_ = NULL;
	this->m_ = 0;
	return;
      }
#endif

  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, previous non-zero entries may be
    lost.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::ReallocateVector(size_t i)
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

	    index_
	      = reinterpret_cast<int*>(AllocatorInt::
				       reallocate(index_, i, this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->data_ = NULL;
	    this->index_ = NULL;
	    return;
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    this->index_ = NULL;
	    return;
	  }
#endif

      }
  }


  //! Changes the number of non-zero entries of the vector.
  /*! Changes the number of non-zero entries to \a n. If \a n non-zero entries
    are available before resizing, they are all kept. Otherwise, only the
    first \n non-zero entries are kept.
    \param n new number of non-zero entries of the vector.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Resize(size_t n)
  {
    ResizeVector(n);
  }


  //! Changes the number of non-zero entries of the vector.
  /*! Changes the number of non-zero entries to \a n. If \a n non-zero entries
    are available before resizing, they are all kept. Otherwise, only the
    first \n non-zero entries are kept.
    \param n new number of non-zero entries of the vector.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::ResizeVector(size_t n)
  {
    // function implemented in the aim that explicit specialization
    // of Resize can call ResizeVector
    if (n == this->m_)
      return;

    Vector<T, VectFull, Allocator> new_value(n);
    Vector<int> new_index(n);
    size_t Nmin = min(this->m_, n);
    for (size_t i = 0; i < Nmin; i++)
      {
	new_value(i) = this->data_[i];
	new_index(i) = index_[i];
      }

    SetData(new_value, new_index);
  }


  /*! \brief Changes the length of the vector and sets its data array (low
    level method). */
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param i new length of the vector.
    \param data the new data array. \a data contains the new elements of the
    vector and must therefore contain \a i elements.
    \param index the new index array. \a index contains the new indices of the
    non-zero entries and it must therefore contain \a i elements.
    \warning \a data has to be used carefully outside the object.  Unless you
    use 'Nullify', \a data will be freed by the destructor, which means that
    \a data must have been allocated carefully. The vector allocator should be
    compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>
  ::SetData(size_t i, T* data, int* index)
  {
    this->Clear();

    this->m_ = i;

    this->data_ = data;
    this->index_ = index;
  }


  /*! \brief Changes the length of the vector and sets its data array (low
    level method). */
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param data the new data array. \a data contains the new elements of the
    vector and must therefore contain \a i elements.
    \param index the new index array. \a index contains the new indices of the
    non-zero entries and it must therefore contain \a i elements.
    \note Vectors \a data and \a index are empty vector on exit.
  */
  template <class T, class Allocator>
  template<class Allocator2>
  void Vector<T, VectSparse, Allocator>
  ::SetData(Vector<T, VectFull, Allocator2>& data, Vector<int>& index)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (data.GetM() != index.GetM())
      throw WrongDim("Vector<VectSparse>::SetData ",
		     string("The data vector and the index vector should")
		     + " have the same size.\n  Size of the data vector: "
		     + to_str(data.GetM()) + "\n  Size of index vector: "
		     + to_str(index.GetM()));
#endif

    SetData(data.GetM(), data.GetData(), index.GetData());
    data.Nullify();
    index.Nullify();
  }


  /*! \brief Lets the current vector point to the data of a second vector (low
    level method). */
  /*! Reallocates the current vector and lets its data point to those of \a V.
    \param V the vector to which the current vector points to (on exit).
    \warning On exit, both \a V and the current vector point to the same
    arrays in memory. Only one of them should eventually deallocate the memory
    blocks. The other one should be nullified by the user. In case the current
    vector is responsible for the deallocations, its allocator should be
    compatible with the allocator that created the memory blocks (which is
    probably the allocator of \a V).
  */
  template <class T, class Allocator>
  template<class Allocator2>
  void Vector<T, VectSparse, Allocator>
  ::SetData(const Vector<T, VectSparse, Allocator2>& V)
  {
    SetData(V.GetM(), V.GetData(), V.GetIndex());
  }


  //! Clears the vector without releasing memory.
  /*!
    On exit, the vector is empty and the memory has not been released.
    It is useful for low level manipulations on a Vector instance.
    \warning Memory is not released.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->data_ = NULL;
    this->index_ = NULL;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at \a i.
  */
  template <class T, class Allocator>
  typename Vector<T, VectSparse, Allocator>::value_type
  Vector<T, VectSparse, Allocator>::operator() (int i) const
  {
    size_t k = 0;
    T zero;
    SetComplexZero(zero);
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      // The entry does not exist, a zero is returned.
      return zero;

    return this->data_[k];
  }


  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the vector).
  */
  template <class T, class Allocator>
  typename Vector<T, VectSparse, Allocator>::value_type
  Vector<T, VectSparse, Allocator>::operator() (int i)
  {
    size_t k = 0;
    T zero;
    SetComplexZero(zero);
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      // The entry does not exist, a zero is returned.
      return zero;
    
    return this->data_[k];
  }
  
  
  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the vector).
  */
  template <class T, class Allocator>
  typename Vector<T, VectSparse, Allocator>::reference
  Vector<T, VectSparse, Allocator>::Get(int i)
  {
    size_t k = 0;
    T zero;
    SetComplexZero(zero);
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      // The entry does not exist yet, so a zero entry is introduced.
      AddInteraction(i, zero);
    
    return this->data_[k];
  }


  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the vector).
  */
  template <class T, class Allocator>
  typename Vector<T, VectSparse, Allocator>::const_reference
  Vector<T, VectSparse, Allocator>::Get(int i) const
  {
    size_t k = 0;
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      throw WrongArgument("Vector<VectSparse>::Val(int)",
                          "No reference to element " + to_str(i)
                          + " can be returned: it is a zero entry.");
    
    return this->data_[k];
  }


  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument if i does not belong to the sparsity pattern
  */
  template <class T, class Allocator>
  typename Vector<T, VectSparse, Allocator>::reference
  Vector<T, VectSparse, Allocator>::Val(int i)
  {
    size_t k = 0;
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      throw WrongArgument("Vector<VectSparse>::Val(int)",
                          "the entry " + to_str(i) +
                          " does not belong to the sparsity pattern.");

    
    return this->data_[k];
  }


  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument if i does not belong to the sparsity pattern
  */
  template <class T, class Allocator>
  typename Vector<T, VectSparse, Allocator>::const_reference
  Vector<T, VectSparse, Allocator>::Val(int i) const
  {
    size_t k = 0;
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      throw WrongArgument("Vector<VectSparse>::Val(int)",
                          "the entry " + to_str(i) +
                          " does not belong to the sparsity pattern.");
    
    return this->data_[k];
  }

  
  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: \a X is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator> template <class T2, class Alloc2>
  void Vector<T, VectSparse, Allocator>
  ::Copy(const Vector<T2, VectSparse, Alloc2>& X)
  {
    this->Reallocate(X.GetLength());
    AllocatorInt::memorycpy(this->index_, X.GetIndex(), this->m_);
    for (size_t i = 0; i < X.GetLength(); i++)
      this->Value(i) = X.Value(i);
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  Vector<T, VectSparse, Allocator>&
  Vector<T, VectSparse, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Displays the vector.
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Print() const
  {
    for (int i = 0; i < this->GetM(); i++)
      cout << (Index(i) + 1) << ' ' << Value(i) << '\n';
  }


  //! Assembles the vector.
  /*!
    \warning If you use the method AddInteraction, you don't need to call
    that method.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Assemble()
  {
    int new_size = this->m_;
    Vector<T, VectFull, Allocator> values(new_size);
    Vector<int> index(new_size);
    for (int i = 0; i < new_size; i++)
      {
	values(i) = this->data_[i];
	index(i) = index_[i];
      }

    Seldon::Assemble(new_size, index, values);
    index.Resize(new_size);
    values.Resize(new_size);
    SetData(values, index);
  }


  //! Removes small entries.
  /*! Any number whose absolute value is below (or equal) to \a epsilon is
    removed.
    \param epsilon the threshold value.
  */
  template <class T, class Allocator> template<class T0>
  void Vector<T, VectSparse, Allocator>::RemoveSmallEntry(const T0& epsilon)
  {
    size_t new_size = this->m_;
    Vector<T, VectFull, Allocator> values(new_size);
    Vector<int> index(new_size);
    new_size = 0;
    for (size_t i = 0; i < this->m_; i++)
      if (abs(this->data_[i]) > epsilon)
	{
	  values(new_size) = this->data_[i];
	  index(new_size) = index_[i];
	  new_size++;
	}

    index.Resize(new_size);
    values.Resize(new_size);
    SetData(values, index);
  }


  //! Adds \a val to the vector component #\a i.
  /*! If the vector has no entry at \a i, a new entry with value \a val is
    introduced. Otherwise, this method sums the existing value and \a val.
    \param[in] i index of the component.
    \param[in] val value to be added to the vector component \a i.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::AddInteraction(int i, const T& val)
  {
    // Searching for the position where the entry may be.
    size_t pos = 0;
    while (pos < this->m_ && index_[pos] < i)
      pos++;

    // If the entry already exists, adds 'val'.
    if (pos < this->m_ && index_[pos] == i)
      {
	this->data_[pos] += val;
	return;
      }

    size_t k;

    // If the entry does not exist, the vector is reallocated.
    Resize(this->m_ + 1);
    
    if (this->m_ >= 2)
      for (k = this->m_-1; k > pos; k--)
        {
          this->data_[k] = this->data_[k-1];
          this->index_[k] = this->index_[k-1];
        }
    
    // The new entry.
    this->index_[pos] = i;
    this->data_[pos] = val;
  }


  //! Adds given values to several components of the vector.
  /*! This method sorts the values to be added (according to their indices)
    and adds them with the vector values. For every component, if the vector
    has no entry, a new entry is introduced. Otherwise, the method sums the
    existing value and the corresponsing value in \a value.
    \param[in] n number of values to be added.
    \param[in] index indices of the values to be added.
    \param[in] value values to be added.
    \param[in] already_sorted true if the indices are already sorted.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::
  AddInteractionRow(size_t n, int* index, T* value, bool already_sorted)
  {
    Vector<int> ind;
    Vector<T, VectFull, Allocator> val;
    ind.SetData(n, index);
    val.SetData(n, value);
    AddInteractionRow(n, ind, val, already_sorted);
    ind.Nullify();
    val.Nullify();
  }


  //! Adds given values to several components of the vector.
  /*! This method sorts the values to be added (according to their indices)
    and adds them with the vector values. For every component, if the vector
    has no entry, a new entry is introduced. Otherwise, the method sums the
    existing value and the corresponsing value in \a value.
    \param[in] n number of values to be added.
    \param[in] index indices of the values to be added.
    \param[in] value values to be added.
    \param[in] already_sorted true if the indices are already sorted.
  */
  template <class T, class Allocator>
  template<class Allocator0>
  void Vector<T, VectSparse, Allocator>::
  AddInteractionRow(size_t n, const Vector<int>& index2,
		    const Vector<T, VectFull, Allocator0>& value2,
		    bool already_sorted)
  {
    Vector<int> index;
    Vector<T, VectFull, Allocator0> value;

    if (!already_sorted)
      {
        // checking if numbers are sorted
        already_sorted = true;
        if (n >= 2)
          for (size_t i = 0; i < n-1; i++)
            if (index2(i+1) <= index2(i))
              already_sorted = false;
      }
    
    if (!already_sorted)
      {
        index.Reallocate(n);
        value.Reallocate(n);
        for (size_t i = 0; i < n; i++)
          {
            index(i) = index2(i);
            value(i) = value2(i);
          }
        
        // Sorts the values to be added according to their indices.
        int ni = n; Seldon::Assemble(ni, index, value); n = ni;
      }
    else
      {
        index.SetData(n, index2.GetData());
        value.SetData(n, value2.GetData());
      }
    
    //  Values that already have an entry 
    
    // Number of values to be added without entry.
    size_t Nnew = 0;
    Vector<bool> new_index(n);
    new_index.Fill(true);
    size_t k = 0;
    for (size_t j = 0; j < n; j++)
      {
	while (k < this->m_ && index_[k] < index(j))
	  k++;

	if (k < this->m_ && index(j) == index_[k])
	  {
	    new_index(j) = false;
	    this->data_[k] += value(j);
	  }
	else
	  Nnew++;
      }

    if (Nnew > 0)
      {
	// Some values to be added have no entry yet.
	Vector<T> new_val(this->m_ + Nnew);
	Vector<int> new_ind(this->m_ + Nnew);
	int nb = 0;
	k = 0;
	for (size_t j = 0; j < n; j++)
	  if (new_index(j))
	    {
	      while (k < this->m_ && index_[k] < index(j))
		{
		  new_ind(nb) = index_[k];
		  new_val(nb) = this->data_[k];
		  k++;
		  nb++;
		}

	      // The new entry.
	      new_ind(nb) = index(j);
	      new_val(nb) = value(j);
	      nb++;
	    }

	// Last entries.
	while (k < this->m_)
	  {
	    new_ind(nb) = index_[k];
	    new_val(nb) = this->data_[k];
	    k++;
	    nb++;
	  }

	SetData(new_val, new_ind);
      }

    if (already_sorted)
      {
        index.Nullify();
        value.Nullify();
      }
  }

  
  //! Returns the infinite norm.
  /*!
    \return The infinite norm.
  */
  template <class T, class Allocator>
  typename ClassComplexType<T>::Treal
  Vector<T, VectSparse, Allocator>::GetNormInf() const
  {
    typename ClassComplexType<T>::Treal res(0);
    for (size_t i = 0; i < this->m_; i++)
      res = max(res, abs(this->data_[i]));
    
    return res;
  }


  //! Returns the index of the highest absolute value.
  /*!
    \return The index of the element that has the highest absolute value.
  */
  template <class T, class Allocator>
  int Vector<T, VectSparse, Allocator>::GetNormInfIndex() const
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (this->GetLength() == 0)
      throw WrongDim("Vector<VectSparse>::GetNormInfIndex()",
		     "Vector is null.");
#endif

    typename ClassComplexType<T>::Treal res(0), temp;
    size_t j = 0;
    for (size_t i = 0; i < this->GetLength(); i++)
      {
	temp = res;
	res = max(res, abs(this->data_[i]));
	if (temp != res) j = i;
      }

    return this->index_[j];
  }

  
  /**************************
   * OUTPUT/INPUT FUNCTIONS *
   **************************/


  //! Writes the vector in a file.
  /*! It stores in binary format: (1) the number of non-zero entries in the
    vector (integer), (2) the indices of the non-zero entries (integers), and
    (3) the non-zero values of the vector.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectSparse>::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();
  }


  //! Writes the vector in a stream, in binary format.
  /*! It writes in binary format: (1) the number of non-zero entries in the
    vector (integer), (2) the indices of the non-zero entries (integers), and
    (3) the non-zero values of the vector.
    \param stream stream in which the vector is to be written.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Write(ostream& stream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!stream.good())
      throw IOError("Vector<VectSparse>::Write(ostream& stream)",
                    "Stream is not ready.");
#endif

    int n = this->m_;
    stream.write(reinterpret_cast<char*>(&n),
		 sizeof(int));

    stream.write(reinterpret_cast<char*>(this->index_),
		 this->m_ * sizeof(int));

    stream.write(reinterpret_cast<char*>(this->data_),
		 this->m_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!stream.good())
      throw IOError("Vector<VectSparse>::Write(ostream& stream)",
                    "Output operation failed.");
#endif

  }


  //! Writes the vector in a text file.
  /*! All non-zero elements of the vector are stored in text format: every
    line of the text file contains one index and one value. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectSparse>::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the vector in a stream, in text format.
  /*! All non-zero elements of the vector are stored in text format: every
    line of the text file contains one index and one value. The length is not
    stored.
    \param stream stream in which the vector is to be written.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::WriteText(ostream& stream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!stream.good())
      throw IOError("Vector<VectSparse>::WriteText(ostream& stream)",
                    "Stream is not ready.");
#endif

    // First entries.
    for (size_t i = 0; i < this->m_; i++)
      stream << (Index(i) + 1) << " " << Value(i) << '\n';
    
#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!stream.good())
      throw IOError("Vector<VectSparse>::WriteText(ostream& stream)",
                    "Output operation failed.");
#endif

  }


  //! Sets the vector from a file in binary format.
  /*! Sets the vector according to a binary file that stores the data like
    method Write(string).
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectSparse>::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Sets the vector from a file stream, in binary format.
  /*! Sets the vector according to a binary stream that stores the data like
    method Write(ostream&).
    \param stream stream from which to read the vector values.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::Read(istream& stream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!stream.good())
      throw IOError("Vector<VectSparse>::Read(istream& stream)",
                    "Stream is not ready.");
#endif

    int m;
    stream.read(reinterpret_cast<char*>(&m), sizeof(int));
    this->Reallocate(m);

    stream.read(reinterpret_cast<char*>(this->index_), m * sizeof(int));

    stream.read(reinterpret_cast<char*>(this->data_),
		m * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!stream.good())
      throw IOError("Vector<VectSparse>::Read(istream& stream)",
                    "Input operation failed.");
#endif

  }


  //! Sets the vector from a text file.
  /*! Sets the vector according to a text file that stores the data like
    method WriteText(string).
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectSparse>::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream);

    FileStream.close();
  }


  //! Sets the vector from a file stream, in text format.
  /*! Sets the vector according to a stream, in text format, that stores the
    data like method WriteText(ostream&).
    \param stream stream from which to read the vector values.
  */
  template <class T, class Allocator>
  void Vector<T, VectSparse, Allocator>::ReadText(istream& stream)
  {
    Clear();

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!stream.good())
      throw IOError("Vector<VectSparse>::ReadText(istream& stream)",
                    "Stream is not ready.");
#endif

    Vector<T, VectFull, Allocator> values;
    Vector<int> index;
    T entry;
    int ind = 0;
    int nb_elt = 0;
    while (!stream.eof())
      {
	// New entry is read.
	stream >> ind >> entry;

	if (stream.fail())
	  break;
	else
	  {
#ifdef SELDON_CHECK_IO
	    if (ind < 1)
	      throw IOError(string("Vector<VectSparse>::ReadText") +
			    "(istream& stream)",
			    string("Index should be greater ")
			    + "than 0 but is equal to " + to_str(ind) + ".");
#endif

	    nb_elt++;
	    ind--;

	    // Inserting a new element.
	    if (nb_elt > values.GetM())
	      {
		values.Resize(2 * nb_elt);
		index.Resize(2 * nb_elt);
	      }

	    values(nb_elt - 1) = entry;
	    index(nb_elt - 1) = ind;
	  }
      }

    if (nb_elt > 0)
      {
	// Allocating to the right size.
	this->Reallocate(nb_elt);
	for (int i = 0; i < nb_elt; i++)
	  {
	    Index(i) = index(i);
	    Value(i) = values(i);
	  }
      }
  }


  //! operator<< overloaded for sparse vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, VectSparse, Allocator>& V)
  {
    for (size_t i = 0; i < V.GetLength(); i++)
      out << (V.Index(i) + 1) << ' ' << V.Value(i) << '\n';

    return out;
  }
  
} // namespace Seldon.

#define SELDON_FILE_SPARSE_VECTOR_CXX
#endif
