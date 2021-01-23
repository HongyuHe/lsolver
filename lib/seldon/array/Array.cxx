// Copyright (C) 2010 Lin Wu
// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_ARRAY_CXX


namespace Seldon
{

  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the array is an empty array.
  */
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>::Array()
  {
    if (N < ARRAY_MINRANK || N > ARRAY_MAXRANK)
      {
	string msg = string("Array dimension should be in [") +
	  to_str(ARRAY_MINRANK) + string(", ") + to_str(ARRAY_MAXRANK) + "].";
	throw WrongDim("Array<T, N, Allocator>::Array(int, ...)", msg);
      }

    data_ = NULL;
    for (int k = 0; k < N; k++)
      {
        length_[k] = 0;
        offset_[k] = 0;
      }
  }


  //! Main constructor.
  /*! Builds a i x j x k 3D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
  */
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>::Array(int i, int j, int k)
  {
    if (N != 3)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 3.");

    long nb_elt = long(i) * long(j) * long(k);

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	data_ = Allocator::allocate(nb_elt, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    offset_[0] = long(j) * long(k);
    offset_[1] = k;
    offset_[2] = nb_elt;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l 4D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
  */
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>::Array(int i, int j, int k, int l)
  {
    if (N != 4)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 4.");

    long nb_elt = long(i) * long(j) * long(k) * long(l);

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	data_ = Allocator::allocate(nb_elt, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    offset_[0] = long(j) * long(k) * long(l);
    offset_[1] = long(k) * long(l);
    offset_[2] = long(l);
    offset_[3] = nb_elt;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m 5D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
  */
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>::Array(int i, int j, int k, int l, int m)
  {
    if (N != 5)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 5.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m);

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	data_ = Allocator::allocate(nb_elt, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x "
		     + to_str(m) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    offset_[0] = long(j) * long(k) * long(l) * long(m);
    offset_[1] = long(k) * long(l) * long(m);
    offset_[2] = long(l) * long(m);
    offset_[3] = m;
    offset_[4] = nb_elt;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m x n 6D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
  */
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>
  ::Array(int i, int j, int k, int l, int m, int n)
  {
    if (N != 6)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 6.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m) * long(n);

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	data_ = Allocator::allocate(nb_elt, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }

    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	n != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(n)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x " + 
		     to_str(m) + " x " + to_str(n) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    length_[5] = n;

    offset_[4] = n;
    offset_[3] = long(m) * long(n);
    offset_[2] = long(l) * offset_[3];
    offset_[1] = long(k) * offset_[2];
    offset_[0] = long(j) * offset_[1];
    offset_[5] = nb_elt;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m x n x o 7D array, but data is not
    initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
  */
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>
  ::Array(int i, int j, int k, int l, int m, int n, int o)
  {
    if (N != 7)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 7.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m)
      * long(n) * long(o);
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	data_ = Allocator::allocate(nb_elt, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	n != 0 && o != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(n)
			      * static_cast<long int>(o)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x " + 
		     to_str(m) + " x " + to_str(n) + " x " + to_str(o) +
		     " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    length_[5] = n;
    length_[6] = o;

    offset_[5] = o;
    offset_[4] = long(n) * long(o);
    offset_[3] = long(m) * offset_[4];
    offset_[2] = long(l) * offset_[3];
    offset_[1] = long(k) * offset_[2];
    offset_[0] = long(j) * offset_[1];
    offset_[6] = nb_elt;
  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m x n x o x p 8D array, but data is not
    initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \param p length in dimension #8.
  */
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>
  ::Array(int i, int j, int k, int l, int m, int n, int o, int p)
  {
    if (N != 8)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 8.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m)
      * long(n) * long(o) * long(p);

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	data_ = Allocator::allocate(nb_elt, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	n != 0 && o != 0 && p != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(n)
			      * static_cast<long int>(o)
			      * static_cast<long int>(p)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x " + 
		     to_str(m) + " x " + to_str(n) + " x " + to_str(o)
		     + " x " + to_str(p) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    length_[5] = n;
    length_[6] = o;
    length_[7] = p;

    offset_[6] = p;
    offset_[5] = long(o) * long(p);
    offset_[4] = long(n) * offset_[5];
    offset_[3] = long(m) * offset_[4];
    offset_[2] = long(l) * offset_[3];
    offset_[1] = long(k) * offset_[2];
    offset_[0] = long(j) * offset_[1];
    offset_[7] = nb_elt;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m x n x o x p x q 9D array, but data is not
    initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \param p length in dimension #8.
    \param q length in dimension #9.
  */
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>
  ::Array(int i, int j, int k, int l, int m, int n, int o, int p, int q)
  {
    if (N != 9)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 9.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m)
      * long(n) * long(o) * long(p) * long(q);

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	data_ = Allocator::allocate(nb_elt, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	n != 0 && o != 0 && p != 0 && q != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(n)
			      * static_cast<long int>(o)
			      * static_cast<long int>(p)
			      * static_cast<long int>(q)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x " + 
		     to_str(m) + " x " + to_str(n) + " x " + to_str(o) + " x "
		     + to_str(p) + " x " + to_str(q) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    length_[5] = n;
    length_[6] = o;
    length_[7] = p;
    length_[8] = q;

    offset_[7] = q;
    offset_[6] = long(p) * long(q);
    offset_[5] = long(o) * offset_[6];
    offset_[4] = long(n) * offset_[5];
    offset_[3] = long(m) * offset_[4];
    offset_[2] = long(l) * offset_[3];
    offset_[1] = long(k) * offset_[2];
    offset_[0] = long(j) * offset_[1];
    offset_[8] = nb_elt;

  }


  //! Copy constructor.
  template <class T, int N, class Allocator>
  Array<T, N, Allocator>::Array(const Array<T, N, Allocator>& A)
  {
    data_ = NULL;

    Copy(A);
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the 3D array.
  /*!
    On exit, the array is a i x j x k 3D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Reallocate(int i, int j, int k)
  {
    if (N != 3)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 3.");

    long nb_elt = long(i) * long(j) * long(k);

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    data_ =
	      reinterpret_cast<pointer>(Allocator::reallocate(data_,
							      nb_elt,
							      this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " elements).");
#endif

        length_[0] = i;
        length_[1] = j;
        length_[2] = k;
        offset_[0] = long(j) * long(k);
        offset_[1] = k;
        offset_[2] = nb_elt;
      }
  }


  //! Reallocates memory to resize the 4D array.
  /*!
    On exit, the array is a i x j x k x l 4D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Reallocate(int i, int j, int k, int l)
  {
    if (N != 4)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 4.");

    long nb_elt = long(i) * long(j) * long(k) * long(l);

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    data_ =
	      reinterpret_cast<pointer>(Allocator::reallocate(data_,
							      nb_elt,
							      this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " elements).");
#endif

        length_[0] = i;
        length_[1] = j;
        length_[2] = k;
        length_[3] = l;
        offset_[0] = long(j) * long(k) * long(l);
        offset_[1] = long(k) * long(l);
        offset_[2] = long(l);
        offset_[3] = nb_elt;
      }
  }


  //! Reallocates memory to resize the 5D array.
  /*!
    On exit, the array is a i x j x k x l x m 5D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m)
  {
    if (N != 5)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 5.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m);

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    data_ =
	      reinterpret_cast<pointer>(Allocator::reallocate(data_,
							      nb_elt,
							      this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " elements).");
#endif

        length_[0] = i;
        length_[1] = j;
        length_[2] = k;
        length_[3] = l;
        length_[4] = m;
        offset_[0] = long(j) * long(k) * long(l) * long(m);
        offset_[1] = long(k) * long(l) * long(m);
        offset_[2] = long(l) * long(m);
        offset_[3] = m;
        offset_[4] = nb_elt;
      }
  }


  //! Reallocates memory to resize the 6D array.
  /*!
    On exit, the array is a i x j x k x l x m x n 6D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m, int n)
  {
    if (N != 6)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 6.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m) * long(n);

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4]
	|| n != length_[5])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    data_ =
	      reinterpret_cast<pointer>(Allocator::reallocate(data_,
                                                              nb_elt,
							      this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	    n != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(n)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " x " + to_str(n)
			 + " elements).");
#endif
        
        length_[0] = i;
        length_[1] = j;
        length_[2] = k;
        length_[3] = l;
        length_[4] = m;
        length_[5] = n;
        
        offset_[4] = n;
        offset_[3] = long(m) * long(n);
        offset_[2] = long(l) * offset_[3];
        offset_[1] = long(k) * offset_[2];
        offset_[0] = long(j) * offset_[1];
        offset_[5] = nb_elt;
      }
  }


  //! Reallocates memory to resize the 7D array.
  /*!
    On exit, the array is a i x j x k x l x m x n x o 7D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m, int n, int o)
  {
    if (N != 7)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 7.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m)
      * long(n) * long(o);

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4]
	|| n != length_[5] || o != length_[6])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    data_ =
	      reinterpret_cast<pointer>(Allocator::reallocate(data_, nb_elt, this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	    n != 0 && o != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(n)
				  * static_cast<long int>(o)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " x " + to_str(n)
			 + " x " + to_str(o) + " elements).");
#endif

        length_[0] = i;
        length_[1] = j;
        length_[2] = k;
        length_[3] = l;
        length_[4] = m;
        length_[5] = n;
        length_[6] = o;
        
        offset_[5] = o;
        offset_[4] = long(n) * long(o);
        offset_[3] = long(m) * offset_[4];
        offset_[2] = long(l) * offset_[3];
        offset_[1] = long(k) * offset_[2];
        offset_[0] = long(j) * offset_[1];
        offset_[6] = nb_elt;
      }
  }


  //! Reallocates memory to resize the 8D array.
  /*!
    On exit, the array is a i x j x k x l x m x n x o x p 8D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \param p length in dimension #8.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m, int n, int o, int p)
  {
    if (N != 8)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 8.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m)
      * long(n) * long(o) * long(p);

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4]
	|| n != length_[5] || o != length_[6] || p != length_[7])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif
	    data_ =
	      reinterpret_cast<pointer>(Allocator::reallocate(data_, nb_elt, this));
            

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	    n != 0 && o != 0 && p != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(n)
				  * static_cast<long int>(o)
				  * static_cast<long int>(p)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " x " + to_str(n)
			 + " x " + to_str(o) + " x " + to_str(p)
			 + " elements).");
#endif

        length_[0] = i;
        length_[1] = j;
        length_[2] = k;
        length_[3] = l;
        length_[4] = m;
        length_[5] = n;
        length_[6] = o;
        length_[7] = p;
        
        offset_[6] = p;
        offset_[5] = long(o) * long(p);
        offset_[4] = long(n) * offset_[5];
        offset_[3] = long(m) * offset_[4];
        offset_[2] = long(l) * offset_[3];
        offset_[1] = long(k) * offset_[2];
        offset_[0] = long(j) * offset_[1];
        offset_[7] = nb_elt;
      }
  }


  //! Reallocates memory to resize the 9D array.
  /*!
    On exit, the array is a i x j x k x l x m x n x o x p x q 9D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \param p length in dimension #8.
    \param q length in dimension #9.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m, int n, int o, int p, int q)
  {
    if (N != 9)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 9.");

    long nb_elt = long(i) * long(j) * long(k) * long(l) * long(m)
      * long(n) * long(o) * long(p) * long(q);
    
    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4]
	|| n != length_[5] || o != length_[6] || p != length_[7]
	|| q != length_[8])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    data_ =
	      reinterpret_cast<pointer>(Allocator::reallocate(data_, nb_elt, this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	    n != 0 && o != 0 && p != 0 && q != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(n)
				  * static_cast<long int>(o)
				  * static_cast<long int>(p)
				  * static_cast<long int>(q)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " x " + to_str(n)
			 + " x " + to_str(o) + " x " + to_str(p)
			 + " x " + to_str(q) + " elements).");
#endif

        length_[0] = i;
        length_[1] = j;
        length_[2] = k;
        length_[3] = l;
        length_[4] = m;
        length_[5] = n;
        length_[6] = o;
        length_[7] = p;
        length_[8] = q;
        
        offset_[7] = q;
        offset_[6] = long(p) * long(q);
        offset_[5] = long(o) * offset_[6];
        offset_[4] = long(n) * offset_[5];
        offset_[3] = long(m) * offset_[4];
        offset_[2] = long(l) * offset_[3];
        offset_[1] = long(k) * offset_[2];
        offset_[0] = long(j) * offset_[1];
        offset_[8] = nb_elt;
      }
  }


  //! Clears the array.
  /*!
    Destructs the array.
    \warning On exit, the array is empty.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Clear()
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (data_ != NULL)
	  {
	    Allocator::deallocate(data_, offset_[N-1]);
	    data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }
#endif
    
  }


  //! Duplicates an array.
  /*!
    \param A array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Copy(const Array<T, N, Allocator>& A)
  {
    if (N == 3)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2));
    if (N == 4)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3));
    if (N == 5)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4));
    if (N == 6)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5));
    if (N == 7)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6));
    if (N == 8)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6), A.GetLength(7));
    if (N == 9)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6), A.GetLength(7), A.GetLength(8));

    Allocator::memorycpy(data_, A.GetData(), GetDataSize());
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/

  
  //! Returns the memory used by the object in bytes.
  /*!
    In this method, the type T is assumed to be "static"
    such that sizeof(T) provides the correct size
  */
  template <class T, int N, class Allocator>
  size_t Array<T, N, Allocator>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += size_t(sizeof(T))*GetDataSize();
    return taille;
  }


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the array stores complex
    structures, use 'Fill' instead.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Zero()
  {
    Allocator::memoryset(data_, char(0),
			 GetDataSize()*sizeof(value_type));
  }


  //! Fills the array.
  /*!
    On exit, the array is filled with 1, 2, 3, 4, ... The order of
    those numbers depends on the storage.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Fill()
  {
    for (long i = 0; i < GetDataSize(); i++)
      SetComplexReal(i, data_[i]);
  }


  //! Fills the array with a given value.
  /*!
    On exit, the array is filled with 'x'.
    \param x the value to fill the array with.
  */
  template <class T, int N, class Allocator>
  template <class T0>
  void Array<T, N, Allocator>::Fill(const T0& x)
  {
    T x_;
    SetComplexReal(x, x_);
    for (long i = 0; i < GetDataSize(); i++)
      data_[i] = x_;
  }


  //! Fills the array randomly.
  /*!
    On exit, the array is filled with random values.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::FillRand()
  {
#ifndef SELDON_WITHOUT_REINIT_RANDOM
    srand(time(NULL));
#endif
    for (long i = 0; i < GetDataSize(); i++)
      SetComplexReal(rand(), this->data_[i]);
  }


  //! Displays the array on the standard output.
  /*!
    Displays elements on the standard output, in text format.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Print() const
  {
    if (N == 3)
      {
	int i, j, k;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  cout << (*this)(i, j, k) << '\t';
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 4)
      {
	int i, j, k, l;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      cout << (*this)(i, j, k, l) << '\t';
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 5)
      {
	int i, j, k, l, m;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  cout << (*this)(i, j, k, l, m) << '\t';
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 6)
      {
	int i, j, k, l, m, n;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      cout << (*this)(i, j, k, l, m, n) << '\t';
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 7)
      {
	int i, j, k, l, m, n, o;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  cout << (*this)(i, j, k, l, m, n, o)
				       << '\t';
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 8)
      {
	int i, j, k, l, m, n, o, p;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  {
				    for (p = 0; p < GetLength(7); p++)
				      cout << (*this)(i, j, k, l, m, n, o, p)
					   << '\t';
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 9)
      {
	int i, j, k, l, m, n, o, p, q;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  {
				    for (p = 0; p < GetLength(7); p++)
				      {
					for (q = 0; q < GetLength(8); q++)
					  cout << (*this)(i, j, k, l, m, n,
							  o, p, q)
					       << '\t';
					cout << endl;
				      }
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the array in a file.
  /*!
    Stores the array in a file in binary format. The dimensions (integer)
    are written, and array elements are then written in the same order as in
    memory.
    \param FileName output file name.
  */
  template <class T, int N, class Allocator> void Array<T, N, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Array::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the array to an output stream.
  /*!
    Writes the array to an output stream in binary format. The dimensions
    (integerS) are written, and array elements are then written in the same
    order as in memory.
    \param FileStream output stream.
  */
  template <class T, int N, class Allocator> void Array<T, N, Allocator>
  ::Write(ofstream& FileStream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Array::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    if (with_size)
      {
	for (int i = 0; i < N; i++)
	  FileStream.write(reinterpret_cast<char*>
			   (const_cast<int*>(&length_[i])),
			   sizeof(int));
      }

    FileStream.write(reinterpret_cast<char*>(data_),
		     offset_[N - 1] * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Array::Write(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string("  The output file may have been removed")
		    + " or there is no space left on device.");
#endif

  }


  //! Reads the array from a file.
  /*!
    Reads a array stored in binary format in a file.  The dimensions
    (integers) of the array are read, and array elements are then read in the
    same order as it should be in memory
    \param FileName input file name.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Array::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    Read(FileStream, with_size);

    FileStream.close();
  }


  //! Reads the array from an input stream.
  /*!
    Reads a array in binary format from an input stream.  The dimensions
    (integers) of the array are read, and array elements are then read in the
    same order as it should be in memory.
    \param FileStream input stream.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>
  ::Read(ifstream& FileStream, bool with_size)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Read(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    if (with_size)
      {
	if (N == 3)
	  {
	    int new_l1, new_l2, new_l3;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3);
	  }
	if (N == 4)
	  {
	    int new_l1, new_l2, new_l3, new_l4;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4);
	  }
	if (N == 5)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5);
	  }
	if (N == 6)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6);
	  }
	if (N == 7)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7);
	  }
	if (N == 8)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7,
	      new_l8;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l8), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7, new_l8);
	  }
	if (N == 9)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7,
	      new_l8, new_l9;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l8), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l9), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7, new_l8, new_l9);
	  }
      }

    FileStream.read(reinterpret_cast<char*>(data_), 
		    offset_[N - 1] * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Array::Read(ifstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The intput file may have been removed")
		    + " or may not contain enough data.");
#endif

  }

  //! operator<< overloaded for a 3D array.
  /*!
    \param out output stream.
    \param A the 3D array.
    \return The updated stream.
  */
  template <class T, int N, class Allocator>
  ostream& operator << (ostream& out,
			const Array<T, N, Allocator>& A)
  {
    if (N == 3)
      {
	int i, j, k;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  out << A(i, j, k) << '\t';
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 4)
      {
	int i, j, k, l;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      out << A(i, j, k, l) << '\t';
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 5)
      {
	int i, j, k, l, m;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  out << A(i, j, k, l, m) << '\t';
			out << endl;
		      }
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 6)
      {
	int i, j, k, l, m, n;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      out << A(i, j, k, l, m, n) << '\t';
			    out << endl;
			  }
			out << endl;
		      }
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 7)
      {
	int i, j, k, l, m, n, o;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  cout << A(i, j, k, l, m, n, o)
				       << '\t';
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 8)
      {
	int i, j, k, l, m, n, o, p;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  {
				    for (p = 0; p < A.GetLength(7); p++)
				      cout << A(i, j, k, l, m, n, o, p)
					   << '\t';
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 9)
      {
	int i, j, k, l, m, n, o, p, q;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  {
				    for (p = 0; p < A.GetLength(7); p++)
				      {
					for (q = 0; q < A.GetLength(8); q++)
					  cout << A(i, j, k, l, m, n, o, p, q)
					       << '\t';
					cout << endl;
				      }
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    return out;
  }


} // namespace Seldon.

#define SELDON_FILE_ARRAY_CXX
#endif
