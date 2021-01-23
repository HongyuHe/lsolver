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


#ifndef SELDON_FILE_ARRAY_INLINE_CXX


namespace Seldon
{
  
  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>::~Array()
  {
    Clear();
  }


  /*****************
   * BASIC METHODS *
   *****************/

  //! Returns the length in dimension #1.
  /*!
    \param dimension index for dimension.
    \return The length in dimension #1.
  */
  template <class T, int N, class Allocator>
  inline int Array<T, N, Allocator>::GetLength(int dimension) const
  {
    return length_[dimension];
  }


  //! Returns the number of elements in the 3D array.
  /*!
    Returns the number of elements stored by the 3D array, i.e.
    the product of the lengths in the three dimensions.
    \return The number of elements in the 3D array.
  */
  template <class T, int N, class Allocator>
  inline long Array<T, N, Allocator>::GetSize() const
  {
    return offset_[N - 1];
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory by
    the array, i.e. the product of lengths in the three
    dimensions.
    \return The number of elements stored in the array.
  */
  template <class T, int N, class Allocator>
  inline long Array<T, N, Allocator>::GetDataSize() const
  {
    return offset_[N - 1];
  }


  //! Returns a pointer to the data array.
  /*!
    Returns a pointer to data, i.e. the data array 'data_' which stores the
    values.
    \return A pointer to the data array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::pointer Array<T, N, Allocator>
  ::GetData() const
  {
    return data_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>::operator() (int i, int j, int k)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 3)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 3.");

    CheckBounds(i, j, k, length_[0], length_[1], length_[2], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>::operator() (int i, int j, int k) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 3)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 3.");

    CheckBounds(i, j, k, length_[0], length_[1], length_[2], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \return Element (i, j, k, l) of the 4D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 4)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 4.");

    CheckBounds(i, j, k, l, length_[0], length_[1],
		length_[2], length_[3], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] + l];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \return Element (i, j, k, l) of the 4D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 4)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 4.");

    CheckBounds(i, j, k, l, length_[0], length_[1],
		length_[2], length_[3], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] + l];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \return Element (i, j, k, l, m) of the 5D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l, int m)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 5)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 5.");

    CheckBounds(i, j, k, l, m, length_[0], length_[1], length_[2],
		length_[3], length_[4], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \return Element (i, j, k, l, m) of the 5D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l, int m) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 5)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 5.");

    CheckBounds(i, j, k, l, m, length_[0], length_[1], length_[2],
		length_[3], length_[4], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \return Element (i, j, k, l, m, n) of the 6D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 6)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 6.");

    CheckBounds(i, j, k, l, m, n, length_[0], length_[1], length_[2],
		length_[3], length_[4], length_[5], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \return Element (i, j, k, l, m, n) of the 6D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 6)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 6.");

    CheckBounds(i, j, k, l, m, n, length_[0], length_[1], length_[2],
		length_[3], length_[4], length_[5], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \return Element (i, j, k, l, m, n, o) of the 7D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n, int o)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 7)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 7.");

    CheckBounds(i, j, k, l, m, n, o, length_[0], length_[1], length_[2],
		length_[3], length_[4], length_[5], length_[6], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] + o];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \return Element (i, j, k, l, m, n, o) of the 7D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n, int o) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 7)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 7.");

    CheckBounds(i, j, k, l, m, n, o, length_[0], length_[1], length_[2],
		length_[3], length_[4], length_[5], length_[6], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] + o];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o, p).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \param p index along dimension #8.
    \return Element (i, j, k, l, m, n, o, p) of the 8D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n, int o, int p)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 8)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 8.");

    CheckBounds(i, j, k, l, m, n, o, p, length_[0], length_[1], length_[2],
		length_[3], length_[4], length_[5], length_[6],
		length_[7], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] +
		 o * offset_[6] + p];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o, p).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \param p index along dimension #8.
    \return Element (i, j, k, l, m, n, o, p) of the 8D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n, int o, int p) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 8)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 8.");

    CheckBounds(i, j, k, l, m, n, o, p, length_[0], length_[1], length_[2],
		length_[3], length_[4], length_[5], length_[6],
		length_[7], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] +
		 o * offset_[6] + p];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o, p, q).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \param p index along dimension #8.
    \param q index along dimension #9.
    \return Element (i, j, k, l, m, n, o, p, q) of the 9D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l, int m,
				      int n, int o, int p, int q)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 9)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 9.");

    CheckBounds(i, j, k, l, m, n, o, p, q, length_[0], length_[1], length_[2],
		length_[3], length_[4], length_[5], length_[6],
		length_[7], length_[8], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] +
		 o * offset_[6] + p * offset_[7] + q];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o, p, q).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \param p index along dimension #8.
    \param q index along dimension #9.
    \return Element (i, j, k, l, m, n, o, p, q) of the 9D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l, int m,
				      int n, int o, int p, int q) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (N != 9)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 9.");

    CheckBounds(i, j, k, l, m, n, o, p, q, length_[0], length_[1], length_[2],
		length_[3], length_[4], length_[5], length_[6],
		length_[7], length_[8], "Array");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] +
		 o * offset_[6] + p * offset_[7] + q];
  }


  //! Duplicates an array (assignment operator).
  /*!
    \param A array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>& Array<T, N, Allocator>::operator=
  (const Array<T, N, Allocator>& A)
  {
    Copy(A);

    return *this;
  }
  
} // namespace Seldon.

#define SELDON_FILE_ARRAY_INLINE_CXX
#endif
