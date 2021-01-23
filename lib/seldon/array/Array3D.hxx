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


// To be included by Seldon.hxx


#ifndef SELDON_FILE_ARRAY3D_HXX

#include "../share/Common.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! 3D array.
  /*!
    This class implements 3D arrays.
  */
  template <class T, class Allocator>
  class Array3D
  {
    // typdef declarations.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    // Attributes.
  protected:
    // Length along dimension #1.
    int length1_;
    // Length along dimension #2.
    int length2_;
    // Length along dimension #3.
    int length3_;

    // Size of a slice (i.e. length1_ by length2_).
    long length23_;

    // Pointer to stored elements.
    pointer data_;

    // Methods.
  public:
    // Constructors.
    Array3D();
    Array3D(int i, int j, int k);
    Array3D(const Array3D<T, Allocator>& A);

    // Destructor. (inline)
    ~Array3D();

    // Basic methods. (inline)
    int GetLength1() const;
    int GetLength2() const;
    int GetLength3() const;
    long GetSize() const;
    long GetDataSize() const;
    pointer GetData() const;
    pointer GetDataPointer(int i, int j, int k) const;

    // Memory management.
    void Reallocate(int i, int j, int k);
    void SetData(int i, int j, int k, pointer data);
    void Nullify();
    void Clear();

    // Element access and affectation. (inline)
    reference operator() (int i, int j, int k);
#ifndef SWIG
    const_reference operator() (int i, int j, int k) const;
    Array3D<T, Allocator>& operator= (const Array3D<T, Allocator>& A);
#endif
    void Copy(const Array3D<T, Allocator>& A);

    // Convenient functions.
    size_t GetMemorySize() const;
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    void FillRand();
    void Print() const;

    // Input/output functions
    void Write(string FileName, bool with_size = true) const;
    void Write(ofstream& FileStream, bool with_size = true) const;
    void Read(string FileName, bool with_size = true);
    void Read(ifstream& FileStream, bool with_size = true);
  };


#ifndef SWIG
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
			const Array3D<T, Allocator>& A);
#endif

  template <class T0, class T, class Allocator>
  void MltScalar(const T0& alpha, Array3D<T, Allocator>& A);
  
} // namespace Seldon.


#define SELDON_FILE_ARRAY3D_HXX
#endif
