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


#ifndef SELDON_FILE_ARRAY_HXX


namespace Seldon
{


  //! Multi-dimensional Array.
  /*!
    This class implements multi-dimensional arrays.
  */
  template <class T, int N, class Allocator>
  class Array
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

    // Length along dimensions.
    int length_[N];
    // Sizes of slices for dimension 1, 2, ..., N.
    long offset_[N];
    // Pointer to stored elements.
    pointer data_;

    // Methods.
  public:
    // Constructors.
    Array();
    Array(int i);
    Array(int i, int j, int k);
    Array(int i, int j, int k, int l);
    Array(int i, int j, int k, int l, int m);
    Array(int i, int j, int k, int l, int m, int n);
    Array(int i, int j, int k, int l, int m, int n, int o);
    Array(int i, int j, int k, int l, int m, int n, int o, int p);
    Array(int i, int j, int k, int l, int m, int n, int o, int p, int q);
    // Array(int* extent);
    Array(const Array<T, N, Allocator>& A);

    // Destructor. (inline)
    ~Array();

    // Inline methods.
    int GetLength(int dimension) const;
    long GetSize() const;
    long GetDataSize() const;
    pointer GetData() const;

    // Memory management.
    void Reallocate(int i, int j, int k);
    void Reallocate(int i, int j, int k, int l);
    void Reallocate(int i, int j, int k, int l, int m);
    void Reallocate(int i, int j, int k, int l, int m, int n);
    void Reallocate(int i, int j, int k, int l, int m, int n, int o);
    void Reallocate(int i, int j, int k, int l, int m, int n, int o, int p);
    void Reallocate(int i, int j, int k, int l, int m, int n, int o, int p,
		    int q);
    void Clear();

    // Element access and affectation. (inline)
    reference operator() (int i, int j, int k);
    reference operator() (int i, int j, int k, int l);
    reference operator() (int i, int j, int k, int l, int m);
    reference operator() (int i, int j, int k, int l, int m, int n);
    reference operator() (int i, int j, int k, int l, int m, int n, int o);
    reference operator() (int i, int j, int k, int l, int m, int n, int o,
			  int p);
    reference operator() (int i, int j, int k, int l, int m, int n, int o,
			  int p, int q);
#ifndef SWIG
    const_reference operator() (int i, int j, int k) const;
    const_reference operator() (int i, int j, int k, int l) const;
    const_reference operator() (int i, int j, int k, int l, int m) const;
    const_reference operator()
    (int i, int j, int k, int l, int m, int n) const;
    const_reference operator()
    (int i, int j, int k, int l, int m, int n, int o) const;
    const_reference operator()
    (int i, int j, int k, int l, int m, int n, int o, int p) const;
    const_reference operator()
    (int i, int j, int k, int l, int m, int n, int o, int p, int q) const;
    Array<T, N, Allocator>& operator= (const Array<T, N, Allocator>& A);
#endif
    void Copy(const Array<T, N, Allocator>& A);

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
  template <class T, int N, class Allocator>
  ostream& operator << (ostream& out,
			const Array<T, N, Allocator>& A);
#endif


} // namespace Seldon.


#define SELDON_FILE_ARRAY_HXX
#endif
