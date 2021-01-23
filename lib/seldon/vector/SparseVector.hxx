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


#ifndef SELDON_FILE_SPARSE_VECTOR_HXX

namespace Seldon
{


  //! Sparse vector class.
  template <class T, class Allocator>
  class Vector<T, VectSparse, Allocator>:
    public Vector<T, VectFull, Allocator>
  {
    // typedef declarations.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename SeldonDefaultAllocator<VectFull, int>::allocator AllocatorInt;
    
    typedef VectSparse storage;

    // Attributes.
  private:
    //! Indices of the non-zero entries.
    int* index_;
    
    // Methods.
  public:
    // Constructor.
    explicit Vector();
    explicit Vector(size_t i);
    Vector(const Vector<T, VectSparse, Allocator>& A);

    // Destructor.
    ~Vector();
    void Clear();

    // Memory management.
    void Reallocate(size_t i);
    void ReallocateVector(size_t i);
    void Resize(size_t i);
    void ResizeVector(size_t i);
    void SetData(size_t nz, T* data, int* index);
    template<class Allocator2>
    void SetData(Vector<T, VectFull, Allocator2>& data,
		 Vector<int>& index);
    template<class Allocator2>
    void SetData(const Vector<T, VectSparse, Allocator2>& V);
    void Nullify();

    // Element access and affectation.
    reference Value(int i);
#ifndef SWIG
    const_reference Value(int i) const;
#endif
    int& Index(int i);
#ifndef SWIG
    int Index(int i) const;
#endif
    
    value_type operator() (int i);
    reference Get(int i);
    reference Val(int i);
#ifndef SWIG
    value_type operator() (int i) const;
    const_reference Get(int i) const;
    const_reference Val(int i) const;
    Vector<T, VectSparse, Allocator>& operator= (const Vector<T, VectSparse,
						 Allocator>& X);
#endif
    template<class T2, class Alloc2>
    void Copy(const Vector<T2, VectSparse, Alloc2>& X);

    // Basic functions.
    int* GetIndex() const;
    size_t GetMemorySize() const;
    
    // Convenient functions.
    template <class T0>
#ifndef SWIG
    Vector<T, VectSparse, Allocator>& operator= (const T0& X);
#endif
    void Print() const;
    void Assemble();
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);
    void AddInteraction(int i, const T& val);
    void AddInteractionRow(size_t, int*, T*, bool already_sorted = false);
    template<class Allocator0>
    void AddInteractionRow(size_t nb, const Vector<int>& col,
			   const Vector<T, VectFull, Allocator0>& val,
                           bool already_sorted = false);
    
    typename ClassComplexType<T>::Treal GetNormInf() const;
    int GetNormInfIndex() const;
    
    // Input/output functions.
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);

  };

#ifndef SWIG
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, VectSparse, Allocator>& V);
#endif


} // namespace Seldon.

#define SELDON_FILE_SPARSE_VECTOR_HXX
#endif
