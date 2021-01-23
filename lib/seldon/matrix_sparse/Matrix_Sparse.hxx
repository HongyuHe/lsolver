// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2001-2011 Marc Duruflé
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

#ifndef SELDON_FILE_MATRIX_SPARSE_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Sparse-matrix class.
  /*!
    Sparse matrices are defined by: (1) the number of rows and columns;
    (2) the number of non-zero entries; (3) an array 'ptr_' of start indices
    (i.e. indices of the first element of each row or column, depending
    on the storage); (4) an array 'ind_' of column or row indices of each
    non-zero entry; (5) values of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class Matrix_Sparse: public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename SeldonDefaultAllocator<VectFull, int>::allocator AllocatorInt;
    typedef typename SeldonDefaultAllocator<VectFull, long>::allocator AllocatorLong;
    typedef value_type entry_type;
    typedef value_type access_type;
    typedef value_type const_access_type;

    // Attributes.
  protected:
    // Number of non-zero elements.
    long nz_;
    // Index (in data_) of first element stored for each row or column.
    long* ptr_;
    // Column or row index (in the matrix) each element.
    int* ind_;

    // Methods.
  public:
    // Constructors.
    Matrix_Sparse();
    explicit Matrix_Sparse(int i, int j);
    explicit Matrix_Sparse(int i, int j, long nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix_Sparse(int i, int j, Vector<T, Storage0, Allocator0>& values,
		  Vector<long, Storage1, Allocator1>& ptr,
		  Vector<int, Storage2, Allocator2>& ind);
    Matrix_Sparse(const Matrix_Sparse<T, Prop, Storage, Allocator>& A);

    // Destructor.
    ~Matrix_Sparse();
    void Clear();

    // Memory management.
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    void SetData(int i, int j,
		 Vector<T, Storage0, Allocator0>& values,
		 Vector<long, Storage1, Allocator1>& ptr,
		 Vector<int, Storage2, Allocator2>& ind);
    void SetData(int i, int j, long nz, pointer values, long* ptr, int* ind);
    void Nullify();
    void Reallocate(int i, int j);
    void Reallocate(int i, int j, long nz);
    void Resize(int i, int j);
    void Resize(int i, int j, long nz);
    void Copy(const Matrix_Sparse<T, Prop, Storage, Allocator>& A);

    // Basic methods.
    long GetNonZeros() const;
    long GetDataSize() const;
    size_t GetMemorySize() const;
    long* GetPtr() const;
    int* GetInd() const;
    int GetPtrSize() const;
    long GetIndSize() const;
    void FillPtrInt(Vector<int>&) const;

    // Element acess and affectation.
    const value_type operator() (int i, int j) const;
    value_type& Val(int i, int j);
    value_type& Get(int i, int j);
#ifndef SWIG
    const value_type& Val(int i, int j) const;
    const value_type& Get(int i, int j) const;
#endif
    void AddInteraction(int i, int j, const T& val);
    
    void AddInteractionRow(int i, int nb, const Vector<int>& col,
			   const Vector<T>& val, bool sorted = false);

    void Set(int i, int j, const T& x);
    void SetEntry(int i, int j, const T& x);

#ifndef SWIG
    Matrix_Sparse<T, Prop, Storage, Allocator>&
    operator= (const Matrix_Sparse<T, Prop, Storage, Allocator>& A);
#endif

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    void FillRand();
    void FillRand(long Nelement);
    void FillRand(long Nelement, const T& x);
    
    void Print() const;
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName, bool cplx = false) const;
    void WriteText(ostream& FileStream, bool cplx = false) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName, bool cplx = false);
    void ReadText(istream& FileStream, bool cplx = false);

#ifdef SELDON_WITH_VIRTUAL
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

    virtual void ApplySor(const SeldonTranspose&, Vector<Treal>& x, const Vector<Treal>& r,
			  const typename ClassComplexType<T>::Treal& omega,
			  int nb_iter, int stage_ssor) const;

    virtual void ApplySor(const SeldonTranspose&, Vector<Tcplx>& x, const Vector<Tcplx>& r,
			  const typename ClassComplexType<T>::Treal& omega,
			  int nb_iter, int stage_ssor) const;

    virtual void MltAddVector(const Treal& alpha, const Vector<Treal>& x,
			      const Treal& beta, Vector<Treal>& y) const;

    virtual void MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
			      const Tcplx& beta, Vector<Tcplx>& y) const;

    virtual void MltAddVector(const Treal& alpha, const SeldonTranspose&,
			      const Vector<Treal>& x,
			      const Treal& beta, Vector<Treal>& y) const;

    virtual void MltAddVector(const Tcplx& alpha, const SeldonTranspose&,
			      const Vector<Tcplx>& x,
			      const Tcplx& beta, Vector<Tcplx>& y) const;
    
    virtual void MltVector(const Vector<Treal>& x, Vector<Treal>& y) const;
    virtual void MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const;
    
    virtual void MltVector(const SeldonTranspose&,
			   const Vector<Treal>& x, Vector<Treal>& y) const;

    virtual void MltVector(const SeldonTranspose&,
			   const Vector<Tcplx>& x, Vector<Tcplx>& y) const;
    
    virtual bool IsSymmetric() const;
#endif
    
  };


  //! Column-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSparse, Allocator>:
    public Matrix_Sparse<T, Prop, ColSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ColSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(int i, int j);
    explicit Matrix(int i, int j, long nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(int i, int j,
	   Vector<T, Storage0, Allocator0>& values,
	   Vector<long, Storage1, Allocator1>& ptr,
	   Vector<int, Storage2, Allocator2>& ind);
  };


  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSparse, Allocator>:
    public Matrix_Sparse<T, Prop, RowSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef RowSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(int i, int j);
    explicit Matrix(int i, int j, long nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(int i, int j,
	   Vector<T, Storage0, Allocator0>& values,
	   Vector<long, Storage1, Allocator1>& ptr,
	   Vector<int, Storage2, Allocator2>& ind);

  };

  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_SPARSE_HXX
#endif
