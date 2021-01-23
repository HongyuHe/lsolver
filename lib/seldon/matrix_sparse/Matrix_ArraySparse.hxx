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


// To be included by Seldon.hxx

#ifndef SELDON_FILE_MATRIX_ARRAY_SPARSE_HXX

namespace Seldon
{

  //! Sparse Array-matrix class.
  /*!
    Sparse matrices are defined by: (1) the number of rows and columns;
    (2) the number of non-zero entries; (3) an array of vectors ind
    ind(i) is a vector, which contains indices of columns of the row i
    (4) an array of vectors val : val(i) is a vector, which contains values of
    the row i
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class Matrix_ArraySparse : public VirtualMatrix<T>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef T entry_type;
    typedef T& access_type;
    typedef T const_access_type;

    // Attributes.
  protected:
    //! rows or columns
    Vector<Vector<T, VectSparse, Allocator>, VectFull,
	   NewAlloc<Vector<T, VectSparse, Allocator> > > val_;

  public:
    // Constructors. (inline)
    Matrix_ArraySparse();
    explicit Matrix_ArraySparse(int i, int j);

    // Destructor. (inline)
    ~Matrix_ArraySparse();
    void Clear();

    // Memory management.
    void Reallocate(int i, int j);
    void Resize(int i, int j);

    // Basic methods.
    long GetNonZeros() const;
    size_t GetMemorySize() const;

    // Inline methods.
    long GetDataSize() const;
    int* GetIndex(int i) const;
    T* GetData(int i) const;

    Vector<T, VectSparse, Allocator>* GetData() const;

    // Element acess and affectation. (inline)
#ifdef SELDON_WITH_MODIFIABLE_PARENTHESIS_OPERATOR
    T& operator() (int i, int j);
#endif
    const T operator() (int i, int j) const;    
    T& Get(int i, int j);
    const T& Get(int i, int j) const;
    T& Val(int i, int j);
    const T& Val(int i, int j) const;
    void Set(int i, int j, const T& x);
    void SetEntry(int i, int j, const T& x);

    const T& Value(int num_row, int i) const;
    T& Value(int num_row, int i);
    int Index(int num_row, int i) const;
    int& Index(int num_row, int i);

    void SetData(int, int, Vector<T, VectSparse, Allocator>*);
    void SetData(int, int, T*, int*);
    void Nullify(int i);
    void Nullify();

    // Convenient functions.
    void Print() const;
    void Assemble();
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);

    void SetIdentity();
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    template <class T0>
    Matrix_ArraySparse<T, Prop, Storage, Allocator>& operator= (const T0& x);
    void FillRand();

    // Input/output functions.
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName, bool cplx = false) const;
    void WriteText(ostream& FileStream, bool cplx = false) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName, bool cplx = false);
    void ReadText(istream& FileStream, bool cplx = false);

#ifdef SELDON_WITH_VIRTUAL
    // Inline virtual methods
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
  class Matrix<T, Prop, ArrayColSparse, Allocator> :
    public Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayColSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(int i, int j);

    // Memory management.
    void ClearColumn(int i);
    void ReallocateColumn(int i, int j);
    void ResizeColumn(int i, int j);
    void SwapColumn(int i, int i_);
    void ReplaceIndexColumn(int i, IVect& new_index);

    int GetColumnSize(int i) const;
    void PrintColumn(int i) const;
    void AssembleColumn(int i);

    void AddInteraction(int i, int j, const T& val);

    void AddInteractionRow(int, int, int*, T*, bool already_sorted = false);
    void AddInteractionColumn(int, int, int*, T*, bool already_sorted = false);

    void AddInteractionRow(int i, int nb, const Vector<int>& col,
			   const Vector<T>& val, bool already_sorted = false);
    
    void AddInteractionColumn(int i, int nb, const Vector<int>& row,
			      const Vector<T>& val,
                              bool already_sorted = false);
    
  };


  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSparse, Allocator> :
    public Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayRowSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(int i, int j);

    // Memory management.
    void ClearRow(int i);
    void ReallocateRow(int i, int j);
    void ResizeRow(int i, int j);
    void SwapRow(int i, int i_);
    void ReplaceIndexRow(int i, IVect& new_index);

    int GetRowSize(int i) const;
    void PrintRow(int i) const;
    void AssembleRow(int i);

    void AddInteraction(int i, int j, const T& val);

    void AddInteractionRow(int, int, int*, T*, bool already_sorted = false);
    void AddInteractionColumn(int, int, int*, T*, bool already_sorted = false);

    void AddInteractionRow(int i, int nb, const Vector<int>& col,
			   const Vector<T>& val, bool already_sorted = false);

    void AddInteractionColumn(int i, int nb, const Vector<int>& row,
			      const Vector<T>& val, bool already_sorted = false);
    
  };

  //! Column-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSymSparse, Allocator>:
    public Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayColSymSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(int i, int j);
    
    // access operator
#ifdef SELDON_WITH_MODIFIABLE_PARENTHESIS_OPERATOR
    T& operator() (int i, int j);
#endif
    const T operator() (int i, int j) const;
    
    T& Get(int i, int j);
    const T& Get(int i, int j) const;
    T& Val(int i, int j);
    const T& Val(int i, int j) const;
    void Set(int i, int j, const T& x);
    void SetEntry(int i, int j, const T& x);

    // Memory management.
    void ClearColumn(int i);
    void ReallocateColumn(int i, int j);
    void ResizeColumn(int i, int j);
    void SwapColumn(int i, int i_);
    void ReplaceIndexColumn(int i, IVect& new_index);

    int GetColumnSize(int i) const;
    void PrintColumn(int i) const;
    void AssembleColumn(int i);

    void AddInteraction(int i, int j, const T& val);

    void AddInteractionRow(int, int, int*, T*, bool already_sorted = false);
    void AddInteractionColumn(int, int, int*, T*,
                              bool already_sorted = false);

    void AddInteractionRow(int i, int nb, const Vector<int>& col,
			   const Vector<T>& val,
			   bool already_sorted = false);

    void AddInteractionColumn(int i, int nb, const Vector<int>& row,
			      const Vector<T>& val,
                              bool already_sorted = false);
    
  };


  //! Row-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymSparse, Allocator>:
    public Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayRowSymSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(int i, int j);

    // access operator
#ifdef SELDON_WITH_MODIFIABLE_PARENTHESIS_OPERATOR
    T& operator() (int i, int j);
#endif
    const T operator() (int i, int j) const;

    T& Get(int i, int j);
    const T& Get(int i, int j) const;
    T& Val(int i, int j);
    const T& Val(int i, int j) const;
    void Set(int i, int j, const T& x);
    void SetEntry(int i, int j, const T& x);

    // Memory management.
    void ClearRow(int i);
    void ReallocateRow(int i, int j);
    void ResizeRow(int i, int j);
    void SwapRow(int i, int i_);
    void ReplaceIndexRow(int i, IVect& new_index);

    int GetRowSize(int i) const;
    void PrintRow(int i) const;
    void AssembleRow(int i);

    void AddInteraction(int i, int j, const T& val);

    void AddInteractionRow(int, int, int*, T*, bool already_sorted = false);
    void AddInteractionColumn(int, int, int*, T*, bool sorted = false);

    void AddInteractionRow(int i, int nb, const Vector<int>& col,
			   const Vector<T>& val, bool already_sorted = false);
    
    void AddInteractionColumn(int i, int nb, const Vector<int>& row,
			      const Vector<T>& val, bool already_sorted = false);
    
  };


  template <class T, class Prop, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix<T, Prop, ArrayRowSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix<T, Prop, ArrayColSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A);
  
} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_SPARSE_HXX
#endif
