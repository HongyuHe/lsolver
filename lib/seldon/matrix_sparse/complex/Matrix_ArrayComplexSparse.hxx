// Copyright (C) 2003-2011 Marc Duruflé
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

#ifndef SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX

namespace Seldon
{
  
  //! for complex sparse matrix, the allocator involves real numbers
  template<class T>
  class SeldonDefaultAllocator<ArrayColComplexSparse, T>
  {
  public :
    typedef typename 
    SeldonDefaultAllocator<VectFull, typename ClassComplexType<T>::Treal>
    ::allocator allocator;    
  };

  template<>
  class SeldonDefaultAllocator<ArrayColComplexSparse, complex<float> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, float>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayColComplexSparse, complex<double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, double>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayColComplexSparse, complex<long double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, long double>::allocator allocator;
  };

  //! for complex sparse matrix, the allocator involves real numbers
  template<class T>
  class SeldonDefaultAllocator<ArrayRowComplexSparse, T>
  {
  public :
    typedef typename 
    SeldonDefaultAllocator<VectFull, typename ClassComplexType<T>::Treal>
    ::allocator allocator;    
  };

  template<>
  class SeldonDefaultAllocator<ArrayRowComplexSparse, complex<float> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, float>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayRowComplexSparse, complex<double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, double>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayRowComplexSparse, complex<long double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, long double>::allocator allocator;
  };


  //! for complex sparse matrix, the allocator involves real numbers
  template<class T>
  class SeldonDefaultAllocator<ArrayColSymComplexSparse, T>
  {
  public :
    typedef typename 
    SeldonDefaultAllocator<VectFull, typename ClassComplexType<T>::Treal>
    ::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayColSymComplexSparse, complex<float> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, float>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayColSymComplexSparse, complex<double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, double>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayColSymComplexSparse, complex<long double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, long double>::allocator allocator;
  };


  //! for complex sparse matrix, the allocator involves real numbers
  template<class T>
  class SeldonDefaultAllocator<ArrayRowSymComplexSparse, T>
  {
  public :
    typedef typename 
    SeldonDefaultAllocator<VectFull, typename ClassComplexType<T>::Treal>
    ::allocator allocator;    
  };

  template<>
  class SeldonDefaultAllocator<ArrayRowSymComplexSparse, complex<float> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, float>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayRowSymComplexSparse, complex<double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, double>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ArrayRowSymComplexSparse, complex<long double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, long double>::allocator allocator;
  };

  
  //! Sparse Array-matrix class.
  /*!
    Sparse matrices are defined by: (1) the number of rows and columns;
    (2) the number of non-zero entries; (3) an array of vectors ind
    ind(i) is a vector, which contains indices of columns of the row i;
    (4) an array of vectors val : val(i) is a vector, which contains values of
    the row i
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class Matrix_ArrayComplexSparse : public VirtualMatrix<T>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef complex<value_type> entry_type;
    typedef complex<value_type> access_type;
    typedef complex<value_type> const_access_type;

    // Attributes.
  protected:
    //! real part rows or columns
    Vector<Vector<value_type, VectSparse, Allocator>, VectFull,
	   NewAlloc<Vector<value_type, VectSparse, Allocator> > > val_real_;
    //! imaginary part rows or columns
    Vector<Vector<value_type, VectSparse, Allocator>, VectFull,
	   NewAlloc<Vector<value_type, VectSparse, Allocator> > > val_imag_;

    // Methods.
  public:
    // Constructors.
    Matrix_ArrayComplexSparse();
    explicit Matrix_ArrayComplexSparse(int i, int j);

    // Destructor.
    ~Matrix_ArrayComplexSparse();
    void Clear();

    // Memory management.
    void Reallocate(int i, int j);
    void Resize(int i, int j);

    // Basic methods.
    long GetRealNonZeros() const;
    long GetImagNonZeros() const;
    long GetNonZeros() const;
    long GetRealDataSize() const;
    long GetImagDataSize() const;
    long GetDataSize() const;
    size_t GetMemorySize() const;
    int* GetRealInd(int i) const;
    int* GetImagInd(int i) const;
    value_type* GetRealData(int i) const;
    value_type* GetImagData(int i) const;
    Vector<value_type, VectSparse, Allocator>* GetRealData() const;
    Vector<value_type, VectSparse, Allocator>* GetImagData() const;

    // Element acess and affectation.
    const entry_type operator() (int i, int j) const;
    entry_type& Val(int i, int j);
    const entry_type& Val(int i, int j) const;
    entry_type& Get(int i, int j);
    const entry_type& Get(int i, int j) const;

    value_type& ValReal(int i, int j);
    const value_type& ValReal(int i, int j) const;
    value_type& ValImag(int i, int j);
    const value_type& ValImag(int i, int j) const;
    value_type& GetReal(int i, int j);
    const value_type& GetReal(int i, int j) const;
    value_type& GetImag(int i, int j);
    const value_type& GetImag(int i, int j) const;    
    
    void Set(int i, int j, const entry_type& x);
    
    const value_type& ValueReal(int num_row,int i) const;
    value_type& ValueReal(int num_row,int i);
    int IndexReal(int num_row,int i) const;
    int& IndexReal(int num_row,int i);
    const value_type& ValueImag(int num_row,int i) const;
    value_type& ValueImag(int num_row,int i);
    int IndexImag(int num_row,int i) const;
    int& IndexImag(int num_row,int i);

    void SetRealData(int, int, Vector<value_type, VectSparse, Allocator>*);
    void SetImagData(int, int, Vector<value_type, VectSparse, Allocator>*);
    void SetRealData(int, int, value_type*, int*);
    void SetImagData(int, int, value_type*, int*);
    void NullifyReal(int i);
    void NullifyImag(int i);
    void NullifyReal();
    void NullifyImag();

    // Convenient functions.
    void Print() const;
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName, bool cplx = false) const;
    void WriteText(ostream& FileStream, bool cplx = false) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName, bool cplx = false);
    void ReadText(istream& FileStream, bool cplx = false);
    
    void Assemble();
    
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);

    void SetIdentity();
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const complex<T0>& x);
    template <class T0>
    Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>& operator=
    (const complex<T0>& x);
    void FillRand();

#ifdef SELDON_WITH_VIRTUAL
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

    virtual void ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
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
  class Matrix<T, Prop, ArrayColComplexSparse, Allocator> :
    public Matrix_ArrayComplexSparse<T, Prop, ArrayColComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayColComplexSparse storage;
    typedef Allocator allocator;
    typedef complex<value_type> entry_type;
    
  public:
    Matrix();
    explicit Matrix(int i, int j);

    // Memory management.
    void ClearRealColumn(int i);
    void ClearImagColumn(int i);
    void ReallocateRealColumn(int i, int j);
    void ReallocateImagColumn(int i, int j);
    void ResizeRealColumn(int i, int j);
    void ResizeImagColumn(int i, int j);
    void SwapRealColumn(int i, int i_);
    void SwapImagColumn(int i, int i_);
    void ReplaceRealIndexColumn(int i, IVect& new_index);
    void ReplaceImagIndexColumn(int i, IVect& new_index);

    int GetRealColumnSize(int i) const;
    int GetImagColumnSize(int i) const;
    void PrintRealColumn(int i) const;
    void PrintImagColumn(int i) const;
    void AssembleRealColumn(int i);
    void AssembleImagColumn(int i);

    void AddInteraction(int i, int j, const entry_type& val);

    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<entry_type>& val,
			   bool sorted = false);

    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<entry_type>& val,
			      bool sorted = false);
  };


  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowComplexSparse, Allocator> :
    public Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayRowComplexSparse storage;
    typedef Allocator allocator;
    typedef complex<value_type> entry_type;
    
  public:
    Matrix();
    explicit Matrix(int i, int j);

    // Memory management.
    void ClearRealRow(int i);
    void ClearImagRow(int i);
    void ClearRow(int i);
    void ReallocateRealRow(int i, int j);
    void ReallocateImagRow(int i, int j);
    void ResizeRealRow(int i, int j);
    void ResizeImagRow(int i, int j);
    void SwapRealRow(int i, int i_);
    void SwapImagRow(int i, int i_);
    void ReplaceRealIndexRow(int i, IVect& new_index);
    void ReplaceImagIndexRow(int i, IVect& new_index);

    int GetRealRowSize(int i) const;
    int GetImagRowSize(int i) const;
    void PrintRealRow(int i) const;
    void PrintImagRow(int i) const;
    void AssembleRealRow(int i);
    void AssembleImagRow(int i);

    void AddInteraction(int i, int j, const entry_type& val);

    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<entry_type>& val,
			   bool sorted = false);

    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<entry_type>& val,
			      bool sorted = false);
    
  };


  //! Column-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>:
    public Matrix_ArrayComplexSparse<T, Prop, ArrayColSymComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayColSymComplexSparse storage;
    typedef Allocator allocator;
    typedef complex<value_type> entry_type;
    
  public:
    Matrix();
    explicit Matrix(int i, int j);

    const entry_type operator() (int i, int j) const;
    
    value_type& ValReal(int i, int j);
    const value_type& ValReal(int i, int j) const;
    value_type& ValImag(int i, int j);
    const value_type& ValImag(int i, int j) const;
    value_type& GetReal(int i, int j);
    const value_type& GetReal(int i, int j) const;
    value_type& GetImag(int i, int j);
    const value_type& GetImag(int i, int j) const;    
    
    void Set(int i, int j, const entry_type& x);

    // Memory management.
    void ClearRealColumn(int i);
    void ClearImagColumn(int i);
    void ReallocateRealColumn(int i, int j);
    void ReallocateImagColumn(int i, int j);
    void ResizeRealColumn(int i, int j);
    void ResizeImagColumn(int i, int j);
    void SwapRealColumn(int i, int i_);
    void SwapImagColumn(int i, int i_);
    void ReplaceRealIndexColumn(int i, IVect& new_index);
    void ReplaceImagIndexColumn(int i, IVect& new_index);

    int GetRealColumnSize(int i) const;
    int GetImagColumnSize(int i) const;
    void PrintRealColumn(int i) const;
    void PrintImagColumn(int i) const;
    void AssembleRealColumn(int i);
    void AssembleImagColumn(int i);

    void AddInteraction(int i, int j, const entry_type& val);

    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<entry_type>& val,
			   bool sorted = false);
    
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<entry_type>& val,
			      bool sorted = false);
    
  };


  //! Row-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>:
    public Matrix_ArrayComplexSparse<T, Prop, ArrayRowSymComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ArrayRowSymComplexSparse storage;
    typedef Allocator allocator;
    typedef complex<value_type> entry_type;
    
  public:
    Matrix();
    explicit Matrix(int i, int j);

    const entry_type operator() (int i, int j) const;
    
    value_type& ValReal(int i, int j);
    const value_type& ValReal(int i, int j) const;
    value_type& ValImag(int i, int j);
    const value_type& ValImag(int i, int j) const;
    value_type& GetReal(int i, int j);
    const value_type& GetReal(int i, int j) const;
    value_type& GetImag(int i, int j);
    const value_type& GetImag(int i, int j) const;    
    
    void Set(int i, int j, const entry_type& x);

    // Memory management.
    void ClearRealRow(int i);
    void ClearImagRow(int i);
    void ClearRow(int i);
    void ReallocateRealRow(int i, int j);
    void ReallocateImagRow(int i, int j);
    void ResizeRealRow(int i, int j);
    void ResizeImagRow(int i, int j);
    void SwapRealRow(int i, int i_);
    void SwapImagRow(int i, int i_);
    void ReplaceRealIndexRow(int i, IVect& new_index);
    void ReplaceImagIndexRow(int i, IVect& new_index);

    int GetRealRowSize(int i) const;
    int GetImagRowSize(int i) const;
    void PrintRealRow(int i) const;
    void PrintImagRow(int i) const;
    void AssembleRealRow(int i);
    void AssembleImagRow(int i);

    void AddInteraction(int i, int j, const entry_type& val);

    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<entry_type>& val,
			   bool sorted = false);
    
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<entry_type>& val,
			      bool sorted= false);
  };
  
} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX
#endif
