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

#ifndef SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_HXX

namespace Seldon
{

  //! for complex sparse matrix, the allocator involves real numbers
  template<class T>
  class SeldonDefaultAllocator<ColSymComplexSparse, T>
  {
  public :
    typedef typename
    SeldonDefaultAllocator<VectFull, typename ClassComplexType<T>::Treal>
    ::allocator allocator;    
  };

  template<>
  class SeldonDefaultAllocator<ColSymComplexSparse, complex<float> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, float>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ColSymComplexSparse, complex<double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, double>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<ColSymComplexSparse, complex<long double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, long double>::allocator allocator;
  };


  //! for complex sparse matrix, the allocator involves real numbers
  template<class T>
  class SeldonDefaultAllocator<RowSymComplexSparse, T>
  {
  public :
    typedef typename
    SeldonDefaultAllocator<VectFull, typename ClassComplexType<T>::Treal>
    ::allocator allocator;    
  };

  template<>
  class SeldonDefaultAllocator<RowSymComplexSparse, complex<float> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, float>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<RowSymComplexSparse, complex<double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, double>::allocator allocator;
  };

  template<>
  class SeldonDefaultAllocator<RowSymComplexSparse, complex<long double> >
  {
  public:
    typedef
    SeldonDefaultAllocator<VectFull, long double>::allocator allocator;
  };

  
  //! Symmetric complex sparse-matrix class.
  /*!
    Symmetric sparse matrices are defined by: (1) the number of rows
    and columns; (2) the number of non-zero entries; (3) an array 'ptr_' of
    start indices (i.e. indices of the first element of each row or column,
    depending on the storage); (4) an array 'ind_' of column or row indices
    of each non-zero entry; (5) values of non-zero entries.\par
    Complex sparse matrices are defined in the same way except that real
    and imaginary parts are splitted. It is as if two matrices were stored.
    There are therefore 6 arrays: 'real_ptr_', 'real_ind_', 'real_data_',
    'imag_ptr_', 'imag_ind_' and 'imag_data_'.\par
    Finally, since the matrix is symmetric, only its upper part is stored.
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class Matrix_SymComplexSparse: public Matrix_Base<T, Allocator>
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
    typedef complex<value_type> entry_type;
    typedef complex<value_type> access_type;
    typedef complex<value_type> const_access_type;

    // Attributes.
  protected:
    // Number of non-zero (stored) elements.
    long real_nz_;
    long imag_nz_;
    // Index (in data_) of first element stored for each row or column.
    long* real_ptr_;
    long* imag_ptr_;
    // Column or row index (in the matrix) each element.
    int* real_ind_;
    int* imag_ind_;

    // Data.
    value_type* real_data_;
    value_type* imag_data_;

    // Methods.
  public:
    // Constructors.
    Matrix_SymComplexSparse();
    explicit Matrix_SymComplexSparse(int i, int j);
    Matrix_SymComplexSparse(int i, int j, long real_nz, long imag_nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix_SymComplexSparse(int i, int j,
			    Vector<value_type, Storage0, Allocator0>& real_values,
			    Vector<long, Storage1, Allocator1>& real_ptr,
			    Vector<int, Storage2, Allocator2>& real_ind,
			    Vector<value_type, Storage0, Allocator0>& imag_values,
			    Vector<long, Storage1, Allocator1>& imag_ptr,
			    Vector<int, Storage2, Allocator2>& imag_ind);
    Matrix_SymComplexSparse(const Matrix_SymComplexSparse<T, Prop, Storage,
			    Allocator>& A);

    // Destructor.
    ~Matrix_SymComplexSparse();
    void Clear();

    // Memory management.
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    void SetData(int i, int j,
		 Vector<value_type, Storage0, Allocator0>& real_values,
		 Vector<long, Storage1, Allocator1>& real_ptr,
		 Vector<int, Storage2, Allocator2>& real_ind,
		 Vector<value_type, Storage0, Allocator0>& imag_values,
		 Vector<long, Storage1, Allocator1>& imag_ptr,
		 Vector<int, Storage2, Allocator2>& imag_ind);
    void SetData(int i, int j,
		 long real_nz, pointer real_values, long* real_ptr,
		 int* real_ind,
		 long imag_nz, pointer imag_values, long* imag_ptr,
		 int* imag_ind);
    void Nullify();
    void Reallocate(int i, int j);
    void Reallocate(int i, int j, long real_nz, long imag_nz);
    void Resize(int i, int j);
    void Resize(int i, int j, long real_nz, long imag_nz);
    void Copy(const Matrix_SymComplexSparse<T, Prop, Storage, Allocator>& A);

    // Basic methods.
    long GetNonZeros() const;
    long GetDataSize() const;
    size_t GetMemorySize() const;
    long* GetRealPtr() const;
    long* GetImagPtr() const;
    int* GetRealInd() const;
    int* GetImagInd() const;
    int GetRealPtrSize() const;
    int GetImagPtrSize() const;
    long GetRealIndSize() const;
    long GetImagIndSize() const;
    long GetRealDataSize() const;
    long GetImagDataSize() const;
    value_type* GetRealData() const;
    value_type* GetImagData() const;

    // Element acess and affectation.
    const complex<value_type> operator() (int i, int j) const;
    value_type& ValReal(int i, int j);
    const value_type& ValReal(int i, int j) const;
    value_type& ValImag(int i, int j);
    const value_type& ValImag(int i, int j) const;
    value_type& GetReal(int i, int j);
    const value_type& GetReal(int i, int j) const;
    value_type& GetImag(int i, int j);
    const value_type& GetImag(int i, int j) const;
    void Set(int i, int j, const entry_type& x);
    void AddInteraction(int i, int j, const entry_type& x);

    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<entry_type>& val,
			   bool sorted = false);
    
    Matrix_SymComplexSparse<T, Prop, Storage, Allocator>&
    operator= (const Matrix_SymComplexSparse<T, Prop, Storage, Allocator>& A);

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    void Fill(const entry_type& x);
    void FillRand();

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

    // methods used for iterative solvers
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


  //! Column-major complex sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymComplexSparse, Allocator>:
    public Matrix_SymComplexSparse<T, Prop, ColSymComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ColSymComplexSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(int i, int j);
    Matrix(int i, int j, long real_nz, long imag_nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(int i, int j,
	   Vector<value_type, Storage0, Allocator0>& real_values,
	   Vector<long, Storage1, Allocator1>& real_ptr,
	   Vector<int, Storage2, Allocator2>& real_ind,
	   Vector<value_type, Storage0, Allocator0>& imag_values,
	   Vector<long, Storage1, Allocator1>& imag_ptr,
	   Vector<int, Storage2, Allocator2>& imag_ind);
  };


  //! Row-major complex sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymComplexSparse, Allocator>:
    public Matrix_SymComplexSparse<T, Prop, RowSymComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef RowSymComplexSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(int i, int j);
    Matrix(int i, int j, long real_nz, long imag_nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(int i, int j,
	   Vector<value_type, Storage0, Allocator0>& values,
	   Vector<long, Storage1, Allocator1>& ptr,
	   Vector<int, Storage2, Allocator2>& ind,
	   Vector<value_type, Storage0, Allocator0>& imag_values,
	   Vector<long, Storage1, Allocator1>& imag_ptr,
	   Vector<int, Storage2, Allocator2>& imag_ind);
  };

} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_HXX
#endif
