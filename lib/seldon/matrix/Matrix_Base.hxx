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

#ifndef SELDON_FILE_MATRIX_BASE_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{

  
  //! Abstract base class for all matrices.
  /*!
    This class is abstract if SELDON_WITH_VIRTUAL is defined. 
    In that case, it can be used for iterative solvers or eigenvalue
    solvers such that a general algorithm can be written and compiled
    with this abstract definition of a matrix.
   */
  template<class T>
  class VirtualMatrix
  {
    // Attributes.
  protected:
    // Number of rows.
    int m_;
    // Number of columns.
    int n_;

    // Methods.
  public:
    // Constructors.
    VirtualMatrix();
    explicit VirtualMatrix(int i, int j);

    // Basic methods.
    int GetM() const;
    int GetN() const;
    int GetM(const Seldon::SeldonTranspose& status) const;
    int GetN(const Seldon::SeldonTranspose& status) const;
    int GetM(const CBLAS_TRANSPOSE& status) const;
    int GetN(const CBLAS_TRANSPOSE& status) const;
    long GetSize() const;

#ifdef SELDON_WITH_VIRTUAL
    virtual ~VirtualMatrix();
    
  protected:
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

    public:
    // basic manipulation of matrix
    virtual void Reallocate(int, int);
    virtual void Zero();
    virtual void SetEntry(int, int, const T&);
    virtual void AddInteraction(int, int, const T&);
    virtual void AddInteractionRow(int, int, const Vector<int>&,
				   const Vector<T>& val, bool sorted = false);

    virtual void AddInteractionColumn(int, int, const Vector<int>&,
				      const Vector<T>& val,
				      bool sorted = false);
    
    virtual void AddDistantInteraction(int i, int jglob, int proc,
				       const T& val);
    
    virtual void AddRowDistantInteraction(int iglob, int j, int proc,
					  const T& val);

    virtual size_t GetMemorySize() const;
    virtual void Clear();
    virtual void ClearRow(int i);

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

    // methods for eigensolvers
    virtual bool IsSymmetric() const;
    bool IsComplex() const;
#endif
    
  };

  
  //! Base class for all matrices.
  /*!
    It stores some data and matrix dimensions. It defines basic
    methods as well.
  */
  template <class T, class Allocator 
	    = typename SeldonDefaultAllocator<VectFull, T>::allocator>
  class Matrix_Base : public VirtualMatrix<T>
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
    // Pointer to stored elements.
    pointer data_;

    // Methods.
  public:
    // Constructors.
    Matrix_Base();
    explicit Matrix_Base(int i, int j);
    Matrix_Base(const Matrix_Base<T, Allocator>& A);

    // Destructor.
    ~Matrix_Base();

    pointer GetData() const;
    const_pointer GetDataConst() const;
    void* GetDataVoid() const;
    const void* GetDataConstVoid() const;

  };


  template <class T, class Prop, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Matrix<T, Prop, Storage, Allocator>& A);


} // namespace Seldon.

#define SELDON_FILE_MATRIX_BASE_HXX
#endif
