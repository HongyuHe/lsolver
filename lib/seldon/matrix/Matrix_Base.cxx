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


#ifndef SELDON_FILE_MATRIX_BASE_CXX

#include "Matrix_Base.hxx"

namespace Seldon
{

#ifdef SELDON_WITH_VIRTUAL
  //! reallocates the matrix
  template <class T>
  void VirtualMatrix<T>::Reallocate(int m, int n)
  {
    // this method should be overloaded
    // if not overloaded, an exception is raised:
    throw Undefined("Reallocate", "Not implemented"); 
  }


  //! sets the stored entries to zero
  template <class T>
  void VirtualMatrix<T>::Zero()
  {
    // this method should be overloaded
    // if not overloaded, an exception is raised:
    throw Undefined("Zero", "Not implemented"); 
  }


  //! Sets x to element (i, j) of the matrix
  template <class T>
  void VirtualMatrix<T>::SetEntry(int i, int j, const T& x)
  {
    // this method should be overloaded in sparse matrices
    // if not overloaded, an exception is raised:
    throw Undefined("SetEntry", "Not implemented");
  }
  
  
  //! Adds x to element (i, j) of the matrix
  template <class T>
  void VirtualMatrix<T>::AddInteraction(int i, int j, const T& x)
  {
    // this method should be overloaded in sparse matrices
    // if not overloaded, an exception is raised:
    throw Undefined("AddInteraction", "Not implemented");
  }

  
  //! Adds values to several non-zero entries of a sparse matrix.
  template <class T>
  void VirtualMatrix<T>::AddInteractionRow(int, int, const Vector<int>&,
					   const Vector<T>& val, bool sorted)
  {
    // this method should be overloaded in sparse matrices
    // if not overloaded, an exception is raised:
    throw Undefined("AddInteractionRow", "Not implemented");
  }


  //! Adds values to several non-zero entries of a sparse matrix.
  template <class T>
  void VirtualMatrix<T>::AddInteractionColumn(int, int, const Vector<int>&,
					      const Vector<T>& val, bool sorted)
  {
    // this method should be overloaded in sparse matrices
    // if not overloaded, an exception is raised:
    throw Undefined("AddInteractionColumn", "Not implemented");
  }


  //! Adds a distant value to the matrix
  /*!
    Here the row is local, and the column is global (associated with 
    another processor), this method is overloaded in DistributedMatrix
    \param[in] i local row number
    \param[in] jglob global column number
    \param[in] proc processor associated with the column jglob
    \param[in] val value to add in entry A(i, j)
   */
  template <class T>
  void VirtualMatrix<T>::AddDistantInteraction(int i, int jglob, int proc,
					       const T& val)
  {
    // this method should be overloaded in distributed sparse matrices
    // if not overloaded, an exception is raised:
    throw Undefined("AddDistantInteraction", "Not implemented");
  }


  //! Adds a distant value to the matrix
  /*!
    Here the column is local, and the row is global (associated with 
    another processor), this method is overloaded in DistributedMatrix
    \param[in] iglob global row number
    \param[in] j local column number
    \param[in] proc processor associated with the column jglob
    \param[in] val value to add in entry A(i, j)
   */    
  template <class T>
  void VirtualMatrix<T>::AddRowDistantInteraction(int iglob, int j, int proc,
						  const T& val)
  {
    // this method should be overloaded in distributed sparse matrices
    // if not overloaded, an exception is raised:
    throw Undefined("AddRowDistantInteraction", "Not implemented");
  }
  

  //! returns the memory used by the object in bytes
  template <class T>
  size_t VirtualMatrix<T>::GetMemorySize() const
  {
    // this method should be overloaded
    // if not overloaded, an exception is raised:
    throw Undefined("GetMemorySize", "Not implemented"); 
  }
  
  
  //! Clears the matrix (the memory is released)
  template <class T>
  void VirtualMatrix<T>::Clear()
  {
    // this method should be overloaded
    // if not overloaded, an exception is raised:
    throw Undefined("Clear", "Not implemented"); 
  }

  
  //! Clears row i
  template <class T>
  inline void VirtualMatrix<T>::ClearRow(int i)
  {
    // this method should be overloaded in sparse matrices
    // if not overloaded, an exception is raised:
    throw Undefined("ClearRow", "Not implemented");
  }


  //! Applies S.O.R method to solve A x = r.
  template<class T> void VirtualMatrix<T>
  ::ApplySor(const SeldonTranspose&, Vector<Treal>& x, const Vector<Treal>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    throw Undefined("ApplySOR", "Not implemented");
  }
  

  //! Applies S.O.R method to solve A x = r or A^T x = r.
  template<class T> void VirtualMatrix<T>
  ::ApplySor(const SeldonTranspose&, Vector<Tcplx>& x, const Vector<Tcplx>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    throw Undefined("ApplySOR with transpose", "Not implemented");
  }


  //! Computes y = beta + alpha A x.
  template<class T> void VirtualMatrix<T>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = beta + alpha A x.
  template<class T> void VirtualMatrix<T>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = beta + alpha A x.
  template<class T> void VirtualMatrix<T>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose&,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = beta + alpha A x.
  template<class T> void VirtualMatrix<T>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose&,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = A x
  template<class T> void VirtualMatrix<T>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    throw Undefined("MltVector", "Not implemented");
  }


  //! Computes y = A x
  template<class T> void VirtualMatrix<T>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    throw Undefined("MltVector", "Not implemented");
  }
  

  //! Computes y = A^T x
  template<class T> void VirtualMatrix<T>
  ::MltVector(const SeldonTranspose&,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    throw Undefined("MltVector with transpose", "Not implemented");
  }


  //! Computes y = A^T x
  template<class T> void VirtualMatrix<T>
  ::MltVector(const SeldonTranspose&,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    throw Undefined("MltVector with transpose", "Not implemented");
  }
  
  
  //! returns true if the matrix is symmetric
  template<class T> bool VirtualMatrix<T>::IsSymmetric() const
  {
    return false;
  }
  
  
  //! returns true if the matrix is complex
  template<class T> bool VirtualMatrix<T>::IsComplex() const
  {
    return IsComplexNumber(T());
  }
#endif


  // operator<< overloaded for matrices.
  /*!
    \param out output stream.
    \param A matrix to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Matrix<T, Prop, Storage, Allocator>& A)
  {
    
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetN(); j++)
	  out << A(i, j) << '\t';
	out << endl;
      }

    return out;

  }

} // namespace Seldon.

#define SELDON_FILE_MATRIX_BASE_CXX
#endif
