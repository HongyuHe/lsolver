// Copyright (C) 2014 INRIA
// Author(s): Marc Duruflé
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

#ifndef SELDON_FILE_DISTRIBUTED_VECTOR_INLINE_CXX

#include "DistributedVector.hxx"

namespace Seldon
{
  
  //! Constructor taking overlapped rows
  /*!
    In the array rows, you store all the row numbers
    which are already counted in another processor.
    For example :
    if processor 0 contains rows (0, 3, 5, 7, 8)
    if processor 1 contains rows (0, 1, 2, 4, 5, 6)
    Then you can set OverlapRowNumbers empty on processor 0,
    and equal to (0, 5) on processor 1.
    If there is no row shared by processors, the array rows
    will be empty for each processor
  */
  template<class T, class Allocator>
  inline DistributedVector<T, Allocator>
  ::DistributedVector(const IVect& rows, const MPI_Comm& comm)
    : OverlapRowNumbers(rows), comm_(comm)
  {
  }
  
  
  //! Copy constructor
  template<class T, class Allocator>
  inline DistributedVector<T, Allocator>::
  DistributedVector(const DistributedVector<T, Allocator>& V)
    : Vector<T, VectFull, Allocator>(V),
      OverlapRowNumbers(V.OverlapRowNumbers), comm_(V.comm_)
  {
  }
  
  
  //! setting pointer
  template<class T, class Allocator>
  inline void DistributedVector<T, Allocator>::
  SetData(Vector<T, Vect_Full, Allocator>& x)
  {
    Vector<T, Vect_Full, Allocator>::SetData(x.GetM(), x.GetData());
    x.Nullify();
  }
  
  
  //! setting pointer
  template<class T, class Allocator>
  inline void DistributedVector<T, Allocator>::SetData(int n, T* data)
  {
    Vector<T, Vect_Full, Allocator>::SetData(n, data);
  }
  
  
  //! operator =
  template<class T, class Allocator>
  inline DistributedVector<T, Allocator>& DistributedVector<T, Allocator>::
  operator=(const DistributedVector<T, Allocator>& X)
  {
    Vector<T>::Copy(static_cast<const Vector<T, Vect_Full, Allocator>& >(X));
    return *this;
  }
  
  
  //! returns the number of rows already counted 
  template<class T, class Allocator>
  inline int DistributedVector<T, Allocator>::GetNbOverlap() const
  {
    return OverlapRowNumbers.GetM();
  }
  
  
  //! returns an overlapped row number
  template<class T, class Allocator>
  inline int DistributedVector<T, Allocator>::GetOverlapRow(int i) const
  {
    return OverlapRowNumbers(i);
  }
  
  
  //! returns communicator
  template<class T, class Allocator>
  inline const MPI_Comm& DistributedVector<T, Allocator>
  ::GetCommunicator() const
  {
    return comm_;
  }


  template<class T>
  inline T DotProd(const DistributedVector<T>& X,
		   const DistributedVector<T>& Y)
  {
    return DotProdVector(X, Y);
  }

  template<class T>
  inline complex<T> DotProd(const DistributedVector<complex<T> >& X,
                            const DistributedVector<T>& Y)
  {
    abort();
    return complex<T>(0, 0);
  }

  template<class T>
  inline complex<T> DotProd(const DistributedVector<T>& X,
                            const DistributedVector<complex<T> >& Y)
  {
    abort();
    return complex<T>(0, 0);
  }

  template<class T>
  inline T DotProdConj(const DistributedVector<T>& X,
		       const DistributedVector<T>& Y)
  {
    return DotProdVector(X, Y);
  }

  template<class T>
  inline complex<T> DotProdConj(const DistributedVector<complex<T> >& X,
				const DistributedVector<complex<T> >& Y)
  {
    return DotProdConjVector(X, Y);
  }

  template<class T>
  inline complex<T> DotProdConj(const DistributedVector<complex<T> >& X,
				const DistributedVector<T>& Y)
  {
    abort();
    return complex<T>(0, 0);
  }

  template<class T>
  inline complex<T> DotProdConj(const DistributedVector<T>& X,
				const DistributedVector<complex<T> >& Y)
  {
    abort();
    return complex<T>(0, 0);
  }
  
}

#define SELDON_FILE_DISTRIBUTED_VECTOR_INLINE_CXX
#endif

