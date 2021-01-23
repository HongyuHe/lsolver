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

#ifndef SELDON_FILE_DISTRIBUTED_VECTOR_HXX

namespace Seldon
{
  
  //! class storing a vector distributed over all the processors
  /*!
    This class is useful when some rows of the vector
    are shared by several processors. Functions DotProd, DotProdConj
    and Norm2 are overloaded in order to take into account this overlap.
   */
  template<class T, class Allocator
	   = typename SeldonDefaultAllocator<VectFull, T>::allocator>
  class DistributedVector : public Vector<T, Vect_Full, Allocator>
  {
  protected :
    //! row numbers shared with other processors
    /*!
      In this array, you store all the row numbers
      which are already counted in another processor.
      For example :
      if processor 0 contains rows (0, 3, 5, 7, 8)
      if processor 1 contains rows (0, 1, 2, 4, 5, 6)
      Then you can set OverlapRowNumbers empty on processor 0,
      and equal to (0, 5) on processor 1
     */
    const IVect& OverlapRowNumbers;    
    //! MPI communicator grouping processors involved in the computation
    const MPI_Comm& comm_;
        
  public :
    // constructors
    DistributedVector(const IVect& rows, const MPI_Comm& comm);
    DistributedVector(const DistributedVector<T, Allocator>& V);
    
    // memory management
    void SetData(Vector<T, Vect_Full, Allocator>& x);
    void SetData(int n, T* data);
    
    // basic functions
    DistributedVector<T, Allocator>&
    operator=(const DistributedVector<T, Allocator>& X);
    
    int GetNbOverlap() const;
    int GetOverlapRow(int i) const;
    const MPI_Comm& GetCommunicator() const;
    
  };
  
  
  // returns X.Y
  template<class T1, class Allocator1>
  T1 DotProdVector(const DistributedVector<T1, Allocator1>& X,
		   const DistributedVector<T1, Allocator1>& Y);
  
  // returns X' . Y
  template<class T1, class Allocator1>
  T1 DotProdConjVector(const DistributedVector<T1, Allocator1>& X,
		       const DistributedVector<T1, Allocator1>& Y);  
  
  // returns euclidian norm of x
  template<class T, class Allocator>
  T Norm2(const DistributedVector<complex<T>, Allocator>& x);
  
  // returns euclidian norm of x
  template<class T, class Allocator>
  T Norm2(const DistributedVector<T, Allocator>& x);

  template<class T>
  T minComplex(const T& x, const T& y);

  template<class T>
  complex<T> minComplex(const complex<T>& x, const complex<T>& y);

  template<class T>
  T maxComplex(const T& x, const T& y);

  template<class T>
  complex<T> maxComplex(const complex<T>& x, const complex<T>& y);

  void ExtractDistributedSharedNumbers(const IVect& MatchingNumber_Subdomain,
				       const Vector<IVect>& MatchingDofOrig_Subdomain,
				       const Vector<bool>& is_local, int nb_ddl_local, 
				       IVect& MatchingNumber_Local,
				       Vector<IVect>& MatchingDof_Local);
  
  void AssembleVectorMin(Vector<int>& X, Vector<int>& Xproc,
                         const IVect& ProcNumber,
			 const Vector<IVect>& DofNumber,
                         const MPI_Comm& comm, int Nvol, int nb_u, int tag);

  template<class T>
  void AssembleVector(Vector<T>& X, const MPI_Op& oper,
                      const IVect& ProcNumber, const Vector<IVect>& DofNumber,
                      const MPI_Comm& comm, int Nvol, int nb_u, int tag);

  template<class T>
  void ExchangeVector(Vector<T>& X,
                      const IVect& ProcNumber, const Vector<IVect>& DofNumber,
                      const MPI_Comm& comm, int Nvol, int nb_u, int tag);

  template<class T, class Treal>
  void ExchangeRelaxVector(Vector<T>& X, const Treal& omega, int proc,
			   const IVect& ProcNumber,
			   const Vector<IVect>& DofNumber,
			   const MPI_Comm& comm, int Nvol, int nb_u, int tag);

  template<class T>
  T DotProd(const DistributedVector<T>& X, const DistributedVector<T>& Y);

  template<class T>
  complex<T> DotProd(const DistributedVector<complex<T> >& X,
                     const DistributedVector<T>& Y);

  template<class T>
  complex<T> DotProd(const DistributedVector<T>& X,
                     const DistributedVector<complex<T> >& Y);

  template<class T>
  T DotProdConj(const DistributedVector<T>& X, const DistributedVector<T>& Y);

  template<class T>
  complex<T> DotProdConj(const DistributedVector<complex<T> >& X,
			 const DistributedVector<complex<T> >& Y);

  template<class T>
  complex<T> DotProdConj(const DistributedVector<complex<T> >& X,
			 const DistributedVector<T>& Y);

  template<class T>
  complex<T> DotProdConj(const DistributedVector<T>& X,
			 const DistributedVector<complex<T> >& Y);
  
}

#define SELDON_FILE_DISTRIBUTED_VECTOR_HXX
#endif

