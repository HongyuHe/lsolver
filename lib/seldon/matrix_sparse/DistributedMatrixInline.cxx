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

#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_INLINE_CXX

#include "DistributedMatrix.hxx"

namespace Seldon
{

  /**************************
   * DistributedMatrix_Base *
   **************************/

  
  //! returns MPI communicator (processors that will share the matrix)
  template<class T>
  inline MPI_Comm& DistributedMatrix_Base<T>::GetCommunicator()
  {
    return comm_;
  }
  
  
  //! returns MPI communicator (processors that will share the matrix)
  template<class T>
  inline const MPI_Comm& DistributedMatrix_Base<T>::GetCommunicator() const
  {
    return comm_;
  }
  

  //! returns the local number of rows
  template<class T>
  inline int DistributedMatrix_Base<T>::GetLocalM() const
  {
    return dist_col.GetM();
  }


  //! returns the local number of columns
  template<class T>
  inline int DistributedMatrix_Base<T>::GetLocalN() const
  {
    return dist_row.GetM();
  }

    
  //! returns the global number of rows
  template<class T>
  inline int DistributedMatrix_Base<T>::GetGlobalM() const
  {
    return nglob_;
  }
  
  
  //! returns the number of scalar unknowns
  template<class T>
  inline int DistributedMatrix_Base<T>::GetNodlScalar() const
  {
    return nodl_scalar_;
  }
  

  //! returns the number of scalar unknowns
  template<class T>
  inline int DistributedMatrix_Base<T>::GetNbScalarUnknowns() const
  {
    return nb_unknowns_scal_;
  }

  
  //! adding a non-zero entry, between a local column
  //! and a non-local column
  /*!
    We assume here that the row number is local to the current processor
    \param[in] i local row number
    \param[in] jglob global column number
    \param[in] proc2 distant processor containing the global column
    \param[in] val value of the non-zero entry
  */
  template<class T>
  inline void DistributedMatrix_Base<T>::
  AddDistantInteraction(int i, int jglob, int proc2, const T& val)
  {
    if (local_number_distant_values)
      SwitchToGlobalNumbers();
    
    AddDistantValue(dist_col(i), proc_col(i), jglob, proc2, val);
  }
  
  
  //! adding a non-zero entry, between a local row and a non-local row
  /*!
    We assume here that the column number is local to the current processor
    \param[in] iglob global row number
    \param[in] j local column number
    \param[in] proc2 distant processor containing the global row
    \param[in] val value of the non-zero entry
  */
  template<class T>
  inline void DistributedMatrix_Base<T>
  ::AddRowDistantInteraction(int iglob, int j,
			     int proc2, const T& val)
  {
    if (local_number_distant_values)
      SwitchToGlobalNumbers();
    
    AddDistantValue(dist_row(j), proc_row(j), iglob, proc2, val);
  }
  
  
  //! returns the maximum number of values in dist_col to exchange
  template<class T>
  inline long DistributedMatrix_Base<T>
  ::GetMaxDataSizeDistantCol() const
  {
    return size_max_distant_col;
  }
  
  
  //! returns the maximum number of values in dist_row to exchange
  template<class T>
  inline long DistributedMatrix_Base<T>
  ::GetMaxDataSizeDistantRow() const
  { 
    return size_max_distant_row;
  }
  
  
  //! returns true if the matrix is ready to perform a matrix-vector product
  template<class T>
  inline bool DistributedMatrix_Base<T>
  ::IsReadyForMltAdd() const
  {
    return local_number_distant_values;
  }


  //! tells that there are no distant rows (mltadd ready)
  template<class T>
  inline void DistributedMatrix_Base<T>::SetReadyForMltAdd(bool all_zero)
  {
    local_number_distant_values = all_zero;
  }


  //! returns the number of distant non-zero entries for the local row i
  template<class T>
  inline int DistributedMatrix_Base<T>::GetDistantColSize(int i) const
  {
    return dist_col(i).GetM();
  }


  //! returns the global column number of distant non-zero entry j for the local row i
  template<class T>
  inline int DistributedMatrix_Base<T>::IndexGlobalCol(int i, int j) const
  {
    if (this->local_number_distant_values)      
      return global_col_to_recv(dist_col(i).Index(j));  

    return dist_col(i).Index(j);
  }


  //! returns the processor associated with distant non-zero entry j for the local row i
  template<class T>
  inline int DistributedMatrix_Base<T>::ProcessorDistantCol(int i, int j) const
  {
    return proc_col(i)(j);  
  }


  //! returns the value of distant non-zero entry j for the local row i
  template<class T>
  inline const T& DistributedMatrix_Base<T>::ValueDistantCol(int i, int j) const
  {
    return dist_col(i).Value(j);  
  }


  //! returns the number of distant non-zero entries for the local column i
  template<class T>
  inline int DistributedMatrix_Base<T>::GetDistantRowSize(int i) const
  {
    return dist_row(i).GetM();
  }


  //! returns the global row number of distant non-zero entry j for the local column i
  template<class T>
  inline int DistributedMatrix_Base<T>::IndexGlobalRow(int i, int j) const
  {
    if (this->local_number_distant_values)
      return global_row_to_recv(dist_row(i).Index(j));  
    
    return dist_row(i).Index(j);
  }


  //! returns the processor associated with distant non-zero entry j for the local column i
  template<class T>
  inline int DistributedMatrix_Base<T>::ProcessorDistantRow(int i, int j) const
  {
    return proc_row(i)(j);  
  }


  //! returns the value of distant non-zero entry j for the local column i
  template<class T>
  inline const T& DistributedMatrix_Base<T>::ValueDistantRow(int i, int j) const
  {
    return dist_row(i).Value(j);  
  }
  

  /*********************
   * DistributedMatrix *
   *********************/


  //! Adds val for local row i and global column jglob (proc owns this column)
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AddDistantInteraction(int i, int jglob, int proc, const T& val)
  {
    DistributedMatrix_Base<T>::AddDistantInteraction(i, jglob, proc, val);
  }

    
  //! Adds val for global row iglob and local column j (proc owns the row)
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AddRowDistantInteraction(int iglob, int j, int proc, const T& val)
  {
    DistributedMatrix_Base<T>::AddRowDistantInteraction(iglob, j, proc, val);
  }
  

#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, *this, x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ApplySor(const SeldonTranspose& trans, Vector<Tcplx>& x, const Vector<Tcplx>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, *this, x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha, trans, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, trans, *this, x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(*this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(*this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(trans, *this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans, *this, x, y);
  }
#endif


  //! converts a distributed matrix to another one
  template<class T, class Prop1, class Storage1, class Allocator1,
           class Prop2, class Storage2, class Allocator2>
  inline void Copy(const DistributedMatrix<complex<T>, Prop1, Storage1, Allocator1>& A,
		   DistributedMatrix<T, Prop2, Storage2, Allocator2>& B)
  {
    throw WrongArgument("Copy", "incompatible types");    
  }


  /****************************************
   * Mlt, MltAdd for distributed matrices *
   ****************************************/
  

  template<class T, class Prop, class Storage, class Allocator>
  inline void Mlt(const T& alpha,
                  DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }

  template<class T, class Prop, class Storage, class Allocator>
  inline void Mlt(const T& alpha,
                  DistributedMatrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }


  //! Performs the multiplication of a matrix with a vector.
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  inline void MltVector(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& M,
			const Vector<T1, Storage1, Allocator1>& X,
			Vector<T2, Storage2, Allocator2>& Y, bool assemble)
  {
    T2 one, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAddVector(one, M, X, zero, Y, assemble);
  }

  
  //! Performs the multiplication of a matrix with a vector.
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  inline void MltVector(const T3& alpha,
			const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
			const Vector<T2, Storage2, Allocator2>& X,
			Vector<T3, Storage3, Allocator3>& Y, bool assemble)
  {
    T3 zero;
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAddVector(alpha, M, X, zero, Y, assemble);
  }

  
  //! Performs the multiplication of a matrix with a vector.
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  inline void MltVector(const SeldonTranspose& Trans,
			const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
			const Vector<T2, Storage2, Allocator2>& X,
			Vector<T3, Storage3, Allocator3>& Y, bool assemble)
  {
    T3 one, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAddVector(one, Trans, M, X, zero, Y, assemble);
  }


  //! multiplication by a scalar
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1>
  inline void MltScalar(const T0& alpha,
			DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    A *= alpha;
  }
  
  
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
		  DistributedMatrix<T, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const complex<T>& alpha,
		  const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
		  DistributedMatrix<complex<T>, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }


  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const DistributedMatrix<complex<T>, Prop1, Storage1, Allocator1>& A,
		  DistributedMatrix<complex<T>, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const DistributedMatrix<complex<T>, Prop1, Storage1, Allocator1>& A,
		  DistributedMatrix<T, Prop2, Storage2, Allocator2>& B)
  {
    throw WrongArgument("Add", "incompatible types");    
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    MltVector(A, X, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y,
		  bool assemble)
  {
    MltVector(A, X, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    MltVector(trans, A, X, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y,
		  bool assemble)
  {
    MltVector(trans, A, X, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAddVector(alpha, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const complex<T0>& alpha,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAddVector(alpha, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltAddComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAddVector(alpha, trans, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAddVector(alpha, trans, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltAddComplex", "Incompatible matrix-vector product");			
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const T& omega, int iter, int type_ssor)
  {
    SorVector(M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const SeldonTranspose& trans,
		  const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const T& omega, int iter, int type_ssor)
  {
    SorVector(trans, M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<complex<T>, Storage2, Allocator2>& Y,
		  const Vector<complex<T>, Storage1, Allocator1>& X,
		  const T& omega, int iter, int type_ssor)
  {
    SorVector(M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const SeldonTranspose& trans,
		  const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<complex<T>, Storage2, Allocator2>& Y,
		  const Vector<complex<T>, Storage1, Allocator1>& X,
		  const T& omega, int iter, int type_ssor)
  {
    SorVector(trans, M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const DistributedMatrix<complex<T>, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const T& omega, int iter, int type_ssor)
  {
    throw WrongArgument("SorVector", "Incompatible matrix-vector product");
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const SeldonTranspose&,
		  const DistributedMatrix<complex<T>, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const T& omega, int iter, int type_ssor)
  {
    throw WrongArgument("SorVector", "Incompatible matrix-vector product");
  }

}

#define SELDON_FILE_DISTRIBUTED_MATRIX_INLINE_CXX
#endif

