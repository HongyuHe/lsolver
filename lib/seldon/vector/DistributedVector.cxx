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

#ifndef SELDON_FILE_DISTRIBUTED_VECTOR_CXX

#include "DistributedVector.hxx"

namespace Seldon
{
  
  //! computation of scalar product X.Y
  /*!
    Overlapped rows are correctly handled so that
    the result is the scalar product of the global vectors.
    The result should be the same on each processor.
   */
  template<class T1, class Allocator1>
  T1 DotProdVector(const DistributedVector<T1, Allocator1>& X,
		   const DistributedVector<T1, Allocator1>& Y)
  {
    T1 value;
    SetComplexZero(value);
    for (int i = 0; i < X.GetNbOverlap(); i++)
      value += X(X.GetOverlapRow(i)) * Y(Y.GetOverlapRow(i));
    
    value = -value + DotProd(static_cast<const Vector<T1, VectFull,
			     Allocator1>& >(X),
			     static_cast<const Vector<T1, VectFull,
			     Allocator1>& >(Y));
    
    const MPI_Comm& comm = X.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc > 1)
      {
	T1 sum; SetComplexZero(sum);
	Vector<int64_t> xtmp;
        MpiAllreduce(comm, &value, xtmp, &sum, 1, MPI_SUM);
	return sum;
      }
    
    return value;
  }
  
  
  //! computation of scalar product X'.Y
  /*!
    Overlapped rows are correctly handled so that
    the result is the scalar product of the global vectors.
    The result should be the same on each processor.
   */
  template<class T1, class Allocator1>
  T1 DotProdConjVector(const DistributedVector<T1, Allocator1>& X,
		       const DistributedVector<T1, Allocator1>& Y)
  {
    T1 value;
    SetComplexZero(value);
    for (int i = 0; i < X.GetNbOverlap(); i++)
      value += conjugate(X(X.GetOverlapRow(i))) * Y(Y.GetOverlapRow(i));
    
    value = -value + 
      DotProdConj(static_cast<const Vector<T1, VectFull, Allocator1>& >(X),
		  static_cast<const Vector<T1, VectFull, Allocator1>& >(Y));
    
    const MPI_Comm& comm = X.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc > 1)
      {
	T1 sum; SetComplexZero(sum);
	Vector<int64_t> xtmp;
        MpiAllreduce(comm, &value, xtmp, &sum, 1, MPI_SUM);
	return sum;
      }
    
    return value;
  }
  
  
  //! returns euclidian norm of vector x
  template<class T, class Allocator>
  T Norm2(const DistributedVector<complex<T>, Allocator>& x)
  {
    T scal = abs(DotProdConj(x, x));
    return sqrt(scal);
  }
  
  
  //! returns euclidian norm of vector x
  template<class T, class Allocator>
  T Norm2(const DistributedVector<T, Allocator>& x)
  {
    T scal = DotProd(x, x);
    return sqrt(scal);
  }
  
  
  //! minimum for real numbers
  template<class T>
  T minComplex(const T& x, const T& y)
  {
    return min(x, y);
  }


  //! minimum for complex numbers
  template<class T>
  complex<T> minComplex(const complex<T>& x, const complex<T>& y)
  {
    if (abs(x) > abs(y))
      return y;
    
    return x;
  }


  //! maximum for real numbers
  template<class T>
  T maxComplex(const T& x, const T& y)
  {
    return max(x, y);
  }


  //! maximum for complex numbers
  template<class T>
  complex<T> maxComplex(const complex<T>& x, const complex<T>& y)
  {
    if (abs(x) > abs(y))
      return x;
    
    return y;
  }
  

  //! extracts MatchingNumber_Local / MatchingDof_Local from arrays
  //! MatchingNumber_Subdomain / MatchingDofOrig_Subdomain
  void ExtractDistributedSharedNumbers(const IVect& MatchingNumber_Subdomain,
				       const Vector<IVect>& MatchingDofOrig_Subdomain,
				       const Vector<bool>& is_local, int nb_ddl_local, 
				       IVect& MatchingNumber_Local,
				       Vector<IVect>& MatchingDof_Local)
  {
    MatchingNumber_Local.Clear();
    MatchingDof_Local.Clear();
    for (int i = 0; i < MatchingNumber_Subdomain.GetM(); i++)
      {
	int proc = MatchingNumber_Subdomain(i);
	int nb_dof = 0;
	for (int j = 0; j < MatchingDofOrig_Subdomain(i).GetM(); j++)
	  if ((nb_ddl_local > 0) &&
	      (is_local(MatchingDofOrig_Subdomain(i)(j))))
	    nb_dof++;
	
	if (nb_dof > 0)
	  {
	    IVect dof_num(nb_dof);
	    nb_dof = 0;
	    for (int j = 0; j < MatchingDofOrig_Subdomain(i).GetM(); j++)
	      if (is_local(MatchingDofOrig_Subdomain(i)(j)))
		{
		  dof_num(nb_dof) = MatchingDofOrig_Subdomain(i)(j);
		  nb_dof++;
		}
	    
	    MatchingNumber_Local.PushBack(proc);
	    MatchingDof_Local.PushBack(dof_num);	    
	  }
      }
  }
  
  
  //! assembles minimums of two vectors
  void AssembleVectorMin(Vector<int>& X, Vector<int>& Xproc,
                         const IVect& ProcNumber,
			 const Vector<IVect>& DofNumber,
                         const MPI_Comm& comm, int Nvol, int nb_u, int tag)
  {
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    // number of procs interacting with the current processor
    int nb_dom = DofNumber.GetM();
    Vector<MPI_Request> request_send(nb_dom), request_recv(nb_dom);
    MPI_Status status;
    Vector<Vector<int> > xsend(nb_dom);
    Vector<Vector<int64_t> > xsend_tmp(nb_dom);
    
    // sending informations to other processors
    for (int i = 0; i < nb_dom; i++)
      {
	int j = ProcNumber(i);
	int nb = DofNumber(i).GetM();
        if (nb > 0)
          {
            xsend(i).Reallocate(2*nb_u*nb);
            for (int m = 0; m < nb_u; m++)
              for (int k = 0; k < nb; k++)
                {
                  xsend(i)(nb_u*k+m) = X(DofNumber(i)(k) + m*Nvol);
                  xsend(i)(nb_u*k+m+nb*nb_u)
		    = Xproc(DofNumber(i)(k) + m*Nvol);
                }
            
            // sending the value to the corresponding processor
            MPI_Isend(&xsend(i)(0), 2*nb*nb_u,
                      MPI_INTEGER, j, tag, comm, &request_send(i));
          }
      }
    
    // receiving the informations
    Vector<Vector<int> > xdom(nb_dom);
    for (int i = 0; i < nb_dom; i++)
      {
	int j = ProcNumber(i);
	int nb = DofNumber(i).GetM();
	if (nb > 0)
          {
            xdom(i).Reallocate(2*nb_u*nb);
            xdom(i).Fill(0);
            
            // receiving the values of domain j
            MPI_Irecv(&xdom(i)(0), 2*nb*nb_u,
                      MPI_INTEGER, j, tag, comm, &request_recv(i));
          }
      }
    
    // now waiting all communications are effective
    for (int i = 0; i < nb_dom; i++)
      if (DofNumber(i).GetM() > 0)
        MPI_Wait(&request_recv(i), &status);
    
    for (int i = 0; i < nb_dom; i++)
      if (DofNumber(i).GetM() > 0)
        MPI_Wait(&request_send(i), &status);
    
    // final reduction step
    for (int i = 0; i < nb_dom; i++)
      {
        for (int m = 0; m < nb_u; m++)
          for (int k = 0; k < DofNumber(i).GetM(); k++)
            {
              int nb = DofNumber(i).GetM();
              int p = DofNumber(i)(k) + m*Nvol;
              int proc = xdom(i)(k*nb_u+m+nb*nb_u);
              int col = xdom(i)(k*nb_u+m);
              if (proc < Xproc(p))
                {
                  Xproc(p) = proc;
                  X(p) = col;
                }
              else if (proc == Xproc(p))
                {
                  if (col < X(p))
                    X(p) = col;
                }
            }
      }
  }

  
  //! assembles a distributed vector (some rows are shared between processors)
  /*!
    \param[inout] X vector to assemble
    \param[in] oper reduction operator (MPI_SUM, MPI_MIN or MPI_MAX)
    \param[in] ProcNumber processors interacting with the current processor
    \param[in] DofNumber for each processor, rows to send/receive
    \param[in] comm MPI communicator
    \param[in] Nvol DofNumber will be included between <0, Nvol-1>
    \param[in] nb_u number of unknowns
    (X is supposed to be of size >= Nvol*nb_u)
    \param[in] tag MPI tag
    For each unknown, it is assumed that the same row numbers are affected,
    that's why you can provide DofNumber for only the first unknown
    and consider that rows of other unknowns are obtained
    with the operation m*Nvol + i, m in <0, nb_u-1> and
    m in <0, Nvol-1>
   */
  template<class T>
  void AssembleVector(Vector<T>& X, const MPI_Op& oper,
                      const IVect& ProcNumber, const Vector<IVect>& DofNumber,
                      const MPI_Comm& comm, int Nvol, int nb_u, int tag)
  {
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;

    // number of procs interacting with the current processor
    int nb_dom = DofNumber.GetM();
    Vector<MPI_Request> request_send(nb_dom), request_recv(nb_dom);
    MPI_Status status;
    Vector<Vector<T> > xsend(nb_dom);
    Vector<Vector<int64_t> > xsend_tmp(nb_dom);
    
    // sending informations to other processors
    for (int i = 0; i < nb_dom; i++)
      {
	int j = ProcNumber(i);
	int nb = DofNumber(i).GetM();
        if (nb > 0)
          {
            xsend(i).Reallocate(nb_u*nb);
            for (int m = 0; m < nb_u; m++)
              for (int k = 0; k < nb; k++)
                xsend(i)(nb_u*k+m) = X(DofNumber(i)(k) + m*Nvol);

            // sending the value to the corresponding processor
            request_send(i) = MpiIsend(comm, xsend(i), xsend_tmp(i),
				       nb*nb_u, j, tag);
          }
      }
    
    // receiving the informations
    T zero; SetComplexZero(zero);
    Vector<Vector<T> > xdom(nb_dom);
    Vector<Vector<int64_t> > xdom_tmp(nb_dom);
    for (int i = 0; i < nb_dom; i++)
      {
	int j = ProcNumber(i);
	int nb = DofNumber(i).GetM();
	if (nb > 0)
          {
            xdom(i).Reallocate(nb_u*nb);
            xdom(i).Fill(zero);

            // receiving the values of domain j
            request_recv(i) = MpiIrecv(comm, xdom(i), xdom_tmp(i),
				       nb*nb_u, j, tag);
          }
      }
    
    // now waiting all communications are effective
    for (int i = 0; i < nb_dom; i++)
      if (DofNumber(i).GetM() > 0)
        MPI_Wait(&request_send(i), &status);
    
    for (int i = 0; i < nb_dom; i++)
      if (DofNumber(i).GetM() > 0)
        MPI_Wait(&request_recv(i), &status);
    
    // completing writing operations for receive
    for (int i = 0; i < nb_dom; i++)
      if (DofNumber(i).GetM() > 0)
        MpiCompleteIrecv(xdom(i), xdom_tmp(i), DofNumber(i).GetM()*nb_u);
    
    // final reduction step
    for (int i = 0; i < nb_dom; i++)
      {
        if (oper == MPI_SUM)
          {
            for (int m = 0; m < nb_u; m++)
              for (int k = 0; k < DofNumber(i).GetM(); k++)
                X(DofNumber(i)(k) + m*Nvol) += xdom(i)(k*nb_u+m);
          }
        else if (oper == MPI_MIN)
          {
            for (int m = 0; m < nb_u; m++)
              for (int k = 0; k < DofNumber(i).GetM(); k++)
                X(DofNumber(i)(k) + m*Nvol) 
                  = minComplex(X(DofNumber(i)(k) + m*Nvol),
			       xdom(i)(k*nb_u+m));
          }
        else if (oper == MPI_MAX)
          {
            for (int m = 0; m < nb_u; m++)
              for (int k = 0; k < DofNumber(i).GetM(); k++)
                X(DofNumber(i)(k) + m*Nvol) 
                  = maxComplex(X(DofNumber(i)(k) + m*Nvol),
			       xdom(i)(k*nb_u+m));
          }
        else
          {
            cout << "Operation not implemented" << endl;
            abort();
          }
      }
  }


  //! assembles a distributed vector (some rows are shared between processors)
  /*!
    \param[inout] X vector to assemble
    \param[in] oper reduction operator (MPI_SUM, MPI_MIN or MPI_MAX)
    \param[in] ProcNumber processors interacting with the current processor
    \param[in] DofNumber for each processor, rows to send/receive
    \param[in] comm MPI communicator
    \param[in] Nvol DofNumber will be included between <0, Nvol-1>
    \param[in] nb_u number of unknowns
    (X is supposed to be of size >= Nvol*nb_u)
    \param[in] tag MPI tag
    For each unknown, it is assumed that the same row numbers are affected,
    that's why you can provide DofNumber for only the first unknown
    and consider that rows of other unknowns are obtained with
    the operation m*Nvol + i, m in <0, nb_u-1> and
    m in <0, Nvol-1>
   */
  template<class T>
  void ExchangeVector(Vector<T>& X,
                      const IVect& ProcNumber, const Vector<IVect>& DofNumber,
                      const MPI_Comm& comm, int Nvol, int nb_u, int tag)
  {
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    // number of procs interacting with the current processor
    int nb_dom = DofNumber.GetM();
    Vector<MPI_Request> request_send(nb_dom), request_recv(nb_dom);
    MPI_Status status;
    Vector<Vector<T> > xsend(nb_dom);
    Vector<Vector<int64_t> > xsend_tmp(nb_dom);
    
    // sending informations to other processors
    for (int i = 0; i < nb_dom; i++)
      {
	int j = ProcNumber(i);
	int nb = DofNumber(i).GetM();
        if (nb > 0)
          {
            xsend(i).Reallocate(nb_u*nb);
            for (int m = 0; m < nb_u; m++)
              for (int k = 0; k < nb; k++)
                xsend(i)(nb_u*k+m) = X(DofNumber(i)(k) + m*Nvol);
            
            // sending the value to the corresponding processor
            request_send(i) = MpiIsend(comm, xsend(i), xsend_tmp(i),
				       nb*nb_u, j, tag);
          }
      }
    
    // receiving the informations
    T zero; SetComplexZero(zero);
    Vector<Vector<T> > xdom(nb_dom);
    Vector<Vector<int64_t> > xdom_tmp(nb_dom);
    for (int i = 0; i < nb_dom; i++)
      {
	int j = ProcNumber(i);
	int nb = DofNumber(i).GetM();
	if (nb > 0)
          {
            xdom(i).Reallocate(nb_u*nb);
            xdom(i).Fill(zero);
            
            // receiving the values of domain j
            request_recv(i) = MpiIrecv(comm, xdom(i), xdom_tmp(i),
				       nb*nb_u, j, tag);
          }
      }
    
    // now waiting all communications are effective
    for (int i = 0; i < nb_dom; i++)
      if (DofNumber(i).GetM() > 0)
        MPI_Wait(&request_send(i), &status);
    
    for (int i = 0; i < nb_dom; i++)
      if (DofNumber(i).GetM() > 0)
        MPI_Wait(&request_recv(i), &status);
    
    // completing writing operations for receive
    for (int i = 0; i < nb_dom; i++)
      if (DofNumber(i).GetM() > 0)
        MpiCompleteIrecv(xdom(i), xdom_tmp(i), DofNumber(i).GetM()*nb_u);
    
    // final reduction step
    for (int i = 0; i < nb_dom; i++)
      {
	for (int m = 0; m < nb_u; m++)
	  for (int k = 0; k < DofNumber(i).GetM(); k++)
	    X(DofNumber(i)(k) + m*Nvol) = xdom(i)(k*nb_u+m);
      }
  }

  
  //! exchanges a distributed vector (some rows are shared between processors)
  //! only values on shared rows are exchanged
  /*!
    \param[inout] X vector to assemble
    \param[in] ProcNumber processors interacting with the current processor
    \param[in] DofNumber for each processor, rows to send/receive
    \param[in] comm MPI communicator
    \param[in] Nvol DofNumber will be included between <0, Nvol-1>
    \param[in] nb_u number of unknowns
    (X is supposed to be of size >= Nvol*nb_u)
    \param[in] tag MPI tag
    For each unknown, it is assumed that the same row numbers are affected,
    that's why you can provide DofNumber for only the first unknown
    and consider that rows of other unknowns are obtained with 
    the operation m*Nvol + i, m in <0, nb_u-1> and
    m in <0, Nvol-1>
   */
  template<class T, class Treal>
  void ExchangeRelaxVector(Vector<T>& X, const Treal& omega, int proc,
			   const IVect& ProcNumber,
			   const Vector<IVect>& DofNumber,
			   const MPI_Comm& comm, int Nvol, int nb_u, int tag)
  {
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    if (nb_proc == 1)
      return;
    
    // number of procs interacting with the current processor
    int nb_dom = DofNumber.GetM();
    MPI_Status status;
    Vector<T> xsend(nb_dom);
    Vector<int64_t> xsend_tmp(nb_dom);
        
    if (rank_proc == proc)
      {
	// processor proc is sending informations to other processors
	for (int i = 0; i < nb_dom; i++)
	  {
	    int j = ProcNumber(i);
	    int nb = DofNumber(i).GetM();
	    if (nb > 0)
	      {
		xsend.Reallocate(nb_u*nb);
		for (int m = 0; m < nb_u; m++)
		  for (int k = 0; k < nb; k++)
		    xsend(nb_u*k+m) = X(DofNumber(i)(k) + m*Nvol);
            
		// sending the value to the corresponding processor
		MpiSsend(comm, xsend, xsend_tmp, nb*nb_u, j, tag);
	      }
	  }
      }
    else
      {
	// others processors are receiving the informations
	Vector<T> xdom(nb_dom);
	Vector<int64_t> xdom_tmp(nb_dom);
	T zero; SetComplexZero(zero);
	for (int i = 0; i < nb_dom; i++)
	  {
	    int j = ProcNumber(i);
	    if (j == proc)
	      {
		int nb = DofNumber(i).GetM();
		if (nb > 0)
		  {
		    xdom.Reallocate(nb_u*nb);
		    xdom.Fill(zero);
		    
		    // receiving the values of domain j
		    MpiRecv(comm, xdom, xdom_tmp, nb*nb_u, j, tag, status);

		    for (int m = 0; m < nb_u; m++)
		      for (int k = 0; k < DofNumber(i).GetM(); k++)
			X(DofNumber(i)(k) + m*Nvol) = 
			  (1.0-omega)* X(DofNumber(i)(k) + m*Nvol) +
			  omega*xdom(k*nb_u+m);
		  }
	      }
	  }
      }
  }
  
}

#define SELDON_FILE_DISTRIBUTED_VECTOR_CXX
#endif

