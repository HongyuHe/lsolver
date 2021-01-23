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

#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_CXX

#include "DistributedMatrix.hxx"

namespace Seldon
{
  
  ///////////////////////////////////
  // DistributedMatrixIntegerArray //


  void DistributedMatrixIntegerArray::Clear()
  {
    OverlapRowNumbers.Clear();
    OverlapProcNumbers.Clear();
    GlobalRowNumbers.Clear();
    ProcSharingRows.Clear();
    SharingRowNumbers.Clear();
  }

  
  void DistributedMatrixIntegerArray
  ::SetData(int nl, int ng, int nodl, int nb_u, const MPI_Comm& comm_,
	    const IVect& overlap_row, const IVect& overlap_proc,
	    const IVect& global_rows, const IVect& proc_rows, const Vector<IVect>& sharing_rows)
  {
    nloc = nl;
    nglob = ng;
    nodl_scalar = nodl;
    nb_unknowns_scal = nb_u;
    comm = comm_;

    OverlapRowNumbers.SetData(overlap_row.GetM(), overlap_row.GetData());
    OverlapProcNumbers.SetData(overlap_proc.GetM(), overlap_proc.GetData());
    GlobalRowNumbers.SetData(global_rows.GetM(), global_rows.GetData());
    ProcSharingRows.SetData(proc_rows.GetM(), proc_rows.GetData());
    SharingRowNumbers.SetData(sharing_rows.GetM(), sharing_rows.GetData());
  }
  

  void DistributedMatrixIntegerArray::Nullify()
  {
    OverlapRowNumbers.Nullify();
    OverlapProcNumbers.Nullify();
    GlobalRowNumbers.Nullify();
    ProcSharingRows.Nullify();
    SharingRowNumbers.Nullify();    
  }

  
  // DistributedMatrixIntegerArray //
  ///////////////////////////////////

  
  ////////////////////////////
  // DistributedMatrix_Base //
  

  /********************
   * Internal methods *
   ********************/


  //! erases any array used for MltAdd
  template<class T>
  void DistributedMatrix_Base<T>::EraseArrayForMltAdd()
  {
    global_row_to_recv.Clear(); global_col_to_recv.Clear();
    ptr_global_row_to_recv.Clear(); ptr_global_col_to_recv.Clear();
    local_row_to_send.Clear(); local_col_to_send.Clear();
    
    proc_col_to_recv.Clear(); proc_col_to_send.Clear();
    proc_row_to_recv.Clear(); proc_row_to_send.Clear();
    
    local_number_distant_values = false;
    size_max_distant_row = 0;
    size_max_distant_col = 0;
  }
  

  //! erases informations for matrix-vector product
  //! and reverts dist_row/dist_col to global numbers
  template<class T>
  void DistributedMatrix_Base<T>
  ::SwitchToGlobalNumbers()
  {
    if (!local_number_distant_values)
      return;
    
    // changing row numbers
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        dist_row(i).Index(j) = global_row_to_recv(dist_row(i).Index(j));
    
    // then col numbers
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        dist_col(i).Index(j) = global_col_to_recv(dist_col(i).Index(j));
    
    // erasing datas needed for matrix-vector product
    EraseArrayForMltAdd();
  }
  

  //! changes global numbers in proc_row to local numbers
  /*!
    \param[in] dist_val list of distant non-zero entries. the distant row
    numbers are contained in dist_val(i).Index(j) for each column i
    \param[in] dist_proc processor rank for each distant non-zero entry
    \param[out] glob_num list of distant row numbers, sorted by processor rank
    it will contain row numbers for processor proc_glob(0),
    then processor proc_glob(1), etc
    \param[out] ptr_glob_num beginning of index for each processor
    in array glob_num the row numbers of processor proc_glob(i)
    in array glob_num are stored 
    in glob_num(ptr_glob_num(i):ptr_glob_num(i+1))
    \param[out] proc_glob list of distant processors involved
    in array dist_num
    \param[out] local_num for each processor proc_local(i)
    local_num(i) will store the local row numbers that should be sent
    to processor proc_local(i) (since these local
    rows will be distant in distant processor proc_local(i)),
    somehow glob_num stores the distants rows that will be received,
    whereas local_num stores the local rows that will be sent
    \param[out] proc_local list of distant processors that will receive
    local row numbers
   */
  template<class T> template<class TypeDist>
  void DistributedMatrix_Base<T>
  ::SortAndAssembleDistantInteractions(TypeDist& dist_val,
				       Vector<IVect>& dist_proc,
                                       IVect& glob_num,
				       IVect& ptr_glob_num, IVect& proc_glob,
                                       Vector<IVect>& local_num,
				       IVect& proc_local)
  {
    MPI_Comm& comm = comm_;
    // counting the number of distant interactions
    long N = 0;
    for (int i = 0; i < dist_proc.GetM(); i++)
      N += dist_proc(i).GetM();
    
    // sorting distant interactions by processor number
    IVect all_proc(N), all_num(N); Vector<long> permut(N);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    Vector<long> nb_inter_per_proc(nb_proc);
    permut.Fill(); long nb = 0;
    nb_inter_per_proc.Fill(0);
    for (int i = 0; i < dist_proc.GetM(); i++)
      for (int j = 0; j < dist_proc(i).GetM(); j++)
        {
          all_proc(nb) = dist_proc(i)(j);
          all_num(nb) = dist_val(i).Index(j);
          nb_inter_per_proc(dist_proc(i)(j))++;
          nb++;
        }
    
    Sort(all_proc, all_num, permut);
    
    // number of processors involved ?
    int nb_global_proc = 0;
    for (int i = 0; i < nb_inter_per_proc.GetM(); i++)
      if (nb_inter_per_proc(i) > 0)
        nb_global_proc++;
    
    proc_glob.Reallocate(nb_global_proc);
    
    // for each processor, sorting numbers,
    // and counting how many different numbers there are
    long offset = 0;
    int nb_glob = 0;
    Vector<int> nb_num_per_proc(nb_proc);
    nb_num_per_proc.Fill(0);
    nb_global_proc = 0;
    for (int p = 0; p < nb_inter_per_proc.GetM(); p++)
      if (nb_inter_per_proc(p) > 0)
        {
          long size = nb_inter_per_proc(p);
          Sort(offset, offset+size-1, all_num, permut);
          int prec = all_num(offset);
          int nb_glob_p = 1;
          for (long j = 1; j < size; j++)
            if (all_num(offset+j) != prec)
              {
                prec = all_num(offset+j);
                nb_glob_p++;
              }
          
          nb_num_per_proc(p) = nb_glob_p;
          nb_glob += nb_glob_p;
          proc_glob(nb_global_proc) = p;
          nb_global_proc++;
          offset += size;
        }
    
    // grouping global numbers    
    all_proc.Clear(); IVect local(N);
    offset = 0; glob_num.Reallocate(nb_glob);
    ptr_glob_num.Reallocate(nb_global_proc+1);
    nb_glob = 0; nb_global_proc = 0; ptr_glob_num(0) = 0;
    for (int p = 0; p < nb_inter_per_proc.GetM(); p++)
      if (nb_inter_per_proc(p) > 0)
        {
          long size = nb_inter_per_proc(p);
          int prec = all_num(offset);
          ptr_glob_num(nb_global_proc+1) = ptr_glob_num(nb_global_proc)
	    + nb_num_per_proc(p);
          
	  glob_num(nb_glob) = prec; nb_glob++;
          int nb_glob_p = 1;
          for (long j = 0; j < size; j++)
            {
              if (all_num(offset+j) != prec)
                {
                  prec = all_num(offset+j);
                  glob_num(nb_glob) = prec;
                  nb_glob_p++; nb_glob++;
                }
              
              local(offset+j) = nb_glob-1;
            }
          
          nb_global_proc++;
          offset += size;
        }
    
    // changing numbers in dist_val
    Vector<long> inv_permut(N);
    for (int i = 0; i < N; i++)
      inv_permut(permut(i)) = i;
    
    long nb_glob_d = 0;
    for (int i = 0; i < dist_val.GetM(); i++)
      for (int j = 0; j < dist_val(i).GetM(); j++)
        {
          int n = inv_permut(nb_glob_d);
          dist_val(i).Index(j) = local(n);
          nb_glob_d++;
        }
    
    // exchanging nb_num_per_proc
    IVect nb_num_send(nb_proc);
    MPI_Alltoall(nb_num_per_proc.GetData(), 1, MPI_INTEGER,
                 nb_num_send.GetData(), 1, MPI_INTEGER, comm);
    
    // sending numbers
    Vector<MPI_Request> request_send(nb_proc), request_recv(nb_proc);
    
    nb_global_proc = 0;
    for (int p = 0; p < nb_proc; p++)
      if (nb_num_per_proc(p) > 0)
        {
          int size = nb_num_per_proc(p);
          MPI_Isend(&glob_num(ptr_glob_num(nb_global_proc)), size,
                    MPI_INTEGER, p, 17, comm, &request_send(p));
          
          nb_global_proc++;
        }
        
    int nb_local_proc = 0;
    for (int p = 0; p < nb_proc; p++)
      if (nb_num_send(p) > 0)
        nb_local_proc++;
    
    local_num.Reallocate(nb_local_proc);
    proc_local.Reallocate(nb_local_proc);
    
    // receiving numbers
    MPI_Status status; nb_local_proc = 0;
    for (int p = 0; p < nb_proc; p++)
      if (nb_num_send(p) > 0)
        {
          proc_local(nb_local_proc) = p;
          local_num(nb_local_proc).Reallocate(nb_num_send(p));
          MPI_Recv(local_num(nb_local_proc).GetData(), nb_num_send(p),
                   MPI_INTEGER, p, 17, comm, &status);
          
          nb_local_proc++;
        }

    for (int i = 0; i < request_send.GetM(); i++)
      if (nb_num_per_proc(i) > 0)
        MPI_Wait(&request_send(i), &status);
    
    // global to local conversion
    IVect Glob_to_local(this->GetGlobalM());
    const IVect& RowNumber = this->GetGlobalRowNumber();
    Glob_to_local.Fill(-1);
    for (int i = 0; i < RowNumber.GetM(); i++)
      Glob_to_local(RowNumber(i)) = i;
    
    // replacing global numbers with local numbers
    for (int i = 0; i < local_num.GetM(); i++)
      for (int j = 0; j < local_num(i).GetM(); j++)
        local_num(i)(j) = Glob_to_local(local_num(i)(j));    
  }
  

  //! internal function
  template<class T> template<class T2>
  void DistributedMatrix_Base<T>
  ::ScatterValues(const Vector<T2>& X, const IVect& num_recv,
                  const IVect& ptr_num_recv, const IVect& proc_recv,
                  const Vector<IVect>& num_send, const IVect& proc_send,
                  Vector<T2>& Xcol) const
  {
    // sending datas
    const MPI_Comm& comm = comm_;
    Vector<Vector<T2> > xsend(proc_send.GetM()), xrecv(proc_recv.GetM());
    Vector<Vector<int64_t> > xsend_tmp(proc_send.GetM()),
      xrecv_tmp(proc_recv.GetM());
    
    int tag = 30;
    Vector<MPI_Request> request_send(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xsend(i).Reallocate(nb);
        for (int j = 0; j < nb; j++)
          xsend(i)(j) = X(num_send(i)(j));
        
        request_send(i) =
          MpiIsend(comm, xsend(i), xsend_tmp(i), nb, proc_send(i), tag);
      }
    
    // receiving datas
    Vector<MPI_Request> request_recv(proc_recv.GetM());
    int N = 0;
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i); N += nb;
        xrecv(i).Reallocate(nb);
        request_recv(i) =
          MpiIrecv(comm, xrecv(i), xrecv_tmp(i), nb, proc_recv(i), tag);
      }
    
    // waiting that transfers are effective
    MPI_Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      MPI_Wait(&request_send(i), &status);

    for (int i = 0; i < request_recv.GetM(); i++)
      MPI_Wait(&request_recv(i), &status);
    
    xsend.Clear();
    // completing receives
    for (int i = 0; i < request_recv.GetM(); i++)
      MpiCompleteIrecv(xrecv(i), xrecv_tmp(i), xrecv(i).GetM());
    
    // values are stored in Xcol
    Xcol.Reallocate(N); N = 0;
    for (int i = 0; i < proc_recv.GetM(); i++)
      for (int j = 0; j < ptr_num_recv(i+1) - ptr_num_recv(i); j++)
        Xcol(N++) = xrecv(i)(j);
  }


  //! internal function
  template<class T> template<class T2>
  void DistributedMatrix_Base<T>::
  AssembleValues(const Vector<T2>& Xcol, const IVect& num_recv,
                 const IVect& ptr_num_recv, const IVect& proc_recv,
                 const Vector<IVect>& num_send, const IVect& proc_send,
		 Vector<T2>& X) const
  {
    // sending datas
    const MPI_Comm& comm = comm_;
    Vector<Vector<T2> > xsend(proc_recv.GetM()), xrecv(proc_send.GetM());    
    Vector<Vector<int64_t> > xsend_tmp(proc_recv.GetM()),
      xrecv_tmp(proc_send.GetM());
    
    int tag = 32, N = 0;
    Vector<MPI_Request> request_send(proc_recv.GetM());
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i);
        xsend(i).Reallocate(nb);
        for (int j = 0; j < nb; j++)
          xsend(i)(j) = Xcol(N++);
        
        request_send(i) =
          MpiIsend(comm, xsend(i), xsend_tmp(i), nb, proc_recv(i), tag);
      }
    
    // receiving datas
    Vector<MPI_Request> request_recv(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xrecv(i).Reallocate(nb);
        request_recv(i) =
          MpiIrecv(comm, xrecv(i), xrecv_tmp(i), nb, proc_send(i), tag);
      }
    
    // waiting that transfers are effective
    MPI_Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      MPI_Wait(&request_send(i), &status);

    for (int i = 0; i < request_recv.GetM(); i++)
      MPI_Wait(&request_recv(i), &status);
    
    xsend.Clear();
    // completing receives
    for (int i = 0; i < request_recv.GetM(); i++)
      MpiCompleteIrecv(xrecv(i), xrecv_tmp(i), xrecv(i).GetM());
    
    // values are added to X
    for (int i = 0; i < num_send.GetM(); i++)
      for (int j = 0; j < num_send(i).GetM(); j++)
        X(num_send(i)(j)) += xrecv(i)(j);
  }

  
  //! assembles the results for each row, by taking the minimum of Yproc
  //! then the minimum of Y
  /*!
    Instead of summing values as in function AssembleRowValues,
    the minimum is taken for Yproc, and if there is equality in Yproc,
    the minimum is taken for Y
   */
  template<class T>
  void DistributedMatrix_Base<T>::
  AssembleValuesMin(const IVect& Xcol, const IVect& Xcol_proc,
                    const IVect& num_recv, const IVect& ptr_num_recv,
		    const IVect& proc_recv,
                    const Vector<IVect>& num_send, const IVect& proc_send,
                    IVect& Y, IVect& Yproc) const
  {
    // sending datas
    const MPI_Comm& comm = comm_;
    Vector<Vector<int> > xsend(proc_recv.GetM()), xrecv(proc_send.GetM());    
    int tag = 35, N = 0;
    Vector<MPI_Request> request_send(proc_recv.GetM());
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i);
        xsend(i).Reallocate(2*nb);
        for (int j = 0; j < nb; j++)
          {
            xsend(i)(j) = Xcol(N);
            xsend(i)(nb+j) = Xcol_proc(N);
            N++;
          }
        
        MPI_Isend(xsend(i).GetDataVoid(), 2*nb,
                  GetMpiDataType(Xcol), proc_recv(i), tag, comm, &request_send(i));
      }
    
    // receiving datas
    Vector<MPI_Request> request_recv(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xrecv(i).Reallocate(2*nb);
        MPI_Irecv(xrecv(i).GetDataVoid(), 2*nb,
                  GetMpiDataType(Xcol), proc_send(i), tag, comm, &request_recv(i));
      }
    
    // waiting that transfers are effective
    MPI_Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      MPI_Wait(&request_send(i), &status);

    for (int i = 0; i < request_recv.GetM(); i++)
      MPI_Wait(&request_recv(i), &status);
    
    xsend.Clear();
    // values are assembled in X
    for (int i = 0; i < num_send.GetM(); i++)
      for (int j = 0; j < num_send(i).GetM(); j++)
        {
          int nb = num_send(i).GetM();
          int proc = xrecv(i)(nb+j);
          int col = xrecv(i)(j);
          if (proc < Yproc(num_send(i)(j)))
            {
              Yproc(num_send(i)(j)) = proc;
              Y(num_send(i)(j)) = col;
            }
          else if (proc == Yproc(num_send(i)(j)))
            {
              if (col < Y(num_send(i)(j)))
                Y(num_send(i)(j)) = col;
            }
        }
  }


  //! assembles the vector (by taking the minimum instead 
  //! of summing for AssembleVec
  /*!
    The minimal Xproc is search, if equality the minimal X is searched
   */
  template<class T>
  void DistributedMatrix_Base<T>
  ::AssembleVecMin(Vector<int>& X, Vector<int>& Xproc) const
  {
    AssembleVectorMin(X, Xproc, *ProcSharingRows, *SharingRowNumbers,
                      comm_, nodl_scalar_, nb_unknowns_scal_, 13);
  }
  
    
  //! removes non-zero entries contained in dist_vec below epsilon
  //! dist_proc is modified in the same way
  template<class T> template<class T0, class TypeDist>
  void DistributedMatrix_Base<T>
  ::RemoveSmallEntryDistant(const T0& epsilon,
                            TypeDist& dist_vec, Vector<IVect>& dist_proc)
  {
    for (int i = 0; i < dist_vec.GetM(); i++)
      {
        int nb = 0, size = dist_vec(i).GetM();
        for (int j = 0; j < size; j++)
          if (abs(dist_vec(i).Value(j)) > epsilon)
            nb++;
        
        if (nb < size)
          {
            IVect num(size), proc(size); Vector<T> val(size);
            for (int j = 0; j < size; j++)
              {
                num(j) = dist_vec(i).Index(j);
                val(j) = dist_vec(i).Value(j);
                proc(j) = dist_proc(i)(j);
              }
            
            dist_vec(i).Reallocate(nb);
            dist_proc(i).Reallocate(nb);
            nb = 0;
            for (int j = 0; j < size; j++)
              if (abs(val(j)) > epsilon)
                {
                  dist_vec(i).Index(nb) = num(j);
                  dist_vec(i).Value(nb) = val(j);
                  dist_proc(i)(nb) = proc(j);
                  nb++;
                }
          }
      }
  }
  

  //! adds contribution of dist_col for GetRowSum
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>
  ::GetRowSumDistantCol(Vector<T0>& vec_sum) const
  {
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        vec_sum(i) += abs(dist_col(i).Value(j));
  }
  
  
  //! adds contribution of dist_row for GetRowSum
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>
  ::GetRowSumDistantRow(Vector<T0>& vec_sum) const
  {
    T0 zero; SetComplexZero(zero);
    Vector<T0> Y(global_row_to_recv.GetM());
    Y.Fill(zero);
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          int jrow = dist_row(i).Index(j);
          Y(jrow) += abs(dist_row(i).Value(j));
        }
    
    AssembleRowValues(Y, vec_sum);
  }
  
  
  //! adds contribution of dist_col for GetColSum
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>
  ::GetColSumDistantCol(Vector<T0>& vec_sum) const
  {
    T0 zero; SetComplexZero(zero);
    Vector<T0> Y(global_col_to_recv.GetM());
    Y.Fill(zero);
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          int jrow = dist_col(i).Index(j);
          Y(jrow) += abs(dist_col(i).Value(j));
        }
    
    AssembleColValues(Y, vec_sum);
  }
  
  
  //! adds contribution of dist_row for GetColSum
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>
  ::GetColSumDistantRow(Vector<T0>& vec_sum) const
  {
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        vec_sum(i) += abs(dist_row(i).Value(j));
  }


  /**********************************************
   * Internal Methods for matrix-vector product *
   **********************************************/
  

  //! prepares a matrix vector product with the distributed matrix
  /*!
    This method prepares the matrix-vector, such that a call to MltAdd
    with a distributed will perform exchanges of values between values
    in order to produce the correct result.
   */
  template<class T>
  void DistributedMatrix_Base<T>::PrepareMltAdd()
  {
    if (local_number_distant_values)
      return;
    
    local_number_distant_values = true;
    const MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;

    long N = 0;
    for (int i = 0; i < dist_col.GetM(); i++)
      N += dist_col(i).GetM();
    
    MPI_Allreduce(&N, &size_max_distant_col, 1, MPI_LONG, MPI_MAX, comm);
    N = 0;
    for (int i = 0; i < dist_row.GetM(); i++)
      N += dist_row(i).GetM();
    
    MPI_Allreduce(&N, &size_max_distant_row, 1, MPI_INTEGER, MPI_MAX, comm);
    
    if (size_max_distant_col > 0)
      SortAndAssembleDistantInteractions(dist_col, proc_col,
					 global_col_to_recv,
                                         ptr_global_col_to_recv,
					 proc_col_to_recv,
                                         local_col_to_send, proc_col_to_send);
    
    if (size_max_distant_row > 0)
      SortAndAssembleDistantInteractions(dist_row, proc_row,
					 global_row_to_recv,
                                         ptr_global_row_to_recv,
					 proc_row_to_recv,
                                         local_row_to_send, proc_row_to_send);
  }
  
  
  //! Sends/receives values of X on distant rows (similar to ScatterColValues)
  template<class T> template<class T2>
  void DistributedMatrix_Base<T>
  ::ScatterRowValues(const Vector<T2>& X, Vector<T2>& Xcol) const
  {
    ScatterValues(X, global_row_to_recv, ptr_global_row_to_recv,
		  proc_row_to_recv,
                  local_row_to_send, proc_row_to_send, Xcol);
  }

  
  //! Sends/receives values of X on distant columns
  /*!
    this method exchanges values of X between processors so that Xcol will
    contain the values of X on distant columns dist_col
   */
  template<class T> template<class T2>
  void DistributedMatrix_Base<T>
  ::ScatterColValues(const Vector<T2>& X, Vector<T2>& Xcol) const
  {
    ScatterValues(X, global_col_to_recv, ptr_global_col_to_recv,
		  proc_col_to_recv,
                  local_col_to_send, proc_col_to_send, Xcol);
  }


  //! Sends/receives values of Xrow on distant rows and adds them to X
  /*!
    this method exchanges values of Xrow between processors,
    and adds them to X
    so that X will sum the results of distant rows dist_row
   */
  template<class T> template<class T2>
  void DistributedMatrix_Base<T>
  ::AssembleRowValues(const Vector<T2>& Xrow, Vector<T2>& X) const
  {
    AssembleValues(Xrow, global_row_to_recv, ptr_global_row_to_recv,
		   proc_row_to_recv,
                   local_row_to_send, proc_row_to_send, X);
  }


  //! Sends/receives values of Xcol on distant columns and adds them to X
  template<class T> template<class T2>
  void DistributedMatrix_Base<T>
  ::AssembleColValues(const Vector<T2>& Xcol, Vector<T2>& X) const
  {
    AssembleValues(Xcol, global_col_to_recv, ptr_global_col_to_recv,
		   proc_col_to_recv,
                   local_col_to_send, proc_col_to_send, X);
  }
  
  
  //! assembles the vector (adds values of shared rows)
  template<class T> template<class T2>
  void DistributedMatrix_Base<T>
  ::AssembleVec(Vector<T2>& X) const
  {
    AssembleVector(X, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
                   comm_, nodl_scalar_, nb_unknowns_scal_, 14);
  }


  //! Y = Y + alpha A X with only distant columns of A
  template<class T>
  template<class T2, class Storage2, class Allocator2,
           class T4, class Storage4, class Allocator4>
  void DistributedMatrix_Base<T>
  ::MltAddCol(const SeldonTranspose& Trans,
              const Vector<T2, Storage2, Allocator2>& X,
              Vector<T4, Storage4, Allocator4>& Y) const
  {
    if (Trans.NoTrans())
      {
	for (int i = 0; i < dist_col.GetM(); i++)
	  for (int j = 0; j < dist_col(i).GetM(); j++)
	    {
	      int jloc = dist_col(i).Index(j);
	      Y(i) += dist_col(i).Value(j)*X(jloc);
	    }
      }
    else
      {	
	T4 zero; SetComplexZero(zero);
	Y.Reallocate(global_col_to_recv.GetM());
	Y.Fill(zero);
	if (Trans.Trans())
	  {
	    for (int i = 0; i < dist_col.GetM(); i++)
	      for (int j = 0; j < dist_col(i).GetM(); j++)
		{
		  int jrow = dist_col(i).Index(j);
		  Y(jrow) += dist_col(i).Value(j)*X(i);
		}
	  }
	else
	  {	
	    for (int i = 0; i < dist_col.GetM(); i++)
	      for (int j = 0; j < dist_col(i).GetM(); j++)
		{
		  int jrow = dist_col(i).Index(j);
		  Y(jrow) += conjugate(dist_col(i).Value(j))*X(i);
		}
	  }
      }
  }
  

  //! Y = Y + alpha A^T X with only distant rows of A
  template<class T>
  template<class T2, class Storage2, class Allocator2,
           class T4, class Storage4, class Allocator4>
  void DistributedMatrix_Base<T>
  ::MltAddRow(const SeldonTranspose& Trans,
              const Vector<T2, Storage2, Allocator2>& X,
              Vector<T4, Storage4, Allocator4>& Y) const
  {
    if (Trans.NoTrans())
      {
	T4 zero; SetComplexZero(zero);
	Y.Reallocate(global_row_to_recv.GetM());
	Y.Fill(zero);
	for (int i = 0; i < dist_row.GetM(); i++)
	  for (int j = 0; j < dist_row(i).GetM(); j++)
	    {
	      int jrow = dist_row(i).Index(j);
	      Y(jrow) += dist_row(i).Value(j)*X(i);
	    }
      }
    else
      {
	if (Trans.Trans())
	  {
	    for (int i = 0; i < dist_row.GetM(); i++)
	      for (int j = 0; j < dist_row(i).GetM(); j++)
		{
		  int jloc = dist_row(i).Index(j);
		  Y(i) += dist_row(i).Value(j)*X(jloc);
		}
	  }
	else
	  {
	    for (int i = 0; i < dist_row.GetM(); i++)
	      for (int j = 0; j < dist_row(i).GetM(); j++)
		{
		  int jloc = dist_row(i).Index(j);
		  Y(i) += conjugate(dist_row(i).Value(j))*X(jloc);
		}
	  }
      }
  }
  
  
  /******************
   * Static methods *
   ******************/
  

  //! adding a non-zero entry, between a local row/column
  //! and a non-local row/column
  /*!
    This method is used only internally
    \param[in] dist_col array to which entry is added
    \param[in] proc_col processors associated with dist_col
    \param[in] jglob global row/column number
    \param[in] proc2 distant processor containing the global row/column
    \param[in] val value of the non-zero entry
  */
  template<class T>
  void DistributedMatrix_Base<T>
  ::AddDistantValue(Vector<T, VectSparse>& dist_col_,
                    IVect& proc_col_, int jglob, int proc2, const T& val)
  {
    int pos = 0;
    int size_row = dist_col_.GetM();
    while ((pos < size_row) && (dist_col_.Index(pos) < jglob))
      pos++;
    
    if ((pos < size_row) && (dist_col_.Index(pos) == jglob))
      {
        // already existing entry
        dist_col_.Value(pos) += val;
      }
    else
      {
        // new entry
        Vector<T> value(size_row);
        IVect index(size_row), proc(size_row);
        for (int k = 0; k < size_row; k++)
          {
            index(k) = dist_col_.Index(k);
            value(k) = dist_col_.Value(k);
            proc(k) = proc_col_(k);
          }
        
        dist_col_.Reallocate(size_row+1);
        proc_col_.Reallocate(size_row+1);
        for (int k = 0; k < pos; k++)
          {
            dist_col_.Index(k) = index(k);
            dist_col_.Value(k) = value(k);
            proc_col_(k) = proc(k);
          }
        
        dist_col_.Index(pos) = jglob;
        dist_col_.Value(pos) = val;
        proc_col_(pos) = proc2;
        for (int k = pos+1; k <= size_row; k++)
          {
            dist_col_.Index(k) = index(k-1);
            dist_col_.Value(k) = value(k-1);
            proc_col_(k) = proc(k-1);
          }
        
      }
  }
  

  //! Sends EntierToSend and FloatToSend to the required processors
  /*!
    nsend_int(i) contains the number of integers to send to processor i
    The targets retrieves the results in EntierToRecv and FloatToRecv
   */
  template<class T>
  void DistributedMatrix_Base<T>
  ::SendAndReceiveDistributed(const MPI_Comm& comm, IVect& nsend_int, Vector<IVect>& EntierToSend, 
                              Vector<Vector<T> >& FloatToSend, IVect& nrecv_int,
                              Vector<IVect>& EntierToRecv, Vector<Vector<T> >& FloatToRecv)
  {
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank; MPI_Comm_rank(comm, &rank);

    Vector<MPI_Request> request(3*nb_proc);
    Vector<Vector<int64_t> > FloatToSend_tmp(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        {
          MPI_Isend(&nsend_int(i), 1, MPI_INTEGER, i, 4, comm, &request(i));
          
          // sending all the values and indices stored to the processor i
          if (nsend_int(i) > 0)
            {               
              MPI_Isend(EntierToSend(i).GetData(), nsend_int(i),
                        MPI_INTEGER, i, 5, comm, &request(i+nb_proc));
              
              if (EntierToSend(i)(0) > 0)
                request(i+2*nb_proc) = 
                  MpiIsend(comm, FloatToSend(i), FloatToSend_tmp(i),
                           FloatToSend(i).GetM(), i, 6);
            }
        }
    
    // receiving the number of entries
    MPI_Status status; 
    nrecv_int.Zero();
    Vector<int64_t> FloatToRecv_tmp;
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        MPI_Recv(&nrecv_int(i), 1, MPI_INTEGER, i, 4, comm, &status);
    
    // waiting for sending of nsend_int effective
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        MPI_Wait(&request(i), &status);
    
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            EntierToRecv(i).Reallocate(nrecv_int(i));
            MPI_Recv(EntierToRecv(i).GetData(), nrecv_int(i),
                     MPI_INTEGER, i, 5, comm, &status);
          }
        else
          EntierToRecv(i).Clear();
      }
    
    // waiting for sending of EntierToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (nsend_int(i) > 0)
	MPI_Wait(&request(i+nb_proc), &status);
    
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            int nb_float = EntierToRecv(i)(0);
            if (nb_float > 0)
              {
                FloatToRecv(i).Reallocate(nb_float);
                MpiRecv(comm, FloatToRecv(i), FloatToRecv_tmp,
                        nb_float, i, 6, status);
              }
            else
              FloatToRecv(i).Clear();
          }
        else
          FloatToRecv(i).Clear();
      }
    
    // waiting for sending of FloatToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (nsend_int(i) > 0)
        if (EntierToSend(i)(0) > 0)
          MPI_Wait(&request(i + 2*nb_proc), &status);
    
    // deleting sending arrays
    for (int i = 0; i < nb_proc; i++)
      {
        nsend_int(i) = 0;
        EntierToSend(i).Clear();
        FloatToSend(i).Clear();
      }
  }
  

  //! Adds received non-zero entries to the matrix B
  template<class T>
  void DistributedMatrix_Base<T>
  ::AddReceivedInteractions(const MPI_Comm& comm, Matrix<T, General, ArrayRowSparse>& B,
                            Vector<IVect>& EntierToRecv, Vector<Vector<T> >& FloatToRecv,
                            IVect& nrecv_int, Vector<IVect>& EntierToSend,
                            Vector<Vector<T> >& FloatToSend, IVect& nsend_int,
                            IVect& Glob_to_local, const IVect& OverlappedCol,
                            const IVect& OverlapProcNumber,
                            Vector<IVect>& procB, bool reorder)
  {
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int nfac = 1;
    if (reorder)
      nfac = 2;
    
    for (int i = 0; i < nb_proc; i++)
      {       
        if (FloatToRecv(i).GetM() > 0)
          {	
            int nrow = EntierToRecv(i)(1);    
            // loop over rows
            nrecv_int(i) = 2; int nrecv_float = 0;
            IVect proc_num;
            for (int j = 0; j < nrow; j++)
              {
                int iglob = EntierToRecv(i)(nrecv_int(i)++);
                int irow = Glob_to_local(iglob);
                int size_row = EntierToRecv(i)(nrecv_int(i)++);
                IVect index(size_row); Vector<T> values(size_row);
                if (reorder)
                  {
                    proc_num.Reallocate(size_row);
                    for (int k = 0; k < size_row; k++)
                      {
                        index(k) = EntierToRecv(i)(nrecv_int(i)++);
                        proc_num(k) = EntierToRecv(i)(nrecv_int(i)++);
                        values(k) = FloatToRecv(i)(nrecv_float++);
                      }
                  }
                else
                  {
                    for (int k = 0; k < size_row; k++)
                      {
                        index(k) = EntierToRecv(i)(nrecv_int(i)++);
                        values(k) = FloatToRecv(i)(nrecv_float++);
                      }
                  }
                
                // adding to matrix B if the row is not shared
                // otherwise we send the row to the original processor
                if (OverlappedCol(irow) == -1)
                  {
                    if (reorder)
                      {
                        // checking if index is sorted
                        bool index_sorted = true;
                        for (int k = 1; k < size_row; k++)
                          if (index(k) < index(k-1))
                            index_sorted = false;

                        if (!index_sorted)
                          Sort(size_row, index, values, proc_num);
                        
                        int size_rowB = B.GetRowSize(irow);
                        IVect new_col(size_rowB + size_row), new_proc(size_rowB + size_row);
                        Vector<T> new_val(size_rowB + size_row);
                        int p = 0, nb = 0;
                        for (int k = 0; k < size_row; k++)
                          {
                            while ((p < size_rowB) && (B.Index(irow, p) < index(k)))
                              {
                                new_col(nb) = B.Index(irow, p);
                                new_val(nb) = B.Value(irow, p);
                                new_proc(nb) = procB(irow)(p);
                                p++; nb++;
                              }
                            
                            if ((p < size_rowB) && (B.Index(irow, p) == index(k)))
                              {
                                // the value is added
                                new_col(nb) = index(k);
                                new_val(nb) = values(k) + B.Value(irow, p);
                                new_proc(nb) = procB(irow)(p);
                                nb++; p++;
                              }
                            else
                              {
                                // the value is created
                                new_col(nb) = index(k);
                                new_val(nb) = values(k);
                                new_proc(nb) = proc_num(k);
                                nb++;                                
                              }
                          }

                        while (p < size_rowB)
                          {
                            new_col(nb) = B.Index(irow, p);
                            new_val(nb) = B.Value(irow, p);
                            new_proc(nb) = procB(irow)(p);
                            p++; nb++;
                          }
                        
                        B.ReallocateRow(irow, nb);
                        procB(irow).Reallocate(nb);
                        for (int k = 0; k < nb; k++)
                          {
                            B.Index(irow, k) = new_col(k);
                            B.Value(irow, k) = new_val(k);
                            procB(irow)(k) = new_proc(k);
                          }
                      }
                    else
                      {
                        if (size_row == 1)
                          B.AddInteraction(irow, index(0), values(0));
                        else
                          B.AddInteractionRow(irow, size_row, index, values);
                      }
                  }
                else
                  {
                    int irow_ = OverlappedCol(irow);
                    int proc = OverlapProcNumber(irow_);
                    
                    int offset_int(2), offset_float(0);
                    if (nsend_int(proc) == 0)
                      {
                        nsend_int(proc) = 2;
                        EntierToSend(proc).Reallocate(nfac*size_row+4);
                        FloatToSend(proc).Reallocate(size_row);
                        EntierToSend(proc)(0) = 0;
                        EntierToSend(proc)(1) = 0;
                      }
                    else
                      {
                        offset_int = EntierToSend(proc).GetM();
                        offset_float = FloatToSend(proc).GetM();
                        EntierToSend(proc).Resize(nfac*size_row+2+offset_int);
                        FloatToSend(proc).Resize(size_row+offset_float);
                      }
                    
                    nsend_int(proc) += nfac*size_row+2;
                    EntierToSend(proc)(0) += size_row;
                    EntierToSend(proc)(1)++;
                    EntierToSend(proc)(offset_int++) = iglob;
                    EntierToSend(proc)(offset_int++) = size_row;
                    for (int k = 0; k < size_row; k++)
                      {
                        EntierToSend(proc)(offset_int++) = index(k);
                        if (reorder)
                          EntierToSend(proc)(offset_int++) = proc_num(k);
                        
                        FloatToSend(proc)(offset_float++) = values(k);
                      }                    
                  }
              }
          }
      }
  }


  //! Adds received non-zero entries to the matrix B
  template<class T>
  void DistributedMatrix_Base<T>
  ::AddReceivedInteractions(const MPI_Comm& comm, Matrix<T, General, ArrayColSparse>& B,
                            Vector<IVect>& EntierToRecv, Vector<Vector<T> >& FloatToRecv,
                            IVect& nrecv_int, Vector<IVect>& EntierToSend,
                            Vector<Vector<T> >& FloatToSend, IVect& nsend_int,
                            IVect& Glob_to_local, const IVect& OverlappedCol,
                            const IVect& OverlapProcNumber,
                            Vector<IVect>& procB, bool reorder)
  {
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int nfac = 1;
    if (reorder)
      nfac = 2;
    
    for (int i = 0; i < nb_proc; i++)
      {       
        if (FloatToRecv(i).GetM() > 0)
          {	
            int ncol = EntierToRecv(i)(1);    
            // loop over column
            nrecv_int(i) = 2; int nrecv_float = 0;
            IVect proc_num;
            for (int j = 0; j < ncol; j++)
              {
                int iglob = EntierToRecv(i)(nrecv_int(i)++);
                int irow = Glob_to_local(iglob);
                int size_col = EntierToRecv(i)(nrecv_int(i)++);
                IVect index(size_col); Vector<T> values(size_col);
                if (reorder)
                  {
                    proc_num.Reallocate(size_col);
                    for (int k = 0; k < size_col; k++)
                      {
                        index(k) = EntierToRecv(i)(nrecv_int(i)++);
                        proc_num(k) = EntierToRecv(i)(nrecv_int(i)++);
                        values(k) = FloatToRecv(i)(nrecv_float++);
                      }
                  }
                else
                  {
                    for (int k = 0; k < size_col; k++)
                      {
                        index(k) = EntierToRecv(i)(nrecv_int(i)++);
                        values(k) = FloatToRecv(i)(nrecv_float++);
                      }
                  }
                
                // adding to matrix B if the row is not shared
                // otherwise we send the row to the original processor
                if (OverlappedCol(irow) == -1)
                  {
                    if (reorder)
                      {
                        // checking if index is sorted
                        bool index_sorted = true;
                        for (int k = 1; k < size_col; k++)
                          if (index(k) < index(k-1))
                            index_sorted = false;

                        if (!index_sorted)
                          Sort(size_col, index, values, proc_num);
                        
                        int size_colB = B.GetColumnSize(irow);
                        IVect new_col(size_colB + size_col), new_proc(size_colB + size_col);
                        Vector<T> new_val(size_colB + size_col);
                        int p = 0, nb = 0;
                        for (int k = 0; k < size_col; k++)
                          {
                            while ((p < size_colB) && (B.Index(irow, p) < index(k)))
                              {
                                new_col(nb) = B.Index(irow, p);
                                new_val(nb) = B.Value(irow, p);
                                new_proc(nb) = procB(irow)(p);
                                p++; nb++;
                              }
                            
                            if ((p < size_colB) && (B.Index(irow, p) == index(k)))
                              {
                                // the value is added
                                new_col(nb) = index(k);
                                new_val(nb) = values(k) + B.Value(irow, p);
                                new_proc(nb) = procB(irow)(p);
                                nb++; p++;
                              }
                            else
                              {
                                // the value is created
                                new_col(nb) = index(k);
                                new_val(nb) = values(k);
                                new_proc(nb) = proc_num(k);
                                nb++;                                
                              }
                          }

                        while (p < size_colB)
                          {
                            new_col(nb) = B.Index(irow, p);
                            new_val(nb) = B.Value(irow, p);
                            new_proc(nb) = procB(irow)(p);
                            p++; nb++;
                          }
                        
                        B.ReallocateColumn(irow, nb);
                        procB(irow).Reallocate(nb);
                        for (int k = 0; k < nb; k++)
                          {
                            B.Index(irow, k) = new_col(k);
                            B.Value(irow, k) = new_val(k);
                            procB(irow)(k) = new_proc(k);
                          }
                      }
                    else
                      {
                        if (size_col == 1)
                          B.AddInteraction(index(0), irow, values(0));
                        else
                          B.AddInteractionColumn(irow, size_col, index, values);
                      }
                  }
                else
                  {
                    int irow_ = OverlappedCol(irow);
                    int proc = OverlapProcNumber(irow_);
                    
                    int offset_int(2), offset_float(0);
                    if (nsend_int(proc) == 0)
                      {
                        nsend_int(proc) = 2;
                        EntierToSend(proc).Reallocate(nfac*size_col+4);
                        FloatToSend(proc).Reallocate(size_col);
                        EntierToSend(proc)(0) = 0;
                        EntierToSend(proc)(1) = 0;
                      }
                    else
                      {
                        offset_int = EntierToSend(proc).GetM();
                        offset_float = FloatToSend(proc).GetM();
                        EntierToSend(proc).Resize(nfac*size_col+2+offset_int);
                        FloatToSend(proc).Resize(size_col+offset_float);
                      }
                    
                    nsend_int(proc) += nfac*size_col+2;
                    EntierToSend(proc)(0) += size_col;
                    EntierToSend(proc)(1)++;
                    EntierToSend(proc)(offset_int++) = iglob;
                    EntierToSend(proc)(offset_int++) = size_col;
                    for (int k = 0; k < size_col; k++)
                      {
                        EntierToSend(proc)(offset_int++) = index(k);
                        if (reorder)
                          EntierToSend(proc)(offset_int++) = proc_num(k);
                        
                        FloatToSend(proc)(offset_float++) = values(k);
                      }                    
                  }
              }
          }
      }
  }


  //! clears distant numbers in a sparse matrix
  template<class T> template<class TypeDist>
  void DistributedMatrix_Base<T>
  ::EraseDistantEntries(MPI_Comm& comm, const Vector<bool>& IsRowDropped,
                        const Vector<bool>& IsRowDroppedDistant,
                        TypeDist& dist_row_, Vector<IVect>& proc_row_,
                        TypeDist& dist_col_, Vector<IVect>& proc_col_)
  {
    typedef typename TypeDist::value_type Vect1;
    typedef typename Vect1::value_type T1;
    for (int i = 0; i < dist_col_.GetM(); i++)
      if (IsRowDropped(i))
        {
          dist_col_(i).Clear();
          proc_col_(i).Clear();
        }

    for (int j = 0; j < dist_row_.GetM(); j++)
      {
        int nb = 0;
        for (int iloc = 0; iloc < dist_row_(j).GetM(); iloc++)
          if (IsRowDroppedDistant(dist_row_(j).Index(iloc)))
            nb++;
        
        if (nb > 0)
          {
            int size = dist_row_(j).GetM();
            IVect row_num(size), proc(size); Vector<T1> val(size);
            for (int iloc = 0; iloc < dist_row_(j).GetM(); iloc++)
              {
                row_num(iloc) = dist_row_(j).Index(iloc);
                val(iloc) = dist_row_(j).Value(iloc);
                proc(iloc) = proc_row_(j)(iloc);
              }
            
            dist_row_(j).Reallocate(size-nb);
            proc_row_(j).Reallocate(size-nb);
            nb = 0;
            for (int iloc = 0; iloc < size; iloc++)
              if (!IsRowDroppedDistant(row_num(iloc)))
                {
                  dist_row_(j).Index(nb) = row_num(iloc);
                  dist_row_(j).Value(nb) = val(iloc);
                  proc_row_(j)(nb) = proc(iloc);
                  nb++;
                }
          }
      }   
  }


  /****************
   * Constructors *
   ****************/

  
  //! default constructor
  template<class T>
  DistributedMatrix_Base<T>::DistributedMatrix_Base()
  {
    GlobalRowNumbers = NULL;
    OverlapProcNumbers = NULL;
    OverlapRowNumbers = NULL;
    ProcSharingRows = NULL;
    SharingRowNumbers = NULL;
    nodl_scalar_ = 0;
    nb_unknowns_scal_ = 1;
    nglob_ = 0;
    comm_ = MPI_COMM_SELF;
    
    local_number_distant_values = false;
    size_max_distant_row = 0;
    size_max_distant_col = 0;
  }
  

  //! construction of an m by n matrix
  /*!
    Here m and n are the number of rows and columns of the local matrix
  */
  template<class T>
  DistributedMatrix_Base<T>::DistributedMatrix_Base(int m, int n)
  {
    GlobalRowNumbers = NULL;
    OverlapProcNumbers = NULL;
    OverlapRowNumbers = NULL;
    ProcSharingRows = NULL;
    SharingRowNumbers = NULL;
    nglob_ = m;
    nodl_scalar_ = 0;
    nb_unknowns_scal_ = 1;
    comm_ = MPI_COMM_SELF;
    
    dist_col.Reallocate(m);
    dist_row.Reallocate(n);
    proc_col.Reallocate(m);
    proc_row.Reallocate(n);

    local_number_distant_values = false;
    size_max_distant_row = 0;
    size_max_distant_col = 0;
  }


  //! Initialisation of pointers
  /*!
    This method is mandatory, otherwise pointers are set to NULL
    \param[in] n global number of rows (as if it was on a single processor)
    \param[in] row_num local to global numbering
    \param[in] overlap_num rows already counted in another processor
    \param[in] proc_num original processor for overlapped rows
  */
  template<class T>
  void DistributedMatrix_Base<T>::
  Init(int n, IVect* row_num, IVect* overlap_num, IVect* proc_num,
       int Nvol, int nb_u, IVect* MatchingProc,
       Vector<IVect>* MatchingDofNumber, const MPI_Comm& comm)
  {
    nglob_ = n;
    GlobalRowNumbers = row_num;
    OverlapRowNumbers = overlap_num;
    OverlapProcNumbers = proc_num;
    
    nodl_scalar_ = Nvol;
    nb_unknowns_scal_ = nb_u;
    ProcSharingRows = MatchingProc;
    SharingRowNumbers = MatchingDofNumber;
    
    comm_ = comm;    
  }
  

  //! Initialisation of pointers
  template<class T>
  void DistributedMatrix_Base<T>::
  Init(DistributedMatrixIntegerArray& Ainfo)
  {
    nglob_ = Ainfo.nglob;
    GlobalRowNumbers = &Ainfo.GlobalRowNumbers;
    OverlapRowNumbers = &Ainfo.OverlapRowNumbers;
    OverlapProcNumbers = &Ainfo.OverlapProcNumbers;
    
    nodl_scalar_ = Ainfo.nodl_scalar;
    nb_unknowns_scal_ = Ainfo.nb_unknowns_scal;
    ProcSharingRows = &Ainfo.ProcSharingRows;
    SharingRowNumbers = &Ainfo.SharingRowNumbers;
    
    comm_ = Ainfo.comm;    
  }

  
  //! inits pointers with those of A
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>
  ::Init(const DistributedMatrix_Base<T0>& A)
  {
    OverlapRowNumbers = A.OverlapRowNumbers;
    OverlapProcNumbers = A.OverlapProcNumbers;
    GlobalRowNumbers = A.GlobalRowNumbers;
    ProcSharingRows = A.ProcSharingRows;
    SharingRowNumbers = A.SharingRowNumbers;
    nodl_scalar_ = A.nodl_scalar_;
    nb_unknowns_scal_ = A.nb_unknowns_scal_;
    nglob_ = A.nglob_;
    comm_ = A.comm_;
  }
  

  //! initialisation of a distributed matrix with global row numbers
  template<class T>
  void DistributedMatrix_Base<T>
  ::Init(Vector<IVect>& all_rows,
         IVect& row_num, IVect& overlap_num, IVect& proc_num,
         IVect& MatchingProc, Vector<IVect>& MatchingDofNumber,
	 const MPI_Comm& comm, bool distribute_row)
  {
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    // throwing an error if the rows are not sorted
    if ((rank_proc == 0) && (distribute_row))
      for (int i = 0; i < nb_proc; i++)
        {
          bool sorted = true;
          for (int j = 1; j < all_rows(i).GetM(); j++)
            if (all_rows(i)(j) < all_rows(i)(j-1))
              sorted = false;
          
          if (!sorted)
            {
              cout << "Row numbers must be sorted in this function" << endl;
              abort();
              //Sort(all_rows(i));
            }
        }
    
    // finding the global number of rows
    int nglob = 0;
    if (rank_proc == 0)
      for (int i = 0; i < all_rows.GetM(); i++)
        for (int j = 0; j < all_rows(i).GetM(); j++)
          if (all_rows(i)(j)+1 > nglob)
            nglob = all_rows(i)(j)+1;
    
    if (rank_proc == 0)
      {
        // sending row numbers to all processors
        for (int i = 1; i < nb_proc; i++)
          {
            int n = all_rows(i).GetM();
            MPI_Send(&nglob, 1, MPI_INTEGER, i, 101, comm);
            if (distribute_row)
              {
                MPI_Send(&n, 1, MPI_INTEGER, i, 102, comm);
                MPI_Send(all_rows(i).GetData(), all_rows(i).GetM(),
                         MPI_INTEGER, i, 103, comm);
              }
          }
        
        if (distribute_row)
          row_num = all_rows(0);
        
        // minimal processor for each row
        Vector<int> MinProc(nglob);
        MinProc.Fill(nb_proc);
        for (int i = 0; i < all_rows.GetM(); i++)
          for (int j = 0; j < all_rows(i).GetM(); j++)
            if (MinProc(all_rows(i)(j)) > i)
              MinProc(all_rows(i)(j)) = i;
        
        // counting all the rows that are overlapped 
	// (i.e. shared with a processor of lower rank)
        for (int i = 0; i < all_rows.GetM(); i++)
          {
            int nb_overlap = 0;
            for (int j = 0; j < all_rows(i).GetM(); j++)
              if (MinProc(all_rows(i)(j)) < i)
                nb_overlap++;
            
            IVect num(nb_overlap), proc(nb_overlap);
            nb_overlap = 0;
            for (int j = 0; j < all_rows(i).GetM(); j++)
              if (MinProc(all_rows(i)(j)) < i)
                {
                  num(nb_overlap) = j;            
                  proc(nb_overlap) = MinProc(all_rows(i)(j));
                  nb_overlap++;
                }
            
            // sending overlapped rows
            if (i == 0)
              {
                overlap_num = num;
                proc_num = proc;
              }
            else
              {
                MPI_Send(&nb_overlap, 1, MPI_INTEGER, i, 104, comm);
                if (nb_overlap > 0)
                  {
                    MPI_Send(num.GetData(), nb_overlap,
                             MPI_INTEGER, i, 105, comm);
                    MPI_Send(proc.GetData(), nb_overlap,
                             MPI_INTEGER, i, 106, comm);
                  }
              }
          }
        
        // counting the number of processors sharing each global row
        MinProc.Fill(0);
        for (int i = 0; i < all_rows.GetM(); i++)
          for (int j = 0; j < all_rows(i).GetM(); j++)
            MinProc(all_rows(i)(j))++;
        
        int nb_shared_row = 0;
        for (int i = 0; i < nglob; i++)
          {
            if (MinProc(i) > 1)
              {
                MinProc(i) = nb_shared_row;
                nb_shared_row++;
              }
            else
              MinProc(i) = -1;
          }
        
        // for each shared row, storing the list of processors
        Vector<IVect> ListProcRow(nb_shared_row);
        for (int i = 0; i < all_rows.GetM(); i++)
          for (int j = 0; j < all_rows(i).GetM(); j++)
            {
              int n = all_rows(i)(j);
              if (MinProc(n) >= 0)
                {
                  int nloc = MinProc(n);
                  ListProcRow(nloc).PushBack(i);
                }
            }
        
        // then arrays MatchingDofNumber are constructed
        for (int i = 0; i < all_rows.GetM(); i++)
          {
            // searching all the processors interacting
            Vector<bool> ProcUsed(nb_proc);
            ProcUsed.Fill(false);
            for (int j = 0; j < all_rows(i).GetM(); j++)
              {
                int n = all_rows(i)(j);
                if (MinProc(n) >= 0)
                  {
                    int nloc = MinProc(n);
                    for (int k = 0; k < ListProcRow(nloc).GetM(); k++)
                      ProcUsed(ListProcRow(nloc)(k)) = true;
                  }
              }
            
            int nb_proc_interac = 0;
            for (int k = 0; k < ProcUsed.GetM(); k++)
              if ((k != i) && (ProcUsed(k)))
                nb_proc_interac++;
            
            Vector<int> matching_proc(nb_proc_interac);
            Vector<Vector<int> > matching_row(nb_proc_interac);
            nb_proc_interac = 0;
            for (int k = 0; k < ProcUsed.GetM(); k++)
              if ((k != i) && (ProcUsed(k)))
                {
                  // counting rows shared with processor k
                  int nb_row_interac = 0;
                  for (int j = 0; j < all_rows(i).GetM(); j++)
                    {
                      int n = all_rows(i)(j);
                      if (MinProc(n) >= 0)
                        {
                          int nloc = MinProc(n);
                          for (int k2 = 0;
			       k2 < ListProcRow(nloc).GetM(); k2++)
                            if (ListProcRow(nloc)(k2) == k)
                              nb_row_interac++;
                        }
                    }
                  
                  // filling arrays
                  matching_proc(nb_proc_interac) = k;
                  matching_row(nb_proc_interac).Reallocate(nb_row_interac);
                  nb_row_interac = 0;
                  for (int j = 0; j < all_rows(i).GetM(); j++)
                    {
                      int n = all_rows(i)(j);
                      if (MinProc(n) >= 0)
                        {
                          int nloc = MinProc(n);
                          for (int k2 = 0;
			       k2 < ListProcRow(nloc).GetM(); k2++)
                            if (ListProcRow(nloc)(k2) == k)
                              {
                                matching_row(nb_proc_interac)
				  (nb_row_interac) = j;
                                
				nb_row_interac++;
                              }
                        }
                    }
                                   
                  nb_proc_interac++;
                }
            
            // sending arrays to the processor i
            if (i == 0)
              {
                MatchingProc = matching_proc;
                MatchingDofNumber = matching_row;
              }
            else
              {
                MPI_Send(&nb_proc_interac, 1, MPI_INTEGER, i, 107, comm);
                if (nb_proc_interac > 0)
                  {
                    MPI_Send(matching_proc.GetData(), nb_proc_interac,
                             MPI_INTEGER, i, 108, comm);
                    
                    for (int k = 0; k < nb_proc_interac; k++)
                      {
                        int nb_row = matching_row(k).GetM();
                        MPI_Send(&nb_row, 1, MPI_INTEGER, i, 109, comm);
                        MPI_Send(matching_row(k).GetData(), nb_row,
                                 MPI_INTEGER, i, 110, comm);
                      }
                  }
              }
          }
      }
    else
      {
        // receiving row numbers
        MPI_Status status; int n;
        MPI_Recv(&nglob, 1, MPI_INTEGER, 0, 101, comm, &status);
        if (distribute_row)
          {
            MPI_Recv(&n, 1, MPI_INTEGER, 0, 102, comm, &status);
            row_num.Reallocate(n);
            MPI_Recv(row_num.GetData(), n, MPI_INTEGER, 0, 103, comm, &status);
          }
        
        // receiving overlapped numbers
        MPI_Recv(&n, 1, MPI_INTEGER, 0, 104, comm, &status);
        if (n > 0)
          {
            overlap_num.Reallocate(n);
            proc_num.Reallocate(n);
            MPI_Recv(overlap_num.GetData(), n, MPI_INTEGER,0, 105, comm, &status);
            MPI_Recv(proc_num.GetData(), n, MPI_INTEGER, 0, 106, comm, &status);
          }
        else
          {
            overlap_num.Clear();
            proc_num.Clear();
          }
        
        // receiving local numbers of rows shared with other processors
        MPI_Recv(&n, 1, MPI_INTEGER, 0, 107, comm, &status);
        if (n > 0)
          {
            MatchingProc.Reallocate(n);
            MatchingDofNumber.Reallocate(n);
            MPI_Recv(MatchingProc.GetData(), n,
                     MPI_INTEGER, 0, 108, comm, &status);
            
	    for (int k = 0; k < MatchingProc.GetM(); k++)
              {
                MPI_Recv(&n, 1, MPI_INTEGER, 0, 109, comm, &status);
                MatchingDofNumber(k).Reallocate(n);
                MPI_Recv(MatchingDofNumber(k).GetData(), n, MPI_INTEGER,
                         0, 110, comm, &status);
              }
          }
        else
          {
            MatchingProc.Clear();
            MatchingDofNumber.Clear();
          }
      }
    
    // once the arrays are constructed, we store the pointers
    nglob_ = nglob;
    GlobalRowNumbers = &row_num;
    OverlapRowNumbers = &overlap_num;
    OverlapProcNumbers = &proc_num;
    
    nodl_scalar_ = row_num.GetM();
    nb_unknowns_scal_ = 1;
    ProcSharingRows = &MatchingProc;
    SharingRowNumbers = &MatchingDofNumber;
    
    comm_ = comm;
  }
  

  //! initialisation of a distributed matrix with global row numbers
  template<class T>
  void DistributedMatrix_Base<T>
  ::Init(IVect& row_num, IVect& overlap_num, IVect& proc_num,
         IVect& MatchingProc,
	 Vector<IVect>& MatchingDofNumber, const MPI_Comm& comm)
  {
    // throwing an error if the rows are not sorted
    bool sorted = true;
    for (int j = 1; j < row_num.GetM(); j++)
      if (row_num(j) < row_num(j-1))
        sorted = false;
    
    if (!sorted)
      {
        cout << "Row numbers must be sorted in this function" << endl;
        abort();        
        //Sort(row_num);
      }
    
    // gathering row numbers to the processor 0
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    Vector<IVect> all_rows(nb_proc);
    if (rank_proc == 0)
      {
        all_rows(0) = row_num;
        int nodl_par; MPI_Status status;
        for (int i = 1; i < nb_proc; i++)
          {
            MPI_Recv(&nodl_par, 1, MPI_INTEGER, i, 13, comm, &status);
            all_rows(i).Reallocate(nodl_par);
            MPI_Recv(all_rows(i).GetData(), nodl_par,
                     MPI_INTEGER, i, 14, comm, &status);
          }
      }
    else
      {
        int nodl = row_num.GetM();
        MPI_Send(&nodl, 1, MPI_INTEGER, 0, 13, comm);
        MPI_Send(row_num.GetData(), nodl, MPI_INTEGER, 0, 14, comm);
      }
    
    // then calling Init with all_rows
    Init(all_rows, row_num, overlap_num, proc_num,
         MatchingProc, MatchingDofNumber, comm, false);
  }



  /*********************
   * Memory Management *
   *********************/
    

  //! changing the size of the local matrix, previous values are lost
  template<class T>
  void DistributedMatrix_Base<T>
  ::ReallocateDist(int m, int n)
  {
    // previous values are erased for simplicity
    Clear();
        
    dist_col.Reallocate(m);
    dist_row.Reallocate(n);
    proc_col.Reallocate(m);
    proc_row.Reallocate(n);
  }
  
  
  //! changing the size of the local matrix, previous values are kept
  template<class T>
  void DistributedMatrix_Base<T>::Resize(int m, int n)
  {
    if (dist_col.GetM() != m)
      {
        dist_col.Resize(m);
        proc_col.Resize(m);
      }

    if (dist_row.GetM() != n)
      {
        dist_row.Resize(n);
        proc_row.Resize(n);
      }
  }
  

  //! matrix is cleared
  template<class T>
  void DistributedMatrix_Base<T>::Clear()
  {
    dist_col.Clear();
    proc_col.Clear();
    dist_row.Clear();
    proc_row.Clear();
    
    global_row_to_recv.Clear(); global_col_to_recv.Clear();
    ptr_global_row_to_recv.Clear(); ptr_global_col_to_recv.Clear();
    local_row_to_send.Clear(); local_col_to_send.Clear();
    proc_col_to_recv.Clear(); proc_col_to_send.Clear();
    proc_row_to_recv.Clear(); proc_row_to_send.Clear();
    size_max_distant_row = 0;
    size_max_distant_col = 0;
    local_number_distant_values = false;
  }
  
    
  /*******************
   * Basic functions *
   *******************/
  

  //! equality *this = X 
  template<class T>
  DistributedMatrix_Base<T>& DistributedMatrix_Base<T>
  ::operator=(const DistributedMatrix_Base<T>& X)
  {
    Copy(X);
    return *this;
  }


  //! Copies content of X in the current object
  template<class T> template<class T2>
  void DistributedMatrix_Base<T>
  ::Copy(const DistributedMatrix_Base<T2>& X)
  {
    dist_col.Reallocate(X.dist_col.GetM());
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i).Copy(X.dist_col(i));
 
    dist_row.Reallocate(X.dist_row.GetM());
    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i).Copy(X.dist_row(i));
 
    proc_col = X.proc_col;
    proc_row = X.proc_row;
    
    GlobalRowNumbers = X.GlobalRowNumbers;
    OverlapProcNumbers = X.OverlapProcNumbers;
    OverlapRowNumbers = X.OverlapRowNumbers;

    ProcSharingRows = X.ProcSharingRows;
    SharingRowNumbers = X.SharingRowNumbers;
    
    nodl_scalar_ = X.nodl_scalar_;
    nb_unknowns_scal_ = X.nb_unknowns_scal_;
    nglob_ = X.nglob_;
    
    comm_ = X.comm_;
    
    global_row_to_recv = X.global_row_to_recv;
    global_col_to_recv = X.global_col_to_recv;
    ptr_global_row_to_recv = X.ptr_global_row_to_recv;
    ptr_global_col_to_recv = X.ptr_global_col_to_recv;
    
    local_row_to_send = X.local_row_to_send;
    local_col_to_send = X.local_col_to_send;
    proc_col_to_recv = X.proc_col_to_recv;
    proc_col_to_send = X.proc_col_to_send;
    proc_row_to_recv = X.proc_row_to_recv;
    proc_row_to_send = X.proc_row_to_send;
    local_number_distant_values = X.local_number_distant_values;
    
    size_max_distant_row = X.size_max_distant_row;
    size_max_distant_col = X.size_max_distant_col;
  }
  
    
  //! multiplication by a scalar
  template<class T> template<class T0>
  DistributedMatrix_Base<T>&
  DistributedMatrix_Base<T>::operator *=(const T0& x)
  {
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i) *= x;

    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i) *= x;
    
    return *this;
  }
  
  
  //! returns local to global numbering
  template<class T>
  const IVect& DistributedMatrix_Base<T>::
  GetGlobalRowNumber() const
  {
    if (this->GlobalRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *GlobalRowNumbers;
  }
  
  
  //! returns local to global numbering
  template<class T>
  IVect& DistributedMatrix_Base<T>
  ::GetGlobalRowNumber()
  {
    if (this->GlobalRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *GlobalRowNumbers;
  }


  //! returns rows already counted in another processor
  template<class T>
  const IVect& DistributedMatrix_Base<T>::
  GetOverlapRowNumber() const
  {
    if (this->OverlapRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapRowNumbers;
  }
  

  //! returns rows already counted in another processor
  template<class T>
  IVect& DistributedMatrix_Base<T>
  ::GetOverlapRowNumber()
  {
    if (this->OverlapRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapRowNumbers;
  }

  
  //! returns processor numbers of the original rows
  template<class T>
  const IVect& DistributedMatrix_Base<T>::GetOverlapProcNumber() const
  {
    if (this->OverlapProcNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapProcNumbers;
  }
  

  //! returns processor numbers of the original rows
  template<class T>
  IVect& DistributedMatrix_Base<T>
  ::GetOverlapProcNumber()
  {
    if (this->OverlapProcNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapProcNumbers;
  }


  //! returns processor numbers for each set of shared rows
  template<class T>
  IVect& DistributedMatrix_Base<T>
  ::GetProcessorSharingRows()
  {
    if (this->ProcSharingRows == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *ProcSharingRows;
  }
  

  //! returns processor numbers for each set of shared rows
  template<class T>
  const IVect& DistributedMatrix_Base<T>
  ::GetProcessorSharingRows() const
  {
    if (this->ProcSharingRows == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *ProcSharingRows;
  }

  
  //! returns row numbers for each set of shared rows
  template<class T>
  Vector<IVect>& DistributedMatrix_Base<T>
  ::GetSharingRowNumbers()
  {
    if (this->SharingRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *SharingRowNumbers;
  }


  //! returns row numbers for each set of shared rows
  template<class T>
  const Vector<IVect>& DistributedMatrix_Base<T>
  ::GetSharingRowNumbers() const
  {
    if (this->SharingRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *SharingRowNumbers;
  }
  
      
  //! returns the size of memory used to store the matrix
  template<class T>
  size_t DistributedMatrix_Base<T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + sizeof(int*)*(proc_row.GetM()+proc_col.GetM());
    taille += sizeof(int*)*(local_row_to_send.GetM() + local_col_to_send.GetM()
                            + dist_row.GetM() + dist_col.GetM());
    
    taille += global_row_to_recv.GetMemorySize() + global_col_to_recv.GetMemorySize() +
      ptr_global_row_to_recv.GetMemorySize() + ptr_global_col_to_recv.GetMemorySize() +
      proc_col_to_recv.GetMemorySize() + proc_col_to_send.GetMemorySize() +
      proc_row_to_recv.GetMemorySize() + proc_row_to_send.GetMemorySize();
    
    for (int i = 0; i < proc_row.GetM(); i++)
      taille += proc_row(i).GetMemorySize();

    for (int i = 0; i < proc_col.GetM(); i++)
      taille += proc_col(i).GetMemorySize();
    
    for (int i = 0; i < local_row_to_send.GetM(); i++)
      taille += local_row_to_send(i).GetMemorySize();

    for (int i = 0; i < local_col_to_send.GetM(); i++)
      taille += local_col_to_send(i).GetMemorySize();
    
    for (int i = 0; i < dist_row.GetM(); i++)
      taille += dist_row(i).GetMemorySize();

    for (int i = 0; i < dist_col.GetM(); i++)
      taille += dist_col(i).GetMemorySize();
    
    return taille;
  }

  
  /**********************
   * Convenient methods *
   **********************/
  
  
  //! returns the number of non-zero entries stored in the matrix
  /*!
    It returns the number of non-zero entries stored on the current processor.
    To obtain an "estimation" of the global number of non-zero entries,
    MPI_Reduce should be called
   */
  template<class T>
  long DistributedMatrix_Base<T>::GetNonZeros() const
  {
    long nnz = 0;
    for (int i = 0; i < dist_col.GetM(); i++)
      nnz += dist_col(i).GetM();

    for (int i = 0; i < dist_row.GetM(); i++)
      nnz += dist_row(i).GetM();
    
    return nnz;
  }
  
  
  //! returns the number of elements stored
  /*!
    It returns the number of elements (of type T)
    stored on the current processor.
    To obtain an "estimation" of the global number of elements
    MPI_Reduce should be called
  */
  template<class T>
  long DistributedMatrix_Base<T>::GetDataSize() const
  {    
    long size_int = 0, size = 0;
    for (int i = 0; i < dist_col.GetM(); i++)
      {
        size_int += 2*dist_col(i).GetM();
        size += dist_col(i).GetM();
      }

    for (int i = 0; i < dist_row.GetM(); i++)
      {
        size_int += 2*dist_row(i).GetM();
        size += dist_row(i).GetM();
      }
    
    size_int += 2*(global_row_to_recv.GetM() + global_col_to_recv.GetM());
    
    for (int i = 0; i < local_row_to_send.GetM(); i++)
      size_int += local_row_to_send(i).GetM();

    for (int i = 0; i < local_col_to_send.GetM(); i++)
      size_int += local_col_to_send(i).GetM();
    
    size_int += 2*(dist_col.GetM() + dist_row.GetM()) +
      2*(local_row_to_send.GetM() + local_col_to_send.GetM()) + 25;
    
    int ratio = sizeof(T)/sizeof(int);

    return size + size_int/ratio;
  }
  

  //! removes small entries of the matrix
  /*!
    The result of this function will be different from a sequential
    execution since values stored on different processors
    may add up and be eliminated or not
   */
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>
  ::RemoveSmallEntry(const T0& epsilon)
  {
    RemoveSmallEntryDistant(epsilon, dist_col, proc_col);
    RemoveSmallEntryDistant(epsilon, dist_row, proc_row);
  }
  
  
  //! sets matrix to identity matrix
  template<class T>
  void DistributedMatrix_Base<T>::SetIdentity()
  {
    for (int i = 0; i < dist_col.GetM(); i++)
      {
        dist_col(i).Clear();
        proc_col(i).Clear();
      }
    
    for (int i = 0; i < dist_row.GetM(); i++)
      {
        dist_row(i).Clear();
        proc_row(i).Clear();
      }      
    
    EraseArrayForMltAdd();
  }
  
    
  //! sets values of non-zero entries to 0
  template<class T>
  void DistributedMatrix_Base<T>::Zero()
  {
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i).Zero();
    
    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i).Zero();    
  }
  
  
  //! sets values of non-zero entries to 0, 1, 2, etc
  /*!
    The result will be different from a sequential execution
    since values can be stored in different processors
   */
  template<class T>
  void DistributedMatrix_Base<T>::Fill()
  {
    long value(0);
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          SetComplexReal(value, dist_col(i).Value(j));
          value++;
        }
    
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          SetComplexReal(value, dist_row(i).Value(j));
          value++;
        }    
  }
  
  
  //! sets values of non-zero entries to x
  /*!
    The result will be different from a sequential execution
    since values can be stored in different processors
   */
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>::Fill(const T0& x)
  {
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i).Fill(x);

    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i).Fill(x);
  }
  
  
  //! sets values of non-zero entries to random values
  template<class T>
  void DistributedMatrix_Base<T>::FillRand()
  {
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i).FillRand();

    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i).FillRand();
  }
  
  
  //! writes the matrix in text format
  template<class T>
  void DistributedMatrix_Base<T>
  ::WriteText(ostream& FileStream, Vector<int>& IndRow, Vector<int>& IndCol,
              Vector<T>& Value, bool cplx) const
  {        
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("DistributedMatrix::WriteText(ofstream& FileStream)",
		    "Stream is not ready.");
#endif
        
    // extending arrays in order to contain non-local part
    long N = 0;
    for (int i = 0; i < dist_col.GetM(); i++)
      N += dist_col(i).GetM();

    for (int i = 0; i < dist_row.GetM(); i++)
      N += dist_row(i).GetM();

    long old_size = IndRow.GetM();
    IndRow.Resize(old_size+N);
    IndCol.Resize(old_size+N);
    Value.Resize(old_size+N);
    
    // filling non-local part
    const IVect& global = this->GetGlobalRowNumber();
    N = old_size;
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          IndRow(N) = global(i) + 1;
          IndCol(N) = dist_col(i).Index(j) + 1;
          if (local_number_distant_values)
            IndCol(N) = global_col_to_recv(dist_col(i).Index(j)) + 1;
          
          Value(N) = dist_col(i).Value(j);
          N++;
        }

    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          IndCol(N) = global(i) + 1;
          IndRow(N) = dist_row(i).Index(j) + 1;
          if (local_number_distant_values)
            IndRow(N) = global_row_to_recv(dist_row(i).Index(j)) + 1;
          
          Value(N) = dist_row(i).Value(j);
          N++;
        }
    
    // changing numbers of local part
    for (int i = 0; i < old_size; i++)
      {
        IndRow(i) = global(IndRow(i)) + 1;
        IndCol(i) = global(IndCol(i)) + 1;
      }

    // writing values on the stream
    WriteCoordinateMatrix(FileStream, IndRow, IndCol, Value, cplx);
  }
  
  
  /********************************************
   * Methods called for matrix-vector product *
   ********************************************/
    

  //! Initializes the matrix-vector product  
  template<class T> 
  template<class T2, class T3, class T4, class Storage4, class Allocator4>
  void DistributedMatrix_Base<T>::
  InitMltAdd(bool& proceed_distant_row, bool& proceed_distant_col,
             const Vector<T2>& X, Vector<T2>& Xcol,
             const T3& beta, Vector<T4, Storage4, Allocator4>& Y,
             Vector<T4, Storage4, Allocator4>& Yres) const
  {
    const MPI_Comm& comm = this->GetCommunicator();
    proceed_distant_row = true;
    proceed_distant_col = true;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        proceed_distant_row = false;
        proceed_distant_col = false;
      }
    else
      {
        if (!this->IsReadyForMltAdd())
          {
            // preparing the matrix vector product
            // this method will be called once for
	    // the first matrix-vector product
            const_cast<DistributedMatrix_Base<T>& >(*this)
              .PrepareMltAdd();
          }
        
        if (this->GetMaxDataSizeDistantCol() == 0)
          proceed_distant_col = false;
        
        if (this->GetMaxDataSizeDistantRow() == 0)
          proceed_distant_row = false;
      }
    
    T3 zero;
    SetComplexZero(zero);
    
    if (beta == zero)
      Y.SetData(Yres.GetM(), Yres.GetData());
    else
      Y.Reallocate(Yres.GetM());
    
    Y.Fill(zero);        

    // scattering column values
    if (proceed_distant_col)
      this->ScatterColValues(X, Xcol);    
  }
  
  
  //! Finalizes the matrix-vector product
  template<class T> 
  template<class T2, class T3, class T4, class Storage4, class Allocator4>
  void DistributedMatrix_Base<T>::
  FinalizeMltAdd(bool proceed_distant_row, bool proceed_distant_col,
                 const Vector<T2>& X, Vector<T2>& Xcol, const T3& alpha,
                 const T3& beta, Vector<T4, Storage4, Allocator4>& Y,
                 Vector<T4, Storage4, Allocator4>& Yres, bool assemble) const
  {
    // adding contributions of distant columns
    if (proceed_distant_col)
      this->MltAddCol(SeldonNoTrans, Xcol, Y);
    
    // contributions of distant rows
    Vector<T4> Yrow;
    if (proceed_distant_row)
      this->MltAddRow(SeldonNoTrans, X, Yrow);

    // assembling row values
    if (proceed_distant_row)
      this->AssembleRowValues(Yrow, Y);
    
    // assembling rows shared between processors
    if (assemble)
      this->AssembleVec(Y);

    T3 zero;
    SetComplexZero(zero);

    if (beta == zero)
      {
        Mlt(alpha, Y);
        Y.Nullify();
      }
    else
      {
        Mlt(beta, Yres);
        Add(alpha, Y, Yres);
      }
  }


  //! Initializes the matrix-vector product  
  template<class T> 
  template<class T2, class T3, class T4, class Storage4, class Allocator4>
  void DistributedMatrix_Base<T>::
  InitMltAdd(bool& proceed_distant_row, bool& proceed_distant_col,
             const SeldonTranspose& trans, const Vector<T2>& X, Vector<T2>& Xrow,
             const T3& beta, Vector<T4, Storage4, Allocator4>& Y,
             Vector<T4, Storage4, Allocator4>& Yres) const
  {
    const MPI_Comm& comm = this->GetCommunicator();
    proceed_distant_row = true;
    proceed_distant_col = true;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        proceed_distant_row = false;
        proceed_distant_col = false;
      }
    else
      {
        if (!this->IsReadyForMltAdd())
          {
            // preparing the matrix vector product
            // this method will be called once for 
	    // the first matrix-vector product
            const_cast<DistributedMatrix_Base<T>& >(*this)
              .PrepareMltAdd();
          }
        
        if (this->GetMaxDataSizeDistantCol() == 0)
          proceed_distant_col = false;
        
        if (this->GetMaxDataSizeDistantRow() == 0)
          proceed_distant_row = false;
      }
    
    T3 zero;
    SetComplexZero(zero);
    
    if (beta == zero)
      Y.SetData(Yres.GetM(), Yres.GetData());
    else
      Y.Reallocate(Yres.GetM());
    
    Y.Fill(zero);
    
    // scattering row values
    if (proceed_distant_row)
      this->ScatterRowValues(X, Xrow);
  }


  //! Finalizes the matrix-vector product
  template<class T> 
  template<class T2, class T3, class T4, class Storage4, class Allocator4>
  void DistributedMatrix_Base<T>::
  FinalizeMltAdd(bool proceed_distant_row, bool proceed_distant_col,
                 const SeldonTranspose& trans, const Vector<T2>& X, Vector<T2>& Xrow,
                 const T3& alpha, const T3& beta, Vector<T4, Storage4, Allocator4>& Y,
                 Vector<T4, Storage4, Allocator4>& Yres, bool assemble) const
  {
    // adding contributions of distant rows
    if (proceed_distant_row)
      this->MltAddRow(trans, Xrow, Y);
    
    // contributions of distant columns
    Vector<T4> Ycol;
    if (proceed_distant_col)
      this->MltAddCol(trans, X, Ycol);

    // assembling row values
    if (proceed_distant_col)
      this->AssembleColValues(Ycol, Y);

    // assembling rows shared between processors
    if (assemble)
      this->AssembleVec(Y);

    T3 zero;
    SetComplexZero(zero);

    if (beta == zero)
      {
        Mlt(alpha, Y);
        Y.Nullify();
      }
    else
      {
        Mlt(beta, Yres);
        Add(alpha, Y, Yres);
      }
  }


  //! Initializes the matrix-vector product MltMin  
  template<class T>
  void DistributedMatrix_Base<T>
  ::InitMltMin(Vector<int>& Y, Vector<int>& Yproc,
               Vector<int>& Xcol, Vector<int>& Xcol_proc) const
  {
    const MPI_Comm& comm = this->GetCommunicator();
    if (!this->IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_Base<T>& >(*this)
          .PrepareMltAdd();
      }
    
    // scattering column values
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc > 1)
      {
        this->ScatterColValues(Y, Xcol);
        this->ScatterColValues(Yproc, Xcol_proc);
      }
  }


  //! Finalizes the matrix-vector product MltMin
  template<class T>
  void DistributedMatrix_Base<T>
  ::FinalizeMltMin(Vector<int>& Y, Vector<int>& Yproc,
                   Vector<int>& Xcol, Vector<int>& Xcol_proc) const
  {
    const MPI_Comm& comm = this->GetCommunicator();

    int N = this->global_col_to_recv.GetM();  
    Vector<int> Yrow(N), Yrow_proc(N);
  
    // contributions of distant columns
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    Yrow.Fill(0); Yrow_proc.Fill(nb_proc);
    for (int i = 0; i < this->dist_col.GetM(); i++)
      for (int j = 0; j < this->dist_col(i).GetM(); j++)
        {
          int k = this->dist_col(i).Index(j);
          if (Xcol_proc(k) < Yproc(i))
            {
              Yproc(i) = Xcol_proc(k);
              Y(i) = Xcol(k);
            }
          else if (Xcol_proc(k) == Yproc(i))
            {
              if (Xcol(k) < Y(i))
                Y(i) = Xcol(k);
            }

          if (Yproc(i) < Yrow_proc(k))
            {
              Yrow_proc(k) = Yproc(i);
              Yrow(k) = Y(i);
            }
          else if (Yproc(i) == Yrow_proc(k))
            {
              if (Y(i) < Yrow(k))
                Yrow(k) = Y(i);
            }
        }
    
    this->AssembleValuesMin(Yrow, Yrow_proc,
                            global_col_to_recv, ptr_global_col_to_recv, 
                            proc_col_to_recv, local_col_to_send,
                            proc_col_to_send, Y, Yproc);
    
    // contributions of distant rows
    if (nb_proc > 1)
      {
        ScatterRowValues(Y, Xcol);
        ScatterRowValues(Yproc, Xcol_proc);
      }

    N = global_row_to_recv.GetM();
    Yrow.Reallocate(N), Yrow_proc.Reallocate(N);
    Yrow.Fill(0); Yrow_proc.Fill(nb_proc);
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          int jrow = dist_row(i).Index(j);
          if (Xcol_proc(jrow) < Yproc(i))
            {
              Yproc(i) = Xcol_proc(jrow);
              Y(i) = Xcol(jrow);
            }
          else if (Xcol_proc(jrow) == Yproc(i))
            {
              if (Xcol(jrow) < Y(i))
                Y(i) = Xcol(jrow);
            }

          if (Yproc(i) < Yrow_proc(jrow))
            {
              Yrow_proc(jrow) = Yproc(i);
              Yrow(jrow) = Y(i);
            }
          else if (Yproc(i) == Yrow_proc(jrow))
            {
              if (Y(i) < Yrow(jrow))
                Yrow(jrow) = Y(i);
            }
        }

    // assembling row values
    if (nb_proc > 1)
      {
        AssembleValuesMin(Yrow, Yrow_proc,
                          global_row_to_recv, ptr_global_row_to_recv,
                          proc_row_to_recv, local_row_to_send,
                          proc_row_to_send, Y, Yproc);
        
        // assembling rows shared between processors
        AssembleVecMin(Y, Yproc);
      }
  }


  /****************************************************
   * Methods called for various functions on matrices *
   ****************************************************/

  
  //! Adds alpha A to the current matrix (only distant values)
  template<class T> template<class T0, class T1>
  void DistributedMatrix_Base<T>
  ::AddDistributedMatrix(const T0& alpha,
                         const DistributedMatrix_Base<T1>& A)
  {
    const_cast<DistributedMatrix_Base<T1>& >(A)
      .SwitchToGlobalNumbers();
    
    this->SwitchToGlobalNumbers();
    
    // adding distant interactions
    for (int i = 0; i < A.dist_row.GetM(); i++)
      for (int j = 0; j < A.dist_row(i).GetM(); j++)
        this->AddRowDistantInteraction(A.dist_row(i).Index(j), i,
                                       A.proc_row(i)(j),
                                       alpha*A.dist_row(i).Value(j));
    
    for (int i = 0; i < A.dist_col.GetM(); i++)
      for (int j = 0; j < A.dist_col(i).GetM(); j++)
        this->AddDistantInteraction(i, A.dist_col(i).Index(j),
                                    A.proc_col(i)(j),
                                    alpha*A.dist_col(i).Value(j));
  }
  

  //! Computes res = max(res, maximum absolute value of non-zero entries)
  template<class T>
  void DistributedMatrix_Base<T>
  ::GetMaxAbsDistant(typename ClassComplexType<T>::Treal& res) const
  {
    const MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;

    for (int i = 0; i < this->dist_row.GetM(); i++)
      for (int j = 0; j < this->dist_row(i).GetM(); j++)
        res = max(res, abs(this->dist_row(i).Value(j)));

    for (int i = 0; i < this->dist_col.GetM(); i++)
      for (int j = 0; j < this->dist_col(i).GetM(); j++)
        res = max(res, abs(this->dist_col(i).Value(j)));
    
    // selecting the maximum between processors    
    Vector<int64_t> xtmp;
    typename ClassComplexType<T>::Treal amax(0);
    MpiAllreduce(comm, &res, xtmp, &amax, 1, MPI_MAX);
    
    res = amax;
  }
  

  //! Adds \sum |a_ij| to vec_sum(i) for distant non-zero entries 
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>::AddRowSumDistant(Vector<T0>& vec_sum) const
  {
    const MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    if (!this->IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_Base<T>& >(*this)
          .PrepareMltAdd();
      }
    
    this->GetRowSumDistantCol(vec_sum);    
    this->GetRowSumDistantRow(vec_sum);
    
    this->AssembleVec(vec_sum);
  }


  //! Adds \sum |a_ij| to vec_sum(j) for distant non-zero entries 
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>::AddColSumDistant(Vector<T0>& vec_sum) const
  {
    const MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    if (!this->IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_Base<T>& >(*this)
          .PrepareMltAdd();
      }
    
    this->GetColSumDistantCol(vec_sum);    
    this->GetColSumDistantRow(vec_sum);
    
    this->AssembleVec(vec_sum);
  }


  //! Adds \sum |a_ij| to sum_row(i) and sum_col(j) for distant non-zero entries 
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>::AddRowColSumDistant(Vector<T0>& sum_row,
                                                      Vector<T0>& sum_col) const
  {
    const MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    if (!this->IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_Base<T>& >(*this)
          .PrepareMltAdd();
      }
    
    this->GetRowSumDistantCol(sum_row);
    this->GetRowSumDistantRow(sum_row);
    
    this->GetColSumDistantCol(sum_col);    
    this->GetColSumDistantRow(sum_col);
    
    this->AssembleVec(sum_row);
    this->AssembleVec(sum_col);
  }


  //! Swaps arrays contained in the current structure with arrays given as parameters
  template<class T>
  void DistributedMatrix_Base<T>
  ::ExchangeParallelData(int& smax_row, int& smax_col, bool& local_number,
                         Vector<Vector<T, VectSparse>, VectFull,
                         NewAlloc<Vector<T, VectSparse> > >& dist_row_,
                         Vector<Vector<T, VectSparse>, VectFull,
                         NewAlloc<Vector<T, VectSparse> > >& dist_col_,
                         Vector<IVect>& proc_row_, Vector<IVect>& proc_col_,
                         IVect& global_row_to_recv_, IVect& global_col_to_recv_,
                         IVect& ptr_global_row_to_recv_, IVect& ptr_global_col_to_recv_,
                         Vector<IVect>& local_row_to_send_, Vector<IVect>& local_col_to_send_,
                         IVect& proc_row_to_recv_, IVect& proc_col_to_recv_,
                         IVect& proc_row_to_send_, IVect& proc_col_to_send_)
  {             
    const MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
           
    long stmp = smax_row;
    smax_row = size_max_distant_row;
    size_max_distant_row = stmp;
    
    stmp = smax_col;
    smax_col = size_max_distant_col;
    size_max_distant_col = stmp;
    
    bool btmp = local_number;
    local_number = local_number_distant_values;
    local_number_distant_values = btmp;

    SwapPointer(dist_row, dist_row_);
    SwapPointer(dist_col, dist_col_);
    SwapPointer(proc_row, proc_row_);
    SwapPointer(proc_col, proc_col_);
    SwapPointer(global_row_to_recv, global_row_to_recv_);
    SwapPointer(global_col_to_recv, global_col_to_recv_);
    SwapPointer(ptr_global_row_to_recv, ptr_global_row_to_recv_);
    SwapPointer(ptr_global_col_to_recv, ptr_global_col_to_recv_);
    SwapPointer(local_row_to_send, local_row_to_send_);
    SwapPointer(local_col_to_send, local_col_to_send_);
    SwapPointer(proc_row_to_recv, proc_row_to_recv_);
    SwapPointer(proc_col_to_recv, proc_col_to_recv_);
    SwapPointer(proc_row_to_send, proc_row_to_send_);
    SwapPointer(proc_col_to_send, proc_col_to_send_);
  }


  //! Conjugates distant non-zero entries
  template<class T>
  void DistributedMatrix_Base<T>::ConjugateDistant()
  {
    const MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    for (int i = 0; i < this->dist_row.GetM(); i++)
      for (int j = 0; j < this->dist_row(i).GetM(); j++)
        this->dist_row(i).Value(j) = conjugate(this->dist_row(i).Value(j));
    
    for (int i = 0; i < this->dist_col.GetM(); i++)
      for (int j = 0; j < this->dist_col(i).GetM(); j++)
        this->dist_col(i).Value(j) = conjugate(this->dist_col(i).Value(j));
  }
  

  //! Computes *this = A^T  for non-zero entries
  template<class T>
  void DistributedMatrix_Base<T>::TransposeDistant(const DistributedMatrix_Base<T>& A)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    this->Init(A);
    this->dist_col = A.dist_row;
    this->proc_col = A.proc_row;
    this->global_col_to_recv = A.global_row_to_recv;
    this->ptr_global_col_to_recv = A.ptr_global_row_to_recv;
    this->local_col_to_send = A.local_row_to_send;
    this->proc_col_to_recv = A.proc_row_to_recv;
    this->proc_col_to_send = A.proc_row_to_send;

    this->dist_row = A.dist_col;
    this->proc_row = A.proc_col;
    this->global_row_to_recv = A.global_col_to_recv;
    this->ptr_global_row_to_recv = A.ptr_global_col_to_recv;
    this->local_row_to_send = A.local_col_to_send;
    this->proc_row_to_recv = A.proc_col_to_recv;
    this->proc_row_to_send = A.proc_col_to_send;
    this->local_number_distant_values = A.local_number_distant_values;
    
    this->size_max_distant_row = A.size_max_distant_col;
    this->size_max_distant_col = A.size_max_distant_row;
  }


  //! Multiplies a_ij by Drow(i) for distant non-zero entries
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>::ScaleLeftDistant(const Vector<T0>& Drow)
  {
    MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    if (!this->IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_Base<T>& >(*this)
          .PrepareMltAdd();
      }
    
    // scaling distant columns
    for (int i = 0; i < this->dist_col.GetM(); i++)
      this->dist_col(i) *= Drow(i);
    
    // scaling distant rows
    Vector<T0> Drow_glob;
    this->ScatterRowValues(Drow, Drow_glob);
    
    for (int i = 0; i < this->dist_row.GetM(); i++)
      for (int j = 0; j < this->dist_row(i).GetM(); j++)
        this->dist_row(i).Value(j) *= Drow_glob(this->dist_row(i).Index(j));
  }
  

  //! Multiplies a_ij by Dcol(j) for distant non-zero entries
  template<class T> template<class T0>
  void DistributedMatrix_Base<T>::ScaleRightDistant(const Vector<T0>& Dcol)
  {
    MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;

    if (!this->IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_Base<T>& >(*this)
          .PrepareMltAdd();
      }
    
    // scaling distant rows
    for (int i = 0; i < this->dist_row.GetM(); i++)
      this->dist_row(i) *= Dcol(i);
    
    // scaling distant columns
    Vector<T0>  Dcol_glob;
    this->ScatterColValues(Dcol, Dcol_glob);
    
    for (int i = 0; i < this->dist_col.GetM(); i++)
      for (int j = 0; j < this->dist_col(i).GetM(); j++)
        this->dist_col(i).Value(j) *= Dcol_glob(this->dist_col(i).Index(j));    
  }


  //! Multiplies a_ij by Drow(i) Dcol(i) for distant non-zero entries
  template<class T> template<class T0, class T1>
  void DistributedMatrix_Base<T>
  ::ScaleDistant(const Vector<T0>& Drow, const Vector<T1>& Dcol)
  {
    MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;

    if (!this->IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_Base<T>& >(*this)
          .PrepareMltAdd();
      }
    
    // scaling distant columns
    for (int i = 0; i < this->dist_col.GetM(); i++)
      this->dist_col(i) *= Drow(i);
    
    // scaling distant rows
    Vector<T0> Drow_glob;
    this->ScatterRowValues(Drow, Drow_glob);
    
    for (int i = 0; i < this->dist_row.GetM(); i++)
      for (int j = 0; j < this->dist_row(i).GetM(); j++)
        this->dist_row(i).Value(j) *= Drow_glob(this->dist_row(i).Index(j));

    // scaling distant rows
    for (int i = 0; i < this->dist_row.GetM(); i++)
      this->dist_row(i) *= Dcol(i);
    
    // scaling distant columns
    Vector<T1> Dcol_glob;
    this->ScatterColValues(Dcol, Dcol_glob);
    
    for (int i = 0; i < this->dist_col.GetM(); i++)
      for (int j = 0; j < this->dist_col(i).GetM(); j++)
        this->dist_col(i).Value(j) *= Dcol_glob(this->dist_col(i).Index(j));    
  }


  //! erases columns num(i) 
  template<class T>
  void DistributedMatrix_Base<T>::EraseColDistant(const IVect& num, bool sym)
  {
    MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    if (!this->IsReadyForMltAdd())
      this->PrepareMltAdd();
    
    Vector<bool> IsColDropped(this->GetLocalM());
    IsColDropped.Fill(false);
    for (int i = 0; i < num.GetM(); i++)
      IsColDropped(num(i)) = true;
    
    Vector<bool> IsColDroppedDistant;
    this->ScatterColValues(IsColDropped, IsColDroppedDistant);
    
    EraseDistantEntries(comm, IsColDropped, IsColDroppedDistant,
                        this->dist_col, this->proc_col, this->dist_row, this->proc_row); 

    if (sym)
      {
        this->ScatterRowValues(IsColDropped, IsColDroppedDistant);
        
        EraseDistantEntries(comm, IsColDropped, IsColDroppedDistant,
                            this->dist_row, this->proc_row, this->dist_col, this->proc_col); 
      }
  }


  //! erases rows num(i) 
  template<class T>
  void DistributedMatrix_Base<T>::EraseRowDistant(const IVect& num, bool sym)
  {
    MPI_Comm& comm = this->GetCommunicator();    
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    if (!this->IsReadyForMltAdd())
      this->PrepareMltAdd();
    
    Vector<bool> IsRowDropped(this->GetLocalM());
    IsRowDropped.Fill(false);
    for (int i = 0; i < num.GetM(); i++)
      IsRowDropped(num(i)) = true;
    
    Vector<bool> IsRowDroppedDistant;
    this->ScatterRowValues(IsRowDropped, IsRowDroppedDistant);
    
    EraseDistantEntries(comm, IsRowDropped, IsRowDroppedDistant,
                        this->dist_row, this->proc_row, this->dist_col, this->proc_col); 

    if (sym)
      {
        this->ScatterColValues(IsRowDropped, IsRowDroppedDistant);
        
        EraseDistantEntries(comm, IsRowDropped, IsRowDroppedDistant,
                            this->dist_col, this->proc_col, this->dist_row, this->proc_row); 
      }
  }


  //! Copies A(row, col) to the current matrix
  template<class T> template<class T1>
  void DistributedMatrix_Base<T>
  ::CopySubDistant(const DistributedMatrix_Base<T1>& A,
                   const IVect& row, const IVect& col, bool sym)
  {
    MPI_Comm& comm = this->GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    if (!A.IsReadyForMltAdd())
      const_cast<DistributedMatrix_Base<T1>& >(A)
        .PrepareMltAdd();
    
    this->Init(A);
    int m = A.GetLocalM(), n = A.GetLocalN();    
    Vector<bool> RowKept(m), ColKept(n);
    RowKept.Fill(false); ColKept.Fill(false);
    for (int i = 0; i < row.GetM(); i++)
      RowKept(row(i)) = true;
    
    for (int i = 0; i < col.GetM(); i++)
      ColKept(col(i)) = true;
    
    // checking consistency of row/col with symmetry of matrix B
    if (sym)
      {
        if (m != n)
          {
            cout << "A is non-symmetric while B is symmetric" << endl;
            abort();
          }
        
        for (int i = 0; i < m; i++)
          if (RowKept(i) ^ ColKept(i))
            {
              cout << "row and col must be identic to obtain "
		   << "a symmetric matrix" << endl;
              abort();
            }
      }

    Vector<bool> RowKeptDistant, ColKeptDistant;
    A.ScatterRowValues(RowKept, RowKeptDistant);
    A.ScatterColValues(ColKept, ColKeptDistant);
    
    // using global numbers for B
    this->SwitchToGlobalNumbers();
    
    // extracting values of dist_col
    this->dist_col.Reallocate(m);
    this->proc_col.Reallocate(m);
    for (int i = 0; i < A.dist_col.GetM(); i++)
      {
        if (RowKept(i))
          {
            int size_row = 0;
            for (int j = 0; j < A.dist_col(i).GetM(); j++)
              if (ColKeptDistant(A.dist_col(i).Index(j)))
                size_row++;
            
            this->dist_col(i).Reallocate(size_row);
            this->proc_col(i).Reallocate(size_row);
            size_row = 0;
            for (int j = 0; j < A.dist_col(i).GetM(); j++)
              if (ColKeptDistant(A.dist_col(i).Index(j)))
                {
                  this->dist_col(i).Value(size_row) = A.dist_col(i).Value(j);
                  this->dist_col(i).Index(size_row)
                    = A.global_col_to_recv(A.dist_col(i).Index(j));
                  
                  this->proc_col(i)(size_row) = A.proc_col(i)(j);
                  size_row++;
                }
          }
        else
          {
            this->dist_col(i).Clear();
            this->proc_col(i).Clear();
          }
      }
    
    // same stuff for distant rows
    this->dist_row.Reallocate(n);
    this->proc_row.Reallocate(n);
    for (int i = 0; i < A.dist_row.GetM(); i++)
      {
        if (ColKept(i))
          {
            int size_col = 0;
            for (int j = 0; j < A.dist_row(i).GetM(); j++)
              if (RowKeptDistant(A.dist_row(i).Index(j)))
                size_col++;
            
            this->dist_row(i).Reallocate(size_col);
            this->proc_row(i).Reallocate(size_col);
            size_col = 0;
            for (int j = 0; j < A.dist_row(i).GetM(); j++)
              if (RowKeptDistant(A.dist_row(i).Index(j)))
                {
                  this->dist_row(i).Value(size_col) = A.dist_row(i).Value(j);
                  this->dist_row(i).Index(size_col)
                    = A.global_row_to_recv(A.dist_row(i).Index(j));
                  
                  this->proc_row(i)(size_col) = A.proc_row(i)(j);
                  size_col++;
                }
          }
        else
          {
            this->dist_row(i).Clear();
            this->proc_row(i).Clear();
          }

      }
  }


  /*************************************************
   * Methods for parallel assembling of the matrix *
   *************************************************/
  

  //! Method called by AssembleDistributed
  template<class T>
  void DistributedMatrix_Base<T>
  ::AssembleParallel(Matrix<T, General, ArrayRowSparse>& B, Vector<IVect>& procB,
                     Symmetric& sym, IVect& row_numbers, IVect& local_row_numbers,
		     IVect& OverlappedCol, bool sym_pattern, bool reorder)
  {
    int m = this->GetLocalM();
    int n = this->GetGlobalM();
    
    MPI_Comm& comm = this->GetCommunicator();    
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank; MPI_Comm_rank(comm, &rank);
    
    IVect RowNumber(this->GetGlobalRowNumber());
    const IVect& OverlapRowNumber = this->GetOverlapRowNumber();
    const IVect& OverlapProcNumber = this->GetOverlapProcNumber();

    // constructing array to know if a column is overlapped 
    // (treated by another processor)
    OverlappedCol.Reallocate(m); OverlappedCol.Fill(-1);
    for (int i = 0; i < OverlapRowNumber.GetM(); i++)
      OverlappedCol(OverlapRowNumber(i)) = i;

    /*********************************
     * Parallel assembling of matrix *
     *********************************/
    
    // in this section, we will send/receive overlapped rows and distant rows
    Vector<int> nb_row_sent(nb_proc);
    Vector<int> nsend_int(nb_proc), nsend_float(nb_proc);
    Vector<IVect> EntierToSend(nb_proc);
    Vector<Vector<T> > FloatToSend(nb_proc);
    
    // counting the size of arrays to send
    // nsend_int : the number of integers
    // nsend_float : the number of floats
    for (int i = 0; i < nb_proc; i++)
      {
        if (i != rank)
          nsend_int(i) = 2;
        else
          nsend_int(i) = 0;
        
        nsend_float(i) = 0;
        nb_row_sent(i) = 0;
      }

    // overlapped rows
    for (int j = 0; j < OverlapRowNumber.GetM(); j++)
      {
        int i = OverlapProcNumber(j);
        if (i != rank)
          {
            nsend_int(i) += 2;
            nb_row_sent(i)++;
            int irow = OverlapRowNumber(j);
            nsend_int(i) += B.GetRowSize(irow);
            nsend_float(i) += B.GetRowSize(irow);
            if (reorder)
              nsend_int(i) += B.GetRowSize(irow);
          }
      }
    
    // distant rows
    for (int j = 0; j < this->dist_row.GetM(); j++)
      for (int k = 0; k < this->dist_row(j).GetM(); k++)
        {
          int i = this->proc_row(j)(k);
          if (i != rank)
            {
              int irow = this->dist_row(j).Index(k);
              if (this->local_number_distant_values)
                irow = this->global_row_to_recv(irow);
              
              if (irow <= RowNumber(j))
                {
                  nb_row_sent(i)++;
                  nsend_int(i) += 3;
                  nsend_float(i)++;
                  if (reorder)
                    nsend_int(i)++;
                }
            }
        }
    
    // allocating arrays EntierToSend and FloatToSend
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        { 
          if (nb_row_sent(i) == 0)
            nsend_int(i) = 0;
          
          if (nb_row_sent(i) > 0)
            {
              EntierToSend(i).Reallocate(nsend_int(i));
              FloatToSend(i).Reallocate(nsend_float(i));
              EntierToSend(i)(0) = nsend_float(i);
              EntierToSend(i)(1) = nb_row_sent(i);
              nsend_int(i) = 2; nsend_float(i) = 0;
            }
        }

    // then arrays EntierToSend and FloatToSend are filled

    // storing values and indices of a row shared with processor i
    // processor i is the owner of this row
    for (int j = 0; j < OverlapRowNumber.GetM(); j++)
      {
        int i = OverlapProcNumber(j);
        if (i != rank)
          {
            int irow = OverlapRowNumber(j);
            EntierToSend(i)(nsend_int(i)++) = RowNumber(irow);
            EntierToSend(i)(nsend_int(i)++) = B.GetRowSize(irow);
            for (int k = 0; k < B.GetRowSize(irow); k++)
              {
                EntierToSend(i)(nsend_int(i)++) = B.Index(irow, k);
                if (reorder)
                  EntierToSend(i)(nsend_int(i)++) = procB(irow)(k);
                
                FloatToSend(i)(nsend_float(i)++) = B.Value(irow, k);
              }

            // the corresponding values of B are cleared
            // they are no longer needed since they will be present in the distant processor
            // after the exchange of datas
            B.ClearRow(irow);
            if (reorder)
              procB(irow).Clear();
          }
      }
    
    // storing values of row associated with processor rank
    for (int j = 0; j < this->dist_row.GetM(); j++)
      for (int k = 0; k < this->dist_row(j).GetM(); k++)
        {
          int i = this->proc_row(j)(k);
          if (i != rank)
            {
              int irow = this->dist_row(j).Index(k);
              if (this->local_number_distant_values)
                irow = this->global_row_to_recv(irow);
              
              if (irow <= RowNumber(j))
                {
                  EntierToSend(i)(nsend_int(i)++) = irow;
                  EntierToSend(i)(nsend_int(i)++) = 1;
                  EntierToSend(i)(nsend_int(i)++) = RowNumber(j);
                  if (reorder)
                    EntierToSend(i)(nsend_int(i)++) = rank;
                  
                  FloatToSend(i)(nsend_float(i)++)
                    = this->dist_row(j).Value(k);
                }
            }
        }
    
    // now the initial matrix can be cleared
    this->Clear();

    // Datas for receiving EntierToSend and FloatToSend
    Vector<IVect> EntierToRecv(nb_proc);
    Vector<Vector<T> > FloatToRecv(nb_proc);
    IVect nrecv_int(nb_proc);
    
    // exchanging datas
    SendAndReceiveDistributed(comm, nsend_int, EntierToSend, FloatToSend,
                              nrecv_int, EntierToRecv, FloatToRecv);
    
    // constructing local row numbers 
    int nloc = m - OverlapRowNumber.GetM();
    local_row_numbers.Reallocate(nloc);
    int nrow = 0;
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
	local_row_numbers(nrow++) = i;
    
    // index array to obtain local numbers in array local_row_numbers
    // from local numbers of the matrix
    IVect inv_local_row_numbers(m);
    inv_local_row_numbers.Fill(-1);
    nrow = 0;
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
        inv_local_row_numbers(i) = nrow++;
    
    // global to local conversion
    IVect Glob_to_local(n);
    Glob_to_local.Fill(-1);
    for (int i = 0; i < m; i++)
      Glob_to_local(RowNumber(i)) = i;

    // assembling matrix B with interactions coming from other processors
    AddReceivedInteractions(comm, B, EntierToRecv, FloatToRecv, nrecv_int,
                            EntierToSend, FloatToSend, nsend_int, Glob_to_local,
                            OverlappedCol, OverlapProcNumber, procB, reorder);

    // exchanging datas
    SendAndReceiveDistributed(comm, nsend_int, EntierToSend, FloatToSend,
                              nrecv_int, EntierToRecv, FloatToRecv);

    // assembling matrix B with last interactions coming from other processors
    AddReceivedInteractions(comm, B, EntierToRecv, FloatToRecv, nrecv_int,
                            EntierToSend, FloatToSend, nsend_int, Glob_to_local,
                            OverlappedCol, OverlapProcNumber, procB, reorder);
    
    /****************************
     * Reordering of the matrix *
     ****************************/

    if (reorder)
      {
        // in this section, global rows are renumbered such that 
        // each processor has consecutive row numbers (mandatory for SuperLU)
        // processor 0 will be affected with rows 0..nloc0
        // processor 1 with rows nloc0 ... nloc0 + nloc1
        // ...
        IVect OverlapLocalNumber(OverlapRowNumber.GetM());
        IVect offset_global(nb_proc+1);
        
        IVect nb_col_sent(nb_proc);
        IVect nb_row_overlap(nb_proc);
        
        Vector<IVect> col_number_sorted(nb_proc);
        Vector<IVect> col_number_to_send(nb_proc);
        nb_col_sent.Zero();

        // counting the number of rows to send
        nb_row_overlap.Zero();
        for (int j = 0; j < OverlapRowNumber.GetM(); j++)
          {
            int i = OverlapProcNumber(j);
            if (i != rank)
              nb_row_overlap(i)++;
          }

        Vector<IVect> row_send_overlap(nb_proc);
        for (int i = 0; i < nb_proc; i++)
          if (nb_row_overlap(i) > 0)
            row_send_overlap(i).Reallocate(nb_row_overlap(i));
        
        nb_row_overlap.Zero();
        for (int j = 0; j < OverlapRowNumber.GetM(); j++)
          {
            int i = OverlapProcNumber(j);
            if (i != rank)
              row_send_overlap(i)(nb_row_overlap(i)++) = RowNumber(OverlapRowNumber(j));
          }
        
        // counting the number of columns to send
        for (int j = 0; j < B.GetM(); j++)
          for (int k = 0; k < B.GetRowSize(j); k++)
            {
              int i = procB(j)(k);
              if (i != rank)
                nb_col_sent(i)++;
            }
        
        for (int i = 0; i < nb_proc; i++)
          if (i != rank)
            col_number_sorted(i).Reallocate(nb_col_sent(i));
        
        // storing all the column numbers with distant processors
        nb_col_sent.Zero();
        for (int j = 0; j < B.GetM(); j++)
          for (int k = 0; k < B.GetRowSize(j); k++)
            {
              int i = procB(j)(k);
              if (i != rank)
                col_number_sorted(i)(nb_col_sent(i)++) = B.Index(j, k);
            }

        // duplicates are removed in order to send a few numbers
        for (int i = 0; i < nb_proc; i++)
          if (i != rank)
            {
              IVect permut(nb_col_sent(i)), col_number_sort(col_number_sorted(i));
              permut.Fill();
              Sort(nb_col_sent(i), col_number_sort, permut);
              
              // counting the number of unique numbers
              int prec = -1; nb_col_sent(i) = 0;
              for (int j = 0; j < col_number_sort.GetM(); j++)
                {
                  if (col_number_sort(j) != prec)
                    nb_col_sent(i)++;
                  
                  prec = col_number_sort(j);
                }
              
              // filling col_number_to_send
              col_number_to_send(i).Reallocate(nb_col_sent(i));
              nb_col_sent(i) = 0;
              prec = -1;
              for (int j = 0; j < col_number_sort.GetM(); j++)
                {
                  if (col_number_sort(j) != prec)
                    {
                      col_number_to_send(i)(nb_col_sent(i)) = col_number_sort(j);
                      nb_col_sent(i)++;
                    }
                  
                  col_number_sorted(i)(permut(j)) = nb_col_sent(i)-1;
                  prec = col_number_sort(j);
                }               
            }
        
        // allocating the array EntierToSend
        nsend_int.Zero();
        for (int i = 0; i < nb_proc; i++)
          if ((i != rank) && (nb_col_sent(i)+nb_row_overlap(i) > 0))
            {
              // column numbers will be sent
              nsend_int(i) = 2+nb_col_sent(i) + nb_row_overlap(i);
              EntierToSend(i).Reallocate(nsend_int(i));
              EntierToSend(i)(0) = 0;
              EntierToSend(i)(1) = nb_col_sent(i) + nb_row_overlap(i);
              nsend_int(i) = 2;
            }
        
        // storing columns numbers associated with processor i
        for (int i = 0; i < nb_proc; i++)
          {
            for (int j = 0; j < nb_row_overlap(i); j++)
              EntierToSend(i)(nsend_int(i)++) = row_send_overlap(i)(j);
            
            for (int j = 0; j < nb_col_sent(i); j++)
              EntierToSend(i)(nsend_int(i)++) = col_number_to_send(i)(j);
          }

        // exchanging datas
        SendAndReceiveDistributed(comm, nsend_int, EntierToSend, FloatToSend,
                                  nrecv_int, EntierToRecv, FloatToRecv);
        
        IVect nsend_intB(nb_proc), nrecv_intB(nb_proc);
        nsend_intB.Zero(); nrecv_intB.Zero();
        Vector<IVect> EntierToSendB(nb_proc), EntierToRecvB(nb_proc);
        // detecting if there are some numbers that do not belong to the current processor
        for (int i = 0; i < nb_proc; i++)
          if ((i != rank) && (nrecv_int(i) > 0))
            {
              int nb_col = EntierToRecv(i)(1);
              for (int j = 0; j < nb_col; j++)
                {
                  int iglob = EntierToRecv(i)(2+j);
                  int irow = Glob_to_local(iglob);
                  if (inv_local_row_numbers(irow) == -1)
                    {
                      int p = OverlappedCol(irow);
                      int nproc = OverlapProcNumber(p);
                      if (nsend_intB(nproc) == 0)
                        {
                          nsend_intB(nproc) = 3;
                          EntierToSendB(nproc).Reallocate(3);
                          EntierToSendB(nproc)(0) = 0;
                          EntierToSendB(nproc)(1) = 1;
                          EntierToSendB(nproc)(2) = iglob;
                        }
                      else
                        {
                          nsend_intB(nproc)++;
                          EntierToSendB(nproc)(1)++;
                          EntierToSendB(nproc).PushBack(iglob);
                        }
                    }
                }
            }
        
        // exchanging non-original dofs
        SendAndReceiveDistributed(comm, nsend_intB, EntierToSendB, FloatToSend,
                                  nrecv_intB, EntierToRecvB, FloatToRecv);        
        
        nsend_intB.Zero();
        for (int i = 0; i < nb_proc; i++)
          if ((i != rank) && (nrecv_intB(i) > 0))
            {
              int nb_col = EntierToRecvB(i)(1);
              if (nb_col > 0)
                {
                  nsend_intB(i) = nb_col+2;
                  EntierToSendB(i).Reallocate(nb_col+2);
                  EntierToSendB(i)(0) = 0;
                  EntierToSendB(i)(1) = nb_col;
                  
                  for (int j = 0; j < nb_col; j++)
                    {
                      int iglob = EntierToRecvB(i)(2+j);
                      int irow = Glob_to_local(iglob);
                      if (inv_local_row_numbers(irow) == -1)
                        {
                          cout << "impossible case" << endl;
                          abort();
                        }
                      else
                        EntierToSendB(i)(2+j) = inv_local_row_numbers(irow);
                    }
                }
            }
        
        // returning the local numbers of non-original dofs
        SendAndReceiveDistributed(comm, nsend_intB, EntierToSendB, FloatToSend,
                                  nrecv_intB, EntierToRecvB, FloatToRecv);        

        // filling local row numbers that need to be sent back
        nsend_intB.Zero();
        for (int i = 0; i < nb_proc; i++)
          if ((i != rank) && (nrecv_int(i) > 0))
            {
              int nb_col = EntierToRecv(i)(1);
              if (nb_col > 0)
                {
                  nsend_int(i) = 2*nb_col+2;
                  EntierToSend(i).Reallocate(2*nb_col+2);
                  EntierToSend(i)(0) = 0;
                  EntierToSend(i)(1) = 2*nb_col;
                  
                  for (int j = 0; j < nb_col; j++)
                    {
                      int iglob = EntierToRecv(i)(2+j);
                      int irow = Glob_to_local(iglob);
                      if (inv_local_row_numbers(irow) == -1)
                        {
                          int p = OverlappedCol(irow);
                          int nproc = OverlapProcNumber(p);
                          int num = EntierToRecvB(nproc)(2+nsend_intB(nproc));
                          nsend_intB(nproc)++;
                          
                          EntierToSend(i)(2+2*j) = num;
                          EntierToSend(i)(3+2*j) = nproc;
                        }
                      else
                        {
                          EntierToSend(i)(2+2*j) = inv_local_row_numbers(irow);
                          EntierToSend(i)(3+2*j) = rank;
                        }
                    }
                }
            }
        
        // exchanging datas
        SendAndReceiveDistributed(comm, nsend_int, EntierToSend, FloatToSend,
                                  nrecv_int, EntierToRecv, FloatToRecv);
        
        // receiving local numbers
        Vector<IVect> proc_number_to_send(nb_proc);
        for (int i = 0; i < nb_proc; i++)
          {     
            if (EntierToRecv(i).GetM() > 0)
              {	
                nrecv_int(i) = 2;
                proc_number_to_send(i).Reallocate(nb_col_sent(i));
                for (int j = 0; j < nb_row_overlap(i); j++)
                  {
                    row_send_overlap(i)(j) = EntierToRecv(i)(nrecv_int(i));
                    if (EntierToRecv(i)(nrecv_int(i)+1) != i)
                      {
                        cout << "Impossible case" << endl;
                        abort();
                      }
                    
                    nrecv_int(i) += 2;
                  }

                for (int j = 0; j < nb_col_sent(i); j++)
                  {
                    col_number_to_send(i)(j) = EntierToRecv(i)(nrecv_int(i));
                    proc_number_to_send(i)(j) = EntierToRecv(i)(nrecv_int(i)+1);
                    nrecv_int(i) += 2;
                  }
              }        
          }
        
        // filling OverlapLocalNumber
        nb_row_overlap.Zero();
        
        OverlapLocalNumber.Fill(-1);
        for (int j = 0; j < OverlapRowNumber.GetM(); j++)
          {
            int i = OverlapProcNumber(j);
            if (i != rank)
              {
                OverlapLocalNumber(j) = row_send_overlap(i)(nb_row_overlap(i));
                nb_row_overlap(i)++;
              }
          }
        
        // now in order to compute the global row numbers for any processor
        // we need to retrieve the offsets (i.e. nloc cumulated)
        offset_global.Zero();
        
        MPI_Allgather(&nloc, 1, MPI_INTEGER, &offset_global(1), 1, MPI_INTEGER, comm);
        
        for (int i = 1; i < nb_proc; i++)
          offset_global(i+1) += offset_global(i);
       
        // RowNumber is modified
        nrow = 0;
        for (int i = 0; i < m; i++)
          {
            if (OverlappedCol(i) == -1)
              {
                RowNumber(i) = offset_global(rank) + nrow;
                nrow++;
              }
            else
              RowNumber(i) = -1;
          }
        
        // then numbers in B are modified with the new numbering
        nb_col_sent.Zero();
        for (int j = 0; j < m; j++)
          if (OverlappedCol(j) == -1)
            {
              for (int k = 0; k < B.GetRowSize(j); k++)
                {
                  int jglob = B.Index(j, k);
                  int i = procB(j)(k);
                  int ireal = i;
                  int iloc = -1;
                  if (i == rank)
                    {
                      int irow = Glob_to_local(jglob);
                      if (OverlappedCol(irow) == -1)
                        iloc = inv_local_row_numbers(irow);
                      else
                        {
                          int p = OverlappedCol(irow);
                          i = OverlapProcNumber(p);
                          ireal = i;
                          iloc = OverlapLocalNumber(p);
                        }
                    }
                  else
                    {
                      int ilocC = col_number_sorted(i)(nb_col_sent(i));
                      ireal = proc_number_to_send(i)(ilocC);
                      iloc = col_number_to_send(i)(ilocC);
                      nb_col_sent(i)++;
                    }
                  
                  if (iloc >= 0)
                    {
                      B.Index(j, k) = offset_global(ireal) + iloc;
                    }
                  else
                    {
                      cout << "Impossible case" << endl;
                      abort();
                    }
                }
            }    
      }

    Glob_to_local.Clear();
    
    nrow = 0;
    row_numbers.Reallocate(nloc);
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
	{
	  row_numbers(nrow) = RowNumber(i);
	  nrow++;
	}    
  }

  
  //! Fills PtrA, IndA and ValA from values contained in B
  template<class T> template<class Tint0, class Tint1>
  void DistributedMatrix_Base<T>
  ::ConvertToCSR(Matrix<T, General, ArrayRowSparse>& B, IVect& OverlappedCol,
		 Vector<Tint0>& PtrA, Vector<Tint1>& IndA, Vector<T>& ValA)
  {
    int m = B.GetM();
    int nloc = 0;
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
	nloc++;
    
    /***************************
     * Final conversion in CSR *
     ***************************/
    
    // now we convert Bh to RowSparse while removing overlapped rows
    PtrA.Reallocate(nloc+1);
    int nrow = 0; long nnz = 0;
    PtrA(nrow) = 0;
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
	{
	  PtrA(nrow+1) = PtrA(nrow) + B.GetRowSize(i);
	  nrow++;
          nnz += B.GetRowSize(i);
	}
    
    IndA.Reallocate(nnz);
    ValA.Reallocate(nnz); nrow = 0; nnz = 0;
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
	{
	  for (int j = 0; j < B.GetRowSize(i); j++)
	    {
              IndA(nnz) = B.Index(i, j);
	      ValA(nnz) = B.Value(i, j);
	      nnz++;
	    }
          
	  nrow++;
	}    
  }


  //! Method called by AssembleDistributed  
  template<class T>
  void DistributedMatrix_Base<T>
  ::AssembleParallel(Matrix<T, General, ArrayColSparse>& B, Vector<IVect>& procB,
                     General& prop, IVect& col_numbers, IVect& local_col_numbers,
                     IVect& OverlappedCol, bool sym_pattern, bool reorder)
  {
    int n = this->GetLocalN();
    int m = this->GetGlobalM();
    
    const MPI_Comm& comm = this->GetCommunicator();  
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank; MPI_Comm_rank(comm, &rank);
    
    IVect RowNumber = this->GetGlobalRowNumber();
    const IVect& OverlapRowNumber = this->GetOverlapRowNumber();
    const IVect& OverlapProcNumber = this->GetOverlapProcNumber();

    // constructing array to know if a column is overlapped 
    // (treated by another processor)
    OverlappedCol.Reallocate(n); OverlappedCol.Fill(-1);
    for (int i = 0; i < OverlapRowNumber.GetM(); i++)
      OverlappedCol(OverlapRowNumber(i)) = i;

    /*********************************
     * Parallel assembling of matrix *
     *********************************/
    
    // we send to each processor additional entries due to overlapping
    // distant columns or distant rows (because of symmetrisation of patterns
    IVect nsend_int(nb_proc), nsend_float(nb_proc), nb_col_sent(nb_proc);
    Vector<IVect> EntierToSend(nb_proc);
    Vector<Vector<T> > FloatToSend(nb_proc);


    // counting the size of arrays to send
    // nsend_int : the number of integers
    // nsend_float : the number of floats
    for (int i = 0; i < nb_proc; i++)
      {
        if (i != rank)
          nsend_int(i) = 2;
        else
          nsend_int(i) = 0;
        
        nsend_float(i) = 0;
        nb_col_sent(i) = 0;
      }

    // overlapped columns
    for (int j = 0; j < OverlapRowNumber.GetM(); j++)
      {
        int i = OverlapProcNumber(j);
        if (i != rank)
          {
            nsend_int(i) += 2;
            nb_col_sent(i)++;
            int irow = OverlapRowNumber(j);
            nsend_int(i) += B.GetColumnSize(irow);
            nsend_float(i) += B.GetColumnSize(irow);
            if (reorder)
              nsend_int(i) += B.GetColumnSize(irow);
          }
      }
    
    // distant rows
    if (sym_pattern)
      for (int j = 0; j < this->dist_row.GetM(); j++)
        for (int k = 0; k < this->dist_row(j).GetM(); k++)
          {
            int i = this->proc_row(j)(k);
            if (i != rank)
              {
                nb_col_sent(i)++;
                nsend_int(i) += 3;
                nsend_float(i)++;
                if (reorder)
                  nsend_int(i)++;
              }
          }
    
    // distant columns
    for (int j = 0; j < this->dist_col.GetM(); j++)
      for (int k = 0; k < this->dist_col(j).GetM(); k++)
        {
          int i = this->proc_col(j)(k);
          if (i != rank)
            {
              nb_col_sent(i)++;
              nsend_int(i) += 3;
              nsend_float(i)++;
              if (reorder)
                nsend_int(i)++;
            }
        }


    // allocating arrays EntierToSend and FloatToSend
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        { 
          if (nb_col_sent(i) == 0)
            nsend_int(i) = 0;
          
          if (nb_col_sent(i) > 0)
            {
              EntierToSend(i).Reallocate(nsend_int(i));
              FloatToSend(i).Reallocate(nsend_float(i));
              EntierToSend(i)(0) = nsend_float(i);
              EntierToSend(i)(1) = nb_col_sent(i);
              nsend_int(i) = 2; nsend_float(i) = 0;
            }
        }

    // then arrays EntierToSend and FloatToSend are filled

    // storing values and indices of a column shared with processor i
    // processor i is the owner of this column
    for (int j = 0; j < OverlapRowNumber.GetM(); j++)
      {
        int i = OverlapProcNumber(j);
        if (i != rank)
          {
            int irow = OverlapRowNumber(j);
            EntierToSend(i)(nsend_int(i)++) = RowNumber(irow);
            EntierToSend(i)(nsend_int(i)++) = B.GetColumnSize(irow);
            for (int k = 0; k < B.GetColumnSize(irow); k++)
              {
                EntierToSend(i)(nsend_int(i)++) = B.Index(irow, k);
                if (reorder)
                  EntierToSend(i)(nsend_int(i)++) = procB(irow)(k);
                
                FloatToSend(i)(nsend_float(i)++) = B.Value(irow, k);
              }
            
            // the corresponding values of B are cleared
            // they are no longer needed since they will be present in the distant processor
            // after the exchange of datas
            B.ClearColumn(irow);
            if (reorder)
              procB(irow).Clear();
          }
      }

    // storing values to enforce a symmetric pattern
    if (sym_pattern)
      for (int j = 0; j < this->dist_row.GetM(); j++)
        for (int k = 0; k < this->dist_row(j).GetM(); k++)
          {
            int i = this->proc_row(j)(k);
            if (i != rank)
              {
                int irow = this->dist_row(j).Index(k);
                if (this->local_number_distant_values)
                  irow = this->global_row_to_recv(irow);
                
                EntierToSend(i)(nsend_int(i)++) = irow;
                EntierToSend(i)(nsend_int(i)++) = 1;
                EntierToSend(i)(nsend_int(i)++) = RowNumber(j);
                if (reorder)
                  EntierToSend(i)(nsend_int(i)++) = rank;
                
                FloatToSend(i)(nsend_float(i)++) = 0;
              }
          }
            
    // storing values of row associated with processor rank
    for (int j = 0; j < this->dist_col.GetM(); j++)
      for (int k = 0; k < this->dist_col(j).GetM(); k++)
        {
          int i = this->proc_col(j)(k);
          if (i != rank)
            {
              int irow = this->dist_col(j).Index(k);
              if (this->local_number_distant_values)
                irow = this->global_col_to_recv(irow);
              
              EntierToSend(i)(nsend_int(i)++) = irow;
              EntierToSend(i)(nsend_int(i)++) = 1;
              EntierToSend(i)(nsend_int(i)++) = RowNumber(j);
              if (reorder)
                EntierToSend(i)(nsend_int(i)++) = rank;
              
              FloatToSend(i)(nsend_float(i)++)
                = this->dist_col(j).Value(k);
            }
        }

    // now the initial matrix can be cleared    
    this->Clear();
    
    // Datas for receiving EntierToSend and FloatToSend
    IVect nrecv_int(nb_proc);
    Vector<IVect> EntierToRecv(nb_proc);
    Vector<Vector<T> > FloatToRecv(nb_proc);

    // exchanging datas
    SendAndReceiveDistributed(comm, nsend_int, EntierToSend, FloatToSend,
                              nrecv_int, EntierToRecv, FloatToRecv);

    // constructing local column numbers 
    int nloc = n - OverlapRowNumber.GetM();
    local_col_numbers.Reallocate(nloc);
    int ncol = 0;
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
	local_col_numbers(ncol++) = i;
    
    // index array to obtain local numbers in array local_col_numbers
    // from local numbers of the matrix
    IVect inv_local_col_numbers(n);
    inv_local_col_numbers.Fill(-1);
    ncol = 0;
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
        inv_local_col_numbers(i) = ncol++;
    
    // global to local conversion
    IVect Glob_to_local(m);
    Glob_to_local.Fill(-1);
    for (int i = 0; i < n; i++)
      Glob_to_local(RowNumber(i)) = i;

    // assembling matrix B with interactions coming from other processors
    AddReceivedInteractions(comm, B, EntierToRecv, FloatToRecv, nrecv_int,
                            EntierToSend, FloatToSend, nsend_int, Glob_to_local,
                            OverlappedCol, OverlapProcNumber, procB, reorder);

    // exchanging datas
    SendAndReceiveDistributed(comm, nsend_int, EntierToSend, FloatToSend,
                              nrecv_int, EntierToRecv, FloatToRecv);

    // assembling matrix B with last interactions coming from other processors
    AddReceivedInteractions(comm, B, EntierToRecv, FloatToRecv, nrecv_int,
                            EntierToSend, FloatToSend, nsend_int, Glob_to_local,
                            OverlappedCol, OverlapProcNumber, procB, reorder);

    /****************************
     * Reordering of the matrix *
     ****************************/

    if (reorder)
      {
        // in this section, global rows are renumbered such that 
        // each processor has consecutive row numbers (mandatory for SuperLU)
        // processor 0 will be affected with rows 0..nloc0
        // processor 1 with rows nloc0 ... nloc0 + nloc1
        // ...
        IVect OverlapLocalNumber(OverlapRowNumber.GetM());
        IVect offset_global(nb_proc+1);
        
        IVect nb_col_overlap(nb_proc);
        
        Vector<IVect> col_number_sorted(nb_proc);
        Vector<IVect> col_number_to_send(nb_proc);
        nb_col_sent.Zero();

        // counting the number of columns to send
        nb_col_overlap.Zero();
        for (int j = 0; j < OverlapRowNumber.GetM(); j++)
          {
            int i = OverlapProcNumber(j);
            if (i != rank)
              nb_col_overlap(i)++;
          }

        Vector<IVect> col_send_overlap(nb_proc);
        for (int i = 0; i < nb_proc; i++)
          if (nb_col_overlap(i) > 0)
            col_send_overlap(i).Reallocate(nb_col_overlap(i));
        
        nb_col_overlap.Zero();
        for (int j = 0; j < OverlapRowNumber.GetM(); j++)
          {
            int i = OverlapProcNumber(j);
            if (i != rank)
              col_send_overlap(i)(nb_col_overlap(i)++) = RowNumber(OverlapRowNumber(j));
          }
        
        // counting the number of columns to send
        for (int j = 0; j < B.GetN(); j++)
          for (int k = 0; k < B.GetColumnSize(j); k++)
            {
              int i = procB(j)(k);
              if (i != rank)
                nb_col_sent(i)++;
            }
        
        for (int i = 0; i < nb_proc; i++)
          if (i != rank)
            col_number_sorted(i).Reallocate(nb_col_sent(i));
        
        // storing all the column numbers with distant processors
        nb_col_sent.Zero();
        for (int j = 0; j < B.GetN(); j++)
          for (int k = 0; k < B.GetColumnSize(j); k++)
            {
              int i = procB(j)(k);
              if (i != rank)
                col_number_sorted(i)(nb_col_sent(i)++) = B.Index(j, k);
            }

        // duplicates are removed in order to send a few numbers
        for (int i = 0; i < nb_proc; i++)
          if (i != rank)
            {
              IVect permut(nb_col_sent(i)), col_number_sort(col_number_sorted(i));
              permut.Fill();
              Sort(nb_col_sent(i), col_number_sort, permut);
              
              // counting the number of unique numbers
              int prec = -1; nb_col_sent(i) = 0;
              for (int j = 0; j < col_number_sort.GetM(); j++)
                {
                  if (col_number_sort(j) != prec)
                    nb_col_sent(i)++;
                  
                  prec = col_number_sort(j);
                }
              
              // filling col_number_to_send
              col_number_to_send(i).Reallocate(nb_col_sent(i));
              nb_col_sent(i) = 0;
              prec = -1;
              for (int j = 0; j < col_number_sort.GetM(); j++)
                {
                  if (col_number_sort(j) != prec)
                    {
                      col_number_to_send(i)(nb_col_sent(i)) = col_number_sort(j);
                      nb_col_sent(i)++;
                    }
                  
                  col_number_sorted(i)(permut(j)) = nb_col_sent(i)-1;
                  prec = col_number_sort(j);
                }               
            }
        
        // allocating the array EntierToSend
        nsend_int.Zero();
        for (int i = 0; i < nb_proc; i++)
          if ((i != rank) && (nb_col_sent(i)+nb_col_overlap(i) > 0))
            {
              // column numbers will be sent
              nsend_int(i) = 2+nb_col_sent(i) + nb_col_overlap(i);
              EntierToSend(i).Reallocate(nsend_int(i));
              EntierToSend(i)(0) = 0;
              EntierToSend(i)(1) = nb_col_sent(i) + nb_col_overlap(i);
              nsend_int(i) = 2;
            }
        
        // storing columns numbers associated with processor i
        for (int i = 0; i < nb_proc; i++)
          {
            for (int j = 0; j < nb_col_overlap(i); j++)
              EntierToSend(i)(nsend_int(i)++) = col_send_overlap(i)(j);
            
            for (int j = 0; j < nb_col_sent(i); j++)
              EntierToSend(i)(nsend_int(i)++) = col_number_to_send(i)(j);
          }

        // exchanging datas
        SendAndReceiveDistributed(comm, nsend_int, EntierToSend, FloatToSend,
                                  nrecv_int, EntierToRecv, FloatToRecv);
        
        IVect nsend_intB(nb_proc), nrecv_intB(nb_proc);
        nsend_intB.Zero(); nrecv_intB.Zero();
        Vector<IVect> EntierToSendB(nb_proc), EntierToRecvB(nb_proc);
        // detecting if there are some numbers that do not belong to the current processor
        for (int i = 0; i < nb_proc; i++)
          if ((i != rank) && (nrecv_int(i) > 0))
            {
              int nb_col = EntierToRecv(i)(1);
              for (int j = 0; j < nb_col; j++)
                {
                  int iglob = EntierToRecv(i)(2+j);
                  int irow = Glob_to_local(iglob);
                  if (inv_local_col_numbers(irow) == -1)
                    {
                      int p = OverlappedCol(irow);
                      int nproc = OverlapProcNumber(p);
                      if (nsend_intB(nproc) == 0)
                        {
                          nsend_intB(nproc) = 3;
                          EntierToSendB(nproc).Reallocate(3);
                          EntierToSendB(nproc)(0) = 0;
                          EntierToSendB(nproc)(1) = 1;
                          EntierToSendB(nproc)(2) = iglob;
                        }
                      else
                        {
                          nsend_intB(nproc)++;
                          EntierToSendB(nproc)(1)++;
                          EntierToSendB(nproc).PushBack(iglob);
                        }
                    }
                }
            }
        
        // exchanging non-original dofs
        SendAndReceiveDistributed(comm, nsend_intB, EntierToSendB, FloatToSend,
                                  nrecv_intB, EntierToRecvB, FloatToRecv);        
        
        nsend_intB.Zero();
        for (int i = 0; i < nb_proc; i++)
          if ((i != rank) && (nrecv_intB(i) > 0))
            {
              int nb_col = EntierToRecvB(i)(1);
              if (nb_col > 0)
                {
                  nsend_intB(i) = nb_col+2;
                  EntierToSendB(i).Reallocate(nb_col+2);
                  EntierToSendB(i)(0) = 0;
                  EntierToSendB(i)(1) = nb_col;
                  
                  for (int j = 0; j < nb_col; j++)
                    {
                      int iglob = EntierToRecvB(i)(2+j);
                      int irow = Glob_to_local(iglob);
                      if (inv_local_col_numbers(irow) == -1)
                        {
                          cout << "impossible case" << endl;
                          abort();
                        }
                      else
                        EntierToSendB(i)(2+j) = inv_local_col_numbers(irow);
                    }
                }
            }
        
        // returning the local numbers of non-original dofs
        SendAndReceiveDistributed(comm, nsend_intB, EntierToSendB, FloatToSend,
                                  nrecv_intB, EntierToRecvB, FloatToRecv);        

        // filling local row numbers that need to be sent back
        nsend_intB.Zero();
        for (int i = 0; i < nb_proc; i++)
          if ((i != rank) && (nrecv_int(i) > 0))
            {
              int nb_col = EntierToRecv(i)(1);
              if (nb_col > 0)
                {
                  nsend_int(i) = 2*nb_col+2;
                  EntierToSend(i).Reallocate(2*nb_col+2);
                  EntierToSend(i)(0) = 0;
                  EntierToSend(i)(1) = 2*nb_col;
                  
                  for (int j = 0; j < nb_col; j++)
                    {
                      int iglob = EntierToRecv(i)(2+j);
                      int irow = Glob_to_local(iglob);
                      if (inv_local_col_numbers(irow) == -1)
                        {
                          int p = OverlappedCol(irow);
                          int nproc = OverlapProcNumber(p);
                          int num = EntierToRecvB(nproc)(2+nsend_intB(nproc));
                          nsend_intB(nproc)++;
                          
                          EntierToSend(i)(2+2*j) = num;
                          EntierToSend(i)(3+2*j) = nproc;
                        }
                      else
                        {
                          EntierToSend(i)(2+2*j) = inv_local_col_numbers(irow);
                          EntierToSend(i)(3+2*j) = rank;
                        }
                    }
                }
            }
        
        // exchanging datas
        SendAndReceiveDistributed(comm, nsend_int, EntierToSend, FloatToSend,
                                  nrecv_int, EntierToRecv, FloatToRecv);
        
        // receiving local numbers
        Vector<IVect> proc_number_to_send(nb_proc);
        for (int i = 0; i < nb_proc; i++)
          {     
            if (EntierToRecv(i).GetM() > 0)
              {	
                nrecv_int(i) = 2;
                proc_number_to_send(i).Reallocate(nb_col_sent(i));
                for (int j = 0; j < nb_col_overlap(i); j++)
                  {
                    col_send_overlap(i)(j) = EntierToRecv(i)(nrecv_int(i));
                    if (EntierToRecv(i)(nrecv_int(i)+1) != i)
                      {
                        cout << "Impossible case" << endl;
                        abort();
                      }
                    
                    nrecv_int(i) += 2;
                  }

                for (int j = 0; j < nb_col_sent(i); j++)
                  {
                    col_number_to_send(i)(j) = EntierToRecv(i)(nrecv_int(i));
                    proc_number_to_send(i)(j) = EntierToRecv(i)(nrecv_int(i)+1);
                    nrecv_int(i) += 2;
                  }
              }        
          }
        
        // filling OverlapLocalNumber
        nb_col_overlap.Zero();
        
        OverlapLocalNumber.Fill(-1);
        for (int j = 0; j < OverlapRowNumber.GetM(); j++)
          {
            int i = OverlapProcNumber(j);
            if (i != rank)
              {
                OverlapLocalNumber(j) = col_send_overlap(i)(nb_col_overlap(i));
                nb_col_overlap(i)++;
              }
          }
        
        // now in order to compute the global row numbers for any processor
        // we need to retrieve the offsets (i.e. nloc cumulated)
        offset_global.Zero();
        
        MPI_Allgather(&nloc, 1, MPI_INTEGER, &offset_global(1), 1, MPI_INTEGER, comm);
        
        for (int i = 1; i < nb_proc; i++)
          offset_global(i+1) += offset_global(i);
       
        // RowNumber is modified
        ncol = 0;
        for (int i = 0; i < n; i++)
          {
            if (OverlappedCol(i) == -1)
              {
                RowNumber(i) = offset_global(rank) + ncol;
                ncol++;
              }
            else
              RowNumber(i) = -1;
          }
        
        // then numbers in B are modified with the new numbering
        nb_col_sent.Zero();
        for (int j = 0; j < n; j++)
          if (OverlappedCol(j) == -1)
            {
              for (int k = 0; k < B.GetColumnSize(j); k++)
                {
                  int jglob = B.Index(j, k);
                  int i = procB(j)(k);
                  int ireal = i;
                  int iloc = -1;
                  if (i == rank)
                    {
                      int irow = Glob_to_local(jglob);
                      if (OverlappedCol(irow) == -1)
                        iloc = inv_local_col_numbers(irow);
                      else
                        {
                          int p = OverlappedCol(irow);
                          i = OverlapProcNumber(p);
                          ireal = i;
                          iloc = OverlapLocalNumber(p);
                        }
                    }
                  else
                    {
                      int ilocC = col_number_sorted(i)(nb_col_sent(i));
                      ireal = proc_number_to_send(i)(ilocC);
                      iloc = col_number_to_send(i)(ilocC);
                      nb_col_sent(i)++;
                    }
                  
                  if (iloc >= 0)
                    {
                      B.Index(j, k) = offset_global(ireal) + iloc;
                    }
                  else
                    {
                      cout << "Impossible case" << endl;
                      abort();
                    }
                }
              
              B.AssembleColumn(j);
            }    
      }

    Glob_to_local.Clear();

    ncol = 0;
    col_numbers.Reallocate(nloc);
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
	{
	  col_numbers(ncol) = RowNumber(i);
	  ncol++;
	}
  }


  //! Fills PtrA, IndA and ValA from values contained in B
  template<class T> template<class Tint0, class Tint1>
  void DistributedMatrix_Base<T>
  ::ConvertToCSC(Matrix<T, General, ArrayColSparse>& B, IVect& OverlappedCol,
		 Vector<Tint0>& PtrA, Vector<Tint1>& IndA, Vector<T>& ValA)
  {  
    int n = B.GetN();
    int nloc = 0;
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
	nloc++;

    /***************************
     * Final conversion in CSC *
     ***************************/
    
    PtrA.Reallocate(nloc+1);
    int ncol = 0; long nnz = 0;
    PtrA(ncol) = 0;
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
        {
          PtrA(ncol+1) = PtrA(ncol) + B.GetColumnSize(i);
          ncol++;
          nnz += B.GetColumnSize(i);
        }
    
    IndA.Reallocate(nnz);
    ValA.Reallocate(nnz);
    ncol = 0; nnz = 0;
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
        {
          for (int j = 0; j < B.GetColumnSize(i); j++)
            {
              IndA(nnz) = B.Index(i, j);
              ValA(nnz) = B.Value(i, j);
              nnz++;
            }
          
          ncol++;
        }
  }
  

  //! grouping all the local rows of the matrix into CSR form
  /*!
    Creation of a sparse matrix B, that will contain the local rows
    of the current matrix. The column numbers are global.
    Values placed on non-local rows are ignored
  */
  template<class T>
  template<class T0, class Allocator0>
  void DistributedMatrix_Base<T>
  ::GetDistributedRows(Matrix<T0, General, ArrayRowSparse,
		       Allocator0>& B, Vector<IVect>& procB, bool sym) const
  {    
    int m = this->GetLocalM();
    int n = this->GetGlobalM();
    int rank; MPI_Comm_rank(this->comm_, &rank);
    
    bool retrieve_proc = false;
    if (procB.GetM() == m)
      retrieve_proc = true;
    
    // now, we are using global numbers
    // and removing lower part of the matrix if symmetric
    B.Resize(m, n);
    const IVect& RowNumber = this->GetGlobalRowNumber();
    for (int i = 0; i < m; i++)
      {
	int size_row = B.GetRowSize(i);
        size_row += dist_col(i).GetM();
	IVect index(size_row), proc_loc(size_row);
	Vector<T0> value(size_row);
	int nb = 0;
	int num_row = RowNumber(i);
	if (sym)
          {
            // local values
            for (int j = 0; j < B.GetRowSize(i); j++)
              {
                int num_col = RowNumber(B.Index(i, j));
                if (num_row <= num_col)
                  {
                    index(nb) = num_col;
                    value(nb) = B.Value(i, j);
                    proc_loc(nb) = rank;
                    nb++;
                  }
              }
            
            // distant values
            for (int j = 0; j < dist_col(i).GetM(); j++)
              {
                int num_col = dist_col(i).Index(j);
                if (this->local_number_distant_values)
                  num_col = global_col_to_recv(num_col);
                
                if (num_row <= num_col)
                  {
                    index(nb) = num_col;
                    value(nb) = dist_col(i).Value(j);
                    proc_loc(nb) = proc_col(i)(j);
                    nb++;
                  }
              }
          }
        else
          {
            // local values
            for (int j = 0; j < B.GetRowSize(i); j++)
              {
                int num_col = RowNumber(B.Index(i, j));
                index(nb) = num_col;
                value(nb) = B.Value(i, j);
                proc_loc(nb) = rank;
                nb++;
              }
            
            // distant values
            for (int j = 0; j < dist_col(i).GetM(); j++)
              {
                int num_col = dist_col(i).Index(j);
                if (this->local_number_distant_values)
                  num_col = global_col_to_recv(num_col);
                
                index(nb) = num_col;
                value(nb) = dist_col(i).Value(j);
                proc_loc(nb) = proc_col(i)(j);
                nb++;
              }
          }
        
        Sort(nb, index, value, proc_loc);
        size_row = 0;
        int prec = -1;
        for (int j = 0; j < nb; j++)
          {
            if (index(j) == prec)
              value(size_row-1) += value(j);
            else
              {
                index(size_row) = index(j);
                value(size_row) = value(j);
                proc_loc(size_row) = proc_loc(j);
                size_row++;
              }
            
            prec = index(j);
          }
        
	B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = index(j);
            B.Value(i, j) = value(j);
          }
        
        if (retrieve_proc)
          {
            procB(i).Reallocate(size_row);
            for (int j = 0; j < size_row; j++)
              procB(i)(j) = proc_loc(j);
          }
      }
  }
  
  
  //! grouping all the local columns of the matrix into CSC form
  /*!
    Creation of a sparse matrix B, that will contain the local columns
    of the current matrix. The row numbers are global.
    Values placed on non-local columns are ignored
  */
  template<class T>
  template<class T0, class Allocator0>
  void DistributedMatrix_Base<T>
  ::GetDistributedColumns(Matrix<T0, General, ArrayColSparse, Allocator0>& B,
                          Vector<IVect>& procB, Vector<long>& Ptr, IVect& Ind,
                          Vector<T0>& Val, bool sym_pattern) const
  {
    int m = this->GetGlobalM();
    int n = this->GetLocalN();
    int rank; MPI_Comm_rank(this->comm_, &rank);
    
    bool retrieve_proc = false;
    if (procB.GetM() == n)
      retrieve_proc = true;
    
    // for row numbers, we put global numbers and we add some distant entries
    // (i.e entries with local columns,
    // and null values by symmetry of local rows )
    B.Clear(); B.Reallocate(m, n);
    const IVect& RowNumber = this->GetGlobalRowNumber();
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
        size_col += dist_row(i).GetM();
        if (sym_pattern)
          size_col += dist_col(i).GetM();
        
	IVect index(size_col), proc_loc(size_col);
	Vector<T0> value(size_col);
	int nb = 0;
        // local values
        for (int j = Ptr(i); j < Ptr(i+1); j++)
          {
            index(nb) = RowNumber(Ind(j));
            value(nb) = Val(j);
            proc_loc(nb) = rank;
            nb++;
          }
        
        // distant values
        for (int j = 0; j < dist_row(i).GetM(); j++)
          {
            index(nb) = dist_row(i).Index(j);
            if (this->local_number_distant_values)    
              index(nb) = global_row_to_recv(index(nb));
            
            proc_loc(nb) = proc_row(i)(j);
            value(nb) = dist_row(i).Value(j);
            nb++;
          }
        
        // values due to symmetrisation of pattern
        if (sym_pattern)
          for (int j = 0; j < dist_col(i).GetM(); j++)
            {
              index(nb) = dist_col(i).Index(j);
              if (this->local_number_distant_values)    
                index(nb) = global_col_to_recv(index(nb));
            
              proc_loc(nb) = proc_col(i)(j);
              value(nb) = 0;
              nb++;
            }
        
        Sort(nb, index, value, proc_loc);
        size_col = 0;
        int prec = -1;
        for (int j = 0; j < nb; j++)
          {
            if (index(j) == prec)
              value(size_col-1) += value(j);
            else
              {
                index(size_col) = index(j);
                value(size_col) = value(j);
                proc_loc(size_col) = proc_loc(j);
                size_col++;
              }
            
            prec = index(j);
          }
        
	B.ReallocateColumn(i, size_col);
        for (int j = 0; j < size_col; j++)
          {
            B.Index(i, j) = index(j);
            B.Value(i, j) = value(j);
          }

        if (retrieve_proc)
          {
            procB(i).Reallocate(size_col);
            for (int j = 0; j < size_col; j++)
              procB(i)(j) = proc_loc(j);
          }
      }
  }


  // DistributedMatrix_Base //
  ////////////////////////////


  ///////////////////////
  // DistributedMatrix //


  //! Default constructor
  template<class T, class Prop, class Storage, class Allocator>
  DistributedMatrix<T, Prop, Storage, Allocator>::DistributedMatrix()
    : Matrix<T, Prop, Storage, Allocator>(),
      DistributedMatrix_Base<T>()
  {
  }

  
  //! Allocates a matrix of size i x j 
  /*!
    Here i and j are the local number of rows and colums
    they are different for each processor
   */
  template<class T, class Prop, class Storage, class Allocator>
  DistributedMatrix<T, Prop, Storage, Allocator>::DistributedMatrix(int i, int j)
    : Matrix<T, Prop, Storage, Allocator>(i, j),
      DistributedMatrix_Base<T>(i, j)
  {
  }
  

  //! changing the size of the local matrix, previous values are lost
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Reallocate(int m, int n)
  {
    Matrix<T, Prop, Storage, Allocator>::Reallocate(m, n);
    DistributedMatrix_Base<T>::ReallocateDist(m, n);
  }


  //! changing the size of the local matrix, previous values are kept
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Resize(int m, int n)
  {
    Matrix<T, Prop, Storage, Allocator>::Resize(m, n);
    DistributedMatrix_Base<T>::Resize(m, n);
  }

  
  //! Clears the matrix
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Clear()
  {
    Matrix<T, Prop, Storage, Allocator>::Clear();
    DistributedMatrix_Base<T>::Clear();
  }


  //! Clears the local matrix
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::ClearLocal()
  {
    Matrix<T, Prop, Storage, Allocator>::Clear();
  }


  //! equality *this = X 
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Prop2, class Storage2, class Allocator2>
  DistributedMatrix<T, Prop, Storage, Allocator>&
  DistributedMatrix<T, Prop, Storage, Allocator>::
  operator=(const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& X)
  {
    Seldon::Copy(static_cast<const Matrix<T2, Prop2, Storage2, Allocator2>& >(X),
		 static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this));

    DistributedMatrix_Base<T>::Copy(X);
    return *this;
  }


  //! multiplication by a scalar
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  DistributedMatrix<T, Prop, Storage, Allocator>&
  DistributedMatrix<T, Prop, Storage, Allocator>::operator *=(const T0& x)
  {
    Mlt(x, static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this));
    
    static_cast<DistributedMatrix_Base<T>& >(*this) *= x;
    return *this;
  }


  //! returns the size of memory used to store the matrix
  template<class T, class Prop, class Storage, class Allocator>
  size_t DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetMemorySize() const
  {
    size_t taille = Matrix<T, Prop, Storage, Allocator>::GetMemorySize();
    taille += DistributedMatrix_Base<T>::GetMemorySize();
    return taille;
  }


  //! returns the number of non-zero entries stored in the matrix
  /*!
    It returns the number of non-zero entries stored on the current processor.
    To obtain an "estimation" of the global number of non-zero entries,
    MPI_Reduce should be called
   */
  template<class T, class Prop, class Storage, class Allocator>
  long DistributedMatrix<T, Prop, Storage, Allocator>::GetNonZeros() const
  {
    long nnz = Matrix<T, Prop, Storage, Allocator>::GetNonZeros();
    nnz += DistributedMatrix_Base<T>::GetNonZeros();
    return nnz;
  }
  

  //! returns the number of elements stored
  /*!
    It returns the number of elements (of type T)
    stored on the current processor.
    To obtain an "estimation" of the global number of elements
    MPI_Reduce should be called
  */
  template<class T, class Prop, class Storage, class Allocator>
  long DistributedMatrix<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    long size = Matrix<T, Prop, Storage, Allocator>::GetDataSize();
    size += DistributedMatrix_Base<T>::GetDataSize();
    return size;
  }
  

  //! removes small entries of the matrix
  /*!
    The result of this function will be different from a sequential
    execution since values stored on different processors
    may add up and be eliminated or not
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::RemoveSmallEntry(const T0& epsilon)
  {
    Seldon::RemoveSmallEntry(static_cast<Matrix<T, Prop, Storage,
			     Allocator>& >(*this),
                             epsilon);

    DistributedMatrix_Base<T>::RemoveSmallEntry(epsilon);
  }


  //! sets matrix to identity matrix
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::SetIdentity()
  {
    Matrix<T, Prop, Storage, Allocator>::SetIdentity();

    MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;

    // setting to 0 diagonal of overlapped rows
    const IVect& overlap = this->GetOverlapRowNumber();
    T zero; SetComplexZero(zero);
    for (int i = 0; i < overlap.GetM(); i++)
      this->Set(overlap(i), overlap(i), zero);
    
    DistributedMatrix_Base<T>::SetIdentity();
  }


  //! sets values of non-zero entries to 0
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Zero()
  {
    Matrix<T, Prop, Storage, Allocator>::Zero();
    DistributedMatrix_Base<T>::Zero();
  }


  //! sets values of non-zero entries to 0, 1, 2, etc
  /*!
    The result will be different from a sequential execution
    since values can be stored in different processors
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Fill()
  {
    Matrix<T, Prop, Storage, Allocator>::Fill();
    DistributedMatrix_Base<T>::Fill();
  }


  //! sets values of non-zero entries to x
  /*!
    The result will be different from a sequential execution
    since values can be stored in different processors
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    Matrix<T, Prop, Storage, Allocator>::Fill(x);
    DistributedMatrix_Base<T>::Fill(x);
  }


  //! sets values of non-zero entries to random values
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::FillRand()
  {
    Matrix<T, Prop, Storage, Allocator>::FillRand();
  }
  

  //! writes the matrix in binary format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    const MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::Write(FileName);
    
    cout << "Write not implemented for distributed matrices" << endl;
    abort();
  }
  

  //! writes the matrix in binary format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {
    const MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::Write(FileStream);
    
    cout << "Write not implemented for distributed matrices" << endl;
    abort();
  }
  

  //! writes the matrix in text format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::WriteText(string FileName, bool cplx) const
  {
    const MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::WriteText(FileName, cplx);

    // a different file name for each processor
    string name = GetBaseString(FileName) + "_P" + to_str(rank_proc)
      + "." + GetExtension(FileName);
    
    // opening the stream
    ofstream FileStream(name.c_str());
    FileStream.precision(15);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("DistributedMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + name + "\".");
#endif
    
    // then writing datas
    WriteText(FileStream, cplx);
    
    // closing files
    FileStream.close();
  }
  

  //! writes the matrix in text format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream, bool cplx) const
  {    
    const MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::WriteText(FileStream, cplx);

    // converting local part into coordinate form
    Vector<int> IndRow, IndCol;
    Vector<T> Value;
    ConvertMatrix_to_Coordinates(*this, IndRow, IndCol, Value, 0, true);
    
    DistributedMatrix_Base<T>::WriteText(FileStream,
                                         IndRow, IndCol, Value, cplx);
  }


  //! reads the matrix in binary format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Read(string FileName)
  {
    MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::Read(FileName);

    cout << "Read not implemented for distributed matrices" << endl;
    abort();
  }
  

  //! reads the matrix in binary format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream)
  {
    MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::Read(FileStream);

    cout << "Read not implemented for distributed matrices" << endl;
    abort();
  }
  

  //! reads the matrix in text format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ReadText(string FileName, bool cplx)
  {
    MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::ReadText(FileName, cplx);

    cout << "ReadText not implemented for distributed matrices" << endl;
    abort();
  }
  

  //! reads the matrix in text format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ReadText(istream& FileStream, bool cplx)
  {
    MPI_Comm& comm = this->comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::ReadText(FileStream, cplx);
    
    cout << "ReadText not implemented for distributed matrices" << endl;
    abort();
  }


  //! grouping all the local rows of the matrix into CSR form
  /*!
    Creation of a sparse matrix B, that will contain the local rows
    of the current matrix. The column numbers are global.
    Values placed on distant rows are ignored
  */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Allocator0>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetDistributedRows(Matrix<T0, General, ArrayRowSparse,
		       Allocator0>& B, Vector<IVect>& procB) const
  {
    Copy(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), B);
    DistributedMatrix_Base<T>::GetDistributedRows(B, procB, IsSymmetricMatrix(*this));
  }
  

  //! grouping all the local columns of the matrix into CSC form
  /*!
    Creation of a sparse matrix B, that will contain the local columns
    of the current matrix. The row numbers are global.
    Values placed on distant columns are ignored
  */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Allocator0>
  void DistributedMatrix<T, Prop, Storage, Allocator>::
  GetDistributedColumns(Matrix<T0, General, ArrayColSparse, Allocator0>& B,
                        Vector<IVect>& procB, bool sym_pattern) const
  {
    // conversion to CSC format of local part and symmetrisation of pattern
    Vector<long> Ptr; Vector<int> Ind; Vector<T0> Val;
    General sym;
    ConvertToCSC(*this, sym, Ptr, Ind, Val, sym_pattern);

    DistributedMatrix_Base<T>::
      GetDistributedColumns(B, procB, Ptr, Ind, Val, sym_pattern);
    
  }

  
  // DistributedMatrix //
  ///////////////////////
  
  
  /*************************
   * Matrix vector product *
   *************************/

    
  //! matrix vector product with a distributed matrix
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Storage4, Allocator4>& Yres, bool assemble)
  {    
    bool proceed_distant_row, proceed_distant_col;
    Vector<T4, Storage4, Allocator4> Y; Vector<T2> Xcol;
    M.InitMltAdd(proceed_distant_row, proceed_distant_col,
                 X, Xcol, beta, Y, Yres);
    
    // local matrix
    MltVector(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(M),
	      X, Y);
    
    // distributed contribution
    M.FinalizeMltAdd(proceed_distant_row, proceed_distant_col,
                     X, Xcol, alpha, beta, Y, Yres, assemble);    
  }
  

  //! matrix vector product with a distributed matrix
  template<class T0, class T1, class Prop1, class Storage1,
	   class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Storage4, Allocator4>& Yres, bool assemble)
  {
    if (Trans.NoTrans())
      {
	MltAddVector(alpha, M, X, beta, Yres, assemble);
	return;
      }

    bool proceed_distant_row, proceed_distant_col;    
    Vector<T4, Storage4, Allocator4> Y; Vector<T2> Xrow;
    M.InitMltAdd(proceed_distant_row, proceed_distant_col,
                 Trans, X, Xrow, beta, Y, Yres);
    
    // local matrix
    MltVector(Trans, static_cast<const Matrix<T1, Prop1,
	      Storage1, Allocator1>& >(M), X, Y);
    
    M.FinalizeMltAdd(proceed_distant_row, proceed_distant_col,
                     Trans, X, Xrow, alpha, beta, Y, Yres, assemble);
  }
  

  //! matrix-vector product selecting the minimal processor number
  //! and column number for each row
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          int col = A.Index(i, j);
          if (Yproc(col) < Yproc(i))
            {
              Yproc(i) = Yproc(col);
              Y(i) = Y(col);
            }
          else if (Yproc(col) == Yproc(i))
            {
              if (Y(col) < Y(i))
                Y(i) = Y(col);
              else
                Y(col) = Y(i);
            }
          else
            {
              Yproc(col) = Yproc(i);
              Y(col) = Y(i);
            }
        }
  }


  //! matrix-vector product selecting the minimal processor number
  //! and column number for each row
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          int col = A.Index(i, j);
          if (Yproc(col) < Yproc(i))
            {
              Yproc(i) = Yproc(col);
              Y(i) = Y(col);
            }
          else if (Yproc(col) == Yproc(i))
            {
              if (Y(col) < Y(i))
                Y(i) = Y(col);
              else
                Y(col) = Y(i);
            }
          else
            {
              Yproc(col) = Yproc(i);
              Y(col) = Y(i);
            }
        }
  }


#ifdef SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX
  //! matrix-vector product selecting the minimal processor number
  //! and column number for each row
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int col = A.IndexReal(i, j);
            if (Yproc(col) < Yproc(i))
              {
                Yproc(i) = Yproc(col);
                Y(i) = Y(col);
              }
            else if (Yproc(col) == Yproc(i))
              {
                if (Y(col) < Y(i))
                  Y(i) = Y(col);
                else
                  Y(col) = Y(i);
              }
            else
              {
                Yproc(col) = Yproc(i);
                Y(col) = Y(i);
              }
          }

        for (int j = 0; j < A.GetImagRowSize(i); j++)
          {
            int col = A.IndexImag(i, j);
            if (Yproc(col) < Yproc(i))
              {
                Yproc(i) = Yproc(col);
                Y(i) = Y(col);
              }
            else if (Yproc(col) == Yproc(i))
              {
                if (Y(col) < Y(i))
                  Y(i) = Y(col);
                else
                  Y(col) = Y(i);
              }
            else
              {
                Yproc(col) = Yproc(i);
                Y(col) = Y(i);
              }
          }
      }
  }


  //! matrix-vector product selecting the minimal processor number
  //! and column number for each row
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1,
	      ArrayRowSymComplexSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int col = A.IndexReal(i, j);
            if (Yproc(col) < Yproc(i))
              {
                Yproc(i) = Yproc(col);
                Y(i) = Y(col);
              }
            else if (Yproc(col) == Yproc(i))
              {
                if (Y(col) < Y(i))
                  Y(i) = Y(col);
                else
                  Y(col) = Y(i);
              }
            else
              {
                Yproc(col) = Yproc(i);
                Y(col) = Y(i);
              }
          }

        for (int j = 0; j < A.GetImagRowSize(i); j++)
          {
            int col = A.IndexImag(i, j);
            if (Yproc(col) < Yproc(i))
              {
                Yproc(i) = Yproc(col);
                Y(i) = Y(col);
              }
            else if (Yproc(col) == Yproc(i))
              {
                if (Y(col) < Y(i))
                  Y(i) = Y(col);
                else
                  Y(col) = Y(i);
              }
            else
              {
                Yproc(col) = Yproc(i);
                Y(col) = Y(i);
              }
          }
      }
  }
#endif
  
              
  //! matrix vector product with a distributed matrix (minimal column numbers)
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void MltMin(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              IVect& Y, IVect& Yproc)
  {
    IVect Xcol, Xcol_proc;
    M.InitMltMin(Y, Yproc, Xcol, Xcol_proc);

    // local matrix
    const IVect& global = M.GetGlobalRowNumber();
    MltMin(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(M),
           global, Y, Yproc);
    
    M.FinalizeMltMin(Y, Yproc, Xcol, Xcol_proc);
  }

  
  /**************************
   * Functions for matrices *
   **************************/

  
  //! Adds two distributed matrices (B = B + alpha A)
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A, 
		 DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    // adding local part
    AddMatrix(alpha,
	      static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(A),
	      static_cast<Matrix<T2, Prop2, Storage2, Allocator2>& >(B));
    
    B.AddDistributedMatrix(alpha, A);
  }
  
  
  //! returns the maximal element (in absolute value) of the matrix
  /*!
    The result will be different from a sequential execution
    since A(i, j) is obtained as a sum of values stored
    in different processors
   */
  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  MaxAbs(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    typename ClassComplexType<T1>::Treal res;
    res = MaxAbs(static_cast<const Matrix<T1, Prop1,
		 Storage1, Allocator1>& >(A));
    
    A.GetMaxAbsDistant(res);
    return res;
  }

  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetRowSum(Vector<T0>& vec_sum,
                 const DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    GetRowSum(vec_sum,
	      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(A));
    
    A.AddRowSumDistant(vec_sum);
  }
  

  //! For each column of the matrix, computation of 
  //! the sum of absolute values
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetColSum(Vector<T0>& vec_sum,
                 const DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    GetColSum(vec_sum,
	      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(A));
    
    A.AddColSumDistant(vec_sum);
  }

  
  //! gets sum of absolute value for each row and each colum
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetRowColSum(Vector<T0>& sum_row, Vector<T0>& sum_col,
                    const DistributedMatrix<T, Prop, Storage, Allocator> & A)
  {    
    GetRowColSum(sum_row, sum_col,
                 static_cast<const Matrix<T, Prop, Storage, Allocator>& >(A));
    
    A.AddRowColSumDistant(sum_row, sum_col);
  }
  
  
  //! returns the 1-norm of A
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  Norm1(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Norm1(static_cast<const Matrix<T1, Prop1,
		   Storage1, Allocator1>& >(A));
    
    Vector<typename ClassComplexType<T1>::Treal> sum_col;
    GetColSum(sum_col, A);

    typename ClassComplexType<T1>::Treal res, amax;
    amax = 0;
    for (int i = 0; i < sum_col.GetM(); i++)
      amax = max(amax, abs(sum_col(i)));
    
    Vector<int64_t> xtmp;
    MpiAllreduce(comm, &amax, xtmp, &res, 1, MPI_MAX);
    return res;
  }


  //! returns the infinity norm of A
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  NormInf(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return NormInf(static_cast<const Matrix<T1, Prop1,
		     Storage1, Allocator1>& >(A));
    
    Vector<typename ClassComplexType<T1>::Treal> sum_row;
    GetRowSum(sum_row, A);

    typename ClassComplexType<T1>::Treal res, amax;
    amax = 0;
    for (int i = 0; i < sum_row.GetM(); i++)
      amax = max(amax, abs(sum_row(i)));
    
    Vector<int64_t> xtmp;
    MpiAllreduce(comm, &amax, xtmp, &res, 1, MPI_MAX);
    return res;
  }


  //! replaces A by its transpose
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Transpose(DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    
    // storing the distributed datas
    Vector<Vector<T1, VectSparse>, VectFull,
           NewAlloc<Vector<T1, VectSparse> > > dist_row, dist_col;
    
    bool local_number_distant_values;
    int smax_row, smax_col;
    Vector<IVect> proc_row, proc_col;
    
    IVect global_row_to_recv, global_col_to_recv;
    IVect ptr_global_row_to_recv, ptr_global_col_to_recv;
    Vector<IVect> local_row_to_send, local_col_to_send;
    IVect proc_col_to_recv, proc_col_to_send,
      proc_row_to_recv, proc_row_to_send;

    // distributed datas are retrieved without copy
    A.ExchangeParallelData(smax_row, smax_col, local_number_distant_values,
                           dist_row, dist_col, proc_row, proc_col,
                           global_row_to_recv, global_col_to_recv,
                           ptr_global_row_to_recv, ptr_global_col_to_recv,
                           local_row_to_send, local_col_to_send,
                           proc_row_to_recv, proc_col_to_recv,
                           proc_row_to_send, proc_col_to_send);

    // local matrix is transposed (A may be erased during the process)
    Transpose(static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(A));

    // transposing distributed datas    
    A.ExchangeParallelData(smax_col, smax_row, local_number_distant_values,
                           dist_col, dist_row, proc_col, proc_row,
                           global_col_to_recv, global_row_to_recv,
                           ptr_global_col_to_recv, ptr_global_row_to_recv,
                           local_col_to_send, local_row_to_send,
                           proc_col_to_recv, proc_row_to_recv,
                           proc_col_to_send, proc_row_to_send);
  }
  

  //! replaces A by its conjugate
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Conjugate(DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    Conjugate(static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(A));
    
    A.ConjugateDistant();
  }
  
  
  //! computes B = transpose(A)
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Transpose(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
                 DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    Transpose(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(A),
              static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(B));
    
    B.TransposeDistant(A);
  }
  
  
  //! Matrix transposition and conjugation.
  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                     DistributedMatrix<T, Prop, Storage, Allocator>& B)
  {
    Transpose(A, B);
    Conjugate(B);
  }


  //! Matrix transposition and conjugation.
  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    Transpose(A);
    Conjugate(A);
  }

  
  /**************************
   * Matrix-matrix products *
   **************************/
  

  //! computes C = A B
  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltMatrix(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
		 const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
		 DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return
	MltMatrix(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(A),
		  static_cast<const Matrix<T2, Prop2, Storage2, Allocator2>& >(B),
		  static_cast<Matrix<T4, Prop4, Storage4, Allocator4>& >(C));
    
    cout << "Mlt not implemented for distributed matrices" << endl;
    abort();
  }

  
  //! computes C = beta C + alpha A B
  template<class T0,
           class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T3,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltAddMatrix(const T0& alpha,
		    const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
		    const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
		    const T3& beta,
		    DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return
	MltAddMatrix(alpha,
		     static_cast<const Matrix<T1, Prop1,Storage1,Allocator1>& >(A),
		     static_cast<const Matrix<T2, Prop2,Storage2,Allocator2>& >(B),
		     beta,
		     static_cast<Matrix<T4, Prop4, Storage4, Allocator4>& >(C));
    
    cout << "MltAdd not implemented for distributed matrices" << endl;
    abort();
  }

  
  //! computes C = beta C + alpha A B (with transposes)
  template<class T0,
           class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T3,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltAddMatrix(const T0& alpha, const SeldonTranspose& transA,
		    const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
		    const SeldonTranspose& transB,
		    const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
		    const T3& beta,
		    DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return MltAddMatrix(alpha, transA,
			  static_cast<const Matrix<T1, Prop1,
			  Storage1, Allocator1>& >(A),
			  transB,
			  static_cast<const Matrix<T2, Prop2,
			  Storage2, Allocator2>& >(B),
			  beta,
			  static_cast<Matrix<T4, Prop4,
			  Storage4, Allocator4>& >(C));
    
    cout << "MltAdd not implemented for distributed matrices" << endl;
    abort();
  }
  

  /********************
   * Matrix functions *
   ********************/

  
  //! X = M(i, :)
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void GetRow(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
              int i, Vector<T1, VectSparse, Allocator1>& X)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return GetRow(static_cast<const Matrix<T0, Prop0,
		    Storage0, Allocator0>& >(A), i, X);
    
    cout << "GetRow not implemented for distributed matrices" << endl;
    abort();
  }


  //! X = M(:, i)
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void GetCol(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
              int i, Vector<T1, VectSparse, Allocator1>& X)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return GetCol(static_cast<const Matrix<T0, Prop0,
		    Storage0, Allocator0>& >(A), i, X);
    
    cout << "GetCol not implemented for distributed matrices" << endl;
    abort();
  }
  
  
  //! M(i, :) = X
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
              int i, DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return SetRow(X, i, static_cast<Matrix<T0, Prop0,
		    Storage0, Allocator0>& >(A));
    
    cout << "SetRow not implemented for distributed matrices" << endl;
    abort();
  }
  

  //! M(:, i) = X
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
              int i, DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return SetCol(X, i, static_cast<Matrix<T0, Prop0,
		    Storage0, Allocator0>& >(A));
    
    cout << "SetCol not implemented for distributed matrices" << endl;
    abort();
  }

  
  //! applies permutation on rows and columns of A
  template<class T, class Prop, class Storage, class Allocator>
  void ApplyPermutation(DistributedMatrix<T, Prop, Storage, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return ApplyPermutation(static_cast<Matrix<T, Prop,
			      Storage, Allocator>& >(A),
                              row_perm, col_perm);
    
    cout << "ApplyPermutation not implemented for distributed matrices"
	 << endl;
    abort();
  }
  

  //! applies permutation on rows and columns of A
  template<class T, class Prop, class Storage, class Allocator>
  void ApplyInversePermutation(DistributedMatrix<T, Prop,
			       Storage, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return ApplyInversePermutation(static_cast<Matrix<T, Prop, Storage,
				     Allocator>& >(A),
                                     row_perm, col_perm);
    
    cout << "ApplyInversePermutation not implemented for distributed matrices"
	 << endl;
    abort();
  }

  
  //! applies successive over-relaxation steps to solve A X = B
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return SorVector(static_cast<const Matrix<T0, Prop0,
		       Storage0, Allocator0>& >(A),
		       X, B, omega, iter, type_ssor);
    
    cout << "SOR not implemented for distributed matrices" << endl;
    abort();    
  }
  
  
  //! applies successive over-relaxation steps to solve A^T X = B
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);
    
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return SOR(transM,
                 static_cast<const Matrix<T0, Prop0,
		 Storage0, Allocator0>& >(A),
                 X, B, omega, iter, type_ssor);
    
    cout << "SOR not implemented for distributed matrices" << endl;
    abort();
  }
  
  
  //! extracts a set of sparse columns from matrix A
  template<class T1, class Prop, class Storage, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const DistributedMatrix<T1, Prop, Storage, Allocator>& A,
              const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>,
	      VectSparse, Allocator3>& V)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return GetCol(static_cast<const Matrix<T1, Prop,
		    Storage, Allocator>& >(A),
                    col_number, V);
    
    cout << "GetCol not implemented for distributed matrices" << endl;
    abort();
  }
  
  
  //! converts a distributed matrix to another one
  template<class T, class Prop1, class Storage1, class Allocator1,
           class Prop2, class Storage2, class Allocator2>
  void Copy(const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
            DistributedMatrix<T, Prop2, Storage2, Allocator2>& B)
  {
    B = A;
  }


  //! converts a distributed matrix to another one
  template<class T, class Prop1, class Storage1, class Allocator1,
           class Prop2, class Storage2, class Allocator2>
  void Copy(const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
            DistributedMatrix<complex<T>, Prop2, Storage2, Allocator2>& B)
  {
    B = A;
  }

  
  //! extracts real part of a matrix
  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void CopyReal(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
		DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    MPI_Comm& comm = B.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return CopyReal(static_cast<const Matrix<T1, Prop1,
		      Storage1, Allocator1>& >(A),
                      static_cast<Matrix<T2, Prop2,
		      Storage2, Allocator2>& >(B));

    cout << "CopyReal not implemented for distributed matrices" << endl;
    abort();
  }
  
  
  //! returns Froebenius norm of matrix A
  template<class T, class Prop, class Storage, class Allocator>
  typename ClassComplexType<T>::Treal
  NormFro(const DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return NormFro(static_cast<const Matrix<T, Prop,
		     Storage, Allocator>& >(A));

    cout << "NormFro not implemented for distributed matrices" << endl;
    abort();    
  }
  

  //! multiplies rows by coefficients
  template<class T, class Storage, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(DistributedMatrix<T, General, Storage, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    ScaleLeftMatrix(static_cast<Matrix<T, General,
		    Storage, Allocator>& >(A), Drow);
    
    A.ScaleLeftDistant(Drow);
  }


  //! multiplies columns by coefficients
  template<class T, class Storage, class Allocator,
           class T1, class Allocator1>
  void ScaleRightMatrix(DistributedMatrix<T, General, Storage, Allocator>& A,
                        const Vector<T1, VectFull, Allocator1>& Dcol)
  {
    ScaleRightMatrix(static_cast<Matrix<T, General,
		     Storage, Allocator>& >(A), Dcol);
    
    A.ScaleRightDistant(Dcol);
  }


  //! multiplies rows and columns by coefficients
  template<class T, class Prop, class Storage, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(DistributedMatrix<T, Prop, Storage, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    ScaleMatrix(static_cast<Matrix<T, Prop, Storage, Allocator>& >(A),
		Drow, Dcol);
    
    A.ScaleDistant(Drow, Dcol);
  }


  //! assembling a distributed matrix, global rows are given in CSR form
  /*!
    \param[in] A distributed matrix
    \param[in] sym the matrix is assumed to be symmetric,
    only upper part is assembled
    \param[in] comm MPI communicator
    \param[out] row_numbers global row numbers
    \param[out] local_row_numbers local row numbers
    \param[out] PtrA Index of first element stored for each row
    \param[out] IndA global column numbers
    \param[out] ValA values
    \param[in] sym_pattern not used
   */
  template<class Prop, class Storage, class Alloc, class Tint0, class Tint1, class T>
  void AssembleDistributed(DistributedMatrix<T, Prop, Storage, Alloc>& A,
			   Symmetric& sym, const MPI_Comm& comm,
                           IVect& row_numbers, IVect& local_row_numbers,
			   Vector<Tint0>& PtrA, Vector<Tint1>& IndA,
			   Vector<T>& ValA, bool sym_pattern, bool reorder)
  {
    PtrA.Clear(); IndA.Clear(); ValA.Clear();

    // we convert A in ArrayRowSparse format
    Matrix<T, General, ArrayRowSparse> B;
    Vector<IVect> procB;
    if (reorder)
      procB.Reallocate(A.GetM());
    
    A.GetDistributedRows(B, procB);
    
    // local matrix is cleared
    A.ClearLocal();
    
    // then calling AssembleParallel
    IVect OverlappedCol;
    A.AssembleParallel(B, procB, sym, row_numbers, local_row_numbers,
                       OverlappedCol, sym_pattern, reorder);
  
    // final conversion with type Tint
    A.ConvertToCSR(B, OverlappedCol, PtrA, IndA, ValA);
  }
  
  
  //! assembling a distributed matrix, global columns are given in CSC form
  /*!
    \param[in] A distributed matrix
    \param[in] prop the matrix is assumed to be unsymmetric
    \param[in] comm MPI communicator
    \param[out] col_numbers global column numbers
    \param[out] local_col_numbers local column numbers
    \param[out] PtrA Index of first element stored for each column
    \param[out] IndA global row numbers
    \param[out] ValA values
    \param[in] sym_pattern if true,
    the pattern of the unsymmetric matrix is symmetrized
   */
  template<class Prop, class Storage, class Alloc, class Tint0, class Tint1, class T>
  void AssembleDistributed(DistributedMatrix<T, Prop, Storage, Alloc>& A,
			   General& prop, const MPI_Comm& comm,
                           IVect& col_numbers, IVect& local_col_numbers,
                           Vector<Tint0>& PtrA, Vector<Tint1>& IndA,
			   Vector<T>& ValA, bool sym_pattern, bool reorder)
  {
    PtrA.Clear(); IndA.Clear(); ValA.Clear();
    
    // we convert A in CSC format
    Matrix<T, General, ArrayColSparse> B;
    Vector<IVect> procB;
    if (reorder)
      procB.Reallocate(A.GetN());

    A.GetDistributedColumns(B, procB, sym_pattern);

    // local matrix is cleared
    A.ClearLocal();
        
    // then AssembleParallel is called
    IVect OverlappedCol;
    A.AssembleParallel(B, procB, prop, col_numbers, local_col_numbers,
                       OverlappedCol, sym_pattern, reorder);

    // final conversion with type Tint
    A.ConvertToCSC(B, OverlappedCol, PtrA, IndA, ValA);
  }

  
  //! clears several columns of a sparse matrix
  /*!
    \param[in] num numbers of the columns to be cleared
    \param[in] comm MPI communicator
    \param[inout] A sparse matrix where columns are erased
    In this functions the column numbers are assumed to be increasing
   */  
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void EraseCol(const IVect& num,
		DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    // erasing columns of the local entries
    EraseCol(num, static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(A));
    
    A.EraseColDistant(num, IsSymmetricMatrix(A));
  }
  

  //! clears several rows of a sparse matrix
  /*!
    \param[in] num numbers of the rows to be cleared
    \param[in] comm MPI communicator
    \param[inout] A sparse matrix where rows are erased
    In this functions the row numbers are assumed to be increasing
   */  
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void EraseRow(const IVect& num,
		DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    // erasing rows of the local entries
    EraseRow(num, static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(A));
    
    A.EraseRowDistant(num, IsSymmetricMatrix(A));
  }


  //! extracts a sub-matrix B from matrix A
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void
  CopySubMatrix(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		const IVect& row, const IVect& col,
		DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    CopySubMatrix(static_cast<const Matrix<T0, Prop0,
		  Storage0, Allocator0>& >(A),
                  row, col, static_cast<Matrix<T1, Prop1,
		  Storage1, Allocator1>& >(B));
    
    B.CopySubDistant(A, row, col, IsSymmetricMatrix(B));
  }
  
}

#define SELDON_FILE_DISTRIBUTED_MATRIX_CXX
#endif

