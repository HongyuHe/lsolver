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

#ifndef SELDON_FILE_MPI_COMMUNICATION_CXX

namespace Seldon
{

  template<class T>
  MPI_Request MpiIsend(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    MPI_Request request;
    MPI_Isend(x, n*GetRatioMpiDataType(*x), 
              GetMpiDataType(*x), proc, tag, comm, &request);
    return request;
  }
  
  template<class T>
  MPI_Request MpiIsend(const MPI_Comm& comm, Vector<T>& x,
			Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    MPI_Request request;
    MPI_Isend(x.GetData(), n*GetRatioMpiDataType(x),
              GetMpiDataType(x), proc, tag, comm, &request); 
    return request;
  }

  template<class T>
  MPI_Request MpiIrecv(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    MPI_Request request;
    MPI_Irecv(x, n*GetRatioMpiDataType(*x), 
              GetMpiDataType(*x), proc, tag, comm, &request);
    return request;
  }

  template<class T>
  MPI_Request MpiIrecv(const MPI_Comm& comm, Vector<T>& x,
			Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    MPI_Request request;
    MPI_Irecv(x.GetData(), n*GetRatioMpiDataType(x),
              GetMpiDataType(x), proc, tag, comm, &request); 
    return request;
  }
  
  template<class T>
  void MpiCompleteIrecv(T* x, Vector<int64_t>& xtmp, int n)
  {
  }
  
  template<class T>
  void MpiCompleteIrecv(Vector<T>& x, Vector<int64_t>& xtmp, int n)
  {    
  }

  template<class T>
  void MpiSsend(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    MPI_Send(x, n*GetRatioMpiDataType(*x),
             GetMpiDataType(*x), proc, tag, comm);
  }
  
  template<class T>
  void MpiSsend(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    MPI_Send(x.GetData(), n*GetRatioMpiDataType(x),
             GetMpiDataType(x), proc, tag, comm);
  }
                       
  template<class T>
  void MpiSend(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    MPI_Send(x, n*GetRatioMpiDataType(*x),
             GetMpiDataType(*x), proc, tag, comm);
  }
  
  template<class T>
  void MpiSend(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag)
  {
    MPI_Send(x.GetData(), n*GetRatioMpiDataType(x),
             GetMpiDataType(x), proc, tag, comm);
  }
                     
  template<class T>
  void MpiGather(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                 T* y, int n, int proc)
  {
    MPI_Gather(x, n*GetRatioMpiDataType(*x), GetMpiDataType(*x),
               y, n*GetRatioMpiDataType(*x), GetMpiDataType(*x),
               proc, comm);
  }
  
  template<class T>
  void MpiGather(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                 Vector<T>& y, int n, int proc)
  {
    MPI_Gather(x.GetData(), n*GetRatioMpiDataType(x), GetMpiDataType(x),
               y.GetData(), n*GetRatioMpiDataType(x), GetMpiDataType(x),
               proc, comm);
  }
  
  template<class T>
  void MpiAllreduce(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                    T* y, int n, const MPI_Op& op)
  {
    MPI_Allreduce(x, y, n*GetRatioMpiDataType(*x), GetMpiDataType(*x), op, comm);
  }
  
  template<class T>
  void MpiAllreduce(const MPI_Comm& comm, Vector<T>& x,
		    Vector<int64_t>& xtmp,
                    Vector<T>& y, int n, const MPI_Op& op)
  {
    MPI_Allreduce(x.GetData(), y.GetData(), n*GetRatioMpiDataType(x),
                  GetMpiDataType(x), op, comm);
  }
  
  template<class T>
  void MpiReduce(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                 T* y, int n, const MPI_Op& op, int proc)
  {
    MPI_Reduce(x, y, n*GetRatioMpiDataType(*x),
               GetMpiDataType(*x), op, proc, comm);
  }
  
  template<class T>
  void MpiReduce(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                 Vector<T>& y, int n, const MPI_Op& op, int proc)
  {
    MPI_Reduce(x.GetData(), y.GetData(), n*GetRatioMpiDataType(x),
               GetMpiDataType(x), op, proc, comm);
  }
  
  template<class T>
  void MpiRecv(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status)
  {
    MPI_Recv(x, n*GetRatioMpiDataType(*x), GetMpiDataType(*x),
             proc, tag, comm, &status);
  }
  
  template<class T>
  void MpiRecv(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status)
  {
    MPI_Recv(x.GetData(), n*GetRatioMpiDataType(x), GetMpiDataType(x),
             proc, tag, comm, &status);
  }

  template<class T>
  void MpiBcast(const MPI_Comm& comm, T* x,
		Vector<int64_t>& xtmp, int n, int proc)
  {
    MPI_Bcast(x, n*GetRatioMpiDataType(*x), GetMpiDataType(*x), proc, comm);
  }
  
  template<class T>
  void MpiBcast(const MPI_Comm& comm, Vector<T>& x,
		Vector<int64_t>& xtmp, int n, int proc)
  {
    MPI_Bcast(x.GetData(), n*GetRatioMpiDataType(x), GetMpiDataType(x), proc, comm);
  }
  
}

#define SELDON_FILE_MPI_COMMUNICATION_CXX
#endif
