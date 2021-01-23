#ifndef SELDON_FILE_MPI_COMMUNICATION_HXX

namespace Seldon
{
  
  MPI_Datatype GetMpiDataType(const Vector<bool>&);
  MPI_Datatype GetMpiDataType(const Vector<int>&);
  MPI_Datatype GetMpiDataType(const Vector<long>&);
  MPI_Datatype GetMpiDataType(const Vector<float>&);
  MPI_Datatype GetMpiDataType(const Vector<complex<float> >&);
  MPI_Datatype GetMpiDataType(const Vector<double>&);
  MPI_Datatype GetMpiDataType(const Vector<complex<double> >&);
  MPI_Datatype GetMpiDataType(const Vector<long double>&);
  MPI_Datatype GetMpiDataType(const Vector<complex<long double> >&);

  MPI_Datatype GetMpiDataType(const bool&);
  MPI_Datatype GetMpiDataType(const int&);
  MPI_Datatype GetMpiDataType(const long&);
  MPI_Datatype GetMpiDataType(const float&);
  MPI_Datatype GetMpiDataType(const complex<float>&);
  MPI_Datatype GetMpiDataType(const double&);
  MPI_Datatype GetMpiDataType(const complex<double>&);
  MPI_Datatype GetMpiDataType(const long double&);
  MPI_Datatype GetMpiDataType(const complex<long double>&);
  
  template<class T>
  int GetRatioMpiDataType(const T&);
 
  template<class T>
  int GetRatioMpiDataType(const Vector<T>&);
 
  template<class T>
  int GetRatioMpiDataType(const complex<T>&);

  template<class T>
  int GetRatioMpiDataType(const Vector<complex<T> >&);

  template<class T>
  MPI_Request MpiIsend(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);
  
  template<class T>
  MPI_Request MpiIsend(const MPI_Comm& comm, Vector<T>& x,
			Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  template<class T>
  MPI_Request MpiIrecv(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  template<class T>
  MPI_Request MpiIrecv(const MPI_Comm& comm, Vector<T>& x,
			Vector<int64_t>& xtmp,
                        int n, int proc, int tag);
  
  template<class T>
  void MpiCompleteIrecv(T* x, Vector<int64_t>& xtmp, int n);
  
  template<class T>
  void MpiCompleteIrecv(Vector<T>& x, Vector<int64_t>& xtmp, int n);

  template<class T>
  void MpiSsend(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag);
  
  template<class T>
  void MpiSsend(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                int n, int proc, int tag);
                       
  template<class T>
  void MpiSend(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag);
  
  template<class T>
  void MpiSend(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag);
  
  template<class T>
  void MpiGather(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                 T* y, int n, int proc);
  
  template<class T>
  void MpiGather(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                 Vector<T>& y, int n, int proc);
  
  template<class T>
  void MpiAllreduce(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                    T* y, int n, const MPI_Op& op);
  
  template<class T>
  void MpiAllreduce(const MPI_Comm& comm, Vector<T>& x,
		    Vector<int64_t>& xtmp,
                    Vector<T>& y, int n, const MPI_Op& op);
  
  template<class T>
  void MpiReduce(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
                 T* y, int n, const MPI_Op& op, int proc);
  
  template<class T>
  void MpiReduce(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                 Vector<T>& y, int n, const MPI_Op& op, int proc);
  
  template<class T>
  void MpiRecv(const MPI_Comm& comm, T* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status);
  
  template<class T>
  void MpiRecv(const MPI_Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status);

  template<class T>
  void MpiBcast(const MPI_Comm& comm, T* x,
		Vector<int64_t>& xtmp, int n, int proc);
  
  template<class T>
  void MpiBcast(const MPI_Comm& comm, Vector<T>& x,
		Vector<int64_t>& xtmp, int n, int proc);
  
}

#define SELDON_FILE_MPI_COMMUNICATION_HXX
#endif
