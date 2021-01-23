#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifdef SELDON_WITH_MPI
#include "SeldonDistributedHeader.hxx"
#include "SeldonDistributedInline.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "vector/DistributedVector.cxx"
#include "share/MpiCommunication.cxx"
#endif

namespace Seldon
{
  
  // dense vectors
  SELDON_EXTERN template class DistributedVector<Real_wp>;
  SELDON_EXTERN template class DistributedVector<Complex_wp>;

  // fonctions DotProd, Norm2
  SELDON_EXTERN template Real_wp DotProdVector(const DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&);
  SELDON_EXTERN template Complex_wp DotProdVector(const DistributedVector<Complex_wp>&, const DistributedVector<Complex_wp>&);
  SELDON_EXTERN template Real_wp DotProdConjVector(const DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&);
  SELDON_EXTERN template Complex_wp DotProdConjVector(const DistributedVector<Complex_wp>&, const DistributedVector<Complex_wp>&);
  SELDON_EXTERN template Real_wp Norm2(const DistributedVector<Real_wp>&);
  SELDON_EXTERN template Real_wp Norm2(const DistributedVector<Complex_wp>&);
  
  // fonctions AssembleVector
  SELDON_EXTERN template void AssembleVector(Vector<bool>&, const MPI_Op&, const IVect&, const Vector<IVect>&, const MPI_Comm&, int, int, int);
  SELDON_EXTERN template void AssembleVector(Vector<int>&, const MPI_Op&, const IVect&, const Vector<IVect>&, const MPI_Comm&, int, int, int);
  SELDON_EXTERN template void AssembleVector(Vector<Real_wp>&, const MPI_Op&, const IVect&, const Vector<IVect>&, const MPI_Comm&, int, int, int);
  SELDON_EXTERN template void AssembleVector(Vector<Complex_wp>&, const MPI_Op&, const IVect&, const Vector<IVect>&, const MPI_Comm&, int, int, int);
  
  SELDON_EXTERN template void ExchangeVector(Vector<Real_wp>&, const IVect&, const Vector<IVect>&, const MPI_Comm&, int, int, int);
  SELDON_EXTERN template void ExchangeVector(Vector<Complex_wp>&, const IVect&, const Vector<IVect>&, const MPI_Comm&, int, int, int);
  SELDON_EXTERN template void ExchangeRelaxVector(Vector<Real_wp>&, const Real_wp&, int, const IVect&, const Vector<IVect>&, const MPI_Comm&, int, int, int);
  SELDON_EXTERN template void ExchangeRelaxVector(Vector<Complex_wp>&, const Real_wp&, int, const IVect&, const Vector<IVect>&, const MPI_Comm&, int, int, int);
  
  // MPI communications
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, bool*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, int*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, long*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, double*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, complex<double>*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Vector<bool>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Vector<int>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Vector<long>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, int, int, int);

  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, bool*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, int*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, long*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, double*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, complex<double>*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Vector<bool>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Vector<int>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Vector<long>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, int, int, int);
  
  SELDON_EXTERN template void MpiCompleteIrecv(bool*, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(int*, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(long*, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(double*, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(complex<double>*, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(Vector<bool>&, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(Vector<int>&, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(Vector<long>&, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(Vector<double>&, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(Vector<complex<double> >&, Vector<int64_t>&, int);
  
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, int*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, long*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, double*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, complex<double>*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, Vector<int>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, Vector<long>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, int, int, int);
  
  SELDON_EXTERN template void MpiSend(const MPI_Comm&, double*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSend(const MPI_Comm&, complex<double>*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSend(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSend(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, int, int, int);
  
  SELDON_EXTERN template void MpiGather(const MPI_Comm&, double*, Vector<int64_t>&, double*, int, int);
  SELDON_EXTERN template void MpiGather(const MPI_Comm&, complex<double>*, Vector<int64_t>&, complex<double>*, int, int);
  SELDON_EXTERN template void MpiGather(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, Vector<double>&, int, int);
  SELDON_EXTERN template void MpiGather(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, Vector<complex<double> >&, int, int);

  SELDON_EXTERN template void MpiAllreduce(const MPI_Comm&, double*, Vector<int64_t>&, double*, int, const MPI_Op&);
  SELDON_EXTERN template void MpiAllreduce(const MPI_Comm&, complex<double>*, Vector<int64_t>&, complex<double>*, int, const MPI_Op&);
  SELDON_EXTERN template void MpiAllreduce(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, Vector<double>&, int, const MPI_Op&);
  SELDON_EXTERN template void MpiAllreduce(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, Vector<complex<double> >&, int, const MPI_Op&);
  
  SELDON_EXTERN template void MpiReduce(const MPI_Comm&, double*, Vector<int64_t>&, double*, int, const MPI_Op&, int);
  SELDON_EXTERN template void MpiReduce(const MPI_Comm&, complex<double>*, Vector<int64_t>&, complex<double>*, int, const MPI_Op&, int);
  SELDON_EXTERN template void MpiReduce(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, Vector<double>&, int, const MPI_Op&, int);
  SELDON_EXTERN template void MpiReduce(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, Vector<complex<double> >&, int, const MPI_Op&, int);
  
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, int*, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, long*, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, double*, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, complex<double>*, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, Vector<int>&, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, Vector<long>&, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, int, int, int, MPI_Status&);
  
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, double*, Vector<int64_t>&, int, int);
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, complex<double>*, Vector<int64_t>&, int, int);
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, Vector<bool>&, Vector<int64_t>&, int, int);
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, Vector<double>&, Vector<int64_t>&, int, int);
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, Vector<complex<double> >&, Vector<int64_t>&, int, int);
  
#ifdef SELDON_WITH_MULTIPLE
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Real_wp*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIsend(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, int, int, int);

  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Real_wp*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template MPI_Request MpiIrecv(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, int, int, int);

  SELDON_EXTERN template void MpiCompleteIrecv(Real_wp*, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(Complex_wp*, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(Vector<Real_wp>&, Vector<int64_t>&, int);
  SELDON_EXTERN template void MpiCompleteIrecv(Vector<Complex_wp >&, Vector<int64_t>&, int);

  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, Real_wp*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSsend(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, int, int, int);

  SELDON_EXTERN template void MpiSend(const MPI_Comm&, Real_wp*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSend(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSend(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, int, int, int);
  SELDON_EXTERN template void MpiSend(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, int, int, int);

  SELDON_EXTERN template void MpiGather(const MPI_Comm&, Real_wp*, Vector<int64_t>&, Real_wp*, int, int);
  SELDON_EXTERN template void MpiGather(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, Complex_wp*, int, int);
  SELDON_EXTERN template void MpiGather(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, Vector<Real_wp>&, int, int);
  SELDON_EXTERN template void MpiGather(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, Vector<Complex_wp >&, int, int);

  SELDON_EXTERN template void MpiAllreduce(const MPI_Comm&, Real_wp*, Vector<int64_t>&, Real_wp*, int, const MPI_Op&);
  SELDON_EXTERN template void MpiAllreduce(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, Complex_wp*, int, const MPI_Op&);
  SELDON_EXTERN template void MpiAllreduce(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, Vector<Real_wp>&, int, const MPI_Op&);
  SELDON_EXTERN template void MpiAllreduce(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, Vector<Complex_wp >&, int, const MPI_Op&);
  
  SELDON_EXTERN template void MpiReduce(const MPI_Comm&, Real_wp*, Vector<int64_t>&, Real_wp*, int, const MPI_Op&, int);
  SELDON_EXTERN template void MpiReduce(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, Complex_wp*, int, const MPI_Op&, int);
  SELDON_EXTERN template void MpiReduce(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, Vector<Real_wp>&, int, const MPI_Op&, int);
  SELDON_EXTERN template void MpiReduce(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, Vector<Complex_wp >&, int, const MPI_Op&, int);
  
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, Real_wp*, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, int, int, int, MPI_Status&);
  SELDON_EXTERN template void MpiRecv(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, int, int, int, MPI_Status&);
  
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, Real_wp*, Vector<int64_t>&, int, int);
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, Complex_wp*, Vector<int64_t>&, int, int);
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, Vector<Real_wp>&, Vector<int64_t>&, int, int);
  SELDON_EXTERN template void MpiBcast(const MPI_Comm&, Vector<Complex_wp >&, Vector<int64_t>&, int, int);
#endif
  

}
