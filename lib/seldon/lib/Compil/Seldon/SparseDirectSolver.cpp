#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/direct/SparseDirectSolver.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class SparseDirectSolver<Real_wp>;
  SELDON_EXTERN template class SparseDirectSolver<Complex_wp>;
  
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::ComputeOrdering(Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::ComputeOrdering(Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::ComputeOrdering(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::ComputeOrdering(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::ComputeOrdering(Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::ComputeOrdering(Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::ComputeOrdering(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::ComputeOrdering(Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::Factorize(Matrix<Real_wp, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::Factorize(Matrix<Complex_wp, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::Factorize(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::Factorize(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::Factorize(Matrix<Real_wp, General, RowSparse>&, bool);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::Factorize(Matrix<Complex_wp, General, RowSparse>&, bool);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::Factorize(Matrix<Real_wp, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::Factorize(Matrix<Complex_wp, General, ArrayRowSparse>&, bool);

  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformAnalysis(Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformAnalysis(Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformAnalysis(Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformAnalysis(Matrix<Complex_wp, General, RowSparse>&);

    SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformAnalysis(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformAnalysis(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformAnalysis(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformAnalysis(Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformFactorization(Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformFactorization(Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformFactorization(Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformFactorization(Matrix<Complex_wp, General, RowSparse>&);

    SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformFactorization(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformFactorization(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformFactorization(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformFactorization(Matrix<Complex_wp, General, ArrayRowSparse>&);

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::FactorizeDistributed(MPI_Comm&, Vector<long>&, Vector<int>&,
										Vector<Real_wp>&, const IVect&, bool, bool);

  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::FactorizeDistributed(MPI_Comm&, Vector<long>&, Vector<int>&,
										   Vector<Complex_wp>&, const IVect&, bool, bool);    

  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::FactorizeDistributed(MPI_Comm&, Vector<int64_t>&, Vector<int64_t>&,
										Vector<Real_wp>&, const IVect&, bool, bool);

  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::FactorizeDistributed(MPI_Comm&, Vector<int64_t>&, Vector<int64_t>&,
										   Vector<Complex_wp>&, const IVect&, bool, bool);    

  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformAnalysisDistributed(MPI_Comm&, Vector<long>&, Vector<int>&,
											   Vector<Real_wp>&, const IVect&, bool, bool);

  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformAnalysisDistributed(MPI_Comm&, Vector<long>&, Vector<int>&,
											      Vector<Complex_wp>&, const IVect&, bool, bool);    

  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformAnalysisDistributed(MPI_Comm&, Vector<int64_t>&, Vector<int64_t>&,
											   Vector<Real_wp>&, const IVect&, bool, bool);

  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformAnalysisDistributed(MPI_Comm&, Vector<int64_t>&, Vector<int64_t>&,
											      Vector<Complex_wp>&, const IVect&, bool, bool);    

  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformFactorizationDistributed(MPI_Comm&, Vector<long>&, Vector<int>&,
											   Vector<Real_wp>&, const IVect&, bool, bool);

  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformFactorizationDistributed(MPI_Comm&, Vector<long>&, Vector<int>&,
											      Vector<Complex_wp>&, const IVect&, bool, bool);    

  SELDON_EXTERN template void SparseDirectSolver<Real_wp>::PerformFactorizationDistributed(MPI_Comm&, Vector<int64_t>&, Vector<int64_t>&,
											   Vector<Real_wp>&, const IVect&, bool, bool);

  SELDON_EXTERN template void SparseDirectSolver<Complex_wp>::PerformFactorizationDistributed(MPI_Comm&, Vector<int64_t>&, Vector<int64_t>&,
											      Vector<Complex_wp>&, const IVect&, bool, bool);    

  SELDON_EXTERN template void SolveLU_Distributed(MPI_Comm&, const SeldonTranspose&, SparseDirectSolver<Real_wp>&, Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void SolveLU_Distributed(MPI_Comm&, const SeldonTranspose&, SparseDirectSolver<Real_wp>&, Vector<Complex_wp>&, Vector<int>&);
  SELDON_EXTERN template void SolveLU_Distributed(MPI_Comm&, const SeldonTranspose&, SparseDirectSolver<Complex_wp>&, Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void SolveLU_Distributed(MPI_Comm&, const SeldonTranspose&, SparseDirectSolver<Complex_wp>&, Vector<Complex_wp>&, Vector<int>&);

  SELDON_EXTERN template void SolveLU_Distributed(MPI_Comm&, const SeldonTranspose&, SparseDirectSolver<Real_wp>&, Matrix<Real_wp, General, ColMajor>&, Vector<int>&);
  SELDON_EXTERN template void SolveLU_Distributed(MPI_Comm&, const SeldonTranspose&, SparseDirectSolver<Real_wp>&, Matrix<Complex_wp, General, ColMajor>&, Vector<int>&);
  SELDON_EXTERN template void SolveLU_Distributed(MPI_Comm&, const SeldonTranspose&, SparseDirectSolver<Complex_wp>&, Matrix<Real_wp, General, ColMajor>&, Vector<int>&);
  SELDON_EXTERN template void SolveLU_Distributed(MPI_Comm&, const SeldonTranspose&, SparseDirectSolver<Complex_wp>&, Matrix<Complex_wp, General, ColMajor>&, Vector<int>&);
#endif
  
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseDirectSolver<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseDirectSolver<Real_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseDirectSolver<Complex_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseDirectSolver<Complex_wp>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseDirectSolver<Real_wp>&, Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseDirectSolver<Real_wp>&, Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseDirectSolver<Complex_wp>&, Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseDirectSolver<Complex_wp>&, Matrix<Complex_wp, General, ColMajor>&);

}



