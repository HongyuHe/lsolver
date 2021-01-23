#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifdef SELDON_WITH_MPI
#include "SeldonDistributedHeader.hxx"
#include "SeldonDistributedInline.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/solver/DistributedSolver.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class SparseDistributedSolver<Real_wp>;
  SELDON_EXTERN template class SparseDistributedSolver<Complex_wp>;

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Factorize(Matrix<Real_wp, General, RowSparse>&, bool, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Factorize(Matrix<Real_wp, Symmetric, RowSymSparse>&, bool, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Factorize(Matrix<Real_wp, General, ArrayRowSparse>&, bool, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Factorize(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformAnalysis(Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformAnalysis(Matrix<Real_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformAnalysis(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformAnalysis(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformFactorization(Matrix<Real_wp, General, RowSparse>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformFactorization(Matrix<Real_wp, Symmetric, RowSymSparse>&, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformFactorization(Matrix<Real_wp, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformFactorization(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Factorize(DistributedMatrix<Real_wp, General, RowSparse>&, bool, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Factorize(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, bool, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Factorize(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, bool, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Factorize(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool, bool);
  
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformAnalysis(DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformAnalysis(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformAnalysis(DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformAnalysis(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformFactorization(DistributedMatrix<Real_wp, General, RowSparse>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformFactorization(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformFactorization(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::PerformFactorization(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);
#endif

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::TransSolve(Vector<Real_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(const SeldonTranspose&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::TransSolve(Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(const SeldonTranspose&, Matrix<Real_wp, General, ColMajor>&);

  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(Vector<Complex_wp>&, const Vector<Complex_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(Vector<Complex_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::TransSolve(Vector<Complex_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(const SeldonTranspose&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::TransSolve(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void SparseDistributedSolver<Real_wp>::Solve(const SeldonTranspose&, Matrix<Complex_wp, General, ColMajor>&);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Factorize(Matrix<Complex_wp, General, RowSparse>&, bool, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Factorize(Matrix<Complex_wp, Symmetric, RowSymSparse>&, bool, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Factorize(Matrix<Complex_wp, General, ArrayRowSparse>&, bool, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Factorize(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformAnalysis(Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformAnalysis(Matrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformAnalysis(Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformAnalysis(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformFactorization(Matrix<Complex_wp, General, RowSparse>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformFactorization(Matrix<Complex_wp, Symmetric, RowSymSparse>&, bool);
  
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformFactorization(Matrix<Complex_wp, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformFactorization(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool);

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Factorize(DistributedMatrix<Complex_wp, General, RowSparse>&, bool, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Factorize(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, bool, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Factorize(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, bool, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Factorize(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool, bool);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformAnalysis(DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformAnalysis(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformAnalysis(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformAnalysis(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformFactorization(DistributedMatrix<Complex_wp, General, RowSparse>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformFactorization(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, bool);
  
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformFactorization(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::PerformFactorization(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool);
#endif

  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Solve(Vector<Complex_wp>&, const Vector<Complex_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Solve(Vector<Complex_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::TransSolve(Vector<Complex_wp>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Solve(const SeldonTranspose&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Solve(const SeldonTranspose&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Solve(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::TransSolve(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void SparseDistributedSolver<Complex_wp>::Solve(const SeldonTranspose&, Matrix<Complex_wp, General, ColMajor>&);
  
}



