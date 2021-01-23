#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifdef SELDON_WITH_MPI
#include "SeldonDistributedHeader.hxx"
#include "SeldonDistributedInline.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/direct/Mumps.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class MatrixMumps<double>;
  SELDON_EXTERN template class MatrixMumps<complex<double> >;
    
  SELDON_EXTERN template void MatrixMumps<double>::FindOrdering(Matrix<double, Symmetric, RowSymSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::FindOrdering(Matrix<complex<double> , Symmetric, RowSymSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::FindOrdering(Matrix<double, Symmetric, ArrayRowSymSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::FindOrdering(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::FindOrdering(Matrix<double, General, RowSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::FindOrdering(Matrix<complex<double> , General, RowSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::FindOrdering(Matrix<double, General, ArrayRowSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::FindOrdering(Matrix<complex<double> , General, ArrayRowSparse>&, IVect&, bool);

  SELDON_EXTERN template void MatrixMumps<double>::FactorizeMatrix(Matrix<double, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::FactorizeMatrix(Matrix<double, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::FactorizeMatrix(Matrix<double, General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::FactorizeMatrix(Matrix<double, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, ArrayRowSparse>&, bool);

  SELDON_EXTERN template void MatrixMumps<double>::PerformAnalysis(Matrix<double, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::PerformAnalysis(Matrix<complex<double> , Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void MatrixMumps<double>::PerformAnalysis(Matrix<double, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::PerformAnalysis(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void MatrixMumps<double>::PerformAnalysis(Matrix<double, General, RowSparse>&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::PerformAnalysis(Matrix<complex<double> , General, RowSparse>&);
  SELDON_EXTERN template void MatrixMumps<double>::PerformAnalysis(Matrix<double, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::PerformAnalysis(Matrix<complex<double> , General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void MatrixMumps<double>::PerformFactorization(Matrix<double, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::PerformFactorization(Matrix<complex<double> , Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void MatrixMumps<double>::PerformFactorization(Matrix<double, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::PerformFactorization(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void MatrixMumps<double>::PerformFactorization(Matrix<double, General, RowSparse>&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::PerformFactorization(Matrix<complex<double> , General, RowSparse>&);
  SELDON_EXTERN template void MatrixMumps<double>::PerformFactorization(Matrix<double, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::PerformFactorization(Matrix<complex<double> , General, ArrayRowSparse>&);

  SELDON_EXTERN template void MatrixMumps<double>::GetSchurMatrix(Matrix<double, Symmetric, RowSymSparse>&, const IVect&, Matrix<double, General, ColMajor>&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::GetSchurMatrix(Matrix<complex<double> , Symmetric, RowSymSparse>&, const IVect&, Matrix<complex<double> , General, ColMajor>&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::GetSchurMatrix(Matrix<double, Symmetric, ArrayRowSymSparse>&, const IVect&, Matrix<double, General, ColMajor>&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::GetSchurMatrix(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, const IVect&, Matrix<complex<double> , General, ColMajor>&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::GetSchurMatrix(Matrix<double, General, RowSparse>&, const IVect&, Matrix<double, General, ColMajor>&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::GetSchurMatrix(Matrix<complex<double> , General, RowSparse>&, const IVect&, Matrix<complex<double> , General, ColMajor>&, bool);
  SELDON_EXTERN template void MatrixMumps<double>::GetSchurMatrix(Matrix<double, General, ArrayRowSparse>&, const IVect&, Matrix<double, General, ColMajor>&, bool);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::GetSchurMatrix(Matrix<complex<double> , General, ArrayRowSparse>&, const IVect&, Matrix<complex<double> , General, ColMajor>&, bool);
  
  SELDON_EXTERN template void MatrixMumps<double>::Solve(Vector<double>&);
  SELDON_EXTERN template void MatrixMumps<double>::Solve(const SeldonTranspose&, Vector<double>&);
  SELDON_EXTERN template void MatrixMumps<double>::Solve(const SeldonTranspose&, Matrix<double, General, ColMajor>&);

  SELDON_EXTERN template void MatrixMumps<complex<double> >::Solve(Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::Solve(const SeldonTranspose&, Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixMumps<complex<double> >::Solve(const SeldonTranspose&, Matrix<complex<double> , General, ColMajor>&);

  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, RowSymSparse>&, MatrixMumps<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, RowSymSparse>&, MatrixMumps<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, ArrayRowSymSparse>&, MatrixMumps<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, MatrixMumps<complex<double> >&, bool);

  SELDON_EXTERN template void GetLU(Matrix<double, General, RowSparse>&, MatrixMumps<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, RowSparse>&, MatrixMumps<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, General, ArrayRowSparse>&, MatrixMumps<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, ArrayRowSparse>&, MatrixMumps<complex<double> >&, bool);

  SELDON_EXTERN template void GetSchurMatrix(Matrix<double, Symmetric, RowSymSparse>&, MatrixMumps<double>&, const IVect&, Matrix<double, General, ColMajor>&, bool);
  SELDON_EXTERN template void GetSchurMatrix(Matrix<complex<double> , Symmetric, RowSymSparse>&, MatrixMumps<complex<double> >&, const IVect&, Matrix<complex<double> , General, ColMajor>&, bool);
  SELDON_EXTERN template void GetSchurMatrix(Matrix<double, Symmetric, ArrayRowSymSparse>&, MatrixMumps<double>&, const IVect&, Matrix<double, General, ColMajor>&, bool);
  SELDON_EXTERN template void GetSchurMatrix(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, MatrixMumps<complex<double> >&, const IVect&, Matrix<complex<double> , General, ColMajor>&, bool);

  SELDON_EXTERN template void GetSchurMatrix(Matrix<double, General, RowSparse>&, MatrixMumps<double>&, const IVect&, Matrix<double, General, ColMajor>&, bool);
  SELDON_EXTERN template void GetSchurMatrix(Matrix<complex<double> , General, RowSparse>&, MatrixMumps<complex<double> >&, const IVect&, Matrix<complex<double> , General, ColMajor>&, bool);
  SELDON_EXTERN template void GetSchurMatrix(Matrix<double, General, ArrayRowSparse>&, MatrixMumps<double>&, const IVect&, Matrix<double, General, ColMajor>&, bool);
  SELDON_EXTERN template void GetSchurMatrix(Matrix<complex<double> , General, ArrayRowSparse>&, MatrixMumps<complex<double> >&, const IVect&, Matrix<complex<double> , General, ColMajor>&, bool);
  
  SELDON_EXTERN template void SolveLU(MatrixMumps<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixMumps<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(MatrixMumps<double>&, Matrix<double, General, ColMajor>&);
  
  SELDON_EXTERN template void SolveLU(MatrixMumps<double>&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixMumps<double>&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(MatrixMumps<double>&, Matrix<complex<double> , General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixMumps<double>&, Matrix<complex<double> , General, ColMajor>&);
  
  SELDON_EXTERN template void SolveLU(MatrixMumps<complex<double> >&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixMumps<complex<double> >&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(MatrixMumps<complex<double> >&, Matrix<double, General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixMumps<complex<double> >&, Matrix<double, General, ColMajor>&);
  
  SELDON_EXTERN template void SolveLU(MatrixMumps<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixMumps<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(MatrixMumps<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixMumps<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template void MatrixMumps<double >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&, Vector<double >&, const IVect&);

  SELDON_EXTERN template void MatrixMumps<double >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&,
		     Matrix<double, General, ColMajor>&, const Vector<int>&);

  SELDON_EXTERN template void MatrixMumps<complex<double> >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&, Vector<complex<double> >&, const IVect&);

  SELDON_EXTERN template void MatrixMumps<complex<double> >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&,
		     Matrix<complex<double>, General, ColMajor>&, const Vector<int>&);

#endif
  
}



