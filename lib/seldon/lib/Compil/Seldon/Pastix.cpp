#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/direct/Pastix.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class MatrixPastix<double>;
  SELDON_EXTERN template class MatrixPastix<complex<double> >;
    
  SELDON_EXTERN template void MatrixPastix<double>::FindOrdering(Matrix<double, Symmetric, RowSymSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::FindOrdering(Matrix<complex<double> , Symmetric, RowSymSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixPastix<double>::FindOrdering(Matrix<double, Symmetric, ArrayRowSymSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::FindOrdering(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixPastix<double>::FindOrdering(Matrix<double, General, RowSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::FindOrdering(Matrix<complex<double> , General, RowSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixPastix<double>::FindOrdering(Matrix<double, General, ArrayRowSparse>&, IVect&, bool);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::FindOrdering(Matrix<complex<double> , General, ArrayRowSparse>&, IVect&, bool);

  SELDON_EXTERN template void MatrixPastix<double>::FactorizeMatrix(Matrix<double, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixPastix<double>::FactorizeMatrix(Matrix<double, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixPastix<double>::FactorizeMatrix(Matrix<double, General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixPastix<double>::FactorizeMatrix(Matrix<double, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, ArrayRowSparse>&, bool);

  SELDON_EXTERN template void MatrixPastix<double>::Solve(Vector<double>&);
  SELDON_EXTERN template void MatrixPastix<double>::Solve(const SeldonTranspose&, Vector<double>&);
  SELDON_EXTERN template void MatrixPastix<double>::Solve(const SeldonTranspose&, Matrix<double, General, ColMajor>&);

  SELDON_EXTERN template void MatrixPastix<complex<double> >::Solve(Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::Solve(const SeldonTranspose&, Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixPastix<complex<double> >::Solve(const SeldonTranspose&, Matrix<complex<double> , General, ColMajor>&);

  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, RowSymSparse>&, MatrixPastix<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, RowSymSparse>&, MatrixPastix<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, ArrayRowSymSparse>&, MatrixPastix<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, MatrixPastix<complex<double> >&, bool);

  SELDON_EXTERN template void GetLU(Matrix<double, General, RowSparse>&, MatrixPastix<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, RowSparse>&, MatrixPastix<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, General, ArrayRowSparse>&, MatrixPastix<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, ArrayRowSparse>&, MatrixPastix<complex<double> >&, bool);

  SELDON_EXTERN template void SolveLU(MatrixPastix<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPastix<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(MatrixPastix<double>&, Matrix<double, General, ColMajor>&);
  
  SELDON_EXTERN template void SolveLU(MatrixPastix<double>&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPastix<double>&, Vector<complex<double> >&);
    
  SELDON_EXTERN template void SolveLU(MatrixPastix<complex<double> >&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPastix<complex<double> >&, Vector<double>&);
    
  SELDON_EXTERN template void SolveLU(MatrixPastix<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPastix<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(MatrixPastix<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPastix<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);
  

  SELDON_EXTERN template void GetCholesky(Matrix<double, Symmetric, ArrayRowSymSparse>&, MatrixPastix<double>&, bool);
  SELDON_EXTERN template void GetCholesky(Matrix<double, Symmetric, RowSymSparse>&, MatrixPastix<double>&, bool);
  
  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, MatrixPastix<double>&, Vector<double>&);
  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, MatrixPastix<double>&, Vector<double>&);
  

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template void MatrixPastix<double >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&, Vector<double >&, const IVect&);

  SELDON_EXTERN template void MatrixPastix<double >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&,
		     Matrix<double, General, ColMajor>&, const Vector<int>&);

  SELDON_EXTERN template void MatrixPastix<complex<double> >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&, Vector<complex<double> >&, const IVect&);

  SELDON_EXTERN template void MatrixPastix<complex<double> >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&,
		     Matrix<complex<double>, General, ColMajor>&, const Vector<int>&);
  
  SELDON_EXTERN template void MatrixPastix<double >
  ::MltDistributed(MPI_Comm&, const SeldonTranspose&, Vector<double >&, const IVect&);

  SELDON_EXTERN template void MatrixPastix<double >
  ::MltDistributed(MPI_Comm&, const SeldonTranspose&,
		   Matrix<double, General, ColMajor>&, const Vector<int>&);

  SELDON_EXTERN template void MatrixPastix<complex<double> >
  ::MltDistributed(MPI_Comm&, const SeldonTranspose&, Vector<complex<double> >&, const IVect&);

  SELDON_EXTERN template void MatrixPastix<complex<double> >
  ::MltDistributed(MPI_Comm&, const SeldonTranspose&,
		     Matrix<complex<double>, General, ColMajor>&, const Vector<int>&);
  
#endif

}



