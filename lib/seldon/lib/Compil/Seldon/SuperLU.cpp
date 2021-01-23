#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/direct/SuperLU.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class MatrixSuperLU<double>;
  SELDON_EXTERN template class MatrixSuperLU<complex<double> >;

  SELDON_EXTERN template class MatrixSuperLU_Base<double>;
  SELDON_EXTERN template class MatrixSuperLU_Base<complex<double> >;
    
  SELDON_EXTERN template void MatrixSuperLU<double>::FactorizeMatrix(Matrix<double, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixSuperLU<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixSuperLU<double>::FactorizeMatrix(Matrix<double, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixSuperLU<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixSuperLU<double>::FactorizeMatrix(Matrix<double, General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixSuperLU<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixSuperLU<double>::FactorizeMatrix(Matrix<double, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void MatrixSuperLU<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, ArrayRowSparse>&, bool);
  
  SELDON_EXTERN template void MatrixSuperLU<double>::Solve(Vector<double>&);
  SELDON_EXTERN template void MatrixSuperLU<double>::Solve(const SeldonTranspose&, Vector<double>&);
  SELDON_EXTERN template void MatrixSuperLU<double>::Solve(const SeldonTranspose&, Matrix<double, General, ColMajor>&);

  SELDON_EXTERN template void MatrixSuperLU<complex<double> >::Solve(Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixSuperLU<complex<double> >::Solve(const SeldonTranspose&, Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixSuperLU<complex<double> >::Solve(const SeldonTranspose&, Matrix<complex<double> , General, ColMajor>&);

  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, RowSymSparse>&, MatrixSuperLU<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, RowSymSparse>&, MatrixSuperLU<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, ArrayRowSymSparse>&, MatrixSuperLU<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, MatrixSuperLU<complex<double> >&, bool);

  SELDON_EXTERN template void GetLU(Matrix<double, General, RowSparse>&, MatrixSuperLU<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, RowSparse>&, MatrixSuperLU<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, General, ArrayRowSparse>&, MatrixSuperLU<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, ArrayRowSparse>&, MatrixSuperLU<complex<double> >&, bool);
  
  SELDON_EXTERN template void SolveLU(MatrixSuperLU<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixSuperLU<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(MatrixSuperLU<double>&, Matrix<double, General, ColMajor>&);
  
  SELDON_EXTERN template void SolveLU(MatrixSuperLU<double>&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixSuperLU<double>&, Vector<complex<double> >&);
  
  SELDON_EXTERN template void SolveLU(MatrixSuperLU<complex<double> >&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixSuperLU<complex<double> >&, Vector<double>&);
  
  SELDON_EXTERN template void SolveLU(MatrixSuperLU<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixSuperLU<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(MatrixSuperLU<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixSuperLU<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);


#ifdef SELDON_WITH_SUPERLU_DIST
  SELDON_EXTERN template void MatrixSuperLU<double >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&, Vector<double >&, const IVect&);

  SELDON_EXTERN template void MatrixSuperLU<double >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&, Matrix<double, General, ColMajor>&, const IVect&);

  SELDON_EXTERN template void MatrixSuperLU<complex<double> >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&, Vector<complex<double> >&, const IVect&);

  SELDON_EXTERN template void MatrixSuperLU<complex<double> >
  ::SolveDistributed(MPI_Comm&, const SeldonTranspose&, Matrix<complex<double>, General, ColMajor>&, const IVect&);
#endif
  
}



