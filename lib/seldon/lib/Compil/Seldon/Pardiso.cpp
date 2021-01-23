#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/direct/Pardiso.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class MatrixPardiso<double>;
  SELDON_EXTERN template class MatrixPardiso<complex<double> >;
    
  SELDON_EXTERN template void MatrixPardiso<double>::FactorizeMatrix(Matrix<double, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixPardiso<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixPardiso<double>::FactorizeMatrix(Matrix<double, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixPardiso<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixPardiso<double>::FactorizeMatrix(Matrix<double, General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixPardiso<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixPardiso<double>::FactorizeMatrix(Matrix<double, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void MatrixPardiso<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, ArrayRowSparse>&, bool);
  
  SELDON_EXTERN template void MatrixPardiso<double>::Solve(Vector<double>&);
  SELDON_EXTERN template void MatrixPardiso<double>::Solve(const SeldonTranspose&, Vector<double>&);
  SELDON_EXTERN template void MatrixPardiso<double>::Solve(const SeldonTranspose&, Matrix<double, General, ColMajor>&);

  SELDON_EXTERN template void MatrixPardiso<complex<double> >::Solve(Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixPardiso<complex<double> >::Solve(const SeldonTranspose&, Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixPardiso<complex<double> >::Solve(const SeldonTranspose&, Matrix<complex<double> , General, ColMajor>&);

  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, RowSymSparse>&, MatrixPardiso<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, RowSymSparse>&, MatrixPardiso<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, ArrayRowSymSparse>&, MatrixPardiso<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, MatrixPardiso<complex<double> >&, bool);

  SELDON_EXTERN template void GetLU(Matrix<double, General, RowSparse>&, MatrixPardiso<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, RowSparse>&, MatrixPardiso<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, General, ArrayRowSparse>&, MatrixPardiso<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, ArrayRowSparse>&, MatrixPardiso<complex<double> >&, bool);
  
  SELDON_EXTERN template void SolveLU(MatrixPardiso<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPardiso<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(MatrixPardiso<double>&, Matrix<double, General, ColMajor>&);
  
  SELDON_EXTERN template void SolveLU(MatrixPardiso<double>&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPardiso<double>&, Vector<complex<double> >&);
  
  SELDON_EXTERN template void SolveLU(MatrixPardiso<complex<double> >&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPardiso<complex<double> >&, Vector<double>&);
  
  SELDON_EXTERN template void SolveLU(MatrixPardiso<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPardiso<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(MatrixPardiso<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixPardiso<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);
  
}



