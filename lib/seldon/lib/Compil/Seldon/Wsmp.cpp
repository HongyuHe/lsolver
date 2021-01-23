#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/direct/Wsmp.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class MatrixWsmp<double>;
  SELDON_EXTERN template class MatrixWsmp<complex<double> >;
    
  SELDON_EXTERN template void MatrixWsmp<double>::FactorizeMatrix(Matrix<double, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixWsmp<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixWsmp<double>::FactorizeMatrix(Matrix<double, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixWsmp<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixWsmp<double>::FactorizeMatrix(Matrix<double, General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixWsmp<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixWsmp<double>::FactorizeMatrix(Matrix<double, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void MatrixWsmp<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, ArrayRowSparse>&, bool);
  
  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, RowSymSparse>&, MatrixWsmp<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, RowSymSparse>&, MatrixWsmp<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, ArrayRowSymSparse>&, MatrixWsmp<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, MatrixWsmp<complex<double> >&, bool);

  SELDON_EXTERN template void GetLU(Matrix<double, General, RowSparse>&, MatrixWsmp<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, RowSparse>&, MatrixWsmp<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, General, ArrayRowSparse>&, MatrixWsmp<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, ArrayRowSparse>&, MatrixWsmp<complex<double> >&, bool);
  
  SELDON_EXTERN template void SolveLU(MatrixWsmp<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixWsmp<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(MatrixWsmp<double>&, Matrix<double, General, ColMajor>&);
  
  SELDON_EXTERN template void SolveLU(MatrixWsmp<double>&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixWsmp<double>&, Vector<complex<double> >&);
  
  SELDON_EXTERN template void SolveLU(MatrixWsmp<complex<double> >&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixWsmp<complex<double> >&, Vector<double>&);
  
  SELDON_EXTERN template void SolveLU(MatrixWsmp<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixWsmp<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(MatrixWsmp<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixWsmp<complex<double> >&, Matrix<complex<double> , General, ColMajor>&);
  
}



