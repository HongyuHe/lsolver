#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/direct/UmfPack.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class MatrixUmfPack<double>;
  SELDON_EXTERN template class MatrixUmfPack<complex<double> >;

  SELDON_EXTERN template class MatrixUmfPack_Base<double>;
  SELDON_EXTERN template class MatrixUmfPack_Base<complex<double> >;
    
  SELDON_EXTERN template void MatrixUmfPack<double>::FactorizeMatrix(Matrix<double, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixUmfPack<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixUmfPack<double>::FactorizeMatrix(Matrix<double, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixUmfPack<complex<double> >::FactorizeMatrix(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void MatrixUmfPack<double>::FactorizeMatrix(Matrix<double, General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixUmfPack<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, RowSparse>&, bool);
  SELDON_EXTERN template void MatrixUmfPack<double>::FactorizeMatrix(Matrix<double, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void MatrixUmfPack<complex<double> >::FactorizeMatrix(Matrix<complex<double> , General, ArrayRowSparse>&, bool);

  SELDON_EXTERN template void MatrixUmfPack<double>::Solve(Vector<double>&);
  SELDON_EXTERN template void MatrixUmfPack<double>::Solve(const SeldonTranspose&, Vector<double>&);

  SELDON_EXTERN template void MatrixUmfPack<complex<double> >::Solve(Vector<complex<double> >&);
  SELDON_EXTERN template void MatrixUmfPack<complex<double> >::Solve(const SeldonTranspose&, Vector<complex<double> >&);

  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, RowSymSparse>&, MatrixUmfPack<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, RowSymSparse>&, MatrixUmfPack<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, Symmetric, ArrayRowSymSparse>&, MatrixUmfPack<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , Symmetric, ArrayRowSymSparse>&, MatrixUmfPack<complex<double> >&, bool);

  SELDON_EXTERN template void GetLU(Matrix<double, General, RowSparse>&, MatrixUmfPack<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, RowSparse>&, MatrixUmfPack<complex<double> >&, bool);
  SELDON_EXTERN template void GetLU(Matrix<double, General, ArrayRowSparse>&, MatrixUmfPack<double>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<complex<double> , General, ArrayRowSparse>&, MatrixUmfPack<complex<double> >&, bool);

  SELDON_EXTERN template void SolveLU(MatrixUmfPack<double>&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(MatrixUmfPack<double>&, Matrix<double, General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixUmfPack<double>&, Vector<double>&);
  
  SELDON_EXTERN template void SolveLU(MatrixUmfPack<double>&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixUmfPack<double>&, Vector<complex<double> >&);
  
  SELDON_EXTERN template void SolveLU(MatrixUmfPack<complex<double> >&, Vector<double>&);
  SELDON_EXTERN template void SolveLU(MatrixUmfPack<complex<double> >&, Matrix<complex<double>, General, ColMajor>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixUmfPack<complex<double> >&, Vector<double>&);
  
  SELDON_EXTERN template void SolveLU(MatrixUmfPack<complex<double> >&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, MatrixUmfPack<complex<double> >&, Vector<complex<double> >&);
  
}



