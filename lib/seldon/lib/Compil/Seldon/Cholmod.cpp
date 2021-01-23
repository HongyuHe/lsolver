#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/direct/Cholmod.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void GetCholesky(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, MatrixCholmod&, bool);
  SELDON_EXTERN template void GetCholesky(Matrix<Real_wp, Symmetric, RowSymSparse>&, MatrixCholmod&, bool);
  
  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, MatrixCholmod&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, MatrixCholmod&, Vector<Real_wp>&);
  
  SELDON_EXTERN template void MatrixCholmod::Solve(const SeldonTranspose&, Vector<Real_wp>&);
  SELDON_EXTERN template void MatrixCholmod::Mlt(const SeldonTranspose&, Vector<Real_wp>&);

  SELDON_EXTERN template void MatrixCholmod::FactorizeMatrix(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);
  
}



