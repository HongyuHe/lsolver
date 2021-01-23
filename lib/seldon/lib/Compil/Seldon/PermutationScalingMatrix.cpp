#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/Permutation_ScalingMatrix.cxx"
#endif

namespace Seldon
{
  
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, General, RowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, General, RowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, General, ColSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, General, ColSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, const Vector<int>&);

  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, General, RowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, General, RowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, General, ColSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, General, ColSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, RowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, ColSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, ColSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, ColSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, ColSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, ColSparse>&, const Vector<Complex_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&);
  
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, RowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, ColSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, ColSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&);
    
}
