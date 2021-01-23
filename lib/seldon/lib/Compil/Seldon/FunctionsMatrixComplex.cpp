#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "SeldonComplexMatrixHeader.hxx"
#include "SeldonComplexMatrixInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/complex/Functions_MatrixComplex.cxx"
#include "matrix_sparse/Functions_MatrixArray.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#endif

namespace Seldon
{
  
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);

  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);

  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, RowComplexSparse>&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  

  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, General, RowComplexSparse>&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);

  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const IVect&, const IVect&);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const IVect&, const IVect&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const IVect&, const IVect&);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const IVect&, const IVect&);
  
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);

  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);

  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);  
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);  
  
  SELDON_EXTERN template bool IsComplexMatrix(const Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template bool IsComplexMatrix(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template bool IsComplexMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template bool IsComplexMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Complex_wp, General, RowComplexSparse>&, const Real_wp&);
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Real_wp&);
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Real_wp&);
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Real_wp&);

  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);

  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);

  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, General, RowComplexSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  
}
