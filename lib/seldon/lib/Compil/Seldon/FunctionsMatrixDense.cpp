#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix/Functions.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, RowSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, RowSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, ColSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, ColSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, Symmetric, RowSymSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, int, Vector<Complex_wp, VectSparse>&);
  
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, RowSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, RowSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, ColSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, ColSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, Symmetric, RowSymSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, int, Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, RowMajor>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, RowMajor>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, ColMajor>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, ColMajor>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, Symmetric, RowSymPacked>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, Symmetric, RowSym>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, Symmetric, RowSym>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, Hermitian, RowHerm>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, RowLoTriang>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, RowLoTriang>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, RowUpTriang>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, RowUpTriang>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, RowLoTriangPacked>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, RowLoTriangPacked>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, RowUpTriangPacked>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, RowUpTriangPacked>&, int, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, RowMajor>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, RowMajor>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, ColMajor>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, ColMajor>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, Symmetric, RowSymPacked>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, Symmetric, RowSym>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Symmetric, RowSym>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Hermitian, RowHerm>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, RowLoTriang>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, RowLoTriang>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, RowUpTriang>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, RowUpTriang>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, RowLoTriangPacked>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, RowLoTriangPacked>&, int, Vector<Complex_wp >&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, RowUpTriangPacked>&, int, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, RowUpTriangPacked>&, int, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void SetRow(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp>&, int, Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp>&, int, Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp>&, int, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowUpTriangPacked>&);

  SELDON_EXTERN template void SetCol(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp>&, int, Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp>&, int, Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp>&, int, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp>&, int, Matrix<Real_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp >&, int, Matrix<Complex_wp, General, RowUpTriangPacked>&);
  
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, General, RowMajor>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, General, RowMajor>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, General, ColMajor>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, General, ColMajor>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Real_wp, Symmetric, RowSym>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyInversePermutation(Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<int>&, const Vector<int>&, int);

  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, General, RowMajor>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, General, RowMajor>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, General, ColMajor>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, General, ColMajor>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Real_wp, Symmetric, RowSym>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<int>&, const Vector<int>&, int);
  SELDON_EXTERN template void ApplyPermutation(Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<int>&, const Vector<int>&, int);

  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, RowMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, ColMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>&);

  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowMajor>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, ColMajor>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&);

  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp >&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp >&);

  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, RowMajor>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, ColMajor>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&);

  
}
