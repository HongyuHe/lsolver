#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/Functions_MatrixArray.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#endif

namespace Seldon
{
  /* Matrix-vector functions */
  
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, General, ArrayRowSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, General, ArrayRowSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetRow(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetRow(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int, Vector<Complex_wp, VectSparse>&);
  
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, ArrayRowSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, ArrayRowSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int, Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  

  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);  
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp , General, ArrayRowSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp , General, ArrayRowSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&,
					   const Matrix<Real_wp, General, ArrayRowSparse>&,
					   const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&,
					   const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&,
					   const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);  
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&,
					   const Matrix<Complex_wp , General, ArrayRowSparse>&,
					   const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp , General, ArrayRowSparse>&,
					   const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&
);


  /* ArrayColSparse */

  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ArrayColSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);  
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ArrayColSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, ArrayColSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayColSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayColSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayColSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ArrayColSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ArrayColSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, ArrayColSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp , General, ArrayColSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp , General, ArrayColSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&,
					   const Matrix<Real_wp, General, ArrayColSparse>&,
					   const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ArrayColSparse>&,
					   const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ArrayColSparse>&,
					   const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);  
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&,
					   const Matrix<Complex_wp , General, ArrayColSparse>&,
					   const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp , General, ArrayColSparse>&,
					   const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  /* ArrayRowSymSparse */
  
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);  
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
    
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp , Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp , Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&,
                                     const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
                                     const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
                                     const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
                                     const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&,
                                     const Matrix<Complex_wp , Symmetric, ArrayRowSymSparse>&,
                                     const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp , Symmetric, ArrayRowSymSparse>&,
                                     const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  
  /* Matrix functions */
  
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
    
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
    
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
    
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  

  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
                                  const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&,
					const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Matrix<Real_wp, General, RowMajor>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, ArrayRowSparse>&, const Matrix<Complex_wp, General, RowMajor>&, const Complex_wp&, Matrix<Complex_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, ArrayRowSparse>&, const Complex_wp&, Matrix<Complex_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, ArrayRowSparse>&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Complex_wp&, Matrix<Complex_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Matrix<Real_wp, General, ArrayRowSparse>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Matrix<Complex_wp, General, ArrayRowSparse>&, const Complex_wp&, Matrix<Complex_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Real_wp, General, RowSparse>&, const Real_wp&);
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Real_wp, Symmetric, RowSymSparse>&, const Real_wp&); 
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Real_wp, General, ArrayRowSparse>&, const Real_wp&);
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Real_wp&);

  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Complex_wp, General, RowSparse>&, const Real_wp&);
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Real_wp&); 
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Complex_wp, General, ArrayRowSparse>&, const Real_wp&);
  SELDON_EXTERN template void RemoveSmallEntry(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Real_wp&);

  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Real_wp, General, RowSparse>&, const IVect&, const IVect&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const IVect&, const IVect&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const IVect&, const IVect&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const IVect&, const IVect&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, General, RowSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const IVect&, const IVect&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const IVect&, const IVect&, Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopySubMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const IVect&, const IVect&, Matrix<Complex_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void MltMatrix(Matrix<Real_wp, General, ArrayRowSparse> const&, Matrix<Real_wp, General, ArrayRowSparse> const&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MltAddMatrix(Real_wp const&, Matrix<Real_wp, General, ArrayRowSparse> const&, Matrix<Real_wp, General, ArrayRowSparse> const&, Real_wp const&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MltAddMatrix(Real_wp const&, SeldonTranspose const&, Matrix<Real_wp, General, ArrayRowSparse> const&, SeldonTranspose const&, Matrix<Real_wp, General, ArrayRowSparse> const&, Real_wp const&, Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void MltAddMatrix(Complex_wp const&, SeldonTranspose const&, Matrix<Complex_wp, General, ArrayRowSparse> const&, SeldonTranspose const&, Matrix<Complex_wp, General, ArrayRowSparse> const&, Complex_wp const&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MltAddMatrix(Complex_wp const&, Matrix<Complex_wp, General, ArrayRowSparse> const&, Matrix<Complex_wp, General, ArrayRowSparse> const&, Complex_wp const&, Matrix<Complex_wp, General, ArrayRowSparse>&);


  // for integer matrices

  SELDON_EXTERN template void MltVector(const Matrix<int, General, RowSparse>&,
					const Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<int, General, RowSparse>&,
					const Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void MltVector(const Matrix<int, General, RowSparse>&,
					const Vector<complex<int> >&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<int, General, RowSparse>&,
					const Vector<complex<int> >&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltAddVector(const int&, const Matrix<int, General, RowSparse>&,
					   const Vector<int>&, const int&, Vector<int>&);

  SELDON_EXTERN template void MltAddVector(const int&, const SeldonTranspose&, const Matrix<int, General, RowSparse>&,
					   const Vector<int>&, const int&, Vector<int>&);

  SELDON_EXTERN template void MltAddVector(const complex<int>&, const Matrix<int, General, RowSparse>&,
					   const Vector<complex<int> >&, const complex<int>&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltAddVector(const complex<int>&, const SeldonTranspose&, const Matrix<int, General, RowSparse>&,
					   const Vector<complex<int> >&, const complex<int>&, Vector<complex<int> >&);


  SELDON_EXTERN template void MltVector(const Matrix<int, General, ColSparse>&,
					const Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<int, General, ColSparse>&,
					const Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void MltVector(const Matrix<int, General, ColSparse>&,
					const Vector<complex<int> >&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<int, General, ColSparse>&,
					const Vector<complex<int> >&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltAddVector(const int&, const Matrix<int, General, ColSparse>&,
					   const Vector<int>&, const int&, Vector<int>&);

  SELDON_EXTERN template void MltAddVector(const int&, const SeldonTranspose&, const Matrix<int, General, ColSparse>&,
					   const Vector<int>&, const int&, Vector<int>&);

  SELDON_EXTERN template void MltAddVector(const complex<int>&, const Matrix<int, General, ColSparse>&,
					   const Vector<complex<int> >&, const complex<int>&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltAddVector(const complex<int>&, const SeldonTranspose&, const Matrix<int, General, ColSparse>&,
					   const Vector<complex<int> >&, const complex<int>&, Vector<complex<int> >&);


  SELDON_EXTERN template void MltVector(const Matrix<int, Symmetric, RowSymSparse>&,
					const Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<int, Symmetric, RowSymSparse>&,
					const Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void MltVector(const Matrix<int, Symmetric, RowSymSparse>&,
					const Vector<complex<int> >&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<int, Symmetric, RowSymSparse>&,
					const Vector<complex<int> >&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltAddVector(const int&, const Matrix<int, Symmetric, RowSymSparse>&,
					   const Vector<int>&, const int&, Vector<int>&);

  SELDON_EXTERN template void MltAddVector(const int&, const SeldonTranspose&, const Matrix<int, Symmetric, RowSymSparse>&,
					   const Vector<int>&, const int&, Vector<int>&);

  SELDON_EXTERN template void MltAddVector(const complex<int>&, const Matrix<int, Symmetric, RowSymSparse>&,
					   const Vector<complex<int> >&, const complex<int>&, Vector<complex<int> >&);

  SELDON_EXTERN template void MltAddVector(const complex<int>&, const SeldonTranspose&, const Matrix<int, Symmetric, RowSymSparse>&,
					   const Vector<complex<int> >&, const complex<int>&, Vector<complex<int> >&);
  
}
