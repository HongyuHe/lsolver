#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "SeldonComplexMatrixHeader.hxx"
#include "SeldonComplexMatrixInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/complex/Matrix_ComplexConversions.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, General, RowComplexSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);
  
  //SELDON_EXTERN template void
  //ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, General, ArrayColComplexSparse>&,
  //                             Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, General, RowComplexSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, General, ArrayRowComplexSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, RowComplexSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Symmetric&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Symmetric&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, RowComplexSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Symmetric&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Symmetric&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Symmetric&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Symmetric&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ColSparse>&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ColSparse>&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, RowComplexSparse>&, Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, RowComplexSparse>&, Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, RowComplexSparse>&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, RowComplexSparse>&, Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  
}
