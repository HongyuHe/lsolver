#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/Matrix_Conversions.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Real_wp, General, RowSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Real_wp>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, General, RowSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<int, General, RowSparse>&,
                               Vector<int>&, Vector<int>&, Vector<int>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Real_wp, General, ColSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Real_wp>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, General, ColSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<int, General, ColSparse>&,
                               Vector<int>&, Vector<int>&, Vector<int>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Real_wp>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, Symmetric, RowSymSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<int, Symmetric, RowSymSparse>&,
                               Vector<int>&, Vector<int>&, Vector<int>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Real_wp, General, ArrayRowSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Real_wp>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Real_wp, General, ArrayColSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Real_wp>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, General, ArrayRowSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, General, ArrayColSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Real_wp>&, int, bool);
  
  SELDON_EXTERN template void
  ConvertMatrix_to_Coordinates(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
                               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, int, bool);

  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Real_wp>&,
                                 Matrix<Real_wp, General, RowSparse>&, int);

  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<int>&,
                                 Matrix<int, General, RowSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, General, RowSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Real_wp>&,
                                 Matrix<Real_wp, General, ColSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, General, ColSparse>&, int);

  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<int>&,
                                 Matrix<int, General, ColSparse>&, int);

  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Real_wp>&,
                                 Matrix<Real_wp, Symmetric, RowSymSparse>&, int);

  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<int>&,
                                 Matrix<int, Symmetric, RowSymSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, Symmetric, RowSymSparse>&, int);

  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Real_wp>&,
                                 Matrix<Real_wp, General, ArrayRowSparse>&, int);

  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Real_wp>&,
                                 Matrix<Real_wp, General, ArrayColSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, General, ArrayRowSparse>&, int);

  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, General, ArrayColSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Real_wp>&,
                                 Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int);
  
  SELDON_EXTERN template void
  ConvertMatrix_from_Coordinates(const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&,
                                 Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, General, RowSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, RowSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, General, ArrayRowSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, ArrayRowSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, bool);
  
  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, General, ColSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, General, ColSparse>&, General&,
               Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, ColSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, ColSparse>&, General&,
               Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, RowSymSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, RowSymSparse>&, General&,
               Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Symmetric&,
               Vector<int>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Symmetric&,
               Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool);
  
  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, General&,
               Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Symmetric&,
               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Symmetric&,
               Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, General&,
               Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&,
               Vector<int>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&,
               Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, General&,
               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, General&,
               Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&,
               Vector<int>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&,
               Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, General, RowSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, RowSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, General, ArrayRowSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, ArrayRowSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool);
  
  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, General, ColSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, General, ColSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, RowSymSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Symmetric&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool);
  
  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Symmetric&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, General&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void
  ConvertToCSC(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&,
               Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool);


  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, RowSparse>&,
       Matrix<Real_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, RowSparse>&,
       Matrix<Complex_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&,
       Matrix<Real_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&,
       Matrix<Complex_wp, General, ColSparse>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
       Matrix<Real_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&,
       Matrix<Complex_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Real_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Real_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Complex_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&,
       Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, RowSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, RowSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, RowSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, RowSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, ColSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, ColSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, ColSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, ColSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, ArrayRowSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, ArrayRowSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, ArrayRowSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, ArrayRowSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, RowSymSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, RowSymSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Symmetric&, Vector<int>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Symmetric&, Vector<long>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Symmetric&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Symmetric&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, General&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, General&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, Vector<int>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, Vector<long>&,
               Vector<int>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, Vector<int>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, Vector<long>&,
               Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, RowSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, RowSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, ColSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, ColSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, General, ArrayRowSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, General, ArrayRowSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, RowSymSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Symmetric&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Symmetric&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, General&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  ConvertToCSR(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, Vector<int64_t>&,
               Vector<int64_t>&, Vector<Complex_wp>&);


  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
       Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&,
       Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ColSparse>&,
       Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, ColSparse>&,
       Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&,
       Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&,
       Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
       Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&,
       Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, RowSparse>&,
       Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, RowSparse>&,
       Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
       Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&,
       Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ColSparse>&,
       Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, ColSparse>&,
       Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&,
       Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&,
       Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&,
       Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, RowSparse>&,
       Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, RowSparse>&,
       Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ColSparse>&,
       Matrix<Real_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, ColSparse>&,
       Matrix<Complex_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
       Matrix<Real_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&,
       Matrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
       Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, RowSparse>&,
	     Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, RowSparse>&,
	     Matrix<Complex_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&,
	     Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&,
	     Matrix<Complex_wp, General, RowMajor>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, RowMajor>&,
	     Matrix<Real_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, General, RowMajor>&,
	     Matrix<Complex_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, RowSymPacked>&,
	     Matrix<Real_wp, Symmetric, RowSymPacked>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, RowSymPacked>&,
	     Matrix<Complex_wp, Symmetric, RowSymPacked>&);

  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
	     Matrix<Real_wp, Symmetric, RowSymPacked>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
	     Matrix<Real_wp, Symmetric, RowSymPacked>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
	     Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
	     Matrix<Real_wp, Symmetric, RowSym>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
	     Matrix<Complex_wp, Symmetric, RowSym>&);
  
  SELDON_EXTERN template void
  CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&,
	     Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, General, RowSparse>&,
                      Matrix<int, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<int, General, RowSparse>&,
                      Matrix<int, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
                      Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
                      Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
                      Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Complex_wp, Symmetric, RowSymSparse>&,
                      Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, General, RowSparse>&,
                      Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Complex_wp, General, RowSparse>&,
                      Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, General, ArrayRowSparse>&,
                      Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Complex_wp, General, ArrayRowSparse>&,
                      Vector<int>&, Vector<int>&);


  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
                      Vector<int64_t>&, Vector<int64_t>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
                      Vector<int64_t>&, Vector<int64_t>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, Symmetric, RowSymSparse>&,
                      Vector<int64_t>&, Vector<int64_t>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Complex_wp, Symmetric, RowSymSparse>&,
                      Vector<int64_t>&, Vector<int64_t>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, General, RowSparse>&,
                      Vector<int64_t>&, Vector<int64_t>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Complex_wp, General, RowSparse>&,
                      Vector<int64_t>&, Vector<int64_t>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Real_wp, General, ArrayRowSparse>&,
                      Vector<int64_t>&, Vector<int64_t>&);

  SELDON_EXTERN template void
  GetSymmetricPattern(const Matrix<Complex_wp, General, ArrayRowSparse>&,
                      Vector<int64_t>&, Vector<int64_t>&);

  SELDON_EXTERN template void 
  ConvertToSparse(const Matrix<Real_wp, Symmetric, RowSymPacked>&,
                  Matrix<Real_wp, Symmetric, RowSymSparse>&, const Real_wp&);
  
  SELDON_EXTERN template void 
  ConvertToSparse(const Matrix<Real_wp, General, RowMajor>&,
                  Matrix<Real_wp, General, ArrayRowSparse>&, const Real_wp&);

  SELDON_EXTERN template void 
  ConvertToSparse(const Matrix<Real_wp, General, RowMajor>&,
                  Matrix<Real_wp, General, RowSparse>&, const Real_wp&);



  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&,
					 Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
}
