#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "SeldonComplexMatrixHeader.hxx"
#include "SeldonComplexMatrixInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/basic_functions/Functions_Vector.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#endif

namespace Seldon
{
  
  SELDON_EXTERN template void CheckDim(const Vector<Real_wp>&, const Vector<Real_wp>&, string, string);
  SELDON_EXTERN template void CheckDim(const Vector<Complex_wp>&, const Vector<Complex_wp>&, string, string);
  SELDON_EXTERN template void CheckDim(const Vector<Real_wp, VectSparse>&, const Vector<Real_wp, VectSparse>&, string, string);
  SELDON_EXTERN template void CheckDim(const Vector<Complex_wp, VectSparse>&, const Vector<Complex_wp, VectSparse>&, string, string);

  // matrix-vector check

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, ColMajor>&, const Vector<Real_wp>&, const Vector<Real_wp>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, string, string); 

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp>&, string, string); 

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, string, string); 

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp>&, string, string); 

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp>&, string, string); 

  // sparse matrix with vectors

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, string, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, const Vector<Complex_wp>& Y, string, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, string, string);


  // matrix-matrix check
  

  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);


  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, ColMajor>&, const Matrix<Real_wp, General, ColMajor>&, const Matrix<Real_wp, General, ColMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, ColMajor>&, const Matrix<Real_wp, General, ColMajor>&, const Matrix<Real_wp, General, ColMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, ColMajor>&, const Matrix<Real_wp, General, ColMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, ColMajor>&, const Matrix<Real_wp, General, ColMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Matrix<Real_wp, General, ColMajor>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, ColMajor>&, const Matrix<Complex_wp, General, ColMajor>&, const Matrix<Complex_wp, General, ColMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, ColMajor>&, const Matrix<Complex_wp, General, ColMajor>&, const Matrix<Complex_wp, General, ColMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, ColMajor>&, const Matrix<Complex_wp, General, ColMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, ColMajor>&, const Matrix<Complex_wp, General, ColMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Matrix<Complex_wp, General, ColMajor>&, string);


  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, Hermitian, RowHerm>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, Hermitian, RowHerm>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Matrix<Complex_wp, General, RowMajor>&, string);


  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSym>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, Symmetric, RowSym>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, Symmetric, RowSym>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, Symmetric, RowSym>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSym>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, Symmetric, RowSym>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, Symmetric, RowSym>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, Symmetric, RowSym>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSym>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSym>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);


  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriang>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowLoTriang>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriang>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowLoTriang>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowLoTriang>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowLoTriang>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowLoTriang>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowLoTriang>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriang>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriang>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);


  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriangPacked>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowLoTriangPacked>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);



  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriang>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowUpTriang>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriang>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowUpTriang>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowUpTriang>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowUpTriang>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowUpTriang>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowUpTriang>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriang>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriang>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);


  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriangPacked>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Real_wp, General, RowUpTriangPacked>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, string);

  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>& Y, string);
  SELDON_EXTERN template void CheckDim(const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonSide&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, string);
  SELDON_EXTERN template void CheckDim(const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, string);

}
