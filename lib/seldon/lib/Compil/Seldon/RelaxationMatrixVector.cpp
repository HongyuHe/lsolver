#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "SeldonComplexMatrixHeader.hxx"
#include "SeldonComplexMatrixInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/Relaxation_MatVect.cxx"
#include "matrix_sparse/complex/Functions_MatVectComplex.cxx"
#include "computation/basic_functions/Functions_Base.cxx"
#endif

namespace Seldon
{
  
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, RowSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, RowSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<int, General, RowSparse>&, Vector<int>&, const Vector<int>&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<int, General, ColSparse>&, Vector<int>&, const Vector<int>&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, ColSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, ColSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, ArrayColSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, ArrayColSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<int, Symmetric, RowSymSparse>&, Vector<int>&, const Vector<int>&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<int, General, RowSparse>&, Vector<int>&, const Vector<int>&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<int, General, RowSparse>&, Vector<complex<int> >&, const Vector<complex<int> >&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<int, General, ColSparse>&, Vector<int>&, const Vector<int>&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<int, General, ColSparse>&, Vector<complex<int> >&, const Vector<complex<int> >&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ColSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayColSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayColSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayColSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<int, Symmetric, RowSymSparse>&, Vector<int>&, const Vector<int>&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<int, Symmetric, RowSymSparse>&, Vector<complex<int> >&, const Vector<complex<int> >&, const int&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, RowComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  
  //SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  //SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  //SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  //SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  
}
