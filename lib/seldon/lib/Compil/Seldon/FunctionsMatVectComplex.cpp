#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "SeldonComplexMatrixHeader.hxx"
#include "SeldonComplexMatrixInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/complex/Functions_MatVectComplex.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, RowComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, const Real_wp&, int, int);
  
}
