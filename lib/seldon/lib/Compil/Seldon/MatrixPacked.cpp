#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix/Matrix_HermPacked.cxx"
#include "matrix/Matrix_SymPacked.cxx"
#include "matrix/Matrix_TriangPacked.cxx"
#include "matrix/Matrix_Base.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template class Matrix_HermPacked<Complex_wp, Hermitian, RowHermPacked>;
  SELDON_EXTERN template class Matrix<Complex_wp, Hermitian, RowHermPacked>;
  SELDON_EXTERN template class Matrix_SymPacked<Real_wp, Symmetric, RowSymPacked>;
  SELDON_EXTERN template class Matrix_SymPacked<Complex_wp, Symmetric, RowSymPacked>;
  SELDON_EXTERN template class Matrix<Real_wp, Symmetric, RowSymPacked>;
  SELDON_EXTERN template class Matrix<Complex_wp, Symmetric, RowSymPacked>;
  SELDON_EXTERN template class Matrix_TriangPacked<Real_wp, General, RowUpTriangPacked>;
  SELDON_EXTERN template class Matrix_TriangPacked<Complex_wp, General, RowUpTriangPacked>;
  SELDON_EXTERN template class Matrix_TriangPacked<Real_wp, General, RowLoTriangPacked>;
  SELDON_EXTERN template class Matrix_TriangPacked<Complex_wp, General, RowLoTriangPacked>;
  SELDON_EXTERN template class Matrix<Real_wp, General, RowUpTriangPacked>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, RowUpTriangPacked>;
  SELDON_EXTERN template class Matrix<Real_wp, General, RowLoTriangPacked>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, RowLoTriangPacked>;

  SELDON_EXTERN template void Matrix_HermPacked<Complex_wp, Hermitian, RowHermPacked>::Fill(const int&);
  SELDON_EXTERN template void Matrix_HermPacked<Complex_wp, Hermitian, RowHermPacked>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_HermPacked<Complex_wp, Hermitian, RowHermPacked>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_SymPacked<Real_wp, Symmetric, RowSymPacked>::Fill(const int&);
  SELDON_EXTERN template void Matrix_SymPacked<Real_wp, Symmetric, RowSymPacked>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_SymPacked<Complex_wp, Symmetric, RowSymPacked>::Fill(const int&);
  SELDON_EXTERN template void Matrix_SymPacked<Complex_wp, Symmetric, RowSymPacked>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_SymPacked<Complex_wp, Symmetric, RowSymPacked>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_TriangPacked<Real_wp, General, RowUpTriangPacked>::Fill(const int&);
  SELDON_EXTERN template void Matrix_TriangPacked<Real_wp, General, RowUpTriangPacked>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_TriangPacked<Complex_wp, General, RowUpTriangPacked>::Fill(const int&);
  SELDON_EXTERN template void Matrix_TriangPacked<Complex_wp, General, RowUpTriangPacked>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_TriangPacked<Complex_wp, General, RowUpTriangPacked>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_TriangPacked<Real_wp, General, RowLoTriangPacked>::Fill(const int&);
  SELDON_EXTERN template void Matrix_TriangPacked<Real_wp, General, RowLoTriangPacked>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_TriangPacked<Complex_wp, General, RowLoTriangPacked>::Fill(const int&);
  SELDON_EXTERN template void Matrix_TriangPacked<Complex_wp, General, RowLoTriangPacked>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_TriangPacked<Complex_wp, General, RowLoTriangPacked>::Fill(const Complex_wp&);
  
  SELDON_EXTERN template Matrix<Complex_wp, Hermitian, RowHermPacked>& Matrix<Complex_wp, Hermitian, RowHermPacked>::operator *=(const Complex_wp&);
  
#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Real_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Real_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, General, RowLoTriangPacked>&);
#endif

}
