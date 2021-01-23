#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix/Matrix_Hermitian.cxx"
#include "matrix/Matrix_Symmetric.cxx"
#include "matrix/Matrix_Triangular.cxx"
#include "matrix/Matrix_Base.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template class Matrix_Hermitian<Complex_wp, Hermitian, RowHerm>;
  SELDON_EXTERN template class Matrix<Complex_wp, Hermitian, RowHerm>;
  SELDON_EXTERN template class Matrix_Symmetric<Real_wp, Symmetric, RowSym>;
  SELDON_EXTERN template class Matrix_Symmetric<Complex_wp, Symmetric, RowSym>;
  SELDON_EXTERN template class Matrix<Real_wp, Symmetric, RowSym>;
  SELDON_EXTERN template class Matrix<Complex_wp, Symmetric, RowSym>;
  SELDON_EXTERN template class Matrix_Triangular<Real_wp, General, RowUpTriang>;
  SELDON_EXTERN template class Matrix_Triangular<Complex_wp, General, RowUpTriang>;
  SELDON_EXTERN template class Matrix_Triangular<Real_wp, General, RowLoTriang>;
  SELDON_EXTERN template class Matrix_Triangular<Complex_wp, General, RowLoTriang>;
  SELDON_EXTERN template class Matrix<Real_wp, General, RowUpTriang>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, RowUpTriang>;
  SELDON_EXTERN template class Matrix<Real_wp, General, RowLoTriang>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, RowLoTriang>;

  SELDON_EXTERN template void Matrix_Hermitian<Complex_wp, Hermitian, RowHerm>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Hermitian<Complex_wp, Hermitian, RowHerm>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Hermitian<Complex_wp, Hermitian, RowHerm>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_Symmetric<Real_wp, Symmetric, RowSym>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Symmetric<Real_wp, Symmetric, RowSym>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Symmetric<Complex_wp, Symmetric, RowSym>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Symmetric<Complex_wp, Symmetric, RowSym>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Symmetric<Complex_wp, Symmetric, RowSym>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_Triangular<Real_wp, General, RowUpTriang>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Triangular<Real_wp, General, RowUpTriang>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Triangular<Complex_wp, General, RowUpTriang>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Triangular<Complex_wp, General, RowUpTriang>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Triangular<Complex_wp, General, RowUpTriang>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_Triangular<Real_wp, General, RowLoTriang>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Triangular<Real_wp, General, RowLoTriang>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Triangular<Complex_wp, General, RowLoTriang>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Triangular<Complex_wp, General, RowLoTriang>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Triangular<Complex_wp, General, RowLoTriang>::Fill(const Complex_wp&);

#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Real_wp, General, RowUpTriang>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Real_wp, General, RowLoTriang>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, General, RowLoTriang>&);
#endif

  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&,
					Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Hermitian, RowHermPacked>&,
					Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);


  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&,
					Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, Symmetric, RowSymPacked>&,
					Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&,
					Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, RowSymPacked>&,
					Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  
}
