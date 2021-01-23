#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix/Matrix_Pointers.cxx"
#include "matrix/Matrix_Base.cxx"
#endif


namespace Seldon
{

  SELDON_EXTERN template class VirtualMatrix<bool>;
  SELDON_EXTERN template class VirtualMatrix<int>;
  SELDON_EXTERN template class VirtualMatrix<float>;
  SELDON_EXTERN template class VirtualMatrix<Real_wp>;
  SELDON_EXTERN template class VirtualMatrix<complex<float> >;
  SELDON_EXTERN template class VirtualMatrix<Complex_wp>;
  SELDON_EXTERN template class VirtualMatrix<Matrix<Real_wp> >;
  
  SELDON_EXTERN template class Matrix_Pointers<int, General, ColMajor>;
  SELDON_EXTERN template class Matrix_Pointers<Real_wp, General, ColMajor>;
  SELDON_EXTERN template class Matrix_Pointers<Complex_wp, General, ColMajor>;
  SELDON_EXTERN template class Matrix_Pointers<bool, General, RowMajor>;
  SELDON_EXTERN template class Matrix_Pointers<int, General, RowMajor>;
  SELDON_EXTERN template class Matrix_Pointers<float, General, RowMajor>;
  SELDON_EXTERN template class Matrix_Pointers<Real_wp, General, RowMajor>;
  SELDON_EXTERN template class Matrix_Pointers<complex<float>, General, RowMajor>;
  SELDON_EXTERN template class Matrix_Pointers<Complex_wp, General, RowMajor>;
  SELDON_EXTERN template class Matrix<int, General, ColMajor>;
  SELDON_EXTERN template class Matrix<Real_wp, General, ColMajor>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, ColMajor>;
  SELDON_EXTERN template class Matrix<int, General, RowMajor>;
  SELDON_EXTERN template class Matrix<float, General, RowMajor>;
  SELDON_EXTERN template class Matrix<Real_wp, General, RowMajor>;
  SELDON_EXTERN template class Matrix<complex<float>, General, RowMajor>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, RowMajor>;
  
  SELDON_EXTERN template void Matrix_Pointers<Real_wp, General, ColMajor>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Pointers<Complex_wp, General, ColMajor>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Pointers<Complex_wp, General, ColMajor>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_Pointers<Real_wp, General, RowMajor>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Pointers<Complex_wp, General, RowMajor>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Pointers<Complex_wp, General, RowMajor>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_Pointers<int, General, ColMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<Real_wp, General, ColMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<Complex_wp, General, ColMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<bool, General, RowMajor>::Fill(const bool&);
  SELDON_EXTERN template void Matrix_Pointers<int, General, RowMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<Real_wp, General, RowMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<Complex_wp, General, RowMajor>::Fill(const int&);
  
#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<int, General, ColMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<int, General, RowMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<bool, General, RowMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<Complex_wp, General, RowMajor>&);
#endif

#ifdef SELDON_WITH_MULTIPLE
  SELDON_EXTERN template class VirtualMatrix<double>;
  SELDON_EXTERN template class VirtualMatrix<complex<double> >;

  SELDON_EXTERN template class Matrix_Pointers<double, General, ColMajor>;
  SELDON_EXTERN template class Matrix_Pointers<complex<double>, General, ColMajor>;
  SELDON_EXTERN template class Matrix_Pointers<double, General, RowMajor>;
  SELDON_EXTERN template class Matrix_Pointers<complex<double>, General, RowMajor>;

  SELDON_EXTERN template void Matrix_Pointers<double, General, ColMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<double, General, ColMajor>::Fill(const double&);
  SELDON_EXTERN template void Matrix_Pointers<complex<double>, General, ColMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<complex<double>, General, ColMajor>::Fill(const double&);
  SELDON_EXTERN template void Matrix_Pointers<complex<double>, General, ColMajor>::Fill(const complex<double>&);

  SELDON_EXTERN template void Matrix_Pointers<double, General, RowMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<double, General, RowMajor>::Fill(const double&);
  SELDON_EXTERN template void Matrix_Pointers<complex<double>, General, RowMajor>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Pointers<complex<double>, General, RowMajor>::Fill(const double&);
  SELDON_EXTERN template void Matrix_Pointers<complex<double>, General, RowMajor>::Fill(const complex<double>&);

  SELDON_EXTERN template void Matrix_Pointers<Real_wp, General, RowMajor>::Fill(const double&);

  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<double, General, ColMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<complex<double>, General, ColMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<double, General, RowMajor>&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<complex<double>, General, RowMajor>&);
  
#endif
  
}
