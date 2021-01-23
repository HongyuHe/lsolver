#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "SeldonComplexMatrixHeader.hxx"
#include "SeldonComplexMatrixInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/complex/Matrix_ComplexSparse.cxx"
#include "matrix_sparse/complex/Matrix_SymComplexSparse.cxx"
#include "matrix_sparse/complex/Matrix_ArrayComplexSparse.cxx"
#include "matrix_sparse/Matrix_Sparse.cxx"
#include "matrix_sparse/IOMatrixMarket.cxx"
#include "matrix/Matrix_Base.cxx"
#endif

namespace Seldon
{
  
  SELDON_EXTERN template class Matrix_ComplexSparse<Complex_wp, General, RowComplexSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, RowComplexSparse>;

  SELDON_EXTERN template void Matrix_ComplexSparse<Complex_wp, General, RowComplexSparse>::SetData(int, int, Vector<Real_wp>&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, Vector<long>&, Vector<int>&);

  SELDON_EXTERN template class Matrix_SymComplexSparse<Complex_wp, Symmetric, RowSymComplexSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, Symmetric, RowSymComplexSparse>;

  SELDON_EXTERN template void Matrix_SymComplexSparse<Complex_wp, Symmetric, RowSymComplexSparse>::SetData(int, int, Vector<Real_wp>&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, Vector<long>&, Vector<int>&);

  SELDON_EXTERN template class Matrix_ArrayComplexSparse<Complex_wp, General, ArrayRowComplexSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, ArrayRowComplexSparse>;

  SELDON_EXTERN template void Matrix_ArrayComplexSparse<Complex_wp, General, ArrayRowComplexSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_ArrayComplexSparse<Complex_wp, General, ArrayRowComplexSparse>::RemoveSmallEntry(const Real_wp&);

  SELDON_EXTERN template class Matrix_ArrayComplexSparse<Complex_wp, Symmetric, ArrayRowSymComplexSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>;

  SELDON_EXTERN template void Matrix_ArrayComplexSparse<Complex_wp, Symmetric, ArrayRowSymComplexSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_ArrayComplexSparse<Complex_wp, Symmetric, ArrayRowSymComplexSparse>::RemoveSmallEntry(const Real_wp&);

  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, General, RowComplexSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, General, RowComplexSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, General, ArrayRowComplexSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, ostream&, Complex_wp&, int, bool);
  

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, General, RowComplexSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, General, RowComplexSparse>&, const string&);

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const string&);

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const string&);

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const string&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, General, RowComplexSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);


  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, General, RowComplexSparse>&, const string&);

  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const string&);

  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const string&);

  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const string&);

}
