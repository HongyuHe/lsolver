#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "matrix_sparse/IOMatrixMarket.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/IOMatrixMarket.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template void ReadCoordinateMatrix(istream&, Vector<int>&, Vector<int>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(istream&, Vector<int>&, Vector<int>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(istream&, Vector<int>&, Vector<int>&, Vector<int>&, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(ostream&, const Vector<int>&, const Vector<int>&, const Vector<Real_wp>&, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(ostream&, const Vector<int>&, const Vector<int>&, const Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(ostream&, const Vector<int>&, const Vector<int>&, const Vector<int>&, bool);
  
  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Real_wp, General, RowSparse>&, istream&, Real_wp&, int, long, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, General, RowSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Real_wp, General, RowSparse>&, ostream&, Real_wp&, int, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, General, RowSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<int, General, RowSparse>&, istream&, int&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<int, General, RowSparse>&, ostream&, int&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<int, General, ColSparse>&, istream&, int&, int, long, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Real_wp, General, ColSparse>&, istream&, Real_wp&, int, long, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, General, ColSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<int, General, ColSparse>&, ostream&, int&, int, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Real_wp, General, ColSparse>&, ostream&, Real_wp&, int, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, General, ColSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Real_wp, General, ArrayRowSparse>&, istream&, Real_wp&, int, long, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, General, ArrayRowSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, ostream&, Real_wp&, int, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Real_wp, General, ArrayColSparse>&, istream&, Real_wp&, int, long, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, General, ArrayColSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Real_wp, General, ArrayColSparse>&, ostream&, Real_wp&, int, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, General, ArrayColSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Real_wp, Symmetric, RowSymSparse>&, istream&, Real_wp&, int, long, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, Symmetric, RowSymSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Real_wp, Symmetric, RowSymSparse>&, ostream&, Real_wp&, int, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<int, Symmetric, RowSymSparse>&, istream&, int&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<int, Symmetric, RowSymSparse>&, ostream&, int&, int, bool);

  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, istream&, Real_wp&, int, long, bool);
  SELDON_EXTERN template void ReadCoordinateMatrix(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, istream&, Complex_wp&, int, long, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, ostream&, Real_wp&, int, bool);
  SELDON_EXTERN template void WriteCoordinateMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, ostream&, Complex_wp&, int, bool);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Real_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Real_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  ReadHarwellBoeing(string, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Real_wp, General, RowSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, General, RowSparse>&, const string&);

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Real_wp, General, ColSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, General, ColSparse>&, const string&);

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Real_wp, General, ArrayRowSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, General, ArrayRowSparse>&, const string&);

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const string&);

  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteHarwellBoeing(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const string&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Real_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, General, ColSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Real_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  ReadMatrixMarket(string, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Real_wp, General, RowSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, General, RowSparse>&, const string&);

  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Real_wp, General, ColSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, General, ColSparse>&, const string&);

  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Real_wp, General, ArrayRowSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, General, ArrayRowSparse>&, const string&);

  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const string&);

  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const string&);
  
  SELDON_EXTERN template void
  WriteMatrixMarket(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const string&);
  
}
