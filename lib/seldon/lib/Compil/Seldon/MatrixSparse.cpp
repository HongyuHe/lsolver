#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/Matrix_Sparse.cxx"
#include "matrix_sparse/Matrix_SymSparse.cxx"
#include "matrix_sparse/Matrix_ArraySparse.cxx"
#include "matrix/Matrix_Base.cxx"
#endif

namespace Seldon
{
  
  SELDON_EXTERN template class Matrix_Sparse<int, General, RowSparse>;
  SELDON_EXTERN template class Matrix_Sparse<int, General, ColSparse>;
  SELDON_EXTERN template class Matrix_Sparse<Real_wp, General, RowSparse>;
  SELDON_EXTERN template class Matrix_Sparse<Real_wp, General, ColSparse>;
  SELDON_EXTERN template class Matrix_Sparse<Complex_wp, General, RowSparse>;
  SELDON_EXTERN template class Matrix_Sparse<Complex_wp, General, ColSparse>;
  SELDON_EXTERN template class Matrix<int, General, RowSparse>;
  SELDON_EXTERN template class Matrix<int, General, ColSparse>;
  SELDON_EXTERN template class Matrix<Real_wp, General, RowSparse>;
  SELDON_EXTERN template class Matrix<Real_wp, General, ColSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, RowSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, ColSparse>;

  SELDON_EXTERN template void Matrix_Sparse<int, General, RowSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Sparse<int, General, ColSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Sparse<Real_wp, General, RowSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Sparse<Real_wp, General, RowSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Sparse<Complex_wp, General, RowSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Sparse<Complex_wp, General, RowSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Sparse<Complex_wp, General, RowSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_Sparse<Real_wp, General, ColSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Sparse<Real_wp, General, ColSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Sparse<Complex_wp, General, ColSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Sparse<Complex_wp, General, ColSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Sparse<Complex_wp, General, ColSparse>::Fill(const Complex_wp&);
  
  SELDON_EXTERN template void Matrix_Sparse<int, General, RowSparse>::SetData(int, int, Vector<int>&, Vector<long>&, Vector<int>&);
  SELDON_EXTERN template void Matrix_Sparse<int, General, ColSparse>::SetData(int, int, Vector<int>&, Vector<long>&, Vector<int>&);
  SELDON_EXTERN template void Matrix_Sparse<Real_wp, General, RowSparse>::SetData(int, int, Vector<Real_wp>&, Vector<long>&, Vector<int>&);
  SELDON_EXTERN template void Matrix_Sparse<Complex_wp, General, RowSparse>::SetData(int, int, Vector<Complex_wp>&, Vector<long>&, Vector<int>&);
  SELDON_EXTERN template void Matrix_Sparse<Real_wp, General, ColSparse>::SetData(int, int, Vector<Real_wp>&, Vector<long>&, Vector<int>&);
  SELDON_EXTERN template void Matrix_Sparse<Complex_wp, General, ColSparse>::SetData(int, int, Vector<Complex_wp>&, Vector<long>&, Vector<int>&);

  SELDON_EXTERN template class Matrix_SymSparse<int, Symmetric, RowSymSparse>;
  SELDON_EXTERN template class Matrix_SymSparse<Real_wp, Symmetric, RowSymSparse>;
  SELDON_EXTERN template class Matrix_SymSparse<Complex_wp, Symmetric, RowSymSparse>;
  SELDON_EXTERN template class Matrix<Real_wp, Symmetric, RowSymSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, Symmetric, RowSymSparse>;

  SELDON_EXTERN template void Matrix_SymSparse<int, Symmetric, RowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_SymSparse<Real_wp, Symmetric, RowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_SymSparse<Real_wp, Symmetric, RowSymSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_SymSparse<Complex_wp, Symmetric, RowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_SymSparse<Complex_wp, Symmetric, RowSymSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_SymSparse<Complex_wp, Symmetric, RowSymSparse>::Fill(const Complex_wp&);

  SELDON_EXTERN template void Matrix_SymSparse<int, Symmetric, RowSymSparse>::SetData(int, int, Vector<int>&, Vector<long>&, Vector<int>&);
  SELDON_EXTERN template void Matrix_SymSparse<Real_wp, Symmetric, RowSymSparse>::SetData(int, int, Vector<Real_wp>&, Vector<long>&, Vector<int>&);
  SELDON_EXTERN template void Matrix_SymSparse<Complex_wp, Symmetric, RowSymSparse>::SetData(int, int, Vector<Complex_wp>&, Vector<long>&, Vector<int>&);
  
  SELDON_EXTERN template class Matrix_ArraySparse<Real_wp, General, ArrayRowSparse>;
  SELDON_EXTERN template class Matrix_ArraySparse<Complex_wp, General, ArrayRowSparse>;
  SELDON_EXTERN template class Matrix<Real_wp, General, ArrayRowSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, ArrayRowSparse>;

  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, General, ArrayRowSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, General, ArrayRowSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, General, ArrayRowSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, General, ArrayRowSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, General, ArrayRowSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, General, ArrayRowSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, General, ArrayRowSparse>::RemoveSmallEntry(const Real_wp&);


  SELDON_EXTERN template class Matrix_ArraySparse<Real_wp, General, ArrayColSparse>;
  SELDON_EXTERN template class Matrix_ArraySparse<Complex_wp, General, ArrayColSparse>;
  SELDON_EXTERN template class Matrix<Real_wp, General, ArrayColSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, ArrayColSparse>;

  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, General, ArrayColSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, General, ArrayColSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, General, ArrayColSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, General, ArrayColSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, General, ArrayColSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, General, ArrayColSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, General, ArrayColSparse>::RemoveSmallEntry(const Real_wp&);

  SELDON_EXTERN template class Matrix_ArraySparse<Real_wp, Symmetric, ArrayRowSymSparse>;
  SELDON_EXTERN template class Matrix_ArraySparse<Complex_wp, Symmetric, ArrayRowSymSparse>;
  SELDON_EXTERN template class Matrix<Real_wp, Symmetric, ArrayRowSymSparse>;
  SELDON_EXTERN template class Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>;

  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, Symmetric, ArrayRowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, Symmetric, ArrayRowSymSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Real_wp, Symmetric, ArrayRowSymSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, Symmetric, ArrayRowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, Symmetric, ArrayRowSymSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, Symmetric, ArrayRowSymSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_ArraySparse<Complex_wp, Symmetric, ArrayRowSymSparse>::RemoveSmallEntry(const Real_wp&);
  

  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
}
