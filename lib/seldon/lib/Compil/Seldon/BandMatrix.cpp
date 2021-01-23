#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "matrix_sparse/BandMatrix.hxx"
#include "matrix_sparse/BandMatrixInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/BandMatrix.cxx"
#include "computation/basic_functions/Functions_Base.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"
#endif

namespace Seldon
{
  /* BandedCol */

  SELDON_EXTERN template class Matrix_Band<Real_wp, General, BandedCol>;
  SELDON_EXTERN template class Matrix_Band<Complex_wp, General, BandedCol>;
  
  SELDON_EXTERN template class Matrix<Real_wp, General, BandedCol>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, BandedCol>;

  SELDON_EXTERN template void Matrix_Band<Real_wp, General, BandedCol>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Band<Real_wp, General, BandedCol>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Band<Real_wp, General, BandedCol>::Solve(Vector<Real_wp>&) const;
  SELDON_EXTERN template void Matrix_Band<Real_wp, General, BandedCol>::Solve(Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, BandedCol>::Solve(Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Matrix_Band<Real_wp, General, BandedCol>::Solve(const Vector<int>&, Vector<Real_wp>&) const;
  SELDON_EXTERN template void Matrix_Band<Real_wp, General, BandedCol>::Solve(const Vector<int>&, Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, BandedCol>::Solve(const Vector<int>&, Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Matrix_Band<Real_wp, General, BandedCol>::Add_(const Real_wp&, const Matrix<Real_wp, General, BandedCol>&);


  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, BandedCol>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, BandedCol>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, BandedCol>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, BandedCol>::Add_(const Complex_wp&, const Matrix<Complex_wp, General, BandedCol>&);

  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, BandedCol>&, Matrix<Real_wp, General, BandedCol>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, BandedCol>&, Matrix<Complex_wp, General, BandedCol>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, BandedCol>&);
  
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, BandedCol>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, BandedCol>&, const Vector<int>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, BandedCol>&, const Vector<int>&, Vector<Complex_wp>&, LapackInfo&);

  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, BandedCol>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const Matrix<Complex_wp, General, BandedCol>&, const Vector<int>&, Vector<Complex_wp>&, LapackInfo&);

  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void Copy(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, BandedCol>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, BandedCol>&, Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, BandedCol>&, Matrix<Complex_wp, General, BandedCol>&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, BandedCol>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, BandedCol>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, BandedCol>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, BandedCol>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, BandedCol>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
    
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, BandedCol>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, BandedCol>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, BandedCol>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, BandedCol>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, General, BandedCol>&);

  SELDON_EXTERN template void ScaleMatrix(Matrix<Real_wp, General, BandedCol>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(Matrix<Complex_wp, General, BandedCol>&, const Vector<Real_wp>&, const Vector<Real_wp>&);

  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, BandedCol>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, BandedCol>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, BandedCol>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, BandedCol>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, BandedCol>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, BandedCol>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  
  /* ArrowCol */  

  SELDON_EXTERN template class Matrix_Band<Real_wp, General, ArrowCol>;
  SELDON_EXTERN template class Matrix_Band<Complex_wp, General, ArrowCol>;

  SELDON_EXTERN template class Matrix_Arrow<Real_wp, General, ArrowCol>;
  SELDON_EXTERN template class Matrix_Arrow<Complex_wp, General, ArrowCol>;

  SELDON_EXTERN template class Matrix<Real_wp, General, ArrowCol>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, ArrowCol>;

  SELDON_EXTERN template void Matrix_Band<Real_wp, General, ArrowCol>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Band<Real_wp, General, ArrowCol>::Fill(const Real_wp&);

  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, ArrowCol>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, ArrowCol>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Band<Complex_wp, General, ArrowCol>::Fill(const Complex_wp&);

  SELDON_EXTERN template void Matrix_Arrow<Real_wp, General, ArrowCol>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Arrow<Real_wp, General, ArrowCol>::Fill(const Real_wp&);

  SELDON_EXTERN template void Matrix_Arrow<Complex_wp, General, ArrowCol>::Fill(const int&);
  SELDON_EXTERN template void Matrix_Arrow<Complex_wp, General, ArrowCol>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_Arrow<Complex_wp, General, ArrowCol>::Fill(const Complex_wp&);
  SELDON_EXTERN template void Matrix_Arrow<Real_wp, General, ArrowCol>::Solve(Vector<Real_wp>&) const;
  SELDON_EXTERN template void Matrix_Arrow<Real_wp, General, ArrowCol>::Solve(Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Matrix_Arrow<Complex_wp, General, ArrowCol>::Solve(Vector<Complex_wp>&) const;
  
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, ArrowCol>&, Matrix<Real_wp, General, ArrowCol>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, ArrowCol>&, Matrix<Complex_wp, General, ArrowCol>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, ArrowCol>&);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, ArrowCol>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, ArrowCol>&, Matrix<Real_wp, General, ArrowCol>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, ArrowCol>&, Matrix<Complex_wp, General, ArrowCol>&);

  SELDON_EXTERN template void Matrix_Arrow<Real_wp, General, ArrowCol>::MltAdd(const Real_wp&, const SeldonTranspose&, const Vector<Real_wp>&, Vector<Real_wp>&) const;
  SELDON_EXTERN template void Matrix_Arrow<Real_wp, General, ArrowCol>::MltAdd(const Complex_wp&, const SeldonTranspose&, const Vector<Complex_wp>&, Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Matrix_Arrow<Complex_wp, General, ArrowCol>::MltAdd(const Complex_wp&, const SeldonTranspose&, const Vector<Complex_wp>&, Vector<Complex_wp>&) const;

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ArrowCol>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, ArrowCol>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ArrowCol>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, ArrowCol>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ArrowCol>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrowCol>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, ArrowCol>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrowCol>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
    
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, General, ArrowCol>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, General, ArrowCol>&);

  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, ArrowCol>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrowCol>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);

  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, ArrowCol>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrowCol>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);

}
