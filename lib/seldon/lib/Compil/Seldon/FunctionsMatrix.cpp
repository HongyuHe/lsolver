#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/basic_functions/Functions_Matrix.cxx"
#ifdef SELDON_WITH_MKL
#include "computation/interfaces/Mkl_Sparse.cxx"
#endif
#endif

namespace Seldon
{
  /* Dense matrices */
  
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, General, ColMajor>&); 
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, General, RowMajor>&); 
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, Symmetric, RowSym>&);
  
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, General, ColMajor>&); 
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, General, RowMajor>&);   
  SELDON_EXTERN template void MltScalar(const Complex_wp &, Matrix<Complex_wp , General, RowMajor>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp &, Matrix<Complex_wp , General, ColMajor>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp &, Matrix<Complex_wp , Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Complex_wp , Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp &, Matrix<Complex_wp , Symmetric, RowSym>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Complex_wp , Symmetric, RowSym>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Complex_wp , Hermitian, RowHerm>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Complex_wp , Hermitian, RowHermPacked>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Real_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp &, Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Real_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp &, Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Real_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp &, Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Real_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template void MltScalar(const Real_wp &, Matrix<Complex_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp &, Matrix<Complex_wp, General, RowUpTriangPacked>&);
  
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, Matrix<Real_wp, General, RowMajor> const&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowMajor>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowMajor>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, ColMajor>&, Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSym>&, Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSym>&, Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, Hermitian, RowHerm>&, Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowLoTriang>&, Matrix<Real_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, RowLoTriang>&, Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriang>&, Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowLoTriangPacked>&, Matrix<Real_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowUpTriang>&, Matrix<Real_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, RowUpTriang>&, Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriang>&, Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowUpTriangPacked>&, Matrix<Real_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, Matrix<Complex_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, Matrix<Complex_wp, General, RowUpTriangPacked>&);

  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, RowLoTriang>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, RowUpTriang>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, RowUpTriangPacked>&);

  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, RowLoTriang>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, RowUpTriang>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, RowUpTriangPacked>&);

  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, RowLoTriang>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, RowLoTriang>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, RowLoTriangPacked>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, RowUpTriang>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, RowUpTriang>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, RowUpTriangPacked>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, RowUpTriangPacked>&);

  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, General, RowMajor>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, General, RowMajor>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, General, ColMajor>&, Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, Symmetric, RowSym>&, Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, Symmetric, RowSym>&, Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, Hermitian, RowHerm>&, Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, Matrix<Complex_wp, Hermitian, RowHermPacked>&);

  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, General, RowMajor>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, General, RowMajor>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, General, ColMajor>&, Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, Symmetric, RowSym>&, Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, Symmetric, RowSym>&, Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, Hermitian, RowHerm>&, Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, ColMajor>&);

  SELDON_EXTERN template void MltAddMatrix(Real_wp const&, SeldonTranspose const&, Matrix<Real_wp, General, ColMajor> const&, SeldonTranspose const&, Matrix<Real_wp, General, ColMajor> const&, Real_wp const&, Matrix<Real_wp, Symmetric, RowSymPacked>&);


  SELDON_EXTERN template void MltAddMatrix(Real_wp const&, SeldonTranspose const&, Matrix<Real_wp, General, RowMajor> const&, SeldonTranspose const&, Matrix<Real_wp, General, RowMajor> const&, Real_wp const&, Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void MltAddMatrix(Real_wp const&, SeldonTranspose const&, Matrix<Real_wp, General, RowMajor> const&, SeldonTranspose const&, Matrix<Real_wp, General, RowMajor> const&, Real_wp const&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  
  SELDON_EXTERN template void MltAddMatrix(Complex_wp const&, SeldonTranspose const&, Matrix<Complex_wp, General, RowMajor> const&, SeldonTranspose const&, Matrix<Complex_wp, General, RowMajor> const&, Complex_wp const&, Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void MltAddMatrix(Complex_wp const&, SeldonTranspose const&, Matrix<Complex_wp, General, RowMajor> const&, SeldonTranspose const&, Matrix<Complex_wp, General, RowMajor> const&, Complex_wp const&, Matrix<Complex_wp, Symmetric, RowSym>&);
  
  SELDON_EXTERN template void MltAddMatrix(Real_wp const&, Matrix<Real_wp, Symmetric, RowSymPacked> const&, Matrix<Real_wp, General, RowMajor> const&, Real_wp const&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void MltAddMatrix(Complex_wp const&, Matrix<Complex_wp, Symmetric, RowSymPacked> const&, Matrix<Complex_wp, General, RowMajor> const&, Complex_wp const&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void MltAddMatrix(Real_wp const&, Matrix<Real_wp, General, RowMajor> const&, Matrix<Real_wp, Symmetric, RowSymPacked> const&, Real_wp const&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void MltAddMatrix(Complex_wp const&, Matrix<Complex_wp, General, RowMajor> const&, Matrix<Complex_wp, Symmetric, RowSymPacked> const&, Complex_wp const&, Matrix<Complex_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowLoTriangPacked>&,
					   const Matrix<Real_wp, General, RowLoTriangPacked>&, const Real_wp&,
					   Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowLoTriang>&,
					   const Matrix<Real_wp, General, RowLoTriang>&, const Real_wp&,
					   Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowUpTriangPacked>&,
					   const Matrix<Real_wp, General, RowUpTriangPacked>&, const Real_wp&,
					   Matrix<Real_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowUpTriang>&,
					   const Matrix<Real_wp, General, RowUpTriang>&, const Real_wp&,
					   Matrix<Real_wp, General, RowMajor>&);


  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriangPacked>&,
					   const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Complex_wp&,
					   Matrix<Complex_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriang>&,
					   const Matrix<Complex_wp, General, RowLoTriang>&, const Complex_wp&,
					   Matrix<Complex_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriangPacked>&,
					   const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Complex_wp&,
					   Matrix<Complex_wp, General, RowMajor>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriang>&,
					   const Matrix<Complex_wp, General, RowUpTriang>&, const Complex_wp&,
					   Matrix<Complex_wp, General, RowMajor>&);


  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowLoTriangPacked>&,
					   const Matrix<Real_wp, General, RowLoTriangPacked>&, const Real_wp&,
					   Matrix<Real_wp, General, RowLoTriangPacked>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowLoTriang>&,
					   const Matrix<Real_wp, General, RowLoTriang>&, const Real_wp&,
					   Matrix<Real_wp, General, RowLoTriang>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowUpTriangPacked>&,
					   const Matrix<Real_wp, General, RowUpTriangPacked>&, const Real_wp&,
					   Matrix<Real_wp, General, RowUpTriangPacked>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowUpTriang>&,
					   const Matrix<Real_wp, General, RowUpTriang>&, const Real_wp&,
					   Matrix<Real_wp, General, RowUpTriang>&);


  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriangPacked>&,
					   const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Complex_wp&,
					   Matrix<Complex_wp, General, RowLoTriangPacked>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriang>&,
					   const Matrix<Complex_wp, General, RowLoTriang>&, const Complex_wp&,
					   Matrix<Complex_wp, General, RowLoTriang>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriangPacked>&,
					   const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Complex_wp&,
					   Matrix<Complex_wp, General, RowUpTriangPacked>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriang>&,
					   const Matrix<Complex_wp, General, RowUpTriang>&, const Complex_wp&,
					   Matrix<Complex_wp, General, RowUpTriang>&);

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Real_wp&, Matrix<Real_wp, Symmetric, RowSymPacked>&);

  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Complex_wp&, Matrix<Complex_wp, Symmetric, RowSymPacked>&);


  /* Sparse matrices */
  
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowSparse>&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, ColSparse>&, Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, RowSparse>&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowSparse>&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, ColSparse>&, Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, ColSparse>&, Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template Real_wp MaxAbs(const Matrix<Complex_wp, Symmetric, RowSymSparse>&);
    
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void Conjugate(Matrix<Complex_wp, General, ColSparse>&);
    
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, General, RowSparse>&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, General, RowSparse>&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void Transpose(Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Real_wp, General, ColSparse>&, Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void Transpose(const Matrix<Complex_wp, General, ColSparse>&, Matrix<Complex_wp, General, ColSparse>&);
    
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, General, RowSparse>&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, General, RowSparse>&, Matrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Matrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void TransposeConj(Matrix<Complex_wp, General, ColSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Real_wp, General, ColSparse>&, Matrix<Real_wp, General, ColSparse>&);
  SELDON_EXTERN template void TransposeConj(const Matrix<Complex_wp, General, ColSparse>&, Matrix<Complex_wp, General, ColSparse>&);
  
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowSparse>&, const Matrix<Real_wp, General, RowMajor>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void MltMatrix(const Matrix<Real_wp, General, RowSparse>&, const Matrix<Real_wp, General, RowSparse>&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowSparse>&, const Matrix<Real_wp, General, RowSparse>&, const Real_wp&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const Real_wp&, Matrix<Real_wp, General, RowSparse>&);
   
}
