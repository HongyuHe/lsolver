#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "SeldonComplexMatrixHeader.hxx"
#include "SeldonComplexMatrixInline.hxx"

#ifdef SELDON_WITH_MPI
#include "SeldonDistributedHeader.hxx"
#include "SeldonDistributedInline.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/DistributedMatrix.cxx"
#include "matrix_sparse/DistributedMatrixFunction.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template class DistributedMatrix_Base<Real_wp>;
  SELDON_EXTERN template class DistributedMatrix_Base<Complex_wp>;

  SELDON_EXTERN template void DistributedMatrix_Base<Real_wp>::AssembleVec(Vector<Real_wp>&) const;
  
  SELDON_EXTERN template void DistributedMatrix_Base<Real_wp>::Init(const DistributedMatrix_Base<Real_wp>&);
  SELDON_EXTERN template void DistributedMatrix_Base<Complex_wp>::Init(const DistributedMatrix_Base<Real_wp>&);
  SELDON_EXTERN template void DistributedMatrix_Base<Real_wp>::Init(const DistributedMatrix_Base<Complex_wp>&);
  SELDON_EXTERN template void DistributedMatrix_Base<Complex_wp>::Init(const DistributedMatrix_Base<Complex_wp>&);

  /* ArrayRowSparse */

  SELDON_EXTERN template class DistributedMatrix<Real_wp, General, ArrayRowSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Real_wp, General, ArrayRowSparse>& DistributedMatrix<Real_wp, General, ArrayRowSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, ArrayRowSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, ArrayRowSparse>::Fill(const int&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, ArrayRowSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, ArrayRowSparse>::GetDistributedRows(Matrix<Real_wp, General, ArrayRowSparse>&, Vector<IVect>& ) const;
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, ArrayRowSparse>::GetDistributedColumns(Matrix<Real_wp, General, ArrayColSparse>&, Vector<IVect>&, bool) const;
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Real_wp&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MltMin(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, IVect&, IVect&);
  SELDON_EXTERN template void MltMin(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, IVect&, IVect&);
  SELDON_EXTERN template void MltMin(const Matrix<Real_wp, General, ArrayRowSparse>&,
				     const IVect&, IVect&, IVect&);
  SELDON_EXTERN template void MltMin(const Matrix<Complex_wp, General, ArrayRowSparse>&,
				     const IVect&, IVect&, IVect&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void GetRow(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void ScaleLeftMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool, bool);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void DistributeCentralizedMatrix(Matrix<Real_wp, General, ArrayRowSparse>&,
							  const MPI_Comm& comm, int, const Vector<int>&, DistributedMatrixIntegerArray&,
							  DistributedMatrix<Real_wp, General, ArrayRowSparse>&, bool);

  SELDON_EXTERN template void DistributeCentralizedMatrix(Matrix<Complex_wp, General, ArrayRowSparse>&,
							  const MPI_Comm& comm, int, const Vector<int>&, DistributedMatrixIntegerArray&,
							  DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, bool);

  /* RowSparse */

  SELDON_EXTERN template class DistributedMatrix<Real_wp, General, RowSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Real_wp, General, RowSparse>& DistributedMatrix<Real_wp, General, RowSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, RowSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, RowSparse>::Fill(const int&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, RowSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, RowSparse>::GetDistributedRows(Matrix<Real_wp, General, ArrayRowSparse>&, Vector<IVect>&) const;
  SELDON_EXTERN template void DistributedMatrix<Real_wp, General, RowSparse>::GetDistributedColumns(Matrix<Real_wp, General, ArrayColSparse>&, Vector<IVect>&, bool) const;
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Real_wp&, const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Real_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Real_wp, General, RowSparse>&, DistributedMatrix<Real_wp, General, RowSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Real_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Real_wp, General, RowSparse>&, DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Real_wp, General, RowSparse>&, DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Real_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void GetRow(const DistributedMatrix<Real_wp, General, RowSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const DistributedMatrix<Real_wp, General, RowSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, DistributedMatrix<Real_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Real_wp, General, RowSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Real_wp, General, RowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, General, RowSparse>&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void ScaleLeftMatrix(DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, General, RowSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, General, RowSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, General, RowSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, General, RowSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool, bool);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Real_wp, General, RowSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void DistributeCentralizedMatrix(Matrix<Real_wp, General, RowSparse>&,
							  const MPI_Comm& comm, int, const Vector<int>&, DistributedMatrixIntegerArray&,
							  DistributedMatrix<Real_wp, General, RowSparse>&, bool);

  SELDON_EXTERN template void DistributeCentralizedMatrix(Matrix<Complex_wp, General, RowSparse>&,
							  const MPI_Comm& comm, int, const Vector<int>&, DistributedMatrixIntegerArray&,
							  DistributedMatrix<Complex_wp, General, RowSparse>&, bool);
  
  /* ArrayRowSymSparse */

  SELDON_EXTERN template class DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>::GetDistributedRows(Matrix<Real_wp, General, ArrayRowSparse>&, Vector<IVect>& ) const;
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>::GetDistributedColumns(Matrix<Real_wp, General, ArrayColSparse>&, Vector<IVect>&, bool) const;
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Real_wp&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void MltMin(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, IVect&, IVect&);
  SELDON_EXTERN template void MltMin(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, IVect&, IVect&);
  SELDON_EXTERN template void MltMin(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
				     const IVect&, IVect&, IVect&);
  SELDON_EXTERN template void MltMin(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
				     const IVect&, IVect&, IVect&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void GetRow(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool, bool);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void DistributeCentralizedMatrix(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
							  const MPI_Comm& comm, int, const Vector<int>&, DistributedMatrixIntegerArray&,
							  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);

  SELDON_EXTERN template void DistributeCentralizedMatrix(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
							  const MPI_Comm& comm, int, const Vector<int>&, DistributedMatrixIntegerArray&,
							  DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool);

  /* RowSymSparse */

  SELDON_EXTERN template class DistributedMatrix<Real_wp, Symmetric, RowSymSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Real_wp, Symmetric, RowSymSparse>& DistributedMatrix<Real_wp, Symmetric, RowSymSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, RowSymSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, RowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, RowSymSparse>::Fill(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, RowSymSparse>::GetDistributedRows(Matrix<Real_wp, General, ArrayRowSparse>&, Vector<IVect>&) const;
  SELDON_EXTERN template void DistributedMatrix<Real_wp, Symmetric, RowSymSparse>::GetDistributedColumns(Matrix<Real_wp, General, ArrayColSparse>&, Vector<IVect>&, bool) const;
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Real_wp&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void GetRow(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, int, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Real_wp, VectSparse>&, int, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Real_wp, VectSparse>&, int, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Real_wp>&, const Vector<Real_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Real_wp>&, bool, bool);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void DistributeCentralizedMatrix(Matrix<Real_wp, Symmetric, RowSymSparse>&,
							  const MPI_Comm& comm, int, const Vector<int>&, DistributedMatrixIntegerArray&,
							  DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, bool);

  SELDON_EXTERN template void DistributeCentralizedMatrix(Matrix<Complex_wp, Symmetric, RowSymSparse>&,
							  const MPI_Comm& comm, int, const Vector<int>&, DistributedMatrixIntegerArray&,
							  DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, bool);

  /* complex ArrayRowSparse */

  SELDON_EXTERN template class DistributedMatrix<Complex_wp, General, ArrayRowSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Complex_wp, General, ArrayRowSparse>& DistributedMatrix<Complex_wp, General, ArrayRowSparse>::operator *=(const Complex_wp&);
  SELDON_EXTERN template DistributedMatrix<Complex_wp, General, ArrayRowSparse>& DistributedMatrix<Complex_wp, General, ArrayRowSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, ArrayRowSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, ArrayRowSparse>::Fill(const int&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, ArrayRowSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, ArrayRowSparse>::GetDistributedRows(Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<IVect>& ) const;
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, ArrayRowSparse>::GetDistributedColumns(Matrix<Complex_wp, General, ArrayColSparse>&, Vector<IVect>&, bool) const;
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Complex_wp&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void GetRow(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);

  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void ScaleLeftMatrix(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Complex_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool, bool);

  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);


  /* complex RowSparse */

  SELDON_EXTERN template class DistributedMatrix<Complex_wp, General, RowSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Complex_wp, General, RowSparse>& DistributedMatrix<Complex_wp, General, RowSparse>::operator *=(const Complex_wp&);
  SELDON_EXTERN template DistributedMatrix<Complex_wp, General, RowSparse>& DistributedMatrix<Complex_wp, General, RowSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, RowSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, RowSparse>::Fill(const int&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, RowSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, RowSparse>::GetDistributedRows(Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<IVect>&) const;
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, RowSparse>::GetDistributedColumns(Matrix<Complex_wp, General, ArrayColSparse>&, Vector<IVect>&, bool) const;
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Complex_wp&, const DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Complex_wp, General, RowSparse>&, DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, General, RowSparse>&, DistributedMatrix<Complex_wp, General, RowSparse>&);  
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Complex_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Complex_wp, General, RowSparse>&, DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Complex_wp, General, RowSparse>&, DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Complex_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void GetRow(const DistributedMatrix<Complex_wp, General, RowSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const DistributedMatrix<Complex_wp, General, RowSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Complex_wp, General, RowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, RowSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, RowSparse>&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void ScaleLeftMatrix(DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Complex_wp, General, RowSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, General, RowSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, General, RowSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, General, RowSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, General, RowSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool, bool);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, General, RowSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  /* complex ArrayRowSymSparse */

  SELDON_EXTERN template class DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>& DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>::operator *=(const Complex_wp&);
  SELDON_EXTERN template DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>& DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>::GetDistributedRows(Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<IVect>& ) const;
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>::GetDistributedColumns(Matrix<Complex_wp, General, ArrayColSparse>&, Vector<IVect>&, bool) const;
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void GetRow(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);

  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool, bool);

  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  /* complex RowSymSparse */

  SELDON_EXTERN template class DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>& DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>::operator *=(const Complex_wp&);
  SELDON_EXTERN template DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>& DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>::Fill(const int&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>::GetDistributedRows(Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<IVect>&) const;
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>::GetDistributedColumns(Matrix<Complex_wp, General, ArrayColSparse>&, Vector<IVect>&, bool) const;
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void GetRow(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void GetCol(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, int, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  
  SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<int>&, const Vector<int>&);
  SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void SorVector(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, const Real_wp&, int, int);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<long>&, Vector<int>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, Symmetric&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void AssembleDistributed(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, General&, const MPI_Comm&, IVect&, IVect&, Vector<int64_t>&, Vector<int64_t>&, Vector<Complex_wp>&, bool, bool);

  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);
  
  /* ArrayRowComplexSparse */

  SELDON_EXTERN template class DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>& DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>::RemoveSmallEntry(const Real_wp&);
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Complex_wp&, const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void MltMin(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&,
				     const IVect&, IVect&, IVect&);
  SELDON_EXTERN template void MltMin(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, IVect&, IVect&);
  
  //SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  
  //SELDON_EXTERN template void GetRow(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, int, Vector<Complex_wp, VectSparse>&);
  //SELDON_EXTERN template void GetCol(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, int, Vector<Complex_wp, VectSparse>&);
  //SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  //SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  
  //SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<int>&, const Vector<int>&);
  //SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  
  SELDON_EXTERN template void ScaleLeftMatrix(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);
  
  /* ArrayRowSymComplexSparse */

  SELDON_EXTERN template class DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>& DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>::operator *=(const Real_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>::Fill(const Complex_wp&);
  SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>::RemoveSmallEntry(const Real_wp&);
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void MltMin(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&,
				     const IVect&, IVect&, IVect&);
  SELDON_EXTERN template void MltMin(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, IVect&, IVect&);
  
  //SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  //SELDON_EXTERN template void GetRow(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, int, Vector<Complex_wp, VectSparse>&);
  //SELDON_EXTERN template void GetCol(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, int, Vector<Complex_wp, VectSparse>&);
  //SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  //SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  //SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<int>&, const Vector<int>&);
  //SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, DistributedMatrix<Complex_wp, General, RowSparse>&);
  
  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  /* RowComplexSparse */

  SELDON_EXTERN template class DistributedMatrix<Complex_wp, General, RowComplexSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Complex_wp, General, RowComplexSparse>& DistributedMatrix<Complex_wp, General, RowComplexSparse>::operator *=(const Real_wp&);
  //SELDON_EXTERN template void DistributedMatrix<Complex_wp, General, RowComplexSparse>::RemoveSmallEntry(const Real_wp&);
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Complex_wp&, const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  
  //SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  
  //SELDON_EXTERN template void GetRow(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, int, Vector<Complex_wp, VectSparse>&);
  //SELDON_EXTERN template void GetCol(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, int, Vector<Complex_wp, VectSparse>&);
  //SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  //SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  
  //SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<int>&, const Vector<int>&);
  //SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  
  SELDON_EXTERN template void ScaleLeftMatrix(DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleRightMatrix(DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, General, RowComplexSparse>&);

  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, General, RowComplexSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);

  /* RowSymComplexSparse */

  SELDON_EXTERN template class DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>;
  
  // methods of distributed matrices
  SELDON_EXTERN template DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>& DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>::operator *=(const Real_wp&);
  //SELDON_EXTERN template void DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>::RemoveSmallEntry(const Real_wp&);
  
  // functions for distributed matrices
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  
  //SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  
  SELDON_EXTERN template Real_wp MaxAbs(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template Real_wp NormInf(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  
  SELDON_EXTERN template void Transpose(DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void Conjugate(DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void Transpose(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void TransposeConj(DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  
  //SELDON_EXTERN template void GetRow(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, int, Vector<Complex_wp, VectSparse>&);
  //SELDON_EXTERN template void GetCol(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, int, Vector<Complex_wp, VectSparse>&);
  //SELDON_EXTERN template void SetRow(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  //SELDON_EXTERN template void SetCol(const Vector<Complex_wp, VectSparse>&, int, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  
  //SELDON_EXTERN template void ApplyPermutation(DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<int>&, const Vector<int>&);
  //SELDON_EXTERN template void ApplyInversePermutation(DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<int>&, const Vector<int>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);

  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void ScaleMatrix(DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const Vector<Real_wp>&, const Vector<Real_wp>&);
  
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  
  SELDON_EXTERN template void CopySubMatrix(const DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&,
					    const IVect&, const IVect&, DistributedMatrix<Complex_wp, General, RowSparse>&);

}
