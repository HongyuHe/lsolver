#ifndef SELDON_WITH_BLAS
#define SELDON_WITH_BLAS
#endif
#ifndef SELDON_WITH_LAPACK
#define SELDON_WITH_LAPACK
#endif

#include "SeldonHeader.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "vector/Vector.cxx"
#include "matrix/Matrix_Base.cxx"
#include "matrix/Matrix_Pointers.cxx"
#include "matrix/Matrix_Symmetric.cxx"
#include "matrix/Matrix_SymPacked.cxx"
#include "matrix/Matrix_Hermitian.cxx"
#include "matrix/Matrix_HermPacked.cxx"
#include "matrix_sparse/Matrix_Sparse.cxx"
#include "share/Allocator.cxx"
#include "share/AllocatorInline.cxx"
#include "share/StorageInline.cxx"
#include "share/MatrixFlagInline.cxx"
#include "computation/interfaces/Blas_1.cxx"
#include "computation/interfaces/Blas_2.cxx"
#include "computation/interfaces/Blas_3.cxx"
#include "computation/basic_functions/Functions_Vector.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#endif

#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;

namespace Seldon
{

  /* Blas Level 1 */

  SELDON_EXTERN template void ApplyRot(Vector<@real>&, Vector<@real>&, const @real&, const @real&);
  SELDON_EXTERN template void ApplyModifRot(Vector<@real>&, Vector<@real>&, const @real*);

  SELDON_EXTERN template void Swap(Vector<@real_complex>&, Vector<@real_complex>&);

  SELDON_EXTERN template void MltScalar(const @scalar&, Vector<@scalar>&);

  SELDON_EXTERN template void CopyVector(const Vector<@real_complex>&, Vector<@real_complex>&);

  SELDON_EXTERN template void AddVector(const @real_complex&, const Vector<@real_complex>&, Vector<@real_complex>&);

  SELDON_EXTERN template @real_complex DotProdVector(const Vector<@real_complex>&, const Vector<@real_complex>&);
  SELDON_EXTERN template @complex DotProdConjVector(const Vector<@complex>&, const Vector<@complex>&);

  SELDON_EXTERN template float Norm1(const Vector<float>&);
  SELDON_EXTERN template double Norm1(const Vector<double>&);
  SELDON_EXTERN template float Norm1(const Vector<complexfloat>&);
  SELDON_EXTERN template double Norm1(const Vector<complexdouble>&);

  SELDON_EXTERN template float Norm2(const Vector<float>&);
  SELDON_EXTERN template double Norm2(const Vector<double>&);
  SELDON_EXTERN template float Norm2(const Vector<complexfloat>&);
  SELDON_EXTERN template double Norm2(const Vector<complexdouble>&);

  SELDON_EXTERN template size_t GetMaxAbsIndex(const Vector<@real_complex>&);

  /* Blas Level 2 */

  SELDON_EXTERN template void Mlt(const Matrix<@real_complex, General, @storage_blasT>&, Vector<@real_complex>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, @storage_blasT>&, Vector<@real_complex>&);
  SELDON_EXTERN template void MltVector(const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, Vector<@real_complex>&);
  SELDON_EXTERN template void MltVector(const Matrix<@real_complex, Symmetric, @storage_blasS>&, const Vector<@real_complex>&, Vector<@real_complex>&);
  SELDON_EXTERN template void MltVector(const Matrix<@complex, Hermitian, @storage_blasH>&, const Vector<@complex>&, Vector<@complex>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, Vector<@real_complex>&);

  SELDON_EXTERN template void MltAddVector(const @real_complex&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, const @real_complex&, Vector<@real_complex>&);
  SELDON_EXTERN template void MltAddVector(const @real_complex&, const SeldonTranspose&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, const @real_complex&, Vector<@real_complex>&);
  SELDON_EXTERN template void MltAddVector(const @real_complex&, const Matrix<@real_complex, Symmetric, @storage_blasS>&, const Vector<@real_complex>&, const @real_complex&, Vector<@real_complex>&);
  SELDON_EXTERN template void MltAddVector(const @complex&, const Matrix<@complex, Hermitian, @storage_blasH>&, const Vector<@complex>&, const @complex&, Vector<@complex>&);

  SELDON_EXTERN template void Rank1Update(const @real_complex&, const Vector<@real_complex>&, const Vector<@real_complex>&, Matrix<@real_complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void Rank1Update(const @complex&, const Vector<@complex>&, const SeldonConjugate&, const Vector<@complex>&, Matrix<@complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void Rank1Update(const @real&, const Vector<@real>&, Matrix<@real, Symmetric, @storage_blasS>&);
  SELDON_EXTERN template void Rank1Update(const float&, const Vector<complexfloat>&, Matrix<complexfloat, Hermitian, @storage_blasH>&);
  SELDON_EXTERN template void Rank1Update(const double&, const Vector<complexdouble>&, Matrix<complexdouble, Hermitian, @storage_blasH>&);

  SELDON_EXTERN template void Rank2Update(const @real&, const Vector<@real>&, const Vector<@real>&, Matrix<@real, Symmetric, @storage_blasS>&);
  SELDON_EXTERN template void Rank2Update(const @complex&, const Vector<@complex>&, const Vector<@complex>&, Matrix<@complex, Hermitian, @storage_blasH>&);

  SELDON_EXTERN template void Solve(const Matrix<@real_complex, General, @storage_blasT>&, Vector<@real_complex>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, @storage_blasT>&, Vector<@real_complex>&);

  /* Blas Level 3 */

  SELDON_EXTERN template void MltAddMatrix(const @real_complex&, const Matrix<@real_complex, General, @storage_blasGE>&, const Matrix<@real_complex, General, @storage_blasGE>&, const @real_complex&, Matrix<@real_complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void MltAddMatrix(const @real_complex&, const SeldonTranspose&, const Matrix<@real_complex, General, @storage_blasGE>&, const SeldonTranspose&, const Matrix<@real_complex, General, @storage_blasGE>&, const @real_complex&, Matrix<@real_complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void MltAdd(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, Symmetric, ColSym>&, const Matrix<@real_complex, General, ColMajor>&, const @real_complex&, Matrix<@real_complex, General, ColMajor>&);
  SELDON_EXTERN template void MltAdd(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, Symmetric, RowSym>&, const Matrix<@real_complex, General, RowMajor>&, const @real_complex&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void MltAdd(const SeldonSide&, const @complex&, const Matrix<@complex, Symmetric, ColHerm>&, const Matrix<@complex, General, ColMajor>&, const @complex&, Matrix<@complex, General, ColMajor>&);
  SELDON_EXTERN template void MltAdd(const SeldonSide&, const @complex&, const Matrix<@complex, Symmetric, RowHerm>&, const Matrix<@complex, General, RowMajor>&, const @complex&, Matrix<@complex, General, RowMajor>&);

  //SELDON_EXTERN template void MltMatrix(const Matrix<@real_complex, General, @storage_blasGE>&, const Matrix<@real_complex, General, @storage_blasGE>&, Matrix<@real_complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, General, RowLoTriang>&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, General, RowUpTriang>&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, General, ColLoTriang>&, Matrix<@real_complex, General, ColMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, General, ColUpTriang>&, Matrix<@real_complex, General, ColMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const @real_complex&, const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, RowLoTriang>&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const @real_complex&, const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, RowUpTriang>&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const @real_complex&, const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, ColLoTriang>&, Matrix<@real_complex, General, ColMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const @real_complex&, const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, ColUpTriang>&, Matrix<@real_complex, General, ColMajor>&);

  SELDON_EXTERN template void Solve(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, General, RowLoTriang>&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, General, RowUpTriang>&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, General, ColLoTriang>&, Matrix<@real_complex, General, ColMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const @real_complex&, const Matrix<@real_complex, General, ColUpTriang>&, Matrix<@real_complex, General, ColMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const @real_complex&, const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, RowLoTriang>&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const @real_complex&, const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, RowUpTriang>&, Matrix<@real_complex, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const @real_complex&, const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, ColLoTriang>&, Matrix<@real_complex, General, ColMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const @real_complex&, const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, ColUpTriang>&, Matrix<@real_complex, General, ColMajor>&);

  /* Other functions */

  SELDON_EXTERN template void Conjugate(Vector<@complex>&);
  SELDON_EXTERN template void Conjugate(Matrix<@complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void Conjugate(Matrix<@complex, Symmetric, @storage_blasS>&);
  SELDON_EXTERN template void Conjugate(Matrix<@complex, Hermitian, @storage_blasH>&);

  SELDON_EXTERN template void Transpose(Matrix<@real_complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void Transpose(Matrix<@real_complex, Symmetric, @storage_blasS>&);
  SELDON_EXTERN template void Transpose(Matrix<@complex, Hermitian, @storage_blasH>&);

  SELDON_EXTERN template void TransposeConj(Matrix<@real_complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void TransposeConj(Matrix<@real_complex, Symmetric, @storage_blasS>&);
  SELDON_EXTERN template void TransposeConj(Matrix<@complex, Hermitian, @storage_blasH>&);

  SELDON_EXTERN template void MltAddVector(const @real_complex&, const Matrix<@real_complex, General, RowSparse>&, const Vector<@real_complex>&, const @real_complex&, Vector<@real_complex>&);
  SELDON_EXTERN template void MltAddVector(const @real_complex&, const SeldonTranspose&, const Matrix<@real_complex, General, RowSparse>&, const Vector<@real_complex>&, const @real_complex&, Vector<@real_complex>&);

}
