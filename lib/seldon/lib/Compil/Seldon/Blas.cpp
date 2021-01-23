#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/Blas_1.cxx"
#include "computation/interfaces/Blas_2.cxx"
#include "computation/interfaces/Blas_3.cxx"
#include "computation/basic_functions/Functions_Vector.cxx"
#include "computation/basic_functions/Functions_Base.cxx"
#ifdef SELDON_WITH_MKL
#include "computation/interfaces/Mkl_Sparse.cxx"
#endif
#endif

namespace Seldon
{
  /*  Blas Level 1 */

  SELDON_EXTERN template void GenRot(Complex_wp&, Complex_wp&, Real_wp&, Complex_wp&);
  SELDON_EXTERN template void ApplyRot(Real_wp&, Real_wp&, const Real_wp&, const Real_wp&);
  SELDON_EXTERN template void ApplyRot(Complex_wp&, Complex_wp&, const Real_wp&, const Complex_wp&);

  SELDON_EXTERN template void ApplyRot(Vector<Real_wp>&, Vector<Real_wp>&, const Real_wp&, const Real_wp&);
  SELDON_EXTERN template void ApplyModifRot(Vector<Real_wp>&, Vector<Real_wp>&, const Real_wp*);
  
  SELDON_EXTERN template void Swap(Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Swap(Vector<Complex_wp>&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void MltScalar(const int&, Vector<int>&);
  SELDON_EXTERN template void MltScalar(const complex<int>&, Vector<complex<int> >&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void CopyVector(const Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void CopyVector(const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void CopyVector(const Vector<Complex_wp>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void AddVector(const Real_wp&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void AddVector(const Complex_wp&, const Vector<Complex_wp>&, Vector<Complex_wp>&);

  SELDON_EXTERN template Real_wp DotProdVector(const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template Complex_wp DotProdVector(const Vector<Complex_wp>&, const Vector<Complex_wp>&);

  SELDON_EXTERN template Real_wp DotProdConjVector(const Vector<Real_wp>&, const Vector<Real_wp>&);
  SELDON_EXTERN template Complex_wp DotProdConjVector(const Vector<Complex_wp>&, const Vector<Complex_wp>&);

  SELDON_EXTERN template Real_wp Norm1(const Vector<Real_wp>&);
  SELDON_EXTERN template Real_wp Norm1(const Vector<Complex_wp>&);

  SELDON_EXTERN template Real_wp Norm2(const Vector<Real_wp>&);
  SELDON_EXTERN template Real_wp Norm2(const Vector<Complex_wp>&);

  SELDON_EXTERN template size_t GetMaxAbsIndex(const Vector<Real_wp>&);
  SELDON_EXTERN template size_t GetMaxAbsIndex(const Vector<Complex_wp>&);

  /* Blas Level 2 */
  
  SELDON_EXTERN template void Mlt(const Matrix<Real_wp, General, RowUpTriang>&, Vector<Real_wp>&); 
  SELDON_EXTERN template void Mlt(const Matrix<Complex_wp, General, RowUpTriang>&, Vector<Complex_wp>&); 
  SELDON_EXTERN template void Mlt(const Matrix<Real_wp, General, RowLoTriang>&, Vector<Real_wp>&); 
  SELDON_EXTERN template void Mlt(const Matrix<Complex_wp, General, RowLoTriang>&, Vector<Complex_wp>&); 
  SELDON_EXTERN template void Mlt(const Matrix<Real_wp, General, RowUpTriangPacked>&, Vector<Real_wp>&); 
  SELDON_EXTERN template void Mlt(const Matrix<Complex_wp, General, RowUpTriangPacked>&, Vector<Complex_wp>&); 
  SELDON_EXTERN template void Mlt(const Matrix<Real_wp, General, RowLoTriangPacked>&, Vector<Real_wp>&); 
  SELDON_EXTERN template void Mlt(const Matrix<Complex_wp, General, RowLoTriangPacked>&, Vector<Complex_wp>&); 
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowUpTriang>&, Vector<Real_wp>&); 
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowUpTriang>&, Vector<Complex_wp>&); 
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowLoTriang>&, Vector<Real_wp>&); 
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowLoTriang>&, Vector<Complex_wp>&); 
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowUpTriangPacked>&, Vector<Real_wp>&); 
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, Vector<Complex_wp>&); 
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowLoTriangPacked>&, Vector<Real_wp>&); 
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, Vector<Complex_wp>&); 
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  
  SELDON_EXTERN template void Rank1Update(const Real_wp&, const Vector<Real_wp>&, const Vector<Real_wp>&, Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void Rank1Update(const Complex_wp&, const Vector<Complex_wp>&, const Vector<Complex_wp>&, Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void Rank1Update(const Real_wp&, const Vector<Real_wp>&, const Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Rank1Update(const Complex_wp&, const Vector<Complex_wp>&, const Vector<Complex_wp>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Rank1Update(const Complex_wp&, const Vector<Complex_wp>&, const SeldonConjugate&, const Vector<Complex_wp>&, Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void Rank1Update(const Complex_wp&, const Vector<Complex_wp>&, const SeldonConjugate&, const Vector<Complex_wp>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Rank1Update(const Real_wp&, const Vector<Real_wp>&, Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void Rank1Update(const Real_wp&, const Vector<Real_wp>&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void Rank1Update(const Real_wp&, const Vector<Complex_wp>&, Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void Rank1Update(const Real_wp&, const Vector<Complex_wp>&, Matrix<Complex_wp, Hermitian, RowHermPacked>&);

  SELDON_EXTERN template void Rank2Update(const Real_wp&, const Vector<Real_wp>&, const Vector<Real_wp>&, Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void Rank2Update(const Real_wp&, const Vector<Real_wp>&, const Vector<Real_wp>&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void Rank2Update(const Complex_wp&, const Vector<Complex_wp>&, const Vector<Complex_wp>&, Matrix<Complex_wp, Hermitian, RowHerm>&);
  SELDON_EXTERN template void Rank2Update(const Complex_wp&, const Vector<Complex_wp>&, const Vector<Complex_wp>&, Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  
  SELDON_EXTERN template void Solve(const Matrix<Real_wp, General, RowUpTriang>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Solve(const Matrix<Complex_wp, General, RowUpTriang>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Solve(const Matrix<Real_wp, General, RowLoTriang>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Solve(const Matrix<Complex_wp, General, RowLoTriang>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Solve(const Matrix<Real_wp, General, RowUpTriangPacked>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Solve(const Matrix<Complex_wp, General, RowUpTriangPacked>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Solve(const Matrix<Real_wp, General, RowLoTriangPacked>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Solve(const Matrix<Complex_wp, General, RowLoTriangPacked>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowUpTriang>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowUpTriang>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowLoTriang>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowLoTriang>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowUpTriangPacked>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowLoTriangPacked>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Solve(const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, Vector<Complex_wp>&);
  
  /* Blas Level 3 */

  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, ColMajor>&, const Matrix<Real_wp, General, ColMajor>&, const Real_wp&, Matrix<Real_wp, General, ColMajor>&); 
  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, ColMajor>&, const Matrix<Complex_wp, General, ColMajor>&, const Complex_wp&, Matrix<Complex_wp, General, ColMajor>&); 
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, General, RowMajor>&, const Matrix<Real_wp, General, RowMajor>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&); 
  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, RowMajor>&, const Matrix<Complex_wp, General, RowMajor>&, const Complex_wp&, Matrix<Complex_wp, General, RowMajor>&); 
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Real_wp&, Matrix<Real_wp, General, ColMajor>&); 
  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Complex_wp&, Matrix<Complex_wp, General, ColMajor>&); 
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&); 
  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Complex_wp&, Matrix<Complex_wp, General, RowMajor>&); 
  SELDON_EXTERN template void MltAdd(const SeldonSide&, const Real_wp&, const Matrix<Real_wp, Symmetric, RowSym>&, const Matrix<Real_wp, General, RowMajor>&, const Real_wp&, Matrix<Real_wp, General, RowMajor>&); 
  SELDON_EXTERN template void MltAdd(const SeldonSide&, const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSym>&, const Matrix<Complex_wp, General, RowMajor>&, const Complex_wp&, Matrix<Complex_wp, General, RowMajor>&); 
  SELDON_EXTERN template void MltAdd(SeldonSide const&, Complex_wp const&, Matrix<Complex_wp, Hermitian, RowHerm> const&, Matrix<Complex_wp, General, RowMajor> const&, Complex_wp const&, Matrix<Complex_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void Mlt(const SeldonSide&, const Real_wp&, const Matrix<Real_wp, General, RowLoTriang>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriang>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const Real_wp&, const Matrix<Real_wp, General, RowUpTriang>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriang>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const Real_wp&, const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowLoTriang>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const Complex_wp&, const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowLoTriang>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const Real_wp&, const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowUpTriang>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Mlt(const SeldonSide&, const Complex_wp&, const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowUpTriang>&, Matrix<Complex_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void Solve(const SeldonSide&, const Real_wp&, const Matrix<Real_wp, General, RowLoTriang>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriang>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const Real_wp&, const Matrix<Real_wp, General, RowUpTriang>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriang>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const Real_wp&, const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowLoTriang>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const Complex_wp&, const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowLoTriang>&, Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const Real_wp&, const SeldonTranspose&, const SeldonDiag&, const Matrix<Real_wp, General, RowUpTriang>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void Solve(const SeldonSide&, const Complex_wp&, const SeldonTranspose&, const SeldonDiag&, const Matrix<Complex_wp, General, RowUpTriang>&, Matrix<Complex_wp, General, RowMajor>&);

  /* functions in Functions_Vector */
  
  SELDON_EXTERN template void AddVector(const Real_wp&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
     
  SELDON_EXTERN template void AddVector(const Real_wp&, const Vector<Real_wp, VectSparse>&, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void AddVector(const Real_wp&, const Vector<Complex_wp, VectSparse>&, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void AddVector(const Complex_wp&, const Vector<Complex_wp, VectSparse>&, Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template void AddVector(const Real_wp&, const Vector<Real_wp, VectSparse>&, Vector<Real_wp, VectFull>&);
  SELDON_EXTERN template void AddVector(const Real_wp&, const Vector<Complex_wp, VectSparse>&, Vector<Complex_wp, VectFull>&);
  SELDON_EXTERN template void AddVector(const Complex_wp&, const Vector<Complex_wp, VectSparse>&, Vector<Complex_wp, VectFull>&);

  SELDON_EXTERN template void Swap(Vector<Real_wp, VectSparse>&, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void Swap(Vector<Complex_wp, VectSparse>&, Vector<Complex_wp, VectSparse>&);
  
  SELDON_EXTERN template void MltScalar(const Real_wp&, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Vector<Complex_wp, VectSparse>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template Real_wp DotProdVector(const Vector<Real_wp, VectSparse>&, const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template Complex_wp DotProdVector(const Vector<Complex_wp, VectSparse>&, const Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template Real_wp DotProdConjVector(const Vector<Real_wp, VectSparse>&, const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template Complex_wp DotProdConjVector(const Vector<Complex_wp, VectSparse>&, const Vector<Complex_wp, VectSparse>&);
  
  SELDON_EXTERN template Real_wp DotProdVector(const Vector<Real_wp, VectSparse>&, const Vector<Real_wp, VectFull>&);
  SELDON_EXTERN template Complex_wp DotProdVector(const Vector<Complex_wp, VectSparse>&, const Vector<Complex_wp, VectFull>&);

  SELDON_EXTERN template Real_wp DotProdVector(const Vector<Real_wp, VectFull>&, const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template Complex_wp DotProdVector(const Vector<Complex_wp, VectFull>&, const Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template void ScatterSparseEntry(const Vector<Real_wp, VectSparse>&, Vector<Real_wp, VectFull>&);
  
  SELDON_EXTERN template void CopyVector(const Vector<Real_wp, VectSparse>&, Vector<Real_wp>&);
  SELDON_EXTERN template void CopyVector(const Vector<Real_wp, VectSparse>&, Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void CopyVector(const Vector<Complex_wp, VectSparse>&, Vector<Complex_wp, VectSparse>&);
  
  SELDON_EXTERN template Real_wp Norm1(const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template Real_wp Norm1(const Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template Real_wp Norm2(const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template Real_wp Norm2(const Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template void Conjugate(Vector<Real_wp>&);
  SELDON_EXTERN template void Conjugate(Vector<Complex_wp>&);
  SELDON_EXTERN template void Conjugate(Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void Conjugate(Vector<Complex_wp, VectSparse>&);
  
  SELDON_EXTERN template long GetMaxAbsIndex(const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template long GetMaxAbsIndex(const Vector<Complex_wp, VectSparse>&);

#ifdef SELDON_WITH_MKL
  SELDON_EXTERN template void Add(const Real_wp&, const Vector<Real_wp, VectSparse>&, Vector<Real_wp, VectFull>&);
  SELDON_EXTERN template void Add(const Real_wp&, const Vector<Complex_wp, VectSparse>&, Vector<Complex_wp, VectFull>&);
#endif

}
