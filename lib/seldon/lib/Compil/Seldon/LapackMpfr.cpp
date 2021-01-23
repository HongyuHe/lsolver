#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/basic_functions/Functions_Base.cxx"
#endif

namespace Seldon
{
  /* Linear equations */
       
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, ColMajor>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, ColMajor>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, RowMajor>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, RowMajor>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, RowSym>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, RowSym>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, RowSymPacked>&, Vector<int>&, LapackInfo&);
  //SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<int>&, LapackInfo&);
  //SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<int>&, LapackInfo&);
  
  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, General, ColMajor>&, const Vector<int>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Complex_wp, General, ColMajor>&, const Vector<int>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, General, RowMajor>&, const Vector<int>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Complex_wp, General, RowMajor>&, const Vector<int>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<int>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<int>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<int>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<int>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<int>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<int>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<int>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<int>&, Vector<Complex_wp>&, LapackInfo&);
  //SELDON_EXTERN template void SolveLuVector(const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<int>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void SolveLuVector(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<int>&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, RowUpTriang>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Complex_wp, General, RowUpTriang>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, RowLoTriang>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Complex_wp, General, RowLoTriang>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, RowUpTriangPacked>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Complex_wp, General, RowUpTriangPacked>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, RowLoTriangPacked>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Complex_wp, General, RowLoTriangPacked>&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, RowMajor>&, const Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, ColMajor>&, const Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, General, RowMajor>&, LapackInfo&);
  //SELDON_EXTERN template void GetInverse(Matrix<Real_wp, Symmetric, RowSym>&);
  //SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, Symmetric, RowSymPacked>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, Symmetric, RowSymPacked>&, LapackInfo&);
  //SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, Hermitian, RowHerm>&);
  //SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, Hermitian, RowHermPacked>&);
  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, General, RowUpTriang>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, General, RowUpTriang>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, General, RowLoTriang>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, General, RowLoTriang>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, General, RowUpTriangPacked>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, General, RowUpTriangPacked>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, General, RowLoTriangPacked>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, General, RowLoTriangPacked>&, LapackInfo&);
  
  SELDON_EXTERN template void GetCholesky(Matrix<Real_wp, Symmetric, RowSym>&, LapackInfo&);
  SELDON_EXTERN template void GetCholesky(Matrix<Real_wp, Symmetric, RowSymPacked>&, LapackInfo&);
  SELDON_EXTERN template void GetCholesky(Matrix<Complex_wp, Hermitian, RowHerm>&, LapackInfo&);
  SELDON_EXTERN template void GetCholesky(Matrix<Complex_wp, Hermitian, RowHermPacked>&, LapackInfo&);
  
  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Complex_wp>&, LapackInfo&);
  
  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Complex_wp>&, LapackInfo&);
    
  /* Least-squares */
  
  SELDON_EXTERN template void GetQR(Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetQR(Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetQR(Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetQR(Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp>&, LapackInfo&);
  //SELDON_EXTERN template void GetLQ(Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetLQ(Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void GetLQ(Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetLQ(Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp>&);
  
  //SELDON_EXTERN template void GetQ_FromQR(Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetQ_FromQR(Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void GetQ_FromQR(Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetQ_FromQR(Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void GetQ_FromLQ(Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetQ_FromLQ(Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void GetQ_FromLQ(Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetQ_FromLQ(Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void MltQ_FromQR(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void MltQ_FromQR(const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void MltQ_FromQR(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void MltQ_FromQR(const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, LapackInfo&);

  //SELDON_EXTERN template void MltQ_FromLQ(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void MltQ_FromLQ(const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void MltQ_FromLQ(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void MltQ_FromLQ(const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void MltQ_FromLQ(const SeldonSide&, const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, Matrix<Real_wp, General, ColMajor>&);
  //SELDON_EXTERN template void MltQ_FromLQ(const SeldonSide&, const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, Matrix<Complex_wp, General, ColMajor>&);
  //SELDON_EXTERN template void MltQ_FromLQ(const SeldonSide&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&);
  //SELDON_EXTERN template void MltQ_FromLQ(const SeldonSide&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, Matrix<Complex_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void SolveQR(const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveQR(const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveQR(const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void SolveQR(const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, LapackInfo&);
  //SELDON_EXTERN template void SolveLQ(const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void SolveLQ(const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void SolveLQ(const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void SolveLQ(const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  
  /* Eigenvalue problems */
  
  SELDON_EXTERN template void GetEigenvalues(Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, Vector<Real_wp>&, Matrix<Real_wp, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&, Matrix<Complex_wp, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp>&, Matrix<Complex_wp, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, Symmetric, RowSym>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, Symmetric, RowSym>&, Vector<Complex_wp>&, Matrix<Complex_wp, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp>&, Matrix<Complex_wp, General, RowMajor>&, LapackInfo&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Real_wp>&, Matrix<Complex_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Real_wp>&, Matrix<Complex_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Real_wp, Symmetric, RowSym>&, Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Real_wp, Symmetric, RowSym>&, Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, Symmetric, RowSym>&, Matrix<Complex_wp, Symmetric, RowSym>&, Vector<Complex_wp >&, Vector<Complex_wp >&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, Symmetric, RowSym>&, Matrix<Complex_wp, Symmetric, RowSym>&, Vector<Complex_wp >&, Vector<Complex_wp >&, Matrix<Complex_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp >&, Vector<Complex_wp >&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp >&, Vector<Complex_wp >&, Matrix<Complex_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, Hermitian, RowHerm>&, Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, Hermitian, RowHerm>&, Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Real_wp>&, Matrix<Complex_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, Hermitian, RowHermPacked>&, Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, Hermitian, RowHermPacked>&, Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Real_wp>&, Matrix<Complex_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Real_wp, General, ColMajor>&, Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Real_wp, General, RowMajor>&, Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Real_wp, General, ColMajor>&, Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&, Matrix<Real_wp, General, ColMajor>&);
  // SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Real_wp, General, RowMajor>&, Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void GetEigenvalues(Matrix<Complex_wp, General, RowMajor>&, Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp>&, Vector<Complex_wp>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&, Vector<Complex_wp>&, Matrix<Complex_wp, General, ColMajor>&);
  //SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<Complex_wp, General, RowMajor>&, Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp>&, Vector<Complex_wp>&, Matrix<Complex_wp, General, RowMajor>&);
  
  //SELDON_EXTERN template void GetSVD(Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, Matrix<Real_wp, General, ColMajor>&, Matrix<Real_wp, General, ColMajor>&);
  //SELDON_EXTERN template void GetSVD(Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, Matrix<Real_wp, General, RowMajor>&, Matrix<Real_wp, General, RowMajor>&);
  //SELDON_EXTERN template void GetSVD(Matrix<Complex_wp, General, ColMajor>&, Vector<Real_wp>&, Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&);
  //SELDON_EXTERN template void GetSVD(Matrix<Complex_wp, General, RowMajor>&, Vector<Real_wp>&, Matrix<Complex_wp, General, RowMajor>&, Matrix<Complex_wp, General, RowMajor>&);
  
  SELDON_EXTERN template void GetPseudoInverse(Matrix<Real_wp, General, ColMajor>&, const Real_wp&, LapackInfo&);
  SELDON_EXTERN template void GetPseudoInverse(Matrix<Real_wp, General, RowMajor>&, const Real_wp&, LapackInfo&);
  SELDON_EXTERN template void GetPseudoInverse(Matrix<Complex_wp, General, ColMajor>&, const Real_wp&, LapackInfo&);
  SELDON_EXTERN template void GetPseudoInverse(Matrix<Complex_wp, General, RowMajor>&, const Real_wp&, LapackInfo&);
  
}

