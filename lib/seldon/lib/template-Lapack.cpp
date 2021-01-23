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
#include "share/Allocator.cxx"
#include "share/StorageInline.cxx"
#include "share/MatrixFlagInline.cxx"
#include "computation/interfaces/Lapack_LinearEquations.cxx"
#include "computation/interfaces/Lapack_LeastSquares.cxx"
#include "computation/interfaces/Lapack_Eigenvalues.cxx"
#endif

#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;

namespace Seldon
{
  /* Linear equations */

  SELDON_EXTERN template void GetLU(Matrix<@real_complex, General, @storage_blasGE>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<@real_complex, Symmetric, @storage_blasS>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<@complex, Hermitian, @storage_blasH>&, Vector<int>&, LapackInfo&);

  SELDON_EXTERN template void SolveLuVector(const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<int>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<int>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<@real_complex, Symmetric, @storage_blasS>&, const Vector<int>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<@complex, Hermitian, @storage_blasH>&, const Vector<int>&, Vector<@complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const Matrix<@real_complex, General, @storage_blasT>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, @storage_blasT>&, Vector<@real_complex>&, LapackInfo&);

  SELDON_EXTERN template @real ReciprocalConditionNumber(const Matrix<@real, General, @storage_blasGE>&, Vector<int>&, SeldonNorm, @real, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const Matrix<complexdouble, General, @storage_blasGE>&, Vector<int>&, SeldonNorm, double, LapackInfo&);
  SELDON_EXTERN template @real ReciprocalConditionNumber(const Matrix<@real, Symmetric, @storage_blasS>&, Vector<int>&, SeldonNorm, @real, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const Matrix<complexdouble, Symmetric, @storage_blasS>&, Vector<int>&, SeldonNorm, double, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const Matrix<complexdouble, Hermitian, @storage_blasH>&, Vector<int>&, SeldonNorm, double, LapackInfo&);
  SELDON_EXTERN template @real ReciprocalConditionNumber(const Matrix<@real, General, @storage_blasT>&, SeldonNorm, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const Matrix<complexdouble, General, @storage_blasT>&, SeldonNorm, LapackInfo&);
  SELDON_EXTERN template @real ReciprocalConditionNumber(const SeldonDiag&, const Matrix<@real, General, @storage_blasT>&, SeldonNorm, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const SeldonDiag&, const Matrix<complexdouble, General, @storage_blasT>&, SeldonNorm, LapackInfo&);

  SELDON_EXTERN template void RefineSolutionLU(const Matrix<@real, General, @storage_blasGE>&, const Matrix<@real, General, @storage_blasGE>&, const Vector<int>&, Vector<@real>&, const Vector<@real>&, @real&, @real&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, General, RowMajor>&, const Matrix<complexdouble, General, RowMajor>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, General, ColMajor>&, const Matrix<complexdouble, General, ColMajor>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const Matrix<@real, General, @storage_blasGE>&, const Matrix<@real, General, @storage_blasGE>&, const Vector<int>&, Vector<@real>&, const Vector<@real>&, @real&, @real&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const Matrix<complexdouble, General, ColMajor>&, const Matrix<complexdouble, General, ColMajor>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const Matrix<complexdouble, General, RowMajor>&, const Matrix<complexdouble, General, RowMajor>&, const Vector<int>&, Vector<complexdouble>&, Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<@real, Symmetric, @storage_blasS>&, const Matrix<@real, Symmetric, @storage_blasS>&, const Vector<int>&, Vector<@real>&, const Vector<@real>&, @real&, @real&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Symmetric, @storage_blasS>&, const Matrix<complexdouble, Symmetric, @storage_blasS>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Hermitian, RowHerm>&, const Matrix<complexdouble, Hermitian, RowHerm>&, const Vector<int>&, Vector<complexdouble>&, Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Hermitian, RowHermPacked>&, const Matrix<complexdouble, Hermitian, RowHermPacked>&, const Vector<int>&, Vector<complexdouble>&, Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Hermitian, ColHerm>&, const Matrix<complexdouble, Hermitian, ColHerm>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Hermitian, ColHermPacked>&, const Matrix<complexdouble, Hermitian, ColHermPacked>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<double, General, @storage_blasT>&, Vector<double>&, const Vector<double>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const SeldonDiag&, const Matrix<double, General, @storage_blasT>&, Vector<double>&, const Vector<double>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, General, @storage_blasT>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const SeldonDiag&, const Matrix<complexdouble, General, @storage_blasTCol>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const SeldonDiag&, const Matrix<complexdouble, General, @storage_blasTRow>&, Vector<complexdouble>&, Vector<complexdouble>&, double&, double&, LapackInfo&);

  SELDON_EXTERN template void GetInverse(Matrix<@real_complex, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<@real_complex, Symmetric, @storage_blasS>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<@complex, Hermitian, @storage_blasH>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<@real_complex, General, @storage_blasT>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(const SeldonDiag&, Matrix<@real_complex, General, @storage_blasT>&, LapackInfo&);

  SELDON_EXTERN template void GetScalingFactors(const Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, @real&, @real&, @real&, LapackInfo&);
  SELDON_EXTERN template void GetScalingFactors(const Matrix<complexdouble, General, @storage_blasGE>&, Vector<double>&, Vector<double>&, double&, double&, double&, LapackInfo&);

  SELDON_EXTERN template void GetCholesky(Matrix<double, Symmetric, @storage_blasS>&, LapackInfo&);
  SELDON_EXTERN template void GetCholesky(Matrix<complexdouble, Hermitian, @storage_blasH>&, LapackInfo&);

  SELDON_EXTERN template void SolveCholesky(const @trans&, const Matrix<double, Symmetric, @storage_blasS>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void SolveCholesky(const @trans&, const Matrix<complexdouble, Hermitian, @storage_blasH>&, Vector<complexdouble>&, LapackInfo&);

  SELDON_EXTERN template void MltCholesky(const @trans&, const Matrix<double, Symmetric, @storage_blasS>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void MltCholesky(const @trans&, const Matrix<complexdouble, Hermitian, @storage_blasH>&, Vector<complexdouble>&, LapackInfo&);

  /* Least-squares */

  SELDON_EXTERN template void GetQR(Matrix<@real_complex, General, @storage_blasGE>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void GetLQ(Matrix<@real_complex, General, @storage_blasGE>&, Vector<@real_complex>&, LapackInfo&);

  SELDON_EXTERN template void GetQ_FromQR(Matrix<@real_complex, General, @storage_blasGE>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void GetQ_FromLQ(Matrix<@real_complex, General, @storage_blasGE>&, Vector<@real_complex>&, LapackInfo&);

  SELDON_EXTERN template void MltQ_FromQR(const @trans&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void MltQ_FromQR(const @side&, const @trans&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, Matrix<@real_complex, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void MltQ_FromLQ(const @trans&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void MltQ_FromLQ(const @side&, const @trans&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, Matrix<@real_complex, General, @storage_blasGE>&, LapackInfo&);

  SELDON_EXTERN template void SolveQR(const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLQ(const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<@real_complex>&, Vector<@real_complex>&, LapackInfo&);

  SELDON_EXTERN template void GetLowerTriangular(const Matrix<@real_complex, General, @storage_blasGE>&, Matrix<@real_complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void GetUpperTriangular(const Matrix<@real_complex, General, @storage_blasGE>&, Matrix<@real_complex, General, @storage_blasGE>&);

  SELDON_EXTERN template void GetQR(const Matrix<@real_complex, General, @storage_blasGE>&, Matrix<@real_complex, General, @storage_blasGE>&, Matrix<@real_complex, General, @storage_blasGE>&);
  SELDON_EXTERN template void GetLQ(const Matrix<@real_complex, General, @storage_blasGE>&, Matrix<@real_complex, General, @storage_blasGE>&, Matrix<@real_complex, General, @storage_blasGE>&);

  /* Eigenvalue problems */

  SELDON_EXTERN template void GetEigenvalues(Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@complex, General, @storage_blasGE>&, Vector<@complex>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, Matrix<@real, General, @storage_blasGE>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@complex, General, @storage_blasGE>&, Vector<@complex>&, Matrix<@complex, General, @storage_blasGE>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@real_complex, Symmetric, RowSym>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@real_complex, Symmetric, RowSym>&, Vector<@real_complex>&, Matrix<@real_complex, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@real_complex, Symmetric, ColSym>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@real_complex, Symmetric, ColSym>&, Vector<@real_complex>&, Matrix<@real_complex, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@real_complex, Symmetric, RowSymPacked>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@real_complex, Symmetric, RowSymPacked>&, Vector<@real_complex>&, Matrix<@real_complex, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@real_complex, Symmetric, ColSymPacked>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@real_complex, Symmetric, ColSymPacked>&, Vector<@real_complex>&, Matrix<@real_complex, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Hermitian, RowHerm>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Hermitian, RowHerm>&, Vector<double>&, Matrix<complex<double>, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Hermitian, ColHerm>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Hermitian, ColHerm>&, Vector<double>&, Matrix<complex<double>, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Hermitian, RowHermPacked>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Hermitian, RowHermPacked>&, Vector<double>&, Matrix<complex<double>, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Hermitian, ColHermPacked>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Hermitian, ColHermPacked>&, Vector<double>&, Matrix<complex<double>, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<double, Symmetric, RowSym>&, Matrix<double, Symmetric, RowSym>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<double, Symmetric, RowSym>&, Matrix<double, Symmetric, RowSym>&, Vector<double>&, Matrix<double, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Symmetric, RowSym>&, Matrix<complex<double>, Symmetric, RowSym>&, Vector<complex<double> >&, Vector<complex<double> >&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Symmetric, RowSym>&, Matrix<complex<double>, Symmetric, RowSym>&, Vector<complex<double> >&, Vector<complex<double> >&, Matrix<complex<double>, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<double, Symmetric, ColSym>&, Matrix<double, Symmetric, ColSym>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<double, Symmetric, ColSym>&, Matrix<double, Symmetric, ColSym>&, Vector<double>&, Matrix<double, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Symmetric, ColSym>&, Matrix<complex<double>, Symmetric, ColSym>&, Vector<complex<double> >&, Vector<complex<double> >&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Symmetric, ColSym>&, Matrix<complex<double>, Symmetric, ColSym>&, Vector<complex<double> >&, Vector<complex<double> >&, Matrix<complex<double>, General, ColMajor>&, LapackInfo&);
    SELDON_EXTERN template void GetEigenvalues(Matrix<double, Symmetric, RowSymPacked>&, Matrix<double, Symmetric, RowSymPacked>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<double, Symmetric, RowSymPacked>&, Matrix<double, Symmetric, RowSymPacked>&, Vector<double>&, Matrix<double, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Symmetric, RowSymPacked>&, Matrix<complex<double>, Symmetric, RowSymPacked>&, Vector<complex<double> >&, Vector<complex<double> >&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Symmetric, RowSymPacked>&, Matrix<complex<double>, Symmetric, RowSymPacked>&, Vector<complex<double> >&, Vector<complex<double> >&, Matrix<complex<double>, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<double, Symmetric, ColSymPacked>&, Matrix<double, Symmetric, ColSymPacked>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<double, Symmetric, ColSymPacked>&, Matrix<double, Symmetric, ColSymPacked>&, Vector<double>&, Matrix<double, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Symmetric, ColSymPacked>&, Matrix<complex<double>, Symmetric, ColSymPacked>&, Vector<complex<double> >&, Vector<complex<double> >&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Symmetric, ColSymPacked>&, Matrix<complex<double>, Symmetric, ColSymPacked>&, Vector<complex<double> >&, Vector<complex<double> >&, Matrix<complex<double>, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Hermitian, RowHerm>&, Matrix<complex<double>, Hermitian, RowHerm>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Hermitian, RowHerm>&, Matrix<complex<double>, Hermitian, RowHerm>&, Vector<double>&, Matrix<complex<double>, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Hermitian, ColHerm>&, Matrix<complex<double>, Hermitian, ColHerm>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Hermitian, ColHerm>&, Matrix<complex<double>, Hermitian, ColHerm>&, Vector<double>&, Matrix<complex<double>, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Hermitian, RowHermPacked>&, Matrix<complex<double>, Hermitian, RowHermPacked>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Hermitian, RowHermPacked>&, Matrix<complex<double>, Hermitian, RowHermPacked>&, Vector<double>&, Matrix<complex<double>, General, RowMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<complex<double>, Hermitian, ColHermPacked>&, Matrix<complex<double>, Hermitian, ColHermPacked>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<complex<double>, Hermitian, ColHermPacked>&, Matrix<complex<double>, Hermitian, ColHermPacked>&, Vector<double>&, Matrix<complex<double>, General, ColMajor>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@real, General, @storage_blasGE>&, Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, Vector<@real>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@real, General, @storage_blasGE>&, Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, Vector<@real>&, Matrix<@real, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@complex, General, @storage_blasGE>&, Matrix<@complex, General, @storage_blasGE>&, Vector<@complex>&, Vector<@complex>&, LapackInfo&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@complex, General, @storage_blasGE>&, Matrix<@complex, General, @storage_blasGE>&, Vector<@complex>&, Vector<@complex>&, Matrix<@complex, General, @storage_blasGE>&, LapackInfo&);

  SELDON_EXTERN template void GetSVD(Matrix<double, General, @storage_blasGE>&, Vector<double>&, Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetSVD(Matrix<complex<double>, General, @storage_blasGE>&, Vector<double>&, Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, LapackInfo&);

  SELDON_EXTERN template void GetPseudoInverse(Matrix<double, General, @storage_blasGE>&, const double&, LapackInfo&);
  SELDON_EXTERN template void GetPseudoInverse(Matrix<complex<double>, General, @storage_blasGE>&, const double&, LapackInfo&);

  SELDON_EXTERN template void GetHessenberg(Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetHessenberg(Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetQZ(Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetHessenberg(Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetHessenberg(Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetQZ(Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, LapackInfo&);

  SELDON_EXTERN template void SolveHessenberg(Matrix<double, General, @storage_blasGE>&, Vector<double>&);
  SELDON_EXTERN template void SolveHessenberg(Matrix<complex<double>, General, @storage_blasGE>&, Vector<complex<double> >&);
  SELDON_EXTERN template void SolveHessenbergTwo(Matrix<double, General, @storage_blasGE>&, Vector<double>&);
  SELDON_EXTERN template void SolveHessenbergTwo(Matrix<complex<double>, General, @storage_blasGE>&, Vector<complex<double> >&);

  SELDON_EXTERN template void SolveSylvester(Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&, Matrix<double, General, @storage_blasGE>&);
  SELDON_EXTERN template void SolveSylvester(Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&, Matrix<complex<double>, General, @storage_blasGE>&);

}
