#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifdef SELDON_WITH_MPI
#include "SeldonDistributedHeader.hxx"
#include "SeldonDistributedInline.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/eigenvalue/VirtualEigenvalueSolver.cxx"
#include "computation/interfaces/eigenvalue/PolynomialEigenvalueSolver.cxx"
#include "computation/interfaces/eigenvalue/NonLinearEigenvalueSolver.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void SortEigenvalues(Vector<Real_wp>&, Vector<Real_wp>&,
                                              Matrix<Real_wp, General, RowMajor>&,
                                              Matrix<Real_wp, General, RowMajor>&,
                                              int, int, const Real_wp&, const Real_wp&);
  
  SELDON_EXTERN template void SortEigenvalues(Vector<Complex_wp >&, Vector<Complex_wp >&,
                                              Matrix<Complex_wp , General, RowMajor>&,
                                              Matrix<Complex_wp , General, RowMajor>&,
                                              int, int, const Complex_wp &, const Complex_wp &);

  SELDON_EXTERN template void SortEigenvalues(Vector<Real_wp>&, Vector<Real_wp>&,
                                              Matrix<Real_wp, General, ColMajor>&,
                                              Matrix<Real_wp, General, ColMajor>&,
                                              int, int, const Real_wp&, const Real_wp&);
  
  SELDON_EXTERN template void SortEigenvalues(Vector<Complex_wp >&, Vector<Complex_wp >&,
                                              Matrix<Complex_wp , General, ColMajor>&,
                                              Matrix<Complex_wp , General, ColMajor>&,
                                              int, int, const Complex_wp &, const Complex_wp &);
  
  SELDON_EXTERN template class EigenProblem_Base<Real_wp>;
  SELDON_EXTERN template class EigenProblem_Base<Complex_wp >;
  
  SELDON_EXTERN template class VirtualEigenProblem<Real_wp, Real_wp, Real_wp>;
  SELDON_EXTERN template class VirtualEigenProblem<Complex_wp, Complex_wp, Real_wp>;
  SELDON_EXTERN template class VirtualEigenProblem<Complex_wp, Real_wp, Complex_wp>;
  SELDON_EXTERN template class VirtualEigenProblem<Complex_wp, Complex_wp, Complex_wp>;
  
  SELDON_EXTERN template void ApplyScalingEigenvec(EigenProblem_Base<Real_wp>&,
                                                   Vector<Real_wp>&, Vector<Real_wp>&, Matrix<Real_wp, General, ColMajor>&,
                                                   const Real_wp&, const Real_wp&);

  SELDON_EXTERN template void ApplyScalingEigenvec(EigenProblem_Base<Complex_wp >&,
                                                   Vector<Complex_wp >&, Vector<Complex_wp >&, Matrix<Complex_wp, General, ColMajor>&,
                                                   const Complex_wp&, const Complex_wp&);

  SELDON_EXTERN template class DenseEigenProblem<Real_wp, Real_wp, General, RowMajor, Real_wp, Symmetric, RowSymPacked>;
  SELDON_EXTERN template class DenseEigenProblem<Complex_wp , Complex_wp, General, RowMajor, Real_wp , Symmetric, RowSymPacked>;
  SELDON_EXTERN template class DenseEigenProblem<Real_wp, Real_wp, Symmetric, RowSymPacked, Real_wp, Symmetric, RowSymPacked>;

  /* Sparse eigenproblem with real symmetric mass matrix */
  
  SELDON_EXTERN template class SparseEigenProblem<Real_wp, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>, Matrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>, Matrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Real_wp, Matrix<Real_wp, General, ArrayRowSparse>, Matrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Complex_wp, General, ArrayRowSparse>, Matrix<Real_wp, Symmetric, ArrayRowSymSparse> >;

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template class SparseEigenProblem<Real_wp, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Real_wp, DistributedMatrix<Real_wp, General, ArrayRowSparse>, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, DistributedMatrix<Complex_wp, General, ArrayRowSparse>, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, DistributedMatrix<Complex_wp, General, ArrayRowSparse>, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> >;
#endif

  SELDON_EXTERN template void GetEigenvaluesEigenvectors(EigenProblem_Base<Real_wp>& var_eig, Vector<Real_wp>& lambda, Vector<Real_wp>& lambda_imag, Matrix<Real_wp, General, ColMajor>& eigen_vec, int type);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(EigenProblem_Base<Real_wp>& var_eig, Vector<Real_wp>& lambda, Vector<Real_wp>& lambda_imag, Matrix<Real_wp, General, RowMajor>& eigen_vec, int type);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(EigenProblem_Base<Complex_wp >& var_eig, Vector<Complex_wp >& lambda, Vector<Complex_wp >& lambda_imag, Matrix<Complex_wp, General, ColMajor>& eigen_vec, int type);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(EigenProblem_Base<Complex_wp >& var_eig, Vector<Complex_wp >& lambda, Vector<Complex_wp >& lambda_imag, Matrix<Complex_wp, General, RowMajor>& eigen_vec, int type);

  /* Sparse eigenproblem with complex symmetric mass matrix */
  
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Real_wp, General, ArrayRowSparse>, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Complex_wp, General, ArrayRowSparse>, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> >;

  /* Sparse eigenproblem with complex unsymmetric mass matrix */
  
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Real_wp, General, ArrayRowSparse>, Matrix<Complex_wp, General, ArrayRowSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>, Matrix<Complex_wp, General, ArrayRowSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Complex_wp, General, ArrayRowSparse>, Matrix<Complex_wp, General, ArrayRowSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>, Matrix<Complex_wp, General, ArrayRowSparse> >;

  /* Sparse eigenproblem with real unsymmetric mass matrix */
  
  SELDON_EXTERN template class SparseEigenProblem<Real_wp, Matrix<Real_wp, General, ArrayRowSparse>, Matrix<Real_wp, General, ArrayRowSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Real_wp, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>, Matrix<Real_wp, General, ArrayRowSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Complex_wp, General, ArrayRowSparse>, Matrix<Real_wp, General, ArrayRowSparse> >;
  SELDON_EXTERN template class SparseEigenProblem<Complex_wp, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>, Matrix<Real_wp, General, ArrayRowSparse> >;

  /* Polynomial Eigenvalue Problems */
  
  SELDON_EXTERN template class PolynomialEigenProblem_Base<Real_wp>;
  SELDON_EXTERN template class PolynomialEigenProblem_Base<Complex_wp >;

  SELDON_EXTERN template class PolynomialEigenProblem<Real_wp>;
  SELDON_EXTERN template class PolynomialEigenProblem<Complex_wp >;

  SELDON_EXTERN template class PolynomialDenseEigenProblem<Real_wp, General, RowMajor>;
#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template class PolynomialSparseEigenProblem<Real_wp, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class PolynomialSparseEigenProblem<Complex_wp, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> >;
#else
  SELDON_EXTERN template class PolynomialSparseEigenProblem<Real_wp, Matrix<Real_wp, Symmetric, ArrayRowSymSparse> >;
  SELDON_EXTERN template class PolynomialSparseEigenProblem<Complex_wp, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> >;
#endif

  SELDON_EXTERN template class NonLinearEigenProblem_Base<Real_wp>;
  SELDON_EXTERN template class NonLinearEigenProblem_Base<Complex_wp >;

}
