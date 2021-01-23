#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/solver/SparseCholeskyFactorisation.cxx"
#ifdef SELDON_WITH_MPI
#include "computation/solver/DistributedCholeskySolver.cxx"
#endif
#endif

namespace Seldon
{
  SELDON_EXTERN template void
  GetCholesky(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int);

  SELDON_EXTERN template void
  SolveCholesky(const SeldonTranspose&,
                const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  SolveCholesky(const SeldonTranspose&,
                const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Real_wp>&);

  SELDON_EXTERN template void
  SolveCholesky(const SeldonTranspose&,
		const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void
  MltCholesky(const SeldonTranspose&,
              const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Real_wp>&);
  
  SELDON_EXTERN template void
  MltCholesky(const SeldonTranspose&,
              const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Real_wp>&);

  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<Complex_wp>&);

  SELDON_EXTERN template class SparseCholeskySolver<Real_wp>;
  
  SELDON_EXTERN template void SparseCholeskySolver<Real_wp>::
  Factorize(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);

  SELDON_EXTERN template void SparseCholeskySolver<Real_wp>::
  Solve(const SeldonTranspose&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void SparseCholeskySolver<Real_wp>::
  Mlt(const SeldonTranspose&, Vector<Real_wp>&, bool);

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template class DistributedCholeskySolver<Real_wp>;

  SELDON_EXTERN template void DistributedCholeskySolver<Real_wp>::
  Factorize(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);

  SELDON_EXTERN template void DistributedCholeskySolver<Real_wp>::
  Factorize(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);

  SELDON_EXTERN template void DistributedCholeskySolver<Real_wp>::
  Solve(const SeldonTranspose&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void DistributedCholeskySolver<Real_wp>::
  Mlt(const SeldonTranspose&, Vector<Real_wp>&, bool);

#endif

}



