#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/solver/preconditioner/IlutPreconditioning.cxx"
#include "computation/solver/preconditioner/SymmetricIlutPreconditioning.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class IlutPreconditioning<Real_wp>;
  SELDON_EXTERN template class IlutPreconditioning<Complex_wp >;
  
  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::FactorizeSymMatrix(const IVect&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::FactorizeUnsymMatrix(const IVect&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::FactorizeUnsymMatrix(const IVect&, Matrix<Real_wp, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void IlutPreconditioning<Complex_wp >::FactorizeSymMatrix(const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void IlutPreconditioning<Complex_wp >::FactorizeUnsymMatrix(const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void IlutPreconditioning<Complex_wp >::FactorizeUnsymMatrix(const IVect&, Matrix<Complex_wp, General, ArrayRowSparse>&, bool);
  
  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::FactorizeMatrix(const IVect&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::FactorizeMatrix(const IVect&, Matrix<Real_wp, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void IlutPreconditioning<Complex_wp >::FactorizeMatrix(const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void IlutPreconditioning<Complex_wp >::FactorizeMatrix(const IVect&, Matrix<Complex_wp, General, ArrayRowSparse>&, bool);

  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::Solve(const SeldonTranspose&, Vector<Real_wp>&);
  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::Solve(Vector<Complex_wp>&);
  SELDON_EXTERN template void IlutPreconditioning<Real_wp>::Solve(const SeldonTranspose&, Vector<Complex_wp>&);

  SELDON_EXTERN template void IlutPreconditioning<Complex_wp >::Solve(Vector<Complex_wp >&);
  SELDON_EXTERN template void IlutPreconditioning<Complex_wp >::Solve(const SeldonTranspose&, Vector<Complex_wp >&);

  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, ArrayRowSparse>&,
				    IlutPreconditioning<Real_wp>&, IVect&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
				    IlutPreconditioning<Real_wp>&, IVect&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, ArrayRowSparse>&,
				    IlutPreconditioning<Complex_wp >&, IVect&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
				    IlutPreconditioning<Complex_wp >&, IVect&, bool);


  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, RowSparse>&,
				    IlutPreconditioning<Real_wp>&, IVect&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, RowSymSparse>&,
				    IlutPreconditioning<Real_wp>&, IVect&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, RowSparse>&,
				    IlutPreconditioning<Complex_wp >&, IVect&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, RowSymSparse>&,
				    IlutPreconditioning<Complex_wp >&, IVect&, bool);
  
}
