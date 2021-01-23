#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/solver/SparseSolver.cxx"
#include "computation/solver/Ordering.cxx"
#include "computation/basic_functions/Functions_Base.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class VirtualSparseDirectSolver<Real_wp>;
  SELDON_EXTERN template class VirtualSparseDirectSolver<Complex_wp>;

  SELDON_EXTERN template class SparseSeldonSolver<Real_wp>;
  SELDON_EXTERN template class SparseSeldonSolver<Complex_wp>;
  
  SELDON_EXTERN template void FindReverseCuthillMcKeeOrdering(const Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<int>&);
  SELDON_EXTERN template void FindReverseCuthillMcKeeOrdering(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Vector<int>&);
  SELDON_EXTERN template void FindReverseCuthillMcKeeOrdering(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<int>&);
  SELDON_EXTERN template void FindReverseCuthillMcKeeOrdering(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Vector<int>&);
  SELDON_EXTERN template void FindReverseCuthillMcKeeOrdering(const Matrix<Real_wp, General, RowSparse>&, Vector<int>&);
  SELDON_EXTERN template void FindReverseCuthillMcKeeOrdering(const Matrix<Complex_wp, General, RowSparse>&, Vector<int>&);
  SELDON_EXTERN template void FindReverseCuthillMcKeeOrdering(const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<int>&);
  SELDON_EXTERN template void FindReverseCuthillMcKeeOrdering(const Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<int>&);

  SELDON_EXTERN template void FindSparseOrdering(Matrix<Real_wp, Symmetric, RowSymSparse>&, Vector<int>&, int);
  SELDON_EXTERN template void FindSparseOrdering(Matrix<Complex_wp, Symmetric, RowSymSparse>&, Vector<int>&, int);
  SELDON_EXTERN template void FindSparseOrdering(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<int>&, int);
  SELDON_EXTERN template void FindSparseOrdering(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Vector<int>&, int);
  SELDON_EXTERN template void FindSparseOrdering(Matrix<Real_wp, General, RowSparse>&, Vector<int>&, int);
  SELDON_EXTERN template void FindSparseOrdering(Matrix<Complex_wp, General, RowSparse>&, Vector<int>&, int);
  SELDON_EXTERN template void FindSparseOrdering(Matrix<Real_wp, General, ArrayRowSparse>&, Vector<int>&, int);
  SELDON_EXTERN template void FindSparseOrdering(Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<int>&, int);
  
  SELDON_EXTERN template void SparseSeldonSolver<Real_wp>::FactorizeMatrix(const IVect&, Matrix<Real_wp, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void SparseSeldonSolver<Complex_wp>::FactorizeMatrix(const IVect&, Matrix<Complex_wp, Symmetric, RowSymSparse>&, bool);
  SELDON_EXTERN template void SparseSeldonSolver<Real_wp>::FactorizeMatrix(const IVect&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void SparseSeldonSolver<Complex_wp>::FactorizeMatrix(const IVect&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, bool);
  SELDON_EXTERN template void SparseSeldonSolver<Real_wp>::FactorizeMatrix(const IVect&, Matrix<Real_wp, General, RowSparse>&, bool);
  SELDON_EXTERN template void SparseSeldonSolver<Complex_wp>::FactorizeMatrix(const IVect&, Matrix<Complex_wp, General, RowSparse>&, bool);
  SELDON_EXTERN template void SparseSeldonSolver<Real_wp>::FactorizeMatrix(const IVect&, Matrix<Real_wp, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void SparseSeldonSolver<Complex_wp>::FactorizeMatrix(const IVect&, Matrix<Complex_wp, General, ArrayRowSparse>&, bool);
  
  SELDON_EXTERN template void SparseSeldonSolver<Real_wp>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void SparseSeldonSolver<Real_wp>::Solve(const SeldonTranspose&, Vector<Real_wp>&);

  SELDON_EXTERN template void SparseSeldonSolver<Real_wp>::Solve(Vector<Complex_wp>&);
  SELDON_EXTERN template void SparseSeldonSolver<Real_wp>::Solve(const SeldonTranspose&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void SparseSeldonSolver<Complex_wp>::Solve(Vector<Complex_wp>&);
  SELDON_EXTERN template void SparseSeldonSolver<Complex_wp>::Solve(const SeldonTranspose&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, ArrayRowSparse>&, IVect&, IVect&, const double&, int);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, ArrayRowSparse>&, IVect&, IVect&, const double&, int);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int);

  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Real_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ArrayRowSparse>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Vector<Complex_wp>&);
  
  
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, RowSymSparse>&, SparseSeldonSolver<Real_wp>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, RowSymSparse>&, SparseSeldonSolver<Complex_wp>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, SparseSeldonSolver<Real_wp>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, SparseSeldonSolver<Complex_wp>&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, RowSparse>&, SparseSeldonSolver<Real_wp>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, RowSparse>&, SparseSeldonSolver<Complex_wp>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, ArrayRowSparse>&, SparseSeldonSolver<Real_wp>&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, ArrayRowSparse>&, SparseSeldonSolver<Complex_wp>&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, RowSymSparse>&, SparseSeldonSolver<Real_wp>&, IVect&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, RowSymSparse>&, SparseSeldonSolver<Complex_wp>&, IVect&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, SparseSeldonSolver<Real_wp>&, IVect&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, SparseSeldonSolver<Complex_wp>&, IVect&, bool);

  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, RowSparse>&, SparseSeldonSolver<Real_wp>&, IVect&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, RowSparse>&, SparseSeldonSolver<Complex_wp>&, IVect&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, ArrayRowSparse>&, SparseSeldonSolver<Real_wp>&, IVect&, bool);
  SELDON_EXTERN template void GetLU(Matrix<Complex_wp, General, ArrayRowSparse>&, SparseSeldonSolver<Complex_wp>&, IVect&, bool);

  
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, RowMajor>&, const Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<int>&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void SolveLU(SparseSeldonSolver<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseSeldonSolver<Real_wp>&, Vector<Real_wp>&);
 
  SELDON_EXTERN template void SolveLU(SparseSeldonSolver<Real_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseSeldonSolver<Real_wp>&, Vector<Complex_wp>&);
   
  //SELDON_EXTERN template void SolveLU(SparseSeldonSolver<Complex_wp>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void SolveLU(const class_SeldonNoTrans&, SparseSeldonSolver<Complex_wp>&, Vector<Real_wp>&);
  //SELDON_EXTERN template void SolveLU(const class_SeldonTrans&, SparseSeldonSolver<Complex_wp>&, Vector<Real_wp>&);
 
  SELDON_EXTERN template void SolveLU(SparseSeldonSolver<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, SparseSeldonSolver<Complex_wp>&, Vector<Complex_wp>&);
  
}



