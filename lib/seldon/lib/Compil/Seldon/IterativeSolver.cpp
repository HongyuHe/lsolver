#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifdef SELDON_WITH_MPI
#include "SeldonDistributedHeader.hxx"
#include "SeldonDistributedInline.hxx"
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/solver/iterative/Iterative.cxx"
#include "computation/solver/preconditioner/Precond_Ssor.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template class Iteration<Real_wp>;  
  SELDON_EXTERN template int Iteration<Real_wp>::Init(const Vector<Real_wp>&);
  SELDON_EXTERN template int Iteration<Real_wp>::Init(const Vector<Complex_wp >&);
  SELDON_EXTERN template bool Iteration<Real_wp>::Finished(const Vector<Real_wp>&) const;
  SELDON_EXTERN template bool Iteration<Real_wp>::Finished(const Vector<Complex_wp >&) const;

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template int Iteration<Real_wp>::Init(const DistributedVector<Real_wp>&);
  SELDON_EXTERN template int Iteration<Real_wp>::Init(const DistributedVector<Complex_wp >&);
  SELDON_EXTERN template bool Iteration<Real_wp>::Finished(const DistributedVector<Real_wp>&) const;
  SELDON_EXTERN template bool Iteration<Real_wp>::Finished(const DistributedVector<Complex_wp >&) const;
#endif

#ifdef SELDON_WITH_PRECONDITIONING
  SELDON_EXTERN template class SorPreconditioner<Real_wp>;
  SELDON_EXTERN template class SorPreconditioner<Complex_wp >;
  
  SELDON_EXTERN template void SorPreconditioner<Real_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Real_wp>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void SorPreconditioner<Real_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Real_wp>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
  SELDON_EXTERN template void SorPreconditioner<Complex_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Complex_wp>&, const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);
#endif

  SELDON_EXTERN template int BiCg(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);

  SELDON_EXTERN template int BiCg(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int BiCgStab(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				      Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int BiCgStab(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				      Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int BiCgStabl(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				       Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int BiCgStabl(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				       Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int BiCgcr(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				      Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int BiCgcr(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				    Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Cg(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Cg(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Cgne(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Cgne(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Cgs(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				 Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Cgs(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				 Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int CoCg(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int CoCg(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Gcr(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				 Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Gcr(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				 Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Gmres(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				   Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Gmres(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				   Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Lsqr(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Lsqr(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int MinRes(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				    Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int MinRes(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				    Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int QCgs(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int QCgs(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Qmr(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				 Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Qmr(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				 Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int QmrSym(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				    Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int QmrSym(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				    Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Symmlq(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				    Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Symmlq(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				    Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int TfQmr(const VirtualMatrix<Real_wp>&, Vector<Real_wp>&, const Vector<Real_wp>&,
				   Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int TfQmr(const VirtualMatrix<Complex_wp >&, Vector<Complex_wp >&, const Vector<Complex_wp >&,
				   Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);


  // distributed vectors
#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template int BiCg(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);

  SELDON_EXTERN template int BiCg(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int BiCgStab(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				      Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int BiCgStab(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				      Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int BiCgStabl(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				       Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int BiCgStabl(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				       Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int BiCgcr(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				      Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int BiCgcr(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				    Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Cg(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Cg(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Cgne(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Cgne(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Cgs(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				 Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Cgs(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				 Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int CoCg(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int CoCg(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Gcr(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				 Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Gcr(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				 Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Gmres(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				   Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Gmres(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				   Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Lsqr(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Lsqr(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int MinRes(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				    Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int MinRes(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				    Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int QCgs(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				  Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int QCgs(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				  Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Qmr(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				 Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Qmr(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				 Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int QmrSym(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				    Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int QmrSym(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				    Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int Symmlq(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				    Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int Symmlq(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				    Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);

  SELDON_EXTERN template int TfQmr(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&,
				   Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);
  
  SELDON_EXTERN template int TfQmr(const VirtualMatrix<Complex_wp >&, DistributedVector<Complex_wp >&, const DistributedVector<Complex_wp >&,
				   Preconditioner_Base<Complex_wp >&, Iteration<Real_wp>&);
#endif

}
