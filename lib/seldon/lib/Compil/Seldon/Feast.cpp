#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/eigenvalue/Feast.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void FindEigenvaluesFeast(EigenProblem_Base<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&,
						   Matrix<Real_wp, General, ColMajor>&);

  SELDON_EXTERN template void FindEigenvaluesFeast(EigenProblem_Base<Complex_wp >&, Vector<Complex_wp >&,
						   Vector<Complex_wp >&,
						   Matrix<Complex_wp, General, ColMajor>&);
  
}
