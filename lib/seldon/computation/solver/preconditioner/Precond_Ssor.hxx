// Copyright (C) 2003-2009 Marc Duruflé
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_FILE_PRECOND_SSOR_HXX

namespace Seldon
{
  
  template <class T>
  class SorPreconditioner : public Preconditioner_Base<T>
  {
  protected :
    bool symmetric_precond; //!< true for Symmetric relaxation
    int nb_iter; //!< number of iterations
    typename ClassComplexType<T>::Treal omega; //!< relaxation parameter

  public :
    SorPreconditioner();
    
    bool IsSymmetric() const;
    void InitSymmetricPreconditioning();
    void InitUnSymmetricPreconditioning();
    
    void SetParameterRelaxation(const typename ClassComplexType<T>::Treal& param);
    void SetNumberIterations(int nb_iterations);
    void CopyParameter(const Preconditioner_Base<T>&);

#ifdef SELDON_WITH_VIRTUAL
    template<class T0>
    void SolveGen(const SeldonTranspose&, const VirtualMatrix<T>&, const Vector<T0>&, Vector<T0>&,
		  bool init_guess_null = true);
    
    void Solve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>& z);
    void TransSolve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>&);

    void Solve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>& z, bool init);
    void TransSolve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>&, bool init);
#else
    template<class Vector1, class Matrix1>
    void Solve(const Matrix1& A, const Vector1& r, Vector1& z,
	       bool init_guess_null = true);

    template<class Vector1, class Matrix1>
    void TransSolve(const Matrix1& A, const Vector1& r, Vector1& z,
		    bool init_guess_null = true);
#endif

  };

}

#define SELDON_FILE_PRECOND_SSOR_HXX
#endif
