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


#ifndef SELDON_FILE_PRECOND_SSOR_CXX

namespace Seldon
{

  //! Default constructor
  template<class T>
  SorPreconditioner<T>::SorPreconditioner()
  {
    nb_iter = 1; omega = 1;
    symmetric_precond = true;
  }

  
  //! returns true if symmetric sor is used
  template<class T>
  bool SorPreconditioner<T>::IsSymmetric() const
  {
    return symmetric_precond;
  }
  
  
  //! if called forward and backward sweep will be applied such that
  //! the preconditioning is symmetric
  template<class T>
  void SorPreconditioner<T>::InitSymmetricPreconditioning()
  {
    symmetric_precond = true;
  }
  
  
  //! if called, forward sweep is applied when calling Solve
  //! backward sweep is applied when calling TransSolve
  template<class T>
  void SorPreconditioner<T>::InitUnSymmetricPreconditioning()
  {
    symmetric_precond = false;
  }
  
  
  //! sets the relaxation parameter omega
  template<class T>
  void SorPreconditioner<T>
  ::SetParameterRelaxation(const typename ClassComplexType<T>::Treal& param)
  {
    omega = param;
  }
  
  
  //! sets the number of SOR sweeps to perform when calling Solve/TransSolve
  template<class T>
  void SorPreconditioner<T>::SetNumberIterations(int nb_iterations)
  {
    nb_iter = nb_iterations;
  }


  //! copies the parameters of another sor preconditioning
  template<class T>
  void SorPreconditioner<T>::CopyParameter(const Preconditioner_Base<T>& prec)
  {
    const SorPreconditioner<T>& prec_sor = dynamic_cast<const SorPreconditioner<T>& >(prec);
    symmetric_precond = prec_sor.symmetric_precond;
    nb_iter = prec_sor.nb_iter;
    omega = prec_sor.omega;
  }
  

#ifdef SELDON_WITH_VIRTUAL
  //! applies preconditioning for real or complex vectors
  template<class T> template<class T0>
  void SorPreconditioner<T>
  ::SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
	     const Vector<T0>& r, Vector<T0>& z, bool init)
  {
    if (init)
      z.Zero();
    
    int stage = 0;
    if (!symmetric_precond)
      {
	if (trans.NoTrans())
	  stage = 2;
	else
	  stage = 3;
      }
    
    A.ApplySor(trans, z, r, omega, nb_iter, stage);
  }

  
  //! performs an SOR iteration
  template<class T>
  void SorPreconditioner<T>
  ::Solve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z, bool init)
  {
    if (init)
      z.Zero();
   
    if (symmetric_precond)
      A.ApplySor(SeldonNoTrans, z, r, omega, nb_iter, 0);
    else
      A.ApplySor(SeldonNoTrans, z, r, omega, nb_iter, 2);
  
  }


  //! performs an SOR iteration (transpose)
  template<class T>
  void SorPreconditioner<T>
  ::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z, bool init)
  {
    if (init)
      z.Zero();
    
    if (symmetric_precond)
      A.ApplySor(SeldonTrans, z, r, omega, nb_iter, 0);
    else
      A.ApplySor(SeldonTrans, z, r, omega, nb_iter, 3);
  }


  //! applies SOR preconditioning
  template<class T>
  void SorPreconditioner<T>
  ::Solve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z)
  {
    Solve(A, r, z, true);
  }


  //! applies transpose SOR preconditioning
  template<class T>
  void SorPreconditioner<T>
  ::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z)
  {
    TransSolve(A, r, z, true);
  }
#else

  //! Solves M z = r
  template<class T> template<class Vector1, class Matrix1>
  void SorPreconditioner<T>::
  Solve(const Matrix1& A, const Vector1& r, Vector1& z, bool init_guess_null)
  {
    if (init_guess_null)
      z.Zero();
    
    if (symmetric_precond)
      SOR(A, z, r, omega, nb_iter, 0);
    else
      SOR(A, z, r, omega, nb_iter, 2);
  }


  //! Solves M^t z = r
  template<class T> template<class Vector1, class Matrix1>
  void SorPreconditioner<T>::
  TransSolve(const Matrix1& A, const Vector1& r,
	     Vector1& z, bool init_guess_null)
  {
    if (init_guess_null)
      z.Zero();
    
    if (symmetric_precond)
      SOR(SeldonTrans, A, z, r, omega, nb_iter, 0);
    else
      SOR(SeldonTrans, A, z, r, omega, nb_iter, 3);
    
  }

#endif

}

#define SELDON_FILE_PRECOND_SSOR_CXX
#endif
