// Copyright (C) 2015 Marc Duruflé
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

#ifndef SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_INLINE_CXX

namespace Seldon
{

  //! returns the number of rows of the factorised matrix
  template<class T>
  inline int SparseDirectSolver<T>::GetM() const
  {
    return n;
  }
  
  
  //! returns the number of rows of the factorised matrix
  template<class T>
  inline int SparseDirectSolver<T>::GetN() const
  {
    return n;
  }
  
  
  //! returns the ordering algorithm to use
  template<class T>
  inline int SparseDirectSolver<T>::GetTypeOrdering() const
  {
    return type_ordering;
  }
  
  
  //! sets directly the new ordering (by giving a permutation vector)
  template<class T>
  inline void SparseDirectSolver<T>::SetPermutation(const IVect& num)
  {
    type_ordering = SparseMatrixOrdering::USER;
    permut = num;
  }
  
  
  //! modifies the ordering algorithm to use
  template<class T>
  inline void SparseDirectSolver<T>::SelectOrdering(int type)
  {
    type_ordering = type;
  }


  //! hiding all messages
  template<class T>
  inline void SparseDirectSolver<T>::HideMessages()
  {
    print_level = 0;
    solver->HideMessages();
  }
  
  
  //! displaying basic messages
  template<class T>
  inline void SparseDirectSolver<T>::ShowMessages()
  {
    print_level = 1;
    solver->ShowMessages();    
  }
  
  
  //! displaying all the messages
  template<class T>
  inline void SparseDirectSolver<T>::ShowFullHistory()
  {
    print_level = 2;
    solver->ShowMessages();        
  }
    
    
  //! sets the threshold used for pivoting
  template<class T>
  void SparseDirectSolver<T>::SetPivotThreshold(const double& eps)
  {
    pivot_threshold = eps;
    solver->SetPivotThreshold(eps);
  }
    
  
  //! modifies the number of threads per node (for Pastix only)
  template<class T>
  inline void SparseDirectSolver<T>::SetNumberOfThreadPerNode(int p)
  {
    nb_threads_per_node = p;
    solver->SetNumberOfThreadPerNode(nb_threads_per_node);
  }
  

  //! modifies the number of threads per node (for Pastix only)
  template<class T>
  inline int SparseDirectSolver<T>::GetNumberOfThreadPerNode() const
  {
    return nb_threads_per_node;
  }
  
  
  //! modifies the direct solver to use
  template<class T>
  inline void SparseDirectSolver<T>::SelectDirectSolver(int type)
  {
    type_solver = type;
    type_ordering = SparseMatrixOrdering::AUTO;
    InitSolver();
  }
  
  
  //! enforces the use of unsymmetric algorithm for ilut solver
  template<class T>
  inline void SparseDirectSolver<T>::SetNonSymmetricIlut()
  {
    enforce_unsym_ilut = true;
  }

  
  //! returns the direct solver to use
  template<class T>
  inline int SparseDirectSolver<T>::GetDirectSolver()
  {
    return type_solver;
  }
    
  
  template<class T>
  inline void SparseDirectSolver<T>::RefineSolution()
  {
    refine_solution = true;
    solver->RefineSolution();
  }
  
  
  template<class T>
  inline void SparseDirectSolver<T>::DoNotRefineSolution()
  {
    refine_solution = false;
    solver->DoNotRefineSolution();
  }
  
  
  template<class T>
  inline void SparseDirectSolver<T>::
  SetCoefficientEstimationNeededMemory(double coef)
  {
    solver->SetCoefficientEstimationNeededMemory(coef);
  }
  
  
  template<class T>
  inline void SparseDirectSolver<T>::
  SetMaximumCoefficientEstimationNeededMemory(double coef)
  {
    solver->SetMaximumCoefficientEstimationNeededMemory(coef);
  }
  
  
  template<class T>
  inline void SparseDirectSolver<T>::
  SetIncreaseCoefficientEstimationNeededMemory(double coef)
  {
    solver->SetIncreaseCoefficientEstimationNeededMemory(coef);
  }
  
  
  //! returns threshold used for ilut (if this solver is selected)
  template<class T>
  inline double SparseDirectSolver<T>::GetThresholdMatrix() const
  {
    return threshold_matrix;
  }

  
  /*************************
   * Solve and SparseSolve *
   *************************/
  
  
  //! Solves a sparse linear system using LU factorization.
  /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
    \f$ X \f$ and \f$ Y \f$ are vectors.
    \param[in] M the sparse matrix of the linear system, to be factorized in
    LU form by UMFPACK, SuperLU or Mumps. On exit, \a M is cleared.
    \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
    solution \f$ X \f$ of the system.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  inline void SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
			  Vector<T1, Storage1, Allocator1>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, Y, "SolveLU(M, Y)");
#endif
    
    typedef typename Matrix<T0, Prop0, Storage0, Allocator0>::entry_type Tmat;
    SparseDirectSolver<Tmat> matrix_lu;
    
    matrix_lu.Factorize(M);
    matrix_lu.Solve(Y);
  }


  //! \copydoc SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
  //! Vector<T1, Storage1, Allocator1>& Y)
  template <class T, class Prop0, class Allocator0, class Allocator1>
  inline void GetAndSolveLU(Matrix<T, Prop0, ColSparse, Allocator0>& M,
			    Vector<T, VectFull, Allocator1>& Y)
  {
    SparseSolve(M, Y);
  }


  //! \copydoc SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
  //! Vector<T1, Storage1, Allocator1>& Y)
  template <class T, class Prop0, class Allocator0, class Allocator1>
  inline void GetAndSolveLU(Matrix<T, Prop0, RowSparse, Allocator0>& M,
			    Vector<T, VectFull, Allocator1>& Y)
  {
    SparseSolve(M, Y);
  }


  //! \copydoc SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
  //! Vector<T1, Storage1, Allocator1>& Y)
  template <class T, class Prop0, class Allocator0, class Allocator1>
  inline void GetAndSolveLU(Matrix<T, Prop0, ColSymSparse, Allocator0>& M,
			    Vector<T, VectFull, Allocator1>& Y)
  {
    SparseSolve(M, Y);
  }


  //! \copydoc SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
  //! Vector<T1, Storage1, Allocator1>& Y)
  template <class T, class Prop0, class Allocator0, class Allocator1>
  inline void GetAndSolveLU(Matrix<T, Prop0, RowSymSparse, Allocator0>& M,
			    Vector<T, VectFull, Allocator1>& Y)
  {
    SparseSolve(M, Y);
  }
  
}

#define SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_INLINE_CXX
#endif

