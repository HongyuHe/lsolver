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

#ifndef SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_CXX

#include "SparseDirectSolver.hxx"

namespace Seldon
{

  //! Default constructor
  template<class T>
  SparseDirectSolver<T>::SparseDirectSolver()
  {
    n = 0;
    type_ordering = SparseMatrixOrdering::AUTO;

    type_solver = SELDON_SOLVER;
    
    // we try to use an other available solver
    // The order of preference is Pastix, Mumps, Pardiso, UmfPack and SuperLU
#ifdef SELDON_WITH_SUPERLU
    type_solver = SUPERLU;
#endif
#ifdef SELDON_WITH_SUPERLU
    type_solver = SUPERLU;
#endif
#ifdef SELDON_WITH_UMFPACK
    type_solver = UMFPACK;
#endif
#ifdef SELDON_WITH_PARDISO
    type_solver = PARDISO;
#endif
#ifdef SELDON_WITH_WSMP
    type_solver = WSMP;
#endif
#ifdef SELDON_WITH_MUMPS
    type_solver = MUMPS;
#endif
#ifdef SELDON_WITH_PASTIX
    type_solver = PASTIX;
#endif
    
    nb_threads_per_node = 1;
    threshold_matrix = 0;
    print_level = 0;
    refine_solution = false;
    enforce_unsym_ilut = false;
    // for this default value, the corresponding method of the
    // sparse solver will not be called
    pivot_threshold = -2.0;
    
    solver = NULL;
    InitSolver();
  }
  
  
  //! Destructor
  template<class T>
  SparseDirectSolver<T>::~SparseDirectSolver()
  {
    if (solver != NULL)
      delete solver;
  }


  //! clearing factorisation
  template<class T>
  void SparseDirectSolver<T>::Clear()
  {
    if (n > 0)
      {
	solver->Clear();
	n = 0;    
      }    
  }
    
  
  //! returns true if the solver type is available
  template<class T> 
  bool SparseDirectSolver<T>::IsAvailableSolver(int type)
  {
    switch (type)
      {
      case SELDON_SOLVER: return true;
      case UMFPACK:
#ifdef SELDON_WITH_UMFPACK
	return true;
#endif
	return false;
      case SUPERLU:
#ifdef SELDON_WITH_SUPERLU
	return true;
#endif
	return false;
      case MUMPS:
#ifdef SELDON_WITH_MUMPS
	return true;
#endif
	return false;
      case PASTIX:
#ifdef SELDON_WITH_PASTIX
	return true;
#endif
	return false;
      case ILUT:
#ifdef SELDON_WITH_PRECONDITIONING
	return true;
#endif
	return false;
      case PARDISO:
#ifdef SELDON_WITH_PARDISO
	return true;
#endif
	return false;
      case WSMP:
#ifdef SELDON_WITH_WSMP
	return true;
#endif
	return false;
      default:
	return false;
      }
  }


  template<class T> 
  bool SparseDirectSolver<T>::AffectOrdering()
  {
    bool user_ordering = false;
    // we set the ordering for each direct solver interfaced
    switch (type_ordering)
      {
      case SparseMatrixOrdering::AUTO :
	{
	  // we choose the default strategy
	  // proposed by the direct solver that will be called
	  if (type_solver == MUMPS)
	    {
	      solver->SelectOrdering(7);
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_PASTIX
	      solver->SelectOrdering(API_ORDER_SCOTCH);
#endif
	    }
	  else if (type_solver == PARDISO)
	    {
	      solver->SelectOrdering(2);
	    }
	  else if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_UMFPACK
	      solver->SelectOrdering(UMFPACK_ORDERING_AMD);
#endif
	    }
	  else if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      solver->SelectOrdering(superlu::COLAMD);
#endif 
	    }
	  else if (type_solver == WSMP)
	    {
	    }
	  else
	    {
              
              type_ordering = SparseMatrixOrdering::IDENTITY;
              
#ifdef SELDON_WITH_UMFPACK
              type_ordering = SparseMatrixOrdering::AMD;
#endif

#ifdef SELDON_WITH_MUMPS
              type_ordering = SparseMatrixOrdering::METIS;
#endif

#ifdef SELDON_WITH_PARDISO
              type_ordering = SparseMatrixOrdering::METIS;
#endif

#ifdef SELDON_WITH_PASTIX
              type_ordering = SparseMatrixOrdering::SCOTCH;
#endif
              
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::IDENTITY :
      case SparseMatrixOrdering::REVERSE_CUTHILL_MCKEE :
      case SparseMatrixOrdering::USER :
	{
	  user_ordering = true;
	}
	break;
      case SparseMatrixOrdering::PORD :
      case SparseMatrixOrdering::AMF :
      case SparseMatrixOrdering::QAMD :
	{
	  // Mumps orderings
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
	      if (type_ordering == SparseMatrixOrdering::PORD)
		solver->SelectOrdering(4);
	      else if (type_ordering == SparseMatrixOrdering::AMF)
		solver->SelectOrdering(2);
	      else
		solver->SelectOrdering(6);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::SCOTCH :
      case SparseMatrixOrdering::PTSCOTCH :
	{
	  // available for Mumps and Pastix
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
              if (type_ordering==SparseMatrixOrdering::PTSCOTCH)
                solver->SelectParallelOrdering(1);
              else
                solver->SelectOrdering(3);
#endif
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_SCOTCH
	      solver->SelectOrdering(API_ORDER_SCOTCH);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::METIS :
      case SparseMatrixOrdering::PARMETIS :
	{
	  // available for Mumps and Pastix
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
              if (type_ordering==SparseMatrixOrdering::PARMETIS)
                solver->SelectParallelOrdering(2);
              else
                solver->SelectOrdering(5);
#endif
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_PASTIX
	      solver->SelectOrdering(API_ORDER_METIS);
#endif
	    }
	  else if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_UMFPACK
	      solver->SelectOrdering(UMFPACK_ORDERING_METIS);
#endif
	    }
          /*
            currently not implemented in SuperLU
	  else if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
              if (type_ordering==SparseMatrixOrdering::PARMETIS)
                solver->SelectOrdering(superlu::PARMETIS);
              else
                solver->SelectOrdering(superlu::METIS_AT_PLUS_A);
#endif
	    }
          */
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::AMD :
      case SparseMatrixOrdering::COLAMD :
	{
	  // default ordering for UmfPack
	  if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_UMFPACK
	      solver->SelectOrdering(UMFPACK_ORDERING_AMD);
#endif
	    }
	  else if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
              solver->SelectOrdering(0);
#endif
            }
	  else if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      solver->SelectOrdering(superlu::COLAMD);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
        break;
      case SparseMatrixOrdering::MMD_AT_PLUS_A:
        {
          if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      solver->SelectOrdering(superlu::MMD_AT_PLUS_A);
#endif
	    }
        }
        break;
      case SparseMatrixOrdering::MMD_ATA:
        {
          if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      solver->SelectOrdering(superlu::MMD_ATA);
#endif
	    }
        }
	break;
      }
    
    return user_ordering;
  }
  

  //! computation of the permutation vector in order to reduce fill-in
  template<class T> template<class T0, class Prop, class Storage, class Alloc>
  void SparseDirectSolver<T>
  ::ComputeOrdering(Matrix<T0, Prop, Storage, Alloc>& A)
  {
    bool user_ordering = AffectOrdering();
    
    if (user_ordering)
      {
	// case where the ordering is not natively available in the direct solver
	// computing separetely the ordering
	FindSparseOrdering(A, permut, type_ordering);
        
	solver->SetPermutation(permut);
      }
	
  }

  
  //! modifies threshold used for ilut
  template<class T>
  void SparseDirectSolver<T>::SetThresholdMatrix(const double& eps)
  {
    threshold_matrix = eps;
#ifdef SELDON_WITH_ILUT
    if (type_solver == ILUT)
      dynamic_cast<IlutPreconditioning<T>& >(*solver)
	.SetDroppingThreshold(eps);
#endif
  }

    
  //! initializes the solver (internal method)
  template<class T>
  void SparseDirectSolver<T>::InitSolver()
  {
    if (solver != NULL)
      delete solver;
    
    switch (type_solver)
      {
      case SUPERLU:
#ifdef SELDON_WITH_SUPERLU
	solver = new MatrixSuperLU<T>();
#else
	cout << "Seldon not compiled with SuperLU" << endl;
	cout << "SELDON_WITH_SUPERLU is not defined" << endl;
	abort();
#endif
	break;
      case UMFPACK:
#ifdef SELDON_WITH_UMFPACK
	solver = new MatrixUmfPack<T>();
#else
	cout << "Seldon not compiled with UmfPack" << endl;
	cout << "SELDON_WITH_UMFPACK is not defined" << endl;
	abort();
#endif
	break;
      case PARDISO:
#ifdef SELDON_WITH_PARDISO
	solver = new MatrixPardiso<T>();
#else
	cout << "Seldon not compiled with Pardiso" << endl;
	cout << "SELDON_WITH_PARDISO is not defined" << endl;
	abort();
#endif
	break;
      case MUMPS:
#ifdef SELDON_WITH_MUMPS
	solver = new MatrixMumps<T>();
#else
	cout << "Seldon not compiled with Mumps" << endl;
	cout << "SELDON_WITH_MUMPS is not defined" << endl;
	abort();
#endif
	break;
      case PASTIX:
#ifdef SELDON_WITH_PASTIX
	solver = new MatrixPastix<T>();
#else
	cout << "Seldon not compiled with Pastix" << endl;
	cout << "SELDON_WITH_PASTIX is not defined" << endl;
	abort();
#endif
	break;
      case WSMP:
#ifdef SELDON_WITH_WSMP
	solver = new MatrixWsmp<T>();
#else
	cout << "Seldon not compiled with Wsmp" << endl;
	cout << "SELDON_WITH_WSMP is not defined" << endl;
	abort();
#endif
	break;
      case ILUT:
	{
#ifdef SELDON_WITH_PRECONDITIONING
	  IlutPreconditioning<T>* ilut_solver = new IlutPreconditioning<T>();
	  ilut_solver->SetDroppingThreshold(threshold_matrix);
	  solver = ilut_solver;
#else
	  cout << "Seldon not compiled with Preconditioning" << endl;
	  cout << "SELDON_WITH_PRECONDITIONING is not defined" << endl;
	  abort();
#endif
	}
	break;
      case SELDON_SOLVER:
	solver = new SparseSeldonSolver<T>();
	break;
      default:
	cout << "Unknown solver" << endl;
	abort();
      }
    
    if (pivot_threshold != -2.0)
      solver->SetPivotThreshold(pivot_threshold);
    
    if (refine_solution)
      solver->RefineSolution();
    else
      solver->DoNotRefineSolution();

    solver->SetNumberOfThreadPerNode(nb_threads_per_node);
    if (print_level > 0)
      solver->ShowMessages();
    else
      solver->HideMessages();
  }
  
  
  //! factorisation of matrix A
  /*!
    LU factorisation is stored in the current object.
    You can ask to clear the matrix given on input (to spare memory)
   */
  template<class T>
  template<class Prop, class Storage, class Allocator>
  void SparseDirectSolver<T>::Factorize(Matrix<T, Prop, Storage, Allocator>& A,
					bool keep_matrix)
  {
    ComputeOrdering(A);
    n = A.GetM();
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
	MatrixUmfPack<T>& mat_umf =
	  dynamic_cast<MatrixUmfPack<T>& >(*solver);
	
	GetLU(A, mat_umf, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with UmfPack support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
	MatrixSuperLU<T>& mat_superlu =
	  dynamic_cast<MatrixSuperLU<T>& >(*solver);

	GetLU(A, mat_superlu, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else if (type_solver == PARDISO)
      {
#ifdef SELDON_WITH_PARDISO
	MatrixPardiso<T>& mat_pardiso =
	  dynamic_cast<MatrixPardiso<T>& >(*solver);

	GetLU(A, mat_pardiso, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Pardiso support.");
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	MatrixMumps<T>& mat_mumps =
	  dynamic_cast<MatrixMumps<T>& >(*solver);
	
	GetLU(A, mat_mumps, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	MatrixPastix<T>& mat_pastix =
	  dynamic_cast<MatrixPastix<T>& >(*solver);

	GetLU(A, mat_pastix, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == WSMP)
      {
#ifdef SELDON_WITH_WSMP
	MatrixWsmp<T>& mat_wsmp =
	  dynamic_cast<MatrixWsmp<T>& >(*solver);

	GetLU(A, mat_wsmp, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Wsmp support.");
#endif
      }
    else if (type_solver == ILUT)
      {
#ifdef SELDON_WITH_PRECONDITIONING        
	IlutPreconditioning<T>& mat_ilut =
	  dynamic_cast<IlutPreconditioning<T>& >(*solver);

        // setting some parameters
        if (enforce_unsym_ilut || (!IsSymmetricMatrix(A)))
          mat_ilut.SetUnsymmetricAlgorithm();
        else
          mat_ilut.SetSymmetricAlgorithm();
        
        // then performing factorization
	GetLU(A, mat_ilut, permut, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with the preconditioners.");
#endif
      }
    else
      {
	SparseSeldonSolver<T>& mat_seldon =
	  static_cast<SparseSeldonSolver<T>& >(*solver);

	GetLU(A, mat_seldon, permut, keep_matrix);
      }
  }
   

  //! Performs the analysis of the matrix before numerical factorization
  template <class T> template<class Prop, class Storage, class Allocator>
  void SparseDirectSolver<T>::PerformAnalysis(Matrix<T, Prop, Storage, Allocator>& A)
  {
    ComputeOrdering(A);
    n = A.GetM();
    if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	MatrixMumps<T>& mat_mumps =
	  dynamic_cast<MatrixMumps<T>& >(*solver);
	
	mat_mumps.PerformAnalysis(A);
#else
        throw Undefined("SparseDirectSolver::PerformAnalysis(MatrixSparse&, bool)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else
      {
	cout << "Not implemented for other solvers" << endl;
	abort();
      }
  }
    
  
  //! Performs the numerical factorization
  template <class T> template<class Prop, class Storage, class Allocator>
  void SparseDirectSolver<T>::PerformFactorization(Matrix<T, Prop, Storage, Allocator>& A)
  {
    if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	MatrixMumps<T>& mat_mumps =
	  dynamic_cast<MatrixMumps<T>& >(*solver);
	
	mat_mumps.PerformFactorization(A);
#else
        throw Undefined("SparseDirectSolver::PerformAnalysis(MatrixSparse&, bool)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else
      {
	cout << "Not implemented for other solvers" << endl;
	abort();
      }
  }


  //! Returns error code of the direct solver
  template <class T>
  int SparseDirectSolver<T>::GetInfoFactorization(int& ierr) const
  {
    ierr = solver->GetInfoFactorization();
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
        switch (ierr)
          {
          case UMFPACK_OK :
            return FACTO_OK;
          case UMFPACK_WARNING_singular_matrix :
            return NUMERICALLY_SINGULAR_MATRIX;
          case UMFPACK_ERROR_out_of_memory :
            return OUT_OF_MEMORY;
          case UMFPACK_ERROR_invalid_Numeric_object :
          case UMFPACK_ERROR_invalid_Symbolic_object :
          case UMFPACK_ERROR_argument_missing :
          case UMFPACK_ERROR_different_pattern :
          case UMFPACK_ERROR_invalid_system :
            return INVALID_ARGUMENT;
          case UMFPACK_ERROR_n_nonpositive :
            return INCORRECT_NUMBER_OF_ROWS;
          case UMFPACK_ERROR_invalid_matrix :
            return MATRIX_INDICES_INCORRECT;
          case UMFPACK_ERROR_invalid_permutation :
            return INVALID_PERMUTATION;
          case UMFPACK_ERROR_ordering_failed :
            return ORDERING_FAILED;
          default :
            return INTERNAL_ERROR;
          }
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
        if (ierr > 0)
          return INTERNAL_ERROR;
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
        switch (ierr)
          {
          case -2 :
            // nz out of range
            return MATRIX_INDICES_INCORRECT;
          case -3 :
            // invalid job number
            return INVALID_ARGUMENT;
          case -4 :
            // invalid permutation
            return INVALID_PERMUTATION;
          case -5 :
            // problem of real workspace allocation
            return OUT_OF_MEMORY;
          case -6 :
            // structurally singular matrix
            return STRUCTURALLY_SINGULAR_MATRIX;
          case -7 :
            // problem of integer workspace allocation
            return OUT_OF_MEMORY;
          case -10 :
            // numerically singular matrix
            return NUMERICALLY_SINGULAR_MATRIX;
          case -13 :
            // allocate failed
            return OUT_OF_MEMORY;
          case -16 :
            // N out of range
            return INCORRECT_NUMBER_OF_ROWS;
          case -22 :
            // invalid pointer
            return INVALID_ARGUMENT;
          case 1 :
            // index out of range
            return MATRIX_INDICES_INCORRECT;
	  case 0 :
	    return FACTO_OK;
          default :
            return INTERNAL_ERROR;
          }
#endif
      }
    else if (type_solver == PARDISO)
      {
#ifdef SELDON_WITH_PARDISO
        switch (ierr)
          {
          case -1 :
            return INVALID_ARGUMENT;
          case -2 :
          case -9 :
            return OUT_OF_MEMORY;
          case -3 :
            return INVALID_PERMUTATION;
          case -4 :
            return NUMERICALLY_SINGULAR_MATRIX;
          case -6 :
            return ORDERING_FAILED;
          case -7 :
            return STRUCTURALLY_SINGULAR_MATRIX;
          case -8 :
            return OVERFLOW_32BIT;            
	  case 0 :
	    return FACTO_OK;
          default :
            return INTERNAL_ERROR;
          }
#endif
      }
    
    return FACTO_OK;
  }
  
    
  //! x_solution is overwritten by solution of A x = b
  /*!
    We assume that Factorize has been called previously
  */
  template<class T>
  void SparseDirectSolver<T>::Solve(Vector<T>& x_solution)
  {
    Solve(SeldonNoTrans, x_solution);
  }
  
  
  //! x_solution is overwritten with solution of A x = b or A^T x = b
  template<class T>
  void SparseDirectSolver<T>::Solve(const SeldonTranspose& trans,
				    Vector<T>& x_solution, bool assemble)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (n != x_solution.GetM())
      throw WrongDim("SparseDirectSolver::Solve", 
		     string("The length of x is equal to ")
		     + to_str(x_solution.GetM())
		     + " while the size of the matrix is equal to "
		     + to_str(n) + ".");
#endif
    
    solver->Solve(trans, x_solution.GetData(), 1);
  }
  

  //! x_solution is overwritten by solution of A x = b
  /*!
    Multiple right hand sides
    We assume that Factorize has been called previously
  */
  template<class T>
  void SparseDirectSolver<T>
  ::Solve(Matrix<T, General, ColMajor>& x_solution)
  {
    Solve(SeldonNoTrans, x_solution);
  }


  //! x_solution is overwritten by solution of A x = b
  /*!
    Multiple right hand sides
    We assume that Factorize has been called previously
  */
  template<class T>
  void SparseDirectSolver<T>
  ::Solve(const SeldonTranspose& trans, Matrix<T, General, ColMajor>& x_sol)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (n != x_sol.GetM())
      throw WrongDim("SparseDirectSolver::Solve", 
		     string("The length of x is equal to ")
		     + to_str(x_sol.GetM())
		     + " while the size of the matrix is equal to "
		     + to_str(n) + ".");
#endif
    
    solver->Solve(trans, x_sol.GetData(), x_sol.GetN());
  }

  
#ifdef SELDON_WITH_MPI
  //! Factorization of a matrix
  /*!
    The matrix is given on each processor of the communicator in CSC
    format. If the matrix is assumed to be symmetric, you provide only
    the lower part of the matrix.
    \param[in] comm_facto communicator grouping processors involved in the
    factorisation.
    \param[inout] Ptr start indices
    \param[inout] Row column indices
    \param[inout] Val values of non-zero entries
    \param[in] glob_num global column numbers, each column of the global
    matrix is associated with one processor and only one
    \param[in] sym if true, the matrix is assumed to be symmetric
    \param[in] keep_matrix if false the input matrix is erased
   */
  template<class T> template<class Tint0, class Tint1>
  void SparseDirectSolver<T>::
  FactorizeDistributed(MPI_Comm& comm_facto,
                       Vector<Tint0>& Ptr, Vector<Tint1>& Row, Vector<T>& Val,
                       const IVect& glob_num, bool sym, bool keep_matrix)
  {
    bool user_ordering = AffectOrdering();
    if (user_ordering)
      {
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI_Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "user orderings not available in parallel");
      }
    
    n = Ptr.GetM()-1;
    solver->FactorizeDistributedMatrix(comm_facto, Ptr, Row, Val,
				       glob_num, sym, keep_matrix);
    
  }

  
  //! Performs analysis before numerical factorization
  template<class T> template<class Tint0, class Tint1>
  void SparseDirectSolver<T>::
  PerformAnalysisDistributed(MPI_Comm& comm_facto,
			     Vector<Tint0>& Ptr, Vector<Tint1>& Row, Vector<T>& Val,
			     const IVect& glob_num, bool sym, bool keep_matrix)
  {
    bool user_ordering = AffectOrdering();
    if (user_ordering)
      {
        throw Undefined("SparseDirectSolver::PerformAnalysisDistributed(MPI_Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "user orderings not available in parallel");
      }
    
    n = Ptr.GetM()-1;
    solver->PerformAnalysisDistributed(comm_facto, Ptr, Row, Val,
				       glob_num, sym, keep_matrix);
    
  }

  
  //! Performs numerical factorization
  template<class T> template<class Tint0, class Tint1>
  void SparseDirectSolver<T>::
  PerformFactorizationDistributed(MPI_Comm& comm_facto,
				  Vector<Tint0>& Ptr, Vector<Tint1>& Row, Vector<T>& Val,
				  const IVect& glob_num, bool sym, bool keep_matrix)
  {
    n = Ptr.GetM()-1;
    solver->PerformFactorizationDistributed(comm_facto, Ptr, Row, Val,
					    glob_num, sym, keep_matrix);    
  }
  
  
  //! solution of linear system Ax = b by using LU factorization
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T>
  void SparseDirectSolver<T>::
  SolveDistributed(MPI_Comm& comm_facto, Vector<T>& x_solution,
		   const IVect& glob_number)
  {
    SolveDistributed(comm_facto, SeldonNoTrans, x_solution, glob_number);
  }


  //! solution of linear system Ax = b by using LU factorization
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T>
  void SparseDirectSolver<T>::
  SolveDistributed(MPI_Comm& comm_facto, const SeldonTranspose& trans,
		   Vector<T>& x_solution, const IVect& glob_number)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (n != x_solution.GetM())
      throw WrongDim("SparseDirectSolver::SolveDistributed", 
		     string("The length of x is equal to ")
		     + to_str(x_solution.GetM())
		     + " while the size of the matrix is equal to "
		     + to_str(n) + ".");
#endif

    solver->SolveDistributed(comm_facto, trans,
			     x_solution.GetData(), 1, glob_number);
  }


  //! solution of linear system Ax = b by using LU factorization
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T>
  void SparseDirectSolver<T>::
  SolveDistributed(MPI_Comm& comm_facto, Matrix<T, General, ColMajor>& x,
		   const IVect& glob_number)
  {
    SolveDistributed(comm_facto, SeldonNoTrans, x, glob_number);
  }


  //! solution of linear system Ax = b by using LU factorization
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T>
  void SparseDirectSolver<T>::
  SolveDistributed(MPI_Comm& comm_facto, const SeldonTranspose& trans,
		   Matrix<T, General, ColMajor>& x, const IVect& glob_number)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (n != x.GetM())
      throw WrongDim("SparseDirectSolver::SolveDistributed", 
		     string("The length of x is equal to ")
		     + to_str(x.GetM())
		     + " while the size of the matrix is equal to "
		     + to_str(n) + ".");
#endif

    solver->SolveDistributed(comm_facto, trans,
			     x.GetData(), x.GetN(), glob_number);
  }


  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
			   SparseDirectSolver<T>& mat_lu,
                           Vector<T>& x, Vector<int>& global_col)
  {
    mat_lu.SolveDistributed(comm, transA, x, global_col);
  }

  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
                           SparseDirectSolver<complex<T> >& mat_lu,
                           Vector<T>& x, Vector<int>& global_col)
  {
    throw WrongArgument("SolveLU(SparseDirectSolver, Vector, Vector)",
			"the result must be complex");
  }
  
  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
			   SparseDirectSolver<T>& mat_lu,
                           Vector<complex<T> >& x, Vector<int>& global_col)
  {
    Vector<T> xreal(x.GetM()), ximag(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
        xreal(i) = real(x(i));
        ximag(i) = imag(x(i));
      }
    
    mat_lu.SolveDistributed(comm, transA, xreal, global_col);
    mat_lu.SolveDistributed(comm, transA, ximag, global_col);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T>(xreal(i), ximag(i));
  }
  
  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
			   SparseDirectSolver<T>& mat_lu,
                           Matrix<T, General, ColMajor>& x, Vector<int>& global_col)
  {
    mat_lu.SolveDistributed(comm, transA, x, global_col);
  }

  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
                           SparseDirectSolver<complex<T> >& mat_lu,
                           Matrix<T, General, ColMajor>& x, Vector<int>& global_col)
  {
    throw WrongArgument("SolveLU(SparseDirectSolver, Vector, Vector)",
			"the result must be complex");
  }
  
  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
			   SparseDirectSolver<T>& mat_lu,
                           Matrix<complex<T>, General, ColMajor>& x,
			   Vector<int>& global_col)
  {
    throw WrongArgument("SolveLU(SparseDirectSolver, Vector, Vector)",
			"This case is currently not handled");
  }
#endif  
  

  // Solve LU with vector
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<T>& mat_lu, Vector<T>& x)
  {
    mat_lu.Solve(transA, x);
  }

  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<complex<T> >& mat_lu, Vector<T>& x)
  {
    throw WrongArgument("SolveLU(SparseDirectSolver, Vector, Vector)",
			"the result must be complex");
  }
  
  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<T>& mat_lu, Vector<complex<T> >& x)
  {
    // using a matrix with two columns
    Matrix<T, General, ColMajor> xm(x.GetM(), 2);
    for (int i = 0; i < x.GetM(); i++)
      {
        xm(i, 0) = real(x(i));
        xm(i, 1) = imag(x(i));
      }
    
    mat_lu.Solve(transA, xm);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T>(xm(i, 0), xm(i, 1));
  }
  

  // Solve LU with matrix  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<T>& mat_lu,
	       Matrix<T, General, ColMajor>& x)
  {
    mat_lu.Solve(transA, x);
  }

  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<complex<T> >& mat_lu,
	       Matrix<T, General, ColMajor>& x)
  {
    throw WrongArgument("SolveLU(SparseDirectSolver, Vector, Vector)",
			"the result must be complex");
  }
  
  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<T>& mat_lu,
	       Matrix<complex<T>, General, ColMajor>& x)
  {
    throw WrongArgument("SolveLU(SparseDirectSolver, Vector, Vector)",
			"This case is currently not handled");
  }
  
}

#define SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_CXX
#endif

