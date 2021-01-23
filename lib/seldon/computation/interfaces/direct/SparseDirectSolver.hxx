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

#ifndef SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_HXX

namespace Seldon
{
  
  //! Class grouping different direct solvers
  template<class T>
  class SparseDirectSolver
  {
  protected :
    //! ordering to use
    int type_ordering;
    //! solver to use
    int type_solver;
    //! number of threads (for Pastix)
    int nb_threads_per_node;
    //! ordering (if supplied by the user)
    IVect permut;
    //! size of factorized linear system 
    int n;
    //! pointer to the used solver
    VirtualSparseDirectSolver<T>* solver;
    //! threshold for ilut solver
    double threshold_matrix;
    double pivot_threshold;
    bool refine_solution;
    int print_level;
    //! use of non-symmetric ilut ?
    bool enforce_unsym_ilut;
        
  public :
    // available solvers
    enum {SELDON_SOLVER, UMFPACK, SUPERLU, MUMPS, PASTIX, ILUT, PARDISO, WSMP};
    
    // error codes
    enum {FACTO_OK, STRUCTURALLY_SINGULAR_MATRIX,
          NUMERICALLY_SINGULAR_MATRIX, OUT_OF_MEMORY, INVALID_ARGUMENT,
          INCORRECT_NUMBER_OF_ROWS, MATRIX_INDICES_INCORRECT,
          INVALID_PERMUTATION, ORDERING_FAILED, INTERNAL_ERROR,
          OVERFLOW_32BIT};
    
    SparseDirectSolver();
    
    // Inline methods
    int GetM() const;
    int GetN() const;
    
    int GetTypeOrdering() const;
    void SetPermutation(const IVect&);
    void SelectOrdering(int);
    
    void HideMessages();
    void ShowMessages();
    void ShowFullHistory();

    void SetPivotThreshold(const double&);
    void SetNumberOfThreadPerNode(int m);
    int GetNumberOfThreadPerNode() const;
    
    void SelectDirectSolver(int);    
    void SetNonSymmetricIlut();

    int GetDirectSolver();
    
    void RefineSolution();
    void DoNotRefineSolution();

    void SetCoefficientEstimationNeededMemory(double);
    void SetMaximumCoefficientEstimationNeededMemory(double);
    void SetIncreaseCoefficientEstimationNeededMemory(double);
    
    double GetThresholdMatrix() const;

    // other methods
    ~SparseDirectSolver();

    void Clear();

    static bool IsAvailableSolver(int type);
    bool AffectOrdering();

    template<class T0, class Prop, class Storage, class Alloc>
    void ComputeOrdering(Matrix<T0, Prop, Storage, Alloc>& A);
    
  protected:
    void InitSolver();
    
  public:
    void SetThresholdMatrix(const double&);
    
    template<class Prop, class Storage, class Allocator>
    void Factorize(Matrix<T, Prop, Storage, Allocator>& A,
		   bool keep_matrix = false);

    template<class Prop, class Storage, class Allocator>
    void PerformAnalysis(Matrix<T, Prop, Storage, Allocator>& A);
    
    template<class Prop, class Storage, class Allocator>
    void PerformFactorization(Matrix<T, Prop, Storage, Allocator>& A);
    
    int GetInfoFactorization(int& ierr) const;
    
    void Solve(Vector<T>& x);
    void Solve(const SeldonTranspose& TransA, Vector<T>& x, bool assemble = true);
    
    void Solve(Matrix<T, General, ColMajor>& x);
    void Solve(const SeldonTranspose&, Matrix<T, General, ColMajor>& x);

#ifdef SELDON_WITH_MPI
    template<class Tint0, class Tint1>
    void FactorizeDistributed(MPI_Comm& comm_facto,
                              Vector<Tint0>& Ptr, Vector<Tint1>& IndRow,
                              Vector<T>& Val, const IVect& glob_num,
                              bool sym, bool keep_matrix = false);

    template<class Tint0, class Tint1>
    void PerformAnalysisDistributed(MPI_Comm& comm_facto,
				    Vector<Tint0>& Ptr, Vector<Tint1>& IndRow,
				    Vector<T>& Val, const IVect& glob_num,
				    bool sym, bool keep_matrix = false);

    template<class Tint0, class Tint1>
    void PerformFactorizationDistributed(MPI_Comm& comm_facto,
					 Vector<Tint0>& Ptr, Vector<Tint1>& IndRow,
					 Vector<T>& Val, const IVect& glob_num,
					 bool sym, bool keep_matrix = false);
    
    void SolveDistributed(MPI_Comm& comm_facto, Vector<T>& x_solution,
                          const IVect& glob_number);
    
    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA, Vector<T>& x_solution,
                          const IVect& glob_number);

    void SolveDistributed(MPI_Comm& comm_facto,
			  Matrix<T, General, ColMajor>& x_solution,
                          const IVect& glob_number);
    
    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA,
			  Matrix<T, General, ColMajor>& x_solution,
                          const IVect& glob_number);
#endif
    
  };

#ifdef SELDON_WITH_MPI
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
			   SparseDirectSolver<T>& mat_lu,
                           Vector<T>& x, Vector<int>& global_col);
  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
                           SparseDirectSolver<complex<T> >& mat_lu,
                           Vector<T>& x, Vector<int>& global_col);
  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
			   SparseDirectSolver<T>& mat_lu,
                           Vector<complex<T> >& x, Vector<int>& global_col);
  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
			   SparseDirectSolver<T>& mat_lu,
                           Matrix<T, General, ColMajor>& x, Vector<int>& global_col);
  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
                           SparseDirectSolver<complex<T> >& mat_lu,
                           Matrix<T, General, ColMajor>& x, Vector<int>& global_col);
  
  template<class T>
  void SolveLU_Distributed(MPI_Comm& comm, const SeldonTranspose& transA,
			   SparseDirectSolver<T>& mat_lu,
                           Matrix<complex<T>, General, ColMajor>& x,
			   Vector<int>& global_col);
#endif  
  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<T>& mat_lu, Vector<T>& x);
  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<complex<T> >& mat_lu, Vector<T>& x);
  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<T>& mat_lu, Vector<complex<T> >& x);

  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<T>& mat_lu,
	       Matrix<T, General, ColMajor>& x);
  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<complex<T> >& mat_lu,
	       Matrix<T, General, ColMajor>& x);
  
  template<class T>
  void SolveLU(const SeldonTranspose& transA,
	       SparseDirectSolver<T>& mat_lu,
	       Matrix<complex<T>, General, ColMajor>& x);

  // SparseSolve and GetAndSolveLU
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                   Vector<T1, Storage1, Allocator1>& Y);


  template <class T, class Prop0, class Allocator0, class Allocator1>
  void GetAndSolveLU(Matrix<T, Prop0, ColSparse, Allocator0>& M,
		     Vector<T, VectFull, Allocator1>& Y);
  
  template <class T, class Prop0, class Allocator0, class Allocator1>
  void GetAndSolveLU(Matrix<T, Prop0, RowSparse, Allocator0>& M,
		     Vector<T, VectFull, Allocator1>& Y);

  template <class T, class Prop0, class Allocator0, class Allocator1>
  void GetAndSolveLU(Matrix<T, Prop0, ColSymSparse, Allocator0>& M,
		     Vector<T, VectFull, Allocator1>& Y);

  template <class T, class Prop0, class Allocator0, class Allocator1>
  void GetAndSolveLU(Matrix<T, Prop0, RowSymSparse, Allocator0>& M,
		     Vector<T, VectFull, Allocator1>& Y);
   
}

#define SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_HXX
#endif
