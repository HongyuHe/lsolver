#ifndef SELDON_FILE_DISTRIBUTED_SOLVER_CXX

#include "DistributedSolver.hxx"

namespace Seldon
{

  //! Default constructor
  template <class T>
  SparseDistributedSolver<T>::SparseDistributedSolver()
    : SparseDirectSolver<T>()
  {
    diagonal_scaling_left = false;
    diagonal_scaling_right = false;
    
    this->refine_solution = true;

#ifdef SELDON_WITH_MPI
    // initialisation for a sequential matrix
    nodl_scalar_ = 1;
    nb_unknowns_scal_ = 1;
    comm_ = MPI_COMM_SELF;
    ProcSharingRows_ = NULL;
    SharingRowNumbers_ = NULL;
#endif
  }


  //! Changes the level of displayed messages
  template <class T>
  void SparseDistributedSolver<T>::SetPrintLevel(int print_level)
  { 
    if (print_level >= 2)
      this->ShowMessages(); 
    else if (print_level >= 6)
      this->ShowFullHistory();
    else
      this->HideMessages();    
  }
  
  
  //! Clears LU matrices
  template <class T>
  void SparseDistributedSolver<T>::Clear()
  {
    SparseDirectSolver<T>::Clear();
    
    diagonal_scale_left.Clear();
    diagonal_scale_right.Clear();
  }
  
    
  //! scales matrix with 1 / \sum a_{i, j} (rows and columns)
  template<class T> template<class MatrixSparse>
  void SparseDistributedSolver<T>::ScaleMatrixRowCol(MatrixSparse& A)
  {
    int n = A.GetM();
    diagonal_scaling_left = true;
    diagonal_scaling_right = true;
    
    GetRowColSum(diagonal_scale_left, diagonal_scale_right, A);
    
    // forming D^-1 A C^-1 where d_i = \sum_j | a_ij |
    // and c_j = \sum_i | a_ij |
    for (int i = 0; i < n; i++)
      {
        if (diagonal_scale_left(i) != 0)
          diagonal_scale_left(i) = abs(1.0/diagonal_scale_left(i));
        else
          diagonal_scale_left(i) = 1.0;
        
        if (diagonal_scale_right(i) != 0)
          diagonal_scale_right(i) = abs(1.0/diagonal_scale_right(i));
        else
          diagonal_scale_right(i) = 1.0;
      }
    
    ScaleMatrix(A, diagonal_scale_left, diagonal_scale_right);
  }
  
  
  //! factorization of matrix in sequential
  template<class T>
  template<class Prop, class Storage, class Allocator>
  void SparseDistributedSolver<T>
  ::Factorize(Matrix<T, Prop, Storage, Allocator>& A,
	      bool keep_matrix, bool scale_matrix)
  {
    diagonal_scaling_left = false;
    diagonal_scaling_right = false;
    if (scale_matrix)
      ScaleMatrixRowCol(A);
    
    // then factorizing the modified matrix
    SparseDirectSolver<T>::Factorize(A, keep_matrix);
  }


  //! Performs the analysis of the matrix A before numerical factorization
  template<class T>
  template<class Prop, class Storage, class Allocator>
  void SparseDistributedSolver<T>
  ::PerformAnalysis(Matrix<T, Prop, Storage, Allocator>& A)
  {
    SparseDirectSolver<T>::PerformAnalysis(A);
  }


  //! Performs the numerical factorization
  template<class T>
  template<class Prop, class Storage, class Allocator>
  void SparseDistributedSolver<T>
  ::PerformFactorization(Matrix<T, Prop, Storage, Allocator>& A,
			 bool scale_matrix)
  {
    diagonal_scaling_left = false;
    diagonal_scaling_right = false;
    if (scale_matrix)
      ScaleMatrixRowCol(A);
    
    // then factorizing the modified matrix
    SparseDirectSolver<T>::PerformFactorization(A);
  }
  

#ifdef SELDON_WITH_MPI  
  //! factorization of distributed matrix
  template<class T>
  template<class Prop0, class Storage0, class Allocator0>
  void SparseDistributedSolver<T>
  ::Factorize(DistributedMatrix<T, Prop0, Storage0, Allocator0>& A,
	      bool keep_matrix, bool scale_matrix)
  {
    diagonal_scaling_left = false;
    diagonal_scaling_right = false;
    if (scale_matrix)
      ScaleMatrixRowCol(A);

    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        comm_ = comm;
        SparseDirectSolver<T>::Factorize(A, keep_matrix);
      }
    else
      {
        // factorisation of a distributed matrix
	this->n = A.GetM();
        comm_ = comm;
        ProcSharingRows_ = &A.GetProcessorSharingRows();
        SharingRowNumbers_ = &A.GetSharingRowNumbers();
        nodl_scalar_ = A.GetNodlScalar();
        nb_unknowns_scal_ = A.GetNbScalarUnknowns();
        bool sym_matrix = IsSymmetricMatrix(A);
	
	DistributedMatrix<T, Prop0, Storage0, Allocator0> Bstore;
	DistributedMatrix<T, Prop0, Storage0, Allocator0>* B;
	B = &A;
	if (keep_matrix)
	  {
	    B = &Bstore;
	    Bstore = A;
	  }

        // assembles distributed matrix
        int rank_proc; MPI_Comm_rank(comm, &rank_proc);
        if ((this->print_level >= 1) &&(rank_proc == 0))
          cout << "Assembling the distributed matrix..." << endl;
        
	if (this->solver->UseInteger8())
	  {
	    Vector<int64_t> Ptr, IndRow;
	    Vector<T> Val;
	    bool sym_pattern = false;
	    if (this->type_solver == this->PASTIX)
	      sym_pattern = true;
	    
	    bool reorder_num = false;
	    if (this->type_solver == this->SUPERLU)
	      {
		General sym;
		reorder_num = true;
		AssembleDistributed(*B, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }
	    else
	      {
		Prop0 sym;
		AssembleDistributed(*B, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }

	    if ((this->print_level >= 1) && (rank_proc == 0))
	      cout << "Factorizing the distributed matrix..." << endl;
	    
	    // factorizes the matrix
	    this->FactorizeDistributed(comm, Ptr, IndRow, Val,
				       global_col_numbers, sym_matrix, reorder_num);
	  }
	else
	  {
	    Vector<long> Ptr; Vector<int> IndRow;
	    Vector<T> Val;
	    bool sym_pattern = false;
	    if (this->type_solver == this->PASTIX)
	      sym_pattern = true;
	    
	    bool reorder_num = false;
	    if (this->type_solver == this->SUPERLU)
	      {
		General sym;
		reorder_num = true;
		AssembleDistributed(*B, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }
	    else
	      {
		Prop0 sym;
		AssembleDistributed(*B, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }

	    if ((this->print_level >= 1) && (rank_proc == 0))
	      cout << "Factorizing the distributed matrix..." << endl;
	    
	    // factorizes the matrix
	    this->FactorizeDistributed(comm, Ptr, IndRow, Val,
				       global_col_numbers, sym_matrix, reorder_num);
	  }        
      }
  }


  //! Performs the analysis of the matrix before numerical factorization
  template<class T>
  template<class Prop0, class Storage0, class Allocator0>
  void SparseDistributedSolver<T>
  ::PerformAnalysis(DistributedMatrix<T, Prop0, Storage0, Allocator0>& A)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        comm_ = comm;
        SparseDirectSolver<T>::PerformAnalysis(A);
      }
    else
      {
        // factorisation of a distributed matrix
	this->n = A.GetM();
        comm_ = comm;
        ProcSharingRows_ = &A.GetProcessorSharingRows();
        SharingRowNumbers_ = &A.GetSharingRowNumbers();
        nodl_scalar_ = A.GetNodlScalar();
        nb_unknowns_scal_ = A.GetNbScalarUnknowns();
        bool sym_matrix = IsSymmetricMatrix(A);
        
        // assembles distributed matrix
        int rank_proc; MPI_Comm_rank(comm, &rank_proc);
        if ((this->print_level >= 1) &&(rank_proc == 0))
          cout << "Assembling the distributed matrix..." << endl;
        
	DistributedMatrix<T, Prop0, Storage0, Allocator0> B; B = A;
	if (this->solver->UseInteger8())
	  {
	    Vector<int64_t> Ptr, IndRow;
	    Vector<T> Val;
	    bool sym_pattern = false;
	    if (this->type_solver == this->PASTIX)
	      sym_pattern = true;
	    
	    bool reorder_num = false;
	    if (this->type_solver == this->SUPERLU)
	      {
		General sym;
		reorder_num = true;
		AssembleDistributed(B, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }
	    else
	      {
		Prop0 sym;
		AssembleDistributed(B, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }

	    if ((this->print_level >= 1) && (rank_proc == 0))
	      cout << "Analysing the distributed matrix..." << endl;
	    
	    // factorizes the matrix
	    this->PerformAnalysisDistributed(comm, Ptr, IndRow, Val,
					     global_col_numbers, sym_matrix, reorder_num);
	  }
	else
	  {
	    Vector<long> Ptr; Vector<int> IndRow;
	    Vector<T> Val;
	    bool sym_pattern = false;
	    if (this->type_solver == this->PASTIX)
	      sym_pattern = true;
	    
	    bool reorder_num = false;
	    if (this->type_solver == this->SUPERLU)
	      {
		General sym;
		reorder_num = true;
		AssembleDistributed(B, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }
	    else
	      {
		Prop0 sym;
		AssembleDistributed(B, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }

	    if ((this->print_level >= 1) && (rank_proc == 0))
	      cout << "Analysing the distributed matrix..." << endl;
	    
	    // analyses the matrix
	    this->PerformAnalysisDistributed(comm, Ptr, IndRow, Val,
					     global_col_numbers, sym_matrix, reorder_num);
	  }        

      }
  }
  

  //! Performs the numerical factorization
  template<class T>
  template<class Prop0, class Storage0, class Allocator0>
  void SparseDistributedSolver<T>
  ::PerformFactorization(DistributedMatrix<T, Prop0, Storage0, Allocator0>& A,
			 bool scale_matrix)
  {
    diagonal_scaling_left = false;
    diagonal_scaling_right = false;
    if (scale_matrix)
      ScaleMatrixRowCol(A);

    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc, rank_proc;
    MPI_Comm_size(comm, &nb_proc); MPI_Comm_rank(comm, &rank_proc);
    if (nb_proc == 1)
      {
        comm_ = comm;
        SparseDirectSolver<T>::PerformFactorization(A);
      }
    else
      {
        // factorisation of a distributed matrix
        bool sym_matrix = IsSymmetricMatrix(A);
        
        // assembles distributed matrix
        if ((this->print_level >= 1) &&(rank_proc == 0))
          cout << "Assembling the distributed matrix..." << endl;
        
	if (this->solver->UseInteger8())
	  {
	    Vector<int64_t> Ptr, IndRow;
	    Vector<T> Val;
	    bool sym_pattern = false;
	    if (this->type_solver == this->PASTIX)
	      sym_pattern = true;
	    
	    bool reorder_num = false;
	    if (this->type_solver == this->SUPERLU)
	      {
		General sym;
		reorder_num = true;
		AssembleDistributed(A, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }
	    else
	      {
		Prop0 sym;
		AssembleDistributed(A, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }

	    if ((this->print_level >= 1) && (rank_proc == 0))
	      cout << "Factorizing the distributed matrix..." << endl;
	    
	    // factorizes the matrix
	    this->PerformFactorizationDistributed(comm, Ptr, IndRow, Val,
						  global_col_numbers, sym_matrix, reorder_num);
	  }
	else
	  {
	    Vector<long> Ptr; Vector<int> IndRow;
	    Vector<T> Val;
	    bool sym_pattern = false;
	    if (this->type_solver == this->PASTIX)
	      sym_pattern = true;
	    
	    bool reorder_num = false;
	    if (this->type_solver == this->SUPERLU)
	      {
		General sym;
		reorder_num = true;
		AssembleDistributed(A, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }
	    else
	      {
		Prop0 sym;
		AssembleDistributed(A, sym, comm, global_col_numbers,
				    local_col_numbers,
				    Ptr, IndRow, Val, sym_pattern, reorder_num);
	      }

	    if ((this->print_level >= 1) && (rank_proc == 0))
	      cout << "Factorizing the distributed matrix..." << endl;
	    
	    // factorizes the matrix
	    this->PerformFactorizationDistributed(comm, Ptr, IndRow, Val,
						  global_col_numbers, sym_matrix, reorder_num);
	  }        

      }
  }

  
  //! assembles the vector (adds values of shared rows)
  template<class T> template<class T2>
  void SparseDistributedSolver<T>::AssembleVec(Vector<T2>& X) const
  {
    AssembleVector(X, MPI_SUM, *ProcSharingRows_, *SharingRowNumbers_,
                   comm_, nodl_scalar_, nb_unknowns_scal_, 20);
  }

  //! assembles the vector (adds values of shared rows)
  template<class T> template<class T2>
  void SparseDistributedSolver<T>::AssembleVec(Matrix<T2, General, ColMajor>& A) const
  {
    int nrhs = A.GetN();
    Vector<T2> X;

    for (int k = 0; k < nrhs; k++)
      {
	X.SetData(A.GetM(), &A(0, k));
	AssembleVector(X, MPI_SUM, *ProcSharingRows_, *SharingRowNumbers_,
		       comm_, nodl_scalar_, nb_unknowns_scal_, 21);
    
	X.Nullify();
      }
  }
#endif
  
  
  //! solution of linear system Ax = b by using LU factorization
  //! (with scaling)
  /*!
    \param[out] x_solution solution
    \param[in] b_rhs right hand side
   */
  template<class T> template<class T1>
  void SparseDistributedSolver<T>::Solve(Vector<T1>& x_solution,
					 const Vector<T1>& b_rhs)
  {
    Copy(b_rhs, x_solution);
    Solve(x_solution);
  }
  
  
  //! solution of linear system Ax = b by using LU factorization
  //! (with scaling)
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T> template<class T1>
  void SparseDistributedSolver<T>::Solve(Vector<T1>& x_solution)
  {
    Solve(SeldonNoTrans, x_solution);
  }
  
  
  //! solution of linear system A^T x = b by using LU factorization
  //! (with scaling)
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T> template<class T1>
  void SparseDistributedSolver<T>::TransSolve(Vector<T1>& x_solution)
  {
    Solve(SeldonTrans, x_solution);
  }
  

  template<class T> template<class T1>
  void SparseDistributedSolver<T>::Solve(const SeldonTranspose& trans,
                                         Vector<T1>& x_solution, bool assemble)
  {
    if (diagonal_scaling_left && trans.NoTrans())
      for (int i = 0; i < x_solution.GetM(); i++)
	x_solution(i) *= diagonal_scale_left(i);
    
    if (diagonal_scaling_right && trans.Trans())
      for (int i = 0; i < x_solution.GetM(); i++)
	x_solution(i) *= diagonal_scale_right(i);
    
    // SolveLU or SolveLU_Distributed is used to handle the case of solving
    // a complex right hand side with a real matrix
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      SolveLU(trans, *this, x_solution);
    else
      {
        // extracting right hand side (we remove overlapped dofs)
        int n = local_col_numbers.GetM();
        Vector<T1> x_sol_extract(n);
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_sol_extract(i) = x_solution(local_col_numbers(i));
        
        SolveLU_Distributed(comm, trans, *this,
                            x_sol_extract, global_col_numbers);
        
        x_solution.Zero();
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_solution(local_col_numbers(i)) = x_sol_extract(i);
        
        // adding overlapped components
        if (assemble)
	  this->AssembleVec(x_solution);
      }
#else
    SolveLU(trans, *this, x_solution);
#endif
    
    if (diagonal_scaling_right && trans.NoTrans())
      for (int i = 0; i < x_solution.GetM(); i++)
        x_solution(i) *= diagonal_scale_right(i);

    if (diagonal_scaling_left && trans.Trans())
      for (int i = 0; i < x_solution.GetM(); i++)
        x_solution(i) *= diagonal_scale_left(i);
  }
  
  
  //! solution of linear system Ax = b by using LU factorization (with scaling)
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T> template<class T1>
  void SparseDistributedSolver<T>::Solve(Matrix<T1, General, ColMajor>& x_solution)
  {
    Solve(SeldonNoTrans, x_solution);
  }


  //! solution of linear system A^T x = b by using LU factorization (with scaling)
  /*!
    \param[in,out] x_solution on input right hand sides, on output solution
   */
  template<class T> template<class T1>
  void SparseDistributedSolver<T>::TransSolve(Matrix<T1, General, ColMajor>& x_solution)
  {
    Solve(SeldonTrans, x_solution);
  }


  template<class T> template<class T1>
  void SparseDistributedSolver<T>::Solve(const SeldonTranspose& trans,
                                         Matrix<T1, General, ColMajor>& x_solution)
  {
    if (diagonal_scaling_left && trans.NoTrans())
      ScaleLeftMatrix(x_solution, diagonal_scale_left);
   
    if (diagonal_scaling_right && trans.Trans())
      ScaleLeftMatrix(x_solution, diagonal_scale_right);
 
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      SolveLU(trans, *this, x_solution);
    else
      {
        // extracting right hand side (we remove overlapped dofs)
	int n = local_col_numbers.GetM();
	int N = x_solution.GetM(), nrhs = x_solution.GetN();
        Matrix<T1, General, ColMajor> x_sol_extract(n, nrhs);
	for (int k = 0; k < nrhs; k++)
	  for (int i = 0; i < local_col_numbers.GetM(); i++)
	    x_sol_extract(i, k) = x_solution(local_col_numbers(i), k);
        
	x_solution.Clear();
        SolveLU_Distributed(comm, trans, *this,
                            x_sol_extract, global_col_numbers);
        
	x_solution.Reallocate(N, nrhs);
        x_solution.Zero();
	for (int k = 0; k < nrhs; k++)
	  for (int i = 0; i < local_col_numbers.GetM(); i++)
	    x_solution(local_col_numbers(i), k) = x_sol_extract(i, k);
        
        // adding overlapped components
        this->AssembleVec(x_solution);
      }
#else
    SolveLU(trans, *this, x_solution);
#endif
    
    if (diagonal_scaling_right && trans.NoTrans())
      ScaleLeftMatrix(x_solution, diagonal_scale_right);

    if (diagonal_scaling_left && trans.Trans())
      ScaleLeftMatrix(x_solution, diagonal_scale_left);
  }
  
  
  //! Computation of a Schur Complement
  /*!
    \param[in] mat_direct initial matrix
    \param[in] num row numbers of the Schur Complement
    \param[out] mat_schur Schur Complement
   */
  template<class T> template<class MatrixSparse, class MatrixFull>
  void SparseDistributedSolver<T>
    ::GetSchurComplement(MatrixSparse& mat_direct, const IVect& num,
                         MatrixFull& mat_schur)
  {
    if (this->type_solver == this->MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	MatrixMumps<T>& mat_mumps =
	  dynamic_cast<MatrixMumps<T>& >(*this->solver);
	
	GetSchurMatrix(mat_direct, this->mat_mumps, num, mat_schur);
#else
	cout<<"Recompile Montjoie with Mumps."
	    << " Schur complement can't be performed otherwise"<<endl;
	abort();
#endif
      }
    else
      {
        cout << "Try to use Mumps, not implemented with other solvers " <<endl;
        abort();
      }
  }

  
  //! returns in bytes the memory used by the direct solver
  template<class T>
  size_t SparseDistributedSolver<T>::GetMemorySize() const
  {
    size_t taille = diagonal_scale_left.GetMemorySize()
      + diagonal_scale_right.GetMemorySize();
    
    taille += this->permut.GetMemorySize();
    taille += this->solver->GetMemorySize();
    
    return taille;
  }
  
} // namespace Seldon

#define SELDON_FILE_DISTRIBUTED_SOLVER_CXX
#endif
