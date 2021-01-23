// Copyright (C) 2010 Marc Duruflé
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


#ifndef SELDON_FILE_DISTRIBUTED_CHOLESKY_SOLVER_CXX


namespace Seldon
{

  /*****************************
   * DistributedCholeskySolver *
   *****************************/


  //! Default constructor
  template <class T>
  DistributedCholeskySolver<T>::DistributedCholeskySolver()
    : SparseCholeskySolver<T>()
  {
#ifdef SELDON_WITH_MPI
    // initialisation for a sequential matrix
    nodl_scalar_ = 1;
    nb_unknowns_scal_ = 1;
    comm_ = MPI_COMM_SELF;
    ProcSharingRows_ = NULL;
    SharingRowNumbers_ = NULL;
#endif
  }


  //! factorization of matrix in sequential
  template<class T>
  template<class Prop, class Storage, class Allocator>
  void DistributedCholeskySolver<T>
  ::Factorize(Matrix<T, Prop, Storage, Allocator>& A, bool keep_matrix)
  {
    // we call the sequential function
    SparseCholeskySolver<T>::Factorize(A, keep_matrix);
  }

  
#ifdef SELDON_WITH_MPI
  //! factorization of distributed matrix
  template<class T>
  template<class Prop0, class Storage0, class Allocator0>
  void DistributedCholeskySolver<T>
  ::Factorize(DistributedMatrix<T, Prop0, Storage0, Allocator0>& A, bool keep_matrix)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        comm_ = comm;
        SparseCholeskySolver<T>::Factorize(A, keep_matrix);
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

	if (this->type_solver != SparseCholeskySolver<T>::PASTIX)
	  {
	    cout << "Only available for Pastix" << endl;
	    abort();
	  }

	DistributedMatrix<T, Prop0, Storage0, Allocator0> Bstore;
	DistributedMatrix<T, Prop0, Storage0, Allocator0>* B;
	B = &A;
	if (keep_matrix)
	  {
	    B = &Bstore;
	    Bstore = A;
	  }

#ifdef SELDON_WITH_PASTIX
	MatrixPastix<T>& mat_pastix = dynamic_cast<MatrixPastix<T>&>(*this->solver);
	mat_pastix.SetCholeskyFacto(true);
#else
	cout << "Recompile Seldon with Pastix" << endl;
	abort();
#endif

	bool sym_pattern = true, reorder_num = false;
	Prop0 sym;
	if (this->solver->UseInteger8())
	  {
	    Vector<int64_t> Ptr, IndRow;
	    Vector<T> Val;
	    
	    AssembleDistributed(*B, sym, comm, global_col_numbers,
				local_col_numbers,
				Ptr, IndRow, Val, sym_pattern, reorder_num);
	    
	    // factorizes the matrix
	    this->solver->FactorizeDistributedMatrix(comm, Ptr, IndRow, Val,
						     global_col_numbers, sym_matrix, reorder_num);
	  }
	else
	  {
	    Vector<long> Ptr; Vector<int> IndRow;
	    Vector<T> Val;

	    AssembleDistributed(*B, sym, comm, global_col_numbers,
				local_col_numbers,
				Ptr, IndRow, Val, sym_pattern, reorder_num);
	    
	    // factorizes the matrix
	    this->solver->FactorizeDistributedMatrix(comm, Ptr, IndRow, Val,
						     global_col_numbers, sym_matrix, reorder_num);
	  }        
      }
  }


  //! assembles the vector (adds values of shared rows)
  template<class T> template<class T2>
  void DistributedCholeskySolver<T>::AssembleVec(Vector<T2>& X) const
  {
    AssembleVector(X, MPI_SUM, *ProcSharingRows_, *SharingRowNumbers_,
                   comm_, nodl_scalar_, nb_unknowns_scal_, 20);
  }

  
  //! assembles the vector (adds values of shared rows)
  template<class T> template<class T2>
  void DistributedCholeskySolver<T>::AssembleVec(Matrix<T2, General, ColMajor>& A) const
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


  template<class T> template<class T1>
  void DistributedCholeskySolver<T>::Solve(const SeldonTranspose& trans,
					   Vector<T1>& x_solution, bool assemble)
  {

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      SparseCholeskySolver<T>::Solve(trans, x_solution);
    else
      {
        // extracting right hand side (we remove overlapped dofs)
        int n = local_col_numbers.GetM();
        Vector<T1> x_sol_extract(n);
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_sol_extract(i) = x_solution(local_col_numbers(i));

#ifdef SELDON_WITH_PASTIX
	MatrixPastix<T>& mat_pastix = dynamic_cast<MatrixPastix<T>&>(*this->solver);
	
	mat_pastix.SolveDistributed(comm, trans, x_sol_extract, global_col_numbers);
#else
	cout << "Recompile Seldon with Pastix" << endl;
	abort();
#endif

        x_solution.Zero();
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_solution(local_col_numbers(i)) = x_sol_extract(i);
        
        // adding overlapped components
        if (assemble)
	  this->AssembleVec(x_solution);
      }
#else
    SparseCholeskySolver<T>::Solve(trans, x_solution);
#endif
  }


  template<class T> template<class T1>
  void DistributedCholeskySolver<T>::Mlt(const SeldonTranspose& trans,
					 Vector<T1>& x_solution, bool assemble)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      SparseCholeskySolver<T>::Mlt(trans, x_solution);
    else
      {
        // extracting right hand side (we remove overlapped dofs)
        int n = local_col_numbers.GetM();
        Vector<T1> x_sol_extract(n);
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_sol_extract(i) = x_solution(local_col_numbers(i));

#ifdef SELDON_WITH_PASTIX
	MatrixPastix<T>& mat_pastix = dynamic_cast<MatrixPastix<T>&>(*this->solver);

	mat_pastix.MltDistributed(comm, trans, x_sol_extract, global_col_numbers);
#else
	cout << "Recompile Seldon with Pastix" << endl;
	abort();
#endif
	
        x_solution.Zero();
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_solution(local_col_numbers(i)) = x_sol_extract(i);
        
        // adding overlapped components
        if (assemble)
	  this->AssembleVec(x_solution);
      }
#else
    SparseCholeskySolver<T>::Mlt(trans, x_solution);
#endif
  }

} // namespace Seldon.


#define SELDON_FILE_DISTRIBUTED_CHOLESKY_SOLVER_CXX
#endif
