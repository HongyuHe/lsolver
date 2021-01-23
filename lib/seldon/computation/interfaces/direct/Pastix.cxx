// Copyright (C) 2001-2010 Marc Duruflé
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

#ifndef SELDON_FILE_PASTIX_CXX

#include "Pastix.hxx"

namespace Seldon
{

  //! Default constructor.
  template<class T>
  MatrixPastix<T>::MatrixPastix()
  {
    comm_ = MPI_COMM_WORLD;
    pastix_data = NULL;
    n = 0;
    for (int i = 0; i < IPARM_SIZE; i++)
      iparm[i] = 0;
    
    for (int i = 0; i < DPARM_SIZE; i++)
      dparm[i] = 0;
    
    // Factorization of a matrix on a single processor.
    distributed = false;

    // No refinement by default.
    refine_solution = false;

    // initializing parameters
    pastix_initParam(iparm, dparm);

    iparm[IPARM_RHS_MAKING] = API_RHS_B;
    iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
    iparm[IPARM_FREE_CSCUSER] = API_CSC_FREE;
    if (refine_solution)
      iparm[IPARM_FREE_CSCPASTIX] = API_CSC_PRESERVE;
    else
      iparm[IPARM_FREE_CSCPASTIX] = API_CSC_FREE;
    
    //iparm[IPARM_BINDTHRD] = API_BIND_NO;
    
    threshold_pivot = 0.0;
    adjust_threshold_pivot = false;
    cholesky = false;
  }


  //! destructor
  template<class T>
  MatrixPastix<T>::~MatrixPastix()
  {
    Clear();
  }


  //! Calling main C-function pastix.
  template<>
  void MatrixPastix<double>::
  CallPastix(const MPI_Comm& comm, pastix_int_t* colptr, pastix_int_t* row,
             double* val, double* b, pastix_int_t nrhs)
  {
    comm_ = comm;
    if (distributed)
      d_dpastix(&pastix_data, comm, n, colptr, row, val,
                col_num.GetData(), perm.GetData(), invp.GetData(),
                b, nrhs, iparm, dparm);
    else
      d_pastix(&pastix_data, comm, n, colptr, row, val,
               perm.GetData(), invp.GetData(), b, nrhs, iparm, dparm);
  }


  //! Calling main C-function pastix.
  template<>
  void MatrixPastix<complex<double> >::
  CallPastix(const MPI_Comm& comm, pastix_int_t* colptr, pastix_int_t* row,
             complex<double>* val, complex<double>* b, pastix_int_t nrhs)
  {
    comm_ = comm;
    if (distributed)
      z_dpastix(&pastix_data, comm, n, colptr, row, val,
                col_num.GetData(), perm.GetData(), invp.GetData(),
                b, nrhs, iparm, dparm);
    else
      z_pastix(&pastix_data, comm, n, colptr, row, val,               
               perm.GetData(), invp.GetData(),
               b, nrhs, iparm, dparm);
  }


  template<class T>
  bool MatrixPastix<T>::UseInteger8() const  
  {
    if (sizeof(pastix_int_t) == 8)
      return true;
    
    return false;
  }


  //! Clearing factorization.
  template<class T>
  void MatrixPastix<T>::Clear()
  {
    if (n > 0)
      {
	pastix_int_t nrhs = 1;
	iparm[IPARM_START_TASK] = API_TASK_CLEAN;
	iparm[IPARM_END_TASK] = API_TASK_CLEAN;

	CallPastix(comm_, NULL, NULL, NULL, NULL, nrhs);

	perm.Clear();
        invp.Clear();
        col_num.Clear();
	n = 0;
        pastix_data = NULL;
        distributed = false;
      }
  }


  //! no message will be displayed
  template<class T>
  void MatrixPastix<T>::HideMessages()
  {
    iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
  }


  //! Low level of display.
  template<class T>
  void MatrixPastix<T>::ShowMessages()
  {
    iparm[IPARM_VERBOSE] = API_VERBOSE_NO;
  }


  //! Displaying all messages.
  template<class T>
  void MatrixPastix<T>::ShowFullHistory()
  {
    iparm[IPARM_VERBOSE] = API_VERBOSE_YES;
  }


  //! selects the algorithm used for reordering
  template<class T>
  void MatrixPastix<T>::SelectOrdering(int type)
  {
    iparm[IPARM_ORDERING] = type;
  }

  
  //! provides a permutation array (instead of using Scotch reordering)
  template<class T>
  void MatrixPastix<T>::SetPermutation(const IVect& permut)
  {
    iparm[IPARM_ORDERING] = API_ORDER_PERSONAL;
    perm.Reallocate(permut.GetM());
    invp.Reallocate(permut.GetM());
    for (int i = 0; i < perm.GetM(); i++)
      {
	perm(i) = permut(i) + 1;
	invp(permut(i)) = i+1;
      }
  }

  
  //! sets Cholesky factorisation
  template<class T>
  void MatrixPastix<T>::SetCholeskyFacto(bool chol)
  {
    cholesky = chol;
  }
  
  
  //! you can change the threshold used for static pivoting
  template<class T>
  void MatrixPastix<T>::SetPivotThreshold(double eps)
  {
    adjust_threshold_pivot = true;
    threshold_pivot = eps;
  }
  
  
  //! You can require that solution is refined after LU resolution.
  template<class T>
  void MatrixPastix<T>::RefineSolution()
  {
    refine_solution = true;
    iparm[IPARM_FREE_CSCPASTIX] = API_CSC_PRESERVE;
  }


  //! You can require that solution is not refined (faster).
  template<class T>
  void MatrixPastix<T>::DoNotRefineSolution()
  {
    refine_solution = false;
    iparm[IPARM_FREE_CSCPASTIX] = API_CSC_FREE;
  }

  
  //! Returns the size of memory used by the factorisation in bytes
  template<class T>
  size_t MatrixPastix<T>::GetMemorySize() const
  {
    size_t taille = sizeof(pastix_int_t)*(perm.GetM()+invp.GetM()+col_num.GetM());
    if (n <= 0)
      return taille;
    
    // assuming that for each term, a value and an index is needed
    taille += (sizeof(T)+sizeof(pastix_int_t))*iparm[IPARM_NNZEROS];
    return taille;
  }
  

  template<class T>
  int MatrixPastix<T>::GetInfoFactorization() const
  {
    return 0;
  }

  
  //! Returning ordering found by Scotch.
  template<class T>
  template<class T0, class Prop, class Storage, class Allocator, class Tint>
  void MatrixPastix<T>::
  FindOrdering(Matrix<T0, Prop, Storage, Allocator> & mat,
               Vector<Tint>& numbers, bool keep_matrix)
  {
    // We clear the previous factorization, if any.
    Clear();

    distributed = false;

    n = mat.GetN();
    if (n <= 0)
      return;

    pastix_int_t nrhs = 1, nnz = 0;
    pastix_int_t* ptr_ = NULL;
    pastix_int_t* ind_ = NULL;
    Vector<pastix_int_t> Ptr, Ind;

    iparm[IPARM_SYM] = API_SYM_YES;
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    GetSymmetricPattern(mat, Ptr, Ind);
    if (!keep_matrix)
      mat.Clear();

    ptr_ = Ptr.GetData();
    // Changing to 1-index notation.
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    nnz = Ind.GetM();
    ind_ = Ind.GetData();
    for (long i = 0; i < nnz; i++)
      ind_[i]++;

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // We get ordering only.
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_ORDERING;

    CallPastix(MPI_COMM_SELF, ptr_, ind_, NULL, NULL, nrhs);

    numbers.Reallocate(perm.GetM());
    for (int i = 0; i < perm.GetM(); i ++)
      numbers(i) = perm(i);
    
    Ptr.Nullify();
    Ind.Nullify();
  }


  //! Factorization of unsymmetric matrix
  template<class T> template<class T0, class Storage, class Allocator>
  void MatrixPastix<T>
  ::FactorizeMatrix(Matrix<T0, General, Storage, Allocator> & mat,
                    bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    Vector<pastix_int_t> Ptr, IndRow;
    Vector<T> Val;

    General prop;
    ConvertToCSC(mat, prop, Ptr, IndRow, Val, true);
    if (!keep_matrix)
      mat.Clear();

    FactorizeCSC(Ptr, IndRow, Val, false);
  }

  
  template<class T>
  void MatrixPastix<T>
  ::FactorizeCSC(Vector<pastix_int_t>& Ptr, Vector<pastix_int_t>& IndRow,
		 Vector<T>& Val, bool sym)
  {
    distributed = false;
    n = Ptr.GetM()-1;
    if (n <= 0)
      return;
    
    pastix_int_t nrhs = 1, nnz = 0;
    pastix_int_t* ptr_ = NULL;
    pastix_int_t* ind_ = NULL;
    T* values_ = NULL;

    ptr_ = Ptr.GetData();
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    nnz = IndRow.GetM();
    ind_ = IndRow.GetData();
    for (long i = 0; i < nnz; i++)
      ind_[i]++;

    values_ = Val.GetData();

    if (sym)
      {
	iparm[IPARM_SYM]           = API_SYM_YES;
	if (cholesky)
	  iparm[IPARM_FACTORIZATION] = API_FACT_LLT;
	else
	  iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
      }
    else
      {
	iparm[IPARM_SYM]           = API_SYM_NO;    
	iparm[IPARM_FACTORIZATION] = API_FACT_LU;
      }

    if (iparm[IPARM_ORDERING] != API_ORDER_PERSONAL)
      {
	perm.Reallocate(n); invp.Reallocate(n);
	perm.Fill(); invp.Fill();
      }

    // pivot threshold
    if (adjust_threshold_pivot)
      dparm[DPARM_EPSILON_MAGN_CTRL] = threshold_pivot;

    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_ANALYSE;

    CallPastix(MPI_COMM_SELF, ptr_, ind_, values_, NULL, nrhs);

    // factorization only
    Vector<pastix_int_t> proc_num(iparm[IPARM_THREAD_NBR]);
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
    proc_num.Fill(rank_proc);
    pastix_bindThreads(pastix_data, iparm[IPARM_THREAD_NBR], proc_num.GetData());

    iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;

    CallPastix(MPI_COMM_SELF, ptr_, ind_, values_, NULL, nrhs);
    
    if (iparm[IPARM_STATIC_PIVOTING] > 0)
      {
        if (!refine_solution)
          {
            cout << "Refining solution is needed when pivoting is used" << endl;
            abort();
          }
      }
    
    if (iparm[IPARM_VERBOSE] != API_VERBOSE_NOT)
      cout << "Factorization successful" << endl;
    
    Ptr.Nullify();
    IndRow.Nullify();
    Val.Nullify();
  }


  //! Factorization of symmetric matrix.
  template<class T> template<class T0, class Storage, class Allocator>
  void MatrixPastix<T>::
  FactorizeMatrix(Matrix<T0, Symmetric, Storage, Allocator> & mat,
                  bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    Vector<pastix_int_t> Ptr, IndRow;
    Vector<T> Val;

    Symmetric prop;
    ConvertToCSR(mat, prop, Ptr, IndRow, Val);

    FactorizeCSC(Ptr, IndRow, Val, true);
  }


  //! solving A x = b (A is already factorized)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! solving A x = b or A^T x = b (A is already factorized)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Solve(const SeldonTranspose& TransA,
                              Vector<T, VectFull, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), 1);
  }


  //! solving A x = b or A^T x = b (A is already factorized)
  template<class T>
  void MatrixPastix<T>::Solve(const SeldonTranspose& TransA, T* x_ptr, int nrhs_)
  {
    pastix_int_t nrhs = nrhs_;

    if (cholesky)
      {
        if (TransA.Trans())
          iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_BACKWARD_ONLY;
        else
          iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_FORWARD_ONLY;
        
        iparm[IPARM_END_TASK] = API_TASK_SOLVE;
      }
    else
      {
        if (TransA.Trans())
          iparm[IPARM_TRANSPOSE_SOLVE] = API_YES;
        else
          iparm[IPARM_TRANSPOSE_SOLVE] = API_NO;
        
	// no refinement for multiple right hand sides
        if ((refine_solution) && (nrhs == 1))
          iparm[IPARM_END_TASK] = API_TASK_REFINE;
        else
          iparm[IPARM_END_TASK] = API_TASK_SOLVE;
      }
    
    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    
    CallPastix(MPI_COMM_SELF, NULL, NULL, NULL, x_ptr, nrhs);
  }
  

  //! solving A x = b or A^T x = b (A is already factorized)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Solve(const SeldonTranspose& TransA,
                              Matrix<T, General, ColMajor, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), x.GetN());
  }


  //! Computes x = L b or x = L^T b (A = L L^T)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Mlt(const SeldonTranspose& TransA,
			    Vector<T, VectFull, Allocator2>& x)
  {
    Mlt(TransA, x.GetData(), 1);
  }


  //! Computes x = L b or x = L^T b (A = L L^T)
  template<class T>
  void MatrixPastix<T>::Mlt(const SeldonTranspose& TransA, T* x_ptr, int nrhs_)
  {
    pastix_int_t nrhs = nrhs_;

    if (cholesky)
      {
        if (TransA.Trans())
          iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_UTRMV;
        else
          iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_LTRMV;
        
        iparm[IPARM_END_TASK] = API_TASK_SOLVE;
      }
    else
      {
        cout << "Mlt is defined for Cholesky only" << endl;
        abort();
      }
    
    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    
    CallPastix(MPI_COMM_SELF, NULL, NULL, NULL, x_ptr, nrhs);
  }


  //! Modifies the number of threads per node.
  template<class T>
  void MatrixPastix<T>::SetNumberOfThreadPerNode(int num_thread)
  {
    iparm[IPARM_THREAD_NBR] = num_thread;
  }


  template<class T>
  void MatrixPastix<T>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
                             Vector<int>& IndRow, Vector<T>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    Vector<int> PtrInt(Ptr.GetM());
    for (int i = 0; i < Ptr.GetM(); i++)
      PtrInt(i) = Ptr(i);

    FactorizeParallel(comm_facto, PtrInt, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }

  
  template<class T>
  void MatrixPastix<T>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
                             Vector<int64_t>& IndRow, Vector<T>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    FactorizeParallel(comm_facto, Ptr, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }
  
  
  template<class T> template<class Tint>
  void MatrixPastix<T>::  
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<Tint>&, Vector<Tint>&,
                    Vector<T>&,
                    const Vector<int>& glob_num,
                    bool sym, bool keep_matrix)
  {
    cout << "Not available for this type of integers " << endl;
    cout << "Size of Tint = " << sizeof(Tint) << endl;
    abort();
  }
  

  //! Distributed factorization (on several nodes).
  template<class T>
  void MatrixPastix<T>::
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<pastix_int_t>& Ptr,
                    Vector<pastix_int_t>& IndRow,
                    Vector<T>& Val, const Vector<int>& glob_number,
                    bool sym, bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    n = Ptr.GetM() - 1;
    if (n <= 0)
      return;
    
    distributed = true;

    if (sym)
      {
        iparm[IPARM_SYM] = API_SYM_YES;
        if (cholesky)
          iparm[IPARM_FACTORIZATION] = API_FACT_LLT;
        else
          iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
      }
    else
      {
        iparm[IPARM_SYM] = API_SYM_NO;
        iparm[IPARM_FACTORIZATION] = API_FACT_LU;
      }
    
    iparm[IPARM_GRAPHDIST] = API_YES;
    iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;
    
    pastix_int_t* ptr_ = Ptr.GetData();
    pastix_int_t nrhs = 1;
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    pastix_int_t nnz = IndRow.GetM();
    pastix_int_t* ind_ = IndRow.GetData();
    for (long i = 0; i < nnz; i++)
      ind_[i]++;

    T* values_ = Val.GetData();

    col_num.Reallocate(n);
    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();
    for (int i = 0; i < n; i++)
      col_num(i) = glob_number(i)+1;

    if (adjust_threshold_pivot)
      dparm[DPARM_EPSILON_MAGN_CTRL] = threshold_pivot;

    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_ANALYSE;

    CallPastix(comm_facto, ptr_, ind_, values_, NULL, nrhs);

    Vector<pastix_int_t> proc_num(iparm[IPARM_THREAD_NBR]);
    int rank_proc; MPI_Comm_rank(comm_facto, &rank_proc);
    proc_num.Fill(rank_proc);
    pastix_bindThreads(pastix_data, iparm[IPARM_THREAD_NBR], proc_num.GetData());

    // factorization only
    iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
    
    //Vector<T> rhs(n); rhs.Zero();
    //T* rhs_ = rhs.GetData();
    //CallPastix(comm_facto, ptr_, ind_, values_, rhs_, nrhs);    
    CallPastix(comm_facto, ptr_, ind_, values_, NULL, nrhs);    

    //iparm[IPARM_FREE_CSCUSER] = API_CSC_FREE;
    //Ptr.Nullify();
    //IndRow.Nullify();
    //Val.Nullify();
  }

  
  //! solves A x = b or A^T x = b in parallel
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::SolveDistributed(MPI_Comm& comm_facto,
                                         const SeldonTranspose& TransA,
                                         Vector<T, VectFull, Allocator2>& x,
                                         const Vector<int>& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), 1, glob_num);
  }


  //! solves A x = b or A^T x = b in parallel
  template<class T>
  void MatrixPastix<T>::SolveDistributed(MPI_Comm& comm_facto,
					 const SeldonTranspose& TransA,
					 T* x_ptr, int nrhs_,
					 const IVect& glob_num)
  {
    pastix_int_t nrhs = nrhs_;
    
    // uncomment to have something working
    for (int k = 0; k < nrhs_; k++)
      {
        // workaround : each linear system is solved separately
        nrhs = 1;
        
        if (cholesky)
          {
            if (TransA.Trans())
              iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_BACKWARD_ONLY;
            else
              iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_FORWARD_ONLY;
            
            iparm[IPARM_END_TASK] = API_TASK_SOLVE;
          }
        else
          {
            if (TransA.Trans())
              iparm[IPARM_TRANSPOSE_SOLVE] = API_YES;
            else
              iparm[IPARM_TRANSPOSE_SOLVE] = API_NO;
            
            if (refine_solution)
              {
                // no refinement for multiple right hand sides
                if (nrhs > 1)
                  iparm[IPARM_END_TASK] = API_TASK_SOLVE;
                else
                  iparm[IPARM_END_TASK] = API_TASK_REFINE;
              }
            else
              iparm[IPARM_END_TASK] = API_TASK_SOLVE;        
          }
        
        iparm[IPARM_START_TASK] = API_TASK_SOLVE;
        CallPastix(comm_facto, NULL, NULL, NULL, x_ptr, nrhs);
        
        // uncomment here
        x_ptr += n;
      }
  }


  //! solves A x = b or A^T x = b in parallel
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::SolveDistributed(MPI_Comm& comm_facto,
                                         const SeldonTranspose& TransA,
                                         Matrix<T, General, ColMajor, Allocator2>& x,
                                         const Vector<int>& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), x.GetN(), glob_num);
  }

  
  //! computes L x or L^T x
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::MltDistributed(MPI_Comm& comm_facto,
				       const SeldonTranspose& TransA,
				       Vector<T, VectFull, Allocator2>& x,
				       const Vector<int>& glob_num)
  {
    MltDistributed(comm_facto, TransA, x.GetData(), 1, glob_num);
  }


  //! computes L x or L^T x
  template<class T>
  void MatrixPastix<T>::MltDistributed(MPI_Comm& comm_facto,
				       const SeldonTranspose& TransA,
				       T* x_ptr, int nrhs_,
				       const IVect& glob_num)
  {
    pastix_int_t nrhs = nrhs_;

    for (int k = 0; k < nrhs_; k++)
      {
	// workaround : each linear system is solved separately
	nrhs = 1;
        
        if (cholesky)
          {
            if (TransA.Trans())
              iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_UTRMV;
            else
              iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_LTRMV;
            
            iparm[IPARM_END_TASK] = API_TASK_SOLVE;
          }
        else
          {
            cout << "Mlt is defined for Cholesky only" << endl;
            abort();
          }
        
        iparm[IPARM_START_TASK] = API_TASK_SOLVE;
        //CallPastix(comm_facto, NULL, NULL, NULL, x_ptr, nrhs);
        
        CallPastix(comm_facto, NULL, NULL, NULL, x_ptr, nrhs);
        x_ptr += n;
      }
  }


  //! solves A x = b or A^T x = b in parallel
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::MltDistributed(MPI_Comm& comm_facto,
				       const SeldonTranspose& TransA,
				       Matrix<T, General, ColMajor, Allocator2>& x,
				       const Vector<int>& glob_num)
  {
    MltDistributed(comm_facto, TransA, x.GetData(), x.GetN(), glob_num);
  }


  //! Factorization of a matrix of same type T as the Pastix object 
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPastix<T>& mat_lu, bool keep_matrix, T& x)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  

  //! Factorization of a complex matrix with a real Pastix object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPastix<T>& mat_lu, bool keep_matrix, complex<T>& x)
  {
    throw WrongArgument("GetLU(Matrix<complex<T> >& A, MatrixPastix<T>& mat_lu, bool)",
			"The LU matrix must be complex");
  }

  
  //! Factorization of a real matrix with a complex Pastix object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPastix<complex<T> >& mat_lu, bool keep_matrix, T& x)
  {
    throw WrongArgument("GetLU(Matrix<T>& A, MatrixMumps<Pastix<T> >& mat_lu, bool)",
			"The sparse matrix must be complex");
  }
  
  
  //! Factorization of a general matrix with Pastix
  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixPastix<T>& mat_lu,
	     bool keep_matrix)
  {
    // we check if the type of non-zero entries of matrix A
    // and of the Pastix object (T) are different
    // we call one of the GetLUs written above
    // such a protection avoids to compile the factorisation of a complex
    // matrix with a real Pastix object
    typename Matrix<T0, Prop, Storage, Allocator>::entry_type x;
    GetLU(A, mat_lu, keep_matrix, x);
  }
  

  //! LU resolution with a vector whose type is the same as Pastix object
  template<class T, class Allocator>
  void SolveLU(MatrixPastix<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  //! LU resolution with a vector whose type is the same as Pastix object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! LU resolution with a matrix whose type is the same as Pastix object
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixPastix<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(SeldonNoTrans, x);
  }


  //! LU resolution with a matrix whose type is the same as UmfPack object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! Solves A x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(MatrixPastix<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x)
  {
    Matrix<double, General, ColMajor> y(x.GetM(), 2);
    
    for (int i = 0; i < x.GetM(); i++)
      {
	y(i, 0) = real(x(i));
	y(i, 1) = imag(x(i));
      }
    
    SolveLU(mat_lu, y);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<double>(y(i, 0), y(i, 1));
  }
  

  //! Solves A x = b or A^T x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<double>& mat_lu, Vector<complex<double>, VectFull, Allocator>& x)
  {
    Matrix<double, General, ColMajor> y(x.GetM(), 2);
    
    for (int i = 0; i < x.GetM(); i++)
      {
	y(i, 0) = real(x(i));
	y(i, 1) = imag(x(i));
      }
    
    SolveLU(TransA, mat_lu, y);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<double>(y(i, 0), y(i, 1));

  }


  //! Solves A x = b, where A is complex and x is real => Forbidden
  template<class Allocator>
  void SolveLU(MatrixPastix<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixPastix<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }

  
  //! Solves A x = b or A^T x = b, where A is complex and x is real => Forbidden  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixPastix<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }


  template<class T, class Prop, class Storage, class Allocator>
  void GetCholesky(Matrix<T, Prop, Storage, Allocator>& A,
                   MatrixPastix<T>& mat_chol, bool keep_matrix)
  {
    mat_chol.SetCholeskyFacto(true);
    //IVect permut(A.GetM()); permut.Fill();
    //mat_chol.SetPermutation(permut);
    mat_chol.FactorizeMatrix(A, keep_matrix);
  }


  template<class T, class Allocator>
  void
  SolveCholesky(const SeldonTranspose& TransA,
                MatrixPastix<T>& mat_chol, Vector<T, VectFull, Allocator>& x)
  {
    mat_chol.Solve(TransA, x);
  }


  template<class T, class Allocator>
  void
  MltCholesky(const SeldonTranspose& TransA,
              MatrixPastix<T>& mat_chol, Vector<T, VectFull, Allocator>& x)
  {
    mat_chol.Mlt(TransA, x);
  }

} // end namespace

#define SELDON_FILE_PASTIX_CXX
#endif
