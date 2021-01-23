// Copyright (C) 2015-2015 Marc Duruflé
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

#ifndef SELDON_FILE_WSMP_CXX

#include "Wsmp.hxx"

namespace Seldon
{

  template<>
  void MatrixWsmp<double>::CallWssmp(int* n_, int* ia, int* ja, double* avals, double* diag,
                                     int* perm, int* invp, double* b, int* ldb, int* nrhs,
                                     double* aux, int* naux, int* mrp, int* iparam, double* dparam)
  {
#ifdef SELDON_WITH_MPI
    pwssmp_(n_, ia, ja, avals, diag, perm, invp,
            b, ldb, nrhs, aux, naux, mrp, iparam, dparam);
#else
    wssmp_(n_, ia, ja, avals, diag, perm, invp,
           b, ldb, nrhs, aux, naux, mrp, iparam, dparam);
#endif
  }


  template<>
  void MatrixWsmp<complex<double> >
  ::CallWssmp(int* n_, int* ia, int* ja, complex<double>* avals, complex<double>* diag,
              int* perm, int* invp, complex<double>* b, int* ldb, int* nrhs,
              complex<double>* aux, int* naux, int* mrp, int* iparam, double* dparam)
  {
#ifdef SELDON_WITH_MPI
    pzssmp_(n_, ia, ja, avals, diag, perm, invp,
            b, ldb, nrhs, aux, naux, mrp, iparam, dparam);
#else
    zssmp_(n_, ia, ja, avals, diag, perm, invp,
           b, ldb, nrhs, aux, naux, mrp, iparam, dparam);
#endif
  }
  
  
  template<>
  void MatrixWsmp<double>::CallWgsmp(int* n_, int* ia, int* ja, double* avals,
                                     double* b, int* ldb, int* nrhs,
                                     double* rmisc, int* iparam, double* dparam)
  {
#ifdef SELDON_WITH_MPI
    pwgsmp_(n_, ia, ja, avals, b, ldb, nrhs, rmisc, iparam, dparam);
#else
    wgsmp_(n_, ia, ja, avals, b, ldb, nrhs, rmisc, iparam, dparam);
#endif
  }
  
  
  template<>
  void MatrixWsmp<complex<double> >
  ::CallWgsmp(int* n_, int* ia, int* ja, complex<double>* avals,
              complex<double>* b, int* ldb, int* nrhs,
              complex<double>* rmisc, int* iparam, double* dparam)
  {
#ifdef SELDON_WITH_MPI
    pzgsmp_(n_, ia, ja, avals, b, ldb, nrhs, rmisc, iparam, dparam);
#else
    zgsmp_(n_, ia, ja, avals, b, ldb, nrhs, rmisc, iparam, dparam);
#endif
  }


  template<class T>
  MatrixWsmp<T>::MatrixWsmp() : iparm(64), dparm(64)
  {
    iparm.Zero();
    dparm.Zero();
    n = 0;
    refine_solution = false;
    use_pivoting = true;
    cholesky = false;
    symmetric = false;
    distributed = false;
    threshold_pivot = 0.001;
    
    int nb_threads = 1;
    wsetmaxthrds_(&nb_threads);
  }
  
  
  template<class T>
  MatrixWsmp<T>::~MatrixWsmp() 
  {
    Clear();
  }
    

  template<class T>
  bool MatrixWsmp<T>::UseInteger8() const  
  {    
    return false;
  }

  
  template<class T>
  void MatrixWsmp<T>::Clear() 
  {
    if (n > 0)
      {
        wsmp_clear_();
        Ptr.Clear(); IndRow.Clear(); Val.Clear();
        permut.Clear(); inverse_permut.Clear();
        n = 0;
      }
  }
  
  
  template<class T>
  void MatrixWsmp<T>::ShowMessages()
  {
  }
  
  
  template<class T>
  void MatrixWsmp<T>::HideMessages()
  {
  }
  
  
  template<class T>
  void MatrixWsmp<T>::ShowFullHistory()
  {
  }
    

  template<class T>
  void MatrixWsmp<T>::SetPivotThreshold(double eps)
  {
    threshold_pivot = eps;
  }
  
  
  template<class T>
  void MatrixWsmp<T>::SetCholeskyFacto(bool chol)
  {
    cholesky = chol; 
  }
  
  
  template<class T>
  int MatrixWsmp<T>::GetInfoFactorization() const
  {
    return iparm(63); 
  }
  
  
  template<class T>
  void MatrixWsmp<T>::RefineSolution()
  {
    refine_solution = true; 
  }
  
  
  template<class T>
  void MatrixWsmp<T>::DoNotRefineSolution()
  {
    refine_solution = false; 
  }
  
  
  template<class T>
  size_t MatrixWsmp<T>::GetMemorySize() const
  {
    if (n > 0)
      return size_t(iparm(22))*1024*sizeof(T);
    
    return 0;
  }
  
  
  template<class T>
  void MatrixWsmp<T>::SetNumberOfThreadPerNode(int nb_threads)
  {
    wsetmaxthrds_(&nb_threads);
  }
  
  
  //! factorizes a symmetric matrix
  template<class T> template<class Storage, class Allocator>
  void MatrixWsmp<T>
  ::FactorizeMatrix(Matrix<T, Symmetric, Storage, Allocator>& mat,
                    bool keep_matrix)
  {
    Clear();

    Symmetric prop;
    ConvertToCSR(mat, prop, Ptr, IndRow, Val);
    if (!keep_matrix)
      mat.Clear();
    
    FactorizeSymmetric();
  }
  

  template<class T>
  void MatrixWsmp<T>::FactorizeCSR(Vector<int>& Ptr_, Vector<int>& IndRow_,
				   Vector<T>& Val_, bool sym)
  {
    Ptr.SetData(Ptr_);
    IndRow.SetData(IndRow_);
    Val.SetData(Val_);
    Ptr_.Nullify(); IndRow_.Nullify(); Val_.Nullify();
    
    if (sym)
      FactorizeSymmetric();
    else
      FactorizeUnsymmetric();
  }


  template<class T>
  void MatrixWsmp<T>::FactorizeSymmetric()
  {
    distributed = false;    
    symmetric = true;
    n = Ptr.GetM()-1;
    
    wsmp_initialize_();
    
    iparm(0) = 0;
    iparm(1) = 0;
    iparm(2) = 0;
    
    // initialisation step
    int nrhs = 0, naux = 0, mrp = 0;
    CallWssmp(&n, NULL, NULL, NULL, NULL, NULL, NULL,
              NULL, &n, &nrhs, NULL, &naux, &mrp, iparm.GetData(), dparm.GetData());      
    
    // CSC/CSR format
    iparm(3) = 0;
    // C-style index
    iparm(4) = 0;
    
    // wsmp can reuse the arrays IndRow and Val
    iparm(13) = 3;
    
    // pivot threshold
    dparm(10) = threshold_pivot;
    
    if (cholesky)
      {
        // L L^T factorization
        iparm(30) = 0;
      }
    else
      {
        // L D L^T factorization
        if (use_pivoting)
          iparm(30) = 2;
        else
          iparm(30) = 1;
        
        if (IsComplexNumber(T()))
          iparm(30) += 2;
      }
    
    // ordering step, symbolic step
    permut.Reallocate(n);
    inverse_permut.Reallocate(n);
    
    iparm(1) = 1;
    iparm(2) = 2;
    CallWssmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
              NULL, permut.GetData(), inverse_permut.GetData(), NULL, &n, &nrhs,
              NULL, &naux, &mrp, iparm.GetData(), dparm.GetData());
    
    // then factorization
    iparm(1) = 3;
    iparm(2) = 3;
    CallWssmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
              NULL, permut.GetData(), inverse_permut.GetData(), NULL, &n, &nrhs,
              NULL, &naux, &mrp, iparm.GetData(), dparm.GetData());
    
  }


  //! factorizes a symmetric matrix
  template<class T> template<class Storage, class Allocator>
  void MatrixWsmp<T>
  ::FactorizeMatrix(Matrix<T, General, Storage, Allocator>& mat,
                    bool keep_matrix)
  {
    Clear();

    General prop;
    ConvertToCSR(mat, prop, Ptr, IndRow, Val);
    if (!keep_matrix)
      mat.Clear();
   
    FactorizeUnsymmetric();
  }


  template<class T>
  void MatrixWsmp<T>::FactorizeUnsymmetric()
  {
    distributed = false;    
    symmetric = false;
    n = Ptr.GetM()-1;
    
    wsmp_initialize_();
    
    iparm(0) = 0;
    iparm(1) = 0;
    iparm(2) = 0;

    // initialisation step
    int nrhs = 0; T darray;
    CallWgsmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
              NULL, &n, &nrhs, &darray, iparm.GetData(), dparm.GetData());
        
    // CSC/CSR format
    iparm(3) = 0;
    // C-style index
    iparm(4) = 0;
    
    // pivot threshold
    dparm(10) = threshold_pivot;
    
    // ordering step, symbolic step
    iparm(1) = 1;
    iparm(2) = 1;
    CallWgsmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
              NULL, &n, &nrhs, &darray, iparm.GetData(), dparm.GetData());
     
    // then factorization
    iparm(1) = 2;
    iparm(2) = 2;
    CallWgsmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
              NULL, &n, &nrhs, &darray, iparm.GetData(), dparm.GetData());
    
  }
  
  
  template<class T>
  void MatrixWsmp<T>::Solve(Vector<T>& b)
  {
    Solve(SeldonNoTrans, b);
  }
  
  
  template<class T>
  void MatrixWsmp<T>::Solve(const SeldonTranspose& trans, Vector<T>& b)
  {
    Solve(trans, b.GetData(), 1);
  }


  template<class T>
  void MatrixWsmp<T>::Solve(const SeldonTranspose& trans, T* x_ptr, int nrhs)
  {
    int naux = 0, mrp = 0;
    
    // solve
    if (symmetric)
      {
        iparm(1) = 4;
        iparm(2) = 4;
        if (cholesky)
          {
            if (trans.Trans())
              iparm(29) = 2;
            else
              iparm(29) = 1;
          }
        else
          {
            if (refine_solution)
              iparm(2) = 5;
          }
        
        CallWssmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
                  NULL, permut.GetData(), inverse_permut.GetData(), x_ptr, &n, &nrhs,
                  NULL, &naux, &mrp, iparm.GetData(), dparm.GetData());        
      }
    else
      {
        iparm(1) = 3;
        iparm(2) = 3;
        if (refine_solution)
          iparm(2) = 4;
        
        if (trans.Trans())
          {
	    for (int i = 0; i < n*nrhs; i++)
	      x_ptr[i] = conjugate(x_ptr[i]);

            iparm(29) = 4;
          }
        else
          iparm(29) = 0;

        CallWgsmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
                  x_ptr, &n, &nrhs,
                  NULL, iparm.GetData(), dparm.GetData());

        if (trans.Trans())
	  for (int i = 0; i < n*nrhs; i++)
	    x_ptr[i] = conjugate(x_ptr[i]);
      }
  }


  template<class T>
  void MatrixWsmp<T>::Solve(const SeldonTranspose& trans, Matrix<T, General, ColMajor>& b)
  {
    Solve(trans, b.GetData(), b.GetN());
  }


#ifdef SELDON_WITH_MPI
  template<class T>
  void MatrixWsmp<T>::
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
  void MatrixWsmp<T>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
                             Vector<int64_t>& IndRow, Vector<T>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    FactorizeParallel(comm_facto, Ptr, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }
  
  
  template<class T> template<class Tint>
  void MatrixWsmp<T>::  
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
  void MatrixWsmp<T>::
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<int>& Ptr_, Vector<int>& IndRow_,
                    Vector<T>& Val_, const Vector<int>& glob_number,
                    bool sym, bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();
    
    Ptr.SetData(Ptr_);
    IndRow.SetData(IndRow_);
    Val.SetData(Val_);
    Ptr_.Nullify(); IndRow_.Nullify(); Val_.Nullify();

    // size of local system
    n = Ptr.GetM() - 1;
    symmetric = sym;
    distributed = true;
    
    wsmp_initialize_();
    
    // finding the size of the overall system
    int nmax = 0, N = 0;
    for (int i = 0; i < glob_number.GetM(); i++)
      nmax = max(glob_number(i)+1, nmax);

    MPI_Allreduce(&nmax, &N, 1, MPI_INTEGER, MPI_MAX, comm_facto);
    
    int nrhs = 0, naux = 0, mrp = 0;
    if (sym)
      {
        CallWssmp(&n, NULL, NULL, NULL, NULL, NULL, NULL,
                  NULL, &n, &nrhs, NULL, &naux, &mrp, iparm.GetData(), dparm.GetData());      
        
        // CSC/CSR format
        iparm(3) = 0;
        // C-style index
        iparm(4) = 0;
        
        // wsmp can reuse the arrays IndRow and Val
        iparm(13) = 3;
        
        // pivot threshold
        dparm(10) = threshold_pivot;
        
        if (cholesky)
          {
            // L L^T factorization
            iparm(30) = 0;
          }
        else
          {
            // L D L^T factorization
            if (use_pivoting)
              iparm(30) = 2;
            else
              iparm(30) = 1;
            
            if (IsComplexNumber(T()))
              iparm(30) += 2;
          }
        
        permut.Reallocate(N);
        inverse_permut.Reallocate(N);
        
        iparm(1) = 1;
        iparm(2) = 3;
        CallWssmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(), NULL,
                  permut.GetData(), inverse_permut.GetData(), NULL,
                  &n, &nrhs, NULL, &naux, &mrp, iparm.GetData(), dparm.GetData());
      }
    else
      {

        iparm(0) = 0;
        iparm(1) = 0;
        iparm(2) = 0;
        
        // initialisation step
        T darray;
        CallWgsmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
                  NULL, &n, &nrhs, &darray, iparm.GetData(), dparm.GetData());

        // CSC/CSR format
        iparm(3) = 1;
        // C-style index
        iparm(4) = 0;
        
        // pivot threshold
        dparm(10) = threshold_pivot;
        
        // ordering step and factorization
        iparm(1) = 1;
        iparm(2) = 2;
        CallWgsmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
                  NULL, &n, &nrhs, &darray, iparm.GetData(), dparm.GetData());
        
      }    
  }

  
  //! solves A x = b or A^T x = b in parallel
  template<class T>
  void MatrixWsmp<T>::SolveDistributed(MPI_Comm& comm_facto,
                                       const SeldonTranspose& TransA,
                                       Vector<T>& x, const Vector<int>& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), 1, glob_num);
  }

  
  template<class T>
  void MatrixWsmp<T>::SolveDistributed(MPI_Comm& comm_facto,
                                       const SeldonTranspose& TransA,
                                       T* x_ptr, int nrhs, const Vector<int>& glob_num)
  {
    int naux = 0, mrp = 0;
    if (symmetric)
      {
        iparm(1) = 4;
        iparm(2) = 4;
        if (cholesky)
          {
            if (TransA.Trans())
              iparm(29) = 2;
            else
              iparm(29) = 1;
          }
        else
          {
            if (refine_solution)
              iparm(2) = 5;
          }

        CallWssmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(), NULL,
                  permut.GetData(), inverse_permut.GetData(), x_ptr,
                  &n, &nrhs, NULL, &naux, &mrp, iparm.GetData(), dparm.GetData());
      }
    else
      {
        iparm(1) = 3;
        iparm(2) = 3;
        if (refine_solution)
          iparm(2) = 4;
        
        if (TransA.Trans())
          {
	    for (int i = 0; i < n*nrhs; i++)
	      x_ptr[i] = conjugate(x_ptr[i]);
	    
            iparm(29) = 4;
          }
        else
          iparm(29) = 0;

        CallWgsmp(&n, Ptr.GetData(), IndRow.GetData(), Val.GetData(),
                  x_ptr, &n, &nrhs,
                  NULL, iparm.GetData(), dparm.GetData());
        
        if (TransA.Trans())
	  for (int i = 0; i < n*nrhs; i++)
	    x_ptr[i] = conjugate(x_ptr[i]);
      }
  }


  //! solves A x = b or A^T x = b in parallel
  template<class T>
  void MatrixWsmp<T>::SolveDistributed(MPI_Comm& comm_facto,
                                       const SeldonTranspose& TransA,
                                       Matrix<T, General, ColMajor>& x,
                                       const Vector<int>& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), x.GetN(), glob_num);
  }
#endif

  
  //! Factorization of a sparse matrix with Wsmp
  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T, Prop, Storage, Allocator>& A, MatrixWsmp<T>& mat_lu,
	     bool keep_matrix)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }


  //! LU solution with a vector whose type is the same as for Wsmp object
  template<class T, class Allocator>
  void SolveLU(MatrixWsmp<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  //! LU solution with a vector whose type is the same as for Wsmp object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixWsmp<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! LU solution with a matrix whose type is the same as for Wsmp object
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixWsmp<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(SeldonNoTrans, x);
  }


  //! LU solution with a matrix whose type is the same as for Wsmp object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixWsmp<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! Solves A x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(MatrixWsmp<double>& mat_lu,
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
	       MatrixWsmp<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x)
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
  void SolveLU(MatrixWsmp<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixSuperLU<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }


  //! Solves A x = b or A^T x = b, where A is complex and x is real => Forbidden  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixWsmp<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixSuperLU<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }

}

#define SELDON_FILE_WSMP_CXX
#endif
