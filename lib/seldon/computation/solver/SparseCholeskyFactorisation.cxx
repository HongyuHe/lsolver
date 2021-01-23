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

#ifndef SELDON_FILE_SPARSE_CHOLESKY_FACTORISATION_CXX


#include "SparseCholeskyFactorisation.hxx"

namespace Seldon
{

  //! Implementation of Cholesky factorization for sparse symmetric matrix
  /*!
    This method may be slow for large matrices. For large matrices,
    it is more efficient to use an external library (Cholmod for example).
    No renumbering is performed by this method.
   */
  template<class T, class Prop, class Allocator>
  void GetCholesky(Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
                   int print_level)
  {
    int n = A.GetN();
    T t, s, fact;
    int j_col, jrow, index_lu, jpos;    
    Vector<T> Row_Val(n);
    IVect Index(n), Row_Ind(n);
    Row_Val.Fill(0);
    Row_Ind.Fill(-1);
    
    Index.Fill(-1);
    
    // conversion to unsymmetric matrix 
    Matrix<T, General, ArrayRowSparse, Allocator> B;
    Copy(A, B);
    
    // A is cleared
    A.Clear();
    A.Reallocate(n, n);

    // Main loop over rows.
    T zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    int new_percent = 0, old_percent = 0;
    for (int i_row = 0; i_row < n; i_row++)
      {
        // Progress bar if print level is high enough.
        if (print_level > 0)
          {
            new_percent = int(double(i_row) / double(n-1) * 78.);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#";
                cout.flush();
              }

            old_percent = new_percent;
          }
        
	int size_row = B.GetRowSize(i_row);
        
        // we are separating lower from upper part
	int length_lower = 0, length_upper = 1; 
        Row_Ind(i_row) = i_row;
	Row_Val(i_row) = zero;
	Index(i_row) = i_row;
	
        for (int j = 0; j < size_row; j++) 
          {
            int k = B.Index(i_row, j);
            t = B.Value(i_row, j);
            if (k < i_row)
              {
                Row_Ind(length_lower) = k; 
                Row_Val(length_lower) = t;
                Index(k) = length_lower;
                length_lower++;
              }
            else if (k == i_row)
              Row_Val(i_row) = t;
            else
              {
                jpos = i_row + length_upper;
		Row_Ind(jpos) = k;
		Row_Val(jpos) = t;
		Index(k) = jpos;
		length_upper++;
              }
          }
        
        B.ClearRow(i_row);
        
	j_col = 0;
	int length = 0;
	
	// previous rows are eliminated	
	while (j_col < length_lower) 
	  {            
	    jrow = Row_Ind(j_col);
            
            // we search first element in lower part
            int k = j_col;
            for (int j = (j_col+1) ; j < length_lower; j++) 
	      {
		if (Row_Ind(j) < jrow) 
		  {
		    jrow = Row_Ind(j);
		    k = j;
		  }
	      }
            
            
	    if (k != j_col) 
	      {
                // if k different from j_col, we are exchanging positions
		int j = Row_Ind(j_col);
		Row_Ind(j_col) = Row_Ind(k);
		Row_Ind(k) = j;
                
		Index(jrow) = j_col;
		Index(j) = k;
		
		s = Row_Val(j_col);
		Row_Val(j_col) = Row_Val(k);
	        Row_Val(k) = s;
	      }
            
	    // zero out element in row
	    Index(jrow) = -1;
	    fact = Row_Val(j_col) * A.Value(jrow, 0);
	    
	    // combine current row and row jrow
	    for (int k = 1; k < A.GetRowSize(jrow); k++) 
	      {
		s = fact * A.Value(jrow, k); 
		int j = A.Index(jrow, k);
		
		jpos = Index(j);
		if (j >= i_row)
		  {
                    // dealing with upper part		    
		    if (jpos == -1)
		      {			
			// this is a fill-in element
                        int i = i_row + length_upper;
			Row_Ind(i) = j;
			Index(j) = i;
			Row_Val(i) = -s;
			length_upper++;
		      } 
		    else 
		      {
                        // this is not a fill-in element
			Row_Val(jpos) -= s;
		      }
		  } 
		else 
		  {		    
		    // dealing  with lower part.		    
		    if (jpos == -1) 
		      {
			//  this is a fill-in element
			Row_Ind(length_lower) = j;
			Index(j) = length_lower;
			Row_Val(length_lower) = -s;
			length_lower++;
		      } 
		    else 
		      {
			// this is not a fill-in element
			Row_Val(jpos) -= s;
		      }
		  }
	      }
            
	    // store this pivot element
            // (from left to right -- no danger of overlap
            // with the working elements in L (pivots).	    
	    Row_Val(length) = fact;
	    Row_Ind(length) = jrow;
	    ++length;
	    j_col++;
	  }
        
        for (int k = 0; k < length_upper; k++) 
	  Index(Row_Ind(i_row+k )) = -1;
        	
	// now we can store the uppert part of row
	size_row = length_upper;
	A.ReallocateRow(i_row, length_upper);

	// We store inverse of square root of diagonal element of u.
        if (Row_Val(i_row) <= 0)
          {
            cout << "Error during Cholesky factorization " << endl;
            cout << "Matrix must be definite positive " << endl;
            cout << "but diagonal element of row " << i_row 
                 << "is equal to " << Row_Val(i_row) << endl;

	    throw WrongArgument("GetCholesky", "Negative diagonal");
          }
        
	A.Value(i_row, 0) = one / Row_Val(i_row);
        A.Index(i_row, 0) = i_row;
        index_lu = 1;
	
	// and extra-diagonal terms
	for (int k = (i_row+1); k < (i_row+length_upper); k++) 
	  {
	    A.Index(i_row, index_lu) = Row_Ind(k);
	    A.Value(i_row, index_lu) = Row_Val(k);
	    index_lu++;
	  }
      }
    
    if (print_level > 0)
      {
        cout << endl;
        cout << "The matrix takes " <<
          int((A.GetDataSize() * (sizeof(T)+4)) / (1024*1024)) << " MB" << endl;
      }

    // Diagonal of A is replaced by its square root.
    for (int i = 0; i < n; i++)
      {
        A.Value(i, 0) = sqrt(A.Value(i,0));
        // and other elements multiplied by this value
        for (int k = 1; k < A.GetRowSize(i); k++)
          A.Value(i, k) *= A.Value(i, 0);
	
	A.Value(i, 0) = one / A.Value(i, 0);
      }
  }
  

  //! Resolution of L x = y or L^T x = y.
  /*!
    \param[in] TransA SeldonTrans or SeldonNoTrans.
    \param[in] A Cholesky factorization obtained after calling "GetCholesky".
    \param[in,out] x on exit, it is overwritten by the solution.
   */
  template<class T0, class Prop, class Allocator0,
           class T1, class Storage, class Allocator1>
  void SolveCholesky(const SeldonTranspose& TransA,
                     const Matrix<T0, Prop, ArrayRowSymSparse, Allocator0>& A,
                     Vector<T1, Storage, Allocator1>& x)
  {
    int n = A.GetM();
    if (n <= 0)
      return;
    
    if (TransA.Trans())
      {
        // we solve L^T x = x
        int j;
        for (int i = n-1; i >= 0; i--) 
          {
            T1 val = x(i);
            for (int k = 1; k < A.GetRowSize(i) ; k++)
              {
                j = A.Index(i, k);
                val -= A.Value(i, k)*x(j);
              }

            x(i) = val / A.Value(i, 0);
          }
      }
    else
      {
        // we solve L x = x
        int j;
        for (int i = 0; i < n; i++) 
          {
            x(i) /= A.Value(i, 0);
            for (int k = 1; k < A.GetRowSize(i) ; k++)
              {
                j = A.Index(i, k);
                x(j) -= A.Value(i, k)*x(i);
              }
          }
      }
  }

  
  //! Resolution of L x = y or L^T x = y.
  /*!
    \param[in] TransA SeldonTrans or SeldonNoTrans.
    \param[in] A Cholesky factorization obtained after calling "GetCholesky".
    \param[in,out] x on exit, it is overwritten by the solution.
   */
  template<class T0, class Prop, class Alloc0,
	   class T1, class Storage, class Allocator1>
  void SolveCholesky(const SeldonTranspose& TransA,
		     const Matrix<T0, Prop, RowSymSparse, Alloc0>& A,
                     Vector<T1, Storage, Allocator1>& X)
  {
    int n = A.GetM();
    if (n <= 0)
      return;
    
    int* ind = A.GetInd();
    long* ptr = A.GetPtr();
    T0* data = A.GetData();
    T1 val;
    
    if (TransA.Trans())
      {
	// Resolution of L^T x = x.
	for (int i = n - 1; i >= 0; i--)
	  {
	    val = X(i);
	    for (long j = ptr[i] + 1; j < ptr[i+1]; j++)
	      val -= data[j] * X(ind[j]);
	    
	    X(i) = val / data[ptr[i]];
	  }
      }
    else
      {
	// Resolution of L x = x.
	for (int i = 0; i < n; i++)
	  {
	    X(i) /= data[ptr[i]];
	    for (long j = ptr[i] + 1; j < ptr[i+1]; j++)
	      X(ind[j]) -= data[j] * X(i);
	  }
      }
  }


  //! Computation of y = L x or y = L^T x.
  /*!
    \param[in] TransA SeldonTrans or SeldonNoTrans.
    \param[in] A Cholesky factorization obtained after calling "GetCholesky".
    \param[in,out] x on exit, it is overwritten by the value of y.
   */
  template<class T0, class Prop, class Allocator0,
           class T1, class Storage, class Allocator1>
  void MltCholesky(const SeldonTranspose& TransA,
                   const Matrix<T0, Prop, ArrayRowSymSparse, Allocator0>& A,
                   Vector<T1, Storage, Allocator1>& x)
  {
    int n = A.GetM();
    if (n <= 0)
      return;
    
    if (TransA.Trans())
      {
        // we overwrite x by L^T x
        int j;
        for (int i = 0; i < n; i++) 
          {
            T1 val = x(i) * A.Value(i, 0);
            for (int k = 1; k < A.GetRowSize(i) ; k++)
              {
                j = A.Index(i, k);
                val += A.Value(i, k)*x(j);
              }
            
            x(i) = val;
          }
      }
    else
      {
        // we overwrite x by L x
        int j;
        for (int i = n-1; i >= 0; i--) 
          {
            for (int k = 1; k < A.GetRowSize(i) ; k++)
              {
                j = A.Index(i, k);
                x(j) += A.Value(i, k)*x(i);
              }
	    
            x(i) *= A.Value(i, 0);
          }
      }
  }
  
  
  //! Computation of y = L x or y = L^T x.
  /*!
    \param[in] TransA SeldonTrans or SeldonNoTrans.
    \param[in] A Cholesky factorization obtained after calling "GetCholesky".
    \param[in,out] x on exit, it is overwritten by the value of y.
   */
  template<class T0, class Prop, class Allocator0,
           class T1, class Storage, class Allocator1>
  void MltCholesky(const SeldonTranspose& TransA,
                   const Matrix<T0, Prop, RowSymSparse, Allocator0>& A,
                   Vector<T1, Storage, Allocator1>& x)
  {
    int n = A.GetM();
    if (n <= 0)
      return;
    
    int* ind = A.GetInd();
    long* ptr = A.GetPtr();
    T0* data = A.GetData();
    T1 val;
    
    if (TransA.Trans())
      {
        // We overwrite x by L^T x
        for (int i = 0; i < n; i++) 
          {
            val = x(i) * data[ptr[i]];
            for (long k = ptr[i] + 1; k < ptr[i+1]; k++)
              val += data[k] * x(ind[k]);
	    
            x(i) = val;
          }
      }
    else
      {
        // We overwrite x by L x.
	for (int i = n-1; i >= 0; i--) 
          {
            for (long k = ptr[i] + 1; k < ptr[i+1]; k++)
	      x(ind[k]) += data[k] * x(i);
	    
            x(i) *= data[ptr[i]];
          }
      }
  }

  
  //! Default constructor.
  template<class T>
  SparseCholeskySolver<T>::SparseCholeskySolver()
  {
    n = 0;
    print_level = -1;
    type_ordering = SparseMatrixOrdering::IDENTITY;

    type_solver = SELDON_SOLVER;
#ifdef SELDON_WITH_CHOLMOD
    type_solver = CHOLMOD;
#endif
#ifdef SELDON_WITH_PASTIX
    type_solver = PASTIX;
#endif
    
    solver = NULL;
    InitSolver();
  }
  
  
  //! Displays no messages.
  template<class T>
  void SparseCholeskySolver<T>::HideMessages()
  {
    print_level = -1;
    if (solver != NULL)
      solver->HideMessages();
  }
  
  
  //! Displays only brief messages.
  template<class T>
  void SparseCholeskySolver<T>::ShowMessages()
  {
    print_level = 1;
    if (solver != NULL)
      solver->ShowMessages();
  }
  
  
  //! Displays a lot of messages.
  template<class T>
  void SparseCholeskySolver<T>::ShowFullHistory()
  {
    print_level = 3;
    if (solver != NULL)
      solver->ShowMessages();
  }
    
  
  //! Clears Cholesky factors.
  template<class T>
  void SparseCholeskySolver<T>::Clear()
  {
    if (n > 0)
      {
	n = 0;
	if (solver != NULL)
	  solver->Clear();

        mat_sym.Clear();
      }    
  }
    
  
  //! Returns the number of rows.
  template<class T>
  int SparseCholeskySolver<T>::GetM() const
  {
    return n;
  }
  
  
  //! Returns the number of rows.
  template<class T>
  int SparseCholeskySolver<T>::GetN() const
  {
    return n;
  }
  
  
  //! returns memory size used by the object in bytes
  template<class T>
  size_t SparseCholeskySolver<T>::GetMemorySize() const
  {
    size_t taille = mat_sym.GetMemorySize();
    taille += xtmp.GetMemorySize();
    if (solver != NULL)
      taille += solver->GetMemorySize();

    return taille;
  }
  
    
  //! Returns the type of ordering used.
  template<class T>
  int SparseCholeskySolver<T>::GetTypeOrdering() const
  {
    return type_ordering;
  }
  
  
  //! Modifies the ordering used.
  template<class T>
  void SparseCholeskySolver<T>::SetPermutation(const IVect& num)
  {
    type_ordering = SparseMatrixOrdering::USER;
    permutation = num;
  }
  
  
  //! Modifies the type of ordering used.
  template<class T>
  void SparseCholeskySolver<T>::SelectOrdering(int type)
  {
    type_ordering = type;
  }
  
  
  //! Modifies the direct solver used.
  template<class T>
  void SparseCholeskySolver<T>::SelectDirectSolver(int type)
  {
    type_solver = type;
    InitSolver();
  }
  
  
  //! returns the type of direct solver used
  template<class T>
  int SparseCholeskySolver<T>::GetDirectSolver()
  {
    return type_solver;
  }


  template<class T>
  void SparseCholeskySolver<T>::InitSolver()
  {
    if (solver != NULL)
      delete solver;
    
    if (type_solver == CHOLMOD)
      {
#ifdef SELDON_WITH_CHOLMOD
	solver = new MatrixCholmod();
#else
	cout << "Seldon compiled without Cholmod support" << endl;
	cout << "SELDON_WITH_CHOLMOD is not defined" << endl;
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	solver = new MatrixPastix<T>();
#else
	cout << "Seldon compiled without Pastix support" << endl;
	cout << "SELDON_WITH_PASTIX is not defined" << endl;
#endif
      }
    else
      {
	// Seldon solver stored in mat_sym
	solver = NULL;	
      }
  }
    
  
  //! Performs Cholesky factorization.
  template<class T>
  template<class Prop, class Storage, class Allocator>
  void SparseCholeskySolver<T>::
  Factorize(Matrix<T, Prop, Storage, Allocator>& A, bool keep_matrix)
  {    
    n = A.GetM();    
    if (type_solver == CHOLMOD)
      {
#ifdef SELDON_WITH_CHOLMOD
        if (print_level >= 1)
          cout << "Calling Cholmod to factorize the matrix" << endl;
	
	MatrixCholmod& mat_chol = 
	  dynamic_cast<MatrixCholmod&>(*this->solver);
	
	mat_chol.FactorizeMatrix(A, keep_matrix);
#else
	throw Error("SparseCholeskySolver::Factorize",
                    "Recompile with Cholmod or change solver type.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
        if (print_level >= 1)
          cout << "Calling Pastix to factorize the matrix" << endl;

	MatrixPastix<T>& mat_pastix = 
	  dynamic_cast<MatrixPastix<T>&>(*this->solver);
        
        GetCholesky(A, mat_pastix, keep_matrix);
#else
	throw Error("SparseCholeskySolver::Factorize",
                    "Recompile with Pastix or change solver type.");
#endif
      }
    else
      {
        FindSparseOrdering(A, permutation, type_ordering);
        Copy(A, mat_sym);
        if (!keep_matrix)
          A.Clear();
        
        ApplyInversePermutation(mat_sym, permutation, permutation);
        
	GetCholesky(mat_sym, print_level);
        xtmp.Reallocate(n);
      }
  }
   
  
  //! Solves L x = b or L^T x = b.
  template<class T> template<class T1>
  void SparseCholeskySolver<T>
  ::Solve(const SeldonTranspose& TransA, Vector<T1>& x_solution, bool assemble)
  {
    if (type_solver == CHOLMOD)
      {
#ifdef SELDON_WITH_CHOLMOD
	MatrixCholmod& mat_chol = 
	  dynamic_cast<MatrixCholmod&>(*this->solver);
	
	mat_chol.Solve(TransA, x_solution);
#else
	throw Error("SparseCholeskySolver::Solve",
                    "Recompile with Cholmod or change solver type.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	MatrixPastix<T>& mat_pastix = 
	  dynamic_cast<MatrixPastix<T>&>(*this->solver);

        SolveCholesky(TransA, mat_pastix, x_solution);
#else
	throw Error("SparseCholeskySolver::Solve",
                    "Recompile with Pastix or change solver type.");
#endif
      }
    else
      {
	if (TransA.NoTrans())
          {
            for (int i = 0; i < x_solution.GetM(); i++)
              xtmp(permutation(i)) = x_solution(i);
            
            SolveCholesky(TransA, mat_sym, xtmp);
            Copy(xtmp, x_solution);
          }
        else
          {	
            Copy(x_solution, xtmp);
            SolveCholesky(TransA, mat_sym, xtmp);
            
            for (int i = 0; i < x_solution.GetM(); i++)
              x_solution(i) = xtmp(permutation(i));
          }
      }
  }
  
  
  //! Computes L x or L^T.
  template<class T> template<class T1>
  void SparseCholeskySolver<T>
  ::Mlt(const SeldonTranspose& TransA, Vector<T1>& x_solution, bool assemble)
  {
    if (type_solver == CHOLMOD)
      {
#ifdef SELDON_WITH_CHOLMOD
	MatrixCholmod& mat_chol = 
	  dynamic_cast<MatrixCholmod&>(*this->solver);

	mat_chol.Mlt(TransA, x_solution);
#else
	throw Error("SparseCholeskySolver::Mlt",
                    "Recompile with Cholmod or change solver type.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	MatrixPastix<T>& mat_pastix = 
	  dynamic_cast<MatrixPastix<T>&>(*this->solver);

        MltCholesky(TransA, mat_pastix, x_solution);
#else
	throw Error("SparseCholeskySolver::Mlt",
                    "Recompile with Pastix or change solver type.");
#endif
      }
    else
      {
	if (TransA.NoTrans())
          {
            Copy(x_solution, xtmp);
            MltCholesky(TransA, mat_sym, xtmp);
            
            for (int i = 0; i < x_solution.GetM(); i++)
              x_solution(i) = xtmp(permutation(i));
          }
        else
          {
            for (int i = 0; i < x_solution.GetM(); i++)
              xtmp(permutation(i)) = x_solution(i);
            
            MltCholesky(TransA, mat_sym, xtmp);
            Copy(xtmp, x_solution);
          }
      }
  }

} // namespace Seldon.


#define SELDON_FILE_SPARSE_CHOLESKY_FACTORISATION_CXX
#endif
