// Copyright (C) 2001-2011 INRIA, Vivien Mallet
// Author(s): Marc Duruflé, Marc Fragu, Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_CXX

#include "Functions.hxx"

#include "../computation/basic_functions/Functions_Vector.cxx"
#include <list>

/*
  Functions defined in this file:
  (storage RowSparse, ColSparse, RowSymSparse, ColSymSparse, and dense)

  X = A(i, :)
  GetRow(A, i, X)

  X = A(:, j)
  GetCol(A, j, X)

  A(i, :) = X
  SetRow(X, i, A)

  A(:, j) = X
  SetCol(X, j, A)
  

  (dense storages only)
  
  A = A(row_perm, col_perm)
  ApplyPermutation(A, row_perm, col_perm)
  ApplyPermutation(A, row_perm, col_perm, starting_index)

  A(row_perm, col_perm) = A
  ApplyInversePermutation(A, row_perm, col_perm)
  ApplyInversePermutation(A, row_perm, col_perm, starting_index)
  
  A = Drow * A * Dcol
  ScaleMatrix(A, Drow, Dcol)
  
  A = Drow * A
  ScaleLeftMatrix(A, Drow)

  A = A * Dcol
  ScaleRightMatrix(A, Dcol)
  
  Implementations of ApplyPermutation, ApplyInversePermutation,
  ScaleMatrix, ScaleLeftMatrix and ScaleRightMatrix for sparse storages
  are present in file matrix_sparse/Permutation_ScalingMatrix.cxx
*/

namespace Seldon
{

  //! Extracts a row from a sparse matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, RowSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    int size_row = ptr[i+1] - ptr[i];

    X.Reallocate(size_row);
    long shift = ptr[i];
    for (int j = 0; j < size_row; j++)
      {
	X.Index(j) = ind[shift + j];
	X.Value(j) = data[shift + j];
      }
  }


  //! Extracts a row from a sparse matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, ColSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    list<pair<int, T0> > vec;
    for (int j = 0; j < M.GetN(); j++)
      for (long k = ptr[j]; k < ptr[j+1]; k++)
	if (ind[k] == i)
	  vec.push_back(make_pair(j, data[k]));
    
    typename list<pair<int, T0> >::iterator it;
    X.Reallocate(vec.size());
    int j = 0;
    for (it = vec.begin(); it != vec.end(); it++)
      {
	X.Index(j) = it->first;
	X.Value(j) = it->second;
	j++;
      }
  }

  
  //! Extracts a row from a sparse matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, Symmetric, RowSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    list<pair<int, T0> > vec;
    for (int j = 0; j < i; j++)
      for (long k = ptr[j]; k < ptr[j+1]; k++)
	if (ind[k] == i)
	  vec.push_back(make_pair(j, data[k]));
    
    for (long k = ptr[i]; k < ptr[i+1]; k++)
      vec.push_back(make_pair(ind[k], data[k]));
    
    typename list<pair<int, T0> >::iterator it;
    X.Reallocate(vec.size());
    int j = 0;
    for (it = vec.begin(); it != vec.end(); it++)
      {
	X.Index(j) = it->first;
	X.Value(j) = it->second;
	j++;
      }
  }

  
  //! Extracts a row from a sparse matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, Symmetric, ColSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    list<pair<int, T0> > vec;
    for (long k = ptr[i]; k < ptr[i+1]; k++)
      vec.push_back(make_pair(ind[k], data[k]));
    
    for (int j = i+1; j < M.GetM(); j++)
      for (long k = ptr[j]; k < ptr[j+1]; k++)
	if (ind[k] == i)
	  vec.push_back(make_pair(j, data[k]));
    
    typename list<pair<int, T0> >::iterator it;
    X.Reallocate(vec.size());
    int j = 0;
    for (it = vec.begin(); it != vec.end(); it++)
      {
	X.Index(j) = it->first;
	X.Value(j) = it->second;
	j++;
      }
  }
  
  
  //! Extracts a row from a matrix
  /*!
    \param M matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetRow(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	      int i, Vector<T1, Storage1, Allocator1>& X)
  {
    if (Storage0::Sparse)
      throw WrongArgument("GetRow", "Function intended to dense matrices");
    
    X.Reallocate(M.GetN());
    for (int j = 0; j < M.GetN(); j++)
      X(j) = M(i, j);
  }

  
  //! Extracts a column from a sparse matrix
  /*!
    \param M sparse matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, ColSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int n = M.GetN();
    if (j < 0 || j >= n)
      throw WrongIndex("GetCol()",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    int size_col = ptr[j+1] - ptr[j];

    X.Reallocate(size_col);
    long shift = ptr[j];
    for (int i = 0; i < size_col; i++)
      {
	X.Index(i) = ind[shift + i];
	X.Value(i) = data[shift + i];
      }
  }

  
  //! Extracts a column from a matrix
  /*!
    \param M matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, RowSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int n = M.GetN();
    if (j < 0 || j >= n)
      throw WrongIndex("GetCol()",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    int m = M.GetM();

    list<pair<int, T0> > vec;
    for (int i = 0; i < m; i++)
      for (long k = ptr[i]; k < ptr[i+1]; k++)
	if (ind[k] == j)
	  vec.push_back(make_pair(i, data[k]));
    
    typename list<pair<int, T0> >::iterator it;
    X.Reallocate(vec.size());
    int i = 0;
    for (it = vec.begin(); it != vec.end(); it++)
      {
	X.Index(i) = it->first;
	X.Value(i) = it->second;
	i++;
      }
  }

  
  //! Extracts a column from a sparse matrix
  /*!
    \param M sparse matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, Symmetric, ColSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
    // symmetric matrix row = col
    GetRow(M, j, X);
  }
  
  
  //! Extracts a column from a sparse matrix
  /*!
    \param M sparse matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, Symmetric, RowSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
    // symmetric matrix row = col
    GetRow(M, j, X);
  }
  
  
  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M,
              int j, Vector<T1, PETScSeq, Allocator1>& X)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      X.SetBuffer(i, M(i, j));
    X.Flush();
  }


  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
              int j, Vector<T1, PETScPar, Allocator1>& X)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      X.SetBuffer(i, M(i, j));
    X.Flush();
  }


  //! Extracts a column from a matrix
  /*!
    \param M matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	      int j, Vector<T1, Storage1, Allocator1>& X)
  {
    if (Storage0::Sparse)
      throw WrongArgument("GetCol", "Function intended to dense matrices");
    
    X.Reallocate(M.GetM());
    for (int i = 0; i < M.GetM(); i++)
      X(i) = M(i, j);
  }

  
  //! Extracts columns of a matrix.
  /*! Columns [\a begin, \a end[ of \a M_in are returned in \a M_out.
    \param[in] M_in input matrix.
    \param[in] begin first column of \a M_in to extract.
    \param[in] end the last column to be extracted from \a M_in is \a end - 1.
    \param[out] M_out on exit, matrix composed of the columns \a begin to
    \a end - 1 of \a M_in. \a M_out is reallocated if necessary.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, Storage0, Allocator0>& M_in,
	      int begin, int end,
              Matrix<T1, Prop1, Storage1, Allocator1>& M_out)
  {
    M_out.Reallocate(M_in.GetM(), end - begin);
    for (int i = 0; i < M_in.GetM(); i++)
      for (int j = begin, k = 0; j < end; j++, k++)
        M_out(i, k) = M_in(i, j);
  }


  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, Storage0, Allocator0>& M)
  {
    if (Storage0::Sparse)
      throw WrongArgument("SetRow", "Function intended to dense matrices");
    
    for (int j = 0; j < M.GetN(); j++)
      M.Set(i, j, X(j));
  }
  
  
  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScSeq, Allocator1>& X,
              int i, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M)
  {
    for (int j = 0; j < M.GetN(); j++)
      M.SetBuffer(i, j, X(j));
    M.Flush();
  }


  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScPar, Allocator1>& X,
              int i, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M)
  {
    for (int j = 0; j < M.GetN(); j++)
      M.SetBuffer(i, j, X(j));
    M.Flush();
  }


  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, RowSparse, Allocator0>& M)
  {
    int m = M.GetM();
    int n = M.GetN();
    long nnz = M.GetDataSize();
    int Nx = X.GetSize();

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= m)
      throw WrongIndex("SetRow(Vector, int, Matrix<RowSparse>)",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    long *ptr_vector =  M.GetPtr();
    long ptr_i0 =  ptr_vector[i], ptr_i1 = ptr_vector[i + 1];
    int row_size_difference = Nx - ptr_i1 + ptr_i0;

    if (row_size_difference == 0)
      {
	for (int k = 0; k < Nx; k++)
	  M.GetInd()[k + ptr_i0] = X.Index(k);
	for (int k = 0; k < Nx; k++)
	  M.GetData()[k + ptr_i0] = X.Value(k);
	return;
      }

    Vector<int>
      new_ind_vector(nnz + row_size_difference);
    for (long k = 0; k <  ptr_i0; k++)
      new_ind_vector(k) = M.GetInd()[k];
    for (int k = 0; k < Nx; k++)
      new_ind_vector(k + ptr_i0) = X.Index(k);
    for (long k = 0; k < nnz - ptr_i1; k++)
      new_ind_vector(k + ptr_i0 + Nx) =  M.GetInd()[k + ptr_i1];

    Vector<T1, VectFull, Allocator0 >
      new_data_vector(nnz + row_size_difference);
    for (long k = 0; k <  ptr_i0; k++)
      new_data_vector(k) = M.GetData()[k];
    for (int k = 0; k < Nx; k++)
      new_data_vector(k + ptr_i0) = X.Value(k);
    for (long k = 0; k < nnz - ptr_i1; k++)
      new_data_vector(k + ptr_i0 + Nx) =  M.GetData()[k + ptr_i1];

    Vector<long> new_ptr_vector(m + 1);
    for (int j = 0; j < i + 1; j++)
      new_ptr_vector(j) = ptr_vector[j];
    for (int j = i + 1; j < m+1; j++)
      new_ptr_vector(j) =  ptr_vector[j] + row_size_difference;

    M.SetData(m, n, new_data_vector, new_ptr_vector, new_ind_vector);
  }


  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, ColSparse, Allocator0>& M)
  {
    int m = M.GetM();
    int n = M.GetN();
    long nnz = M.GetDataSize();
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= m)
      throw WrongIndex("SetRow(Vector, int, Matrix<ColSparse>)",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif
    
    // we retrieve the current row of matrix M
    Vector<T0, VectSparse, Allocator0> row;
    GetRow(M, i, row);
    
    // counting the new number of new elements
    long new_nnz = nnz + X.GetM() - row.GetM();
    
    // if we have the same pattern for old and new row
    // we change only the values
    bool same_pattern = true;
    if (X.GetM() == row.GetM())
      {
	for (int k = 0; k < X.GetM(); k++)
	  if (X.Index(k) != row.Index(k))
	    same_pattern = false;
      }
    else
      same_pattern = false;
    
    if (same_pattern)
      {
	for (int k = 0; k < X.GetM(); k++)
	  {
	    int j = X.Index(k);
	    for (long k2 = ptr[j]; k2 < ptr[j+1]; k2++)
	      if (ind[k2] == i)
		data[k2] = X.Value(k);
	  }
      }
    else
      {
	// the pattern has to be modified, reallocating a new matrix
	Vector<long> Ptr(n+1); Vector<int> Ind(new_nnz);
	Vector<T0, VectFull, Allocator0> Val(new_nnz);
	
	// loop on first rows
	int kx = 0, kr = 0;
	long nb = 0;
	Ptr(0) = 0;
	for (int j = 0; j < n; j++)
	  {
	    // trying to find the corresponding index for X and row
	    bool valX = false, val_row = false;
	    while ( (kx < X.GetM()) && (X.Index(kx) < j))
	      kx++;
	    
	    if ( (kx < X.GetM()) && (X.Index(kx) == j))
	      valX = true;
	    
	    while ( (kr < row.GetM()) && (row.Index(kr) < j))
	      kr++;

	    if ( (kr < row.GetM()) && (row.Index(kr) == j))
	      val_row = true;
	    
	    // filling matrix until index i
	    long k = ptr[j];
	    while ((k < ptr[j+1]) && (ind[k] < i))
	      {
		Ind(nb) = ind[k];
		Val(nb) = data[k];
		nb++;
		k++;
	      }
	    
	    // if there is a value on X for this index, adding it
	    if (valX)
	      {
		Ind(nb) = i;
		Val(nb) = X.Value(kx);
		nb++;
		kx++;
	      }
	    
	    // if there is a value on row for this index
	    // we skip this value
	    if (val_row)
	      k++;
	    
	    // last values of column
	    while (k < ptr[j+1])
	      {
		Ind(nb) = ind[k];
		Val(nb) = data[k];
		nb++;
		k++;
	      }
	    
	    Ptr(j+1) = nb;
	  }
	
	M.SetData(m, n, Val, Ptr, Ind);
      }
  }
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, RowSymSparse, Allocator0>& M)
  {
    int n = M.GetN();
    long nnz = M.GetDataSize();
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= n)
      throw WrongIndex("SetRow(Vector, int, Matrix<RowSymSparse>)",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    // we retrieve the current row of matrix M
    Vector<T0, VectSparse, Allocator0> row;
    GetRow(M, i, row);
    
    // counting the new number of new elements
    long new_nnz = nnz + X.GetM() - row.GetM();
    
    // if we have the same pattern for old and new row
    // we change only the values
    bool same_pattern = true;
    if (X.GetM() == row.GetM())
      {
	for (int k = 0; k < X.GetM(); k++)
	  if (X.Index(k) != row.Index(k))
	    same_pattern = false;
      }
    else
      same_pattern = false;
    
    if (same_pattern)
      {
	int kdiag = 0;
	for (int k = 0; k < X.GetM(); k++)
	  {
	    int j = X.Index(k);
	    if (j < i)
	      {
		for (long k2 = ptr[j]; k2 < ptr[j+1]; k2++)
		  if (ind[k2] == i)
		    data[k2] = X.Value(k);
		
		kdiag = k+1;
	      }
	    else
	      {
		long k2 = ptr[i] + k - kdiag;
		data[k2] = X.Value(k);
	      }
	  }
      }
    else
      {
	// the pattern has to be modified, reallocating a new matrix
	Vector<long> Ptr(n+1); Vector<int> Ind(new_nnz);
	Vector<T0, VectFull, Allocator0> Val(new_nnz);
	
	// loop on first rows
	int kx = 0, kr = 0;
	long nb = 0;
	Ptr(0) = 0;
	for (int j = 0; j < i; j++)
	  {
	    // trying to find the corresponding index for X and row
	    bool valX = false, val_row = false;
	    while ( (kx < X.GetM()) && (X.Index(kx) < j))
	      kx++;
	    
	    if ( (kx < X.GetM()) && (X.Index(kx) == j))
	      valX = true;
	    
	    while ( (kr < row.GetM()) && (row.Index(kr) < j))
	      kr++;

	    if ( (kr < row.GetM()) && (row.Index(kr) == j))
	      val_row = true;
	    
	    // filling matrix until index i
	    long k = ptr[j];
	    while ((k < ptr[j+1]) && (ind[k] < i))
	      {
		Ind(nb) = ind[k];
		Val(nb) = data[k];
		nb++;
		k++;
	      }
	    
	    // if there is a value on X for this index, adding it
	    if (valX)
	      {
		Ind(nb) = i;
		Val(nb) = X.Value(kx);
		nb++;
		kx++;
	      }
	    
	    // if there is a value on row for this index
	    // we go to the next value
	    if (val_row)
	      k++;
	    
	    // last values of row
	    while (k < ptr[j+1])
	      {
		Ind(nb) = ind[k];
		Val(nb) = data[k];
		nb++;
		k++;
	      }
	    
	    Ptr(j+1) = nb;
	  }
	
	// then changing row i
	while (kx < X.GetM())
	  {
	    Ind(nb) = X.Index(kx);
	    Val(nb) = X.Value(kx);
	    nb++;
	    kx++;
	  }
	
	Ptr(i+1) = nb;
	
	// then last rows of M
	for (int j = i+1; j < n; j++)
	  {
	    for (long k = ptr[j]; k < ptr[j+1]; k++)
	      {
		Ind(nb) = ind[k];
		Val(nb) = data[k];
		nb++;
	      }
	    
	    Ptr(j+1) = nb;
	  }
	
	M.SetData(n, n, Val, Ptr, Ind);
      }
  }
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, ColSymSparse, Allocator0>& M)
  {
    int n = M.GetN();
    long nnz = M.GetDataSize();
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();
    
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= n)
      throw WrongIndex("SetRow(Vector, int, Matrix<RowSymSparse>)",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    // we retrieve the current row of matrix M
    Vector<T0, VectSparse, Allocator0> row;
    GetRow(M, i, row);
    
    // counting the new number of new elements
    long new_nnz = nnz + X.GetM() - row.GetM();
    
    // if we have the same pattern for old and new row
    // we change only the values
    bool same_pattern = true;
    if (X.GetM() == row.GetM())
      {
	for (int k = 0; k < X.GetM(); k++)
	  if (X.Index(k) != row.Index(k))
	    same_pattern = false;
      }
    else
      same_pattern = false;
    
    if (same_pattern)
      {
	for (int k = 0; k < X.GetM(); k++)
	  {
	    int j = X.Index(k);
	    if (j <= i)
	      {
		long k2 = ptr[i] + k;
		data[k2] = X.Value(k);
	      }
	    else
	      {
		for (long k2 = ptr[j]; k2 < ptr[j+1]; k2++)
		  if (ind[k2] == i)
		    data[k2] = X.Value(k);
	      }
	  }
      }
    else
      {
	// the pattern has to be modified, reallocating a new matrix
	Vector<long> Ptr(n+1); Vector<int> Ind(new_nnz);
	Vector<T0, VectFull, Allocator0> Val(new_nnz);
	
	// first columns of M
	Ptr(0) = 0;
	long nb = 0;
	for (int j = 0; j < i; j++)
	  {
	    for (long k = ptr[j]; k < ptr[j+1]; k++)
	      {
		Ind(nb) = ind[k];
		Val(nb) = data[k];
		nb++;
	      }
	    
	    Ptr(j+1) = nb;
	  }
	
	// then changing column i
	int kx = 0;
	while ( (kx < X.GetM()) && (X.Index(kx) <= i) )
	  {
	    Ind(nb) = X.Index(kx);
	    Val(nb) = X.Value(kx);
	    nb++;
	    kx++;
	  }
	
	Ptr(i+1) = nb;
	
	// loop on last columns
	for (int j = i+1; j < n; j++)
	  {
	    // trying to find the corresponding index for X
	    bool valX = false;
	    while ( (kx < X.GetM()) && (X.Index(kx) < j))
	      kx++;
	    
	    if ( (kx < X.GetM()) && (X.Index(kx) == j))
	      valX = true;
	    
	    // filling matrix until index i
	    long k = ptr[j];
	    while ((k < ptr[j+1]) && (ind[k] < i))
	      {
		Ind(nb) = ind[k];
		Val(nb) = data[k];
		nb++;
		k++;
	      }
	    
	    // if there is a value on X for this index, adding it
	    if (valX)
	      {
		Ind(nb) = i;
		Val(nb) = X.Value(kx);
		nb++;
		kx++;
	      }
	    
	    // if there is a value on column j for this index
	    // we go to the next value
	    if ( (k < ptr[j+1]) && (ind[k] == i))
	      k++;
	    
	    // last values of column
	    while (k < ptr[j+1])
	      {
		Ind(nb) = ind[k];
		Val(nb) = data[k];
		nb++;
		k++;
	      }
	    
	    Ptr(j+1) = nb;
	  }	
	
	M.SetData(n, n, Val, Ptr, Ind);
      }
  }
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, RowLoTriang, Allocator0>& M)
  {
    for (int j = 0; j <= i; j++)
      M.Set(i, j, X(j));
  }


  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, RowLoTriangPacked, Allocator0>& M)
  {
    for (int j = 0; j <= i; j++)
      M.Set(i, j, X(j));
  }


  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, ColLoTriang, Allocator0>& M)
  {
    for (int j = 0; j <= i; j++)
      M.Set(i, j, X(j));
  }


  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, ColLoTriangPacked, Allocator0>& M)
  {
    for (int j = 0; j <= i; j++)
      M.Set(i, j, X(j));
  }

  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, RowUpTriang, Allocator0>& M)
  {
    for (int j = i; j < M.GetN(); j++)
      M.Set(i, j, X(j));
  }

  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, RowUpTriangPacked, Allocator0>& M)
  {
    for (int j = i; j < M.GetN(); j++)
      M.Set(i, j, X(j));
  }

  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, ColUpTriang, Allocator0>& M)
  {
    for (int j = i; j < M.GetN(); j++)
      M.Set(i, j, X(j));
  }

  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, ColUpTriangPacked, Allocator0>& M)
  {
    for (int j = i; j < M.GetN(); j++)
      M.Set(i, j, X(j));
  }

  
  //! Sets a column of a matrix
  /*!
    \param[in] X new column \a j of \a M.
    \param j column index
    \param M matrix
    M(:, j) = X
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, Storage0, Allocator0>& M)
  {
    if (Storage0::Sparse)
      throw WrongArgument("SetCol", "Function intended to dense matrices");
    
    for (int i = 0; i < M.GetM(); i++)
      M.Set(i, j, X(i));
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScSeq, Allocator1>& X,
              int j, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      M.SetBuffer(i, j, X(i));
    M.Flush();
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScPar, Allocator1>& X,
              int j, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      M.SetBuffer(i, j, X(i));
    M.Flush();
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectFull, Allocator1>& X,
	      int j, Matrix<T0, General, RowSparse, Allocator0>& M)
  {
    Vector<T1, VectSparse, Allocator1> X_sparse;
    for (int k = 0; k < X.GetLength(); k++)
      {
        T1 value = X(k);
        if (value != T1(0.))
          X_sparse.AddInteraction(k, value);
      }

    SetCol(X_sparse, j, M);
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
    M(:, j) = X
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, RowSparse, Allocator0>& M)
  {
    int m = M.GetM();
    int n = M.GetN();
    long nnz = M.GetDataSize();
    int Nx = X.GetSize();

#ifdef SELDON_CHECK_BOUNDS
    if (j < 0 || j >= n)
      throw WrongIndex("SetCol(Vector, int, Matrix<RowSparse>)",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif

    // The column to be changed.
    Vector<T1, VectSparse, Allocator1> column_j;
    GetCol(M, j, column_j);
    int Ncolumn_j = column_j.GetSize();
    int column_size_difference = Nx - Ncolumn_j;

    // Built a vector indexed with the rows of column_j and X.
    Vector<int, VectSparse> column_j_mask;
    Vector<int> index_j(Ncolumn_j);
    Vector<int> value_j(Ncolumn_j);
    for (int p = 0; p < Ncolumn_j; p++)
      index_j(p) = column_j.Index(p);
    value_j.Fill(-1);
    column_j_mask.SetData(value_j, index_j);
    value_j.Nullify();
    index_j.Nullify();
    Vector<int, VectSparse> X_mask;
    Vector<int> index_x(Nx);
    Vector<int> value_x(Nx);
    for (int p = 0; p < Nx; p++)
      index_x(p) = X.Index(p);
    value_x.Fill(1);
    X_mask.SetData(value_x, index_x);
    value_x.Nullify();
    index_x.Nullify();
    X_mask.AddInteractionRow(column_j_mask.GetSize(),
			     column_j_mask.GetIndex(),
			     column_j_mask.GetData(), true);

    // Built the new pointer vector.
    Vector<long> ptr_vector;
    ptr_vector.SetData(m + 1, M.GetPtr());
    Vector<long> new_ptr_vector(m + 1);
    new_ptr_vector.Zero();
    for (int p = 0; p < X_mask.GetM(); p++)
      new_ptr_vector(X_mask.Index(p) + 1) = X_mask.Value(p);

    for (int p = 0; p < m; p++)
      new_ptr_vector(p + 1) += new_ptr_vector(p);

    new_ptr_vector += ptr_vector;

    // Built the new index and the new data vectors row by row.
    Vector<int>
      new_ind_vector(nnz + column_size_difference);
    Vector<T0, VectFull, Allocator0>
      new_data_vector(nnz + column_size_difference);

    Vector<T0, VectSparse, Allocator0> working_vector;
    int Nworking_vector;

    int line = 0;
    for (int interaction = 0; interaction < X_mask.GetM(); interaction++)
      {
	int ind_x =  X_mask.Index(interaction);
	for (int k = 0; k < ptr_vector(ind_x) -  ptr_vector(line); k++)
	  new_ind_vector.GetData()[k + new_ptr_vector(line)] =
	    M.GetInd()[k + ptr_vector(line)];
	for (int k = 0; k < ptr_vector(ind_x) -  ptr_vector(line); k++)
	  new_data_vector.GetData()[k + new_ptr_vector(line)] =
	    M.GetData()[k + ptr_vector(line)];

	int ind_j;
	Nworking_vector = ptr_vector(ind_x + 1) - ptr_vector(ind_x);
	working_vector.SetData(Nworking_vector,
			       M.GetData() + ptr_vector(ind_x),
			       M.GetInd() + ptr_vector(ind_x));
	switch(X_mask.Value(interaction))
	  {
	    // Collision.
	  case 0:
	    working_vector.Get(j) = X(ind_x);
	    for (int k = 0; k < Nworking_vector; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetIndex()[k];
	    for (int k = 0; k < Nworking_vector; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetData()[k];
	    break;

	    // Suppression.
	  case -1:
	    ind_j = 0;
	    while (ind_j < Nworking_vector &&
		   working_vector.Index(ind_j) != j)
	      ind_j++;

	    for (int k = 0; k < ind_j; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetIndex()[k];
	    for (int k = 0; k < Nworking_vector - ind_j - 1; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x) + ind_j] =
		working_vector.GetIndex()[k + ind_j + 1];

	    for (int k = 0; k < ind_j; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetData()[k];
	    for (int k = 0; k < Nworking_vector - ind_j - 1; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x) + ind_j] =
		working_vector.GetData()[k + ind_j + 1];
	    break;

	    // Addition.
	  case 1:
	    ind_j = 0;
	    while (ind_j < Nworking_vector &&
		   working_vector.Index(ind_j) < j)
	      ind_j++;
	    for (int k = 0; k < ind_j; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetIndex()[k];
	    new_ind_vector.GetData()[new_ptr_vector(ind_x) + ind_j] = j;
	    for (int k = 0; k < Nworking_vector - ind_j; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x) + ind_j + 1]
		= working_vector.GetIndex()[k + ind_j];

	    for (int k = 0; k < ind_j; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetData()[k];
	    new_data_vector.GetData()[new_ptr_vector(ind_x)  + ind_j]
	      = X(ind_x);
	    for (int k = 0; k < Nworking_vector - ind_j; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x) + ind_j + 1]
		= working_vector.GetData()[k + ind_j];
	  }

	line = ind_x + 1;
	working_vector.Nullify();
      }
    for (int k = 0; k < ptr_vector(m) -  ptr_vector(line); k++)
      new_ind_vector.GetData()[k + new_ptr_vector(line)] =
	M.GetInd()[k + ptr_vector(line)];
    for (int k = 0; k < ptr_vector(m) -  ptr_vector(line); k++)
      new_data_vector.GetData()[k + new_ptr_vector(line)] =
	M.GetData()[k + ptr_vector(line)];

    M.SetData(m, n, new_data_vector, new_ptr_vector, new_ind_vector);
    ptr_vector.Nullify();
    new_data_vector.Nullify();
    new_ind_vector.Nullify();
    new_ptr_vector.Nullify();
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, ColSparse, Allocator0>& M)
  {
    int m = M.GetM();
    int n = M.GetN();
    long nnz = M.GetDataSize();
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T0* data = M.GetData();

#ifdef SELDON_CHECK_BOUNDS
    if (j < 0 || j >= n)
      throw WrongIndex("SetCol(Vector, int, Matrix<ColSparse>)",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif
    
    int size_col = ptr[j+1] - ptr[j];
    if (size_col == X.GetM())
      {
	// no need of reallocating matrix M
	// changing indexes and values
	for (int k = 0; k < X.GetM(); k++)
	  {
	    ind[ptr[j] + k] = X.Index(k);
	    data[ptr[j] + k] = X.Value(k);
	  }
      }
    else
      {
	long new_nnz = nnz + X.GetM() - size_col;
	// new matrix
	Vector<long> Ptr(n+1); Vector<int> Ind(new_nnz);
	Vector<T0, VectFull, Allocator0> Val(new_nnz);
	Ptr(0) = 0;
	for (int i = 0; i < j; i++)
	  {
	    for (long k = ptr[i]; k < ptr[i+1]; k++)
	      {
		Ind(k) = ind[k];
		Val(k) = data[k];
	      }
	    Ptr(i+1) = ptr[i+1];
	  }
	
	new_nnz = Ptr(j);
	for (int k = 0; k < X.GetM(); k++)
	  {
	    Ind(new_nnz) = X.Index(k);
	    Val(new_nnz) = X.Value(k);
	    new_nnz++;
	  }
	
	Ptr(j+1) = new_nnz;
	for (int i = j+1; i < n; i++)
	  {
	    for (long k = ptr[i]; k < ptr[i+1]; k++)
	      {
		Ind(new_nnz) = ind[k];
		Val(new_nnz) = data[k];
		new_nnz++;
	      }
	    Ptr(i+1) = new_nnz;
	  }
	
	M.SetData(m, n, Val, Ptr, Ind);
      }
  }
  

  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, RowSymSparse, Allocator0>& M)
  {
    // symmetric matrix, row = col
    SetRow(X, j, M);
  }
  
  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, ColSymSparse, Allocator0>& M)
  {
    // symmetric matrix, row = col
    SetRow(X, j, M);
  }
  
  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, RowLoTriang, Allocator0>& M)
  {
    for (int i = j; i < M.GetM(); i++)
      M.Set(i, j, X(i));
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, RowLoTriangPacked, Allocator0>& M)
  {
    for (int i = j; i < M.GetM(); i++)
      M.Set(i, j, X(i));
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, ColLoTriang, Allocator0>& M)
  {
    for (int i = j; i < M.GetM(); i++)
      M.Set(i, j, X(i));
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, ColLoTriangPacked, Allocator0>& M)
  {
    for (int i = j; i < M.GetM(); i++)
      M.Set(i, j, X(i));
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, RowUpTriang, Allocator0>& M)
  {
    for (int i = 0; i <= j; i++)
      M.Set(i, j, X(i));
  }

  
    //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, RowUpTriangPacked, Allocator0>& M)
  {
    for (int i = 0; i <= j; i++)
      M.Set(i, j, X(i));
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, ColUpTriang, Allocator0>& M)
  {
    for (int i = 0; i <= j; i++)
      M.Set(i, j, X(i));
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, ColUpTriangPacked, Allocator0>& M)
  {
    for (int i = 0; i <= j; i++)
      M.Set(i, j, X(i));
  }

  
  //! Permutation of a general matrix stored by rows.
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowMajor, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, RowMajor, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(i, j) = A_copy(row_perm(i) - starting_index,
                         col_perm(j) - starting_index);
  }


  //! Permutation of a general matrix stored by columns.
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColMajor, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, ColMajor, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i < A.GetM(); i++)
        A(i, j) = A_copy(row_perm(i) - starting_index,
                         col_perm(j) - starting_index);
  }

  
  //! Permutation of a symmetric matrix stored by rows.
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowSymPacked, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, RowSymPacked, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A(i, j) = A_copy(row_perm(i) - starting_index,
                         row_perm(j) - starting_index);
  }


  //! Permutation of a symmetric matrix stored by columns.
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColSymPacked, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, ColSymPacked, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A(i, j) = A_copy(row_perm(i) - starting_index,
                         row_perm(j) - starting_index);
  }

  
  //! Permutation of a symmetric matrix stored by rows.
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowSym, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, RowSym, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) = A_copy(row_perm(i) - starting_index,
                             row_perm(j) - starting_index);
  }


  //! Permutation of a symmetric matrix stored by columns.
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColSym, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, ColSym, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Val(i, j) = A_copy(row_perm(i) - starting_index,
                             row_perm(j) - starting_index);
  }

  
  //! Permutation of an hermitian matrix stored by rows.
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowHermPacked, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, RowHermPacked, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) = A_copy(row_perm(i) - starting_index,
                             row_perm(j) - starting_index);
  }


  //! Permutation of an hermitian matrix stored by columns.
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColHermPacked, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, ColHermPacked, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Val(i, j) = A_copy(row_perm(i) - starting_index,
                             row_perm(j) - starting_index);
  }

  
  //! Permutation of an hermitian matrix stored by rows.
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowHerm, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, RowHerm, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) = A_copy(row_perm(i) - starting_index,
                             row_perm(j) - starting_index);
  }


  //! Permutation of an hermitian matrix stored by columns.
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColHerm, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, ColHerm, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Val(i, j) = A_copy(row_perm(i) - starting_index,
                             row_perm(j) - starting_index);
  }
  
  
  //! Inverse permutation of a general matrix stored by rows.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowMajor, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, RowMajor, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(row_perm(i) - starting_index, col_perm(j) - starting_index)
          = A_copy(i, j);
  }


  //! Inverse permutation of a general matrix stored by columns.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColMajor, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, ColMajor, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i < A.GetM(); i++)
        A(row_perm(i) - starting_index, col_perm(j) - starting_index)
          = A_copy(i, j);
  }

  
  //! Inverse permutation of a symmetric matrix stored by rows.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowSymPacked, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, RowSymPacked, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Set(row_perm(i) - starting_index, row_perm(j) - starting_index,
              A_copy(i, j));
  }


  //! Inverse permutation of a symmetric matrix stored by columns.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColSymPacked, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, ColSymPacked, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Set(row_perm(i) - starting_index, row_perm(j) - starting_index,
              A_copy(i, j));
  }
  
  
  //! Inverse permutation of a symmetric matrix stored by rows.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowSym, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, RowSym, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Set(row_perm(i) - starting_index, row_perm(j) - starting_index,
              A_copy(i, j));
  }


  //! Inverse permutation of a symmetric matrix stored by columns.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColSym, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, ColSym, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Set(row_perm(i) - starting_index, row_perm(j) - starting_index,
              A_copy(i, j));
  }
  
  
  //! Inverse permutation of an hermitian matrix stored by rows.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowHermPacked, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, RowHermPacked, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Set(row_perm(i) - starting_index, row_perm(j) - starting_index,
              A_copy(i, j));
  }


  //! Inverse permutation of an hermitian matrix stored by columns.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColHermPacked, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, ColHermPacked, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Set(row_perm(i) - starting_index, row_perm(j) - starting_index,
              A_copy(i, j));
  }
  
  
  //! Inverse permutation of an hermitian matrix stored by rows.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowHerm, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, RowHerm, Allocator> A_copy; A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Set(row_perm(i) - starting_index, row_perm(j) - starting_index,
              A_copy(i, j));
  }


  //! Inverse permutation of an hermitian matrix stored by columns.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColHerm, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, ColHerm, Allocator> A_copy; A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Set(row_perm(i) - starting_index, row_perm(j) - starting_index,
              A_copy(i, j));
  }
  
  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowMajor, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(i, j) *= Drow(i)*Dcol(j);
  }


  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColMajor, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i < A.GetM(); i++)
        A(i, j) *= Drow(i)*Dcol(j);
  }
  
  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Drow
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowSymPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A(i, j) *= Drow(i)*Drow(j);
  }


  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Drow
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColSymPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A(i, j) *= Drow(i)*Drow(j);
  }
  
  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Drow
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowSym, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i)*Drow(j);
  }


  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Drow
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColSym, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Val(i, j) *= Drow(i)*Drow(j);
  }
  
  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Drow
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowHermPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i)*Drow(j);
  }


  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Drow
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColHermPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Val(i, j) *= Drow(i)*Drow(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Drow
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowHerm, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i)*Drow(j);
  }


  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Drow
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColHerm, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i <= j; i++)
        A.Val(i, j) *= Drow(i)*Drow(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowLoTriangPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Drow(i)*Dcol(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowLoTriang, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Drow(i)*Dcol(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColLoTriangPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Drow(i)*Dcol(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColLoTriang, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Drow(i)*Dcol(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowUpTriangPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetM(); j++)
        A.Val(i, j) *= Drow(i)*Dcol(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowUpTriang, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i)*Dcol(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColUpTriangPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetM(); j++)
        A.Val(i, j) *= Drow(i)*Dcol(j);
  }

  
  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColUpTriang, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i)*Dcol(j);
  }

  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowMajor, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(i, j) *= Drow(i);
  }
    

  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColMajor, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(i, j) *= Drow(i);
  }

  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowLoTriangPacked, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Drow(i);
  }

  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowLoTriang, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Drow(i);
  }

  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColLoTriangPacked, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Drow(i);
  }

  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColLoTriang, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Drow(i);
  }
  
  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowUpTriangPacked, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i);
  }

  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowUpTriang, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i);
  }
  
  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColUpTriangPacked, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i);
  }

  
  //! Left-scaling of a matrix
  /*!
    A is replaced by Drow A
    where Drow is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColUpTriang, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Drow(i);
  }
  
  
  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowMajor, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(i, j) *= Dcol(j);
  }

  
  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColMajor, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(i, j) *= Dcol(j);
  }
  
  
  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowLoTriangPacked, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Dcol(j);
  }


  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowLoTriang, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Dcol(j);
  }


  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColLoTriangPacked, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Dcol(j);
  }


  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColLoTriang, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j <= i; j++)
        A.Val(i, j) *= Dcol(j);
  }


  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowUpTriangPacked, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Dcol(j);
  }


  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowUpTriang, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Dcol(j);
  }


  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColUpTriangPacked, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Dcol(j);
  }


  //! Right-scaling of a matrix
  /*!
    A is replaced by A Dcol
    where Dcol is diagonal and stored as a dense vector
  */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColUpTriang, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = i; j < A.GetN(); j++)
        A.Val(i, j) *= Dcol(j);
  }

} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_CXX
#endif
