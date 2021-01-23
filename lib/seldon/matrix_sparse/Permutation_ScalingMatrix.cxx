// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_PERMUTATION_SCALING_MATRIX_CXX

/*
  Functions defined in this file:

  A(I, J) = A
  ApplyInversePermutation(A, I, J)
  
  A = A(I, J)
  ApplyPermutation(A, I, J)
  
  A = Drow * A * Dcol
  ScaleMatrix(A, Drow, Dcol)

  A = Drow * A
  ScaleLeftMatrix(A, Drow)

  A = A * Dcol
  ScaleRightMatrix(A, Dcol)
*/

namespace Seldon
{

  /////////////////////////////
  // ApplyInversePermutation //
  
  
  //! Permutation of a general matrix stored by rows.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowSparse, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm)
  {
    int i, j, k; long l, nnz;
    int m = A.GetM();
    int n = A.GetN();
    int* ind = A.GetInd();
    long* ptr = A.GetPtr();
    T* data = A.GetData();

    /*** Permutation of the columns ***/

    Vector<T, VectFull, Allocator> row_data;
    // Column indexes of a given row.
    Vector<int> row_index;
    for (i = 0; i < m; i++)
      if ((nnz = ptr[i + 1] - ptr[i]) != 0)
        {
          row_index.Reallocate(nnz);
          row_data.Reallocate(nnz);
          for (k = 0, l = ptr[i]; k < nnz; k++, l++)
            {
              // Applies the permutation on the columns.
              row_index(k) = col_perm(ind[l]);
              row_data(k) = data[l];
            }

          // The columns must remain in increasing order.
          Sort(row_index, row_data);

          // Putting back the data into the array.
          for (k = 0, l = ptr[i]; k < nnz; k++, l++)
            {
              ind[l] = row_index(k);
              data[l] = row_data(k);
            }
        }
    
    row_index.Clear();
    row_data.Clear();

    /*** Perturbation of the rows ***/

    // Total number of non-zero elements.
    nnz = ptr[m];

    // 'row_perm' is const, so it must be copied.
    Vector<int> row_permutation(row_perm);
    // Row indexes in the origin matrix: prev_row_index(i) should be the
    // location of the i-th row (from the permuted matrix) in the matrix
    // before permutation.
    Vector<int> prev_row_index(m);
    prev_row_index.Fill();

    Sort(row_permutation, prev_row_index);
    row_permutation.Clear();

    // Description of the matrix after permutations.
    Vector<long> new_ptr(m + 1);
    Vector<int> new_ind(nnz);
    Vector<T, VectFull, Allocator> new_data(nnz);

    long ptr_count = 0; int length;
    for (i = 0; i < m; i++)
      {
        length = ptr[prev_row_index(i) + 1] - ptr[prev_row_index(i)];
        for (j = 0; j < length; j++)
          {
            new_data(ptr_count + j) = data[ptr[prev_row_index(i)] + j];
            new_ind(ptr_count + j) = ind[ptr[prev_row_index(i)] + j];
          }
        new_ptr(i) = ptr_count;
        ptr_count += length;
      }
    new_ptr(m) = ptr_count;

    A.SetData(m, n, new_data, new_ptr, new_ind);
  }


  //! Permutation of a general matrix stored by columns.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColSparse, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm)
  {
    int i, j, k; long l, nnz;
    int m = A.GetM();
    int n = A.GetN();
    int* ind = A.GetInd();
    long* ptr = A.GetPtr();
    T* data = A.GetData();

    /*** Permutation of the rows ***/

    Vector<T, VectFull, Allocator> col_data;
    // Row indexes of a given column.
    Vector<int> col_index;
    for (i = 0; i < n; i++)
      if ((nnz = ptr[i + 1] - ptr[i]) != 0)
        {
          col_index.Reallocate(nnz);
          col_data.Reallocate(nnz);
          for (k = 0, l = ptr[i]; k < nnz; k++, l++)
            {
              // Applies the permutation on the rows.
              col_index(k) = row_perm(ind[l]);
              col_data(k) = data[l];
            }

          // The rows must remain in increasing order.
          Sort(col_index, col_data);

          // Putting back the data into the array.
          for (k = 0, l = ptr[i]; k < nnz; k++, l++)
            {
              ind[l] = col_index(k);
              data[l] = col_data(k);
            }
        }
    col_index.Clear();
    col_data.Clear();

    /*** Perturbation of the columns ***/

    // Total number of non-zero elements.
    nnz = ptr[n];

    // 'col_perm' is const, so it must be copied.
    Vector<int> col_permutation(col_perm);
    // Column indexes in the origin matrix: prev_col_index(i) should be the
    // location of the i-th column (from the permuted matrix) in the matrix
    // before permutation.
    Vector<int> prev_col_index(n);
    prev_col_index.Fill();

    Sort(col_permutation, prev_col_index);
    col_permutation.Clear();

    // Description of the matrix after permutations.
    Vector<long> new_ptr(n + 1);
    Vector<int> new_ind(nnz);
    Vector<T, VectFull, Allocator> new_data(nnz);

    int ptr_count = 0, length;
    for (i = 0; i < n; i++)
      {
        length = ptr[prev_col_index(i) + 1] - ptr[prev_col_index(i)];
        for (j = 0; j < length; j++)
          {
            new_data(ptr_count + j) = data[ptr[prev_col_index(i)] + j];
            new_ind(ptr_count + j) = ind[ptr[prev_col_index(i)] + j];
          }
        new_ptr(i) = ptr_count;
        ptr_count += length;
      }
    new_ptr(n) = ptr_count;

    A.SetData(m, n, new_data, new_ptr, new_ind);
  }

  
  //! Permutation of a symmetric matrix stored by rows.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowSymSparse, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm)
  {
    int m = A.GetM(), n = A.GetN();
    long nnz = A.GetDataSize();
    Vector<int> IndRow(nnz), IndCol(nnz);
    Vector<T, VectFull, Allocator> Val(nnz);

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();
    
    // First we convert the matrix in coordinate format and we permute the
    // indices.
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    long k = 0;
    for (int i = 0; i < m; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
	{
	  IndRow(k) = row_perm(i);
	  Val(k) = data[j];
	  IndCol(k) = row_perm(ind[j]);
	  if (IndCol(k) <= IndRow(k))
	    {
	      // We store only the superior part of the symmetric matrix.
	      int ind_tmp = IndRow(k);
	      IndRow(k) = IndCol(k);
	      IndCol(k) = ind_tmp;
	    }
	  k++;
	}
    
    // We sort with respect to row numbers.
    Sort(nnz, IndRow, IndCol, Val);
    
    // then column numbers
    A.Clear();
    Vector<long> Ptr(m+1);
    Ptr(0) = 0; k = 0;
    for (int i = 0; i < m; i++)
      {
	long first_index = k;
	// We get the size of the column i.
	while (k < nnz && IndRow(k) <= i)
	  k++;
        
        Sort(first_index, k-1, IndCol, Val);
        Ptr(i+1) = k;
      }
    
    A.SetData(m, n, Val, Ptr, IndCol);
  }
  
  
  //! Permutation of a symmetric matrix stored by columns.
  /*!
    B(row_perm(i), row_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, row_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColSymSparse, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm)
  {
    int m = A.GetM(), n = A.GetN();
    long nnz = A.GetDataSize();
    Vector<int> IndRow(nnz), IndCol(nnz);
    Vector<T, VectFull, Allocator> Val(nnz);

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();
    
    // First we convert the matrix in coordinate format and we permute the
    // indices.
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    long k = 0;
    for (int i = 0; i < m; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
	{
	  IndCol(k) = row_perm(i);
	  Val(k) = data[j];
	  IndRow(k) = row_perm(ind[j]);
	  if (IndCol(k) <= IndRow(k))
	    {
	      // We store only the superior part of the symmetric matrix.
	      int ind_tmp = IndRow(k);
	      IndRow(k) = IndCol(k);
	      IndCol(k) = ind_tmp;
	    }
	  k++;
	}
    
    // We sort with respect to column numbers.
    Sort(nnz, IndCol, IndRow, Val);
    
    // then row numbers
    A.Clear();
    Vector<long> Ptr(n+1);
    Ptr(0) = 0; k = 0;
    for (int i = 0; i < n; i++)
      {
	long first_index = k;
	// We get the size of the column i.
	while (k < nnz && IndCol(k) <= i)
	  k++;
        
        Sort(first_index, k-1, IndRow, Val);
        Ptr(i+1) = k;
      }
    
    A.SetData(m, n, Val, Ptr, IndRow);
  }

  
  //! Permutation of a general matrix stored by rows.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ArrayRowSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm)
  {
    int m = A.GetM(), n, i, i_, j, i2;
    IVect ind_tmp, iperm(m), rperm(m);
    for (i = 0; i < m; i++)
      {
	iperm(i) = i;
	rperm(i) = i;
      }
    // A(rperm(i),:) will be the place where is the initial row i.

    // Algorithm avoiding the allocation of another matrix.
    for (i = 0; i < m; i++)
      {
	// We get the index of row where the row initially placed on row i is.
	i2 = rperm(i);
	// We get the new index of this row.
	i_ = row_perm(i);

	// We fill ind_tmp of the permuted indices of columns of row i.
	n = A.GetRowSize(i2);
	ind_tmp.Reallocate(n);
	for (j = 0; j < n; j++)
	  ind_tmp(j) = col_perm(A.Index(i2,j));

	// We swap the two rows i and its destination row_perm(i).
	A.SwapRow(i2, i_);
	A.ReplaceIndexRow(i_, ind_tmp);

	// We update the indices iperm and rperm in order to keep in memory
	// the place where the row row_perm(i) is.
	int i_tmp = iperm(i_);
	iperm(i_) = iperm(i2);
	iperm(i2) = i_tmp;
	rperm(iperm(i_)) = i_;
	rperm(iperm(i2)) = i2;

	// We assemble the row i.
	A.AssembleRow(i_);
      }
  }


  //! Permutation of a general matrix stored by columns.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ArrayColSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm)
  {
    int n = A.GetN();
    IVect ind_tmp, iperm(n), rperm(n);
    for (int i = 0; i < n; i++)
      {
	iperm(i) = i;
	rperm(i) = i;
      }
    // A(:, rperm(i)) will be the place where is the initial column i.

    // Algorithm avoiding the allocation of another matrix.
    for (int i = 0; i < n; i++)
      {
	// We get the index of column where the column initially placed on column i is.
	int i2 = rperm(i);
	// We get the new index of this column.
	int i_ = col_perm(i);

	// We fill ind_tmp of the permuted indices of columns of row i.
	int p = A.GetColumnSize(i2);
	ind_tmp.Reallocate(p);
	for (int j = 0; j < p; j++)
	  ind_tmp(j) = row_perm(A.Index(i2, j));

	// We swap the two rows i and its destination col_perm(i).
	A.SwapColumn(i2, i_);
	A.ReplaceIndexColumn(i_, ind_tmp);

	// We update the indices iperm and rperm in order to keep in memory
	// the place where the column col_perm(i) is.
	int i_tmp = iperm(i_);
	iperm(i_) = iperm(i2);
	iperm(i2) = i_tmp;
	rperm(iperm(i_)) = i_;
	rperm(iperm(i2)) = i2;

	// We assemble the column i (to sort row numbers)
	A.AssembleColumn(i_);
      }
  }
  
  
  //! Permutation of a symmetric matrix stored by columns.
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void
  ApplyInversePermutation(Matrix<T, Prop, ArrayColSymSparse, Allocator>& A,
                          const IVect& row_perm, const IVect& col_perm)
  {
    // It is assumed that the permuted matrix is still symmetric! For example,
    // the user can provide row_perm = col_perm.
    int m = A.GetM();
    long nnz = A.GetDataSize();
    Vector<int> IndRow(nnz), IndCol(nnz);
    Vector<T, VectFull, Allocator> Val(nnz);

    // First we convert the matrix in coordinate format and we permute the
    // indices.
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    long k = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	{
	  IndCol(k) = col_perm(i);
	  Val(k) = A.Value(i,j);
	  IndRow(k) = row_perm(A.Index(i, j));
	  if (IndCol(k) <= IndRow(k))
	    {
	      // We store only the superior part of the symmetric matrix.
	      int ind_tmp = IndRow(k);
	      IndRow(k) = IndCol(k);
	      IndCol(k) = ind_tmp;
	    }
	  k++;
	}

    // We sort with respect to column numbers.
    Sort(nnz, IndCol, IndRow, Val);

    // A is filled.
    k = 0;
    for (int i = 0; i < m; i++)
      {
	long first_index = k;
	// We get the size of the column i.
	while (k < nnz && IndCol(k) <= i)
	  k++;
        
	int size_column = k - first_index;
	// If column not empty.
	if (size_column > 0)
	  {
	    A.ReallocateColumn(i, size_column);
	    k = first_index;
	    Sort(k, k+size_column-1, IndRow, Val);
	    for (int j = 0; j < size_column; j++)
	      {
		A.Index(i,j) = IndRow(k);
		A.Value(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearColumn(i);
      }
  }


  //! Permutation of a symmetric matrix stored by rows.
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void
  ApplyInversePermutation(Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
                          const IVect& row_perm, const IVect& col_perm)
  {
    // It is assumed that the permuted matrix is still symmetric! For example,
    // the user can provide row_perm = col_perm.
    int m = A.GetM();
    long nnz = A.GetDataSize();
    Vector<int> IndRow(nnz), IndCol(nnz);
    Vector<T,VectFull,Allocator> Val(nnz);

    // First we convert the matrix in coordinate format and we permute the
    // indices.
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    long k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    IndRow(k) = row_perm(i);
	    Val(k) = A.Value(i, j);
	    IndCol(k) = col_perm(A.Index(i, j));
	    if (IndCol(k) <= IndRow(k))
	      {
		// We store only the superior part of the symmetric matrix.
		int ind_tmp = IndRow(k);
		IndRow(k) = IndCol(k);
		IndCol(k) = ind_tmp;
	      }
	    k++;
	  }
      }
    
    // We sort with respect to row numbers.
    Sort(nnz, IndRow, IndCol, Val);

    // A is filled.
    k = 0;
    for (int i = 0; i < m; i++)
      {
	long first_index = k;
	// We get the size of the row i.
	while (k < nnz && IndRow(k) <= i)
	  k++;
	
        int size_row = k - first_index;
	// If row not empty.
	if (size_row > 0)
	  {
	    A.ReallocateRow(i, size_row);
	    k = first_index;
	    Sort(k, k+size_row-1, IndCol, Val);
	    for (int j = 0; j < size_row; j++)
	      {
		A.Index(i,j) = IndCol(k);
		A.Value(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearRow(i);
      }
  }

  
  // ApplyInversePermutation //  
  /////////////////////////////

  
  //////////////////////
  // ApplyPermutation //  
    
  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, col_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    Vector<int> inv_col_perm(col_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    for (int i = 0; i < col_perm.GetM(); i++)
      inv_col_perm(col_perm(i)) = i;
    
    ApplyInversePermutation(A, inv_row_perm, inv_col_perm);
  }


  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, col_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    Vector<int> inv_col_perm(col_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    for (int i = 0; i < col_perm.GetM(); i++)
      inv_col_perm(col_perm(i)) = i;
    
    ApplyInversePermutation(A, inv_row_perm, inv_col_perm);
  }
  
  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, row_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowSymSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    ApplyInversePermutation(A, inv_row_perm, inv_row_perm);
  }

  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, row_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColSymSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    ApplyInversePermutation(A, inv_row_perm, inv_row_perm);
  }
  
  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, col_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    Vector<int> inv_col_perm(col_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    for (int i = 0; i < col_perm.GetM(); i++)
      inv_col_perm(col_perm(i)) = i;
    
    ApplyInversePermutation(A, inv_row_perm, inv_col_perm);
  }
  

  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, col_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayColSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    Vector<int> inv_col_perm(col_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    for (int i = 0; i < col_perm.GetM(); i++)
      inv_col_perm(col_perm(i)) = i;
    
    ApplyInversePermutation(A, inv_row_perm, inv_col_perm);
  }
  
  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, row_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    ApplyInversePermutation(A, inv_row_perm, inv_row_perm);
  }
  
  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, row_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayColSymSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    ApplyInversePermutation(A, inv_row_perm, inv_row_perm);
  }
  
  
  // ApplyPermutation //  
  //////////////////////
  
  
  /////////////////
  // ScaleMatrix //    
  
  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, RowSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    T1* data = A.GetData();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (long j = ptr[i]; j < ptr[i+1]; j++ )
        data[j] *= scale_left(i) * scale_right(ind[j]);
  }
  
  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ColSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    T1* data = A.GetData();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    
    for (int i = 0; i < A.GetN(); i++ )
      for (long j = ptr[i]; j < ptr[i+1]; j++ )
        data[j] *= scale_left(ind[j]) * scale_right(i);
  }

  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_left).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, RowSymSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    T1* data = A.GetData();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (long j = ptr[i]; j < ptr[i+1]; j++ )
        data[j] *= scale_left(i) * scale_left(ind[j]);
  }
  
  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_left).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ColSymSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    T1* data = A.GetData();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    
    for (int i = 0; i < A.GetM(); i++ )
      for (long j = ptr[i]; j < ptr[i+1]; j++ )
        data[j] *= scale_left(i) * scale_left(ind[j]);
  }

  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetRowSize(i); j++ )
	A.Value(i,j) *= scale_left(i) * scale_right(A.Index(i, j));

  }

  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayColSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    for (int i = 0; i < A.GetN(); i++ )
      for (int j = 0; j < A.GetColumnSize(i); j++ )
	A.Value(i, j) *= scale_right(i) * scale_left(A.Index(i, j));

  }

  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowSymSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetRowSize(i); j++ )
	A.Value(i,j) *= scale_left(i) * scale_right(A.Index(i, j));

  }

  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayColSymSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetColumnSize(i); j++ )
	A.Value(i, j) *= scale_left(i) * scale_right(A.Index(i, j));

  }
  
  
  // ScaleMatrix //
  /////////////////
  
  
  /////////////////////
  // ScaleLeftMatrix //
  

  //! Each row is scaled.
  /*!
    We compute diag(S)*A where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, RowSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    T1* data = A.GetData();
    long* ptr = A.GetPtr();
    
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (long j = ptr[i]; j < ptr[i+1]; j++ )
        data[j] *= scale(i);
  }

  
  //! Each row is scaled.
  /*!
    We compute diag(S)*A where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ColSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    T1* data = A.GetData();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    
    for (int i = 0; i < A.GetN(); i++ )
      for (long j = ptr[i]; j < ptr[i+1]; j++ )
        data[j] *= scale(ind[j]);
  }
  
  
  //! Each row is scaled.
  /*!
    We compute diag(S)*A where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ArrayRowSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetRowSize(i); j++ )
	A.Value(i,j) *= scale(i);
  }

  
  //! Each row is scaled.
  /*!
    We compute diag(S)*A where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ArrayColSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    for (int i = 0; i < A.GetN(); i++ )
      for (int j = 0; j < A.GetColumnSize(i); j++ )
	A.Value(i, j) *= scale(A.Index(i, j));
  }

  
  // ScaleLeftMatrix //
  /////////////////////
  
  
  //////////////////////
  // ScaleRightMatrix //
  

  //! Each column is scaled.
  /*!
    We compute A*diag(S) where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, RowSparse, Allocator1>& A,
                        const Vector<T2, VectFull, Allocator2>& scale)
  {
    T1* data = A.GetData();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (long j = ptr[i]; j < ptr[i+1]; j++ )
        data[j] *= scale(ind[j]);
  }
  
  
  //! Each column is scaled.
  /*!
    We compute A*diag(S) where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, ColSparse, Allocator1>& A,
                        const Vector<T2, VectFull, Allocator2>& scale)
  {
    T1* data = A.GetData();
    long* ptr = A.GetPtr();
    
    for (int i = 0; i < A.GetN(); i++ )
      for (long j = ptr[i]; j < ptr[i+1]; j++ )
        data[j] *= scale(i);
  }

  
  //! Each column is scaled.
  /*!
    We compute A*diag(S) where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, ArrayRowSparse, Allocator1>& A,
                        const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetRowSize(i); j++ )
	A.Value(i,j) *= scale(A.Index(i, j));
  }
  
  
  //! Each column is scaled.
  /*!
    We compute A*diag(S) where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, ArrayColSparse, Allocator1>& A,
                        const Vector<T2, VectFull, Allocator2>& scale)
  {
    for (int i = 0; i < A.GetN(); i++ )
      for (int j = 0; j < A.GetColumnSize(i); j++ )
	A.Value(i, j) *= scale(i);
  }
  
  
  // ScaleRightMatrix //
  //////////////////////
  
  
} // end namespace

#define SELDON_FILE_PERMUTATION_SCALING_MATRIX_CXX
#endif
