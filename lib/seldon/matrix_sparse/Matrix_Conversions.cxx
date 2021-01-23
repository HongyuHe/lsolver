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


#ifndef SELDON_FILE_MATRIX_CONVERSIONS_CXX


#include "Matrix_Conversions.hxx"


namespace Seldon
{

  /*
    From CSR formats to "Matlab" coordinate format.
    index => starting index (usually 0 or 1)
    sym = true => the upper part and lower part are both generated
  */
  
  
  //! Conversion from RowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i; long j;
    int m = A.GetM();
    long nnz = A.GetDataSize();
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i < m; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndRow(j) = i + index;
	  IndCol(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }


  //! Conversion from ColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i; long j;
    int n = A.GetN();
    long nnz = A.GetDataSize();
    IndCol.Reallocate(nnz);
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i < n; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndCol(j) = i + index;
	  IndRow(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }


  //! Conversion from RowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i; long j;
    int m = A.GetM();
    long nnz = A.GetDataSize();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	// first we count the number of non-zero entries
 	// added by symmetry for each row
	Vector<long> NumNnzRow(m+1);
	NumNnzRow.Zero();
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    if (ind[j] != i)
	      NumNnzRow(ind[j] + 1)++;
	
	// number of non-zero entries in total
	nnz = 0;
	for (i = 0; i < m; i++)
	  nnz += NumNnzRow(i+1) + ptr[i+1] - ptr[i];
	
	// we construct cumulated index
	for (i = 0; i < m; i++)
	  NumNnzRow(i+1) += NumNnzRow(i);
	
	// arrays are filled already sorted (by rows)
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      // values located in upper part of the matrix
	      long num = NumNnzRow(i+1) + j;
	      IndRow(num) = i + index; 
	      IndCol(num) = ind[j] + index;
	      Val(num) = val[j];
	      
	      if (ind[j] != i)
		{
		  // values located in lower part of the matrix (by symmetry)
		  num = NumNnzRow(ind[j]) + ptr[ind[j]] + Ptr(ind[j]);
		  IndRow(num) = ind[j] + index;
		  IndCol(num) = i + index;
		  Val(num) = val[j];
		  Ptr(ind[j])++;
		}
	    }
      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(j) = i + index;
	      IndCol(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }


  //! Conversion from ColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i; long j;
    int m = A.GetM();
    long nnz = A.GetDataSize();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	// first we count the number of non-zero entries
	// added by symmetry for each column
	Vector<long> NumNnzCol(m+1);
	NumNnzCol.Zero();
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    if (ind[j] != i)
	      NumNnzCol(ind[j] + 1)++;
	
	// number of non-zero entries in total
	nnz = 0;
	for (i = 0; i < m; i++)
	  nnz += NumNnzCol(i+1) + ptr[i+1] - ptr[i];
	
	// we construct cumulated index
	for (i = 0; i < m; i++)
	  NumNnzCol(i+1) += NumNnzCol(i);
	
	// arrays are filled already sorted
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      // values located in the upper part of the matrix
	      long num = NumNnzCol(i) + j;
	      IndRow(num) = ind[j] + index;
	      IndCol(num) = i + index;
	      Val(num) = val[j];
	      
	      if (ind[j] != i)
		{
		  // values located in the lower part of the matrix
		  num = NumNnzCol(ind[j]) + ptr[ind[j]+1] + Ptr(ind[j]);
		  IndRow(num) = i + index;
		  IndCol(num) = ind[j] + index;
		  Val(num) = val[j];
		  Ptr(ind[j])++;
		}
	    }
      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndCol(j) = i + index;
	      IndRow(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }

  
  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from ArrayRowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    long nnz = A.GetDataSize();
    
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    long nb = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	{
	  IndRow(nb) = i + index;
	  IndCol(nb) = A.Index(i, j) + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int n = A.GetN();
    long nnz = A.GetDataSize();
    
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    long nb = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetColumnSize(i); j++)
	{
	  IndRow(nb) = A.Index(i, j) + index;
	  IndCol(nb) = i + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayRowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    long nnz = A.GetDataSize();
    if (sym)
      {
	// first we count the number of non-zero entries
	// added by symmetry for each row
	Vector<long> NumNnzRow(m+1);
	NumNnzRow.Zero();
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i, j) != i)
	      NumNnzRow(A.Index(i, j) + 1)++;
	
	// number of non-zero entries in total
	nnz = 0;
	for (i = 0; i < m; i++)
	  nnz += NumNnzRow(i+1) + A.GetRowSize(i);
	
	// we construct cumulated index
	Vector<long> NumUpper(m+1);
	NumUpper(0) = 0;
	for (i = 0; i < m; i++)
	  {
	    NumNnzRow(i+1) += NumNnzRow(i);
	    NumUpper(i+1) = NumUpper(i) + A.GetRowSize(i);
	  }
	
	// arrays are filled already sorted (by rows)
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      // values located in upper part of the matrix
	      long num = NumNnzRow(i+1) + NumUpper(i) + j;
	      int jcol = A.Index(i, j);
	      IndRow(num) = i + index; 
	      IndCol(num) = jcol + index;
	      Val(num) = A.Value(i, j);
	      
	      if (jcol != i)
		{
		  // values located in lower part of the matrix (by symmetry)
		  num = NumNnzRow(jcol) + NumUpper(jcol) + Ptr(jcol);
		  IndRow(num) = jcol + index;
		  IndCol(num) = i + index;
		  Val(num) = A.Value(i, j);
		  Ptr(jcol)++;
		}
	    }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	long nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }


  //! Conversion from ArrayColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    long nnz = A.GetDataSize();
    if (sym)
      {
	// first we count the number of non-zero entries
	// added by symmetry for each column
	Vector<long> NumNnzCol(m+1);
	NumNnzCol.Zero();
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    if (A.Index(i, j) != i)
	      NumNnzCol(A.Index(i, j) + 1)++;
	
	// number of non-zero entries in total
	nnz = 0;
	for (int i = 0; i < m; i++)
	  nnz += NumNnzCol(i+1) + A.GetColumnSize(i);
	
	// we construct cumulated index
	Vector<long> NumUpper(m+1);
	NumUpper(0) = 0;
	for (int i = 0; i < m; i++)
	  {
	    NumNnzCol(i+1) += NumNnzCol(i);
	    NumUpper(i+1) = NumUpper(i) + A.GetColumnSize(i);
	  }
	
	// arrays are filled already sorted
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      // values located in the upper part of the matrix
	      long num = NumNnzCol(i) + NumUpper(i) + j;
	      int jrow = A.Index(i, j);
	      IndRow(num) = jrow + index;
	      IndCol(num) = i + index;
	      Val(num) = A.Value(i, j);
	      
	      if (jrow != i)
		{
		  // values located in the lower part of the matrix
		  num = NumNnzCol(jrow) + NumUpper(jrow+1) + Ptr(jrow);
		  IndRow(num) = i + index;
		  IndCol(num) = jrow + index;
		  Val(num) = A.Value(i, j);
		  Ptr(jrow)++;
		}
	    }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	long nb = 0;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndRow(nb) = A.Index(i, j) + index;
	      IndCol(nb) = i + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }

  
  /*
    From "Matlab" coordinate format to CSR formats.
  */


  //! Conversion from coordinate format to RowSparse.
  /*! Contrary to the other conversion functions
    ConvertMatrix_from_Coordinates, this one accepts duplicates.
    \param[in] IndRow row indexes of the non-zero elements.
    \param[in] IndCol column indexes of the non-zero elements.
    \param[in] Val values of the non-zero elements.
    \param[inout] A matrix defined by \a IndRow, \a IndCol and \a Val.
    \param[in] index index of the first column and of the first row.
    In input, the matrix A can be allocated at the correct size,
    in that case the number of rows and columns are given by the size
    of the matrix A. If you want to take exactly the numbers of rows
    and columns contained in the triplet (IndRow, IndCol, Val)
    you should call A.Clear() before calling this function.
    \warning This function assumes that there are no duplicate values
    in the input arrays IndRow, IndCol
  */
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Values,
				 Matrix<T, Prop, RowSparse, Allocator3>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      return;

    long Nelement = IndRow.GetLength();    
    
    // detecting the size of the matrix
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we don't change the size given by the user
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // counting the number of elements for each row
    Vector<int> NumNnzRow(m+1);
    NumNnzRow.Zero();
    for (int i = 0; i < Nelement; i++)
      NumNnzRow(IndRow(i)-index)++;
    
    // the array Ptr can be constructed
    Vector<long> Ptr(m+1);
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzRow(i);
    
    // then the arrays Ind and Val are filled
    Vector<int> Ind(Nelement);
    Vector<T> Val(Nelement);
    NumNnzRow.Zero();
    for (long i = 0; i < Nelement; i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	long num = Ptr(irow) + NumNnzRow(irow);
	Ind(num) = icol;
	Val(num) = Values(i);
	NumNnzRow(irow)++;
      }

    // column numbers are sorted
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i+1)-1, Ind, Val);
    
    // A is set with arrays Val, Ptr, Ind
    A.SetData(m, n, Val, Ptr, Ind);
  }


#ifndef SWIG

  //! Conversion from coordinate format to ColSparse.
  /*!
    \warning This function assumes that there are no duplicate values
    in the input arrays IndRow, IndCol
  */
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Values,
				 Matrix<T, Prop, ColSparse, Allocator3>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      return;

    long Nelement = IndRow.GetLength();    
    
    // detecting the size of the matrix
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we don't change the size given by the user
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // counting the number of elements for each column
    Vector<int> NumNnzCol(n+1);
    NumNnzCol.Zero();
    for (int i = 0; i < Nelement; i++)
      NumNnzCol(IndCol(i)-index)++;
    
    // the array Ptr can be constructed
    Vector<long> Ptr(n+1);
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);

    // then the arrays Ind and Val are filled
    Vector<int> Ind(Nelement);
    Vector<T> Val(Nelement);
    NumNnzCol.Zero();
    for (long i = 0; i < Nelement; i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	int num = Ptr(icol) + NumNnzCol(icol);
	Ind(num) = irow;
	Val(num) = Values(i);
	NumNnzCol(icol)++;
      }

    // row numbers are sorted
    for (int i = 0; i < n; i++)
      Sort(Ptr(i), Ptr(i+1)-1, Ind, Val);
        
    // A is set with arrays Ptr, Ind and Val
    A.SetData(m, n, Val, Ptr, Ind);
  }


  //! Conversion from coordinate format to RowSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Values,
				 Matrix<T, Prop, RowSymSparse, Allocator3>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      return;

    long Nelement = IndRow.GetM();

    // detecting the size of the matrix
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we take the size of the input matrix
    m = max(m, n); n = m;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // only the upper part of the matrix is stored in A
    
    // counting the number of elements for each row
    Vector<int> NumNnzRow(m+1);
    NumNnzRow.Zero();
    for (int i = 0; i < Nelement; i++)
      if (IndRow(i) <= IndCol(i))
	NumNnzRow(IndRow(i)-index)++;
    
    // the array Ptr can be constructed
    Vector<long> Ptr(m+1);
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzRow(i);
    
    // then the arrays Ind and Val are filled
    long nnz = Ptr(m);
    Vector<int> Ind(nnz);
    Vector<T> Val(nnz);
    NumNnzRow.Zero();
    for (long i = 0; i < Nelement; i++)
      if (IndRow(i) <= IndCol(i))
	{
	  int irow = IndRow(i) - index;
	  int icol = IndCol(i) - index;
	  long num = Ptr(irow) + NumNnzRow(irow);
	  Ind(num) = icol;
	  Val(num) = Values(i);
	  NumNnzRow(irow)++;
	}
    
    // column numbers are sorted
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i+1)-1, Ind, Val);
    
    A.SetData(m, n, Val, Ptr, Ind);
  }


  //! Conversion from coordinate format to ColSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Values,
				 Matrix<T, Prop, ColSymSparse, Allocator3>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      return;

    long Nelement = IndRow.GetM();

    // detecting the size of the matrix
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we take the size of the input matrix
    m = max(m, n); n = m;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // only the upper part of the matrix is stored in A
    
    // counting the number of elements for each column
    Vector<int> NumNnzCol(n+1);
    NumNnzCol.Zero();
    for (int i = 0; i < Nelement; i++)
      if (IndRow(i) <= IndCol(i))
	NumNnzCol(IndCol(i)-index)++;
    
    // the array Ptr can be constructed
    Vector<long> Ptr(n+1);
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);

    // then the arrays Ind and Val are filled
    long nnz = Ptr(n);
    Vector<int> Ind(nnz);
    Vector<T> Val(nnz);
    NumNnzCol.Zero();
    for (long i = 0; i < Nelement; i++)
      if (IndRow(i) <= IndCol(i))
	{
	  int irow = IndRow(i) - index;
	  int icol = IndCol(i) - index;
	  long num = Ptr(icol) + NumNnzCol(icol);
	  Ind(num) = irow;
	  Val(num) = Values(i);
	  NumNnzCol(icol)++;
	}
    
    // row numbers are sorted
    for (int i = 0; i < n; i++)
      Sort(Ptr(i), Ptr(i+1)-1, Ind, Val);
    
    A.SetData(m, n, Val, Ptr, Ind);
  }

  
  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from coordinate format to ArrayRowSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Values,
				 Matrix<T, Prop, ArrayRowSparse,
				 Allocator3>& A, int index)
  {
    if (IndRow.GetM() <= 0)
      return;
    
    long Nelement = IndRow.GetLength();    
    
    // detecting the size of the matrix
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we don't change the size given by the user
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    A.Reallocate(m, n);
    
    // counting the number of elements for each row
    Vector<int> NumNnzRow(m+1);
    NumNnzRow.Zero();
    for (int i = 0; i < Nelement; i++)
      NumNnzRow(IndRow(i)-index)++;
    
    // rows are allocated
    for (int i = 0; i < m; i++)
      A.ReallocateRow(i, NumNnzRow(i));
    
    // then the matrix is filled
    NumNnzRow.Zero();
    for (long i = 0; i < Nelement; i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	int num = NumNnzRow(irow);
	A.Index(irow, num) = icol;
	A.Value(irow, num) = Values(i);
	NumNnzRow(irow)++;
      }

    // the matrix A is assembled to sort column numbers
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayColSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Values,
				 Matrix<T, Prop, ArrayColSparse,
				 Allocator3>& A, int index)
  {
    if (IndRow.GetM() <= 0)
      return;
    
    long Nelement = IndRow.GetLength();    
    
    // detecting the size of the matrix
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we don't change the size given by the user
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    A.Reallocate(m, n);
    
    // counting the number of elements for each column
    Vector<int> NumNnzCol(n+1);
    NumNnzCol.Zero();
    for (int i = 0; i < Nelement; i++)
      NumNnzCol(IndCol(i)-index)++;
    
    // columns are allocated
    for (int i = 0; i < n; i++)
      A.ReallocateColumn(i, NumNnzCol(i));
    
    // then the matrix is filled
    NumNnzCol.Zero();
    for (long i = 0; i < Nelement; i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	int num = NumNnzCol(icol);
	A.Index(icol, num) = irow;
	A.Value(icol, num) = Values(i);
	NumNnzCol(icol)++;
      }
    
    // the matrix A is assembled to sort row numbers
    A.Assemble();
  }
  
  
  //! Conversion from coordinate format to ArrayRowSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Values,
				 Matrix<T, Prop, ArrayRowSymSparse,
				 Allocator3>& A, int index)
  {
    if (IndRow.GetM() <= 0)
      return;
    
    long Nelement = IndRow.GetLength();    
    
    // detecting the size of the matrix
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we don't change the size given by the user
    m = max(m, n); n = m;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    A.Reallocate(m, n);
    
    // only the upper part of the matrix is stored in A
    
    // counting the number of elements for each row
    Vector<int> NumNnzRow(m+1);
    NumNnzRow.Zero();
    for (int i = 0; i < Nelement; i++)
      if (IndRow(i) <= IndCol(i))
	NumNnzRow(IndRow(i)-index)++;
    
    // rows are allocated
    for (int i = 0; i < m; i++)
      A.ReallocateRow(i, NumNnzRow(i));
    
    // then the matrix is filled
    NumNnzRow.Zero();
    for (long i = 0; i < Nelement; i++)
      if (IndRow(i) <= IndCol(i))
	{
	  int irow = IndRow(i) - index;
	  int icol = IndCol(i) - index;
	  int num = NumNnzRow(irow);
	  A.Index(irow, num) = icol;
	  A.Value(irow, num) = Values(i);
	  NumNnzRow(irow)++;
	}
    
    // the matrix A is assembled to sort column numbers
    A.Assemble();
  }
  
  
  //! Conversion from coordinate format to ArrayColSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Values,
				 Matrix<T, Prop, ArrayColSymSparse,
				 Allocator3>& A,
				 int index)
  {
    if (IndRow.GetM() <= 0)
      return;
    
    long Nelement = IndRow.GetLength();    
    
    // detecting the size of the matrix
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we don't change the size given by the user
    m = max(m, n); n = m;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    A.Reallocate(m, n);
    
    // counting the number of elements for each column
    Vector<int> NumNnzCol(n+1);
    NumNnzCol.Zero();
    for (int i = 0; i < Nelement; i++)
      if (IndRow(i) <= IndCol(i))
	NumNnzCol(IndCol(i)-index)++;
    
    // columns are allocated
    for (int i = 0; i < n; i++)
      A.ReallocateColumn(i, NumNnzCol(i));
    
    // then the matrix is filled
    NumNnzCol.Zero();
    for (long i = 0; i < Nelement; i++)
      if (IndRow(i) <= IndCol(i))
	{
	  int irow = IndRow(i) - index;
	  int icol = IndCol(i) - index;
	  int num = NumNnzCol(icol);
	  A.Index(icol, num) = irow;
	  A.Value(icol, num) = Values(i);
	  NumNnzCol(icol)++;
	}
    
    // the matrix A is assembled to sort row numbers
    A.Assemble();
  }
#endif


  /*
    From Sparse formats to CSC format
  */

  
  //! Conversion from RowSparse to CSC format
  /*!
    if sym_pat is equal to true, the pattern is symmetrized
    by adding artificial null entries
   */
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    long nnz = A.GetDataSize();
    int n = A.GetN();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    // we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(n);
    NumNnzCol.Zero();
    for (int i = 0; i < nnz; i++)
      NumNnzCol(ind_[i])++;
    
    // The array Ptr can be constructed
    Ptr.Reallocate(n+1);
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);

    // The arrays IndRow and Val are filled
    // (numbers of IndRow are already sorted by construction)
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzCol.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (long j = ptr_[i]; j < ptr_[i+1]; j++)
	{
	  long num = Ptr(ind_[j]) + NumNnzCol(ind_[j]);
	  IndRow(num) = i;
	  Val(num) = data_[j];
	  NumNnzCol(ind_[j])++;
	}
    
    long nb_new_val = 0;
    if (sym_pat)
      {
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        for (int i = 0; i < n; i++)
          {
	    long k = Ptr(i);
            for (long j = ptr_[i]; j < ptr_[i+1]; j++)
              {
                int irow = ind_[j];
                while (k < Ptr(i+1) && IndRow(k) < irow)
                  k++;
		
                if (k < Ptr(i+1) && IndRow(k) == irow)
                  // Already existing entry.
                  k++;
                else
                  {
                    // New entry.
                    NumNnzCol(i)++;
                    nb_new_val++;
                  }
              }
          }    
      }

    if (sym_pat && (nb_new_val > 0))
      {
        // Changing 'IndRow' and 'Val', and assembling the pattern.
        Vector<Tint1, VectFull, Alloc3> OldInd(IndRow);
        Vector<T, VectFull, Alloc4> OldVal(Val);
	IndRow.Reallocate(nnz + nb_new_val);
        Val.Reallocate(nnz + nb_new_val);
        long k = 0, nb = 0;
        T zero; SetComplexZero(zero);
        for (int i = 0; i < n; i++)
          {
	    // loop over non-zero entries of row i
	    for (long j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		int irow = ind_[j];
		while (k < Ptr(i+1) && OldInd(k) < irow)
		  {
		    IndRow(nb) = OldInd(k);
		    Val(nb) = OldVal(k);
		    nb++;
		    k++;
		  }
		
		if (k < Ptr(i+1) && OldInd(k) == irow)
		  {
		    // Already existing entry.
		    IndRow(nb) = OldInd(k);
		    Val(nb) = OldVal(k);
		    nb++;
		    k++;
		  }
		else
		  {
		    // New entry (null).
		    IndRow(nb) = irow;
		    Val(nb) = zero;
		    nb++;
		  }
	      }
	    
	    // last non-zero entries of column i
	    while (k < Ptr(i+1))
              {
                IndRow(nb) = OldInd(k);
                Val(nb) = OldVal(k);
		nb++;
                k++;
              }
          }
	
	// The array Ptr is updated
	Ptr(0) = 0;
	for (int i = 0; i < n; i++)
	  Ptr(i + 1) = Ptr(i) + NumNnzCol(i);	
      }
  }


  //! Conversion from ArrayRowSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    long nnz = A.GetDataSize();
    int n = A.GetN();

    // we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(n);
    NumNnzCol.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	NumNnzCol(A.Index(i, j))++;
    
    // The array Ptr can be constructed
    Ptr.Reallocate(n+1);
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);

    // The arrays IndRow and Val are filled
    // (numbers of IndRow are already sorted by construction)
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzCol.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	{
	  int jcol = A.Index(i, j);
	  long num = Ptr(jcol) + NumNnzCol(jcol);
	  IndRow(num) = i;
	  Val(num) = A.Value(i, j);
	  NumNnzCol(jcol)++;
	}

    long nb_new_val = 0;
    if (sym_pat)
      {
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        for (int i = 0; i < n; i++)
          {
	    long k = Ptr(i);
	    for (int j = 0; j < A.GetRowSize(i); j++)
              {
                int irow = A.Index(i, j);
                while (k < Ptr(i+1) && IndRow(k) < irow)
                  k++;
		
                if (k < Ptr(i+1) && IndRow(k) == irow)
                  // Already existing entry.
                  k++;
                else
                  {
                    // New entry.
                    NumNnzCol(i)++;
                    nb_new_val++;
                  }
              }
	  }
      }
    
    if (sym_pat && (nb_new_val > 0))
      {
        // Changing 'IndRow' and 'Val', and assembling the pattern.
        Vector<Tint1, VectFull, Alloc3> OldInd(IndRow);
        Vector<T, VectFull, Alloc4> OldVal(Val);
        IndRow.Reallocate(nnz + nb_new_val);
        Val.Reallocate(nnz + nb_new_val);
        long k = 0, nb = 0;
        T zero; SetComplexZero(zero);
        for (int i = 0; i < n; i++)
          {
	    for (int j = 0; j < A.GetRowSize(i); j++)
	      {
		int irow = A.Index(i, j);
		while (k < Ptr(i+1) && OldInd(k) < irow)
		  {
		    IndRow(nb) = OldInd(k);
		    Val(nb) = OldVal(k);
		    nb++;
		    k++;
		  }
		
		if (k < Ptr(i+1) && OldInd(k) == irow)
		  {
		    // Already existing entry.
		    IndRow(nb) = OldInd(k);
		    Val(nb) = OldVal(k);
		    nb++;
		    k++;
		  }
		else
		  {
		    // New entry (null).
		    IndRow(nb) = irow;
		    Val(nb) = zero;
		    nb++;
		  }
	      }
	    
	    while (k < Ptr(i+1))
              {
                IndRow(nb) = OldInd(k);
                Val(nb) = OldVal(k);
                nb++;
                k++;
              }
	  }
	
	// The array Ptr is updated
	Ptr(0) = 0;
	for (int i = 0; i < n; i++)
	  Ptr(i + 1) = Ptr(i) + NumNnzCol(i);	
      }
  }

  
  //! Conversion from ColSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with 'nnz' entries.
    long nnz = A.GetDataSize();
    int n = A.GetN();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    if (!sym_pat)
      {
	// direct conversion
	Ptr.Reallocate(n+1);
	IndRow.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (int i = 0; i <= n; i++)
	  Ptr(i) = ptr_[i];
	
	for (long i = 0; i < nnz; i++)
	  {
	    IndRow(i) = ind_[i];
	    Val(i) = data_[i];
	  }
      }
    else
      {
	// we count the number of non-zero entries for each row
	Vector<long> PtrRow(A.GetM()+1);
	PtrRow.Zero();
	for (int i = 0; i < nnz; i++)
	  PtrRow(ind_[i]+1)++;
    
	// replacing PtrRow with cumulated indexes
	for (int i = 0; i < A.GetM(); i++)
	  PtrRow(i+1) += PtrRow(i);
	
	// The array IndCol is filled (corresponding to CSR format)
	Vector<int> IndCol(nnz);
	for (int i = 0; i < n; i++)
	  for (long j = ptr_[i]; j < ptr_[i+1]; j++)
	    {
	      long num = PtrRow(ind_[j]);
	      IndCol(num) = i;
	      PtrRow(ind_[j])++;
	    }
	
	// reverting PtrRow
	for (int i = A.GetM()-1; i >= 1; i--)
	  PtrRow(i) -= PtrRow(i) - PtrRow(i-1);
	
	PtrRow(0) = 0;
	
	Ptr.Reallocate(n+1);
	for (int i = 0; i <= n; i++)
	  Ptr(i) = ptr_[i];

	long nb_new_val = 0;
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
	for (int i = 0; i < n; i++)
	  {
	    long k = ptr_[i];
	      
	    for (long j = PtrRow(i); j < PtrRow(i+1); j++)
	      {
		int irow = IndCol(j);
		while (k < ptr_[i+1] && ind_[k] < irow)
		  k++;
		
		if (k < ptr_[i+1] && ind_[k] == irow)
		  // Already existing entry.
		  k++;
		else
		  {
		    // New entry.
		    nb_new_val++;
		  }
	      }
	  }

	if (nb_new_val > 0)
	  {
	    // IndRow and Val are filled
	    IndRow.Reallocate(nnz + nb_new_val);
	    Val.Reallocate(nnz + nb_new_val);
	    long k = 0, nb = 0;
            T zero; SetComplexZero(zero);
	    Ptr.Zero();
	    for (int i = 0; i < n; i++)
	      {
		// loop over non-zero entries of row i
		for (long j = PtrRow(i); j < PtrRow(i+1); j++)
		  {
		    int irow = IndCol(j);
		    while (k < ptr_[i+1] && ind_[k] < irow)
		      {
			IndRow(nb) = ind_[k];
			Val(nb) = data_[k];
			nb++;
			k++;
		      }
		    
		    if (k < ptr_[i+1] && ind_[k] == irow)
		      {
			// Already existing entry.
			IndRow(nb) = ind_[k];
			Val(nb) = data_[k];
			nb++;
			k++;
		      }
		    else
		      {
			// New entry (null).
			IndRow(nb) = irow;
			Val(nb) = zero;
			nb++;
		      }
		  }
		
		// last non-zero entries of column i
		while (k < ptr_[i+1])
		  {
		    IndRow(nb) = ind_[k];
		    Val(nb) = data_[k];
		    nb++;
		    k++;
		  }

		Ptr(i+1) = nb;
	      }
	  }
      }
  }
  
  
  //! Conversion from ArrayColSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with 'nnz' entries.
    int n = A.GetN();
        
    if (!sym_pat)
      {
	// direct conversion
	Ptr.Reallocate(n+1);
	long nnz = 0;
	Ptr(0) = 0;
	for (int i = 0; i < n; i++)
	  {
	    nnz += A.GetColumnSize(i);
	    Ptr(i+1) = nnz;
	  }
	
	IndRow.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndRow(Ptr(i) + j) = A.Index(i, j);
	      Val(Ptr(i) + j) = A.Value(i, j);
	    }
      }
    else
      {
	// we count the number of non-zero entries for each row
	Vector<long> PtrRow(A.GetM()+1);
	PtrRow.Zero();
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    PtrRow(A.Index(i, j) + 1)++;
    
	// replacing PtrRow with cumulated indexes
	for (int i = 0; i < A.GetM(); i++)
	  PtrRow(i+1) += PtrRow(i);
	
	// The array IndCol is filled (corresponding to CSR format)
	long nnz = PtrRow(A.GetM());
	Vector<int> IndCol(nnz);
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      long num = PtrRow(A.Index(i, j));
	      IndCol(num) = i;
	      PtrRow(A.Index(i, j))++;
	    }
	
	// reverting PtrRow
	for (long i = A.GetM()-1; i >= 1; i--)
	  PtrRow(i) -= PtrRow(i) - PtrRow(i-1);
	
	PtrRow(0) = 0;
	
	Ptr.Reallocate(n+1);
	Ptr(0) = 0;
	for (int i = 0; i < n; i++)
	  Ptr(i+1) = Ptr(i) + A.GetColumnSize(i);
	
	long nb_new_val = 0;
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
	for (int i = 0; i < n; i++)
	  {
	    int k = 0;
	    for (long j = PtrRow(i); j < PtrRow(i+1); j++)
	      {
		int irow = IndCol(j);
		while (k < A.GetColumnSize(i) && A.Index(i, k) < irow)
		  k++;
		
		if (k < A.GetColumnSize(i) && A.Index(i, k) == irow)
		  // Already existing entry.
		  k++;
		else
		  {
		    // New entry.
		    nb_new_val++;
		  }
	      }
	  }

	if (nb_new_val > 0)
	  {
	    // IndRow and Val are filled
	    IndRow.Reallocate(nnz + nb_new_val);
	    Val.Reallocate(nnz + nb_new_val);
	    long nb = 0;
            T zero; SetComplexZero(zero);
	    Ptr.Zero();
	    for (int i = 0; i < n; i++)
	      {
		int k = 0;
		// loop over non-zero entries of row i
		for (long j = PtrRow(i); j < PtrRow(i+1); j++)
		  {
		    int irow = IndCol(j);
		    while (k < A.GetColumnSize(i) && A.Index(i, k) < irow)
		      {
			IndRow(nb) = A.Index(i, k);
			Val(nb) = A.Value(i, k);
			nb++;
			k++;
		      }
		    
		    if (k < A.GetColumnSize(i) && A.Index(i, k) == irow)
		      {
			// Already existing entry.
			IndRow(nb) = A.Index(i, k);
			Val(nb) = A.Value(i, k);
			nb++;
			k++;
		      }
		    else
		      {
			// New entry (null).
			IndRow(nb) = irow;
			Val(nb) = zero;
			nb++;
		      }
		  }
		
		// last non-zero entries of column i
		while (k < A.GetColumnSize(i))
		  {
		    IndRow(nb) = A.Index(i, k);
		    Val(nb) = A.Value(i, k);
		    nb++;
		    k++;
		  }

		Ptr(i+1) = nb;
	      }
	  }
      }    
  }
  
  
  //! Conversion from ColSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();    
    long nnz = A.GetDataSize();
    
    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz);
    Value.Reallocate(nnz);
    
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    for (int i = 0; i <= n; i++)
      Ptr(i) = ptr_[i];
    
    for (long i = 0; i < nnz; i++)
      {
	Ind(i) = ind_[i];
	Value(i) = data_[i];
      }
  }

  
  //! Conversion from ColSymSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    
    Vector<Tint1, VectFull, Alloc3> IndCol;
    
    // this function returns IndRow, IndCol with sorted column numbers
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    long nnz = 0;
    for (long i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }
    
    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from ArrayColSymSparse to symmetric CSC format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    long nnz = A.GetDataSize();
    
    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz);
    Value.Reallocate(nnz);
    
    Ptr(0) = 0;
    for (int i = 1; i <= n; i++)
      Ptr(i) = Ptr(i-1) + A.GetColumnSize(i-1);
    
    long nb = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	{
	  Ind(nb) = A.Index(i, j);
	  Value(nb) = A.Value(i, j);
	  nb++;
	}
  }

  
  //! Conversion from ArrayColSymSparse to CSC format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    
    Vector<Tint1, VectFull, Alloc3> IndCol;

    // this function returns IndRow, IndCol with sorted column numbers    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    long nnz = 0;
    for (long i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }
    
    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from RowSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    long nnz = A.GetDataSize();
    int n = A.GetN();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    // we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(n);
    NumNnzCol.Zero();
    for (long i = 0; i < nnz; i++)
      NumNnzCol(ind_[i])++;
    
    // The array Ptr can be constructed
    Ptr.Reallocate(n+1);
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);
    
    // The arrays IndRow and Val are filled
    // (numbers of IndRow are already sorted by construction)
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzCol.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (long j = ptr_[i]; j < ptr_[i+1]; j++)
	{
	  long num = Ptr(ind_[j]) + NumNnzCol(ind_[j]);
	  IndRow(num) = i;
	  Val(num) = data_[j];
	  NumNnzCol(ind_[j])++;
	}
  }

  
  //! Conversion from RowSymSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    
    Vector<Tint1, VectFull, Alloc2> IndCol;
    
    // this function sorts by row numbers
    // by symmetry, row numbers = col numbers, we invert the two arguments
    // to have sorted column numbers
    ConvertMatrix_to_Coordinates(A, IndCol, IndRow, Value, 0, true);
    
    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    long nnz = 0;
    for (long i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }
    
    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);
  }
  
  
  //! Conversion from ArrayRowSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    long nnz = A.GetDataSize();
    int n = A.GetN();
    
    // we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(n);
    NumNnzCol.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	NumNnzCol(A.Index(i, j))++;
    
    // The array Ptr can be constructed
    Ptr.Reallocate(n+1);
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);
    
    // The arrays IndRow and Val are filled
    // (numbers of IndRow are already sorted by construction)
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzCol.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	{
	  int icol = A.Index(i, j);
	  long num = Ptr(icol) + NumNnzCol(icol);
	  IndRow(num) = i;
	  Val(num) = A.Value(i, j);
	  NumNnzCol(icol)++;
	}
  }
  
  
  //! Conversion from ArrayRowSymSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetM();

    Vector<Tint1, VectFull, Alloc2> IndCol;
    
    // this function sorts by row numbers
    // by symmetry, row numbers = col numbers, we invert the two arguments
    // to have sorted column numbers
    ConvertMatrix_to_Coordinates(A, IndCol, IndRow, Value, 0, true);
    
    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    long nnz = 0;
    for (long i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }
    
    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowMajor, Allocator>& A,
		  Matrix<T, Prop, RowMajor, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSymPacked, Allocator>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator>& B)
  {
    B = A;
  }

  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColMajor, Allocator>& A,
		  Matrix<T, Prop, ColMajor, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSymPacked, Allocator>& A,
		  Matrix<T, Prop, ColSymPacked, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Allocator>& A,
		  Matrix<T, Prop, RowSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Allocator>& A,
		  Matrix<T, Prop, RowSymSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSparse, Allocator>& A,
		  Matrix<T, Prop, ColSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSymSparse, Allocator>& A,
		  Matrix<T, Prop, ColSymSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayRowSymSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayRowSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayColSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayColSymSparse, Allocator>& B)
  {
    B = A;
  }

  
  //! Conversion from ArrayColSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csc)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int> IndRow;
    Vector<long> IndCol;

    General sym;
    ConvertToCSC(mat_array, sym, IndCol, IndRow, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();

    mat_csc.SetData(m, n, Val, IndCol, IndRow);
  }

  
  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSparse, Alloc2>& B)
  {
    Vector<long> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! Conversion from RowSparse to ArrayColSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Alloc1>& A,
		  Matrix<T, Prop, ArrayColSparse, Alloc2>& B)
  {
    Vector<long> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
	B.ReallocateColumn(i, size_col);
	for (long j = Ptr(i); j < Ptr(i+1); j++)
	  {
	    B.Index(i, j-Ptr(i)) = Ind(j);
	    B.Value(i, j-Ptr(i)) = Val(j);
	  }
      }
  }

  
  //! Conversion from ArrayRowSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Val;

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    General sym;
    ConvertToCSC(mat_array, sym, Ptr, IndRow, Val);

    mat_csr.SetData(m, n, Val, Ptr, IndRow);
  }
  
  
  //! Conversion from RowSymSparse to ColSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B)
  {
    Vector<long> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from ArrayRowSymSparse to ColSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B)
  {
    Vector<long> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from ArrayColSymSparse to ColSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B)
  {
    Vector<long> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from RowSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, RowSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<long> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from ArrayRowSymSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T1, VectFull, Allocator1> AllVal;

    int n = A.GetM();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, AllVal);

    B.SetData(n, n, AllVal, Ptr, Ind);
  }


  //! Conversion from ColSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<long> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from ArrayColSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<long> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  /*
    From Sparse formats to CSR format
  */
  
  
  //! Conversion from RowSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    long nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);
    for (int i = 0; i <= m; i++)
      Ptr(i) = ptr_[i];
    
    for (long i = 0; i < nnz; i++)
      {
        IndCol(i) = ind_[i];
        Value(i) = data_[i];
      }
  }
  
  
  //! Conversion from ColSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Computation of the indexes of the beginning of rows.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    // Counting the number of entries per row.
    for (long i = 0; i < nnz; i++)
      Ptr(ind_[i])++;

    // Incrementing in order to get the row indexes.
    long increment = 0, size; int num_row;
    for (int i = 0; i < m; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }
    
    // Last index.
    Ptr(m) = increment;
    
    // 'Offset' will be used to get current positions of new entries.
    Vector<Tint0, VectFull, Alloc2> Offset(Ptr);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);

    // Loop over the columns.
    for (int j = 0; j < n; j++)
      for (long i = ptr_[j]; i < ptr_[j + 1]; i++)
	{
	  num_row = ind_[i];
	  IndCol(Offset(num_row)) = j;
	  Value(Offset(num_row)) = data_[i];
	  Offset(num_row)++;
	}
  }
  
  
  //! Conversion from ArrayColSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    int n = A.GetN();
    long  nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    // Computation of the indexes of the beginning of rows.
    Ptr.Reallocate(m + 1);
    Ptr.Zero();
    // Counting the number of entries per row.
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	Ptr(A.Index(i, j))++;
    
    // Incrementing in order to get the row indexes.
    long increment = 0, size; int num_row;
    for (int i = 0; i < m; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }
    
    // Last index.
    Ptr(m) = increment;
    
    // 'Offset' will be used to get current positions of new entries.
    Vector<Tint0, VectFull, Alloc2> Offset(Ptr);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);

    // Loop over the columns.
    for (int j = 0; j < n; j++)
      for (int i = 0; i < A.GetColumnSize(j); i++)
	{
	  num_row = A.Index(j, i);
	  IndCol(Offset(num_row)) = j;
	  Value(Offset(num_row)) = A.Value(j, i);
	  Offset(num_row)++;
	}
  }
  
  
  //! Conversion from ArrayRowSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    long nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0; Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    IndCol(nnz + j) = A.Index(i, j);
	    Value(nnz + j) = A.Value(i, j);
	  }
	
	nnz += A.GetRowSize(i);
	Ptr(i+1) = nnz;
      }
  }


  //! Conversion from ArrayRowSparse to symmetric CSR 
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    // only upper part of A is stored
    long  nnz = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) >= i)
	  nnz++;
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0; Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (A.Index(i, j) >= i)
	    {
	      IndCol(nnz) = A.Index(i, j);
	      Value(nnz) = A.Value(i, j);
	      nnz++;
	    }
	
	Ptr(i+1) = nnz;
      }
  }

  
  //! Conversion from ColSymSparse to symmetric CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    long nnz = A.GetDataSize();
    int m = A.GetM();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    // we count the number of non-zero entries for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (long i = 0; i < nnz; i++)
      NumNnzRow(ind_[i])++;
    
    // The array Ptr can be constructed
    Ptr.Reallocate(m+1);
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzRow(i);
    
    // The arrays IndCol and Val are filled
    // (numbers of IndCol are already sorted by construction)
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzRow.Zero();
    for (int i = 0; i < m; i++)
      for (long j = ptr_[i]; j < ptr_[i+1]; j++)
	{
	  long num = Ptr(ind_[j]) + NumNnzRow(ind_[j]);
	  IndCol(num) = i;
	  Val(num) = data_[j];
	  NumNnzRow(ind_[j])++;
	}
  }
  
  
  //! Conversion from ColSymSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint1, VectFull, Alloc3> IndRow;
   
    // this function sorts by column numbers
    // by symmetry, row numbers = col numbers, we invert the two arguments
    // to have sorted column numbers 
    ConvertMatrix_to_Coordinates(A, IndCol, IndRow, Value, 0, true);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from ArrayColSymSparse to symmetric CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    long nnz = A.GetDataSize();
    int m = A.GetM();
    
    // we count the number of non-zero entries for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	NumNnzRow(A.Index(i, j))++;
    
    // The array Ptr can be constructed
    Ptr.Reallocate(m+1);
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzRow(i);
    
    // The arrays IndRow and Val are filled
    // (numbers of IndRow are already sorted by construction)
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzRow.Zero();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	{
	  int icol = A.Index(i, j);
	  long num = Ptr(icol) + NumNnzRow(icol);
	  IndCol(num) = i;
	  Val(num) = A.Value(i, j);
	  NumNnzRow(icol)++;
	}
  }
  
  
  //! Conversion from ArrayColSymSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint1, VectFull, Alloc3> IndRow;
    
    // this function sorts by column numbers
    // by symmetry, row numbers = col numbers, we invert the two arguments
    // to have sorted row numbers
    ConvertMatrix_to_Coordinates(A, IndCol, IndRow, Value, 0, true);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }
  
  
  //! Conversion from RowSymSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& IndRow,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Number of rows and non-zero entries.
    long nnz = A.GetDataSize();
    int m = A.GetM();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Allocation of arrays for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    long ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (long k = ptr_[i]; k < ptr_[i+1]; k++)
	  {
	    IndCol(ind) = ind_[k];
	    Val(ind) = data_[k];
	    ind++;
	  }
	
	IndRow(i + 1) = ind;
      }
  }

  
  //! Conversion from RowSymSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint1, VectFull, Alloc3> IndRow;

    // this function returns IndRow, IndCol with sorted row numbers    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from ArrayRowSymSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& IndRow,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Number of rows and non-zero entries.
    long nnz = A.GetDataSize();
    int m = A.GetM();

    // Allocation of arrays for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    long ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < A.GetRowSize(i); k++)
	  {
	    IndCol(ind) = A.Index(i, k);
	    Val(ind) = A.Value(i, k);
	    ind++;
	  }
	IndRow(i + 1) = ind;
      }
  }

  
  //! Conversion from ArrayRowSymSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint1, VectFull, Alloc3> IndRow;
    
    // this function returns IndRow, IndCol with sorted row numbers    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
    
  }

  
  //! Conversion from column-oriented sparse to row-oriented sparse.
  /*!
    \param[in] A matrix to be converted.
    \param[out] B converted matrix.
  */
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T1, Prop1, ColSparse, Alloc1>& A,
		  Matrix<T2, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T1, VectFull, Alloc2> Value;
    
    General sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }

  
  //! Conversion from ArrayColSparse to RowSparse
  /*!
    \param[in] A matrix to be converted.
    \param[out] B converted matrix.
  */
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T1, Prop1, ArrayColSparse, Alloc1>& A,
		  Matrix<T2, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T1, VectFull, Alloc2> Value;
    
    General sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }
  
  
  //! Conversion from ArrayRowSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<long> IndRow;
    Vector<int> IndCol;

    General sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Val, IndRow, IndCol);
  }


  //! Conversion from ArrayRowSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<long> IndRow;
    Vector<int> IndCol;
    
    Symmetric sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Val, IndRow, IndCol);
  }

  
  //! Conversion from ColSymSparse to RowSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, RowSymSparse, Alloc2>& B)
  {
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Value;
    
    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from ColSymSparse to RowSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Value;
    
    General unsym;
    ConvertToCSR(A, unsym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from ArrayColSymSparse to RowSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, RowSymSparse, Alloc2>& B)
  {
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Value;
    
    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from ArrayColSymSparse to RowSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Value;
    
    General unsym;
    ConvertToCSR(A, unsym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from RowSymSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<long> IndRow;
    Vector<int> IndCol;

    General unsym;
    ConvertToCSR(mat_array, unsym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }
  
  
  //! Conversion from ArrayRowSymSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<long> IndRow;
    Vector<int> IndCol;

    Symmetric sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }


  //! Conversion from ArrayRowSymSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<long> IndRow;
    Vector<int> IndCol;

    General unsym;
    ConvertToCSR(mat_array, unsym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }


  /******************************
   * From Sparse to ArraySparse *
   ******************************/
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
    
  }

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);
    
    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = Ptr(i+1) - Ptr(i);
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = Ind(Ptr(i) + j);
	    B.Value(i, j) = Value(Ptr(i) + j);
	  }
      }    
  }
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Value);
    
    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
	B.ReallocateColumn(i, size_col);
	for (int j = 0; j < size_col; j++)
	  {
	    B.Index(i, j) = Ind(Ptr(i) + j);
	    B.Value(i, j) = Value(Ptr(i) + j);
	  }
      }    
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int i; long j;
    int m = A.GetM();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T0* val = A.GetData();

    // first we count the number of non-zero entries
    // added by symmetry for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (i = 0; i < m; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  NumNnzRow(ind[j])++;
    
    // Allocation of B
    B.Reallocate(m, m);
    for (i = 0; i < m; i++)
      B.ReallocateRow(i, ptr[i+1] - ptr[i] + NumNnzRow(i));
    
    // B is filled with sorted column numbers
    Vector<int> Ptr(m);
    Ptr.Zero();
    for (i = 0; i < m; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	{
	  // values located in upper part of the matrix
	  int num = NumNnzRow(i) + j - ptr[i];
	  B.Index(i, num) = ind[j];
	  B.Value(i, num) = val[j];
	  
	  if (ind[j] != i)
	    {
	      // values located in lower part of the matrix (by symmetry)
	      num = Ptr(ind[j]);
	      B.Index(ind[j], num) = i;
	      B.Value(ind[j], num) = val[j];
	      Ptr(ind[j])++;
	    }
	}
  }

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
    
  }

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
	int size_col = ptr_[i+1] - ptr_[i];
	B.ReallocateColumn(i, size_col);
	for (int j = 0; j < size_col; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }    
  }

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& Acsc,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {    
    int m = Acsc.GetM();
    int n = Acsc.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    // conversion to RowSparse
    Matrix<T0, Prop0, RowSparse, Allocator0> A;
    CopyMatrix(Acsc, A);
    
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
  }
  
  
  /***********************************
   * From ArraySparse to ArraySparse *
   ***********************************/
  
  
  //! From ArrayRowSymSparse to ArrayRowSymSparse (T0 and T1 different)
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRowSize(i);
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = A.Index(i, j);
            B.Value(i, j) = A.Value(i, j);
          }
      }
  }
  
  
  //! From ArrayRowSparse to ArrayRowSparse (T0 and T1 different)
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRowSize(i);
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = A.Index(i, j);
            B.Value(i, j) = A.Value(i, j);
          }
      }
  }
  
  
  //! upper part of A is used to obtain a symmetric matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    B.Reallocate(A.GetM(), A.GetN());
    for (int i = 0; i < A.GetM(); i++)
      {
	int k = 0;
	while ( (k < A.GetRowSize(i)) && (A.Index(i, k) < i))
	  k++;
	
	if (k < A.GetRowSize(i))
	  {
	    int size_row = A.GetRowSize(i) - k;
	    B.ReallocateRow(i, size_row);
	    for (int j = k; j < A.GetRowSize(i); j++)
	      {
		B.Index(i, j-k) = A.Index(i, j);
		B.Value(i, j-k) = A.Value(i, j);
	      }
	  }
      }
  }

  
  //! conversion from ArrayColSymSparse to ArrayRowSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);
    
    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = Ptr(i+1) - Ptr(i);
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = Ind(Ptr(i) + j);
	    B.Value(i, j) = Value(Ptr(i) + j);
	  }
      }    
  }
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& Acsc,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {    
    int m = Acsc.GetM();
    int n = Acsc.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    // conversion to RowSparse
    Matrix<T0, Prop0, RowSparse, Allocator0> A;
    CopyMatrix(Acsc, A);
    
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
  }
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int i, j;
    int m = A.GetM();

    // first we count the number of non-zero entries
    // added by symmetry for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  NumNnzRow(A.Index(i, j))++;
    
    // Allocation of B
    B.Reallocate(m, m);
    for (i = 0; i < m; i++)
      B.ReallocateRow(i, A.GetRowSize(i) + NumNnzRow(i));
    
    // B is filled with sorted column numbers
    Vector<int> Ptr(m);
    Ptr.Zero();
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	{
	  // values located in upper part of the matrix
	  int num = NumNnzRow(i) + j;
	  int jcol = A.Index(i, j);
	  B.Index(i, num) = jcol;
	  B.Value(i, num) = A.Value(i, j);
	  
	  if (jcol != i)
	    {
	      // values located in lower part of the matrix (by symmetry)
	      num = Ptr(jcol);
	      B.Index(jcol, num) = i;
	      B.Value(jcol, num) = A.Value(i, j);
	      Ptr(jcol)++;
	    }
	}
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i, j;
    int m = A.GetM();
    
    // first we count the number of non-zero entries
    // added by symmetry for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  NumNnzRow(A.Index(i, j))++;
    
    // Allocation of B
    B.Reallocate(m, m);
    for (i = 0; i < m; i++)
      B.ReallocateColumn(i, A.GetRowSize(i) + NumNnzRow(i));
    
    // B is filled with sorted row numbers
    Vector<int> Ptr(m);
    Ptr.Zero();
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	{
	  // values located in upper part of the matrix
	  int num = NumNnzRow(i) + j;
	  int jcol = A.Index(i, j);
	  B.Index(i, num) = jcol;
	  B.Value(i, num) = A.Value(i, j);
	  
	  if (jcol != i)
	    {
	      // values located in lower part of the matrix (by symmetry)
	      num = Ptr(jcol);
	      B.Index(jcol, num) = i;
	      B.Value(jcol, num) = A.Value(i, j);
	      Ptr(jcol)++;
	    }
	}
  }

  
  //! Conversion from ArrayRowSparse to ArrayColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    // Matrix (m,n) with nnz entries.
    int n = A.GetN();

    // we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(n);
    NumNnzCol.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	NumNnzCol(A.Index(i, j))++;
    
    // allocating matrix B
    B.Reallocate(A.GetM(), n);
    for (int i = 0; i < n; i++)
      B.ReallocateColumn(i, NumNnzCol(i));
    
    // The matrix B can be filled
    // (numbers of IndRow are already sorted by construction)
    NumNnzCol.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	{
	  int jcol = A.Index(i, j);
	  int num = NumNnzCol(jcol);
	  B.Index(jcol, num) = i;
	  B.Value(jcol, num) = A.Value(i, j);
	  NumNnzCol(jcol)++;
	}
  }
  
  
  /***********************
   * GetSymmetricPattern *
   ***********************/


  //! Returns pattern of A + A' in CSR format
  /*!
    From a sparse matrix, we compute the pattern of A + A'
    so that this pattern is symmetric even if A is non-symmetric    
   */
  template<class T, class Prop, class Storage, class Allocator,
           class Tint0, class Tint1, class Allocator2, class Allocator3>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Vector<Tint0, VectFull, Allocator2>& Ptr,
                           Vector<Tint1, VectFull, Allocator3>& Ind)
  {
    typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type T0;
    int n = A.GetM();

    // Converting to coordinates.
    Vector<int> IndRow, IndCol;
    Vector<T0> Value;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value);

    // clearing values
    Value.Clear();

    Vector<int> IndRow2(IndRow), IndCol2(IndCol);

    // we count the number of non-zero entries for each row
    // and each column
    Vector<int> NumNnzRow(n), NumNnzCol(n);
    NumNnzRow.Zero();
    NumNnzCol.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	NumNnzRow(IndRow(i))++;
	NumNnzCol(IndCol(i))++;
      }

    Vector<long> PtrRow(n+1), PtrCol(n+1);
    PtrRow(0) = 0; PtrCol(0) = 0;
    for (int i = 0; i < n; i++)
      {
	PtrRow(i+1) = PtrRow(i) + NumNnzRow(i);
	PtrCol(i+1) = PtrCol(i) + NumNnzCol(i);
      }

    // row numbers are sorted
    NumNnzRow.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow2(i);
	long num = PtrRow(irow) + NumNnzRow(irow);
	IndRow(num) = irow;
	IndCol(num) = IndCol2(i);
	NumNnzRow(irow)++;
      }

    // column numbers are sorted
    NumNnzCol.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int icol = IndCol(i);
	long num = PtrCol(icol) + NumNnzCol(icol);
	IndRow2(num) = IndRow(i);
	IndCol2(num) = icol;
	NumNnzCol(icol)++;
      }

    // intermediary arrays are cleared
    NumNnzRow.Clear(); NumNnzCol.Clear();
    PtrRow.Clear(); PtrCol.Clear();

    Tint0 max_nnz = 0;
    for (long i = 0; i < IndRow.GetM(); i++)
      if (IndRow(i) <= IndCol(i))
        max_nnz++;

    for (long i = 0; i < IndRow.GetM(); i++)
      if (IndCol2(i) <= IndRow2(i))
        max_nnz++;

    // then symmetrization of pattern and conversion to csr.
    Vector<int> Index(2*n);
    Ptr.Reallocate(n+1);
    Ind.Reallocate(max_nnz);
    Tint0 j_end = 0;
    int size_row = 0;
    Tint0 j2_end = 0;
    Ptr(0) = 0;
    for (int i = 0; i < A.GetM(); i++)
      {
        size_row = 0;
        // We retrieve column numbers.
        while ( (j_end < IndRow.GetM()) && (IndRow(j_end) == i))
          {
            if (IndRow(j_end) <= IndCol(j_end))
              {
                Index(size_row) = IndCol(j_end);
                size_row++;
              }

            j_end++;
          }

        while ( (j2_end < IndCol2.GetM()) && (IndCol2(j2_end) == i))
          {
            if (IndCol2(j2_end) <= IndRow2(j2_end))
              {
                Index(size_row) = IndRow2(j2_end);
                size_row++;
              }

            j2_end++;
          }

        // Sorting indexes.
        Assemble(size_row, Index);

        // Updating Ptr, Ind.
        for (int j = 0; j < size_row; j++)
	  Ind(Ptr(i) + j) = Index(j);

        Ptr(i+1) = Ptr(i) + size_row;
      }

    IndRow2.Clear(); IndCol2.Clear();
    IndRow.Clear(); IndCol.Clear();
    Ind.Resize(Ptr(n));
  }


  template<class T, class Prop, class Storage, class Allocator, class AllocI>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Matrix<int, Symmetric, RowSymSparse, AllocI>& B)
  {
    Vector<long> Ptr; Vector<int> Ind;

    GetSymmetricPattern(A, Ptr, Ind);

    int n = A.GetM();
    Vector<int, VectFull, AllocI> Val(Ptr(n));
    // We put Ptr and Ind into the matrix B.
    B.SetData(n, n, Val, Ptr, Ind);
  }

  
  /*****************************************************
   * Conversion from sparse matrices to dense matrices *
   *****************************************************/
  
  
  //! conversion from RowSparse to RowMajor
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Allocator1>& A,
		  Matrix<T, Prop, RowMajor, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();
    
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
        B(i, ind[j]) = data[j];
    
  }

  
  //! conversion from ArrayRowSparse to RowMajor
  template<class T1, class T2, class Prop1, 
	   class Prop2, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		  Matrix<T2, Prop2, RowMajor, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    
    B.Reallocate(m, n);
    T2 zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        B(i, A.Index(i, j)) = A.Value(i, j);
    
  }


  //! conversion from RowSymSparse to RowSymPacked
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();
    
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
        B(i, ind[j]) = data[j];
    
  }

  
  //! conversion from ArrayRowSymSparse to RowSymPacked
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B(i, A.Index(i, j)) = A.Value(i, j);
    
  }


  //! conversion from ArrayRowSymSparse to ColSymPacked
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, ColSymPacked, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B(i, A.Index(i, j)) = A.Value(i, j);
    
  }


  //! conversion from ArrayRowSymSparse to RowSym
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSym, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B.Val(i, A.Index(i, j)) = A.Value(i, j);
    
  }


  //! conversion from ArrayRowSymSparse to ColSym
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, ColSym, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B.Val(i, A.Index(i, j)) = A.Value(i, j);
    
  }
  

  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, RowMajor, Alloc2>& B)
  {
    int m, n;
    MatGetLocalSize(A.GetPetscMatrix(), &m, &n);
    n = A.GetN();
    B.Reallocate(m, n);
    T *local_a;
    MatGetArray(A.GetPetscMatrix(), &local_a);
    for (int i = 0; i < m; i++)
      for(int j = 0; j < n; j++)
        B(i, j) = local_a[i + j * m];
    MatRestoreArray(A.GetPetscMatrix(), &local_a);
  }


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, RowMajor, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIDense, Alloc2>& B)
  {
    T *local_data;
    MatGetArray(B.GetPetscMatrix(), &local_data);
    int mlocal, nlocal;
    MatGetLocalSize(B.GetPetscMatrix(), &mlocal, &nlocal);
    Matrix<T, Prop1, ColMajor, Alloc1> local_D;
    local_D.SetData(mlocal, B.GetN(), local_data);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        local_D(i, j) = A(i, j);
    local_D.Nullify();
    MatRestoreArray(B.GetPetscMatrix(), &local_data);
  }

  
  /*****************************************************
   * Conversion from dense matrices to sparse matrices *
   *****************************************************/
  
  
  //! conversion from RowSymPacked to RowSymSparse
  template<class T>
  void ConvertToSparse(const Matrix<T, Symmetric, RowSymPacked>& A,
                       Matrix<T, Symmetric, RowSymSparse>& B,
		       const T& threshold)
  {
    long nnz = 0;
    int n = A.GetM();
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
        if (abs(A(i, j)) > threshold)
          nnz++;
    
    Vector<int> IndCol(nnz); Vector<long> IndRow(n+1); 
    Vector<T> Value(nnz);
    nnz = 0; IndRow(0) = 0;
    for (int i = 0; i < n; i++)
      {
        for (int j = i; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              IndCol(nnz) = j;
              Value(nnz) = A(i, j);
              nnz++;
            }
        
        IndRow(i+1) = nnz;
      }
    
    B.SetData(n, n, Value, IndRow, IndCol);
    
  }
  
  
  //! conversion from RowMajor to ArrayRowSparse
  template<class T>
  void ConvertToSparse(const Matrix<T, General, RowMajor>& A,
                       Matrix<T, General, ArrayRowSparse>& B,
		       const T& threshold)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = 0;
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            size_row++;
        
        B.ReallocateRow(i, size_row);
        
        size_row = 0;
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              B.Index(i, size_row) = j;
              B.Value(i, size_row) = A(i, j);
              size_row++;
            }
      }
  }
  
  
  //! conversion from RowMajor to RowSparse
  template<class T>
  void ConvertToSparse(const Matrix<T, General, RowMajor>& A,
                       Matrix<T, General, RowSparse>& B,
		       const T& threshold)
  {
    long nnz = 0;
    int m = A.GetM();
    int n = A.GetN();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        if (abs(A(i, j)) > threshold)
          nnz++;
    
    Vector<int> IndCol(nnz); Vector<long> IndRow(m+1); 
    Vector<T> Value(nnz);
    nnz = 0; IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              IndCol(nnz) = j;
              Value(nnz) = A(i, j);
              nnz++;
            }
        
        IndRow(i+1) = nnz;
      }
    
    B.SetData(m, n, Value, IndRow, IndCol);
    
  }
    
} // namespace Seldon.

#define SELDON_FILE_MATRIX_CONVERSIONS_CXX
#endif
