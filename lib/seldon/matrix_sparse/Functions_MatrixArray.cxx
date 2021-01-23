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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_ARRAY_CXX

/*
  Functions defined in this file:
  (storage ArrayRowSparse, ArrayColSparse, etc)

  X = A(i, :)
  GetRow(A, i, X)

  X = A(:, j)
  GetCol(A, j, X)

  A(i, :) = X
  SetRow(X, i, A)

  A(:, j) = X
  SetCol(X, j, A)

  alpha.M*X + beta.Y -> Y
  MltAdd(alpha, M, X, beta, Y)

  alpha.A + B -> B
  Add(alpha, A, B)

  alpha.M -> M
  Mlt(alpha, M)
  
  Highest absolute value of A.
  MaxAbs(A)

  1-norm of matrix A.
  Norm1(A)

  infinity norm of matrix A.
  NormInf(A)

  transpose of matrix A
  Transpose(A)

  B = transpose(A)
  Transpose(A, B)
  
  conjugate of transpose of matrix A
  TransposeConj(A)

  conjugate of matrix A
  Conjugate(A)

  alpha.A*B + beta.C -> C
  MltAdd(alpha, A, B, beta, C)

*/

namespace Seldon
{

  
  ////////////////////
  // GetRow, SetRow //
  
  
  //! Extracts a row from a sparse matrix
  /*!
    \param A sparse matrix
    \param i row index
    \param X row extracted
    X = A(i, :)
  */
  template<class T0, class Allocator0,
	   class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, ArrayRowSparse, Allocator0>& A,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = A.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif
    
    int size_row = A.GetRowSize(i);
    X.Reallocate(size_row);
    for (int j = 0; j < size_row; j++)
      {
	X.Index(j) = A.Index(i, j);
	X.Value(j) = A.Value(i, j);
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
  void GetRow(const Matrix<T0, General, ArrayColSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    list<pair<int, T0> > vec;
    for (int j = 0; j < M.GetN(); j++)
      for (int k = 0; k < M.GetColumnSize(j); k++)
	if (M.Index(j, k) == i)
	  vec.push_back(make_pair(j, M.Value(j, k)));
    
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
  void GetRow(const Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    list<pair<int, T0> > vec;
    // beginning of row i
    for (int j = 0; j < i; j++)
      for (int k = 0; k < M.GetRowSize(j); k++)
	if (M.Index(j, k) == i)
	  vec.push_back(make_pair(j, M.Value(j, k)));
    
    // end of row i
    for (int k = 0; k < M.GetRowSize(i); k++)
      vec.push_back(make_pair(M.Index(i, k), M.Value(i, k)));
    
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
  void GetRow(const Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    list<pair<int, T0> > vec;
    // beginning of row i
    for (int k = 0; k < M.GetColumnSize(i); k++)
      vec.push_back(make_pair(M.Index(i, k), M.Value(i, k)));
    
    // end of row i
    for (int j = i+1; j < M.GetN(); j++)
      for (int k = 0; k < M.GetColumnSize(j); k++)
	if (M.Index(j, k) == i)
	  vec.push_back(make_pair(j, M.Value(j, k)));
    
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
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, ArrayRowSparse, Allocator0>& M)
  {
    M.ClearRow(i);
    M.ReallocateRow(i, X.GetM());    
    for (int k = 0; k < X.GetM(); k++)
      {
	M.Index(i, k) = X.Index(k);
	M.Value(i, k) = X.Value(k);
      }
  }
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, ArrayColSparse, Allocator0>& M)
  {
    T0 val;
    int p = 0;
    for (int j = 0; j < M.GetN(); j++)
      {
	while ( (p < X.GetM()) && (X.Index(p) < j))
	  p++;
	
	bool present_X = false;
	if ( (p < X.GetM()) && (X.Index(p) == j))	  
	  {
	    present_X = true;
	    val = X.Value(p);
	  }
	
	int size_col = M.GetColumnSize(j);
	bool present_val = false;
	int k = 0;
	while ( (k < size_col) && (M.Index(j, k) < i))
	  k++;
	
	if ( (k < size_col) && (M.Index(j, k) == i))
	  {
	    if (!present_X)
	      {
		// reducing size of column
		if (size_col > 1)
		  {
		    if (k == size_col-1)
		      M.ResizeColumn(j, size_col-1);
		    else
		      {
			int last_col = M.Index(j, size_col-1);
			val = M.Value(j, size_col-1);
			M.ResizeColumn(j, size_col-1);
			for (int q = k; q < size_col-2; q++)
			  {
			    M.Index(j, q) = M.Index(j, q+1);
			    M.Value(j, q) = M.Value(j, q+1);
			  }
			
			M.Index(j, size_col-2) = last_col;
			M.Value(j, size_col-2) = val;
		      }
		  }
		else
		  M.ClearColumn(j);
	      }
	    else
	      M.Value(j, k) = val;
	    
	    present_val = true;
	  }
	
	if (!present_val && present_X)
	  {
	    // increasing size of column
	    M.ResizeColumn(j, size_col+1);
	    for (int q = size_col; q > k; q--)
	      {
		M.Index(j, q) = M.Index(j, q-1);
		M.Value(j, q) = M.Value(j, q-1);
	      }
	    
	    M.Index(j, k) = i;
	    M.Value(j, k) = val;
	  }
      }
  }

  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M)
  {
    T0 val;
    SetComplexZero(val);
    int p = 0;
    for (int j = 0; j < i; j++)
      {
	while ( (p < X.GetM()) && (X.Index(p) < j))
	  p++;
	
	bool present_X = false;
	if ( (p < X.GetM()) && (X.Index(p) == j))	  
	  {
	    present_X = true;
	    val = X.Value(p);
	    p++;
	  }
	
	int size_row = M.GetRowSize(j);
	bool present_val = false;
	int k = 0;
	while ( (k < size_row) && (M.Index(j, k) < i))
	  k++;
	
	if ( (k < size_row) && (M.Index(j, k) == i))
	  {
	    if (!present_X)
	      {
		// reducing size of row
		if (size_row > 1)
		  {
		    if (k == size_row-1)
		      M.ResizeRow(j, size_row-1);
		    else
		      {
			int last_col = M.Index(j, size_row-1);
			val = M.Value(j, size_row-1);
			M.ResizeRow(j, size_row-1);
			for (int q = k; q < size_row-2; q++)
			  {
			    M.Index(j, q) = M.Index(j, q+1);
			    M.Value(j, q) = M.Value(j, q+1);
			  }
			
			M.Index(j, size_row-2) = last_col;
			M.Value(j, size_row-2) = val;
		      }
		  }
		else
		  M.ClearRow(j);
	      }
	    else
	      M.Value(j, k) = val;
	    
	    present_val = true;
	  }
	
	if (!present_val && present_X)
	  {
	    // increasing size of row
	    M.ResizeRow(j, size_row+1);
	    for (int q = size_row; q > k; q--)
	      {
		M.Index(j, q) = M.Index(j, q-1);
		M.Value(j, q) = M.Value(j, q-1);
	      }
	    
	    M.Index(j, k) = i;
	    M.Value(j, k) = val;
	  }
      }
    
    // and changing row i
    M.ClearRow(i);
    if (p < X.GetM())
      {
	M.ReallocateRow(i, X.GetM() - p);
	int k = 0;
	while (p < X.GetM())
	  {
	    M.Index(i, k) = X.Index(p);
	    M.Value(i, k) = X.Value(p);
	    k++;
	    p++;
	  }
      }
  }
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M)
  {
    T1 val;
    int p = 0;
    // setting column i
    M.ClearColumn(i);
    while ( (p < X.GetM()) && (X.Index(p) <= i))
      p++;
    
    if (p > 0)
      {
	M.ReallocateColumn(i, p);
	for (int k = 0; k < p; k++)
	  {
	    M.Index(i, k) = X.Index(k);
	    M.Value(i, k) = X.Value(k);
	  }	
      }
    
    // then modifying last columns
    for (int j = i+1; j < M.GetN(); j++)
      {
	while ( (p < X.GetM()) && (X.Index(p) < j))
	  p++;
	
	bool present_X = false;
	if ( (p < X.GetM()) && (X.Index(p) == j))	  
	  {
	    present_X = true;
	    val = X.Value(p);
	  }
	
	int size_col = M.GetColumnSize(j);
	bool present_val = false;
	int k = 0;
	while ( (k < size_col) && (M.Index(j, k) < i))
	  k++;
	
	if ( (k < size_col) && (M.Index(j, k) == i))
	  {
	    if (!present_X)
	      {
		// reducing size of column
		if (size_col > 1)
		  {
		    if (k == size_col-1)
		      M.ResizeColumn(j, size_col-1);
		    else
		      {
			int last_col = M.Index(j, size_col-1);
			val = M.Value(j, size_col-1);
			M.ResizeColumn(j, size_col-1);
			for (int q = k; q < size_col-2; q++)
			  {
			    M.Index(j, q) = M.Index(j, q+1);
			    M.Value(j, q) = M.Value(j, q+1);
			  }
			
			M.Index(j, size_col-2) = last_col;
			M.Value(j, size_col-2) = val;
		      }
		  }
		else
		  M.ClearColumn(j);
	      }
	    else
	      M.Value(j, k) = val;
	    
	    present_val = true;
	  }
	
	if (!present_val && present_X)
	  {
	    // increasing size of column
	    M.ResizeColumn(j, size_col+1);
	    for (int q = size_col; q > k; q--)
	      {
		M.Index(j, q) = M.Index(j, q-1);
		M.Value(j, q) = M.Value(j, q-1);
	      }
	    
	    M.Index(j, k) = i;
	    M.Value(j, k) = val;
	  }

      }
  }
  
  
  // GetRow, SetRow //
  ////////////////////


  ////////////////////
  // GetCol, SetCol //

  
  //! Extracts a column from a matrix
  /*!
    \param M matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, ArrayRowSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int n = M.GetN();
    if (j < 0 || j >= n)
      throw WrongIndex("GetCol()",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif

    int m = M.GetM();

    list<pair<int, T0> > vec;
    for (int i = 0; i < m; i++)
      for (int k = 0; k < M.GetRowSize(i); k++)
	if (M.Index(i, k) == j)
	  vec.push_back(make_pair(i, M.Value(i, k)));
    
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
  
  
  //! Extracts a column from a matrix
  /*!
    \param M matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template<class T0, class Allocator0,
	   class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, ArrayColSparse, Allocator0>& A,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int n = A.GetN();
    if (j < 0 || j >= n)
      throw WrongIndex("GetCol()",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif
    
    int size_col = A.GetColumnSize(j);
    X.Reallocate(size_col);
    for (int k = 0; k < size_col; k++)
      {
	X.Index(k) = A.Index(j, k);
	X.Value(k) = A.Value(j, k);
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
  void GetCol(const Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M,
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
  void GetCol(const Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
    // symmetric matrix row = col
    GetRow(M, j, X);
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, ArrayRowSparse, Allocator0>& M)
  {
    T0 val;
    SetComplexZero(val);
    int p = 0;
    for (int i = 0; i < M.GetM(); i++)
      {
	while ( (p < X.GetM()) && (X.Index(p) < i))
	  p++;
	
	bool present_X = false;
	if ( (p < X.GetM()) && (X.Index(p) == i))	  
	  {
	    present_X = true;
	    val = X.Value(p);
	  }
	
	int size_row = M.GetRowSize(i);
	bool present_val = false;
	int k = 0;
	while ( (k < size_row) && (M.Index(i, k) < j))
	  k++;
	

	if ( (k < size_row) && (M.Index(i, k) == j))
	  {
	    if (!present_X)
	      {
		// reducing size of row
		if (size_row > 1)
		  {
		    if (k == size_row-1)
		      M.ResizeRow(j, size_row-1);
		    else
		      {
			int last_row = M.Index(i, size_row-1);
			val = M.Value(i, size_row-1);
			M.ResizeRow(i, size_row-1);
			for (int q = k; q < size_row-2; q++)
			  {
			    M.Index(i, q) = M.Index(i, q+1);
			    M.Value(i, q) = M.Value(i, q+1);
			  }
			
			M.Index(i, size_row-2) = last_row;
			M.Value(i, size_row-2) = val;
		      }
		  }
		else
		  M.ClearRow(i);
	      }
	    else
	      M.Value(i, k) = val;
	    
	    present_val = true;
	  }
	
	if (!present_val && present_X)
	  {
	    // increasing size of row
	    M.ResizeRow(i, size_row+1);
	    for (int q = size_row; q > k; q--)
	      {
		M.Index(i, q) = M.Index(i, q-1);
		M.Value(i, q) = M.Value(i, q-1);
	      }
	    
	    M.Index(i, k) = j;
	    M.Value(i, k) = val;
	  }	
      }
  }
  
  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, ArrayColSparse, Allocator0>& M)
  {
    M.ClearColumn(j);
    M.ReallocateColumn(j, X.GetM());    
    for (int k = 0; k < X.GetM(); k++)
      {
	M.Index(j, k) = X.Index(k);
	M.Value(j, k) = X.Value(k);
      }
  }
  
  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M)
  {
    // symmetric matrix, row = column
    SetRow(X, j, M);
  }
   
  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M)
  {
    // symmetric matrix, row = column
    SetRow(X, j, M);
  }
  
  
  // GetCol, SetCol //
  ////////////////////


  /////////
  // Mlt //

  /*** ArrayRowSymSparse ***/


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, Symmetric,
		 ArrayRowSymSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    T1 val;
    int m = A.GetM(), n, p;
    C.Fill(0);
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.Index(i, k);
	    val = A.Value(i, k);
	    
	    if (p == i)
	      C(i) += val * B(i);
	    else
	      {
		C(i) += val * B(p);
		C(p) += val * B(i);
	      }
	  }
      }
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		 ArrayRowSymSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    if (!Trans.ConjTrans())
      {
	MltVector(A, B, C);
	return;
      }

    C.Fill(0);
    
    int m = A.GetM(), n, p;
    T1 val;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.Index(i, k);
	    val = conjugate(A.Value(i, k));
	    
	    if (p == i)
	      C(i) += val * B(i);
	    else
	      {
		C(i) += val * B(p);
		C(p) += val * B(i);
	      }
	  }
      }
  }

  
  /*** ArrayColSymSparse ***/


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, Symmetric,
		 ArrayColSymSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    C.Fill(0);
    int m = A.GetM(), n, p;
    T1 val;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetColumnSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.Index(i, k);
	    val = A.Value(i, k);
	    
	    if (p == i)
	      C(i) += val * B(i);
	    else
	      {
		C(i) += val * B(p);
		C(p) += val * B(i);
	      }
	  }
      }
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		 ArrayColSymSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    if (!Trans.ConjTrans())
      {
	MltVector(A, B, C);
	return;
      }

    C.Fill(0);
    
    int m = A.GetM(), n, p;
    T3 val;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetColumnSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.Index(i, k);
	    val = conjugate(A.Value(i, k));
	    
	    if (p == i)
	      C(i) += val * B(i);
	    else
	      {
		C(i) += val * B(p);
		C(p) += val * B(i);
	      }
	  }
      }
  }
  
  
  /*** ArrayRowSparse ***/


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    T3 zero, temp; SetComplexZero(zero);
    int m = A.GetM(), n;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	temp = zero;
	for (int k = 0; k < n ; k++)
	  temp += A.Value(i, k)*B(A.Index(i, k));
	
	C(i) = temp;
      }
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    if (Trans.NoTrans())
      {
	MltVector(A, B, C);
	return;
      }

    C.Fill(0);        
    int m = A.GetM(), n;

    if (Trans.Trans())
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      C(A.Index(i, k)) += A.Value(i, k)*B(i);
	  }
      }
    else
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      C(A.Index(i, k)) += conjugate(A.Value(i, k))*B(i);
	  }
      }
  }
  
  
  /*** ArrayColSparse ***/
  
  
  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    C.Fill(0);
    for (int i = 0 ; i < A.GetN(); i++)
      for (int k = 0; k < A.GetColumnSize(i); k++)
	C(A.Index(i, k)) += A.Value(i, k) * B(i);
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    if (Trans.NoTrans())
      {
	MltVector(A, B, C);
	return;
      }

    T3 zero, temp;
    SetComplexZero(zero);
    
    if (Trans.Trans())
      {
	for (int i = 0 ; i < A.GetN(); i++)
	  {
	    temp = zero;
	    for (int k = 0; k < A.GetColumnSize(i); k++)
	      temp += A.Value(i, k) * B(A.Index(i, k));
	    
	    C(i) = temp;
	  }
      }
    else
      {
	for (int i = 0 ; i < A.GetN(); i++)
	  {
	    temp = zero;
	    for (int k = 0; k < A.GetColumnSize(i); k++)
	      temp += conjugate(A.Value(i, k)) * B(A.Index(i, k));
	    
	    C(i) = temp;
	  }
      }
  }
  

  // Mlt //
  /////////
  
  
  ////////////
  // MltAdd //


  /*** ArrayRowSymSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric,
		    ArrayRowSymSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T1 val;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);

                if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);

		if (p==i)
		  C(i) += alpha * val * B(i);
		else
		  {
		    C(i) += alpha * val * B(p);
		    C(p) += alpha * val * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		    ArrayRowSymSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    if (!Trans.ConjTrans())
      {
	MltAddVector(alpha, A, B, beta, C);
	return;
      }

    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T1 val;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = conjugate(A.Value(i, k));

		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = conjugate(A.Value(i, k));

		if (p==i)
		  C(i) += alpha * val * B(i);
		else
		  {
		    C(i) += alpha * val * B(p);
		    C(p) += alpha * val * B(i);
		  }
	      }
	  }
      }
  }

  
  /*** ArrayColSymSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric,
		    ArrayColSymSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T1 val;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);

		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);

		if (p==i)
		  C(i) += alpha * val * B(i);
		else
		  {
		    C(i) += alpha * val * B(p);
		    C(p) += alpha * val * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		    ArrayColSymSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    if (!Trans.ConjTrans())
      {
	MltAddVector(alpha, A, B, beta, C);
	return;
      }

    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T3 val;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = conjugate(A.Value(i, k));

		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = alpha * conjugate(A.Value(i, k));

		if (p==i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
  }
  
  
  /*** ArrayRowSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero; SetComplexZero(zero);
    T0 one; SetComplexOne(one);
    
    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);
    
    int m = A.GetM(), n, p;
    T1 val;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(i) += val * B(p);
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(i) += alpha * val * B(p);
              }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    if (Trans.NoTrans())
      {
	MltAddVector(alpha, A, B, beta, C);
	return;
      }

    T4 zero; SetComplexZero(zero);
    T0 one; SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);
    
    int m = A.GetM(), n, p;
    T1 val;

    if (Trans.Trans())
      {
	if (alpha == one)
	  {
	    for (int i = 0 ; i < m ; i++)
	      {
		n = A.GetRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.Index(i, k);
		    val = A.Value(i, k);
		    C(p) += val * B(i);
		  }
	      }
	  }
	else // alpha != 1.
	  {
	    for (int i = 0 ; i < m ; i++)
	      {
		n = A.GetRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.Index(i, k);
		    val = A.Value(i, k);
		    C(p) += alpha * val * B(i);
		  }
	      }
	  }
      }
    else
      {
	if (alpha == one)
	  {
	    for (int i = 0 ; i < m ; i++)
	      {
		n = A.GetRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.Index(i, k);
		    val = conjugate(A.Value(i, k));
		    C(p) += val * B(i);
		  }
	      }
	  }
	else // alpha != 1.
	  {
	    for (int i = 0 ; i < m ; i++)
	      {
		n = A.GetRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.Index(i, k);
		    val = conjugate(A.Value(i, k));
		    C(p) += alpha * val * B(i);
		  }
	      }
	  }
      }
  }
  
  
  /*** ArrayColSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    if (alpha == one)
      {
	for (int i = 0 ; i < A.GetN(); i++)
          for (int k = 0; k < A.GetColumnSize(i); k++)
            C(A.Index(i, k)) += A.Value(i, k) * B(i);
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < A.GetN(); i++)
          for (int k = 0; k < A.GetColumnSize(i); k++)
            C(A.Index(i, k)) += alpha * A.Value(i, k) * B(i);
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    if (Trans.NoTrans())
      {
	MltAddVector(alpha, A, B, beta, C);
	return;
      }

    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    if (Trans.Trans())
      {
	if (alpha == one)
	  {
	    for (int i = 0 ; i < A.GetN(); i++)
	      for (int k = 0; k < A.GetColumnSize(i); k++)
		C(i) += A.Value(i, k) * B(A.Index(i, k));
	  }
	else // alpha != 1.
	  {
	    for (int i = 0 ; i < A.GetN(); i++)
	      for (int k = 0; k < A.GetColumnSize(i); k++)
		C(i) += alpha * A.Value(i, k) * B(A.Index(i, k));
	  }
      }
    else
      {
	if (alpha == one)
	  {
	    for (int i = 0 ; i < A.GetN(); i++)
	      for (int k = 0; k < A.GetColumnSize(i); k++)
		C(i) += conjugate(A.Value(i, k)) * B(A.Index(i, k));
	  }
	else // alpha != 1.
	  {
	    for (int i = 0 ; i < A.GetN(); i++)
	      for (int k = 0; k < A.GetColumnSize(i); k++)
		C(i) += alpha * conjugate(A.Value(i, k)) * B(A.Index(i, k));
	  }
      }
  }

  
  // MltAdd //
  ////////////



  /////////
  // Add //


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		 Matrix<T2, General, ArrayRowSparse, Allocator2>& B)
  {
    int m = B.GetM(), n;
    Vector<T2, VectFull, Allocator2> value(B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	for (int j = 0; j < n; j++)
	  value(j) = alpha*A.Value(i, j);
        
	B.AddInteractionRow(i, n, A.GetIndex(i), value.GetData());
      }
  }

  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		 Matrix<T2, General, ArrayColSparse, Allocator2>& B)
  {
    int m = B.GetN(), n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetColumnSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.Value(i, j);
        
	Mlt(alpha, value);
	B.AddInteractionColumn(i, n, A.GetIndex(i), value.GetData());
      }
  }

  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
		 Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B)
  {
    int m = B.GetM(),n;
    Vector<T2, VectFull, Allocator2> value(B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	for (int j = 0; j < n; j++)
	  value(j) = alpha*A.Value(i, j);

	B.AddInteractionRow(i, n, A.GetIndex(i), value.GetData());
      }
  }
  
  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, ArrayColSymSparse, Allocator1>& A,
		 Matrix<T2, Symmetric, ArrayColSymSparse, Allocator2>& B)
  {
    int m = B.GetN(), n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetColumnSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.Value(i, j);
        
	Mlt(alpha, value);
	B.AddInteractionColumn(i, n, A.GetIndex(i), value.GetData());
      }
  }


  // C = C + alpha complex(A,B)
  template<class T0, class T1, class T2, class T3, class Allocator1,
	   class Allocator2, class Allocator3>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		 const Matrix<T2, General, ArrayRowSparse, Allocator2>& B,
		 Matrix<complex<T3>, General, ArrayRowSparse, Allocator3>& C)
  {
    int m = B.GetM(),n1,n2,size_row;;
    Vector<complex<T3>, VectFull, Allocator3> val_row;
    IVect ind_row;
    for (int i = 0 ; i < m ; i++)
      {
	n1 = A.GetRowSize(i);
	n2 = B.GetRowSize(i);
	size_row = n1 + n2;
	val_row.Reallocate(size_row);
	ind_row.Reallocate(size_row);
	for (int j = 0 ; j < n1 ; j++)
	  {
	    ind_row(j) = A.Index(i, j);
	    val_row(j) = alpha*complex<T3>(A.Value(i, j), 0);
	  }

	for (int j = 0 ; j < n2 ; j++)
	  {
	    ind_row(j+n1) = B.Index(i, j);
	    val_row(j+n1) = alpha * complex<T3>(0, B.Value(i, j));
	  }

	C.AddInteractionRow(i, size_row, ind_row, val_row);
      }
  }


  // C = C + alpha complex(A,B)
  template<class T0, class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
		 const Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B,
		 Matrix<complex<T3>, Symmetric, ArrayRowSymSparse, Allocator3>& C)
  {
    int m = B.GetM(), n1, n2, size_row;
    Vector<complex<T3>, VectFull, Allocator3> val_row;
    IVect ind_row;
    for (int i = 0 ; i < m ; i++)
      {
	n1 = A.GetRowSize(i);
	n2 = B.GetRowSize(i);
	size_row = n1 + n2;
	val_row.Reallocate(size_row);
	ind_row.Reallocate(size_row);
	for (int j = 0 ; j < n1 ; j++)
	  {
	    ind_row(j) = A.Index(i, j);
	    val_row(j) = alpha * complex<T3>(A.Value(i, j), 0);
	  }

	for (int j = 0 ; j < n2 ; j++)
	  {
	    ind_row(j+n1) = B.Index(i, j);
	    val_row(j+n1) = alpha * complex<T3>(0, B.Value(i, j));
	  }

	C.AddInteractionRow(i, size_row, ind_row, val_row);
      }
  }


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
		 Matrix<T2, General, ArrayRowSparse, Allocator2>& B)
  {
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0; i < B.GetM(); i++)
      {
	int n = A.GetRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  {
            value(j) = A.Value(i, j);
            if (A.Index(i, j) != i)
              B.AddInteraction(A.Index(i, j), i, alpha*value(j));
          }
        
	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetIndex(i), value.GetData());
      }
  }
  
  
  // Add //
  /////////



  /////////
  // Mlt //


  template<class T0, class T, class Allocator>
  void MltScalar(const T0& alpha, Matrix<T, General, ArrayRowSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void MltScalar(const T0& alpha, Matrix<T, General, ArrayColSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void MltScalar(const T0& alpha,
		 Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void MltScalar(const T0& alpha,
		 Matrix<T, Symmetric, ArrayColSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  // Matrix-matrix product (sparse matrix against dense matrix)
  template<class T0, class T1, class Prop1, class Allocator1, class T4,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T3, class Prop3, class Storage3, class Allocator3>
  void MltAddMatrix(const T0& alpha,
		    const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		    const Matrix<T2, Prop2, Storage2, Allocator2>& B,
		    const T4& beta,
		    Matrix<T3, Prop3, Storage3, Allocator3>& C)
  {
    if (Storage2::Sparse || Storage3::Sparse)
      throw WrongArgument("Mlt", "Function intended for product "
                          " between a sparse matrix and a dense matrix");

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, "MltAdd(alpha, A, B, beta, C)");
#endif
    
    T4 zero; SetComplexZero(zero);
    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    int m = A.GetM();
    int n = B.GetN();
    T3 val;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	{
	  SetComplexZero(val);
	  for (int ind = 0; ind < A.GetRowSize(i); ind++)
	    {
	      int k = A.Index(i, ind);
	      val += A.Value(i, ind) * B(k, j);
	    }
	  
	  C(i, j) += alpha*val;
	}
  }


  // Matrix-matrix product (sparse matrix against dense matrix)
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
	   class T4, class T2, class Prop2, class Allocator2,
           class T3, class Prop3, class Storage3, class Allocator3>
  void MltAddMatrix(const T0& alpha,
		    const Matrix<T1, Prop1, Storage1, Allocator1>& A,
		    const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		    const T4& beta,
		    Matrix<T3, Prop3, Storage3, Allocator3>& C)
  {
    if (Storage1::Sparse || Storage3::Sparse)
      throw WrongArgument("Mlt", "Function intended for product "
                          " between a sparse matrix and a dense matrix");

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, "MltAdd(alpha, A, B, beta, C)");
#endif
    
    T4 zero; SetComplexZero(zero);
    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    int m = A.GetM();
    T2 val;
    for (int j = 0; j < B.GetM(); j++)
      {
	for (int ind = 0; ind < B.GetRowSize(j); ind++)
	  {
	    int k = B.Index(j, ind);
	    val = B.Value(j, ind);
	    for (int i = 0; i < m; i++)
	      C(i, k) += alpha*val*A(i, j);
	  }
      }
  }


  // Mlt //
  /////////
  
  
  ///////////
  // Norms //
  
  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        res = max(res, ComplexAbs(A.Value(i, j)));
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal> sum(A.GetN());
    sum.Fill(Treal(0));
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        sum(A.Index(i, j)) += ComplexAbs( A.Value(i, j));
    
    return sum.GetNormInf();
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal res(0), sum;
    for (int i = 0; i < A.GetM(); i++)
      {
        sum = Treal(0);
        for (int j = 0; j < A.GetRowSize(i); j++)
          sum += ComplexAbs(A.Value(i, j));
        
        res = max(res, sum);
      }
    
    return res;
  }

  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayColSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        res = max(res, ComplexAbs(A.Value(i, j)));
    
    return res;
  }
  
  
  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayColSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal res(0), sum;
    for (int i = 0; i < A.GetN(); i++)
      {
        sum = Treal(0);
        for (int j = 0; j < A.GetColumnSize(i); j++)
          sum += ComplexAbs(A.Value(i, j));
        
        res = max(res, sum);
      }
    
    return res;
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayColSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal> sum(A.GetM());
    sum.Fill(Treal(0));
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        sum(A.Index(i, j)) += ComplexAbs(A.Value(i, j));
    
    return sum.GetNormInf();
  }
  
  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        res = max(res, ComplexAbs(A.Value(i, j)));
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal> sum(A.GetN());
    sum.Fill(Treal(0));
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          sum(A.Index(i, j)) += ComplexAbs( A.Value(i, j));
          if (A.Index(i, j) != i)
            sum(i) += ComplexAbs(A.Value(i, j));
        }
    
    return sum.GetNormInf();
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    return Norm1(A);
  }
  
  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        res = max(res, ComplexAbs(A.Value(i, j)));
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal> sum(A.GetN());
    sum.Fill(Treal(0));
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        {
          sum(A.Index(i, j)) += ComplexAbs( A.Value(i, j));
          if (A.Index(i, j) != i)
            sum(i) += ComplexAbs(A.Value(i, j));
        }
    
    return sum.GetNormInf();
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A)
  {
    return Norm1(A);
  }
  
  
  // Norms //
  ///////////
  
  
  ///////////////
  // Transpose //
  
  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                 Matrix<T, General, ArrayRowSparse, Allocator>& B)
  {
    B.Clear();
    
    int m = A.GetM();
    int n = A.GetN();
    Vector<int> ptr_T(n);
    
    B.Reallocate(n, m);

    // For each column j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    ptr_T.Zero();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        ptr_T(A.Index(i, j))++;
    
    for (int i = 0; i < n; i++)
      B.ReallocateRow(i, ptr_T(i));
    
    // filling matrix B
    ptr_T.Zero();
    for (int i = 0; i < m; i++)
      for (int jp = 0; jp < A.GetRowSize(i); jp++)
    	{
	  int j = A.Index(i, jp);
	  int k = ptr_T(j);
	  ++ptr_T(j);
          B.Value(j, k) = A.Value(i, jp);
          B.Index(j, k) = i;
    	}
    
    B.Assemble();
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayRowSparse, Allocator>& A)
  {
    Matrix<T, General, ArrayRowSparse, Allocator> Acopy(A);
    Transpose(Acopy, A);
  }
  

  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayColSparse, Allocator>& A,
                 Matrix<T, General, ArrayColSparse, Allocator>& B)
  {
    B.Clear();
    
    int m = A.GetM();
    int n = A.GetN();
    Vector<int> ptr_T(m);
    
    B.Reallocate(n, m);

    // For each row j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    ptr_T.Zero();
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        ptr_T(A.Index(i, j))++;
    
    for (int i = 0; i < m; i++)
      B.ReallocateColumn(i, ptr_T(i));
    
    // filling matrix B
    ptr_T.Zero();
    for (int i = 0; i < n; i++)
      for (int jp = 0; jp < A.GetColumnSize(i); jp++)
    	{
	  int j = A.Index(i, jp);
	  int k = ptr_T(j);
	  ++ptr_T(j);
          B.Value(j, k) = A.Value(i, jp);
          B.Index(j, k) = i;
    	}
    
    B.Assemble();
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayColSparse, Allocator>& A)
  {
    Matrix<T, General, ArrayColSparse, Allocator> Acopy(A);
    Transpose(Acopy, A);
  }


  //! Replaces A by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, General, ArrayRowSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        A.Value(i, j) = conjugate(A.Value(i, j));
  }


  //! Replaces A by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, General, ArrayColSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        A.Value(i, j) = conjugate(A.Value(i, j));
  }


  //! Replaces A by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        A.Value(i, j) = conjugate(A.Value(i, j));
  }


  //! Replaces A by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, Symmetric, ArrayColSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        A.Value(i, j) = conjugate(A.Value(i, j));
  }

  
  // Transpose //
  ///////////////


  ////////////////////////////////////
  // MltAdd (matrix-matrix product) //  


  template<class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2,
           class T4, class Prop4, class Allocator4>
  void MltMatrix(const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		 const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		 Matrix<T4, Prop4, ArrayRowSparse, Allocator4>& C)
  {
    T4 one, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    MltAdd(one, A, B, zero, C);
  }
  
  
  // C = beta * C + alpha * A * B
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2, class T3,
           class T4, class Prop4, class Allocator4>
  void MltAddMatrix(const T0& alpha,
		    const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		    const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		    const T3& beta,
		    Matrix<T4, Prop4, ArrayRowSparse, Allocator4>& C)
  {
    int m = A.GetM();
    int n = B.GetN();
    
    T3 zero;
    SetComplexZero(zero);
    if (beta == zero)
      {
        C.Clear();
        C.Reallocate(m, n);
      }
    else
      Mlt(beta, C);
    
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif
    
    Vector<int> Index(n), IndCol(n);
    Vector<T4, VectFull, Allocator4> Value(n);    
    Index.Fill(-1);
    int col, ind; T4 vloc;
    // loop over rows of matrix A
    for (int i = 0; i < m; i++)
      {
        int nnz = 0;
        for (int j = 0; j < A.GetRowSize(i); j++)
          {
            col = A.Index(i, j);
            // for each non-zero entry of the row i of A
            // loop over non-zeros entries of the corresponding row of B
            for (int k = 0; k < B.GetRowSize(col); k++)
              {
                ind = B.Index(col, k);
                vloc = alpha*A.Value(i, j)*B.Value(col, k);
                if (Index(ind) >= 0)
                  {
                    // already existing entry, we add it
                    Value(Index(ind)) += vloc;
                  }
                else
                  {
                    // new non-zero entry
                    Value(nnz) = vloc;
                    IndCol(nnz) = ind;
                    Index(ind) = nnz++;
                  }
              }
          }
        
        // sorting entries
        Sort(nnz, IndCol, Value);
        
        // adding interactions to matrix C
        C.AddInteractionRow(i, nnz, IndCol, Value, true);
        
        // resetting Index
        for (int j = 0; j < nnz; j++)
          Index(IndCol(j)) = -1;
      }
  }


  // C = beta * C + alpha * A * B
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2, class T3,
           class T4, class Prop4, class Allocator4>
  void MltAddMatrix(const T0& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		    const T3& beta,
		    Matrix<T4, Prop4, ArrayRowSparse, Allocator4>& C)
  {
    int m = A.GetM();
    int n = B.GetN();
    if (!TransA.NoTrans())
      m = A.GetN();
    
    if (!TransB.NoTrans())
      n = B.GetM();
        
    T3 zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    if (beta == zero)
      {
        C.Clear();
        C.Reallocate(m, n);
      }
    else
      {
	if ((!TransA.NoTrans()) || (!TransB.NoTrans()))
	  Mlt(beta, C);
      }

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    int col, ind; T4 vloc;
    
    if (TransA.NoTrans())
      {
	if (TransB.NoTrans())
	  MltAdd(alpha, A, B, beta, C);
	else if (TransB.Trans())
	  {
	    Vector<int> Index(n), IndCol(n);
	    Vector<T4, VectFull, Allocator4> Value(n);    
	    Index.Fill(-1);
	    
	    // loop over rows of matrix A
	    for (int i = 0; i < A.GetM(); i++)
	      {
		int nnz = 0;
		for (int j = 0; j < A.GetRowSize(i); j++)
		  {
		    col = A.Index(i, j);
		    // for each non-zero entry of the row i of A
		    // loop over non-zeros entries of the corresponding row of B
		    for (int ib = 0; ib < B.GetM(); ib++)
		      for (int k = 0; k < B.GetRowSize(ib); k++)
			if (B.Index(ib, k) == col)
			  {
			    ind = ib;
			    vloc = alpha*A.Value(i, j)*B.Value(ib, k);
			    if (Index(ind) >= 0)
			      {
				// already existing entry, we add it
				Value(Index(ind)) += vloc;
			      }
			    else
			      {
				// new non-zero entry
				Value(nnz) = vloc;
				IndCol(nnz) = ind;
				Index(ind) = nnz++;
			      }
			  }
		  }
		
		// sorting entries
		Sort(nnz, IndCol, Value);
		
		// adding interactions to matrix C
		C.AddInteractionRow(i, nnz, IndCol, Value, true);
		
		// resetting Index
		for (int j = 0; j < nnz; j++)
		  Index(IndCol(j)) = -1;
	      }
	  }
	else
	  {
	    Vector<int> Index(n), IndCol(n);
	    Vector<T4, VectFull, Allocator4> Value(n);    
	    Index.Fill(-1);
	    
	    // loop over rows of matrix A
	    for (int i = 0; i < A.GetM(); i++)
	      {
		int nnz = 0;
		for (int j = 0; j < A.GetRowSize(i); j++)
		  {
		    col = A.Index(i, j);
		    // for each non-zero entry of the row i of A
		    // loop over non-zeros entries of the corresponding row of B
		    for (int ib = 0; ib < B.GetM(); ib++)
		      for (int k = 0; k < B.GetRowSize(ib); k++)
			if (B.Index(ib, k) == col)
			  {
			    ind = ib;
			    vloc = alpha*A.Value(i, j)*conjugate(B.Value(ib, k));
			    if (Index(ind) >= 0)
			      {
				// already existing entry, we add it
				Value(Index(ind)) += vloc;
			      }
			    else
			      {
				// new non-zero entry
				Value(nnz) = vloc;
				IndCol(nnz) = ind;
				Index(ind) = nnz++;
			      }
			  }
		  }
		
		// sorting entries
		Sort(nnz, IndCol, Value);
		
		// adding interactions to matrix C
		C.AddInteractionRow(i, nnz, IndCol, Value, true);
		
		// resetting Index
		for (int j = 0; j < nnz; j++)
		  Index(IndCol(j)) = -1;
	      }
	  }
      }
    else if (TransA.Trans())
      {
	if (TransB.NoTrans())
	  {
	    // C_{ind, col} = C_{ind, col} + \sum_i A_{i, ind}  B_{i, col}
	    for (int i = 0; i < B.GetM(); i++)
	      {
		for (int j = 0; j < B.GetRowSize(i); j++)
		  {
		    col = B.Index(i, j);
		    for (int k = 0; k < A.GetRowSize(i); k++)
		      {
			ind = A.Index(i, k);
			vloc = alpha*B.Value(i, j)*A.Value(i, k);
			C.AddInteraction(ind, col, vloc);
		      }
		  }
	      }
	  }
	else if (TransB.Trans())
	  {
	    Vector<int> Index(m), IndCol(m);
	    Vector<T4, VectFull, Allocator4> Value(m);    
	    Index.Fill(-1);
	    
	    // loop over rows of matrix B
	    for (int i = 0; i < B.GetM(); i++)
	      {
		int nnz = 0;
		for (int j = 0; j < B.GetRowSize(i); j++)
		  {
		    col = B.Index(i, j);
		    // for each non-zero entry of the row i of B
		    // loop over non-zeros entries of the corresponding row of A
		    for (int k = 0; k < A.GetRowSize(col); k++)
		      {
			ind = A.Index(col, k);
			vloc = alpha*B.Value(i, j)*A.Value(col, k);
			if (Index(ind) >= 0)
			  {
			    // already existing entry, we add it
			    Value(Index(ind)) += vloc;
			  }
			else
			  {
			    // new non-zero entry
			    Value(nnz) = vloc;
			    IndCol(nnz) = ind;
			    Index(ind) = nnz++;
			  }
		      }
		  }
		
		// sorting entries
		Sort(nnz, IndCol, Value);
		
		// adding interactions to matrix C
		C.AddInteractionColumn(i, nnz, IndCol, Value);
		
		// resetting Index
		for (int j = 0; j < nnz; j++)
		  Index(IndCol(j)) = -1;
	      }
	  }
	else
	  {
	    Vector<int> Index(m), IndCol(m);
	    Vector<T4, VectFull, Allocator4> Value(m);    
	    Index.Fill(-1);

	    // loop over rows of matrix B
	    for (int i = 0; i < B.GetM(); i++)
	      {
		int nnz = 0;
		for (int j = 0; j < B.GetRowSize(i); j++)
		  {
		    col = B.Index(i, j);
		    // for each non-zero entry of the row i of B
		    // loop over non-zeros entries of the corresponding row of A
		    for (int k = 0; k < A.GetRowSize(col); k++)
		      {
			ind = A.Index(col, k);
			vloc = alpha*conjugate(B.Value(i, j))*A.Value(col, k);
			if (Index(ind) >= 0)
			  {
			    // already existing entry, we add it
			    Value(Index(ind)) += vloc;
			  }
			else
			  {
			    // new non-zero entry
			    Value(nnz) = vloc;
			    IndCol(nnz) = ind;
			    Index(ind) = nnz++;
			  }
		      }
		  }
		
		// sorting entries
		Sort(nnz, IndCol, Value);
		
		// adding interactions to matrix C
		C.AddInteractionColumn(i, nnz, IndCol, Value);
		
		// resetting Index
		for (int j = 0; j < nnz; j++)
		  Index(IndCol(j)) = -1;
	      }	    
	  }
      }
    else if (TransA.ConjTrans())
      {
	if (TransB.NoTrans())
	  {
	    // C_{ind, col} = C_{ind, col} + \sum_i A_{i, ind}  B_{i, col}
	    for (int i = 0; i < B.GetM(); i++)
	      {
		for (int j = 0; j < B.GetRowSize(i); j++)
		  {
		    col = B.Index(i, j);
		    for (int k = 0; k < A.GetRowSize(i); k++)
		      {
			ind = A.Index(i, k);
			vloc = alpha*B.Value(i, j)*conjugate(A.Value(i, k));
			C.AddInteraction(ind, col, vloc);
		      }
		  }
	      }
	  }
	else if (TransB.Trans())
	  {
	    Vector<int> Index(m), IndCol(m);
	    Vector<T4, VectFull, Allocator4> Value(m);    
	    Index.Fill(-1);
	    
	    // loop over rows of matrix B
	    for (int i = 0; i < B.GetM(); i++)
	      {
		int nnz = 0;
		for (int j = 0; j < B.GetRowSize(i); j++)
		  {
		    col = B.Index(i, j);
		    // for each non-zero entry of the row i of B
		    // loop over non-zeros entries of the corresponding row of A
		    for (int k = 0; k < A.GetRowSize(col); k++)
		      {
			ind = A.Index(col, k);
			vloc = alpha*B.Value(i, j)*conjugate(A.Value(col, k));
			if (Index(ind) >= 0)
			  {
			    // already existing entry, we add it
			    Value(Index(ind)) += vloc;
			  }
			else
			  {
			    // new non-zero entry
			    Value(nnz) = vloc;
			    IndCol(nnz) = ind;
			    Index(ind) = nnz++;
			  }
		      }
		  }
		
		// sorting entries
		Sort(nnz, IndCol, Value);
		
		// adding interactions to matrix C
		C.AddInteractionColumn(i, nnz, IndCol, Value);
		
		// resetting Index
		for (int j = 0; j < nnz; j++)
		  Index(IndCol(j)) = -1;
	      }	    
	  }
	else
	  {
	    Vector<int> Index(m), IndCol(m);
	    Vector<T4, VectFull, Allocator4> Value(m);    
	    Index.Fill(-1);
	    
	    // loop over rows of matrix B
	    for (int i = 0; i < B.GetM(); i++)
	      {
		int nnz = 0;
		for (int j = 0; j < B.GetRowSize(i); j++)
		  {
		    col = B.Index(i, j);
		    // for each non-zero entry of the row i of B
		    // loop over non-zeros entries of the corresponding row of A
		    for (int k = 0; k < A.GetRowSize(col); k++)
		      {
			ind = A.Index(col, k);
			vloc = alpha*conjugate(B.Value(i, j))*conjugate(A.Value(col, k));
			if (Index(ind) >= 0)
			  {
			    // already existing entry, we add it
			    Value(Index(ind)) += vloc;
			  }
			else
			  {
			    // new non-zero entry
			    Value(nnz) = vloc;
			    IndCol(nnz) = ind;
			    Index(ind) = nnz++;
			  }
		      }
		  }
		
		// sorting entries
		Sort(nnz, IndCol, Value);
		
		// adding interactions to matrix C
		C.AddInteractionColumn(i, nnz, IndCol, Value);
		
		// resetting Index
		for (int j = 0; j < nnz; j++)
		  Index(IndCol(j)) = -1;
	      }
	  }
      }
  }
  
  
  // MltAdd (matrix-matrix product) //  
  ////////////////////////////////////


  //! drops non-zero entries below epsilon
  template<class T, class Prop, class  Storage, class Allocator, class T0>
  void RemoveSmallEntry(Matrix<T, Prop, Storage, Allocator>& A,
                        const T0& epsilon)
  {
    A.RemoveSmallEntry(epsilon);
  }


  //! drops non-zero entries below epsilon
  template<class T, class Prop, class Allocator, class T0>
  void RemoveSmallEntry(Matrix<T, Prop, RowSparse, Allocator>& A,
                        const T0& epsilon)
  {
    // TO BE DONE
  }


  //! drops non-zero entries below epsilon
  template<class T, class Prop, class Allocator, class T0>
  void RemoveSmallEntry(Matrix<T, Prop, RowSymSparse, Allocator>& A,
                        const T0& epsilon)
  {
    // TO BE DONE
  }


  //! clears several columns of a sparse matrix
  /*!
    \param[in] col_number numbers of the columns to be cleared
    \param[inout] A sparse matrix where columns are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    int m = col_number.GetM();
    // index array to know fastly if it is a column to erase
    IVect index(A.GetM()); index.Fill(-1);
    for (int i = 0; i < m; i++)
      index(col_number(i)) = i;
    
    // first, we remove rows
    for (int i = 0; i < A.GetM(); i++)
      {
	if (index(i) != -1)
	  A.ClearRow(i);
      }
    
    // then columns
    for (int i = 0; i < A.GetM(); i++)
      {
	bool something_to_remove = false;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (index(A.Index(i,j)) != -1)
	    something_to_remove = true;
	
	if (something_to_remove)
	  {
	    int nb = 0;
	    for (int j = 0; j < A.GetRowSize(i); j++)
	      if (index(A.Index(i,j)) == -1)
		{
		  A.Index(i, nb) = A.Index(i, j);
		  A.Value(i, nb) = A.Value(i, j);	      
		  nb++;
		}
	    
	    A.ResizeRow(i, nb);
	  }
      }
  }


  //! clears several columns of a sparse matrix
  /*!
    \param[in] col_number numbers of the columns to be cleared
    \param[inout] A sparse matrix where columns are erased
   */  
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSparse, Allocator>& A)
  {
    int m = col_number.GetM();
    // index array to know fastly if it is a column to erase
    IVect index(A.GetM()); index.Fill(-1);
    for (int i = 0; i < m; i++)
      index(col_number(i)) = i;
    
    for (int i = 0; i < A.GetM(); i++)
      {
	bool something_to_remove = false;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (index(A.Index(i,j)) != -1)
	    something_to_remove = true;
	
	if (something_to_remove)
	  {
	    int nb = 0;
	    for (int j = 0; j < A.GetRowSize(i); j++)
	      if (index(A.Index(i,j)) == -1)
		{
		  A.Index(i, nb) = A.Index(i, j);
		  A.Value(i, nb) = A.Value(i, j);	      
		  nb++;
		}
	    
	    A.ResizeRow(i, nb);
	  }
      }
  }
  

  //! clears several columns of a sparse matrix
  /*!
    \param[in] col_number numbers of the columns to be cleared
    \param[inout] A sparse matrix where columns are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, RowSparse, Allocator>& A)
  {
    int m = A.GetM(), n = A.GetN();
    long nnz = A.GetIndSize();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T1* data = A.GetData();
    Vector<bool> ColToKeep(n);
    ColToKeep.Fill(true);
    for (int i = 0; i < col_number.GetM(); i++)
      ColToKeep(col_number(i)) = false;
    
    for (long i = 0; i < A.GetIndSize(); i++)
      if (!ColToKeep(ind[i]))
        nnz--;
    
    if (nnz == A.GetIndSize())
      return;
    
    Vector<long> Ptr(m+1); Vector<int> Ind(nnz);
    Vector<T1, VectFull, Allocator> Val(nnz);
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
        long jA = Ptr(i); int size_row = 0;
        for (long j = ptr[i]; j < ptr[i+1]; j++)
          if (ColToKeep(ind[j]))
            {
              Ind(jA) = ind[j];
              Val(jA) = data[j];
              size_row++; jA++;
            }
        
        Ptr(i+1) = Ptr(i) + size_row;
      }
    
    A.SetData(m, n, Val, Ptr, Ind);
  }


  //! clears several columns of a sparse matrix
  /*!
    \param[in] col_number numbers of the columns to be cleared
    \param[inout] A sparse matrix where columns are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, RowSymSparse, Allocator>& A)
  {
    int m = A.GetM(), n = A.GetN();
    long nnz = A.GetIndSize();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T1* data = A.GetData();
    Vector<bool> ColToKeep(n);
    ColToKeep.Fill(true);
    for (int i = 0; i < col_number.GetM(); i++)
      ColToKeep(col_number(i)) = false;
    
    for (int i = 0; i < m; i++)
      {
        if (!ColToKeep(i))
          nnz -= ptr[i+1] - ptr[i];
        else
          {
            for (long j = ptr[i]; j < ptr[i+1]; j++)
              if (!ColToKeep(ind[j]))
                nnz--;
          }
      }
    
    if (nnz == A.GetIndSize())
      return;
    
    Vector<long> Ptr(m+1); Vector<int> Ind(nnz);
    Vector<T1, VectFull, Allocator> Val(nnz);
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
        long jA = Ptr(i); int size_row = 0;
        if (ColToKeep(i))            
          for (long j = ptr[i]; j < ptr[i+1]; j++)
            if (ColToKeep(ind[j]))
              {
                Ind(jA) = ind[j];
                Val(jA) = data[j];
                size_row++; jA++;
              }
        
        Ptr(i+1) = Ptr(i) + size_row;
      }
    
    A.SetData(m, n, Val, Ptr, Ind);
  }
  

  //! clears several columns of a sparse matrix
  /*!
    \param[in] col_number numbers of the columns to be cleared
    \param[inout] A sparse matrix where columns are erased
   */
  template<class T1, class Prop, class Storage, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, Storage, Allocator>& A)
  {
    cout << "Not implemented for any matrix" << endl;
    abort();
  }
   

  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */
  template<class T1, class Prop, class Storage, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, Storage, Allocator>& A)
  {
    cout << "Not implemented for any matrix" << endl;
    abort();
  }


  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */    
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    EraseCol(col_number, A);
  }
  

  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSparse, Allocator>& A)
  {
    for (int i = 0; i < col_number.GetM(); i++)
      A.ClearRow(col_number(i));
  }
  

  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, RowSparse, Allocator>& A)
  {
    int m = A.GetM(), n = A.GetN();
    long nnz = A.GetIndSize();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T1* data = A.GetData();
    Vector<bool> RowToKeep(m);
    RowToKeep.Fill(true);
    for (int i = 0; i < col_number.GetM(); i++)
      RowToKeep(col_number(i)) = false;
    
    for (int i = 0; i < m; i++)
      if (!RowToKeep(i))
        nnz -= ptr[i+1] - ptr[i];
    
    Vector<long> Ptr(m+1); Vector<int> Ind(nnz);
    Vector<T1, VectFull, Allocator> Val(nnz);
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
        if (RowToKeep(i))
          {
            int size_row = ptr[i+1] - ptr[i];
            for (int j = 0; j < size_row; j++)
              {
                Ind(Ptr(i) + j) = ind[ptr[i] + j];
                Val(Ptr(i) + j) = data[ptr[i] + j];
              }
            
            Ptr(i+1) = Ptr(i) + size_row;
          }
        else
          Ptr(i+1) = Ptr(i);
      }
    
    A.SetData(m, n, Val, Ptr, Ind);
  }
  
  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, RowSymSparse, Allocator>& A)
  {
    EraseCol(col_number, A);
  }
  

  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat_direct given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General,
		 ArrayRowSparse, Allocator> & mat_direct)
  {
    int n = mat_direct.GetM();
    diagonal_scale_left.Reallocate(n);
    diagonal_scale_left.Fill(0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < mat_direct.GetRowSize(i); j++)
	diagonal_scale_left(i) += abs(mat_direct.Value(i,j));
    
  }
  
  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat_direct given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left, 
		 const Matrix<Complexe, Symmetric,
		 ArrayRowSymSparse, Allocator> & mat_direct)
  {
    int n = mat_direct.GetM();
    diagonal_scale_left.Reallocate(n);
    diagonal_scale_left.Fill(0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < mat_direct.GetRowSize(i); j++)
	{
	  diagonal_scale_left(i) += abs(mat_direct.Value(i,j));
	  if (i != mat_direct.Index(i,j))
	    diagonal_scale_left(mat_direct.Index(i,j))
	      += abs(mat_direct.Value(i,j));
	}
  }
  
  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat_direct given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General,
		 RowSparse, Allocator> & mat_direct)
  {
    int n = mat_direct.GetM();
    diagonal_scale_left.Reallocate(n);
    diagonal_scale_left.Fill(0);
    long* ptr = mat_direct.GetPtr();
    Complexe* data = mat_direct.GetData();
    for (int i = 0; i < n; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
	diagonal_scale_left(i) += abs(data[j]);
    
  }
  
  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat_direct given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left, 
		 const Matrix<Complexe, Symmetric,
		 RowSymSparse, Allocator> & mat_direct)
  {
    int n = mat_direct.GetM();
    diagonal_scale_left.Reallocate(n);
    diagonal_scale_left.Fill(0);
    long* ptr = mat_direct.GetPtr();
    int* ind = mat_direct.GetInd();
    Complexe* data = mat_direct.GetData();
    for (int i = 0; i < n; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
	{
	  diagonal_scale_left(i) += abs(data[j]);
	  if (i != ind[j])
	    diagonal_scale_left(ind[j]) += abs(data[j]);
	}
  }
  
  
  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat_direct given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale,
		 const Matrix<Complexe, General,
		 ArrayRowSparse, Allocator> & mat_direct)
  {
    int n = mat_direct.GetM();
    diagonal_scale.Reallocate(mat_direct.GetN());
    diagonal_scale.Fill(0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < mat_direct.GetRowSize(i); j++)
	diagonal_scale(mat_direct.Index(i, j)) += abs(mat_direct.Value(i, j));
    
  }
  
  
  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat_direct given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale, 
		 const Matrix<Complexe, Symmetric,
		 ArrayRowSymSparse, Allocator> & mat_direct)
  {
    GetRowSum(diagonal_scale, mat_direct);
  }
  
  
  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat_direct given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale,
		 const Matrix<Complexe, General,
		 RowSparse, Allocator> & mat_direct)
  {
    int n = mat_direct.GetM();
    diagonal_scale.Reallocate(mat_direct.GetN());
    diagonal_scale.Fill(0);
    long* ptr = mat_direct.GetPtr();
    int* ind = mat_direct.GetInd();
    Complexe* data = mat_direct.GetData();
    for (int i = 0; i < n; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
	diagonal_scale(ind[j]) += abs(data[j]);
    
  }
  
  
  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat_direct given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale, 
		 const Matrix<Complexe, Symmetric,
		 RowSymSparse, Allocator> & mat_direct)
  {
    GetRowSum(diagonal_scale, mat_direct);
  }
  

  //! For each row and column of the matrix,
  //! computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] A given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, General,
                    ArrayRowSparse, Allocator> & A)
  {
    int n = A.GetM();
    sum_row.Reallocate(n);
    sum_col.Reallocate(A.GetN());
    sum_row.Fill(0);
    sum_col.Fill(0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	{
          sum_row(i) += abs(A.Value(i,j));
          sum_col(A.Index(i, j)) += abs(A.Value(i,j));
        }    
  }
  

  //! For each row and column of the matrix
  //! computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] A given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, Symmetric,
                    ArrayRowSymSparse, Allocator> & A)
  {
    GetRowSum(sum_row, A);
    sum_col = sum_row;
  }

  
  //! For each row and column of the matrix,
  //! computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] A given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, General,
                    RowSparse, Allocator> & A)
  {
    int n = A.GetM();
    sum_row.Reallocate(n);
    sum_col.Reallocate(A.GetN());
    sum_row.Fill(0);
    sum_col.Fill(0);
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    Complexe* data = A.GetData();
    for (int i = 0; i < n; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
	{
          sum_row(i) += abs(data[j]);
          sum_col(ind[j]) += abs(data[j]);
        }    
  }
  

  //! For each row and column of the matrix,
  //! computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] A given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, Symmetric,
                    RowSymSparse, Allocator> & A)
  {
    GetRowSum(sum_row, A);
    sum_col = sum_row;
  }


  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<T1, VectFull, Allocator1>& Value)
  {
    int m = A.GetM(), n = A.GetN();
    if ((m <= 0) || (n <= 0) || (row.GetM() <= 0) || (col.GetM() <= 0))
      {
        RowNum.Clear(); ColNum.Clear(); Value.Clear();
        return;
      }
    
    Vector<bool> RowKept(m), ColKept(n);
    RowKept.Fill(false); ColKept.Fill(false);
    for (int i = 0; i < row.GetM(); i++)
      RowKept(row(i)) = true;
    
    for (int i = 0; i < col.GetM(); i++)
      ColKept(col(i)) = true;
    
    long nnz = 0;
    for (int i = 0; i < A.GetM(); i++)
      if (RowKept(i))
        for (int j = 0; j < A.GetRowSize(i); j++)
          if (ColKept(A.Index(i, j)))
            nnz++;
    
    RowNum.Reallocate(nnz);
    ColNum.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < A.GetM(); i++)
      if (RowKept(i))
        for (int j = 0; j < A.GetRowSize(i); j++)
          if (ColKept(A.Index(i, j)))
            {
              RowNum(nnz) = i;
              ColNum(nnz) = A.Index(i, j);
              Value(nnz) = A.Value(i, j);
              nnz++;
            }
  }
  

  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0,
		     ArrayRowSymSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<T1, VectFull, Allocator1>& Value)
  {
    int m = A.GetM(), n = A.GetN();
    if ((m <= 0) || (n <= 0) || (row.GetM() <= 0) || (col.GetM() <= 0))
      {
        RowNum.Clear(); ColNum.Clear(); Value.Clear();
        return;
      }
    
    Vector<bool> RowKept(m), ColKept(n);
    RowKept.Fill(false); ColKept.Fill(false);
    for (int i = 0; i < row.GetM(); i++)
      RowKept(row(i)) = true;
    
    for (int i = 0; i < col.GetM(); i++)
      ColKept(col(i)) = true;
    
    long nnz = 0;
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          if (ColKept(A.Index(i, j)) && RowKept(i))
            nnz++;
          
          if (A.Index(i, j) != i)
            if (RowKept(A.Index(i, j)) && ColKept(i))
              nnz++;
        }
    
    RowNum.Reallocate(nnz);
    ColNum.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          if (ColKept(A.Index(i, j)) && RowKept(i))
            {
              RowNum(nnz) = i;
              ColNum(nnz) = A.Index(i, j);
              Value(nnz) = A.Value(i, j);
              nnz++;
            }
          
          if (A.Index(i, j) != i)
            if (RowKept(A.Index(i, j)) && ColKept(i))
              {
                RowNum(nnz) = A.Index(i, j);
                ColNum(nnz) = i;
                Value(nnz) = A.Value(i, j);
                nnz++;
              }
        }
  }


  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<T1, VectFull, Allocator1>& Value)
  {
    int m = A.GetM(), n = A.GetN();
    if ((m <= 0) || (n <= 0) || (row.GetM() <= 0) || (col.GetM() <= 0))
      {
        RowNum.Clear(); ColNum.Clear(); Value.Clear();
        return;
      }
    
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T0* data = A.GetData();
    Vector<bool> RowKept(m), ColKept(n);
    RowKept.Fill(false); ColKept.Fill(false);
    for (int i = 0; i < row.GetM(); i++)
      RowKept(row(i)) = true;
    
    for (int i = 0; i < col.GetM(); i++)
      ColKept(col(i)) = true;
    
    long nnz = 0;
    for (int i = 0; i < m; i++)
      if (RowKept(i))
        for (long j = ptr[i]; j < ptr[i+1]; j++)
          if (ColKept(ind[j]))
            nnz++;
    
    RowNum.Reallocate(nnz);
    ColNum.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < m; i++)
      if (RowKept(i))
        for (long j = ptr[i]; j < ptr[i+1]; j++)
          if (ColKept(ind[j]))
            {
              RowNum(nnz) = i;
              ColNum(nnz) = ind[j];
              Value(nnz) = data[j];
              nnz++;
            }
  }
  

  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<T1, VectFull, Allocator1>& Value)
  {
    int m = A.GetM(), n = A.GetN();
    if ((m <= 0) || (n <= 0) || (row.GetM() <= 0) || (col.GetM() <= 0))
      {
        RowNum.Clear(); ColNum.Clear(); Value.Clear();
        return;
      }

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T0* data = A.GetData();    
    Vector<bool> RowKept(m), ColKept(n);
    RowKept.Fill(false); ColKept.Fill(false);
    for (int i = 0; i < row.GetM(); i++)
      RowKept(row(i)) = true;
    
    for (int i = 0; i < col.GetM(); i++)
      ColKept(col(i)) = true;
    
    long nnz = 0;
    for (int i = 0; i < m; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
        {
          if (ColKept(ind[j]) && RowKept(i))
            nnz++;
          
          if (ind[j] != i)
            if (RowKept(ind[j]) && ColKept(i))
              nnz++;
        }
    
    RowNum.Reallocate(nnz);
    ColNum.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < m; i++)
      for (long j = ptr[i]; j < ptr[i+1]; j++)
        {
          if (ColKept(ind[j]) && RowKept(i))
            {
              RowNum(nnz) = i;
              ColNum(nnz) = ind[j];
              Value(nnz) = data[j];
              nnz++;
            }
          
          if (ind[j] != i)
            if (RowKept(ind[j]) && ColKept(i))
              {
                RowNum(nnz) = ind[j];
                ColNum(nnz) = i;
                Value(nnz) = data[j];
                nnz++;
              }
        }
  }


  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Matrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    Vector<int> RowNum, ColNum;
    typedef typename Matrix<T0, Prop0, Storage0, Allocator0>::entry_type T;
    Vector<T> Value;
    
    // extracts rows/columns in coordinate format
    CopySubMatrix(A, row, col, RowNum, ColNum, Value);
    
    // converts to the sparse matrix B
    B.Reallocate(A.GetM(), A.GetN());
    ConvertMatrix_from_Coordinates(RowNum, ColNum, Value, B);
  }

} // namespace Seldon

#define SELDON_FILE_FUNCTIONS_MATRIX_ARRAY_CXX
#endif
