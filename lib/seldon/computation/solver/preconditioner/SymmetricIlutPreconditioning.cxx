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


#ifndef SELDON_FILE_SYMMETRIC_ILUT_PRECONDITIONING_CXX

namespace Seldon
{

  //! Incomplete Factorization without pivot for symmetric matrix.
  template<class cplx, class Allocator1, class Allocator2>
  void GetIlut(const IlutPreconditioning<cplx, Allocator1>& param,
               Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator2>& A)
  {
    int size_row;
    int n = A.GetN();
    int lfil = param.GetFillLevel();
    typename ClassComplexType<cplx>::Treal zero = 0.0;
    typename ClassComplexType<cplx>::Treal droptol = param.GetDroppingThreshold();
    typename ClassComplexType<cplx>::Treal alpha = param.GetDiagonalCoefficient();
    bool variable_fill = false;
    bool standard_dropping = true;
    int type_factorization = param.GetFactorisationType();
    int additional_fill = param.GetAdditionalFillNumber();
    int print_level = param.GetPrintLevel();
    if (type_factorization == param.ILUT)
      standard_dropping = false;
    else if (type_factorization == param.ILU_D)
      standard_dropping = true;
    else if (type_factorization == param.ILUT_K)
      {
        variable_fill = true;
	standard_dropping = false;
      }
    else if (type_factorization == param.ILU_0)
      {
        GetIlu0(A);
        return;
      }
    else if (type_factorization == param.MILU_0)
      {
	GetMilu0(A);
        return;
      }
    else if (type_factorization == param.ILU_K)
      {
	GetIluk(lfil, print_level, A);
        return;
      }

    cplx fact, s, t;
    typename ClassComplexType<cplx>::Treal tnorm, one(1);
    int length_lower, length_upper, jpos, jrow;
    int i_row, j_col, index_lu, length;
    int i, j, k;

    lfil = n;
    typedef Vector<cplx, VectFull, Allocator2> VectCplx;
    VectCplx Row_Val(n);
    IVect Index(n), Row_Ind(n);
    Row_Val.Fill(0); Row_Ind.Fill(-1);

    Index.Fill(-1);

    bool element_dropped; cplx dropsum, czero;
    SetComplexZero(czero);

    // We convert A into an unsymmetric matrix.
    Matrix<cplx, General, ArrayRowSparse, Allocator2> B;
    Seldon::Copy(A, B);

    A.Clear();
    A.Reallocate(n, n);

    // Main loop.
    int new_percent = 0, old_percent = 0;
    for (i_row = 0; i_row < n; i_row++)
      {
        // Progress bar if print level is high enough.
        if (print_level > 0)
          {
            new_percent = int(double(i_row)/(n-1)*80);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#"; cout.flush();
              }

            old_percent = new_percent;
          }

        // 1-norm of the row of initial matrix.
	size_row = B.GetRowSize(i_row);
	tnorm = zero;
	dropsum = czero;
	for (k = 0 ; k < size_row; k++)
          tnorm += abs(B.Value(i_row, k));

	if (tnorm == zero)
	  {
            cout << "Structurally singular matrix." << endl;
            cout << "Norm of row " << i_row << " is equal to 0." << endl;
            abort();
          }

        // tnorm is the sum of absolute value of coefficients of row i_row
	tnorm /= typename ClassComplexType<cplx>::Treal(size_row);
	if (variable_fill)
	  lfil = size_row + additional_fill;


        // Separating lower part from upper part for this row.
	length_upper = 1;
	length_lower = 0;
	Row_Ind(i_row) = i_row;
	Row_Val(i_row) = czero;
	Index(i_row) = i_row;

	for (j = 0; j < size_row; j++)
	  {
	    k = B.Index(i_row,j);
            t = B.Value(i_row,j);
	    if (k < i_row)
	      {
		Row_Ind(length_lower) = k;
		Row_Val(length_lower) = t;
		Index(k) = length_lower;
		++length_lower;
	      }
	    else if (k == i_row)
	      {
		Row_Val(i_row) = t;
	      }
	    else
	      {
		jpos = i_row + length_upper;
		Row_Ind(jpos) = k;
		Row_Val(jpos) = t;
		Index(k) = jpos;
		length_upper++;
	      }
          }

        // This row of B is cleared.
        B.ClearRow(i_row);

	j_col = 0;
	length = 0;
	
        // We eliminate previous rows.
	while (j_col <length_lower)
	  {
	    // In order to do the elimination in the correct order, we must
            // select the smallest column index.
	    jrow = Row_Ind(j_col);
	    k = j_col;

	    // We determine smallest column index.
	    for (j = (j_col+1) ; j < length_lower; j++)
	      {
		if (Row_Ind(j) < jrow)
		  {
		    jrow = Row_Ind(j);
		    k = j;
		  }
	      }

            // If needed, we exchange positions of this element in
            // Row_Ind/Row_Val so that it appears first.
	    if (k != j_col)
	      {

		j = Row_Ind(j_col);
		Row_Ind(j_col) = Row_Ind(k);
		Row_Ind(k) = j;

		Index(jrow) = j_col;
		Index(j) = k;

		s = Row_Val(j_col);
		Row_Val(j_col) = Row_Val(k);
	        Row_Val(k) = s;
	      }

            // Zero out element in row by setting Index to -1.
	    Index(jrow) = -1;

	    element_dropped = false;
	    if (standard_dropping)
	      if (abs(Row_Val(j_col)) <= droptol*tnorm)
		{
		  dropsum += Row_Val(j_col);
		  element_dropped = true;
		}

            // Gets the multiplier for row to be eliminated.
	    if (!element_dropped)
	      {
                fact = Row_Val(j_col) * A.Value(jrow, 0);

		if (!standard_dropping)
		  {
		    if (abs(fact) <= droptol)
		      element_dropped = true;
		  }
	      }

	    if (!element_dropped)
	      {
		// Combines current row and row jrow.
		for (k = 1; k < A.GetRowSize(jrow); k++)
	 	  {
		    s = fact * A.Value(jrow,k);
		    j = A.Index(jrow,k);

		    jpos = Index(j);
		    if (j >= i_row)
		      {

			// Dealing with upper part.
			if (jpos == -1)
			  {

			    // This is a fill-in element.
			    i = i_row + length_upper;
			    Row_Ind(i) = j;
			    Index(j) = i;
			    Row_Val(i) = -s;
			    ++length_upper;
			  }
			else
			  {
			    // This is not a fill-in element.
			    Row_Val(jpos) -= s;
			  }
		      }
		    else
		      {
                        // Dealing  with lower part.
			if (jpos == -1)
			  {
                            // This is a fill-in element.
                            Row_Ind(length_lower) = j;
			    Index(j) = length_lower;
			    Row_Val(length_lower) = -s;
			    ++length_lower;
			  }
			else
			  {
			    // This is not a fill-in element.
			    Row_Val(jpos) -= s;
			  }
		      }
		  }

		// We store this pivot element (from left to right -- no
		// danger of overlap with the working elements in L (pivots).
		Row_Val(length) = fact;
		Row_Ind(length) = jrow;
		++length;
	      }

	    j_col++;
	  }

	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          Index(Row_Ind(i_row + k)) = -1;

	// Updating U-matrix -- first apply dropping strategy.
	length = 0;
	for (k = 1; k <= (length_upper-1); k++)
	  {
	    if (abs(Row_Val(i_row+k)) > droptol * tnorm)
	      {
		++length;
		Row_Val(i_row+length) = Row_Val(i_row+k);
		Row_Ind(i_row+length) = Row_Ind(i_row+k);
	      }
	    else
	      dropsum += Row_Val(i_row+k);
	  }


	if (!standard_dropping)
	  {
	    length_upper = length + 1;
	    length = min(length_upper, lfil);

            qsplit_ilut(Row_Val, Row_Ind, i_row+1,
                        i_row+length_upper, i_row+length+1, tnorm);
          }
	else
	  length++;

	// Copies U-part in matrix A.
	A.ReallocateRow(i_row, length);
	index_lu = 1;
	// sorting column numbers
	Sort(i_row+1, i_row + length - 1, Row_Ind, Row_Val);
	for (k = (i_row+1) ; k <= (i_row+length-1) ; k++)
	  {
	    A.Index(i_row,index_lu) = Row_Ind(k);
	    A.Value(i_row,index_lu) = Row_Val(k);
	    ++index_lu;
	  }

	// Stores the inverse of the diagonal element of u.
	if (standard_dropping)
	  Row_Val(i_row) += alpha*dropsum;

	if (Row_Val(i_row) == czero)
          Row_Val(i_row) = (droptol + 1e-4) * tnorm;

	A.Value(i_row,0) = one / Row_Val(i_row);

      } // end main loop.

    if (print_level > 0)
      cout<<endl;

    // for each row of A, we divide by diagonal value
    for (int i = 0; i < n; i++)
      for (int j = 1; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= A.Value(i,0);

    if (print_level > 0)
      cout << "The matrix takes " <<
        int((A.GetDataSize()*(sizeof(cplx)+4))/(1024*1024)) << " MB" << endl;

  }

  
  //! Basic ilu(k) solver 
  /*!
    \param lfil level k
    \param print_level if print_level > 0, messages are displayed
    \param A on input the matrix, on output factors L and U
  */
  template<class cplx, class Allocator>
  void GetIluk(int lfil, int print_level,
               Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    int n = A.GetM();
    cplx fact, s, t;
    int length_lower, length_upper, jpos, jrow;
    int i_row, j_col, index_lu, length;
    int i, j, k;
    typename ClassComplexType<cplx>::Treal tnorm, one(1);
    
    if (lfil < 0)
      {
        cout << "Incorrect fill level." << endl;
        abort();
      }

    typedef Vector<cplx, VectFull, Allocator> VectCplx;
    VectCplx Row_Val(n);
    IVect Index(n), Row_Ind(n), Row_Level(n);
    Row_Val.Fill(0); Row_Ind.Fill(-1);
    Row_Level.Fill(-1);
    Index.Fill(-1);

    bool element_dropped;

    // We convert A into an unsymmetric matrix.
    Matrix<cplx, General, ArrayRowSparse, Allocator> B;
    Seldon::Copy(A, B);
    
    A.Clear();
    A.Reallocate(n, n);
    Vector<IVect, VectFull, NewAlloc<IVect> > levs(n);
    
    // Main loop.
    int new_percent = 0, old_percent = 0;
    for (i_row = 0; i_row < n; i_row++)
      {
        // Progress bar if print level is high enough.
        if (print_level > 0)
          {
            new_percent = int(double(i_row)/(n-1)*80);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#"; cout.flush();
              }

            old_percent = new_percent;
          }

        // 1-norm of the row of initial matrix.
	int size_row = B.GetRowSize(i_row);
	tnorm = 0.0;
	for (k = 0 ; k < size_row; k++)
          tnorm += abs(B.Value(i_row, k));
	
	if (tnorm == 0.0)
	  {
            cout << "Structurally singular matrix." << endl;
            cout << "Norm of row " << i_row << " is equal to 0." << endl;
            abort();
          }

	// current row will be stored in arrays Row_Val, Row_Ind 
	// (respectively values and column indexes)
	// the array Index returns for global column index
	// the local column index used in Row_Ind
	// (Index is the reciprocal array of Row_Ind)	
	
	// Row_Level will store the level associated with each
	// column index (-1 for an entry in the original matrix,
	// 0 for an entry appearing because of a single combination,
	// 1 for an entry appearing because of at least two combinations, etc)
	
        // Separating lower part from upper part for this row.
	length_upper = 1;
	length_lower = 0;
	
	Row_Ind(i_row) = i_row;
	Row_Val(i_row) = 0.0;
	Index(i_row) = i_row;

	for (j = 0; j < size_row; j++)
	  {
	    k = B.Index(i_row, j);
            t = B.Value(i_row, j);
	    if (k < i_row)
	      {
		// part in L
		Row_Ind(length_lower) = k;
		Row_Val(length_lower) = t;
		Index(k) = length_lower;
		Row_Level(length_lower) = -1;
		++length_lower;
	      }
	    else if (k == i_row)
	      {
		// diagonal part
		Row_Val(i_row) = t;
		Row_Level(i_row) = -1;
	      }
	    else
	      {
		// part in U
		jpos = i_row + length_upper;
		Row_Ind(jpos) = k;
		Row_Val(jpos) = t;
		Row_Level(jpos) = -1;
		Index(k) = jpos;
		length_upper++;
	      }
          }

        // This row of B is cleared (since already stored in Row_Ind/Row_Val)
        B.ClearRow(i_row);

	j_col = 0;
	length = 0;

        // We eliminate previous rows.
	while (j_col < length_lower)
	  {
	    // In order to do the elimination in the correct order, we must
            // select the smallest column index.
	    jrow = Row_Ind(j_col);
	    k = j_col;

	    // We determine smallest column index.
	    for (j = (j_col+1) ; j < length_lower; j++)
	      {
		if (Row_Ind(j) < jrow)
		  {
		    jrow = Row_Ind(j);
		    k = j;		    
		  }
	      }

            // If needed, we exchange positions of this element in
            // Row_Ind/Row_Val so that it appears first.
	    if (k != j_col)
	      {

		j = Row_Ind(j_col);
		Row_Ind(j_col) = Row_Ind(k);
		Row_Ind(k) = j;

		Index(jrow) = j_col;
		Index(j) = k;
		
		j = Row_Level(j_col);
		Row_Level(j_col) = Row_Level(k);
		Row_Level(k) = j;

		s = Row_Val(j_col);
		Row_Val(j_col) = Row_Val(k);
	        Row_Val(k) = s;
	      }

            // Zero out element in row by setting Index to -1.
	    Index(jrow) = -1;

	    element_dropped = false;
	    fact = Row_Val(j_col) * A.Value(jrow, 0);
	    
	    int jlev = Row_Level(j_col) + 1;
	    if (jlev > lfil)
	      element_dropped = true;
	    
	    if (!element_dropped)
	      {
		// Combines current row and row jrow.
		for (k = 1; k < A.GetRowSize(jrow); k++)
	 	  {
		    s = fact * A.Value(jrow, k);
		    j = A.Index(jrow, k);

		    jpos = Index(j);
		    if (j >= i_row)
		      {
			// Dealing with upper part.
			if (jpos == -1)
			  {
			    // This is a fill-in element.
			    i = i_row + length_upper;
			    Row_Ind(i) = j;
			    Index(j) = i;
			    Row_Val(i) = -s;
			    Row_Level(i) = jlev + levs(jrow)(k) + 1;
			    ++length_upper;
			  }
			else
			  {
			    // This is not a fill-in element.
			    Row_Val(jpos) -= s;
			    Row_Level(jpos) = min(Row_Level(jpos),
						  jlev + levs(jrow)(k)+1);
			  }
		      }
		    else
		      {
                        // Dealing  with lower part.
			if (jpos == -1)
			  {
                            // This is a fill-in element.
                            Row_Ind(length_lower) = j;
			    Index(j) = length_lower;
			    Row_Val(length_lower) = -s;
			    Row_Level(length_lower) = jlev + levs(jrow)(k) + 1;
			    ++length_lower;
			  }
			else
			  {
			    // This is not a fill-in element.
			    Row_Val(jpos) -= s;
			    Row_Level(jpos) = min(Row_Level(jpos),
						  jlev + levs(jrow)(k)+1);
			  }
		      }
		  }

		// We store this pivot element (from left to right -- no
		// danger of overlap with the working elements in L (pivots).
		Row_Val(length) = fact;
		Row_Ind(length) = jrow;
		++length;
	      }

	    j_col++;
	  }
	
	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          Index(Row_Ind(i_row+k )) = -1;

	// couting size of upper part after dropping
	length = 1;
	for (k = 1; k <= (length_upper-1); k++)
	  if (Row_Level(i_row+k) < lfil)
	    length++;
	
	// sorting column indexes in U
	Sort(i_row+1, i_row+length_upper-1, Row_Ind, Row_Val, Row_Level);
	
	// Copies U-part in matrix A.
	// L-part is not stored since A is symmetric
	A.ReallocateRow(i_row, length);
	levs(i_row).Reallocate(length);
	
	// diagonal element
	A.Index(i_row, 0) = i_row;
	A.Value(i_row, 0) = one / Row_Val(i_row);
	levs(i_row)(0) = -1;
	
	// extra-diagonal elements
	index_lu = 1;
	for (k = (i_row+1) ; k <= (i_row+length_upper-1) ; k++)
	  if (Row_Level(k) < lfil)
	    {
	      A.Index(i_row, index_lu) = Row_Ind(k);
	      A.Value(i_row, index_lu) = Row_Val(k);
	      levs(i_row)(index_lu) = Row_Level(k);
	      ++index_lu;
	    }
	
      } // end main loop.
    
    if (print_level > 0)
      cout<<endl;
    
    // for each row of A, we divide by diagonal value
    for (int i = 0; i < n; i++)
      for (int j = 1; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= A.Value(i,0);
    
    if (print_level > 0)
      cout << "The matrix takes " <<
        int((A.GetDataSize()*(sizeof(cplx)+4))/(1024*1024)) << " MB" << endl;

  }


  template<class cplx, class Allocator>
  void GetIlu0(Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    int n = A.GetM();
    cplx one, invDiag, fact, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    IVect Index(n);
    Index.Fill(-1);
    // loop on rows
    for (int i = 0; i < n; i++)
      {
	if (A.GetRowSize(i) == 0)
	  {
	    cout << "Empty row " << i << endl;
	    cout << "Factorisation cannot be completed" << endl;
	    abort();
	  }
	
	if (A.Index(i, 0) != i)
	  {
	    cout << "No diagonal element on row " << i << endl;
	    cout << "ILU(0) needs one" << endl;
	    abort();
	  }
	
	if (A.Value(i, 0) == zero)
	  {
	    cout << "Factorization fails because we found a null coefficient"
                 << " on diagonal " << i << endl;	    
            abort(); 
	  }

	// updating Index, Index(j) will be the local position of column j
	// in the row i (-1 if there is no non-zero entry (i, j))
	for (int jloc = 1; jloc < A.GetRowSize(i); jloc++)
	  Index(A.Index(i, jloc)) = jloc;
	
	invDiag = one / A.Value(i, 0);
	// loop on each extra-diagonal element of the row
	for (int jloc = 1; jloc < A.GetRowSize(i); jloc++)
	  {
	    int j = A.Index(i, jloc);
	    // combination Lj <- Lj - a_ji / a_ii L_i
	    // for upper part of the row Lj
	    fact = A.Value(i, jloc) * invDiag;
	    
	    for (int kloc = 0; kloc < A.GetRowSize(j); kloc++)
	      {
		// ILU(0) -> combination performed only if a non-zero entry
		// exists at position (j, k)
		int k = A.Index(j, kloc);
		if (Index(k) >= 0)
		  A.Value(j, kloc) -= fact * A.Value(i, Index(k));
	      }
	  }
	
	// storing inverse of diagonal
	A.Value(i, 0) = invDiag;
	
	// reverting Index to -1
	for (int jloc = 1; jloc < A.GetRowSize(i); jloc++)
	  Index(A.Index(i, jloc)) = -1;
	
      }

    // for each row of A, we divide by diagonal value
    for (int i = 0; i < n; i++)
      for (int j = 1; j < A.GetRowSize(i); j++)
	A.Value(i, j) *= A.Value(i, 0);
  }


  template<class cplx, class Allocator>
  void GetMilu0(Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    int n = A.GetM();
    cplx one, invDiag, fact, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    Vector<cplx> SumRow(n);
    SumRow.Fill(zero);
    // loop on rows
    for (int i = 0; i < n; i++)
      {
	if (A.GetRowSize(i) == 0)
	  {
	    cout << "Empty row " << i << endl;
	    cout << "Factorisation cannot be completed" << endl;
	    abort();
	  }
	
	if (A.Index(i, 0) != i)
	  {
	    cout << "No diagonal element on row " << i << endl;
	    cout << "ILU(0) needs one" << endl;
	    abort();
	  }
	
	// adding fill-in values to the diagonal
	A.Value(i, 0) += SumRow(i);
	if (A.Value(i, 0) == zero)
	  {
	    cout << "Factorization fails because we found a null coefficient"
                 << " on diagonal " << i << endl;	    
            abort(); 
	  }
	
	invDiag = one / A.Value(i, 0);
	// loop on each extra-diagonal element of the row
	for (int jloc = 1; jloc < A.GetRowSize(i); jloc++)
	  {
	    int j = A.Index(i, jloc);
	    // combination Lj <- Lj - a_ji / a_ii L_i
	    // for upper part of the row Lj
	    fact = A.Value(i, jloc) * invDiag;
	    
	    int k2 = 0;
	    while ((k2 < A.GetRowSize(i)) && (A.Index(i, k2) < j))
	      k2++;
	    
	    for (int kloc = 0; kloc < A.GetRowSize(j); kloc++)
	      {
		int k = A.Index(j, kloc);
		// MILU(0) -> combination performed only if a non-zero entry
		// exists at position (j, k)
		// Fill-in elements are summed and reported to the diagonal
		while ((k2 < A.GetRowSize(i)) && (A.Index(i, k2) < k))
		  SumRow(j) -= fact * A.Value(i, k2++);
		
		if ((k2 < A.GetRowSize(i)) && (A.Index(i, k2) == k))
		  A.Value(j, kloc) -= fact * A.Value(i, k2++);
	      }
	    
	    while (k2 < A.GetRowSize(i))
	      SumRow(j) -= fact * A.Value(i, k2++);
	  }
	
	// storing inverse of diagonal
	A.Value(i, 0) = invDiag;		
      }

    // for each row of A, we divide by diagonal value
    for (int i = 0; i < n; i++)
      for (int j = 1; j < A.GetRowSize(i); j++)
	A.Value(i, j) *= A.Value(i, 0);
  }


} // end namespace

#define SELDON_FILE_SYMMETRIC_ILUT_PRECONDITIONING_CXX
#endif
