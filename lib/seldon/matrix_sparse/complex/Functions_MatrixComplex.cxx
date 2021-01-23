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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_COMPLEX_CXX

/*
  Functions defined in this file:
  (storage RowComplexSparse, ArrayRowComplexSparse, etc)
  
  alpha.A + B -> B
  Add(alpha, A, B)

  alpha.M -> M
  Mlt(alpha, M)

  A = A(row_perm, col_perm)
  ApplyPermutation(A, row_perm, col_perm)

  A(row_perm, col_perm) = A
  ApplyInversePermutation(A, row_perm, col_perm)
  
  A = Drow * A * Dcol
  ScaleMatrix(A, Drow, Dcol)
  
  A = Drow * A
  ScaleLeftMatrix(A, Drow)

  A = A * Dcol
  ScaleRightMatrix(A, Dcol)

*/

namespace Seldon
{
  
  //! B = B + alpha A
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric,
		 ArrayRowSymComplexSparse, Allocator1>& A,
		 Matrix<T2, Symmetric, ArrayRowSymComplexSparse, Allocator2>& B)
  {
    int m = B.GetM(), n, ni;
    Vector<T2> value(2*B.GetN());
    IVect index(2*B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	ni = A.GetImagRowSize(i);
        for (int j = 0; j < n; j++)
	  {
	    value(j) = alpha*T2(A.ValueReal(i, j), 0);
	    index(j) = A.IndexReal(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(j+n) = alpha*T2(0, A.ValueImag(i, j));
	    index(j+n) = A.IndexImag(i, j);
	  }
	
	B.AddInteractionRow(i, n+ni, index, value);
      }
  }
  

  //! B = B + alpha A  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha, const Matrix<T1, Symmetric,
		 ArrayRowSymComplexSparse, Allocator1>& A,
		 Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B)
  {
    int m = B.GetM(), n;
    Vector<T2, VectFull, Allocator2> value(B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	for (int j = 0; j < n; j++)
	  value(j) = alpha*T1(A.ValueReal(i, j), 0);

	B.AddInteractionRow(i, n, A.GetRealInd(i), value.GetData());
	
        n = A.GetImagRowSize(i);
        for (int j = 0; j < n; j++)
	  value(j) = alpha*T1(0, A.ValueImag(i, j));

	B.AddInteractionRow(i, n, A.GetImagInd(i), value.GetData());
      }
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void AddMatrix(const T0& alpha, const Matrix<T1, General,
		 ArrayRowComplexSparse, Allocator1>& A,
		 Matrix<T2, General, ArrayRowSparse, Allocator2>& B)
  {
    int m = B.GetM(), n;
    Vector<T2, VectFull, Allocator2> value(B.GetN());
    for (int i = 0; i < m; i++)
      {
	n = A.GetRealRowSize(i);
        for (int j = 0; j < n; j++)
          value(j) = alpha*T1(A.ValueReal(i, j), 0);
        
        B.AddInteractionRow(i, n, A.GetRealInd(i), value.GetData());
              
        n = A.GetImagRowSize(i);
        for (int j = 0; j < n; j++)
          value(j) = alpha*T1(0, A.ValueImag(i, j));
            
        B.AddInteractionRow(i, n, A.GetImagInd(i), value.GetData());
      }
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void AddMatrix(const T0& alpha, const Matrix<T1, General,
		 ArrayRowComplexSparse, Allocator1>& A,
		 Matrix<T2, General, ArrayRowComplexSparse, Allocator2>& B)
  {
    int m = B.GetM(), n, ni;
    Vector<T2> value(2*B.GetN()); IVect index(2*B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	ni = A.GetImagRowSize(i);
        for (int j = 0; j < n; j++)
          {
            value(j) = alpha*T2(A.ValueReal(i, j), 0);
            index(j) = A.IndexReal(i, j);
          }
        
        for (int j = 0; j < ni; j++)
          {
            value(n+j) = alpha*T2(0, A.ValueImag(i, j));
            index(n+j) = A.IndexImag(i, j);
          }
            
        B.AddInteractionRow(i, n+ni, index, value);
      }
  }
  
  
  //! C = C + complex(A,B)
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2,
           class T3, class Prop3, class Allocator3>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& A,
		 const Matrix<T2, Prop2, ArrayRowSymSparse, Allocator2>& B,
		 Matrix<T3, Prop3, ArrayRowSymComplexSparse, Allocator3>& C)
  {
    int m = B.GetM(), n, ni;
    Vector<T3> value(2*B.GetN()); IVect index(2*B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	ni = B.GetRowSize(i);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = alpha*T3(A.Value(i, j), 0);
	    index(j) = A.Index(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = alpha*T3(0, B.Value(i, j));
	    index(n+j) = B.Index(i, j);
	  }

	C.AddInteractionRow(i, n+ni, index, value);
      }
  }


  //! C = C + complex(A,B)
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2,
           class T3, class Prop3, class Allocator3>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		 const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		 Matrix<T3, Prop3, ArrayRowComplexSparse, Allocator3>& C)
  {
    int m = B.GetM(), n, ni;
    Vector<T3> value(2*B.GetN());
    IVect index(2*B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	ni = B.GetRowSize(i);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = alpha*T3(A.Value(i, j), 0);
	    index(j) = A.Index(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = alpha*T3(0, B.Value(i, j));
	    index(n+j) = B.Index(i, j);
	  }

	C.AddInteractionRow(i, n+ni, index, value);
      }
  }
  

  template<class T, class Allocator>
  void Add_csr_ptr(const T& alpha, long* ptr_A, int* ind_A, T* data_A,
                   long* ptr_B, int* ind_B, T* data_B, int m,
                   Vector<long>& Ptr,
                   Vector<int>& Ind,
                   Vector<T, VectFull, Allocator>& Val)
  {
    int i = 0;
    long j = 0;
    long k;
    
    // A and B might have the same structure
    // Loop over all non-zeros. If the structures of A and B differ at any
    // time, the loop is broken and a different strategy is undertaken.
    for (i = 0; i < m; i++)
      if (ptr_A[i + 1] == ptr_B[i + 1])
        {
          for (j = ptr_A[i]; j < ptr_A[i + 1]; j++)
            if (ind_A[j] == ind_B[j])
              data_B[j] += alpha * data_A[j];
            else
              break;
          if (j != ptr_A[i + 1])
            break;
        }
      else
        break;
    
    // Success: A and B have the same structure.
    if (i == m)
      return;
    
    // The addition is performed row by row in the following lines. Thus the
    // additions already performed in the current line, if started, should be
    // canceled.
    for (k = ptr_A[i]; k < j; k++)
      if (ind_A[k] == ind_B[k])
        data_B[k] -= alpha * data_A[k];

    // Number of non zero entries currently found.
    long Nnonzero = ptr_A[i];
    
    // counting the number of non-zero entries
    long kb, jb(0), ka, ja(0);
    for (int i2 = i; i2 < m; i2++)
      {
        kb = ptr_B[i2];
        
        for (ka = ptr_A[i2]; ka < ptr_A[i2 + 1]; ka++)
          {
            ja = ind_A[ka];
            while (kb < ptr_B[i2 + 1] && ind_B[kb] < ja)
              {
                kb++;
                Nnonzero++;
              }
            
            if (kb < ptr_B[i2 + 1] && ja == ind_B[kb])
              kb++;
            
            Nnonzero++;
          }

        while (kb < ptr_B[i2 + 1])
          {
            kb++;
            Nnonzero++;
          }
      }
    
    // A and B do not have the same structure. An intermediate matrix will be
    // needed. The first i rows have already been added. These computations
    // are preserved in arrays Ptr, Ind Val.
    Ptr.Reallocate(m+1); Ind.Reallocate(Nnonzero);
    Val.Reallocate(Nnonzero);
    for (int i2 = 0; i2 <= i; i2++)
      Ptr(i2) = ptr_B[i2];
    
    for (j = 0; j < ptr_B[i]; j++)
      {
        Ind(j) = ind_B[j];
        Val(j) = data_B[j];
      }

    // Now deals with the remaining lines.
    Nnonzero = ptr_A[i];
    for (; i < m; i++)
      {
        kb = ptr_B[i];
        if (kb < ptr_B[i + 1])
          jb = ind_B[kb];
        for (ka = ptr_A[i]; ka < ptr_A[i + 1]; ka++)
          {
            ja = ind_A[ka];
            while (kb < ptr_B[i + 1] && jb < ja)
              // For all elements in B that are before the ka-th element of A.
              {
                Ind(Nnonzero) = jb;
                Val(Nnonzero) = data_B[kb];
                kb++;
                if (kb < ptr_B[i + 1])
                  jb = ind_B[kb];
                Nnonzero++;
              }

            if (kb < ptr_B[i + 1] && ja == jb)
              // The element in A is also in B.
              {
                Ind(Nnonzero) = jb;
                Val(Nnonzero) = data_B[kb] + alpha * data_A[ka];
                kb++;
                if (kb < ptr_B[i + 1])
                  jb = ind_B[kb];
              }
            else
              {
                Ind(Nnonzero) = ja;
                Val(Nnonzero) = alpha * data_A[ka];
              }
            Nnonzero++;
          }

        // The remaining elements from B.
        while (kb < ptr_B[i + 1])
          {
            Ind(Nnonzero) = jb;
            Val(Nnonzero) = data_B[kb];
            kb++;
            if (kb < ptr_B[i + 1])
              jb = ind_B[kb];
            Nnonzero++;
          }

        Ptr(i + 1) = Nnonzero;
      }    
  }
  
  
  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, RowSymComplexSparse, Allocator1>& A,
		 Matrix<T2, Symmetric, RowSymComplexSparse, Allocator2>& B)
  {
    Vector<int> IndReal, IndImag;
    Vector<long> PtrReal, PtrImag;
    
    Vector<typename ClassComplexType<T2>::Treal,
	   VectFull, Allocator2> DataReal, DataImag;
    
    Add_csr_ptr(alpha, A.GetRealPtr(), A.GetRealInd(), A.GetRealData(),
                B.GetRealPtr(), B.GetRealInd(), B.GetRealData(), B.GetM(),
                PtrReal, IndReal, DataReal);

    Add_csr_ptr(alpha, A.GetImagPtr(), A.GetImagInd(), A.GetImagData(),
                B.GetImagPtr(), B.GetImagInd(), B.GetImagData(), B.GetM(),
                PtrImag, IndImag, DataImag);

    B.SetData(B.GetM(), B.GetN(), DataReal, PtrReal, IndReal,
              DataImag, PtrImag, IndImag);
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, RowComplexSparse, Allocator1>& A,
		 Matrix<T2, General, RowComplexSparse, Allocator2>& B)
  {
    Vector<int> IndReal, IndImag;
    Vector<long> PtrReal, PtrImag;
    
    Vector<typename ClassComplexType<T2>::Treal,
	   VectFull, Allocator2> DataReal, DataImag;
    
    Add_csr_ptr(alpha, A.GetRealPtr(), A.GetRealInd(), A.GetRealData(),
                B.GetRealPtr(), B.GetRealInd(), B.GetRealData(), B.GetM(),
                PtrReal, IndReal, DataReal);

    Add_csr_ptr(alpha, A.GetImagPtr(), A.GetImagInd(), A.GetImagData(),
                B.GetImagPtr(), B.GetImagInd(), B.GetImagData(), B.GetM(),
                PtrImag, IndImag, DataImag);

    B.SetData(B.GetM(), B.GetN(), DataReal, PtrReal, IndReal,
              DataImag, PtrImag, IndImag);
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void AddMatrix(const complex<T0>& alpha,
		 const Matrix<T1, Symmetric, RowSymComplexSparse, Allocator1>& A,
		 Matrix<T2, Symmetric, RowSymComplexSparse, Allocator2>& B)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Add(Matrix<RowSymComplexSparse>)",
                      "Function not implemented for complex scalars");

    Vector<int> IndReal, IndImag;
    Vector<long> PtrReal, PtrImag;
    
    Vector<typename ClassComplexType<T2>::Treal,
	   VectFull, Allocator2> DataReal, DataImag;
    
    Add_csr_ptr(real(alpha), A.GetRealPtr(), A.GetRealInd(), A.GetRealData(),
                B.GetRealPtr(), B.GetRealInd(), B.GetRealData(), B.GetM(),
                PtrReal, IndReal, DataReal);

    Add_csr_ptr(real(alpha), A.GetImagPtr(), A.GetImagInd(), A.GetImagData(),
                B.GetImagPtr(), B.GetImagInd(), B.GetImagData(), B.GetM(),
                PtrImag, IndImag, DataImag);

    B.SetData(B.GetM(), B.GetN(), DataReal, PtrReal, IndReal,
              DataImag, PtrImag, IndImag);
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void AddMatrix(const complex<T0>& alpha,
		 const Matrix<T1, General, RowComplexSparse, Allocator1>& A,
		 Matrix<T2, General, RowComplexSparse, Allocator2>& B)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Add(Matrix<RowComplexSparse>)",
                      "Function not implemented for complex scalars");

    Vector<int> IndReal, IndImag;
    Vector<long> PtrReal, PtrImag;
    
    Vector<typename ClassComplexType<T2>::Treal,
	   VectFull, Allocator2> DataReal, DataImag;
    
    Add_csr_ptr(real(alpha), A.GetRealPtr(), A.GetRealInd(), A.GetRealData(),
                B.GetRealPtr(), B.GetRealInd(), B.GetRealData(), B.GetM(),
                PtrReal, IndReal, DataReal);

    Add_csr_ptr(real(alpha), A.GetImagPtr(), A.GetImagInd(), A.GetImagData(),
                B.GetImagPtr(), B.GetImagInd(), B.GetImagData(), B.GetM(),
                PtrImag, IndImag, DataImag);

    B.SetData(B.GetM(), B.GetN(), DataReal, PtrReal, IndReal,
              DataImag, PtrImag, IndImag);
  }
  

  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void MltScalar(const T0& alpha,
		 Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i, j) *= alpha;
        
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i, j) *= alpha;
        
      }
  }


  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void MltScalar(const complex<T0>& alpha,
		 Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Mlt(Matrix<ArrayRowComplexSparse>)",
                      "Function not implemented for complex scalars");
    
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i, j) *= real(alpha);
        
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i, j) *= real(alpha);        
      }
  }


  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void MltScalar(const T0& alpha,
		 Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i,j) *= alpha;

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i,j) *= alpha;
      }
  }
  

  //! multiplication by a scalar  
  template<class T0, class T, class Prop, class Allocator>
  void MltScalar(const complex<T0>& alpha,
		 Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Mlt(Matrix<ArrayRowComplexSparse>)",
                      "Function not implemented for complex scalars");
    
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i,j) *= real(alpha);

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i,j) *= real(alpha);
      }
  }


  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void MltScalar(const T0& alpha,
		 Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal* data_A = A.GetRealData();   
    for (long i = 0; i < A.GetRealDataSize(); i++)
      data_A[i] *= alpha;

    data_A = A.GetImagData();
    for (long i = 0; i < A.GetImagDataSize(); i++)
      data_A[i] *= alpha;
  }
  

  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void MltScalar(const complex<T0>& alpha,
		 Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Mlt(Matrix<RowComplexSparse>)",
                      "Function not implemented for complex scalars");
    
    typename ClassComplexType<T>::Treal* data_A = A.GetRealData();
    T0 alpha_r = real(alpha);
    for (long i = 0; i < A.GetRealDataSize(); i++)
      data_A[i] *= alpha_r;

    data_A = A.GetImagData();
    for (long i = 0; i < A.GetImagDataSize(); i++)
      data_A[i] *= alpha_r;
  }


  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void MltScalar(const T0& alpha,
		 Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal* data_A = A.GetRealData();
    for (long i = 0; i < A.GetRealDataSize(); i++)
      data_A[i] *= alpha;

    data_A = A.GetImagData();
    for (long i = 0; i < A.GetImagDataSize(); i++)
      data_A[i] *= alpha;
  }
  

  //! multiplication by a scalar  
  template<class T0, class T, class Prop, class Allocator>
  void MltScalar(const complex<T0>& alpha,
		 Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Mlt(Matrix<RowSymComplexSparse>)",
                      "Function not implemented for complex scalars");
    
    typename ClassComplexType<T>::Treal* data_A = A.GetRealData();
    T0 alpha_r = real(alpha);
    for (long i = 0; i < A.GetRealDataSize(); i++)
      data_A[i] *= alpha_r;

    data_A = A.GetImagData();
    for (long i = 0; i < A.GetImagDataSize(); i++)
      data_A[i] *= alpha_r;
  }

  
  //! Permutation of a general matrix stored by rows.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop,
			       ArrayRowComplexSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm)
  {
    int m = A.GetM();
    IVect ind_tmp, iperm(m), rperm(m);
    for (int i = 0; i < m; i++)
      {
	iperm(i) = i;
	rperm(i) = i;
      }
    
    // A(rperm(i),:) will be the place where is the initial row i.

    // Algorithm avoiding the allocation of another matrix.
    for (int i = 0; i < m; i++)
      {
	// We get the index of row where the row initially placed on row i is.
	int i2 = rperm(i);
	// We get the new index of this row.
	int i_ = row_perm(i);

	// We fill ind_tmp of the permuted indices of columns of row i.
	int nr = A.GetRealRowSize(i2);
	ind_tmp.Reallocate(nr);
	for (int j = 0; j < nr; j++)
	  ind_tmp(j) = col_perm(A.IndexReal(i2, j));

	// We swap the two rows i and its destination row_perm(i).
	A.SwapRealRow(i2, i_);
	A.ReplaceRealIndexRow(i_, ind_tmp);

	int ni = A.GetImagRowSize(i2);
	ind_tmp.Reallocate(ni);
	for (int j = 0; j < ni; j++)
	  ind_tmp(j) = col_perm(A.IndexImag(i2, j));

	A.SwapImagRow(i2, i_);
	A.ReplaceImagIndexRow(i_, ind_tmp);
        
	// We update the indices iperm and rperm in order to keep in memory
	// the place where the row row_perm(i) is.
	int i_tmp = iperm(i_);
	iperm(i_) = iperm(i2);
	iperm(i2) = i_tmp;
	rperm(iperm(i_)) = i_;
	rperm(iperm(i2)) = i2;

	// We assemble the row i.
	A.AssembleRealRow(i_);
        A.AssembleImagRow(i_);
      }
  }

  
  //! Permutation of a symmetric matrix stored by rows.
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop,
                               ArrayRowSymComplexSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm)
  {
    // It is assumed that the permuted matrix is still symmetric! For example,
    // the user can provide row_perm = col_perm.
    int m = A.GetM();
    long nnz_real = A.GetRealDataSize(), nnz_imag = A.GetImagDataSize();
    IVect IndRow(nnz_real), IndCol(nnz_real);
    Vector<typename ClassComplexType<T>::Treal,
	   VectFull, Allocator> Val(nnz_real);

    // First we convert the matrix in coordinate format and we permute the
    // indices.
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    long k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  {
	    IndRow(k) = row_perm(i);
	    Val(k) = A.ValueReal(i,j);
	    IndCol(k) = col_perm(A.IndexReal(i, j));
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

    // We sort by row number.
    Sort(nnz_real, IndRow, IndCol, Val);

    // A is filled.
    k = 0;
    for (int i = 0; i < m; i++)
      {
	long first_index = k;
	// We get the size of the row i.
	while (k < nnz_real && IndRow(k) <= i)
	  k++;

	int size_row = k - first_index;
	// If row not empty.
	if (size_row > 0)
	  {
	    A.ReallocateRealRow(i, size_row);
	    k = first_index;
	    Sort(k, k+size_row-1, IndCol, Val);
	    for (int j = 0; j < size_row; j++)
	      {
		A.IndexReal(i,j) = IndCol(k);
		A.ValueReal(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearRealRow(i);
      }

    // Same procedure for imaginary part.

    IndRow.Reallocate(nnz_imag);
    IndCol.Reallocate(nnz_imag);
    Val.Reallocate(nnz_imag);

    // First we convert the matrix in coordinate format and we permute the
    // indices.
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  {
	    IndRow(k) = row_perm(i);
	    Val(k) = A.ValueImag(i,j);
	    IndCol(k) = col_perm(A.IndexImag(i,j));
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
    // We sort by row number.
    Sort(nnz_imag, IndRow, IndCol, Val);

    // A is filled
    k = 0;
    for (int i = 0; i < m; i++)
      {
	long first_index = k;
	// We get the size of the row i.
	while (k < nnz_imag && IndRow(k) <= i)
	  k++;
	int size_row = k - first_index;
	// If row not empty.
	if (size_row > 0)
	  {
	    A.ReallocateImagRow(i, size_row);
	    k = first_index;
	    Sort(k, k+size_row-1, IndCol, Val);
	    for (int j = 0; j < size_row; j++)
	      {
		A.IndexImag(i,j) = IndCol(k);
		A.ValueImag(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearImagRow(i);
      }
  }

  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, col_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
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
  template<class T, class Prop, class Allocator> void
  ApplyPermutation(Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A,
		   const Vector<int>& row_perm,
		   const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    ApplyInversePermutation(A, inv_row_perm, inv_row_perm);
  }
  
  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop,
		   ArrayRowSymComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale_left(i) * scale_right(A.IndexReal(i, j));

	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale_left(i) * scale_right(A.IndexImag(i, j));
      }
  }


  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale_left(i) * scale_right(A.IndexReal(i, j));

	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale_left(i) * scale_right(A.IndexImag(i, j));
      }
  }

  
  //! Each row is scaled.
  /*!
    We compute diag(S)*A where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ArrayRowComplexSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale(i);

	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale(i);
      }
  }
  
  
  //! Each column is scaled.
  /*!
    We compute A*diag(S) where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop,
			ArrayRowComplexSparse, Allocator1>& A,
			const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale(A.IndexReal(i, j));

	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale(A.IndexImag(i, j));
      }
  }

  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, RowSymComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T1>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T1>::Treal* data_imag = A.GetImagData();
    for (int i = 0; i < m; i++ )
      {
	for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
	  data_real[j] *= scale_left(i) * scale_right(ind_real[j]);

	for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
	  data_imag[j] *= scale_left(i) * scale_right(ind_imag[j]);
      }
  }


  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, RowComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T1>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T1>::Treal* data_imag = A.GetImagData();
    for (int i = 0; i < m; i++ )
      {
	for (long j = ptr_real[i]; j < ptr_real[i+1]; j++ )
	  data_real[j] *= scale_left(i) * scale_right(ind_real[j]);

	for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++ )
	  data_imag[j] *= scale_left(i) * scale_right(ind_imag[j]);
      }
  }

  
  //! Each row is scaled.
  /*!
    We compute diag(S)*A where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, RowComplexSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    long* ptr_real = A.GetRealPtr();
    typename ClassComplexType<T1>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    typename ClassComplexType<T1>::Treal* data_imag = A.GetImagData();
    for (int i = 0; i < m; i++ )
      {
	for (long j = ptr_real[i]; j < ptr_real[i+1]; j++ )
	  data_real[j] *= scale(i);

	for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++ )
	  data_imag[j] *= scale(i);
      }
  }
  
  
  //! Each column is scaled.
  /*!
    We compute A*diag(S) where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, RowComplexSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T1>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T1>::Treal* data_imag = A.GetImagData();
    for (int i = 0; i < m; i++ )
      {
	for (long j = ptr_real[i]; j < ptr_real[i+1]; j++ )
	  data_real[j] *= scale(ind_real[j]);

	for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++ )
	  data_imag[j] *= scale(ind_imag[j]);
      }
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General,
		 ArrayRowComplexSparse, Allocator>& A,
                 Matrix<T, General, ArrayRowComplexSparse, Allocator>& B)
  {
    B.Clear();
    
    int m = A.GetM();
    int n = A.GetN();
    Vector<int> ptr_r(n), ptr_i(n);
    
    B.Reallocate(n, m);

    // For each column j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    ptr_r.Zero();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRealRowSize(i); j++)
        ptr_r(A.IndexReal(i, j))++;

    ptr_i.Zero();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetImagRowSize(i); j++)
        ptr_i(A.IndexImag(i, j))++;
    
    for (int i = 0; i < n; i++)
      {
        B.ReallocateRealRow(i, ptr_r(i));
        B.ReallocateImagRow(i, ptr_i(i));
      }
    
    // filling matrix B
    ptr_r.Zero();
    ptr_i.Zero();
    for (int i = 0; i < m; i++)
      {
        for (int jp = 0; jp < A.GetRealRowSize(i); jp++)
          {
            int j = A.IndexReal(i, jp);
            int k = ptr_r(j);
            ++ptr_r(j);
            B.ValueReal(j, k) = A.ValueReal(i, jp);
            B.IndexReal(j, k) = i;
          }

        for (int jp = 0; jp < A.GetImagRowSize(i); jp++)
          {
            int j = A.IndexImag(i, jp);
            int k = ptr_i(j);
            ++ptr_i(j);
            B.ValueImag(j, k) = A.ValueImag(i, jp);
            B.IndexImag(j, k) = i;
          }
      }
    
    // sorting numbers 
    B.Assemble();
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, RowComplexSparse, Allocator>& A,
                 Matrix<T, General, RowComplexSparse, Allocator>& B)
  {
    B.Clear();
    
    int m = A.GetM();
    int n = A.GetN();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = A.GetImagData();
    Vector<int> ptr_r(n), ptr_i(n);
    
    // For each column j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    ptr_r.Zero();
    for (int i = 0; i < m; i++)
      for (int j = ptr_real[i]; j < ptr_real[i+1]; j++)
        ptr_r(ind_real[j])++;

    ptr_i.Zero();
    for (int i = 0; i < m; i++)
      for (int j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
        ptr_i(ind_imag[j])++;
    
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal, VectFull, Allocator> ValReal(A.GetRealDataSize());
    Vector<Treal, VectFull, Allocator> ValImag(A.GetImagDataSize());
    Vector<long> PtrReal(n+1), PtrImag(n+1);
    Vector<int> IndReal(A.GetRealDataSize()), IndImag(A.GetImagDataSize());
    
    PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < n; i++)
      {
        PtrReal(i+1) = PtrReal(i) + ptr_r(i);
        PtrImag(i+1) = PtrImag(i) + ptr_i(i);
      }

    // filling matrix B
    ptr_r.Zero();
    ptr_i.Zero();
    for (int i = 0; i < m; i++)
      {
        for (long jp = ptr_real[i]; jp < ptr_real[i+1]; jp++)
          {
            int j = ind_real[jp];
            int k = ptr_r(j);
            ValReal(PtrReal(j) + k) = data_real[jp];
            IndReal(PtrReal(j) + k) = i;
            ++ptr_r(j);
          }

        for (long jp = ptr_imag[i]; jp < ptr_imag[i+1]; jp++)
          {
            int j = ind_imag[jp];
            int k = ptr_i(j);
            ValImag(PtrImag(j) + k) = data_imag[jp];
            IndImag(PtrImag(j) + k) = i;
            ++ptr_i(j);
          }
      }
    
    // sorting numbers
    for (int i = 0; i < n; i++)
      {
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
      }
    
    B.SetData(n, m, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                res = max(res, abs(A.ValueImag(i, ji)));
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                res = max(res, ComplexAbs(T(A.ValueReal(i, j), 
					    A.ValueImag(i, ji))));
                ji++;
              }
            else
              res = max(res, abs(A.ValueReal(i, j)));
          }
        
        while (ji < size_imag)
          {
            res = max(res, abs(A.ValueImag(i, ji)));
            ji++;
          }
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
  MaxAbs(const Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = A.GetImagData();
    for (int i = 0; i < A.GetM(); i++)
      {
        long ji = ptr_imag[i];
        for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                res = max(res, abs(data_imag[ji]));
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                res = max(res, ComplexAbs(T(data_real[j], 
					    data_imag[ji])));
                ji++;
              }
            else
              res = max(res, abs(data_real[j]));
          }
        
        while (ji < ptr_imag[i+1])
          {
            res = max(res, abs(data_imag[ji]));
            ji++;
          }
      }
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    Vector<typename ClassComplexType<T>::Treal> sum(A.GetN());
    sum.Fill(0);
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                sum(A.IndexImag(i, ji)) += abs(A.ValueImag(i, ji));
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                sum(k) += ComplexAbs(T(A.ValueReal(i, j), 
				       A.ValueImag(i, ji)));
                ji++;
              }
            else
              sum(k) += abs(A.ValueReal(i, j));
          }
        
        while (ji < size_imag)
          {
            sum(A.IndexImag(i, ji)) += abs(A.ValueImag(i, ji));
            ji++;
          }
      }
    
    return sum.GetNormInf();
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    Vector<typename ClassComplexType<T>::Treal> sum(A.GetN());
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = A.GetImagData();
    sum.Fill(0);
    for (int i = 0; i < A.GetM(); i++)
      {
        long ji = ptr_imag[i];
        for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                sum(ind_imag[ji]) += abs(data_imag[ji]);
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                sum(k) += ComplexAbs(T(data_real[j], 
				       data_imag[ji]));
                ji++;
              }
            else
              sum(k) += abs(data_real[j]);
          }
        
        while (ji < ptr_imag[i+1])
          {
            sum(ind_imag[ji]) += abs(data_imag[ji]);
            ji++;
          }
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
  NormInf(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0), sum;
    for (int i = 0; i < A.GetM(); i++)
      {
        sum = 0;
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                sum += abs(A.ValueImag(i, ji));
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                sum += ComplexAbs(T(A.ValueReal(i, j), 
				    A.ValueImag(i, ji)));
                ji++;
              }
            else
              sum += abs(A.ValueReal(i, j));
          }
        
        while (ji < size_imag)
          {
            sum += abs(A.ValueImag(i, ji));
            ji++;
          }
        
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
  NormInf(const Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0), sum;
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = A.GetImagData();
    for (int i = 0; i < A.GetM(); i++)
      {
        sum = 0;
        long ji = ptr_imag[i];
        for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                sum += abs(data_imag[ji]);
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                sum += ComplexAbs(T(data_real[j], 
				    data_imag[ji]));
                ji++;
              }
            else
              sum += abs(data_real[j]);
          }
        
        while (ji < ptr_imag[i+1])
          {
            sum += abs(data_imag[ji]);
            ji++;
          }
        
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
  MaxAbs(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                res = max(res, abs(A.ValueImag(i, ji)));
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                res = max(res, ComplexAbs(T(A.ValueReal(i, j), 
					    A.ValueImag(i, ji))));
                ji++;
              }
            else
              res = max(res, abs(A.ValueReal(i, j)));
          }
        
        while (ji < size_imag)
          {
            res = max(res, abs(A.ValueImag(i, ji)));
            ji++;
          }
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
  MaxAbs(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = A.GetImagData();
    for (int i = 0; i < A.GetM(); i++)
      {
        long ji = ptr_imag[i];
        for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                res = max(res, abs(data_imag[ji]));
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                res = max(res, ComplexAbs(T(data_real[j], 
					    data_imag[ji])));
                ji++;
              }
            else
              res = max(res, abs(data_real[j]));
          }
        
        while (ji < ptr_imag[i+1])
          {
            res = max(res, abs(data_imag[ji]));
            ji++;
          }
      }
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal val;
    Vector<typename ClassComplexType<T>::Treal> sum(A.GetN());
    sum.Fill(0);
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                val = abs(A.ValueImag(i, ji));
                sum(A.IndexImag(i, ji)) += val;
                if (A.IndexImag(i, ji) != i)
                  sum(i) += val;
                
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                val = ComplexAbs(T(A.ValueReal(i, j), 
				   A.ValueImag(i, ji)));
                sum(k) += val;
                if (k != i)
                  sum(i) += val;
                
                ji++;
              }
            else
              {
                val = abs(A.ValueReal(i, j));
                sum(k) += val;
                if (k != i)
                  sum(i) += val;
              }
          }
        
        while (ji < size_imag)
          {
            val = abs(A.ValueImag(i, ji));
            sum(A.IndexImag(i, ji)) += val;
            if (A.IndexImag(i, ji) != i)
              sum(i) += val;
            
            ji++;
          }
      }
    
    return sum.GetNormInf();
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    Vector<typename ClassComplexType<T>::Treal> sum(A.GetN());
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = A.GetImagData();
    sum.Fill(0);
    typename ClassComplexType<T>::Treal val;
    for (int i = 0; i < A.GetM(); i++)
      {
        long ji = ptr_imag[i];
        for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                val = abs(data_imag[ji]);
                sum(ind_imag[ji]) += val;
                if (ind_imag[ji] != i)
                  sum(i) += val;
                
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                val = ComplexAbs(T(data_real[j], 
				   data_imag[ji]));
                sum(k) += val;
                if (k != i)
                  sum(i) += val;
                
                ji++;
              }
            else
              {
                val = abs(data_real[j]);
                sum(k) += val;
                if (k != i)
                  sum(i) += val;
              }
          }
        
        while (ji < ptr_imag[i+1])
          {
            val = abs(data_imag[ji]);
            sum(ind_imag[ji]) += val;
            if (ind_imag[ji] != i)
              sum(i) += val;
            
            ji++;
          }
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
  NormInf(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    return Norm1(A);
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    return Norm1(A);
  }
  
  
  //! A is replaced by its conjugate
  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetImagRowSize(i); j++)
        A.ValueImag(i, j) = -A.ValueImag(i, j);
  }


  //! A is replaced by its conjugate
  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetImagRowSize(i); j++)
        A.ValueImag(i, j) = -A.ValueImag(i, j); 
  }


  //! A is replaced by its conjugate
  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal* data = A.GetImagData();
    for (long i = 0; i < A.GetImagDataSize(); i++)
      data[i] = -data[i];
  }


  //! A is replaced by its conjugate
  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal* data = A.GetImagData();
    for (long i = 0; i < A.GetImagDataSize(); i++)
      data[i] = -data[i];
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, Symmetric,
                 ArrayRowSymComplexSparse, Allocator>& A,
                 Matrix<T, Symmetric,
                 ArrayRowSymComplexSparse, Allocator>& B)
  {
    B = A;
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayRowComplexSparse, Allocator>& A)
  {
    Matrix<T, General, ArrayRowComplexSparse, Allocator> B(A);
    Transpose(B, A);
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, Symmetric, ArrayRowSymComplexSparse, Allocator>& A)
  {
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, Symmetric,
                 RowSymComplexSparse, Allocator>& A,
                 Matrix<T, Symmetric,
                 RowSymComplexSparse, Allocator>& B)
  {
    B = A;
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, RowComplexSparse, Allocator>& A)
  {
    Matrix<T, General, RowComplexSparse, Allocator> B(A);
    Transpose(B, A);
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, Symmetric, RowSymComplexSparse, Allocator>& A)
  {
  }
  
  
  //! drops non-zero entries below epsilon
  template<class T, class Prop, class Allocator, class T0>
  void RemoveSmallEntry(Matrix<T, Prop, RowComplexSparse, Allocator>& A,
                        const T0& epsilon)
  {
    // TO BE DONE
  }


  //! drops non-zero entries below epsilon
  template<class T, class Prop, class Allocator, class T0>
  void RemoveSmallEntry(Matrix<T, Prop, RowSymComplexSparse, Allocator>& A,
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
		Matrix<T1, Prop, ArrayRowSymComplexSparse, Allocator>& A)
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
	  {
	    A.ClearRealRow(i);
	    A.ClearImagRow(i);
	  }
      }
	
    
    // then columns
    for (int i = 0; i < A.GetM(); i++)
      {
	bool something_to_remove = false;
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  if (index(A.IndexReal(i,j)) != -1)
	    something_to_remove = true;
	
	if (something_to_remove)
	  {
	    int nb = 0;
	    for (int j = 0; j < A.GetRealRowSize(i); j++)
	      if (index(A.IndexReal(i,j)) == -1)
		{
		  A.IndexReal(i, nb) = A.IndexReal(i, j);
		  A.ValueReal(i, nb) = A.ValueReal(i, j);	      
		  nb++;
		}
	    
	    A.ResizeRealRow(i, nb);
	  }

	something_to_remove = false;
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  if (index(A.IndexImag(i,j)) != -1)
	    something_to_remove = true;
	
	if (something_to_remove)
	  {
	    int nb = 0;
	    for (int j = 0; j < A.GetImagRowSize(i); j++)
	      if (index(A.IndexImag(i,j)) == -1)
		{
		  A.IndexImag(i, nb) = A.IndexImag(i, j);
		  A.ValueImag(i, nb) = A.ValueImag(i, j);	      
		  nb++;
		}
	    
	    A.ResizeImagRow(i, nb);
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
		Matrix<T1, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    int m = col_number.GetM();
    // index array to know fastly if it is a column to erase
    IVect index(A.GetM()); index.Fill(-1);
    for (int i = 0; i < m; i++)
      index(col_number(i)) = i;
    
    for (int i = 0; i < A.GetM(); i++)
      {
	bool something_to_remove = false;
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  if (index(A.IndexReal(i,j)) != -1)
	    something_to_remove = true;
	
	if (something_to_remove)
	  {
	    int nb = 0;
	    for (int j = 0; j < A.GetRealRowSize(i); j++)
	      if (index(A.IndexReal(i,j)) == -1)
		{
		  A.IndexReal(i, nb) = A.IndexReal(i, j);
		  A.ValueReal(i, nb) = A.ValueReal(i, j);	      
		  nb++;
		}
	    
	    A.ResizeRealRow(i, nb);
	  }

	something_to_remove = false;
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  if (index(A.IndexImag(i,j)) != -1)
	    something_to_remove = true;
	
	if (something_to_remove)
	  {
	    int nb = 0;
	    for (int j = 0; j < A.GetImagRowSize(i); j++)
	      if (index(A.IndexImag(i,j)) == -1)
		{
		  A.IndexImag(i, nb) = A.IndexImag(i, j);
		  A.ValueImag(i, nb) = A.ValueImag(i, j);	      
		  nb++;
		}
	    
	    A.ResizeImagRow(i, nb);
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
		Matrix<T1, Prop, RowComplexSparse, Allocator>& A)
  {
    int m = A.GetM(), n = A.GetN();
    long nnz_real = A.GetRealIndSize();
    long nnz_imag = A.GetImagIndSize();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T1>::Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T1>::Treal* data_imag = A.GetImagData();
    Vector<bool> ColToKeep(n);
    ColToKeep.Fill(true);
    for (int i = 0; i < col_number.GetM(); i++)
      ColToKeep(col_number(i)) = false;
    
    for (long i = 0; i < A.GetRealIndSize(); i++)
      if (!ColToKeep(ind_real[i]))
        nnz_real--;

    for (long i = 0; i < A.GetImagIndSize(); i++)
      if (!ColToKeep(ind_imag[i]))
        nnz_imag--;
    
    if ((nnz_real == A.GetRealIndSize()) && (nnz_imag == A.GetImagIndSize()))
      return;
    
    Vector<long> PtrReal(m+1), PtrImag(m+1);
    Vector<int> IndReal(m+1), IndImag(nnz_imag);
    Vector<typename ClassComplexType<T1>::Treal,
	   VectFull, Allocator> ValReal(nnz_real), ValImag(nnz_imag);
    
    PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < m; i++)
      {
        long jA = PtrReal(i); int size_row = 0;
        for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
          if (ColToKeep(ind_real[j]))
            {
              IndReal(jA) = ind_real[j];
              ValReal(jA) = data_real[j];
              size_row++; jA++;
            }
        
        PtrReal(i+1) = PtrReal(i) + size_row;

        jA = PtrImag(i); size_row = 0;
        for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
          if (ColToKeep(ind_imag[j]))
            {
              IndImag(jA) = ind_imag[j];
              ValImag(jA) = data_imag[j];
              size_row++; jA++;
            }
        
        PtrImag(i+1) = PtrImag(i) + size_row;
      }
    
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }


  //! clears several columns of a sparse matrix
  /*!
    \param[in] col_number numbers of the columns to be cleared
    \param[inout] A sparse matrix where columns are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, RowSymComplexSparse, Allocator>& A)
  {
    int m = A.GetM(), n = A.GetN();
    long nnz_real = A.GetRealIndSize();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T1>::Treal* data_real = A.GetRealData();
    long nnz_imag = A.GetImagIndSize();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T1>::Treal* data_imag = A.GetImagData();
    Vector<bool> ColToKeep(n);
    ColToKeep.Fill(true);
    for (int i = 0; i < col_number.GetM(); i++)
      ColToKeep(col_number(i)) = false;
    
    for (int i = 0; i < m; i++)
      {
        if (!ColToKeep(i))
          {
            nnz_real -= ptr_real[i+1] - ptr_real[i];
            nnz_imag -= ptr_imag[i+1] - ptr_imag[i];
          }
        else
          {
            for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
              if (!ColToKeep(ind_real[j]))
                nnz_real--;
            
            for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
              if (!ColToKeep(ind_imag[j]))
                nnz_imag--;
          }
      }
    
    if ((nnz_real == A.GetRealIndSize()) && (nnz_imag == A.GetImagIndSize()))
      return;
    
    Vector<long> PtrReal(m+1), PtrImag(m+1);
    Vector<int> IndReal(nnz_real), IndImag(nnz_imag);
    Vector<typename ClassComplexType<T1>::Treal,
	   VectFull, Allocator> ValReal(nnz_real), ValImag(nnz_imag);
    
    PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < m; i++)
      {
        long jA = PtrReal(i), size_row = 0;
        if (ColToKeep(i))            
          for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
            if (ColToKeep(ind_real[j]))
              {
                IndReal(jA) = ind_real[j];
                ValReal(jA) = data_real[j];
                size_row++; jA++;
              }
        
        PtrReal(i+1) = PtrReal(i) + size_row;

        jA = PtrImag(i); size_row = 0;
        if (ColToKeep(i))            
          for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
            if (ColToKeep(ind_imag[j]))
              {
                IndImag(jA) = ind_imag[j];
                ValImag(jA) = data_imag[j];
                size_row++; jA++;
              }
        
        PtrImag(i+1) = PtrImag(i) + size_row;
      }
    
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }


  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSymComplexSparse, Allocator>& A)
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
		Matrix<T1, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < col_number.GetM(); i++)
      {
	A.ClearRealRow(col_number(i));
	A.ClearImagRow(col_number(i));
      }
  }
  

  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, RowComplexSparse, Allocator>& A)
  {
    int m = A.GetM(), n = A.GetN();
    long nnz_real = A.GetRealIndSize();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T1>::Treal* data_real = A.GetRealData();
    long nnz_imag = A.GetImagIndSize();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T1>::Treal* data_imag = A.GetImagData();
    Vector<bool> RowToKeep(m);
    RowToKeep.Fill(true);
    for (int i = 0; i < col_number.GetM(); i++)
      RowToKeep(col_number(i)) = false;
    
    for (int i = 0; i < m; i++)
      if (!RowToKeep(i))
        {
          nnz_real -= ptr_real[i+1] - ptr_real[i];
          nnz_imag -= ptr_imag[i+1] - ptr_imag[i];
        }
    
    Vector<long> PtrReal(m+1), PtrImag(m+1);
    Vector<int> IndReal(nnz_real), IndImag(nnz_imag);
    Vector<typename ClassComplexType<T1>::Treal,
	   VectFull, Allocator> ValReal(nnz_real), ValImag(nnz_imag);
    
    PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < m; i++)
      {
        if (RowToKeep(i))
          {
            int size_row = ptr_real[i+1] - ptr_real[i];
            for (int j = 0; j < size_row; j++)
              {
                IndReal(PtrReal(i) + j) = ind_real[ptr_real[i] + j];
                ValReal(PtrReal(i) + j) = data_real[ptr_real[i] + j];
              }
            
            PtrReal(i+1) = PtrReal(i) + size_row;

            size_row = ptr_imag[i+1] - ptr_imag[i];
            for (int j = 0; j < size_row; j++)
              {
                IndImag(PtrImag(i) + j) = ind_imag[ptr_imag[i] + j];
                ValImag(PtrImag(i) + j) = data_imag[ptr_imag[i] + j];
              }
            
            PtrImag(i+1) = PtrImag(i) + size_row;
          }
        else
          {
            PtrReal(i+1) = PtrReal(i);
            PtrImag(i+1) = PtrImag(i);
          }
      }
    
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }
  
  
  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, RowSymComplexSparse, Allocator>& A)
  {
    EraseCol(col_number, A);
  }

  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowSum(Vector<typename ClassComplexType<T>::Treal>& diagonal_scale_left,
		 const Matrix<T, Symmetric, ArrayRowSymComplexSparse> & mat)
  {
    int n = mat.GetM();
    diagonal_scale_left.Reallocate(n);
    diagonal_scale_left.Fill(0);
    for (int i = 0; i < n; i++)
      {
	for (int j = 0; j < mat.GetRealRowSize(i); j++)
	  {
	    diagonal_scale_left(i) += abs(mat.ValueReal(i,j));
	    if (i != mat.IndexReal(i,j))
	      diagonal_scale_left(mat.IndexReal(i,j))
		+= abs(mat.ValueReal(i,j));
	  }
	
	for (int j = 0; j < mat.GetImagRowSize(i); j++)
	  {
	    diagonal_scale_left(i) += abs(mat.ValueImag(i,j));
	    if (i != mat.IndexImag(i,j))
	      diagonal_scale_left(mat.IndexImag(i,j))
		+= abs(mat.ValueImag(i,j));
	  }
      }
  }
  
  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowSum(Vector<typename ClassComplexType<T>::Treal>& diagonal_scale_left,
		 const Matrix<T, General, ArrayRowComplexSparse>& mat)
  {
    int n = mat.GetM();
    diagonal_scale_left.Reallocate(n);
    diagonal_scale_left.Fill(0);
    for (int i = 0; i < n; i++)
      {
	for (int j = 0; j < mat.GetRealRowSize(i); j++)
	  diagonal_scale_left(i) += abs(mat.ValueReal(i,j));
	
	for (int j = 0; j < mat.GetImagRowSize(i); j++)
	  diagonal_scale_left(i) += abs(mat.ValueImag(i,j));
      }
  }
  

  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowSum(Vector<typename ClassComplexType<T>::Treal>& diagonal_scale_left,
		 const Matrix<T, Symmetric, RowSymComplexSparse> & mat)
  {
    int n = mat.GetM();
    diagonal_scale_left.Reallocate(n);
    diagonal_scale_left.Fill(0);
    long* ptr_real = mat.GetRealPtr();
    int* ind_real = mat.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = mat.GetRealData();
    long* ptr_imag = mat.GetImagPtr();
    int* ind_imag = mat.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = mat.GetImagData();
    for (int i = 0; i < n; i++)
      {
	for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
	  {
	    diagonal_scale_left(i) += abs(data_real[j]);
	    if (i != ind_real[j])
	      diagonal_scale_left(ind_real[j]) += abs(data_real[j]);
	  }
	
	for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
	  {
	    diagonal_scale_left(i) += abs(data_imag[j]);
	    if (i != ind_imag[j])
	      diagonal_scale_left(ind_imag[j]) += abs(data_imag[j]);
	  }
      }
  }
  
  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowSum(Vector<typename ClassComplexType<T>::Treal>& diagonal_scale_left,
		 const Matrix<T, General, RowComplexSparse>& mat)
  {
    int n = mat.GetM();
    diagonal_scale_left.Reallocate(n);
    diagonal_scale_left.Fill(0);
    long* ptr_real = mat.GetRealPtr();
    typename ClassComplexType<T>::Treal* data_real = mat.GetRealData();
    long* ptr_imag = mat.GetImagPtr();
    typename ClassComplexType<T>::Treal* data_imag = mat.GetImagData();
    for (int i = 0; i < n; i++)
      {
	for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
	  diagonal_scale_left(i) += abs(data_real[j]);
	
	for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
	  diagonal_scale_left(i) += abs(data_imag[j]);
      }
  }


  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetColSum(Vector<typename ClassComplexType<T>::Treal>& diagonal_scale,
		 const Matrix<T, Symmetric, ArrayRowSymComplexSparse> & mat)
  {
    GetRowSum(diagonal_scale, mat);
  }
  
  
  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetColSum(Vector<typename ClassComplexType<T>::Treal>& diagonal_scale,
		 const Matrix<T, General, ArrayRowComplexSparse>& mat)
  {
    int n = mat.GetM();
    diagonal_scale.Reallocate(mat.GetN());
    diagonal_scale.Fill(0);
    for (int i = 0; i < n; i++)
      {
	for (int j = 0; j < mat.GetRealRowSize(i); j++)
	  diagonal_scale(mat.IndexReal(i, j)) += abs(mat.ValueReal(i,j));
	
	for (int j = 0; j < mat.GetImagRowSize(i); j++)
	  diagonal_scale(mat.IndexImag(i, j)) += abs(mat.ValueImag(i,j));
      }
  }
  

  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetColSum(Vector<typename ClassComplexType<T>::Treal>& diagonal_scale,
		 const Matrix<T, Symmetric, RowSymComplexSparse> & mat)
  {
    GetRowSum(diagonal_scale, mat);
  }
  
  
  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetColSum(Vector<typename ClassComplexType<T>::Treal>& diagonal_scale,
		 const Matrix<T, General, RowComplexSparse>& mat)
  {
    int n = mat.GetM();
    diagonal_scale.Reallocate(mat.GetN());
    diagonal_scale.Fill(0);
    long* ptr_real = mat.GetRealPtr();
    int* ind_real = mat.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = mat.GetRealData();
    long* ptr_imag = mat.GetImagPtr();
    int* ind_imag = mat.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = mat.GetImagData();
    for (int i = 0; i < n; i++)
      {
	for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
	  diagonal_scale(ind_real[j]) += abs(data_real[j]);
	
	for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
	  diagonal_scale(ind_imag[j]) += abs(data_imag[j]);
      }
  }


  //! For each row and column of the matrix,
  //! computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowColSum(Vector<typename ClassComplexType<T>::Treal>& sum_row,
		    Vector<typename ClassComplexType<T>::Treal>& sum_col,
		    const Matrix<T, Symmetric, ArrayRowSymComplexSparse> & mat)
  {
    GetRowSum(sum_row, mat);
    sum_col = sum_row;
  }
  

  //! For each row and column of the matrix,
  //! computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowColSum(Vector<typename ClassComplexType<T>::Treal>& sum_row,
		    Vector<typename ClassComplexType<T>::Treal>& sum_col,
		    const Matrix<T, General, ArrayRowComplexSparse>& mat)
  {
    int n = mat.GetM();
    sum_row.Reallocate(n);
    sum_col.Reallocate(mat.GetN());
    sum_row.Fill(0);
    sum_col.Fill(0);
    for (int i = 0; i < n; i++)
      {
	for (int j = 0; j < mat.GetRealRowSize(i); j++)
	  {
	    sum_row(i) += abs(mat.ValueReal(i,j));
	    sum_col(mat.IndexReal(i, j)) += abs(mat.ValueReal(i,j));
	  }
	
	for (int j = 0; j < mat.GetImagRowSize(i); j++)
	  {
	    sum_row(i) += abs(mat.ValueImag(i,j));
	    sum_col(mat.IndexImag(i, j)) += abs(mat.ValueImag(i,j));
	  }
      }
  }


  //! For each row and column of the matrix,
  //! computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowColSum(Vector<typename ClassComplexType<T>::Treal>& sum_row,
		    Vector<typename ClassComplexType<T>::Treal>& sum_col,
		    const Matrix<T, Symmetric, RowSymComplexSparse> & mat)
  {
    GetRowSum(sum_row, mat);
    sum_col = sum_row;
  }
  

  //! For each row and column of the matrix,
  //! computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowColSum(Vector<typename ClassComplexType<T>::Treal>& sum_row,
		    Vector<typename ClassComplexType<T>::Treal>& sum_col,
		    const Matrix<T, General, RowComplexSparse>& mat)
  {
    int n = mat.GetM();
    sum_row.Reallocate(n);
    sum_col.Reallocate(mat.GetN());
    sum_row.Fill(0);
    sum_col.Fill(0);
    long* ptr_real = mat.GetRealPtr();
    int* ind_real = mat.GetRealInd();
    typename ClassComplexType<T>::Treal* data_real = mat.GetRealData();
    long* ptr_imag = mat.GetImagPtr();
    int* ind_imag = mat.GetImagInd();
    typename ClassComplexType<T>::Treal* data_imag = mat.GetImagData();
    for (int i = 0; i < n; i++)
      {
	for (long j = ptr_real[i]; j < ptr_real[i+1]; j++)
	  {
	    sum_row(i) += abs(data_real[j]);
	    sum_col(ind_real[j]) += abs(data_real[j]);
	  }
	
	for (long j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
	  {
	    sum_row(i) += abs(data_imag[j]);
	    sum_col(ind_imag[j]) += abs(data_imag[j]);
	  }
      }
  }


  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0,
		     ArrayRowComplexSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<complex<T1>, VectFull, Allocator1>& Value)
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
    
    // counting the number of non-zero elements to keep
    long nnz = 0;
    for (int i = 0; i < A.GetM(); i++)
      if (RowKept(i))
        {
          int ji = 0;
          int size_imag = A.GetImagRowSize(i);        
          for (int jr = 0; jr < A.GetRealRowSize(i); jr++)
            {
              int jcol = A.IndexReal(i, jr);
              while ((ji < size_imag) && (A.IndexImag(i, ji) < jcol))
                {
                  if (ColKept(A.IndexImag(i, ji)))
                    nnz++;
                  
                  ji++;
                }
              
              if ((ji < size_imag) && (A.IndexImag(i, ji) == jcol))
                ji++;
              
              if (ColKept(jcol))
                nnz++;
            }
          
          while (ji < size_imag)
            {
              if (ColKept(A.IndexImag(i, ji)))
                nnz++;
              
              ji++;
            }
        }
    
    RowNum.Reallocate(nnz);
    ColNum.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    // then filling the arrays RowNum, ColNum, Value
    for (int i = 0; i < A.GetM(); i++)
      if (RowKept(i))
        {
          int ji = 0;
          int size_imag = A.GetImagRowSize(i);        
          for (int jr = 0; jr < A.GetRealRowSize(i); jr++)
            {
              int jcol = A.IndexReal(i, jr);
              while ((ji < size_imag) && (A.IndexImag(i, ji) < jcol))
                {
                  if (ColKept(A.IndexImag(i, ji)))
                    {
                      RowNum(nnz) = i;
                      ColNum(nnz) = A.IndexImag(i, ji);
                      Value(nnz) = T0(0, A.ValueImag(i, ji));
                      nnz++;
                    }
                  
                  ji++;
                }
              
              if ((ji < size_imag) && (A.IndexImag(i, ji) == jcol))
                {                  
                  if (ColKept(jcol))
                    {
                      RowNum(nnz) = i;
                      ColNum(nnz) = jcol;
                      Value(nnz)
			= T0(A.ValueReal(i, jr), A.ValueImag(i, ji));
                      nnz++;
                    }
                  
                  ji++;
                }
              else
                {
                  if (ColKept(jcol))
                    {
                      RowNum(nnz) = i;
                      ColNum(nnz) = jcol;
                      Value(nnz) = T0(A.ValueReal(i, jr), 0);
                      nnz++;
                    }
                }
            }
          
          while (ji < size_imag)
            {
              if (ColKept(A.IndexImag(i, ji)))
                {
                  RowNum(nnz) = i;
                  ColNum(nnz) = A.IndexImag(i, ji);
                  Value(nnz) = T0(0, A.ValueImag(i, ji));
                  nnz++;
                }
              
              ji++;
            }
        }
  }
  

  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0,
		     ArrayRowSymComplexSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<complex<T1>, VectFull, Allocator1>& Value)
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
    
    // counting the number of non-zero elements to keep
    long nnz = 0;
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int jr = 0; jr < A.GetRealRowSize(i); jr++)
          {
            int jcol = A.IndexReal(i, jr);
            while ((ji < size_imag) && (A.IndexImag(i, ji) < jcol))
              {                              
                if (ColKept(A.IndexImag(i, ji)) && RowKept(i))
                  nnz++;
                
                if (A.IndexImag(i, ji) != i)
                  if (RowKept(A.IndexImag(i, ji)) && ColKept(i))
                    nnz++;
                
                ji++;
              }
            
            if ((ji < size_imag) && (A.IndexImag(i, ji) == jcol))
              ji++;
            
            if (ColKept(jcol) && RowKept(i))
              nnz++;
            
            if (jcol != i)
              if (RowKept(jcol) && ColKept(i))
                nnz++;
          }
        
        while (ji < size_imag)
          {
            if (ColKept(A.IndexImag(i, ji)) && RowKept(i))
              nnz++;
            
            if (A.IndexImag(i, ji) != i)
              if (RowKept(A.IndexImag(i, ji)) && ColKept(i))
                nnz++;
            
            ji++;
          }
      }
    
    RowNum.Reallocate(nnz);
    ColNum.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int jr = 0; jr < A.GetRealRowSize(i); jr++)
          {
            int jcol = A.IndexReal(i, jr);
            while ((ji < size_imag) && (A.IndexImag(i, ji) < jcol))
              {                              
                if (ColKept(A.IndexImag(i, ji)) && RowKept(i))
                  {
                    RowNum(nnz) = i;
                    ColNum(nnz) = A.IndexImag(i, ji);
                    Value(nnz) = T0(0, A.ValueImag(i, ji));
                    nnz++;
                  }
                
                if (A.IndexImag(i, ji) != i)
                  if (RowKept(A.IndexImag(i, ji)) && ColKept(i))
                    {
                      RowNum(nnz) = A.IndexImag(i, ji);
                      ColNum(nnz) = i;
                      Value(nnz) = T0(0, A.ValueImag(i, ji));
                      nnz++;
                    }
                
                ji++;
              }
            
            if (ColKept(jcol) && RowKept(i))
              {
                RowNum(nnz) = i;
                ColNum(nnz) = jcol;
                if ((ji < size_imag) && (A.IndexImag(i, ji) == jcol))
                  Value(nnz)
		    = T0(A.ValueReal(i, jr), A.ValueImag(i, ji));
                else
                  Value(nnz) = T0(A.ValueReal(i, jr), 0);
                
                nnz++;
              }
            
            if (jcol != i)
              if (RowKept(jcol) && ColKept(i))
                {
                  RowNum(nnz) = jcol;
                  ColNum(nnz) = i;
                  if ((ji < size_imag) && (A.IndexImag(i, ji) == jcol))
                    Value(nnz)
		      = T0(A.ValueReal(i, jr), A.ValueImag(i, ji));
                  else
                    Value(nnz) = T0(A.ValueReal(i, jr), 0);
                  
                  nnz++;
                }
            
            if ((ji < size_imag) && (A.IndexImag(i, ji) == jcol))
              ji++;
          }
        
        while (ji < size_imag)
          {
            if (ColKept(A.IndexImag(i, ji)) && RowKept(i))
              {
                RowNum(nnz) = i;
                ColNum(nnz) = A.IndexImag(i, ji);
                Value(nnz) = T0(0, A.ValueImag(i, ji));
                nnz++;
              }
            
            if (A.IndexImag(i, ji) != i)
              if (RowKept(A.IndexImag(i, ji)) && ColKept(i))
                {
                  RowNum(nnz) = A.IndexImag(i, ji);
                  ColNum(nnz) = i;
                  Value(nnz) = T0(0, A.ValueImag(i, ji));
                  nnz++;
                }
            
            ji++;
          }
      }
  }
  

  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<complex<T1>, VectFull, Allocator1>& Value)
  {
    int m = A.GetM(), n = A.GetN();
    if ((m <= 0) || (n <= 0) || (row.GetM() <= 0) || (col.GetM() <= 0))
      {
        RowNum.Clear(); ColNum.Clear(); Value.Clear();
        return;
      }

    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T0>::Treal* data_real = A.GetRealData();    
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T0>::Treal* data_imag = A.GetImagData();    
    
    Vector<bool> RowKept(m), ColKept(n);
    RowKept.Fill(false); ColKept.Fill(false);
    for (int i = 0; i < row.GetM(); i++)
      RowKept(row(i)) = true;
    
    for (int i = 0; i < col.GetM(); i++)
      ColKept(col(i)) = true;
    
    // counting the number of non-zero elements to keep
    long nnz = 0;
    for (int i = 0; i < m; i++)
      if (RowKept(i))
        {
          long ji = ptr_imag[i];
          for (long jr = ptr_real[i]; jr < ptr_real[i+1]; jr++)
            {
              int jcol = ind_real[jr];
              while ((ji < ptr_imag[i+1]) && (ind_imag[ji] < jcol))
                {
                  if (ColKept(ind_imag[ji]))
                    nnz++;
                  
                  ji++;
                }
              
              if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == jcol))
                ji++;
              
              if (ColKept(jcol))
                nnz++;
            }
          
          while (ji < ptr_imag[i+1])
            {
              if (ColKept(ind_imag[ji]))
                nnz++;
              
              ji++;
            }
        }
    
    RowNum.Reallocate(nnz);
    ColNum.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    // then filling the arrays RowNum, ColNum, Value
    for (int i = 0; i < A.GetM(); i++)
      if (RowKept(i))
        {
          long ji = ptr_imag[i];
          for (long jr = ptr_real[i]; jr < ptr_real[i+1]; jr++)
            {
              int jcol = ind_real[jr];
              while ((ji < ptr_imag[i+1]) && (ind_imag[ji] < jcol))
                {
                  if (ColKept(ind_imag[ji]))
                    {
                      RowNum(nnz) = i;
                      ColNum(nnz) = ind_imag[ji];
                      Value(nnz) = T0(0, data_imag[ji]);
                      nnz++;
                    }
                  
                  ji++;
                }
              
              if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == jcol))
                {                  
                  if (ColKept(jcol))
                    {
                      RowNum(nnz) = i;
                      ColNum(nnz) = jcol;
                      Value(nnz) = T0(data_real[jr], data_imag[ji]);
                      nnz++;
                    }
                  
                  ji++;
                }
              else
                {
                  if (ColKept(jcol))
                    {
                      RowNum(nnz) = i;
                      ColNum(nnz) = jcol;
                      Value(nnz) = T0(data_real[jr], 0);
                      nnz++;
                    }
                }
            }
          
          while (ji < ptr_imag[i+1])
            {
              if (ColKept(ind_imag[ji]))
                {
                  RowNum(nnz) = i;
                  ColNum(nnz) = ind_imag[ji];
                  Value(nnz) = T0(0, data_imag[ji]);
                  nnz++;
                }
              
              ji++;
            }
        }
  }
  

  //! extracts some rows/columns of a matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0,
		     RowSymComplexSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<complex<T1>, VectFull, Allocator1>& Value)
  {
    int m = A.GetM(), n = A.GetN();
    if ((m <= 0) || (n <= 0) || (row.GetM() <= 0) || (col.GetM() <= 0))
      {
        RowNum.Clear(); ColNum.Clear(); Value.Clear();
        return;
      }
    
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T0>::Treal* data_real = A.GetRealData();    
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T0>::Treal* data_imag = A.GetImagData();    
    
    Vector<bool> RowKept(m), ColKept(n);
    RowKept.Fill(false); ColKept.Fill(false);
    for (int i = 0; i < row.GetM(); i++)
      RowKept(row(i)) = true;
    
    for (int i = 0; i < col.GetM(); i++)
      ColKept(col(i)) = true;
    
    // counting the number of non-zero elements to keep
    long nnz = 0;
    for (int i = 0; i < m; i++)
      {
        long ji = ptr_imag[i];
        for (long jr = ptr_real[i]; jr < ptr_real[i+1]; jr++)
          {
            int jcol = ind_real[jr];
            while ((ji < ptr_imag[i+1]) && (ind_imag[ji] < jcol))
              {                              
                if (ColKept(ind_imag[ji]) && RowKept(i))
                  nnz++;
                
                if (ind_imag[ji] != i)
                  if (RowKept(ind_imag[ji]) && ColKept(i))
                    nnz++;
                
                ji++;
              }
            
            if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == jcol))
              ji++;
            
            if (ColKept(jcol) && RowKept(i))
              nnz++;
            
            if (jcol != i)
              if (RowKept(jcol) && ColKept(i))
                nnz++;
          }
        
        while (ji < ptr_imag[i+1])
          {
            if (ColKept(ind_imag[ji]) && RowKept(i))
              nnz++;
            
            if (ind_imag[ji] != i)
              if (RowKept(ind_imag[ji]) && ColKept(i))
                nnz++;
            
            ji++;
          }
      }
    
    RowNum.Reallocate(nnz);
    ColNum.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < m; i++)
      {
        long ji = ptr_imag[i];
        for (long jr = ptr_real[i]; jr < ptr_real[i+1]; jr++)
          {
            int jcol = ind_real[jr];
            while ((ji < ptr_imag[i+1]) && (ind_imag[ji] < jcol))
              {                              
                if (ColKept(ind_imag[ji]) && RowKept(i))
                  {
                    RowNum(nnz) = i;
                    ColNum(nnz) = ind_imag[ji];
                    Value(nnz) = T0(0, data_imag[ji]);
                    nnz++;
                  }
                
                if (ind_imag[ji] != i)
                  if (RowKept(ind_imag[ji]) && ColKept(i))
                    {
                      RowNum(nnz) = ind_imag[ji];
                      ColNum(nnz) = i;
                      Value(nnz) = T0(0, data_imag[ji]);
                      nnz++;
                    }
                
                ji++;
              }
            
            if (ColKept(jcol) && RowKept(i))
              {
                RowNum(nnz) = i;
                ColNum(nnz) = jcol;
                if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == jcol))
                  Value(nnz) = T0(data_real[jr], data_imag[ji]);
                else
                  Value(nnz) = T0(data_real[jr], 0);
                
                nnz++;
              }
            
            if (jcol != i)
              if (RowKept(jcol) && ColKept(i))
                {
                  RowNum(nnz) = jcol;
                  ColNum(nnz) = i;
                  if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == jcol))
                    Value(nnz) = T0(data_real[jr], data_imag[ji]);
                  else
                    Value(nnz) = T0(data_real[jr], 0);
                  
                  nnz++;
                }
            
            if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == jcol))
              ji++;
          }
        
        while (ji < ptr_imag[i+1])
          {
            if (ColKept(ind_imag[ji]) && RowKept(i))
              {
                RowNum(nnz) = i;
                ColNum(nnz) = ind_imag[ji];
                Value(nnz) = T0(0, data_imag[ji]);
                nnz++;
              }
            
            if (ind_imag[ji] != i)
              if (RowKept(ind_imag[ji]) && ColKept(i))
                {
                  RowNum(nnz) = ind_imag[ji];
                  ColNum(nnz) = i;
                  Value(nnz) = T0(0, data_imag[ji]);
                  nnz++;
                }
            
            ji++;
          }
      }
  }
  
}

#define SELDON_FILE_FUNCTIONS_MATRIX_COMPLEX_CXX
#endif

