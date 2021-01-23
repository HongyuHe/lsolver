// Copyright (C) 2003-2011 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_COMPLEX_CONVERSIONS_CXX


#include "Matrix_ComplexConversions.hxx"

/*
  Same functions as in Matrix_Conversions.cxx 
  for complex matrices (RowComplexSparse, etc)
 */

namespace Seldon
{
  
  //! Conversion from RowComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    long nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    long* real_ptr = A.GetRealPtr();
    long* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();
    for (int i = 0; i < m; i++)
      {
        long num_old = nnz;
	long num_i = imag_ptr[i];
        for (long j = real_ptr[i]; j < real_ptr[i+1]; j++)
          {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		// imaginary part alone
		IndCol(nnz) = imag_ind[num_i] + index;
		Val(nnz) = T(0, imag_data[num_i]);
		nnz++; num_i++;
	      }
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = real_ind[j] + index;
		Val(nnz) = T(real_data[j], imag_data[num_i]);
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = real_ind[j] + index;
		Val(nnz) = T(real_data[j], 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
          {
            IndCol(nnz) = imag_ind[j] + index;
            Val(nnz) = T(0, imag_data[j]);
            nnz++;
          }

	// filling row numbers
        for (long j = num_old; j < nnz; j++)
	  IndRow(j) = index + i;
      }
    
    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }

  
  //! Conversion from ColComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int n = A.GetN();
    long nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    long* real_ptr = A.GetRealPtr();
    long* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();
    for (int i = 0; i < n; i++)
      {
	long num_old = nnz;
        long num_i = imag_ptr[i];
        for (long j = real_ptr[i]; j < real_ptr[i+1]; j++)
          {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		// imaginary part alone
		IndRow(nnz) = imag_ind[num_i] + index;
		Val(nnz) = T(0, imag_data[num_i]);
		nnz++; num_i++;
	      }
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		// real and imaginary part are both present
		IndRow(nnz) = real_ind[j] + index;
		Val(nnz) = T(real_data[j], imag_data[num_i]);
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndRow(nnz) = real_ind[j] + index;
		Val(nnz) = T(real_data[j], 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
          {
            IndRow(nnz) = imag_ind[j] + index;
            Val(nnz) = T(0, imag_data[j]);
            nnz++;
          }

	// filling column numbers
        for (long j = num_old; j < nnz; j++)
	  IndCol(j) = index + i;
      }
    
    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }

  
  //! Conversion from ArrayRowSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    long nnz = A.GetDataSize();
    long* real_ptr = A.GetRealPtr();
    long* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();
    
    if (sym)
      {
	// first we count the number of non-zero entries
 	// added by symmetry for each row
	// nnz : total number of non-zero entries
	Vector<int> NumNnzRow(m); Vector<long> PtrLower(m+1), PtrUpper(m);
	NumNnzRow.Zero();
	PtrLower(0) = 0;
	PtrUpper(0) = 0;
	nnz = 0;
	for (int i = 0; i < m; i++)
	  {
	    long num_i = imag_ptr[i], nnz_old = nnz;
	    for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	      {
		while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
		  {
		    // imaginary part alone
		    if (imag_ind[num_i] != i)
		      NumNnzRow(imag_ind[num_i])++;
		    
		    nnz++; num_i++;
		  }

			    
		if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
		  {
		    // real and imaginary part are both present
		    if (real_ind[j] != i)
		      NumNnzRow(real_ind[j])++;
		    
		    num_i++; nnz++;
		  }
		else
		  {
		    // real part alone
		    if (real_ind[j] != i)
		      NumNnzRow(real_ind[j])++;

		    nnz++;
		  }
	      }
	    
	    // last values of imaginary part
	    nnz += imag_ptr[i+1] - num_i;
	    for (long j = num_i; j < imag_ptr[i+1]; j++)
	      if (imag_ind[j] != i)
		NumNnzRow(imag_ind[j])++;
          
	    PtrLower(i+1) = PtrLower(i) + NumNnzRow(i) + nnz - nnz_old;
	    if (i < m-1)
	      PtrUpper(i+1) = PtrLower(i+1) + NumNnzRow(i+1);
          }
	
	// total number of non-zero entries
	nnz = PtrLower(m);
	
	// arrays are filled already sorted by rows
	// PtrLower and PtrUpper are incremented progressively
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (int i = 0; i < m; i++)
	  {
	    long num_i = imag_ptr[i]; int jcol = 0;
	    for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	      {
		while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
		  {
		    jcol = imag_ind[num_i];
		    // imaginary part alone
		    if (jcol != i)
		      {
			IndRow(PtrLower(jcol)) = jcol + index;
			IndCol(PtrLower(jcol)) = i + index;
			Val(PtrLower(jcol)) = T(0, imag_data[num_i]);
			PtrLower(jcol)++;
		      }
		    
		    IndRow(PtrUpper(i)) = i + index;
		    IndCol(PtrUpper(i)) = jcol + index;
		    Val(PtrUpper(i)) = T(0, imag_data[num_i]);
		    PtrUpper(i)++;
		    num_i++;
		  }

			    
		if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
		  {
		    jcol = real_ind[j];
		    // real and imaginary part are both present
		    if (real_ind[j] != i)
		      {
			IndRow(PtrLower(jcol)) = jcol + index;
			IndCol(PtrLower(jcol)) = i + index;
			Val(PtrLower(jcol)) = T(real_data[j], imag_data[num_i]);
			PtrLower(jcol)++;
		      }
		   
		    IndRow(PtrUpper(i)) = i + index;
		    IndCol(PtrUpper(i)) = jcol + index;
		    Val(PtrUpper(i)) = T(real_data[j], imag_data[num_i]);
		    PtrUpper(i)++; 
		    num_i++;
		  }
		else
		  {
		    jcol = real_ind[j];
		    // real part alone
		    if (real_ind[j] != i)
		      {
			IndRow(PtrLower(jcol)) = jcol + index;
			IndCol(PtrLower(jcol)) = i + index;
			Val(PtrLower(jcol)) = T(real_data[j], 0);
			PtrLower(jcol)++;
		      }

		    IndRow(PtrUpper(i)) = i + index;
		    IndCol(PtrUpper(i)) = jcol + index;
		    Val(PtrUpper(i)) = T(real_data[j], 0);
		    PtrUpper(i)++; 
		  }
	      }
	    
	    // last values of imaginary part
	    for (long j = num_i; j < imag_ptr[i+1]; j++)
	      {
		jcol = imag_ind[j];
		if (imag_ind[j] != i)
		  {
		    IndRow(PtrLower(jcol)) = jcol + index;
		    IndCol(PtrLower(jcol)) = i + index;
		    Val(PtrLower(jcol)) = T(0, imag_data[j]);
		    PtrLower(jcol)++;
		  }

		IndRow(PtrUpper(i)) = i + index;
		IndCol(PtrUpper(i)) = jcol + index;
		Val(PtrUpper(i)) = T(0, imag_data[j]);
		PtrUpper(i)++; 
	      }
          }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
        nnz = 0;
	for (int i = 0; i < m; i++)
	  {
	    long num_old = nnz;
	    long num_i = imag_ptr[i];
	    for (long j = real_ptr[i]; j < real_ptr[i+1]; j++)
	      {
		while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
		  {
		    // imaginary part alone
		    IndCol(nnz) = imag_ind[num_i] + index;
		    Val(nnz) = T(0, imag_data[num_i]);
		    nnz++; num_i++;
		  }
		
		if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
		  {
		    // real and imaginary part are both present
		    IndCol(nnz) = real_ind[j] + index;
		    Val(nnz) = T(real_data[j], imag_data[num_i]);
		    num_i++; nnz++;
		  }
		else
		  {
		    // real part alone
		    IndCol(nnz) = real_ind[j] + index;
		    Val(nnz) = T(real_data[j], 0);
		    nnz++;
		  }
	      }
	    
	    // last values of imaginary part
	    for (long j = num_i; j < imag_ptr[i+1]; j++)
	      {
		IndCol(nnz) = imag_ind[j] + index;
		Val(nnz) = T(0, imag_data[j]);
		nnz++;
	      }
	    
	    // filling row numbers
	    for (long j = num_old; j < nnz; j++)
	      IndRow(j) = index + i;
	  }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }

  
  //! Conversion from ColSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    long nnz = A.GetDataSize();
    long* real_ptr = A.GetRealPtr();
    long* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();

    if (sym)
      {
	// first we count the number of non-zero entries
 	// added by symmetry for each column
	// nnz : total number of non-zero entries
	Vector<int> NumNnzCol(m); Vector<long> PtrUpper(m+1), PtrLower(m);
	NumNnzCol.Zero();
	nnz = 0;
	for (int i = 0; i < m; i++)
	  {
	    long num_i = imag_ptr[i], nnz_old = nnz;
	    for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	      {
		while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
		  {
		    // imaginary part alone
		    if (imag_ind[num_i] != i)
		      NumNnzCol(imag_ind[num_i])++;
		    
		    nnz++; num_i++;
		  }
		
		if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
		  {
		    // real and imaginary part are both present
		    if (real_ind[j] != i)
		      NumNnzCol(real_ind[j])++;
		    
		    num_i++; nnz++;
		  }
		else
		  {
		    // real part alone
		    if (real_ind[j] != i)
		      NumNnzCol(real_ind[j])++;

		    nnz++;
		  }
	      }
	    
	    // last values of imaginary part
	    nnz += imag_ptr[i+1] - num_i;
	    for (long j = num_i; j < imag_ptr[i+1]; j++)
	      if (imag_ind[j] != i)
		NumNnzCol(imag_ind[j])++;
          
	    PtrUpper(i+1) = nnz - nnz_old;
          }
	
	// PtrUpper and PtrLower are computed
	PtrUpper(0) = 0;
	PtrLower(0) = PtrUpper(1);
	for (int i = 0; i < m; i++)
	  {
	    long ncol_upper = PtrUpper(i+1);
	    int ncol_lower = NumNnzCol(i);
	    
	    PtrUpper(i+1) = PtrUpper(i) + ncol_upper + ncol_lower;
	    if (i < m-1)
	      PtrLower(i+1) = PtrUpper(i+1) + PtrUpper(i+2);
	  }
	
	// total number of non-zero entries
	nnz = PtrUpper(m);
	
	// arrays are filled already sorted by rows
	// PtrLower and PtrUpper are incremented progressively
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (int i = 0; i < m; i++)
	  {
	    long num_i = imag_ptr[i]; int jrow = 0;
	    for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	      {
		while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
		  {
		    jrow = imag_ind[num_i];
		    // imaginary part alone
		    if (jrow != i)
		      {
			IndRow(PtrLower(jrow)) = i + index;
			IndCol(PtrLower(jrow)) = jrow + index;
			Val(PtrLower(jrow)) = T(0, imag_data[num_i]);
			PtrLower(jrow)++;
		      }
		    
		    IndRow(PtrUpper(i)) = jrow + index;
		    IndCol(PtrUpper(i)) = i + index;
		    Val(PtrUpper(i)) = T(0, imag_data[num_i]);
		    PtrUpper(i)++;
		    num_i++;
		  }

			    
		if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
		  {
		    jrow = real_ind[j];
		    // real and imaginary part are both present
		    if (real_ind[j] != i)
		      {
			IndRow(PtrLower(jrow)) = i + index;
			IndCol(PtrLower(jrow)) = jrow + index;
			Val(PtrLower(jrow)) = T(real_data[j], imag_data[num_i]);
			PtrLower(jrow)++;
		      }
		   
		    IndRow(PtrUpper(i)) = jrow + index;
		    IndCol(PtrUpper(i)) = i + index;
		    Val(PtrUpper(i)) = T(real_data[j], imag_data[num_i]);
		    PtrUpper(i)++; 
		    num_i++;
		  }
		else
		  {
		    jrow = real_ind[j];
		    // real part alone
		    if (real_ind[j] != i)
		      {
			IndRow(PtrLower(jrow)) = i + index;
			IndCol(PtrLower(jrow)) = jrow + index;
			Val(PtrLower(jrow)) = T(real_data[j], 0);
			PtrLower(jrow)++;
		      }

		    IndRow(PtrUpper(i)) = jrow + index;
		    IndCol(PtrUpper(i)) = i + index;
		    Val(PtrUpper(i)) = T(real_data[j], 0);
		    PtrUpper(i)++; 
		  }
	      }
	    
	    // last values of imaginary part
	    for (long j = num_i; j < imag_ptr[i+1]; j++)
	      {
		jrow = imag_ind[j];
		if (imag_ind[j] != i)
		  {
		    IndRow(PtrLower(jrow)) = i + index;
		    IndCol(PtrLower(jrow)) = jrow + index;
		    Val(PtrLower(jrow)) = T(0, imag_data[j]);
		    PtrLower(jrow)++;
		  }

		IndRow(PtrUpper(i)) = jrow + index;
		IndCol(PtrUpper(i)) = i + index;
		Val(PtrUpper(i)) = T(0, imag_data[j]);
		PtrUpper(i)++; 
	      }
          }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	nnz = 0;
	
	for (int i = 0; i < m; i++)
	  {
	    long num_old = nnz;
	    long num_i = imag_ptr[i];
	    for (long j = real_ptr[i]; j < real_ptr[i+1]; j++)
	      {
		while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
		  {
		    // imaginary part alone
		    IndRow(nnz) = imag_ind[num_i] + index;
		    Val(nnz) = T(0, imag_data[num_i]);
		    nnz++; num_i++;
		  }
		
		if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
		  {
		    // real and imaginary part are both present
		    IndRow(nnz) = real_ind[j] + index;
		    Val(nnz) = T(real_data[j], imag_data[num_i]);
		    num_i++; nnz++;
		  }
		else
		  {
		    // real part alone
		    IndRow(nnz) = real_ind[j] + index;
		    Val(nnz) = T(real_data[j], 0);
		    nnz++;
		  }
	      }
	    
	    // last values of imaginary part
	    for (long j = num_i; j < imag_ptr[i+1]; j++)
	      {
		IndRow(nnz) = imag_ind[j] + index;
		Val(nnz) = T(0, imag_data[j]);
		nnz++;
	      }
	    
	    // filling column numbers
	    for (long j = num_old; j < nnz; j++)
	      IndCol(j) = index + i;
	  }
	
	IndRow.Resize(nnz);
	IndCol.Resize(nnz);
	Val.Resize(nnz);	
      }
  }

  
  //! Conversion from ArrayRowComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    long nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < m; i++)
      {
	long num_old = nnz;
	int num_i = 0, size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		IndCol(nnz) = A.IndexImag(i, num_i) + index;
		Val(nnz) = T(0, A.ValueImag(i, num_i));
		nnz++; num_i++;
	      }

            if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = A.IndexReal(i, j) + index;
		Val(nnz) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = A.IndexReal(i, j) + index;
		Val(nnz) = T(A.ValueReal(i, j), 0);
		nnz++;
	      }
          }
	
	for (int j = num_i; j < size_imag; j++)
	  {
	    IndCol(nnz) = A.IndexImag(i, j) + index;
	    Val(nnz) = T(0, A.ValueImag(i, j));
	    nnz++;
	  }
	
	// row numbers for row i
	for (long j = num_old; j < nnz; j++)
	  IndRow(j) = i + index;
      }
    
    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }

  
  //! Conversion from ArrayColComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    long nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < A.GetN(); i++)
      {
	long num_old = nnz;
	int num_i = 0, size_imag = A.GetImagColumnSize(i);
        for (int j = 0; j < A.GetRealColumnSize(i); j++)
          {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		IndRow(nnz) = A.IndexImag(i, num_i) + index;
		Val(nnz) = T(0, A.ValueImag(i, num_i));
		nnz++; num_i++;
	      }
            
            if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		IndRow(nnz) = A.IndexReal(i, j) + index;
		Val(nnz) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndRow(nnz) = A.IndexReal(i, j) + index;
		Val(nnz) = T(A.ValueReal(i, j), 0);
		nnz++;
	      }
          }
	
	for (int j = num_i; j < size_imag; j++)
	  {
	    IndRow(nnz) = A.IndexImag(i, j) + index;
	    Val(nnz) = T(0, A.ValueImag(i, j));
	    nnz++;
	  }
	
	// column numbers for row i
	for (long j = num_old; j < nnz; j++)
	  IndCol(j) = i + index;
      }
    
    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }

  
  //! Conversion from ArrayRowSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymComplexSparse,
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
 	// added by symmetry for each row
	// nnz : total number of non-zero entries
	Vector<int> NumNnzRow(m); Vector<long> PtrLower(m+1), PtrUpper(m);
	NumNnzRow.Zero();
	PtrLower(0) = 0;
	PtrUpper(0) = 0;
	nnz = 0;
	for (int i = 0; i < m; i++)
	  {
	    int num_i = 0; long nnz_old = nnz;
	    int size_imag = A.GetImagRowSize(i);
	    for (int j = 0; j < A.GetRealRowSize(i); j++)
	      {
		while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
		  {
		    // imaginary part alone
		    if (A.IndexImag(i, num_i) != i)
		      NumNnzRow(A.IndexImag(i, num_i))++;
		    
		    nnz++; num_i++;
		  }
					    
		if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
		  {
		    // real and imaginary part are both present
		    if (A.IndexReal(i, j) != i)
		      NumNnzRow(A.IndexReal(i, j))++;
		    
		    num_i++; nnz++;
		  }
		else
		  {
		    // real part alone
		    if (A.IndexReal(i, j) != i)
		      NumNnzRow(A.IndexReal(i, j))++;

		    nnz++;
		  }
	      }
	    
	    // last values of imaginary part
	    nnz += (size_imag-num_i);
	    for (int j = num_i; j < size_imag; j++)
	      if (A.IndexImag(i, j) != i)
		NumNnzRow(A.IndexImag(i, j))++;
          
	    PtrLower(i+1) = PtrLower(i) + NumNnzRow(i) + nnz - nnz_old;
	    if (i < m-1)
	      PtrUpper(i+1) = PtrLower(i+1) + NumNnzRow(i+1);
          }
	
	// total number of non-zero entries
	nnz = PtrLower(m);
	
	// arrays are filled already sorted by rows
	// PtrLower and PtrUpper are incremented progressively
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (int i = 0; i < m; i++)
	  {
	    int num_i = 0, jcol = 0;
	    int size_imag = A.GetImagRowSize(i); 
	    for (int j = 0; j < A.GetRealRowSize(i); j++)
	      {
		while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
		  {
		    jcol = A.IndexImag(i, num_i);
		    // imaginary part alone
		    if (jcol != i)
		      {
			IndRow(PtrLower(jcol)) = jcol + index;
			IndCol(PtrLower(jcol)) = i + index;
			Val(PtrLower(jcol)) = T(0, A.ValueImag(i, num_i));
			PtrLower(jcol)++;
		      }
		    
		    IndRow(PtrUpper(i)) = i + index;
		    IndCol(PtrUpper(i)) = jcol + index;
		    Val(PtrUpper(i)) = T(0, A.ValueImag(i, num_i));
		    PtrUpper(i)++;
		    num_i++;
		  }

			    
		if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
		  {
		    jcol = A.IndexReal(i, j);
		    // real and imaginary part are both present
		    if (jcol != i)
		      {
			IndRow(PtrLower(jcol)) = jcol + index;
			IndCol(PtrLower(jcol)) = i + index;
			Val(PtrLower(jcol)) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
			PtrLower(jcol)++;
		      }
		   
		    IndRow(PtrUpper(i)) = i + index;
		    IndCol(PtrUpper(i)) = jcol + index;
		    Val(PtrUpper(i)) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		    PtrUpper(i)++; 
		    num_i++;
		  }
		else
		  {
		    jcol = A.IndexReal(i, j);
		    // real part alone
		    if (jcol != i)
		      {
			IndRow(PtrLower(jcol)) = jcol + index;
			IndCol(PtrLower(jcol)) = i + index;
			Val(PtrLower(jcol)) = T(A.ValueReal(i, j), 0);
			PtrLower(jcol)++;
		      }
		    
		    IndRow(PtrUpper(i)) = i + index;
		    IndCol(PtrUpper(i)) = jcol + index;
		    Val(PtrUpper(i)) = T(A.ValueReal(i, j), 0);
		    PtrUpper(i)++; 
		  }
	      }
	    
	    // last values of imaginary part
	    for (int j = num_i; j < size_imag; j++)
	      {
		jcol = A.IndexImag(i, j);
		if (jcol != i)
		  {
		    IndRow(PtrLower(jcol)) = jcol + index;
		    IndCol(PtrLower(jcol)) = i + index;
		    Val(PtrLower(jcol)) = T(0, A.ValueImag(i, j));
		    PtrLower(jcol)++;
		  }

		IndRow(PtrUpper(i)) = i + index;
		IndCol(PtrUpper(i)) = jcol + index;
		Val(PtrUpper(i)) = T(0, A.ValueImag(i, j));
		PtrUpper(i)++; 
	      }
          }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
        nnz = 0;
	for (int i = 0; i < m; i++)
	  {
	    long num_old = nnz;
	    int num_i = 0, size_imag = A.GetImagRowSize(i);
	    for (int j = 0; j < A.GetRealRowSize(i); j++)
	      {
		while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
		  {
		    // imaginary part alone
		    IndCol(nnz) = A.IndexImag(i, num_i) + index;
		    Val(nnz) = T(0, A.ValueImag(i, num_i));
		    nnz++; num_i++;
		  }
		
		if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
		  {
		    // real and imaginary part are both present
		    IndCol(nnz) = A.IndexReal(i, j) + index;
		    Val(nnz) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		    num_i++; nnz++;
		  }
		else
		  {
		    // real part alone
		    IndCol(nnz) = A.IndexReal(i, j) + index;
		    Val(nnz) = T(A.ValueReal(i, j), 0);
		    nnz++;
		  }
	      }
	    
	    // last values of imaginary part
	    for (int j = num_i; j < size_imag; j++)
	      {
		IndCol(nnz) = A.IndexImag(i, j) + index;
		Val(nnz) = T(0, A.ValueImag(i, j));
		nnz++;
	      }
	    
	    // filling row numbers
	    for (long j = num_old; j < nnz; j++)
	      IndRow(j) = index + i;
	  }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }
  
  
  //! Conversion from ArrayColSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymComplexSparse,
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
	// nnz : total number of non-zero entries
	Vector<int> NumNnzCol(m); Vector<long> PtrUpper(m+1), PtrLower(m);
	NumNnzCol.Zero();
	nnz = 0;
	for (int i = 0; i < m; i++)
	  {
	    int num_i = 0; long nnz_old = nnz;
	    int size_imag = A.GetImagColumnSize(i);
	    for (int j = 0; j < A.GetRealColumnSize(i); j++)
	      {
		while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
		  {
		    // imaginary part alone
		    if (A.IndexImag(i, num_i) != i)
		      NumNnzCol(A.IndexImag(i, num_i))++;
		    
		    nnz++; num_i++;
		  }
		
		if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
		  {
		    // real and imaginary part are both present
		    if (A.IndexReal(i, j) != i)
		      NumNnzCol(A.IndexReal(i, j))++;
		    
		    num_i++; nnz++;
		  }
		else
		  {
		    // real part alone
		    if (A.IndexReal(i, j) != i)
		      NumNnzCol(A.IndexReal(i, j))++;
		    
		    nnz++;
		  }
	      }
	    
	    // last values of imaginary part
	    nnz += size_imag - num_i;
	    for (int j = num_i; j < size_imag; j++)
	      if (A.IndexImag(i, j) != i)
		NumNnzCol(A.IndexImag(i, j))++;
          
	    PtrUpper(i+1) = nnz - nnz_old;
          }
	
	// PtrUpper and PtrLower are computed
	PtrUpper(0) = 0;
	PtrLower(0) = PtrUpper(1);
	for (int i = 0; i < m; i++)
	  {
	    long ncol_upper = PtrUpper(i+1);
	    int ncol_lower = NumNnzCol(i);
	    
	    PtrUpper(i+1) = PtrUpper(i) + ncol_upper + ncol_lower;
	    if (i < m-1)
	      PtrLower(i+1) = PtrUpper(i+1) + PtrUpper(i+2);
	  }
	
	// total number of non-zero entries
	nnz = PtrUpper(m);
	
	// arrays are filled already sorted by rows
	// PtrLower and PtrUpper are incremented progressively
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (int i = 0; i < m; i++)
	  {
	    int num_i = 0, jrow = 0;
	    int size_imag = A.GetImagColumnSize(i);
	    for (int j = 0; j < A.GetRealColumnSize(i); j++)
	      {
		while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
		  {
		    jrow = A.IndexImag(i, num_i);
		    // imaginary part alone
		    if (jrow != i)
		      {
			IndRow(PtrLower(jrow)) = i + index;
			IndCol(PtrLower(jrow)) = jrow + index;
			Val(PtrLower(jrow)) = T(0, A.ValueImag(i, num_i));
			PtrLower(jrow)++;
		      }
		    
		    IndRow(PtrUpper(i)) = jrow + index;
		    IndCol(PtrUpper(i)) = i + index;
		    Val(PtrUpper(i)) = T(0, A.ValueImag(i, num_i));
		    PtrUpper(i)++;
		    num_i++;
		  }
		
		if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
		  {
		    jrow = A.IndexReal(i, j);
		    // real and imaginary part are both present
		    if (jrow != i)
		      {
			IndRow(PtrLower(jrow)) = i + index;
			IndCol(PtrLower(jrow)) = jrow + index;
			Val(PtrLower(jrow)) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
			PtrLower(jrow)++;
		      }
		   
		    IndRow(PtrUpper(i)) = jrow + index;
		    IndCol(PtrUpper(i)) = i + index;
		    Val(PtrUpper(i)) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		    PtrUpper(i)++; 
		    num_i++;
		  }
		else
		  {
		    jrow = A.IndexReal(i, j);
		    // real part alone
		    if (jrow != i)
		      {
			IndRow(PtrLower(jrow)) = i + index;
			IndCol(PtrLower(jrow)) = jrow + index;
			Val(PtrLower(jrow)) = T(A.ValueReal(i, j), 0);
			PtrLower(jrow)++;
		      }

		    IndRow(PtrUpper(i)) = jrow + index;
		    IndCol(PtrUpper(i)) = i + index;
		    Val(PtrUpper(i)) = T(A.ValueReal(i, j), 0);
		    PtrUpper(i)++; 
		  }
	      }
	    
	    // last values of imaginary part
	    for (int j = num_i; j < size_imag; j++)
	      {
		jrow = A.IndexImag(i, j);
		if (jrow != i)
		  {
		    IndRow(PtrLower(jrow)) = i + index;
		    IndCol(PtrLower(jrow)) = jrow + index;
		    Val(PtrLower(jrow)) = T(0, A.ValueImag(i, j));
		    PtrLower(jrow)++;
		  }

		IndRow(PtrUpper(i)) = jrow + index;
		IndCol(PtrUpper(i)) = i + index;
		Val(PtrUpper(i)) = T(0, A.ValueImag(i, j));
		PtrUpper(i)++; 
	      }
          }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	nnz = 0;

	for (int i = 0; i < m; i++)
	  {
	    long num_old = nnz;
	    int num_i = 0, size_imag = A.GetImagColumnSize(i);
	    for (int j = 0; j < A.GetRealColumnSize(i); j++)
	      {
		while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
		  {
		    // imaginary part alone
		    IndRow(nnz) = A.IndexImag(i, num_i) + index;
		    Val(nnz) = T(0, A.ValueImag(i, num_i));
		    nnz++; num_i++;
		  }
		
		if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
		  {
		    // real and imaginary part are both present
		    IndRow(nnz) = A.IndexReal(i, j) + index;
		    Val(nnz) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		    num_i++; nnz++;
		  }
		else
		  {
		    // real part alone
		    IndRow(nnz) = A.IndexReal(i, j) + index;
		    Val(nnz) = T(A.ValueReal(i, j), 0);
		    nnz++;
		  }
	      }
	    
	    for (int j = num_i; j < size_imag; j++)
	      {
		IndRow(nnz) = A.IndexImag(i, j) + index;
		Val(nnz) = T(0, A.ValueImag(i, j));
		nnz++;
	      }
	    
	    // column numbers for row i
	    for (long j = num_old; j < nnz; j++)
	      IndCol(j) = i + index;
	  }
	
	IndRow.Resize(nnz);
	IndCol.Resize(nnz);
	Val.Resize(nnz);	
      }
  }
  
  
  //! Conversion from coordinate format to RowComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);

    // the size of the matrix is detected from input arrays
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we don't change this size
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // the number of elements for each row is computed
    Vector<int> NumRealNnzRow(m), NumImagNnzRow(m);
    NumRealNnzRow.Zero(); NumImagNnzRow.Zero();

    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	if (real(Val(i)) != zero)
	  NumRealNnzRow(irow)++;

	if (imag(Val(i)) != zero)
	  NumImagNnzRow(irow)++;
      }

    // the arrays PtrReal and PtrImag can be constructed
    Vector<long> PtrReal(m+1), PtrImag(m+1);
    PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < m; i++)
      {
        PtrReal(i+1) = PtrReal(i) + NumRealNnzRow(i);
        PtrImag(i+1) = PtrImag(i) + NumImagNnzRow(i);
      }    
    
    long real_nz = PtrReal(m), imag_nz = PtrImag(m);
    
    // Fills matrix 'A'.
    Vector<int> IndReal(real_nz), IndImag(imag_nz);
    Vector<Treal, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    NumRealNnzRow.Zero(); NumImagNnzRow.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (real(Val(i)) != zero)
	  {
	    long num_r = PtrReal(irow) + NumRealNnzRow(irow);
	    IndReal(num_r) = icol;
	    ValReal(num_r) = real(Val(i));
	    NumRealNnzRow(irow)++;
	  }
	
	if (imag(Val(i)) != zero)
	  {
	    long num_i = PtrImag(irow) + NumImagNnzRow(irow);
	    IndImag(num_i) = icol;
	    ValImag(num_i) = imag(Val(i));
	    NumImagNnzRow(irow)++;
	  }		
      }
    
    // column numbers are sorted
    for (int i = 0; i < m; i++)
      {
	Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
	Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
      }
    
    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }
  
  
  //! Conversion from coordinate format to ColComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);

    // the size of the matrix is detected from input arrays
    int row_max = IndRow.GetNormInf() - index;
    int col_max = IndCol.GetNormInf() - index;
    int m = row_max + 1;
    int n = col_max + 1;

    // if A is already allocated, we don't change this size
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // the number of elements for each column is computed
    Vector<int> NumRealNnzCol(n), NumImagNnzCol(n);
    NumRealNnzCol.Zero(); NumImagNnzCol.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      {
	int irow = IndCol(i) - index;
	if (real(Val(i)) != zero)
	  NumRealNnzCol(irow)++;

	if (imag(Val(i)) != zero)
	  NumImagNnzCol(irow)++;
      }

    // the arrays PtrReal and PtrImag can be constructed
    Vector<long> PtrReal(n+1), PtrImag(n+1);
    PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < n; i++)
      {
        PtrReal(i+1) = PtrReal(i) + NumRealNnzCol(i);
        PtrImag(i+1) = PtrImag(i) + NumImagNnzCol(i);
      }    
    
    long real_nz = PtrReal(n), imag_nz = PtrImag(n);
    
    // Fills matrix 'A'.
    Vector<int> IndReal(real_nz), IndImag(imag_nz);
    Vector<Treal, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    NumRealNnzCol.Zero(); NumImagNnzCol.Zero();
    for (long i = 0; i < IndCol.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (real(Val(i)) != zero)
	  {
	    long num_r = PtrReal(icol) + NumRealNnzCol(icol);
	    IndReal(num_r) = irow;
	    ValReal(num_r) = real(Val(i));
	    NumRealNnzCol(icol)++;
	  }
	
	if (imag(Val(i)) != zero)
	  {
	    long num_i = PtrImag(icol) + NumImagNnzCol(icol);
	    IndImag(num_i) = irow;
	    ValImag(num_i) = imag(Val(i));
	    NumImagNnzCol(icol)++;
	  }		
      }
    
    // row numbers are sorted
    for (int i = 0; i < n; i++)
      {
	Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
	Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
      }
    
    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! Conversion from coordinate format to RowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);

    // the number of row and columns is detected
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    
    // if A is already allocated, we take the number of rows/columns of A
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // Number of elements per row.
    Vector<long> PtrReal(m+1), PtrImag(m+1);
    PtrReal.Zero(); PtrImag.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (irow <= icol)
          {
	    // only upper part of the matrix is kept
            if (real(Val(i)) != zero)
              PtrReal(irow+1)++;
            
            if (imag(Val(i)) != zero)
              PtrImag(irow+1)++;
          }
      }

    // PtrReal and PtrImag are constructed
    for (int i = 0; i < m; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }
    
    long real_nz = PtrReal(m), imag_nz = PtrImag(m);
    
    Vector<int> NumRealNnzRow(m), NumImagNnzRow(m);
    NumRealNnzRow.Zero(); NumImagNnzRow.Zero();
    
    // Fills matrix 'A'.
    Vector<int> IndReal(real_nz), IndImag(imag_nz);
    Vector<Treal, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (irow <= icol)
	  {
	    if (real(Val(i)) != zero)
	      {
		long num_r = PtrReal(irow) + NumRealNnzRow(irow);
		IndReal(num_r) = icol;
		ValReal(num_r) = real(Val(i));
		NumRealNnzRow(irow)++;
	      }
	    
	    if (imag(Val(i)) != zero)
	      {
		long num_i = PtrImag(irow) + NumImagNnzRow(irow);
		IndImag(num_i) = icol;
		ValImag(num_i) = imag(Val(i));
		NumImagNnzRow(irow)++;
	      }		
	  }
      }
    
    // column numbers are sorted
    for (int i = 0; i < m; i++)
      {
	Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
	Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
      }
    
    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! Conversion from coordinate format to ColSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);

    // the number of row and columns is detected
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // if A is already allocated, we take the number of rows/columns of A
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // Number of elements per column.
    Vector<long> PtrReal(n+1), PtrImag(n+1);
    PtrReal.Zero(); PtrImag.Zero();
    for (long i = 0; i < IndCol.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (irow <= icol)
          {
	    // only upper part of the matrix is kept
            if (real(Val(i)) != zero)
              PtrReal(icol + 1)++;
            
            if (imag(Val(i)) != zero)
              PtrImag(icol + 1)++;
          }
      }

    // PtrReal and PtrImag are constructed
    for (int i = 0; i < n; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }
    
    long real_nz = PtrReal(n), imag_nz = PtrImag(n);

    Vector<int> NumRealNnzCol(n), NumImagNnzCol(n);
    NumRealNnzCol.Zero(); NumImagNnzCol.Zero();
    
    // Fills matrix 'A'.
    Vector<int> IndReal(real_nz), IndImag(imag_nz);
    Vector<Treal, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    for (long i = 0; i < IndCol.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (irow <= icol)
	  {
	    if (real(Val(i)) != zero)
	      {
		long num_r = PtrReal(icol) + NumRealNnzCol(icol);
		IndReal(num_r) = irow;
		ValReal(num_r) = real(Val(i));
		NumRealNnzCol(icol)++;
	      }
	    
	    if (imag(Val(i)) != zero)
	      {
		long num_i = PtrImag(icol) + NumImagNnzCol(icol);
		IndImag(num_i) = irow;
		ValImag(num_i) = imag(Val(i));
		NumImagNnzCol(icol)++;
	      }		
	  }
      }
    
    // row numbers are sorted
    for (int i = 0; i < n; i++)
      {
	Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
	Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
      }
        
    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! Conversion from coordinate format to ArrayRowComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);

    // the size of the matrix is detected from input arrays
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // if A is already allocated, we don't change this size
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // the number of elements for each row is computed
    Vector<int> NumRealNnzRow(m), NumImagNnzRow(m);
    NumRealNnzRow.Zero(); NumImagNnzRow.Zero();

    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	if (real(Val(i)) != zero)
	  NumRealNnzRow(irow)++;

	if (imag(Val(i)) != zero)
	  NumImagNnzRow(irow)++;
      }
    
    // The matrix A is allocated
    A.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	if (NumRealNnzRow(i) > 0)
	  A.ReallocateRealRow(i, NumRealNnzRow(i));
	
	if (NumImagNnzRow(i) > 0)
	  A.ReallocateImagRow(i, NumImagNnzRow(i));
      }
    
    // Fills matrix 'A'.
    NumRealNnzRow.Zero(); NumImagNnzRow.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (real(Val(i)) != zero)
	  {
	    int num = NumRealNnzRow(irow);
	    A.IndexReal(irow, num) = icol;
	    A.ValueReal(irow, num) = real(Val(i));
	    NumRealNnzRow(irow)++;
	  }
	
	if (imag(Val(i)) != zero)
	  {
	    int num = NumImagNnzRow(irow);
	    A.IndexImag(irow, num) = icol;
	    A.ValueImag(irow, num) = imag(Val(i));
	    NumImagNnzRow(irow)++;
	  }
      }
    
    // column numbers are sorted
    A.Assemble();
  }

  
  //! Conversion from coordinate format to ArrayColComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);

    // the size of the matrix is detected from input arrays
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // if A is already allocated, we don't change this size
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // the number of elements for each column is computed
    Vector<int> NumRealNnzCol(n), NumImagNnzCol(n);
    NumRealNnzCol.Zero(); NumImagNnzCol.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      {
	int icol = IndCol(i) - index;
	if (real(Val(i)) != zero)
	  NumRealNnzCol(icol)++;

	if (imag(Val(i)) != zero)
	  NumImagNnzCol(icol)++;
      }

    // The matrix A is allocated
    A.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
	if (NumRealNnzCol(i) > 0)
	  A.ReallocateRealColumn(i, NumRealNnzCol(i));
	
	if (NumImagNnzCol(i) > 0)
	  A.ReallocateImagColumn(i, NumImagNnzCol(i));
      }

    // Fills matrix 'A'.
    NumRealNnzCol.Zero(); NumImagNnzCol.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (real(Val(i)) != zero)
	  {
	    int num = NumRealNnzCol(icol);
	    A.IndexReal(icol, num) = irow;
	    A.ValueReal(icol, num) = real(Val(i));
	    NumRealNnzCol(icol)++;
	  }
	
	if (imag(Val(i)) != zero)
	  {
	    int num = NumImagNnzCol(icol);
	    A.IndexImag(icol, num) = irow;
	    A.ValueImag(icol, num) = imag(Val(i));
	    NumImagNnzCol(icol)++;
	  }
      }
    
    // row numbers are sorted
    A.Assemble();
  }

  
  //! Conversion from coordinate format to ArrayRowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);

    // the size of the matrix is detected from input arrays
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // if A is already allocated, we don't change this size
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // the number of elements for each row is computed
    Vector<int> NumRealNnzRow(m), NumImagNnzRow(m);
    NumRealNnzRow.Zero(); NumImagNnzRow.Zero();

    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (irow <= icol)
	  {
	    // only upper part of the matrix is kept
	    if (real(Val(i)) != zero)
	      NumRealNnzRow(irow)++;
	    
	    if (imag(Val(i)) != zero)
	      NumImagNnzRow(irow)++;
	  }
      }
    
    // The matrix A is allocated
    A.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	if (NumRealNnzRow(i) > 0)
	  A.ReallocateRealRow(i, NumRealNnzRow(i));
	
	if (NumImagNnzRow(i) > 0)
	  A.ReallocateImagRow(i, NumImagNnzRow(i));
      }
    
    // Fills matrix 'A'.
    NumRealNnzRow.Zero(); NumImagNnzRow.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (irow <= icol)
	  {
	    if (real(Val(i)) != zero)
	      {
		int num = NumRealNnzRow(irow);
		A.IndexReal(irow, num) = icol;
		A.ValueReal(irow, num) = real(Val(i));
		NumRealNnzRow(irow)++;
	      }
	    
	    if (imag(Val(i)) != zero)
	      {
		int num = NumImagNnzRow(irow);
		A.IndexImag(irow, num) = icol;
		A.ValueImag(irow, num) = imag(Val(i));
		NumImagNnzRow(irow)++;
	      }
	  }
      }
    
    // column numbers are sorted
    A.Assemble();
  }

  
  //! Conversion from coordinate format to ArrayRowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);

    // the size of the matrix is detected from input arrays
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // if A is already allocated, we don't change this size
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // the number of elements for each column is computed
    Vector<int> NumRealNnzCol(n), NumImagNnzCol(n);
    NumRealNnzCol.Zero(); NumImagNnzCol.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (irow <= icol)
	  {
	    // only upper part of the matrix is kept
	    if (real(Val(i)) != zero)
	      NumRealNnzCol(icol)++;
	    
	    if (imag(Val(i)) != zero)
	      NumImagNnzCol(icol)++;
	  }
      }

    // The matrix A is allocated
    A.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
	if (NumRealNnzCol(i) > 0)
	  A.ReallocateRealColumn(i, NumRealNnzCol(i));
	
	if (NumImagNnzCol(i) > 0)
	  A.ReallocateImagColumn(i, NumImagNnzCol(i));
      }

    // Fills matrix 'A'.
    NumRealNnzCol.Zero(); NumImagNnzCol.Zero();
    for (long i = 0; i < IndRow.GetM(); i++)
      {
	int irow = IndRow(i) - index;
	int icol = IndCol(i) - index;
	if (irow <= icol)
	  {
	    if (real(Val(i)) != zero)
	      {
		int num = NumRealNnzCol(icol);
		A.IndexReal(icol, num) = irow;
		A.ValueReal(icol, num) = real(Val(i));
		NumRealNnzCol(icol)++;
	      }
	    
	    if (imag(Val(i)) != zero)
	      {
		int num = NumImagNnzCol(icol);
		A.IndexImag(icol, num) = irow;
		A.ValueImag(icol, num) = imag(Val(i));
		NumImagNnzCol(icol)++;
	      }
	  }
      }
    
    // row numbers are sorted
    A.Assemble();
  }

  
  //! Conversion from RowComplexSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
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

    typedef typename ClassComplexType<T>::Treal Treal;
    long* real_ptr_ = A.GetRealPtr();
    int* real_ind_ = A.GetRealInd();
    Treal* real_data_ = A.GetRealData();

    long* imag_ptr_ = A.GetImagPtr();
    int* imag_ind_ = A.GetImagInd();
    Treal* imag_data_ = A.GetImagData();
    
    long real_nnz = A.GetRealDataSize();
    long imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    long nnz = 0;
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	long num_i = imag_ptr_[i];
        for (long j = real_ptr_[i]; j < real_ptr_[i+1]; j++)
          {
	    while ((num_i < imag_ptr_[i+1]) && (imag_ind_[num_i] < real_ind_[j]))
	      {
		// imaginary part alone
		IndCol(nnz) = imag_ind_[num_i];
		Value(nnz) = T(0, imag_data_[num_i]);
		nnz++; num_i++;
	      }
	    
	    if ((num_i < imag_ptr_[i+1]) && (imag_ind_[num_i] == real_ind_[j]))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = real_ind_[j];
		Value(nnz) = T(real_data_[j], imag_data_[num_i]);
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = real_ind_[j];
		Value(nnz) = T(real_data_[j], 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr_[i+1]; j++)
          {
            IndCol(nnz) = imag_ind_[j];
            Value(nnz) = T(0, imag_data_[j]);
            nnz++;
          }
	
	Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from ColComplexSparse to CSR format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();
    long* real_ptr = A.GetRealPtr();
    long* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();

    // first we count the number of non-zero entries for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (int i = 0; i < n; i++)
      {
	long num_i = imag_ptr[i];
	for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		// imaginary part alone
		NumNnzRow(imag_ind[num_i])++;		
		num_i++;
	      }
	    
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		// real and imaginary part are both present
		NumNnzRow(real_ind[j])++;		
		num_i++;
	      }
	    else
	      NumNnzRow(real_ind[j])++;
	  }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
	  NumNnzRow(imag_ind[j])++;	
      }
    
    // 'Ptr' array can be constructed
    Ptr.Reallocate(m+1);
    Ptr(0) = 0;	
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzRow(i);
    
    // total number of non-zero entries
    nnz = Ptr(m);
	
    // arrays 'IndCol' and 'Val' are filled already sorted by rows
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzRow.Zero();
    for (int i = 0; i < n; i++)
      {
	long num_i = imag_ptr[i]; int jcol = 0;
	for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		jcol = imag_ind[num_i];
		// imaginary part alone
		long num = Ptr(jcol) + NumNnzRow(jcol);
		IndCol(num) = i;
		Val(num) = T(0, imag_data[num_i]);
		NumNnzRow(jcol)++;
		num_i++;
	      }
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		jcol = real_ind[j];
		long num = Ptr(jcol) + NumNnzRow(jcol);
		// real and imaginary part are both present
		IndCol(num) = i;
		Val(num) = T(real_data[j], imag_data[num_i]);
		NumNnzRow(jcol)++; 
		num_i++;
	      }
	    else
	      {
		jcol = real_ind[j];
		long num = Ptr(jcol) + NumNnzRow(jcol);
		// real part alone
		IndCol(num) = i;
		Val(num) = T(real_data[j], 0);
		NumNnzRow(jcol)++; 
	      }
	  }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
	  {
	    jcol = imag_ind[j];
	    long num = Ptr(jcol) + NumNnzRow(jcol);
	    IndCol(num) = i;
	    Val(num) = T(0, imag_data[j]);
	    NumNnzRow(jcol)++;
	  }
      }
  }

  
  //! Conversion from RowSymComplexSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
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

    typedef typename ClassComplexType<T>::Treal Treal;
    long* real_ptr_ = A.GetRealPtr();
    int* real_ind_ = A.GetRealInd();
    Treal* real_data_ = A.GetRealData();

    long* imag_ptr_ = A.GetImagPtr();
    int* imag_ind_ = A.GetImagInd();
    Treal* imag_data_ = A.GetImagData();
    
    long real_nnz = A.GetRealDataSize();
    long imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    long nnz = 0;
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	long num_i = imag_ptr_[i];
        for (long j = real_ptr_[i]; j < real_ptr_[i+1]; j++)
          {
	    while ((num_i < imag_ptr_[i+1]) && (imag_ind_[num_i] < real_ind_[j]))
	      {
		// imaginary part alone
		IndCol(nnz) = imag_ind_[num_i];
		Value(nnz) = T(0, imag_data_[num_i]);
		nnz++; num_i++;
	      }
	    
	    if ((num_i < imag_ptr_[i+1]) && (imag_ind_[num_i] == real_ind_[j]))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = real_ind_[j];
		Value(nnz) = T(real_data_[j], imag_data_[num_i]);
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = real_ind_[j];
		Value(nnz) = T(real_data_[j], 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr_[i+1]; j++)
          {
            IndCol(nnz) = imag_ind_[j];
            Value(nnz) = T(0, imag_data_[j]);
            nnz++;
          }
	
	Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }


  //! Conversion from RowSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    // ConvertMatrix_to_Coordinate provides sorted numbers (by rows here)
    Vector<Tint1> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (long i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);    
  }

  
  //! Conversion from ColSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();
    long* real_ptr = A.GetRealPtr();
    long* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();

    // first we count the number of non-zero entries for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (int i = 0; i < n; i++)
      {
	long num_i = imag_ptr[i];
	for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		// imaginary part alone
		NumNnzRow(imag_ind[num_i])++;		
		num_i++;
	      }
	    
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		// real and imaginary part are both present
		NumNnzRow(real_ind[j])++;		
		num_i++;
	      }
	    else
	      NumNnzRow(real_ind[j])++;
	  }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
	  NumNnzRow(imag_ind[j])++;	
      }
    
    // 'Ptr' array can be constructed
    Ptr.Reallocate(m+1);
    Ptr(0) = 0;	
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzRow(i);
    
    // total number of non-zero entries
    nnz = Ptr(m);
	
    // arrays 'IndCol' and 'Val' are filled already sorted by rows
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzRow.Zero();
    for (int i = 0; i < n; i++)
      {
	long num_i = imag_ptr[i]; int jcol = 0;
	for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		jcol = imag_ind[num_i];
		// imaginary part alone
		long num = Ptr(jcol) + NumNnzRow(jcol);
		IndCol(num) = i;
		Val(num) = T(0, imag_data[num_i]);
		NumNnzRow(jcol)++;
		num_i++;
	      }
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		jcol = real_ind[j];
		long num = Ptr(jcol) + NumNnzRow(jcol);
		// real and imaginary part are both present
		IndCol(num) = i;
		Val(num) = T(real_data[j], imag_data[num_i]);
		NumNnzRow(jcol)++; 
		num_i++;
	      }
	    else
	      {
		jcol = real_ind[j];
		long num = Ptr(jcol) + NumNnzRow(jcol);
		// real part alone
		IndCol(num) = i;
		Val(num) = T(real_data[j], 0);
		NumNnzRow(jcol)++; 
	      }
	  }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
	  {
	    jcol = imag_ind[j];
	    long num = Ptr(jcol) + NumNnzRow(jcol);
	    IndCol(num) = i;
	    Val(num) = T(0, imag_data[j]);
	    NumNnzRow(jcol)++;
	  }
      }    
  }


  //! Conversion from ColSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    Vector<Tint1, VectFull, Alloc3> IndRow;
   
    // this function sorts by column numbers
    // by symmetry, row numbers = col numbers, we invert the two arguments
    // to have sorted column numbers 
    ConvertMatrix_to_Coordinates(A, IndCol, IndRow, Val, 0, true);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }

  
  //! Conversion from ArrayRowComplexSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
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

    long real_nnz = A.GetRealDataSize();
    long imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    long nnz = 0;
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	int num_i = 0;
	int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		IndCol(nnz) = A.IndexImag(i, num_i);
		Value(nnz) = T(0, A.ValueImag(i, num_i));
		nnz++; num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = A.IndexReal(i, j);
		Value(nnz) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = A.IndexReal(i, j);
		Value(nnz) = T(A.ValueReal(i, j), 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
          {
            IndCol(nnz) = A.IndexImag(i, j);
            Value(nnz) = T(0, A.ValueImag(i, j));
            nnz++;
          }
	
	Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from ArrayColComplexSparse to CSR format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();

    // first we count the number of non-zero entries for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (int i = 0; i < n; i++)
      {
	int num_i = 0;
	int size_imag = A.GetImagColumnSize(i);
	for (int j = 0; j < A.GetRealColumnSize(i); j++)
	  {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		NumNnzRow(A.IndexImag(i, num_i))++;		
		num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		NumNnzRow(A.IndexReal(i, j))++;		
		num_i++;
	      }
	    else
	      NumNnzRow(A.IndexReal(i, j))++;
	  }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
	  NumNnzRow(A.IndexImag(i, j))++;
      }
    
    // 'Ptr' array can be constructed
    Ptr.Reallocate(m+1);
    Ptr(0) = 0;	
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzRow(i);
    
    // total number of non-zero entries
    nnz = Ptr(m);
	
    // arrays 'IndCol' and 'Val' are filled already sorted by rows
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzRow.Zero();
    for (int i = 0; i < n; i++)
      {
	int num_i = 0, jcol = 0;
	int size_imag = A.GetImagColumnSize(i);
	for (int j = 0; j < A.GetRealColumnSize(i); j++)
	  {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		jcol = A.IndexImag(i, num_i);
		// imaginary part alone
		long num = Ptr(jcol) + NumNnzRow(jcol);
		IndCol(num) = i;
		Val(num) = T(0, A.ValueImag(i, num_i));
		NumNnzRow(jcol)++;
		num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		jcol = A.IndexReal(i, j);
		long num = Ptr(jcol) + NumNnzRow(jcol);
		// real and imaginary part are both present
		IndCol(num) = i;
		Val(num) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		NumNnzRow(jcol)++; 
		num_i++;
	      }
	    else
	      {
		jcol = A.IndexReal(i, j);
		long num = Ptr(jcol) + NumNnzRow(jcol);
		// real part alone
		IndCol(num) = i;
		Val(num) = T(A.ValueReal(i, j), 0);
		NumNnzRow(jcol)++; 
	      }
	  }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
	  {
	    jcol = A.IndexImag(i, j);
	    long num = Ptr(jcol) + NumNnzRow(jcol);
	    IndCol(num) = i;
	    Val(num) = T(0, A.ValueImag(i, j));
	    NumNnzRow(jcol)++;
	  }
      }
  }

  
  //! Conversion from ArrayRowSymComplexSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
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

    long real_nnz = A.GetRealDataSize();
    long imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    long nnz = 0;
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	int num_i = 0;
	int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		IndCol(nnz) = A.IndexImag(i, num_i);
		Value(nnz) = T(0, A.ValueImag(i, num_i));
		nnz++; num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = A.IndexReal(i, j);
		Value(nnz) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = A.IndexReal(i, j);
		Value(nnz) = T(A.ValueReal(i, j), 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
          {
            IndCol(nnz) = A.IndexImag(i, j);
            Value(nnz) = T(0, A.ValueImag(i, j));
            nnz++;
          }
	
	Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }


  //! Conversion from ArrayRowSymComplexSparse to CSR
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    // ConvertMatrix_to_Coordinate provides sorted numbers (by rows here)
    Vector<Tint1> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (long i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);    
  }
  
    
  //! Conversion from ArrayColSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();

    // first we count the number of non-zero entries for each row
    Vector<int> NumNnzRow(m);
    NumNnzRow.Zero();
    for (int i = 0; i < n; i++)
      {
	int num_i = 0;
	int size_imag = A.GetImagColumnSize(i);
	for (int j = 0; j < A.GetRealColumnSize(i); j++)
	  {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		NumNnzRow(A.IndexImag(i, num_i))++;		
		num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		NumNnzRow(A.IndexReal(i, j))++;		
		num_i++;
	      }
	    else
	      NumNnzRow(A.IndexReal(i, j))++;
	  }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
	  NumNnzRow(A.IndexImag(i, j))++;
      }
    
    // 'Ptr' array can be constructed
    Ptr.Reallocate(m+1);
    Ptr(0) = 0;	
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzRow(i);
    
    // total number of non-zero entries
    nnz = Ptr(m);
	
    // arrays 'IndCol' and 'Val' are filled already sorted by rows
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzRow.Zero();
    for (int i = 0; i < n; i++)
      {
	int num_i = 0, jcol = 0;
	int size_imag = A.GetImagColumnSize(i);
	for (int j = 0; j < A.GetRealColumnSize(i); j++)
	  {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		jcol = A.IndexImag(i, num_i);
		// imaginary part alone
		long num = Ptr(jcol) + NumNnzRow(jcol);
		IndCol(num) = i;
		Val(num) = T(0, A.ValueImag(i, num_i));
		NumNnzRow(jcol)++;
		num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		jcol = A.IndexReal(i, j);
		long num = Ptr(jcol) + NumNnzRow(jcol);
		// real and imaginary part are both present
		IndCol(num) = i;
		Val(num) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		NumNnzRow(jcol)++; 
		num_i++;
	      }
	    else
	      {
		jcol = A.IndexReal(i, j);
		long num = Ptr(jcol) + NumNnzRow(jcol);
		// real part alone
		IndCol(num) = i;
		Val(num) = T(A.ValueReal(i, j), 0);
		NumNnzRow(jcol)++; 
	      }
	  }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
	  {
	    jcol = A.IndexImag(i, j);
	    long num = Ptr(jcol) + NumNnzRow(jcol);
	    IndCol(num) = i;
	    Val(num) = T(0, A.ValueImag(i, j));
	    NumNnzRow(jcol)++;
	  }
      }
  }

  
  //! Conversion from ArrayColSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    Vector<Tint1, VectFull, Alloc3> IndRow;
   
    // this function sorts by column numbers
    // by symmetry, row numbers = col numbers, we invert the two arguments
    // to have sorted column numbers 
    ConvertMatrix_to_Coordinates(A, IndCol, IndRow, Val, 0, true);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }

  
  //! Conversion from RowComplexSparse to CSC format
  /*!
    if sym_pat is equal to true, the pattern is symmetrized
    by adding artificial null entries
   */
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    if (sym_pat)
      throw Undefined("ConvertToCSC", "Symmetrization of pattern not"
                      " implemented for RowComplexSparse storage");

    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();
    long* real_ptr = A.GetRealPtr();
    long* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();

    // first we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(n);
    NumNnzCol.Zero();
    for (int i = 0; i < m; i++)
      {
	long num_i = imag_ptr[i];
	for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		// imaginary part alone
		NumNnzCol(imag_ind[num_i])++;		
		num_i++;
	      }
	    
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		// real and imaginary part are both present
		NumNnzCol(real_ind[j])++;		
		num_i++;
	      }
	    else
	      NumNnzCol(real_ind[j])++;
	  }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
	  NumNnzCol(imag_ind[j])++;	
      }
    
    // 'Ptr' array can be constructed
    Ptr.Reallocate(n+1);
    Ptr(0) = 0;	
    for (int i = 0; i < n; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);
    
    // total number of non-zero entries
    nnz = Ptr(n);
	
    // arrays 'IndRow' and 'Val' are filled already sorted by columns
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzCol.Zero();
    for (int i = 0; i < m; i++)
      {
	long num_i = imag_ptr[i]; int jcol = 0;
	for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		jcol = imag_ind[num_i];
		// imaginary part alone
		long num = Ptr(jcol) + NumNnzCol(jcol);
		IndRow(num) = i;
		Val(num) = T(0, imag_data[num_i]);
		NumNnzCol(jcol)++;
		num_i++;
	      }
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		jcol = real_ind[j];
		long num = Ptr(jcol) + NumNnzCol(jcol);
		// real and imaginary part are both present
		IndRow(num) = i;
		Val(num) = T(real_data[j], imag_data[num_i]);
		NumNnzCol(jcol)++; 
		num_i++;
	      }
	    else
	      {
		jcol = real_ind[j];
		long num = Ptr(jcol) + NumNnzCol(jcol);
		// real part alone
		IndRow(num) = i;
		Val(num) = T(real_data[j], 0);
		NumNnzCol(jcol)++; 
	      }
	  }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
	  {
	    jcol = imag_ind[j];
	    long num = Ptr(jcol) + NumNnzCol(jcol);
	    IndRow(num) = i;
	    Val(num) = T(0, imag_data[j]);
	    NumNnzCol(jcol)++;
	  }
      }
  }
  
  
  //! Conversion from ColComplexSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat)
  {
    if (sym_pat)
      throw Undefined("ConvertToCSC", "Symmetrization of pattern not"
                      " implemented for ColComplexSparse storage");

    int n = A.GetN();
    if (n <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    long* real_ptr_ = A.GetRealPtr();
    int* real_ind_ = A.GetRealInd();
    Treal* real_data_ = A.GetRealData();

    long* imag_ptr_ = A.GetImagPtr();
    int* imag_ind_ = A.GetImagInd();
    Treal* imag_data_ = A.GetImagData();
    
    long real_nnz = A.GetRealDataSize();
    long imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(n+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    long nnz = 0;
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      {
	long num_i = imag_ptr_[i];
        for (long j = real_ptr_[i]; j < real_ptr_[i+1]; j++)
          {
	    while ((num_i < imag_ptr_[i+1]) && (imag_ind_[num_i] < real_ind_[j]))
	      {
		// imaginary part alone
		IndCol(nnz) = imag_ind_[num_i];
		Value(nnz) = T(0, imag_data_[num_i]);
		nnz++; num_i++;
	      }
	    
	    if ((num_i < imag_ptr_[i+1]) && (imag_ind_[num_i] == real_ind_[j]))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = real_ind_[j];
		Value(nnz) = T(real_data_[j], imag_data_[num_i]);
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = real_ind_[j];
		Value(nnz) = T(real_data_[j], 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr_[i+1]; j++)
          {
            IndCol(nnz) = imag_ind_[j];
            Value(nnz) = T(0, imag_data_[j]);
            nnz++;
          }
	
	Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from RowSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();
    long* real_ptr = A.GetRealPtr();
    long* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();

    // first we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(m);
    NumNnzCol.Zero();
    for (int i = 0; i < n; i++)
      {
	long num_i = imag_ptr[i];
	for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		// imaginary part alone
		NumNnzCol(imag_ind[num_i])++;		
		num_i++;
	      }
	    
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		// real and imaginary part are both present
		NumNnzCol(real_ind[j])++;		
		num_i++;
	      }
	    else
	      NumNnzCol(real_ind[j])++;
	  }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
	  NumNnzCol(imag_ind[j])++;	
      }
    
    // 'Ptr' array can be constructed
    Ptr.Reallocate(m+1);
    Ptr(0) = 0;	
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);
    
    // total number of non-zero entries
    nnz = Ptr(m);
	
    // arrays 'IndRow' and 'Val' are filled already sorted by columns
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzCol.Zero();
    for (int i = 0; i < n; i++)
      {
	long num_i = imag_ptr[i]; int jcol = 0;
	for (long j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
	    while ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] < real_ind[j]))
	      {
		jcol = imag_ind[num_i];
		// imaginary part alone
		long num = Ptr(jcol) + NumNnzCol(jcol);
		IndRow(num) = i;
		Val(num) = T(0, imag_data[num_i]);
		NumNnzCol(jcol)++;
		num_i++;
	      }
	    
	    if ((num_i < imag_ptr[i+1]) && (imag_ind[num_i] == real_ind[j]))
	      {
		jcol = real_ind[j];
		long num = Ptr(jcol) + NumNnzCol(jcol);
		// real and imaginary part are both present
		IndRow(num) = i;
		Val(num) = T(real_data[j], imag_data[num_i]);
		NumNnzCol(jcol)++; 
		num_i++;
	      }
	    else
	      {
		jcol = real_ind[j];
		long num = Ptr(jcol) + NumNnzCol(jcol);
		// real part alone
		IndRow(num) = i;
		Val(num) = T(real_data[j], 0);
		NumNnzCol(jcol)++; 
	      }
	  }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr[i+1]; j++)
	  {
	    jcol = imag_ind[j];
	    long num = Ptr(jcol) + NumNnzCol(jcol);
	    IndRow(num) = i;
	    Val(num) = T(0, imag_data[j]);
	    NumNnzCol(jcol)++;
	  }
      }    
  }
  

  //! Conversion from RowSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    Vector<Tint1, VectFull, Alloc3> IndCol;
   
    // this function sorts by row numbers
    // by symmetry, row numbers = col numbers, we invert the two arguments
    // to have sorted column numbers 
    ConvertMatrix_to_Coordinates(A, IndCol, IndRow, Val, 0, true);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }

  
  //! Conversion from ColSymComplexSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat)
  {
    int m = A.GetM();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    long* real_ptr_ = A.GetRealPtr();
    int* real_ind_ = A.GetRealInd();
    Treal* real_data_ = A.GetRealData();

    long* imag_ptr_ = A.GetImagPtr();
    int* imag_ind_ = A.GetImagInd();
    Treal* imag_data_ = A.GetImagData();
    
    long real_nnz = A.GetRealDataSize();
    long imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    long nnz = 0;
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	long num_i = imag_ptr_[i];
        for (long j = real_ptr_[i]; j < real_ptr_[i+1]; j++)
          {
	    while ((num_i < imag_ptr_[i+1]) && (imag_ind_[num_i] < real_ind_[j]))
	      {
		// imaginary part alone
		IndCol(nnz) = imag_ind_[num_i];
		Value(nnz) = T(0, imag_data_[num_i]);
		nnz++; num_i++;
	      }
	    
	    if ((num_i < imag_ptr_[i+1]) && (imag_ind_[num_i] == real_ind_[j]))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = real_ind_[j];
		Value(nnz) = T(real_data_[j], imag_data_[num_i]);
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = real_ind_[j];
		Value(nnz) = T(real_data_[j], 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (long j = num_i; j < imag_ptr_[i+1]; j++)
          {
            IndCol(nnz) = imag_ind_[j];
            Value(nnz) = T(0, imag_data_[j]);
            nnz++;
          }
	
	Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from ColSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    // ConvertMatrix_to_Coordinate provides sorted numbers (by columns here)
    Vector<Tint1> IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per columnx.
    for (long i = 0; i < IndRow.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);    
    
  }

  
  //! Conversion from ArrayRowComplexSparse to CSC format
  /*!
    if sym_pat is equal to true, the pattern is symmetrized
    by adding artificial null entries
   */
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    if (sym_pat)
      throw Undefined("ConvertToCSC", "Symmetrization of pattern not"
                      "implemented for RowComplexSparse storage");
    
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();

    // first we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(n);
    NumNnzCol.Zero();
    for (int i = 0; i < m; i++)
      {
	int num_i = 0;
	int size_imag = A.GetImagRowSize(i);
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		NumNnzCol(A.IndexImag(i, num_i))++;		
		num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		NumNnzCol(A.IndexReal(i, j))++;		
		num_i++;
	      }
	    else
	      NumNnzCol(A.IndexReal(i, j))++;
	  }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
	  NumNnzCol(A.IndexImag(i, j))++;
      }
    
    // 'Ptr' array can be constructed
    Ptr.Reallocate(n+1);
    Ptr(0) = 0;	
    for (int i = 0; i < n; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);
    
    // total number of non-zero entries
    nnz = Ptr(n);
	
    // arrays 'IndRow' and 'Val' are filled already sorted by columns
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzCol.Zero();
    for (int i = 0; i < m; i++)
      {
	int num_i = 0, jcol = 0;
	int size_imag = A.GetImagRowSize(i);
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		jcol = A.IndexImag(i, num_i);
		// imaginary part alone
		long num = Ptr(jcol) + NumNnzCol(jcol);
		IndRow(num) = i;
		Val(num) = T(0, A.ValueImag(i, num_i));
		NumNnzCol(jcol)++;
		num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		jcol = A.IndexReal(i, j);
		long num = Ptr(jcol) + NumNnzCol(jcol);
		// real and imaginary part are both present
		IndRow(num) = i;
		Val(num) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		NumNnzCol(jcol)++; 
		num_i++;
	      }
	    else
	      {
		jcol = A.IndexReal(i, j);
		long num = Ptr(jcol) + NumNnzCol(jcol);
		// real part alone
		IndRow(num) = i;
		Val(num) = T(A.ValueReal(i, j), 0);
		NumNnzCol(jcol)++; 
	      }
	  }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
	  {
	    jcol = A.IndexImag(i, j);
	    long num = Ptr(jcol) + NumNnzCol(jcol);
	    IndRow(num) = i;
	    Val(num) = T(0, A.ValueImag(i, j));
	    NumNnzCol(jcol)++;
	  }
      }
  }
  
  
  //! Conversion from ArrayColComplexSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat)
  {
    if (sym_pat)
      throw Undefined("ConvertToCSC", "Symmetrization of pattern not"
                      "implemented for ColComplexSparse storage");

    int n = A.GetN();
    if (n <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    long real_nnz = A.GetRealDataSize();
    long imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(n+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    long nnz = 0;
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      {
	int num_i = 0;
	int size_imag = A.GetImagColumnSize(i);
        for (int j = 0; j < A.GetRealColumnSize(i); j++)
          {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		IndCol(nnz) = A.IndexImag(i, num_i);
		Value(nnz) = T(0, A.ValueImag(i, num_i));
		nnz++; num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = A.IndexReal(i, j);
		Value(nnz) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = A.IndexReal(i, j);
		Value(nnz) = T(A.ValueReal(i, j), 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
          {
            IndCol(nnz) = A.IndexImag(i, j);
            Value(nnz) = T(0, A.ValueImag(i, j));
            nnz++;
          }
	
	Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from ArrayRowSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();

    // first we count the number of non-zero entries for each column
    Vector<int> NumNnzCol(m);
    NumNnzCol.Zero();
    for (int i = 0; i < n; i++)
      {
	int num_i = 0;
	int size_imag = A.GetImagRowSize(i);
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		NumNnzCol(A.IndexImag(i, num_i))++;		
		num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		NumNnzCol(A.IndexReal(i, j))++;		
		num_i++;
	      }
	    else
	      NumNnzCol(A.IndexReal(i, j))++;
	  }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
	  NumNnzCol(A.IndexImag(i, j))++;
      }
    
    // 'Ptr' array can be constructed
    Ptr.Reallocate(m+1);
    Ptr(0) = 0;	
    for (int i = 0; i < m; i++)
      Ptr(i+1) = Ptr(i) + NumNnzCol(i);
    
    // total number of non-zero entries
    nnz = Ptr(m);
	
    // arrays 'IndRow' and 'Val' are filled already sorted by rows
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    NumNnzCol.Zero();
    for (int i = 0; i < n; i++)
      {
	int num_i = 0, jcol = 0;
	int size_imag = A.GetImagRowSize(i);
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		jcol = A.IndexImag(i, num_i);
		// imaginary part alone
		long num = Ptr(jcol) + NumNnzCol(jcol);
		IndRow(num) = i;
		Val(num) = T(0, A.ValueImag(i, num_i));
		NumNnzCol(jcol)++;
		num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		jcol = A.IndexReal(i, j);
		long num = Ptr(jcol) + NumNnzCol(jcol);
		// real and imaginary part are both present
		IndRow(num) = i;
		Val(num) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		NumNnzCol(jcol)++; 
		num_i++;
	      }
	    else
	      {
		jcol = A.IndexReal(i, j);
		long num = Ptr(jcol) + NumNnzCol(jcol);
		// real part alone
		IndRow(num) = i;
		Val(num) = T(A.ValueReal(i, j), 0);
		NumNnzCol(jcol)++; 
	      }
	  }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
	  {
	    jcol = A.IndexImag(i, j);
	    long num = Ptr(jcol) + NumNnzCol(jcol);
	    IndRow(num) = i;
	    Val(num) = T(0, A.ValueImag(i, j));
	    NumNnzCol(jcol)++;
	  }
      }
  }
  

  //! Conversion from ArrayRowSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    Vector<Tint1, VectFull, Alloc3> IndCol;
   
    // this function sorts by row numbers
    // by symmetry, row numbers = col numbers, we invert the two arguments
    // to have sorted column numbers 
    ConvertMatrix_to_Coordinates(A, IndCol, IndRow, Val, 0, true);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (long i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
    
  }

  
  //! Conversion from ArrayColSymComplexSparse to CSC
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat)
  {
    int n = A.GetN();
    if (n <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    long real_nnz = A.GetRealDataSize();
    long imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(n+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    long nnz = 0;
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      {
	int num_i = 0;
	int size_imag = A.GetImagColumnSize(i);
        for (int j = 0; j < A.GetRealColumnSize(i); j++)
          {
	    while ((num_i < size_imag) && (A.IndexImag(i, num_i) < A.IndexReal(i, j)))
	      {
		// imaginary part alone
		IndCol(nnz) = A.IndexImag(i, num_i);
		Value(nnz) = T(0, A.ValueImag(i, num_i));
		nnz++; num_i++;
	      }
	    
	    if ((num_i < size_imag) && (A.IndexImag(i, num_i) == A.IndexReal(i, j)))
	      {
		// real and imaginary part are both present
		IndCol(nnz) = A.IndexReal(i, j);
		Value(nnz) = T(A.ValueReal(i, j), A.ValueImag(i, num_i));
		num_i++; nnz++;
	      }
	    else
	      {
		// real part alone
		IndCol(nnz) = A.IndexReal(i, j);
		Value(nnz) = T(A.ValueReal(i, j), 0);
		nnz++;
	      }
          }
	
	// last values of imaginary part
	for (int j = num_i; j < size_imag; j++)
          {
            IndCol(nnz) = A.IndexImag(i, j);
            Value(nnz) = T(0, A.ValueImag(i, j));
            nnz++;
          }
	
	Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }


  //! Conversion from ArrayColSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    // ConvertMatrix_to_Coordinate provides sorted numbers (by columns here)
    Vector<Tint1> IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (long i = 0; i < IndRow.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);    
  }


  //! Conversion from ArrayRowSymComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayRowSymComplexSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayRowComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayColSymComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayColSymComplexSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayColComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from RowSymComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from RowSymComplexSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from RowComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ColSymComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ColSymComplexSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ColComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<long> Ptr; Vector<int> IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayRowComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayRowSymComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayRowSymComplexSparse to ColSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSymSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayColComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayColSymComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayColSymComplexSparse to ColSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSymSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from RowComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from RowSymComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }  
  

  //! Conversion from RowSymComplexSparse to ColSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSymSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }  


  //! Conversion from ColComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }
  
  
  //! Conversion from ColSymComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }  
  

  //! Conversion from ColSymComplexSparse to ColSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSymSparse, Allocator1>& B)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }  

  
  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowComplexSparse, Allocator>& A,
		  Matrix<T, Prop, RowComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A,
		  Matrix<T, Prop, RowSymComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColComplexSparse, Allocator>& A,
		  Matrix<T, Prop, ColComplexSparse, Allocator>& B)
  {
    B = A;
  }
  

  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSymComplexSparse, Allocator>& A,
		  Matrix<T, Prop, ColSymComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColComplexSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayColComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>& B)
  {
    B = A;
  }

  
  //! Conversion from ColSparse to RowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
	     Matrix<T1, Prop1, RowComplexSparse, Allocator1>& B)
  {
    int m = A.GetM(), n = A.GetN();
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    General sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);
    
    // counting number of non-zero elements
    typedef typename ClassComplexType<T1>::Treal Treal;
    long imag_nnz = 0, real_nnz = 0;
    Treal zero(0);
    for (long i = 0; i < Value.GetM(); i++)
      {
        if (real(Value(i)) != zero)
          real_nnz++;
        
        if (imag(Value(i)) != zero)
          imag_nnz++;
      }
    
    Vector<long> Ptr_real(m+1), Ptr_imag(m+1);
    Vector<int> Ind_real(real_nnz), Ind_imag(imag_nnz);
    Vector<Treal, VectFull, Allocator1> Val_real(real_nnz), Val_imag(imag_nnz);
    
    // filling arrays Ind_real, Ind_imag, Val_real, Val_imag
    Ptr_real(0) = 0; Ptr_imag(0) = 0;
    real_nnz = 0; imag_nnz = 0;
    for (int i = 0; i < m; i++)
      {
        for (long j = Ptr(i); j < Ptr(i+1); j++)
          {
            if (real(Value(j)) != zero)
              {
                Ind_real(real_nnz) = Ind(j);
                Val_real(real_nnz) = real(Value(j));
                real_nnz++;
              }
            
            if (imag(Value(j)) != zero)
              {
                Ind_imag(imag_nnz) = Ind(j);
                Val_imag(imag_nnz) = imag(Value(j));
                imag_nnz++;
              }
          }
        
        Ptr_real(i+1) = real_nnz;
        Ptr_imag(i+1) = imag_nnz;
      }

    
    // creating the matrix
    B.SetData(m, n, Val_real, Ptr_real, Ind_real,
              Val_imag, Ptr_imag, Ind_imag);
  }
  
  
  //! conversion from ColSparse to ColComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColComplexSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    long nnz = A.GetDataSize();
    
    // counting number of non-zero elements
    long imag_nnz = 0, real_nnz = 0;
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    for (long i = 0; i < nnz; i++)
      {
        if (real(data_[i]) != zero)
          real_nnz++;
        
        if (imag(data_[i]) != zero)
          imag_nnz++;
      }
    
    Vector<long> Ptr_real(n+1), Ptr_imag(n+1);
    Vector<int> Ind_real(real_nnz), Ind_imag(imag_nnz);
    
    Vector<Treal, VectFull, Allocator1> Val_real(real_nnz), Val_imag(imag_nnz);
    
    // filling arrays Ind_real, Ind_imag, Val_real, Val_imag
    Ptr_real(0) = 0; Ptr_imag(0) = 0;
    real_nnz = 0; imag_nnz = 0;
    for (int i = 0; i < n; i++)
      {
        for (long j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              {
                Ind_real(real_nnz) = ind_[j];
                Val_real(real_nnz) = real(data_[j]);
                real_nnz++;
              }
            
            if (imag(data_[j]) != zero)
              {
                Ind_imag(imag_nnz) = ind_[j];
                Val_imag(imag_nnz) = imag(data_[j]);
                imag_nnz++;
              }
          }
        
        Ptr_real(i+1) = real_nnz;
        Ptr_imag(i+1) = imag_nnz;
      }

    
    // creating the matrix
    B.SetData(m, n, Val_real, Ptr_real, Ind_real,
              Val_imag, Ptr_imag, Ind_imag);
  }
  
  
  //! Conversion from RowSymSparse to ColSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& B)
  {
    int m = A.GetM(), n = A.GetN();
    Vector<long> Ptr; Vector<int> Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Value);
    
    // counting number of non-zero elements
    long imag_nnz = 0, real_nnz = 0;
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    for (long i = 0; i < Value.GetM(); i++)
      {
        if (real(Value(i)) != zero)
          real_nnz++;
        
        if (imag(Value(i)) != zero)
          imag_nnz++;
      }
    
    Vector<long> Ptr_real(m+1), Ptr_imag(m+1);
    Vector<int> Ind_real(real_nnz), Ind_imag(imag_nnz);
    
    Vector<Treal, VectFull, Allocator1> Val_real(real_nnz), Val_imag(imag_nnz);
    
    // filling arrays Ind_real, Ind_imag, Val_real, Val_imag
    Ptr_real(0) = 0; Ptr_imag(0) = 0;
    real_nnz = 0; imag_nnz = 0;
    for (int i = 0; i < m; i++)
      {
        for (long j = Ptr(i); j < Ptr(i+1); j++)
          {
            if (real(Value(j)) != zero)
              {
                Ind_real(real_nnz) = Ind(j);
                Val_real(real_nnz) = real(Value(j));
                real_nnz++;
              }
            
            if (imag(Value(j)) != zero)
              {
                Ind_imag(imag_nnz) = Ind(j);
                Val_imag(imag_nnz) = imag(Value(j));
                imag_nnz++;
              }
          }
        
        Ptr_real(i+1) = real_nnz;
        Ptr_imag(i+1) = imag_nnz;
      }

    
    // creating the matrix
    B.SetData(m, n, Val_real, Ptr_real, Ind_real,
              Val_imag, Ptr_imag, Ind_imag);
  }
  
  
  //! conversion from RowSymSparse to RowSymComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
	     Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    long nnz = A.GetDataSize();
    
    // counting number of non-zero elements
    long imag_nnz = 0, real_nnz = 0;
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    for (long i = 0; i < nnz; i++)
      {
        if (real(data_[i]) != zero)
          real_nnz++;
        
        if (imag(data_[i]) != zero)
          imag_nnz++;
      }
    
    Vector<long> Ptr_real(n+1), Ptr_imag(n+1);
    Vector<int> Ind_real(real_nnz), Ind_imag(imag_nnz);
    
    Vector<Treal, VectFull, Allocator1> Val_real(real_nnz), Val_imag(imag_nnz);
    
    // filling arrays Ind_real, Ind_imag, Val_real, Val_imag
    Ptr_real(0) = 0; Ptr_imag(0) = 0;
    real_nnz = 0; imag_nnz = 0;
    for (int i = 0; i < n; i++)
      {
        for (long j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              {
                Ind_real(real_nnz) = ind_[j];
                Val_real(real_nnz) = real(data_[j]);
                real_nnz++;
              }
            
            if (imag(data_[j]) != zero)
              {
                Ind_imag(imag_nnz) = ind_[j];
                Val_imag(imag_nnz) = imag(data_[j]);
                imag_nnz++;
              }
          }
        
        Ptr_real(i+1) = real_nnz;
        Ptr_imag(i+1) = imag_nnz;
      }

    
    // creating the matrix
    B.SetData(m, n, Val_real, Ptr_real, Ind_real,
              Val_imag, Ptr_imag, Ind_imag);
  }

  
  //! Conversion from ColSparse to ArrayRowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& B)
  {
    int m = A.GetM(), n = A.GetN();
    Matrix<T0, Prop0, ArrayRowSparse, Allocator0> Ar;
    CopyMatrix(A, Ar);
    
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_real = 0, size_imag = 0;
        for (int j = 0; j < Ar.GetRowSize(i); j++)
          {
            if (real(Ar.Value(i, j)) != zero)
              size_real++;
            
            if (imag(Ar.Value(i, j)) != zero)
              size_imag++;
          }
        
        B.ReallocateRealRow(i, size_real);
        B.ReallocateImagRow(i, size_imag);
        size_real = 0;
        size_imag = 0;
        for (int j = 0; j < Ar.GetRowSize(i); j++)
          {
            if (real(Ar.Value(i, j)) != zero)
              {
                B.IndexReal(i, size_real) = Ar.Index(i, j);
                B.ValueReal(i, size_real) = real(Ar.Value(i, j));
                size_real++;
              }
            
            if (imag(Ar.Value(i, j)) != zero)
              {
                B.IndexImag(i, size_imag) = Ar.Index(i, j);
                B.ValueImag(i, size_imag) = imag(Ar.Value(i, j));
                size_imag++;
              }
          }
        
        Ar.ClearRow(i);
      }
  }
  
  
  //! conversion from ColSparse to ArrayColComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ArrayColComplexSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    B.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
        int size_real = 0, size_imag = 0;
        for (long j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              size_real++;
            
            if (imag(data_[j]) != zero)
              size_imag++;
          }
        
        B.ReallocateRealColumn(i, size_real);
        B.ReallocateImagColumn(i, size_imag);
        size_real = 0;
        size_imag = 0;
        for (long j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              {
                B.IndexReal(i, size_real) = ind_[j];
                B.ValueReal(i, size_real) = real(data_[j]);
                size_real++;
              }
            
            if (imag(data_[j]) != zero)
              {
                B.IndexImag(i, size_imag) = ind_[j];
                B.ValueImag(i, size_imag) = imag(data_[j]);
                size_imag++;
              }
          }
      }
  }
  
  
  //! Conversion from RowSymSparse to ArrayColSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ArrayColSymComplexSparse, Allocator1>& B)
  {
    int m = A.GetM(), n = A.GetN();
    Matrix<T0, Prop0, ArrayColSymSparse, Allocator0> Ar;
    CopyMatrix(A, Ar);

    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_real = 0, size_imag = 0;
        for (int j = 0; j < Ar.GetColumnSize(i); j++)
          {
            if (real(Ar.Value(i, j)) != zero)
              size_real++;
            
            if (imag(Ar.Value(i, j)) != zero)
              size_imag++;
          }
        
        B.ReallocateRealColumn(i, size_real);
        B.ReallocateImagColumn(i, size_imag);
        size_real = 0;
        size_imag = 0;
        for (int j = 0; j < Ar.GetColumnSize(i); j++)
          {
            if (real(Ar.Value(i, j)) != zero)
              {
                B.IndexReal(i, size_real) = Ar.Index(i, j);
                B.ValueReal(i, size_real) = real(Ar.Value(i, j));
                size_real++;
              }
            
            if (imag(Ar.Value(i, j)) != zero)
              {
                B.IndexImag(i, size_imag) = Ar.Index(i, j);
                B.ValueImag(i, size_imag) = imag(Ar.Value(i, j));
                size_imag++;
              }
          }
        
        Ar.ClearColumn(i);
      }
  }
  
  
  //! conversion from RowSymSparse to ArrayRowSymComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ArrayRowSymComplexSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    B.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
        int size_real = 0, size_imag = 0;
        for (long j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              size_real++;
            
            if (imag(data_[j]) != zero)
              size_imag++;
          }
        
        B.ReallocateRealRow(i, size_real);
        B.ReallocateImagRow(i, size_imag);
        size_real = 0;
        size_imag = 0;
        for (long j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              {
                B.IndexReal(i, size_real) = ind_[j];
                B.ValueReal(i, size_real) = real(data_[j]);
                size_real++;
              }
            
            if (imag(data_[j]) != zero)
              {
                B.IndexImag(i, size_imag) = ind_[j];
                B.ValueImag(i, size_imag) = imag(data_[j]);
                size_imag++;
              }
          }       
      }
  }

  
  //! Conversion from ArrayRowComplexSparse to RowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    long nnz_real = mat_array.GetRealDataSize();
    long nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    
    // Allocation of arrays for CSR.
    typedef typename ClassComplexType<T1>::Treal Treal;
    Vector<Treal, VectFull, Allocator1> Val_real(nnz_real), Val_imag(nnz_imag);
    Vector<long> IndRow_real(m + 1);
    Vector<long> IndRow_imag(m + 1);
    Vector<int> IndCol_real(nnz_real);
    Vector<int> IndCol_imag(nnz_imag);

    long ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }

	IndRow_real(i + 1) = ind_real;
	for (k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }

	IndRow_imag(i + 1) = ind_imag;
      }

    mat_csr.SetData(m, n, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }


  //! Conversion from ArrayRowSymComplexSparse to RowSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0,
	     ArrayRowSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    long nnz_real = mat_array.GetRealDataSize();
    long nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();

    // Allocation of arrays for CSR.
    typedef typename ClassComplexType<T1>::Treal Treal;
    Vector<Treal, VectFull, Allocator1> Val_real(nnz_real), Val_imag(nnz_imag);
    Vector<long> IndRow_real(m + 1);
    Vector<long> IndRow_imag(m + 1);
    Vector<int> IndCol_real(nnz_real);
    Vector<int> IndCol_imag(nnz_imag);

    long ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }

	IndRow_real(i + 1) = ind_real;
	for (int k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }

	IndRow_imag(i + 1) = ind_imag;
      }

    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }
  
  
  //! Conversion from ArrayRowSymComplexSparse to ArrayRowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& mat_csr)
  {
    Vector<int> IndRow, IndCol;
    Vector<T1, VectFull, Allocator1> Value;
    
    // ConvertMatrix_to_Coordinates already sorts by rows
    ConvertMatrix_to_Coordinates(mat_array, IndRow, IndCol, Value, 0, true);

    int m = mat_array.GetM();
    int n = mat_array.GetN();    
    mat_csr.Reallocate(m, n);
    long k = 0;
    for (int i = 0; i < m; i++)
      {
        long k0 = k;
        while ((k < IndRow.GetM()) && (IndRow(k) <= i))
          k++;

        int size_row = k - k0;
        if (size_row > 0)
          {
            mat_csr.ReallocateRow(i, size_row);
            for (int j = 0; j < size_row; j++)
              {
                mat_csr.Index(i, j) = IndCol(k0 + j);
                mat_csr.Value(i, j) = Value(k0 + j);
              }
          }
        else
          mat_csr.ClearRow(i);
      }
  }



  //! Conversion from ArrayRowSymComplexSparse to ArrayRowSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    Vector<T0> value(n);
    Vector<int> col_num(n);
    for (int i = 0; i < m; i++)
      {
        int size_real = A.GetRealRowSize(i);
        int size_imag = A.GetImagRowSize(i);
        int jr = 0, ji = 0;
        int size_row = 0;
        while (jr < size_real)
          {
            int col = A.IndexReal(i, jr);
            while ((ji < size_imag) && (A.IndexImag(i, ji) < col))
              {
                col_num(size_row) = A.IndexImag(i, ji);
                value(size_row) = T0(0, A.ValueImag(i, ji));
                ji++; size_row++;
              }
            
            if ((ji < size_imag) && (A.IndexImag(i, ji) == col))
              {
                col_num(size_row) = col;
                value(size_row) = T0(A.ValueReal(i, jr),
				     A.ValueImag(i, ji));
                ji++; size_row++;
              }
            else
              {
                col_num(size_row) = col;
                value(size_row) = T0(A.ValueReal(i, jr), 0);
                size_row++;
              }
            
            jr++;
          }

        while (ji < size_imag)
          {
            col_num(size_row) = A.IndexImag(i, ji);
            value(size_row) = T0(0, A.ValueImag(i, ji));
            ji++; size_row++;
          }
        
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = col_num(j);
            B.Value(i, j) = value(j);
          }
      }
  }


  //! Conversion from ArrayRowComplexSparse to ArrayRowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    Vector<T0> value(n);
    Vector<int> col_num(n);
    for (int i = 0; i < m; i++)
      {
        int size_real = A.GetRealRowSize(i);
        int size_imag = A.GetImagRowSize(i);
        int jr = 0, ji = 0;
        int size_row = 0;
        while (jr < size_real)
          {
            int col = A.IndexReal(i, jr);
            while ((ji < size_imag) && (A.IndexImag(i, ji) < col))
              {
                col_num(size_row) = A.IndexImag(i, ji);
                value(size_row) = T0(0, A.ValueImag(i, ji));
                ji++; size_row++;
              }
            
            if ((ji < size_imag) && (A.IndexImag(i, ji) == col))
              {
                col_num(size_row) = col;
                value(size_row) = T0(A.ValueReal(i, jr),
				     A.ValueImag(i, ji));
                ji++; size_row++;
              }
            else
              {
                col_num(size_row) = col;
                value(size_row) = T0(A.ValueReal(i, jr), 0);
                size_row++;
              }
            
            jr++;
          }
        
        while (ji < size_imag)
          {
            col_num(size_row) = A.IndexImag(i, ji);
            value(size_row) = T0(0, A.ValueImag(i, ji));
            ji++; size_row++;
          }
        
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = col_num(j);
            B.Value(i, j) = value(j);
          }
      }
  }


  //! Conversion from RowSymComplexSparse to ArrayRowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
	     Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& mat_csr)
  {
    Vector<int> IndRow, IndCol;
    Vector<T1, VectFull, Allocator1> Value;
    
    // ConvertMatrix_to_Coordinates already sorts by rows
    ConvertMatrix_to_Coordinates(mat_array, IndRow, IndCol, Value, 0, true);

    int m = mat_array.GetM();
    int n = mat_array.GetN();    
    mat_csr.Reallocate(m, n);
    long k = 0;
    for (int i = 0; i < m; i++)
      {
        long k0 = k;
        while ((k < IndRow.GetM()) && (IndRow(k) <= i))
          k++;

        int size_row = k - k0;
        if (size_row > 0)
          {
            mat_csr.ReallocateRow(i, size_row);
            for (int j = 0; j < size_row; j++)
              {
                mat_csr.Index(i, j) = IndCol(k0 + j);
                mat_csr.Value(i, j) = Value(k0 + j);
              }
          }
        else
          mat_csr.ClearRow(i);
      }
  }


  //! Conversion from RowComplexSparse to ArrayRowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {  
    typedef typename ClassComplexType<T0>::Treal Treal;  
    int m = A.GetM();
    int n = A.GetN();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    Treal* data_imag = A.GetImagData();
    B.Reallocate(m, n);
    Vector<T0> value(n);
    Vector<int> col_num(n);
    for (int i = 0; i < m; i++)
      {
        long jr = ptr_real[i], ji = ptr_imag[i];
        int size_row = 0;
        while (jr < ptr_real[i+1])
          {
            int col = ind_real[jr];
            while ((ji < ptr_imag[i+1]) && (ind_imag[ji] < col))
              {
                col_num(size_row) = ind_imag[ji];
                value(size_row) = T0(0, data_imag[ji]);
                ji++; size_row++;
              }
            
            if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == col))
              {
                col_num(size_row) = col;
                value(size_row) = T0(data_real[jr],
				     data_imag[ji]);
                ji++; size_row++;
              }
            else
              {
                col_num(size_row) = col;
                value(size_row) = T0(data_real[jr], 0);
                size_row++;
              }
            
            jr++;
          }
        
        while (ji < ptr_imag[i+1])
          {
            col_num(size_row) = ind_imag[ji];
            value(size_row) = T0(0, data_imag[ji]);
            ji++; size_row++;
          }
        
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = col_num(j);
            B.Value(i, j) = value(j);
          }
      }
  }


  //! Conversion from RowSymComplexSparse to ArrayRowSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
	     Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  { 
    typedef typename ClassComplexType<T0>::Treal Treal;   
    int m = A.GetM();
    int n = A.GetN();
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    Treal* data_imag = A.GetImagData();
    B.Reallocate(m, n);
    Vector<T0> value(n);
    Vector<int> col_num(n);
    for (int i = 0; i < m; i++)
      {
        long jr = ptr_real[i], ji = ptr_imag[i];
        int size_row = 0;
        while (jr < ptr_real[i+1])
          {
            int col = ind_real[jr];
            while ((ji < ptr_imag[i+1]) && (ind_imag[ji] < col))
              {
                col_num(size_row) = ind_imag[ji];
                value(size_row) = T0(0, data_imag[ji]);
                ji++; size_row++;
              }
            
            if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == col))
              {
                col_num(size_row) = col;
                value(size_row) = T0(data_real[jr],
				     data_imag[ji]);
                ji++; size_row++;
              }
            else
              {
                col_num(size_row) = col;
                value(size_row) = T0(data_real[jr], 0);
                size_row++;
              }
            
            jr++;
          }
        
        while (ji < ptr_imag[i+1])
          {
            col_num(size_row) = ind_imag[ji];
            value(size_row) = T0(0, data_imag[ji]);
            ji++; size_row++;
          }
        
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = col_num(j);
            B.Value(i, j) = value(j);
          }
      }
  }

}


#define SELDON_FILE_MATRIX_COMPLEX_CONVERSIONS_CXX
#endif
