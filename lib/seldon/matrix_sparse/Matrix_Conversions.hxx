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


#ifndef SELDON_FILE_MATRIX_CONVERSIONS_HXX


/*
  Functions defined in this file:
  
  conversion of a sparse matrix into coordinate format
  (i, j, val) = A  
  ConvertMatrix_to_Coordinates(A, index_row, index_col, val, index, sym)
  
  conversion from coordinate format to other sparse format
  A = (i, j, val)
  ConvertMatrix_from_Coordinates(index_row, index_col, val, A, index)
  
  conversion from sparse matrices to Compressed Sparse Column
  (Ptr, Ind, Val) = A
  ConvertToCSC(A, Ptr, Ind, Val, sym_pat)
  
  conversion from sparse matrices to Compressed Sparse Row
  (Ptr, Ind, Val) = A
  ConvertToCSR(A, Ptr, Ind, Val, sym_pat)
  
  conversion between sparse matrices
  B = A
  CopyMatrix(A, B)
  
  pattern of A + A' in CSR format
  (Ptr, Ind) = A
  GetSymmetricPattern
  
  conversion from a sparse matrix to a dense matrix
  B = A
  CopyMatrix(A, B)
  
  conversion from a dense matrix to a sparse matrix
  B = sparse(A.*(abs(A) > threshold))
  ConvertToSparse(A, B, threshold)
  
*/

namespace Seldon
{


  /*
    From CSR formats to "Matlab" coordinate format.
  */


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


    template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  /*
    From "Matlab" coordinate format to CSR formats.
  */


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSparse, Allocator3>& A,
				 int index = 0);


#ifndef SWIG

  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>&,
				 const Vector<int, VectFull, Allocator2>&,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSparse, Allocator3>& A,
				 int index = 0);


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>&,
				 const Vector<int, VectFull, Allocator2>&,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymSparse, Allocator3>& A,
				 int index = 0);


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>&,
				 const Vector<int, VectFull, Allocator2>&,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymSparse, Allocator3>& A,
				 int index = 0);
  
  
  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>&,
				 const Vector<int, VectFull, Allocator2>&,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSparse,
				 Allocator3>& A, int index = 0);


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>&,
				 const Vector<int, VectFull, Allocator2>&,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSparse,
				 Allocator3>& A, int index = 0);


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>&,
				 const Vector<int, VectFull, Allocator2>&,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymSparse,
				 Allocator3>& A, int index = 0);


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>&,
				 const Vector<int, VectFull, Allocator2>&,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymSparse,
				 Allocator3>& A, int index = 0);
  
#endif // SWIG.


  /*
    From Sparse formats to CSC format
  */


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat = false);
  
    
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat = false);
  
      
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat = false);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& AllVal, bool sym_pat = false);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowMajor, Allocator>& A,
		  Matrix<T, Prop, RowMajor, Allocator>& B);
  

  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSymPacked, Allocator>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColMajor, Allocator>& A,
		  Matrix<T, Prop, ColMajor, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSymPacked, Allocator>& A,
		  Matrix<T, Prop, ColSymPacked, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Allocator>& A,
		  Matrix<T, Prop, RowSparse, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Allocator>& A,
		  Matrix<T, Prop, RowSymSparse, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSparse, Allocator>& A,
		  Matrix<T, Prop, ColSparse, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSymSparse, Allocator>& A,
		  Matrix<T, Prop, ColSymSparse, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayRowSymSparse, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayRowSparse, Allocator>& B);
  
  
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayColSparse, Allocator>& B);
  

  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayColSymSparse, Allocator>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csc);
  
  
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSparse, Alloc2>& B);

  
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Alloc1>& A,
		  Matrix<T, Prop, ArrayColSparse, Alloc2>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B);
  
  
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B);
  

  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B);


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, RowSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B);
  
  
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B);
  
  
  /*
    From Sparse formats to CSR format
  */
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
              class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& IndRow,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& IndRow,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val);
    

  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& IndRow,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& IndRow,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val);

  
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T1, Prop1, ColSparse, Alloc1>& A,
		  Matrix<T2, Prop2, RowSparse, Alloc2>& B);
  
  
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T1, Prop1, ArrayColSparse, Alloc1>& A,
		  Matrix<T2, Prop2, RowSparse, Alloc2>& B);
  
    
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr);


  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, RowSymSparse, Alloc2>& B);


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, RowSparse, Alloc2>& B);
  
  
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, RowSymSparse, Alloc2>& B);


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, RowSparse, Alloc2>& B);

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, RowSymSparse, Allocator1>& B);


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& B);
  
  
  /***********************************
   * From ArraySparse to ArraySparse *
   ***********************************/
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSymSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& Acsc,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& Acsc,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B);
    
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B);
  
  
  /***********************
   * GetSymmetricPattern *
   ***********************/


  template<class T, class Prop, class Storage, class Allocator, class Tint0,
           class Tint1, class Allocator2, class Allocator3>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Vector<Tint0, VectFull, Allocator2>& Ptr,
                           Vector<Tint1, VectFull, Allocator3>& Ind);
  
  
  template<class T, class Prop, class Storage, class Allocator, class AllocI>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Matrix<int, Symmetric, RowSymSparse, AllocI>& B);
  
  
  /*****************************************************
   * Conversion from sparse matrices to dense matrices *
   *****************************************************/
  
  
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Allocator1>& A,
		  Matrix<T, Prop, RowMajor, Allocator2>& B);
  
  
  template<class T1, class T2, class Prop1, class Prop2,
	   class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		  Matrix<T2, Prop2, RowMajor, Allocator2>& B);
  

  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator2>& B);

  
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator2>& B);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, ColSymPacked, Allocator2>& B);


  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSym, Allocator2>& B);
  

  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, ColSym, Allocator2>& B);
  

  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, RowMajor, Alloc2>& B);


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, RowMajor, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIDense, Alloc2>& B);

  
  /*****************************************************
   * Conversion from dense matrices to sparse matrices *
   *****************************************************/
  
  
  template<class T>
  void ConvertToSparse(const Matrix<T, Symmetric, RowSymPacked>& A,
                       Matrix<T, Symmetric, RowSymSparse>& B,
		       const T& threshold);


  template<class T>
  void ConvertToSparse(const Matrix<T, General, RowMajor>& A,
                       Matrix<T, General, ArrayRowSparse>& B,
		       const T& threshold);


  template<class T>
  void ConvertToSparse(const Matrix<T, General, RowMajor>& A,
                       Matrix<T, General, RowSparse>& B,
		       const T& threshold);

} // namespace Seldon.


#define SELDON_FILE_MATRIX_CONVERSIONS_HXX
#endif
