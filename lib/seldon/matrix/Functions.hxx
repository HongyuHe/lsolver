// Copyright (C) 2001-2011 Vivien Mallet, Marc Fragu, Marc Duruflé
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


#ifndef SELDON_FILE_FUNCTIONS_HXX

namespace Seldon
{

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, RowSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, ColSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, Symmetric, RowSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, Symmetric, ColSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetRow(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	      int i, Vector<T1, Storage1, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, RowSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, ColSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, Symmetric, ColSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, Symmetric, RowSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Prop0, class Allocator0,
  class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M,
              int j, Vector<T1, PETScSeq, Allocator1>& X);

  template <class T0, class Prop0, class Allocator0,
  class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
              int j, Vector<T1, PETScPar, Allocator1>& X);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	      int j, Vector<T1, Storage1, Allocator1>& X);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, Storage0, Allocator0>& M_in,
	      int begin, int end,
              Matrix<T1, Prop1, Storage1, Allocator1>& M_out);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, Storage0, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
  class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScSeq, Allocator1>& X,
              int i, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
  class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScPar, Allocator1>& X,
              int i, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, RowSparse, Allocator0>& M);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, ColSparse, Allocator0>& M);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, RowSymSparse, Allocator0>& M);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, ColSymSparse, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, RowLoTriang, Allocator0>& M);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, RowLoTriangPacked, Allocator0>& M);
    
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, ColLoTriang, Allocator0>& M);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, ColLoTriangPacked, Allocator0>& M);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, RowUpTriang, Allocator0>& M);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, RowUpTriangPacked, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, ColUpTriang, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, ColUpTriangPacked, Allocator0>& M);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, Storage0, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
  class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScSeq, Allocator1>& X,
              int j, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
  class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScPar, Allocator1>& X,
              int j, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M);

  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectFull, Allocator1>& X,
	      int j, Matrix<T0, General, RowSparse, Allocator0>& M);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, RowSparse, Allocator0>& M);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, ColSparse, Allocator0>& M);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, RowSymSparse, Allocator0>& M);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, ColSymSparse, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, RowLoTriang, Allocator0>& M);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, RowLoTriangPacked, Allocator0>& M);
  
  template <class T0, class Prop0, class Allocator0,
            class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, ColLoTriang, Allocator0>& M);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, ColLoTriangPacked, Allocator0>& M);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, RowUpTriang, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, RowUpTriangPacked, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, ColUpTriang, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, ColUpTriangPacked, Allocator0>& M);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowMajor, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColMajor, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowSymPacked, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColSymPacked, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowSym, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColSym, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowHermPacked, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColHermPacked, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowHerm, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColHerm, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowMajor, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColMajor, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);
  
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowSymPacked, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColSymPacked, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowSym, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColSym, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowHermPacked, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColHermPacked, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowHerm, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColHerm, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index = 0);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowMajor, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColMajor, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowSymPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColSymPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowSym, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColSym, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowHermPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColHermPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowHerm, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColHerm, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowLoTriangPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowLoTriang, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColLoTriangPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColLoTriang, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowUpTriangPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, RowUpTriang, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColUpTriangPacked, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, ColUpTriang, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowMajor, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColMajor, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowLoTriangPacked, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowLoTriang, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColLoTriangPacked, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColLoTriang, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowUpTriangPacked, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, RowUpTriang, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);

  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColUpTriangPacked, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);
  
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(Matrix<T, Prop, ColUpTriang, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);

  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowMajor, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColMajor, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowLoTriangPacked, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowLoTriang, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColLoTriangPacked, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColLoTriang, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowUpTriangPacked, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, RowUpTriang, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColUpTriangPacked, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T, Prop, ColUpTriang, Allocator>& A,
                        const Vector<T2, VectFull, Allocator2>& Dcol);
  
} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_HXX
#endif
