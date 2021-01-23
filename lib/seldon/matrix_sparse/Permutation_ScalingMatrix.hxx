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

#ifndef SELDON_FILE_PERMUTATION_SCALING_MATRIX_HXX

namespace Seldon
{
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowSparse, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColSparse, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowSymSparse, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColSymSparse, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ArrayRowSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ArrayColSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm);
  
  template<class T, class Prop, class Allocator>
  void
  ApplyInversePermutation(Matrix<T, Prop, ArrayColSymSparse, Allocator>& A,
                          const IVect& row_perm, const IVect& col_perm);
  
  template<class T, class Prop, class Allocator>
  void
  ApplyInversePermutation(Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
                          const IVect& row_perm, const IVect& col_perm);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowSymSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColSymSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayColSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayColSymSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);

  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, RowSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);
  
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ColSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);

  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, RowSymSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);
  
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ColSymSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);

  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);

  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayColSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);

  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowSymSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);

 template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayColSymSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);

  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, RowSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale);

  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ColSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale);

  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ArrayRowSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale);
  
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ArrayColSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale);

  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, RowSparse, Allocator1>& A,
                        const Vector<T2, VectFull, Allocator2>& scale);
  
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, ColSparse, Allocator1>& A,
                        const Vector<T2, VectFull, Allocator2>& scale);

  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, ArrayRowSparse, Allocator1>& A,
                        const Vector<T2, VectFull, Allocator2>& scale);

 template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, ArrayColSparse, Allocator1>& A,
                        const Vector<T2, VectFull, Allocator2>& scale);

}

#define SELDON_FILE_PERMUTATION_SCALING_MATRIX_HXX
#endif
