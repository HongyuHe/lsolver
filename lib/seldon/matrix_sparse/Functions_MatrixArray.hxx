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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_ARRAY_HXX

namespace Seldon
{

  template<class T0, class Allocator0,
	   class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, ArrayRowSparse, Allocator0>& A,
	      int i, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, ArrayColSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X);
  
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X);
  
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X);
  
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, ArrayRowSparse, Allocator0>& M);
  
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, ArrayColSparse, Allocator0>& M);
  
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M);
  
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M);
  
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, ArrayRowSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X);
  
  template<class T0, class Allocator0,
	   class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, ArrayColSparse, Allocator0>& A,
	      int j, Vector<T1, VectSparse, Allocator1>& X);
  
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X);

  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, ArrayRowSparse, Allocator0>& M);
  
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, ArrayColSparse, Allocator0>& M);

  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M);

  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, Symmetric,
		 ArrayRowSymSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		 ArrayRowSymSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, Symmetric,
		 ArrayColSymSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		 ArrayColSymSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric,
		    ArrayRowSymSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);
  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		    ArrayRowSymSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);
  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric,
		    ArrayColSymSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		    ArrayColSymSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);
  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		 Matrix<T2, General, ArrayRowSparse, Allocator2>& B);
  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
		 Matrix<T2, General, ArrayColSparse, Allocator2>& B);
  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
		 Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B);
  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, ArrayColSymSparse, Allocator1>& A,
		 Matrix<T2, Symmetric, ArrayColSymSparse, Allocator2>& B);
  
  template<class T0, class T1, class T2, class T3, class Allocator1,
	   class Allocator2, class Allocator3>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		 const Matrix<T2, General, ArrayRowSparse, Allocator2>& B,
		 Matrix<complex<T3>, General, ArrayRowSparse, Allocator3>& C);
  
  template<class T0, class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
		 const Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B,
		 Matrix<complex<T3>, Symmetric, ArrayRowSymSparse, Allocator3>& C);

  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
		 Matrix<T2, General, ArrayRowSparse, Allocator2>& B);
  
  template<class T0, class T, class Allocator>
  void MltScalar(const T0& alpha, Matrix<T, General, ArrayRowSparse, Allocator>& A);

  template<class T0, class T, class Allocator>
  void MltScalar(const T0& alpha, Matrix<T, General, ArrayColSparse, Allocator>& A);

  template<class T0, class T, class Allocator>
  void MltScalar(const T0& alpha,
		 Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A);

  template<class T0, class T, class Allocator>
  void MltScalar(const T0& alpha,
		 Matrix<T, Symmetric, ArrayColSymSparse, Allocator>& A);
  
  template<class T0, class T1, class Prop1, class Allocator1, class T4,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T3, class Prop3, class Storage3, class Allocator3>
  void MltAddMatrix(const T0& alpha,
		    const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		    const Matrix<T2, Prop2, Storage2, Allocator2>& B,
		    const T4& beta,
		    Matrix<T3, Prop3, Storage3, Allocator3>& C);

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
	   class T4, class T2, class Prop2, class Allocator2,
           class T3, class Prop3, class Storage3, class Allocator3>
  void MltAddMatrix(const T0& alpha,
		    const Matrix<T1, Prop1, Storage1, Allocator1>& A,
		    const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		    const T4& beta,
		    Matrix<T3, Prop3, Storage3, Allocator3>& C);
  
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A);
  
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayColSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayColSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayColSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A);
  
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A);

  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                 Matrix<T, General, ArrayRowSparse, Allocator>& B);

  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayRowSparse, Allocator>& A);

  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayColSparse, Allocator>& A,
                 Matrix<T, General, ArrayColSparse, Allocator>& B);

  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayColSparse, Allocator>& A);

  template<class T, class Allocator>
  void Conjugate(Matrix<T, General, ArrayRowSparse, Allocator>& A);

  template<class T, class Allocator>
  void Conjugate(Matrix<T, General, ArrayColSparse, Allocator>& A);

  template<class T, class Allocator>
  void Conjugate(Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A);

  template<class T, class Allocator>
  void Conjugate(Matrix<T, Symmetric, ArrayColSymSparse, Allocator>& A);

  template<class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2,
           class T4, class Prop4, class Allocator4>
  void MltMatrix(const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		 const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		 Matrix<T4, Prop4, ArrayRowSparse, Allocator4>& C);
  
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2, class T3,
           class T4, class Prop4, class Allocator4>
  void MltAddMatrix(const T0& alpha,
		    const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		    const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		    const T3& beta,
		    Matrix<T4, Prop4, ArrayRowSparse, Allocator4>& C);

  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2, class T3,
           class T4, class Prop4, class Allocator4>
  void MltAddMatrix(const T0& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
		    const T3& beta,
		    Matrix<T4, Prop4, ArrayRowSparse, Allocator4>& C);
  
  template<class T, class Prop, class  Storage, class Allocator, class T0>
  void RemoveSmallEntry(Matrix<T, Prop, Storage, Allocator>& A, const T0&);

  template<class T, class Prop, class Allocator, class T0>
  void RemoveSmallEntry(Matrix<T, Prop, RowSparse, Allocator>& A, const T0&);

  template<class T, class Prop, class Allocator, class T0>
  void RemoveSmallEntry(Matrix<T, Prop, RowSymSparse, Allocator>& A,
			const T0&);

  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSymSparse, Allocator>& A);
  
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSparse, Allocator>& A);
  
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, RowSparse, Allocator>& A);

  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, RowSymSparse, Allocator>& A);

  template<class T1, class Prop, class Storage, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, Storage, Allocator>& A);
  
  template<class T1, class Prop, class Storage, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, Storage, Allocator>& A);
  
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSymSparse, Allocator>& A);

  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, ArrayRowSparse, Allocator>& A);

  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, RowSparse, Allocator>& A);

  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, RowSymSparse, Allocator>& A);

  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General,
                 ArrayRowSparse, Allocator> & mat_direct);
  
  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left, 
		 const Matrix<Complexe, Symmetric,
                 ArrayRowSymSparse, Allocator> & mat_direct);
  
  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General,
                 RowSparse, Allocator> & mat_direct);
  
  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left, 
		 const Matrix<Complexe, Symmetric,
                 RowSymSparse, Allocator> & mat_direct);
  
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General,
                 ArrayRowSparse, Allocator> & mat_direct);
  
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale_left, 
		 const Matrix<Complexe, Symmetric,
                 ArrayRowSymSparse, Allocator> & mat_direct);
  
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General,
                 RowSparse, Allocator> & mat_direct);
  
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale_left, 
		 const Matrix<Complexe, Symmetric,
                 RowSymSparse, Allocator> & mat_direct);

  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, General,
				 ArrayRowSparse, Allocator> & A);
  
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, Symmetric,
				 ArrayRowSymSparse, Allocator> & A);

  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, General,
				 RowSparse, Allocator> & A);
  
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, Symmetric,
				 RowSymSparse, Allocator> & A);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<T1, VectFull, Allocator1>& Value);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1> void
  CopySubMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		const IVect& row, const IVect& col,
		Vector<int>& RowNum,
		Vector<int>& ColNum,
		Vector<T1, VectFull, Allocator1>& Value);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<T1, VectFull, Allocator1>& Value);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Vector<int>& RowNum,
                     Vector<int>& ColNum,
                     Vector<T1, VectFull, Allocator1>& Value);
  
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Matrix<T1, Prop1, Storage1, Allocator1>& B);
  
} // namespace Seldon

#define SELDON_FILE_FUNCTIONS_MATRIX_ARRAY_HXX
#endif
