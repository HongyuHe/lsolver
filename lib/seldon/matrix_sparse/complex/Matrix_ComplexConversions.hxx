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

#ifndef SELDON_FILE_MATRIX_COMPLEX_CONVERSIONS_HXX

namespace Seldon
{
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);

  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowComplexSparse,
				 Allocator4>& A, int index = 0);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColComplexSparse,
				 Allocator4>& A, int index = 0);


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymComplexSparse,
				 Allocator4>& A, int index = 0);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymComplexSparse,
				 Allocator4>& A, int index = 0);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowComplexSparse,
				 Allocator4>& A, int index = 0);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColComplexSparse,
				 Allocator4>& A, int index = 0);

  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymComplexSparse,
				 Allocator4>& A, int index = 0);

  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(const Vector<int, VectFull, Allocator1>& IndRow,
				 const Vector<int, VectFull, Allocator2>& IndCol,
				 const Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymComplexSparse,
				 Allocator4>& A, int index = 0);  
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val,
                    bool sym_pat = false);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat = false);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val,
                    bool sym_pat = false);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val,
                    bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat = false);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat = false);
  
  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val,
                    bool sym_pat = false);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat = false);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val,
                    bool sym_pat = false);

  
  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val,
                    bool sym_pat = false);

  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat = false);


  template<class T, class Prop, class Alloc1, class Tint0,
           class Tint1, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint0, VectFull, Alloc2>& Ptr,
                    Vector<Tint1, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat = false);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSymSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSymSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSymSparse, Allocator1>& B);;
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSymSparse, Allocator1>& B);;
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
             Matrix<T1, Prop1, RowComplexSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColComplexSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
             Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
             Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
             Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
             Matrix<T1, Prop1, ArrayColComplexSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
             Matrix<T1, Prop1, ArrayColSymComplexSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
             Matrix<T1, Prop1, ArrayRowSymComplexSparse, Allocator1>& B);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0,
             ArrayRowSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
             Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& B);
  
  template<class T0, class Prop0, class Allocator0>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
             Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& B);
  
  template<class T0, class Prop0, class Allocator0>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
             Matrix<T0, Prop0, RowComplexSparse, Allocator0>& B);
  
  template<class T0, class Prop0, class Allocator0>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& A,
             Matrix<T0, Prop0, ColComplexSparse, Allocator0>& B);
  
  template<class T0, class Prop0, class Allocator0>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
             Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& B);
  
  template<class T0, class Prop0, class Allocator0>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
             Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& B);
  
  template<class T0, class Prop0, class Allocator0>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
             Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& B);
  
  template<class T0, class Prop0, class Allocator0>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& A,
             Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
             Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& mat_csr);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
             Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);  
  
} 

#define SELDON_FILE_MATRIX_COMPLEX_CONVERSIONS_HXX
#endif 
