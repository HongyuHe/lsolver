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


#ifndef SELDON_FILE_RELAXATION_MATVECT_HXX

namespace Seldon
{

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega,
		 int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ColSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayColSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ColSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayColSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);

} // end namespace

#define SELDON_FILE_RELAXATION_MATVECT_HXX
#endif
