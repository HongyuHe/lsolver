// Copyright (C) 2001-2011 Marc Duruflé
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


#ifndef SELDON_FILE_FUNCTIONS_MATVECT_COMPLEX_HXX

namespace Seldon
{

  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, ColComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y);

  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y);

  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, ColComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y);

  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y);

  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, Symmetric,
		 ArrayRowSymComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		 ArrayRowSymComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, Symmetric,
		 ArrayColSymComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		 ArrayColSymComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, General,
		 ArrayRowComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, General,
		 ArrayRowComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, General,
		 ArrayColComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);

  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, General,
		 ArrayColComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C);
  
  
  /****************
   * MltAddVector *
   ****************/

  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, ColComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, ColComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric,
		    ArrayRowSymComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		    ArrayRowSymComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);
  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric,
		    ArrayColSymComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);
  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		    ArrayColSymComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);
  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, General,
		    ArrayRowComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);
  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, General,
		    ArrayRowComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, General,
		    ArrayColComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);

  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, General,
		    ArrayColComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C);
    
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);
}

#define SELDON_FILE_FUNCTIONS_MATVECT_COMPLEX_HXX
#endif

