// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu
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


#ifndef SELDON_FILE_FUNCTIONS_VECTOR_HXX

#define SELDON_FILE_FUNCTIONS_VECTOR_HXX


/*
  Functions defined in this file:

  alpha X -> X
  Mlt(alpha, X)

  alpha X + Y -> Y
  Add(alpha, X, Y)

  alpha X + beta Y -> Y
  Add(alpha, X, beta, Y)

  X -> Y
  Copy(X, Y)

  X <-> Y
  Swap(X, Y)

  X.Y
  DotProd(X, Y)
  DotProdConj(X, Y)

  ||X||
  Norm1(X)
  Norm2(X)
  GetMaxAbsIndex(X)

  Omega X
  GenRot(x, y, cos, sin)
  ApplyRot(x, y, cos, sin)

*/

namespace Seldon
{

  /////////
  // MLT //
  
  
  template <class T0,
	    class T1, class Storage1, class Allocator1>
  void MltScalar(const T0& alpha,
		 Vector<T1, Storage1, Allocator1>& X)  throw();
  
  
  // MLT //
  /////////
  
  
  /////////
  // ADD //
  
  
  template <class T0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, Storage1, Allocator1>& X,
		 Vector<T2, Storage2, Allocator2>& Y);
  
  template <class T0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, Storage1, Allocator1>& X,
		 const T0& beta,
		 Vector<T2, Storage2, Allocator2>& Y);  
  
  template <class T0,
            class T1, class Allocator1,
            class T2, class Allocator2>
  void AddVector(const T0 alpha,
		 const Vector<T1, PETScSeq, Allocator1>& X,
		 Vector<T2, PETScSeq, Allocator2>& Y);
  
  template <class T0,
            class T1, class Allocator1,
            class T2, class Allocator2>
  void AddVector(const T0 alpha,
		 const Vector<T1, PETScPar, Allocator1>& X,
		 Vector<T2, PETScPar, Allocator2>& Y);
  
  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, VectSparse, Allocator1>& X,
		 Vector<T2, VectSparse, Allocator2>& Y);
  
  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, VectSparse, Allocator1>& X,
		 Vector<T2, VectFull, Allocator2>& Y);
  
  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, Collection, Allocator1>& X,
		 Vector<T2, Collection, Allocator2>& Y);
  
  template <class T0,
	    class T1, template <class U1> class Allocator1,
	    class T2, template <class U2> class Allocator2>
  void AddVector(const T0& alpha,
		 const
		 Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
		 Vector<FloatDouble, DenseSparseCollection, Allocator2<T2> >& Y);
  
  template <class T0,
	    class T1, template <class U1> class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void AddVector(const T0& alpha,
		 const
		 Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
		 Vector<T2, Storage2, Allocator2>& Y);
  
  // ADD //
  /////////
  
  
  //////////
  // COPY //
  
  
  template <class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CopyVector(const Vector<T1, Storage1, Allocator1>& X,
		  Vector<T2, Storage2, Allocator2>& Y);
  
  template <class T1, class Allocator1,
	    class T2, class Allocator2>
  void CopyVector(const Vector<T1, Collection, Allocator1>& X,
		  Vector<T2, VectFull, Allocator2>& Y);
  
  template<class T, class Alloc1, class Alloc2>
  void CopyVector(const Vector<T, PETScPar, Alloc1>& A,
		  Vector<T, VectFull, Alloc2>& B);

  template<class T, class Alloc1, class Alloc2>
  void CopyVector(const Vector<T, VectFull, Alloc1>& A,
		  Vector<T, PETScPar, Alloc2>& B);
  
  
  // COPY //
  //////////
  
  
  //////////
  // SWAP //
  
  
  template <class T, class Storage, class Allocator>
  void Swap(Vector<T, Storage, Allocator>& X,
	    Vector<T, Storage, Allocator>& Y);
  
  template <class T, class Allocator>
  void Swap(Vector<T, VectSparse, Allocator>& X,
	    Vector<T, VectSparse, Allocator>& Y);
  
  template<class T>
  void SwapPointer(Vector<T>& X, Vector<T>& Y);


  // SWAP //
  //////////
  
  
  /////////////
  // DOTPROD //
  
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  T1 DotProdVector(const Vector<T1, Storage1, Allocator1>& X,
		   const Vector<T2, Storage2, Allocator2>& Y);
  
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  typename T1::value_type
  DotProdVector(const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Collection, Allocator2>& Y);
  
  template<class T1, template <class U1> class Allocator1,
	   class T2, template <class U2> class Allocator2>
  double
  DotProdVector(const
		Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
		const
		Vector<FloatDouble, DenseSparseCollection, Allocator2<T2> >& Y);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  T1 DotProdConjVector(const Vector<T1, Storage1, Allocator1>& X,
		       const Vector<T2, Storage2, Allocator2>& Y);
  
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1 DotProdVector(const Vector<T1, VectSparse, Allocator1>& X,
		   const Vector<T2, VectSparse, Allocator2>& Y);

  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1 DotProdVector(const Vector<T1, VectSparse, Allocator1>& X,
		   const Vector<T2, VectFull, Allocator2>& Y);

  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1 DotProdVector(const Vector<T1, VectFull, Allocator1>& X,
		   const Vector<T2, VectSparse, Allocator2>& Y);

  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1
  DotProdConjVector(const Vector<T1, VectSparse, Allocator1>& X,
		    const Vector<T2, VectSparse, Allocator2>& Y);
  
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1 
  DotProdConjVector(const Vector<T1, VectSparse, Allocator1>& X,
		    const Vector<T2, VectFull, Allocator2>& Y);
  
  
  // DOTPROD //
  /////////////


  ///////////
  // NORM1 //


  template<class T1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  Norm1(const Vector<T1, Storage1, Allocator1>& X);
  
  template<class T1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  Norm1(const Vector<T1, VectSparse, Allocator1>& X);
  

  // NORM1 //
  ///////////


  ///////////
  // NORM2 //

  
  template<class T1, class Storage1, class Allocator1>
  T1 Norm2(const Vector<T1, Storage1, Allocator1>& X);
  
  template<class T1, class Storage1, class Allocator1>
  T1 Norm2(const Vector<complex<T1>, Storage1, Allocator1>& X);
  
  template<class T1, class Allocator1>
  T1 Norm2(const Vector<T1, VectSparse, Allocator1>& X);
  
  template<class T1, class Allocator1>
  T1 Norm2(const Vector<complex<T1>, VectSparse, Allocator1>& X);
  

  // NORM2 //
  ///////////


  ////////////////////
  // GETMAXABSINDEX //

  
  template<class T, class Storage, class Allocator>
  long GetMaxAbsIndex(const Vector<T, Storage, Allocator>& X);


  // GETMAXABSINDEX //
  ////////////////////


  //////////////
  // APPLYROT //

  
  template<class T>
  void GenRot(T& a_in, T& b_in, T& c_, T& s_);
  
  template<class T>
  void GenRot(complex<T>& a_in, complex<T>& b_in, T& c_, complex<T>& s_);
  
  template<class T>
  void ApplyRot(T& x, T& y, const T& c_, const T& s_);
  
  template<class T>
  void ApplyRot(complex<T>& x, complex<T>& y,
		const T& c_, const complex<T>& s_);
  
  template<class T, class Allocator1, class Allocator2>
  void ApplyRot(Vector<T, VectFull, Allocator1>& X,
		Vector<T, VectFull, Allocator2>& Y,
		const T& c, const T& s);

  template<class T, class Allocator1, class Allocator2>
  void ApplyRot(Vector<T, VectSparse, Allocator1>& X,
		Vector<T, VectFull, Allocator2>& Y,
		const T& c, const T& s);

  
  // APPLYROT //
  //////////////


  //////////////
  // CHECKDIM //


  template <class T0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void CheckDim(const Vector<T0, Storage0, Allocator0>& X,
		const Vector<T1, Storage1, Allocator1>& Y,
		string function = "", string op = "X + Y");

  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, VectSparse, Allocator0>& X,
		const Vector<T1, VectSparse, Allocator1>& Y,
		string function = "", string op = "X + Y");

  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, Collection, Allocator0>& X,
		const Vector<T1, Collection, Allocator1>& Y,
		string function = "", string op = "X + Y");

  template <class T0, class Allocator0, class Allocator00,
	    class T1, class Allocator1, class Allocator11>
  void CheckDim(const Vector<Vector<T0, VectSparse, Allocator0>,
                Collection, Allocator00>& X,
		const Vector<Vector<T1, VectSparse, Allocator1>,
                Collection, Allocator11>& Y,
		string function = "", string op = "X + Y");

  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, VectFull, Allocator0>& X,
		Vector<T1, Collection, Allocator1>& Y,
		string function = "", string op = "X + Y");

  template <class T0, template <class U0> class Allocator0,
	    class T1, template <class U1> class Allocator1>
  void
  CheckDim(const
           Vector<FloatDouble, DenseSparseCollection, Allocator0<T0> >& X,
           const
           Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& Y,
           string function = "", string op = "X + Y");


  // CHECKDIM //
  //////////////

  
  ////////////////////
  // GATHER/SCATTER //


  template<class T, class Allocator1, class Allocator2>
  void GatherSparseEntry(const Vector<T, VectFull, Allocator1>& y,
			 Vector<T, VectSparse, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void GatherSparseEntryZero(Vector<T, VectFull, Allocator1>& y,
			     Vector<T, VectSparse, Allocator2>& x);
  
  template<class T, class Allocator1, class Allocator2>
  void ScatterSparseEntry(const Vector<T, VectSparse, Allocator1>& x,
			  Vector<T, VectFull, Allocator2>& y);
  

  // GATHER/SCATTER //
  ////////////////////  


  ///////////////
  // CONJUGATE //

  
  template<class T, class Prop, class Allocator>
  void Conjugate(Vector<T, Prop, Allocator>& X);
  
  template<class T, class Allocator>
  void Conjugate(Vector<T, VectSparse, Allocator>& X);
  
  
  // CONJUGATE //
  ///////////////


} // namespace Seldon.

#endif
