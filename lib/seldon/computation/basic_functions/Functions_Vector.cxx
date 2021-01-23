// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_VECTOR_CXX
#define SELDON_FILE_FUNCTIONS_VECTOR_CXX


#include "Functions_Vector.hxx"


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

  
  //! Multiplication of a vector by a scalar
  template <class T0,
	    class T1, class Storage1, class Allocator1>
  void MltScalar(const T0& alpha,
		 Vector<T1, Storage1, Allocator1>& X)  throw()
  {
    X *= alpha;
  }


  // MLT //
  /////////


  /////////
  // ADD //

  
  //! Adds two vectors Y = Y + alpha X
  template <class T0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, Storage1, Allocator1>& X,
		 Vector<T2, Storage2, Allocator2>& Y)
  {
    T0 zero; SetComplexZero(zero);
    if (alpha != zero)
      {
	long ma = X.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
	CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

	for (long i = 0; i < ma; i++)
	  Y(i) += alpha * X(i);
      }
  }


  //! Adds two vectors Y = beta Y + alpha X
  template <class T0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, Storage1, Allocator1>& X,
		 const T0& beta,
		 Vector<T2, Storage2, Allocator2>& Y)
  {
    T0 zero; SetComplexZero(zero);
    if (alpha != zero)
      {
	long ma = X.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
	CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

	for (long i = 0; i < ma; i++)
	  Y(i) = beta*Y(i) + alpha * X(i);
      }
    else
      Mlt(beta, Y);
  }

  
  //! Adds two vectors Y = Y + alpha X
  template <class T0,
            class T1, class Allocator1,
            class T2, class Allocator2>
  void AddVector(const T0 alpha,
		 const Vector<T1, PETScSeq, Allocator1>& X,
		 Vector<T2, PETScSeq, Allocator2>& Y)
  {
    if (alpha != T0(0))
      {
        T1 alpha_ = alpha;
	
#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

        VecAXPY(Y.GetPetscVector(), alpha_, X.GetPetscVector());
      }
  }


  template <class T0,
            class T1, class Allocator1,
            class T2, class Allocator2>
  void AddVector(const T0 alpha,
		 const Vector<T1, PETScPar, Allocator1>& X,
		 Vector<T2, PETScPar, Allocator2>& Y)
  {
    if (alpha != T0(0))
      {
        T1 alpha_ = alpha;
	
#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

        VecAXPY(Y.GetPetscVector(), alpha_, X.GetPetscVector());
      }
  }


  //! Adds two vectors Y = Y + alpha X
  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, VectSparse, Allocator1>& X,
		 Vector<T2, VectSparse, Allocator2>& Y)
  {
    T0 zero;
    SetComplexZero(zero);
    if (alpha != zero)
      {
	Vector<T1, VectSparse, Allocator1> Xalpha = X;
	Xalpha *= alpha;
	Y.AddInteractionRow(Xalpha.GetSize(),
			    Xalpha.GetIndex(), Xalpha.GetData(), true);
      }
  }


  //! Adds two vectors Y = Y + alpha X
  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, VectSparse, Allocator1>& X,
		 Vector<T2, VectFull, Allocator2>& Y)
  {
    T0 zero;
    SetComplexZero(zero);
    if (alpha != zero)
      {
	for (int i = 0; i < X.GetM(); i++)
	  Y(X.Index(i)) += alpha * X.Value(i);
      }
  }


  //! Adds two vectors Y = Y + alpha X
  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void AddVector(const T0& alpha,
		 const Vector<T1, Collection, Allocator1>& X,
		 Vector<T2, Collection, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "Add(X, Y)", "X + Y");
#endif

    for (int i = 0; i < X.GetNvector(); i++)
      Add(alpha, X.GetVector(i), Y.GetVector(i));
  }


  template <class T0,
	    class T1, template <class U1> class Allocator1,
	    class T2, template <class U2> class Allocator2>
  void AddVector(const T0& alpha,
		 const
		 Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
		 Vector<FloatDouble, DenseSparseCollection, Allocator2<T2> >& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "Add(X, Y)");
#endif

    AddVector(alpha, X.GetFloatDense(), Y.GetFloatDense());
    AdVectord(alpha, X.GetFloatSparse(), Y.GetFloatSparse());
    AddVector(alpha, X.GetDoubleDense(), Y.GetDoubleDense());
    AddVector(alpha, X.GetDoubleSparse(), Y.GetDoubleSparse());
  }


  //! Adds two vectors Y = Y + alpha X
  template <class T0,
	    class T1, template <class U1> class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void AddVector(const T0& alpha,
		 const
		 Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
		 Vector<T2, Storage2, Allocator2>& Y)
  {
    if (alpha != T0(0))
      {
	double alpha_ = alpha;

	int ma = X.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
	CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

	for (int i = 0; i < ma; i++)
	  Y(i) += alpha_ * X(i);

      }
  }


  // ADD //
  /////////


  //////////
  // COPY //

  
  //! Y = X
  template <class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CopyVector(const Vector<T1, Storage1, Allocator1>& X,
		  Vector<T2, Storage2, Allocator2>& Y)
  {
    Y.Copy(X);
  }


  //! Y = X
  template <class T1, class Allocator1,
	    class T2, class Allocator2>
  void CopyVector(const Vector<T1, Collection, Allocator1>& X,
		  Vector<T2, VectFull, Allocator2>& Y)
  {
    Y.Clear();
    for (int i = 0; i < X.GetNvector(); i++)
      Y.PushBack(X.GetVector(i));
  }


  template<class T, class Alloc1, class Alloc2>
  void CopyVector(const Vector<T, PETScPar, Alloc1>& A,
		  Vector<T, VectFull, Alloc2>& B)
  {
    B.Reallocate(A.GetLocalM());
    T *local_data;
    VecGetArray(A.GetPetscVector(), &local_data);
    for (int i = 0; i < A.GetLocalM(); i++)
      B(i) = local_data[i];
    VecRestoreArray(A.GetPetscVector(), &local_data);
  }


  template<class T, class Alloc1, class Alloc2>
  void CopyVector(const Vector<T, VectFull, Alloc1>& A,
		  Vector<T, PETScPar, Alloc2>& B)
  {
    T *local_data;
    VecGetArray(B.GetPetscVector(), &local_data);
    for (int i = 0; i < A.GetM(); i++)
      local_data[i] = A(i);
    VecRestoreArray(B.GetPetscVector(), &local_data);
  }



  // COPY //
  //////////


  //////////
  // SWAP //


  //! Swaps X and Y without copying all elements
  template <class T, class Storage, class Allocator>
  void Swap(Vector<T, Storage, Allocator>& X,
	    Vector<T, Storage, Allocator>& Y)
  {
    long nx = X.GetM();
    T* data = X.GetData();
    X.Nullify();
    X.SetData(Y.GetM(), Y.GetData());
    Y.Nullify();
    Y.SetData(nx, data);
  }


  //! Swaps X and Y without copying all elements
  template <class T, class Allocator>
  void Swap(Vector<T, VectSparse, Allocator>& X,
	    Vector<T, VectSparse, Allocator>& Y)
  {
    long nx = X.GetM();
    T* data = X.GetData();
    int* index = X.GetIndex();
    X.Nullify();
    X.SetData(Y.GetM(), Y.GetData(), Y.GetIndex());
    Y.Nullify();
    Y.SetData(nx, data, index);
  }


  //! Swaps X and Y with pointers
  template <class T>
  void SwapPointer(Vector<T>& X, Vector<T>& Y)
  {
    long nx = X.GetM();
    T* data = X.GetData();
    X.Nullify();
    X.SetData(Y.GetM(), Y.GetData());
    Y.Nullify();
    Y.SetData(nx, data);
  }
  

  // SWAP //
  //////////


  /////////////
  // DOTPROD //


  //! Scalar product between two vectors.
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  T1 DotProdVector(const Vector<T1, Storage1, Allocator1>& X,
		   const Vector<T2, Storage2, Allocator2>& Y)
  {
    T1 value;
    SetComplexZero(value);

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "DotProd(X, Y)");
#endif

    for (long i = 0; i < X.GetM(); i++)
      value += X(i) * Y(i);

    return value;
  }


  //! Scalar product between two vector collections.
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  typename T1::value_type
  DotProdVector(const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Collection, Allocator2>& Y)
  {
    typename T1::value_type value(0);

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "DotProd(X, Y)", "<X, Y>");
#endif

    for (int i = 0; i < X.GetNvector(); i++)
      value += DotProdVector(X.GetVector(i), Y.GetVector(i));
    return value;
  }


  //! Scalar product between two vector collections.
  template<class T1, template <class U1> class Allocator1,
	   class T2, template <class U2> class Allocator2>
  double
  DotProdVector(const
		Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& X,
		const
		Vector<FloatDouble, DenseSparseCollection, Allocator2<T2> >& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "DotProd(X, Y)");
#endif

    double value(0.);
    value += DotProdVector(X.GetFloatDense(), Y.GetFloatDense());
    value += DotProdVector(X.GetFloatSparse(), Y.GetFloatSparse());
    value += DotProdVector(X.GetDoubleDense(), Y.GetDoubleDense());
    value += DotProdVector(X.GetDoubleSparse(), Y.GetDoubleSparse());
    return value;
  }


  //! Scalar product between two vectors conj(X).Y .
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  T1 DotProdConjVector(const Vector<T1, Storage1, Allocator1>& X,
		       const Vector<T2, Storage2, Allocator2>& Y)
  {
    T1 value;
    SetComplexZero(value);

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(X, Y, "DotProdConj(X, Y)");
#endif

    for (long i = 0; i < X.GetM(); i++)
      value += conjugate(X(i)) * Y(i);

    return value;
  }


  //! Scalar product between two sparse vectors.
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1 DotProdVector(const Vector<T1, VectSparse, Allocator1>& X,
		   const Vector<T2, VectSparse, Allocator2>& Y)
  {
    T1 value;
    SetComplexZero(value);

    int size_x = X.GetSize();
    int size_y = Y.GetSize();
    int kx = 0, ky = 0, pos_x;
    while (kx < size_x)
      {
	pos_x = X.Index(kx);
	while (ky < size_y && Y.Index(ky) < pos_x)
	  ky++;

	if (ky < size_y && Y.Index(ky) == pos_x)
	  value += X.Value(kx) * Y.Value(ky);

	kx++;
      }

    return value;
  }


  //! Scalar product between a sparse vector and a dense vector
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1 DotProdVector(const Vector<T1, VectSparse, Allocator1>& X,
		   const Vector<T2, VectFull, Allocator2>& Y)
  {
    T1 value;
    SetComplexZero(value);
    for (int i = 0; i < X.GetM(); i++)
      value += X.Value(i)*Y(X.Index(i));
    
    return value;
  }


  //! Scalar product between a sparse vector and a dense vector
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1 DotProdVector(const Vector<T1, VectFull, Allocator1>& X,
		   const Vector<T2, VectSparse, Allocator2>& Y)
  {
    return DotProdVector(Y, X);
  }


  //! Scalar product between two sparse vectors conj(X).Y.
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1
  DotProdConjVector(const Vector<T1, VectSparse, Allocator1>& X,
		    const Vector<T2, VectSparse, Allocator2>& Y)
  {
    T1 value;
    SetComplexZero(value);

    int size_x = X.GetSize();
    int size_y = Y.GetSize();
    int kx = 0, ky = 0, pos_x;
    while (kx < size_x)
      {
	pos_x = X.Index(kx);
	while (ky < size_y && Y.Index(ky) < pos_x)
	  ky++;

	if (ky < size_y && Y.Index(ky) == pos_x)
	  value += conjugate(X.Value(kx)) * Y.Value(ky);

	kx++;
      }

    return value;
  }


  //! Scalar product between a sparse vector and dense one conj(X).Y.
  template<class T1, class Allocator1,
	   class T2, class Allocator2>
  T1
  DotProdConjVector(const Vector<T1, VectSparse, Allocator1>& X,
		    const Vector<T2, VectFull, Allocator2>& Y)
  {
    T1 value;
    SetComplexZero(value);

    for (int i = 0; i < X.GetM(); i++)
      value += conjugate(X.Value(i))*Y(X.Index(i));
    
    return value;
  }


  // DOTPROD //
  /////////////


  ///////////
  // NORM1 //

  
  //! returns 1-norm of X
  /*!
    For complex numbers, we use |z| = |Re(z)| + |Im(z)|
    so that the function is the same as Blas equivalent dzasum
  */
  template<class T1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  Norm1(const Vector<T1, Storage1, Allocator1>& X)
  {
    typename ClassComplexType<T1>::Treal value(0);
    
    for (long i = 0; i < X.GetM(); i++)
      value += ComplexAbs(X(i));

    return value;
  }


  //! returns 1-norm of X
  /*!
    For complex numbers, we use |z| = |Re(z)| + |Im(z)|
    so that the function is the same as Blas equivalent dzasum
  */
  template<class T1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  Norm1(const Vector<T1, VectSparse, Allocator1>& X)
  {
    typename ClassComplexType<T1>::Treal value(0);
    
    for (int i = 0; i < X.GetM(); i++)
      value += ComplexAbs(X.Value(i));

    return value;
  }

  
  // NORM1 //
  ///////////


  ///////////
  // NORM2 //


  //! returns 2-norm of X
  template<class T1, class Storage1, class Allocator1>
  T1 Norm2(const Vector<T1, Storage1, Allocator1>& X)
  {
    T1 value(0);

    for (long i = 0; i < X.GetM(); i++)
      value += X(i) * X(i);

    return sqrt(value);
  }


  //! returns 2-norm of X
  template<class T1, class Storage1, class Allocator1>
  T1 Norm2(const Vector<complex<T1>, Storage1, Allocator1>& X)
  {
    T1 value(0);

    for (size_t i = 0; i < X.GetSize(); i++)
      value += absSquare(X(i));

    return sqrt(value);
  }


  //! returns 2-norm of X
  template<class T1, class Allocator1>
  T1 Norm2(const Vector<T1, VectSparse, Allocator1>& X)
  {
    T1 value(0);

    for (size_t i = 0; i < X.GetSize(); i++)
      value += X.Value(i) * X.Value(i);

    return sqrt(value);
  }


  //! returns 2-norm of X
  template<class T1, class Allocator1>
  T1 Norm2(const Vector<complex<T1>, VectSparse, Allocator1>& X)
  {
    T1 value(0);

    for (size_t i = 0; i < X.GetSize(); i++)
      value += absSquare(X.Value(i));
    
    return sqrt(value);
  }


  // NORM2 //
  ///////////


  ////////////////////
  // GETMAXABSINDEX //


  //! returns index for which X(i) is maximal
  template<class T, class Storage, class Allocator>
  long GetMaxAbsIndex(const Vector<T, Storage, Allocator>& X)
  {
    return X.GetNormInfIndex();
  }


  // GETMAXABSINDEX //
  ////////////////////


  //////////////
  // APPLYROT //


  //! Computation of rotation between two points.
  template<class T>
  void GenRot(T& a_in, T& b_in, T& c_, T& s_)
  {
    // Old BLAS version.
    T roe;
    if (abs(a_in) > abs(b_in))
      roe = a_in;
    else
      roe = b_in;

    T scal = abs(a_in) + abs(b_in);
    T r, z;
    if (scal != T(0))
      {
	T a_scl = a_in / scal;
	T b_scl = b_in / scal;
	r = scal * sqrt(a_scl * a_scl + b_scl * b_scl);
	if (roe < T(0))
	  r *= T(-1);

	c_ = a_in / r;
	s_ = b_in / r;
	z = T(1);
	if (abs(a_in) > abs(b_in))
	  z = s_;
	else if (abs(b_in) >= abs(a_in) && c_ != T(0))
	  z = T(1) / c_;
      }
    else
      {
	c_ = T(1);
	s_ = T(0);
	r = T(0);
	z = T(0);
      }
    a_in = r;
    b_in = z;
  }


  //! Computation of rotation between two points.
  template<class T>
  void GenRot(complex<T>& a_in, complex<T>& b_in, T& c_, complex<T>& s_)
  {

    T a = abs(a_in), b = abs(b_in);
    if (a == T(0))
      {
	c_ = T(0);
	s_ = complex<T>(1, 0);
	a_in = b_in;
      }
    else
      {
	T scale = a + b;
	T a_scal = abs(a_in / scale);
	T b_scal = abs(b_in / scale);
	T norm = sqrt(a_scal * a_scal + b_scal * b_scal) * scale;

	c_ = a / norm;
	complex<T> alpha = a_in / a;
	s_ = alpha * conjugate(b_in) / norm;
	a_in = alpha * norm;
      }
    b_in = complex<T>(0, 0);
  }


  //! Rotation of a point in 2-D.
  template<class T>
  void ApplyRot(T& x, T& y, const T& c_, const T& s_)
  {
    T temp = c_ * x + s_ * y;
    y = c_ * y - s_ * x;
    x = temp;
  }


  //! Rotation of a complex point in 2-D.
  template<class T>
  void ApplyRot(complex<T>& x, complex<T>& y,
		const T& c_, const complex<T>& s_)
  {
    complex<T> temp = s_ * y + c_ * x;
    y = -conjugate(s_) * x + c_ * y;
    x = temp;
  }
  
  
  //! Rotation of a vector of points in 2-D
  template<class T, class Allocator1, class Allocator2>
  void ApplyRot(Vector<T, VectFull, Allocator1>& X,
		Vector<T, VectFull, Allocator2>& Y,
		const T& c, const T& s)
  {
    T tmp;
    for (long i = 0; i < X.GetM(); i++)
      {
	tmp = c*X(i) + s*Y(i);
	Y(i) = c*Y(i) - s*X(i);
	X(i) = tmp;
      }
  }
  
  
  //! Rotation of a vector of points in 2-D
  template<class T, class Allocator1, class Allocator2>
  void ApplyRot(Vector<T, VectSparse, Allocator1>& X,
		Vector<T, VectFull, Allocator2>& Y,
		const T& c, const T& s)
  {
    T tmp;
    for (int i = 0; i < X.GetM(); i++)
      {
	tmp = c*X.Value(i) + s*Y(X.Index(i));
	Y(X.Index(i)) = c*Y(X.Index(i)) - s*X.Value(i);
	X.Value(i) = tmp;
      }    
  }
  

  // APPLYROT //
  //////////////


  //////////////
  // CHECKDIM //


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void CheckDim(const Vector<T0, Storage0, Allocator0>& X,
		const Vector<T1, Storage1, Allocator1>& Y,
		string function, string op)
  {
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
    string Xchar = to_str(&X), Ychar = to_str(&Y);
#else
    string Xchar("X"), Ychar("Y");
#endif

    if (X.GetLength() != Y.GetLength())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + Xchar + string(") is a ")
		     + string("vector of length ") + to_str(X.GetLength())
		     + string(";\n     Y (") + Ychar + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetLength())
		     + string("."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, VectSparse, Allocator0>& X,
		const Vector<T1, VectSparse, Allocator1>& Y,
		string function, string op)
  {
    // The dimension of a Vector<VectSparse> is infinite,
    // so no vector dimension checking has to be done.
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, Collection, Allocator0>& X,
		const Vector<T1, Collection, Allocator1>& Y,
		string function, string op)
  {
    if (X.GetLength() != Y.GetLength())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetLength())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetLength())
		     + string("."));

    if (X.GetNvector() != Y.GetNvector())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetNvector())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetNvector())
		     + string("."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Allocator0, class Allocator00,
	    class T1, class Allocator1, class Allocator11>
  void CheckDim(const Vector<Vector<T0, VectSparse, Allocator0>,
                Collection, Allocator00>& X,
		const Vector<Vector<T1, VectSparse, Allocator1>,
                Collection, Allocator11>& Y,
		string function, string op)
  {
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
    string Xchar = to_str(&X), Ychar = to_str(&Y);
#else
    string Xchar("X"), Ychar("Y");
#endif

    if (X.GetNvector() != Y.GetNvector())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + Xchar + string(") is a ")
		     + string("vector of length ") + to_str(X.GetNvector())
		     + string(";\n     Y (") + Ychar + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetNvector())
		     + string("."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void CheckDim(const Vector<T0, VectFull, Allocator0>& X,
		Vector<T1, Collection, Allocator1>& Y,
		string function, string op)
  {
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
    string Xchar = to_str(&X), Ychar = to_str(&Y);
#else
    string Xchar("X"), Ychar("Y");
#endif

    if (X.GetLength() != Y.GetM())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + Xchar + string(") is a ")
		     + string("vector of length ") + to_str(X.GetLength())
		     + string(";\n     Y (") + Ychar + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetM())
		     + string("."));
  }


   //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, template <class U0> class Allocator0,
	    class T1, template <class U1> class Allocator1>
  void
  CheckDim(const
           Vector<FloatDouble, DenseSparseCollection, Allocator0<T0> >& X,
           const
           Vector<FloatDouble, DenseSparseCollection, Allocator1<T1> >& Y,
           string function, string op)
  {
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
    string Xchar = to_str(&X), Ychar = to_str(&Y);
#else
    string Xchar("X"), Ychar("Y");
#endif

    if (X.GetNvector() != Y.GetNvector())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + Xchar + string(") is a ")
		     + string("vector of length ") + to_str(X.GetNvector())
		     + string(";\n     Y (") + Ychar + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetNvector())
		     + string("."));
  }


  // CHECKDIM //
  //////////////


  ////////////////////
  // GATHER/SCATTER //


  template<class T, class Allocator1, class Allocator2>
  void GatherSparseEntry(const Vector<T, VectFull, Allocator1>& Y,
			 Vector<T, VectSparse, Allocator2>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X.Value(i) = Y(X.Index(i));
  }

  template<class T, class Allocator1, class Allocator2>
  void GatherSparseEntryZero(Vector<T, VectFull, Allocator1>& Y,
			     Vector<T, VectSparse, Allocator2>& X)
  {
    T zero; SetComplexZero(zero);
    for (int i = 0; i < X.GetM(); i++)
      {
	X.Value(i) = Y(X.Index(i));
	Y(X.Index(i)) = zero;
      }
  }
  

  template<class T, class Allocator1, class Allocator2>
  void ScatterSparseEntry(const Vector<T, VectSparse, Allocator1>& X,
			  Vector<T, VectFull, Allocator2>& Y)
  {
    for (int i = 0; i < X.GetM(); i++)
      Y(X.Index(i)) = X.Value(i);
  }
  

  // GATHER/SCATTER //
  ////////////////////  

  
  ///////////////
  // CONJUGATE //


  //! Sets a vector to its conjugate.
  template<class T, class Prop, class Allocator>
  void Conjugate(Vector<T, Prop, Allocator>& X)
  {
    for (long i = 0; i < X.GetM(); i++)
      X(i) = conjugate(X(i));
  }


  //! Sets a vector to its conjugate.
  template<class T, class Allocator>
  void Conjugate(Vector<T, VectSparse, Allocator>& X)
  {
    for (size_t i = 0; i < X.GetSize(); i++)
      X.Value(i) = conjugate(X.Value(i));
  }


  // CONJUGATE //
  ///////////////


} // namespace Seldon.


#endif
