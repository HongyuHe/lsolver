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


#ifndef SELDON_FILE_FUNCTIONS_MATVECT_COMPLEX_CXX

/*
  Functions defined in this file:
  (storage RowComplexSparse, ArrayRowComplexSparse, etc)
  
  alpha.M*X + beta.Y -> Y
  MltAdd(alpha, M, X, beta, Y)
  
  SOR(A, X, B, omega, nb_iter, type_ssor)
  
*/

namespace Seldon
{

  /*************
   * MltVector *
   *************/


  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    int i; long j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    typename ClassComplexType<T1>::Treal rzero(0);
    T1 zero(0, 0);
    T1 temp;

    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * X(real_ind[j]);
	
        Y(i) = temp;
      }

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += T1(rzero, imag_data[j]) * X(imag_ind[j]);
	
        Y(i) += temp;
      }
  }
  
  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, ColComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    int i; long j;

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    typename ClassComplexType<T1>::Treal rzero(0);
    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, ColComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, ColComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    Y.Fill(0);
    for (i = 0; i < M.GetN(); i++)
      {
        for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
          Y(real_ind[j]) += real_data[j] * X(i);

	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  Y(imag_ind[j]) += T1(rzero, imag_data[j]) * X(i);
      }
  }

  
  /*** Symmetric complex sparse matrices ***/
  
  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    int i; long j;
    typename ClassComplexType<T1>::Treal rzero(0);
    T1 zero(0, 0);
    T1 temp;

    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    Y.Fill(0);
    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
            temp += real_data[j] * X(real_ind[j]);
            if (real_ind[j] != i)
              Y(real_ind[j]) += real_data[j] * X(i);
          }
        
        for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  {
            temp += T1(rzero, imag_data[j]) * X(imag_ind[j]);
            if (imag_ind[j] != i)
              Y(imag_ind[j]) += T1(rzero, imag_data[j]) * X(i);
          }
        
        Y(i) += temp;
      }
  }

  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    int i; long j;
    typename ClassComplexType<T1>::Treal rzero(0);
    T1 zero(0, 0);
    T1 temp;

    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    Y.Fill(0);
    for (i = 0; i<ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
            temp += real_data[j] * X(real_ind[j]);
            if (real_ind[j] != i)
              Y(real_ind[j]) += real_data[j] * X(i);
          }

        for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  {
            temp += T1(rzero, imag_data[j]) * X(imag_ind[j]);
            if (imag_ind[j] != i)
              Y(imag_ind[j]) += T1(rzero, imag_data[j]) * X(i);
          }
        
        Y(i) += temp;
      }
  }

  
  /*** Complex sparse matrices, *Trans ***/


  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
	MltVector(M, X, Y);
	return;
      }

    int i; long j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif

    typename ClassComplexType<T1>::Treal rzero(0);
    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    Y.Fill(0);

    if (Trans.Trans())
      {
	for (i = 0; i < ma; i++)
	  for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	    Y(real_ind[j]) += real_data[j] * X(i);
	
	for (i = 0; i < ma; i++)
	  for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	    Y(imag_ind[j]) += T1(rzero, imag_data[j]) * X(i);
      }
    else
      {
	for (i = 0; i < ma; i++)
	  for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	    Y(real_ind[j]) += real_data[j] * X(i);
	
	for (i = 0; i < ma; i++)
	  for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	    Y(imag_ind[j]) += T1(rzero, - imag_data[j]) * X(i);
      }
  }
  
  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, ColComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
	MltVector(M, X, Y);
	return;
      }

    int i; long j;

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif
    
    T4 temp;
    typename ClassComplexType<T1>::Treal rzero(0);
    
    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, ColComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, ColComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    if (Trans.Trans())
      {
	for (i = 0; i < M.GetN(); i++)
	  {
	    SetComplexZero(temp);
	    for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	      temp += real_data[j] * X(real_ind[j]);
	    
	    for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	      temp += T1(rzero, imag_data[j]) * X(imag_ind[j]);
	    
	    Y(i) = temp;
	  }
      }
    else
      {
	for (i = 0; i < M.GetN(); i++)
	  {
	    SetComplexZero(temp);
	    for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	      temp += real_data[j] * X(real_ind[j]);
	    
	    for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	      temp += T1(rzero, -imag_data[j]) * X(imag_ind[j]);
	    
	    Y(i) = temp;
	  }
      }
  }
  
  
  /*** Symmetric complex sparse matrices, *Trans ***/

  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    if (!Trans.ConjTrans())
      {
	MltVector(M, X, Y);
	return;
      }

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    int i; long j;
    T1 zero(0, 0);
    T1 temp;
    typename ClassComplexType<T1>::Treal rzero(0);

    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    Y.Fill(0);
    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
            temp += real_data[j] * X(real_ind[j]);
            if (real_ind[j] != i)
              Y(real_ind[j]) += real_data[j] * X(i);
          }
	
        for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  {
            temp += T1(rzero, - imag_data[j]) * X(imag_ind[j]);
            if (imag_ind[j] != i)
              Y(imag_ind[j]) += T1(rzero, - imag_data[j]) * X(i);
          }
        
	Y(i) += temp;
      }
  }
  
  
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    if (!Trans.ConjTrans())
      {
	MltVector(M, X, Y);
	return;
      }

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    int i; long j;
    T1 zero(0, 0);
    T1 temp;
    typename ClassComplexType<T1>::Treal rzero(0);
    
    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    Y.Fill(0);
    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
            temp += real_data[j] * X(real_ind[j]);
            if (real_ind[j] != i)
              Y(real_ind[j]) += real_data[j] * X(i);
          }
	
        for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  {
            temp += T1(rzero, - imag_data[j]) * X(imag_ind[j]);
            if (imag_ind[j] != i)
              Y(imag_ind[j]) += T1(rzero, - imag_data[j]) * X(i);
          }
        
	Y(i) += temp;
      }
  }
  
  
  /*** ArrayRowSymComplexSparse ***/


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, Symmetric,
		 ArrayRowSymComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    C.Fill(0);
    
    int m = A.GetM(), n, p;
    typename ClassComplexType<T1>::Treal val;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexReal(i, k);
	    val = A.ValueReal(i, k);
	    if (p == i)
	      C(i) += val * B(i);
	    else
	      {
		C(i) += val * B(p);
		C(p) += val * B(i);
	      }
	  }
	
	n = A.GetImagRowSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexImag(i, k);
	    val = A.ValueImag(i, k);
	    if (p == i)
	      C(i) += T1(0, val) * B(i);
	    else
	      {
		C(i) += T1(0, val) * B(p);
		C(p) += T1(0, val) * B(i);
	      }
	  }
      }
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		 ArrayRowSymComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    if (!Trans.ConjTrans())
      {
	MltVector(A, B, C);
	return;
      }

    C.Fill(0);
    int m = A.GetM(),n,p;
    typename ClassComplexType<T1>::Treal val;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexReal(i, k);
	    val = A.ValueReal(i, k);
	    if (p == i)
	      C(i) += val * B(i);
	    else
	      {
		C(i) += val * B(p);
		C(p) += val * B(i);
	      }
	  }
	n = A.GetImagRowSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexImag(i, k);
	    val = A.ValueImag(i, k);
	    if (p == i)
	      C(i) -= T1(0, val) * B(i);
	    else
	      {
		C(i) -= T1(0, val) * B(p);
		C(p) -= T1(0, val) * B(i);
	      }
	  }
      }
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, Symmetric,
		 ArrayColSymComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    C.Fill(0);
    int m = A.GetM(), n, p;
    typename ClassComplexType<T1>::Treal val;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealColumnSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexReal(i, k);
	    val = A.ValueReal(i, k);
	    if (p == i)
	      C(i) += val * B(i);
	    else
	      {
		C(i) += val * B(p);
		C(p) += val * B(i);
	      }
	  }
	
	n = A.GetImagColumnSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexImag(i, k);
	    val = A.ValueImag(i, k);
	    if (p == i)
	      C(i) += T1(0, val) * B(i);
	    else
	      {
		C(i) += T1(0, val) * B(p);
		C(p) += T1(0, val) * B(i);
	      }
	  }
      }
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		 ArrayColSymComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    if (!Trans.ConjTrans())
      {
	MltVector(A, B, C);
	return;
      }

    C.Fill(0);
        
    int m = A.GetM(), n, p;
    typename ClassComplexType<T1>::Treal val;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealColumnSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexReal(i, k);
	    val = A.ValueReal(i, k);
	    if (p == i)
	      C(i) += val * B(i);
	    else
	      {
		C(i) += val * B(p);
		C(p) += val * B(i);
	      }
	  }
	
	n = A.GetImagColumnSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexImag(i, k);
	    val = A.ValueImag(i, k);
	    if (p == i)
	      C(i) -= T1(0, val) * B(i);
	    else
	      {
		C(i) -= T1(0, val) * B(p);
		C(p) -= T1(0, val) * B(i);
	      }
	  }
      }
  }


  /*** ArrayRowComplexSparse ***/


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, General,
		 ArrayRowComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    int m = A.GetM(), n;
    T3 temp;
    for (int i = 0 ; i < m ; i++)
      {
	SetComplexZero(temp);
	n = A.GetRealRowSize(i);
	for (int k = 0; k < n ; k++)
	  temp += A.ValueReal(i, k)*B(A.IndexReal(i, k));
	
	n = A.GetImagRowSize(i);
	for (int k = 0; k < n ; k++)
	  temp += T1(0, A.ValueImag(i, k)) * B(A.IndexImag(i, k));
	
	C(i) = temp;
      }
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, General,
		 ArrayRowComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    if (Trans.NoTrans())
      {
	MltVector(A, B, C);
	return;
      }
    
    C.Fill(0);
        
    int m = A.GetM(), n, p;
    typename ClassComplexType<T1>::Treal val;
      
    if (Trans.Trans())
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(p) += val * B(i);
	      }
	    
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(p) += T1(0, val) * B(i);
	      }
	  }
      }
    else
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(p) += val * B(i);
	      }
	    
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(p) -= T1(0, val) * B(i);
	      }
	  }
      }
  }
  
  
  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const Matrix<T1, General,
		 ArrayColComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    C.Fill(0);
        
    int n, p;
    typename ClassComplexType<T1>::Treal val;
    for (int i = 0 ; i < A.GetN(); i++)
      {
	n = A.GetRealColumnSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexReal(i, k);
	    val = A.ValueReal(i, k);
	    C(p) += val * B(i);
	  }
	
	n = A.GetImagColumnSize(i);
	for (int k = 0; k < n ; k++)
	  {
	    p = A.IndexImag(i, k);
	    val = A.ValueImag(i, k);
	    C(p) += T1(0, val) * B(i);
	  }
      }
  }


  template<class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltVector(const SeldonTranspose& Trans, const Matrix<T1, General,
		 ArrayColComplexSparse, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& B,
		 Vector<T3, VectFull, Allocator3>& C)
  {
    if (Trans.NoTrans())
      {
	MltVector(A, B, C);
	return;
      }

    T3 temp; int n;
    
    if (Trans.Trans())
      {
	for (int i = 0; i < A.GetN(); i++)
	  {
	    SetComplexZero(temp);
	    n = A.GetRealColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      temp += A.ValueReal(i, k) * B(A.IndexReal(i, k));
	    
	    n = A.GetImagColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      temp += T1(0, A.ValueImag(i, k)) * B(A.IndexImag(i, k));
	    
	    C(i) = temp;
	  }
      }
    else
      {
	for (int i = 0; i < A.GetN(); i++)
	  {
	    SetComplexZero(temp);
	    n = A.GetRealColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      temp += A.ValueReal(i, k) * B(A.IndexReal(i, k));
	    
	    n = A.GetImagColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      temp -= T1(0, A.ValueImag(i, k)) * B(A.IndexImag(i, k));
	    
	    C(i) = temp;
	  }
      }
  }


  /****************
   * MltAddVector *
   ****************/
  
  
  /*** Complex sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i; long j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    typename ClassComplexType<T1>::Treal rzero(0);
    T1 zero(0, 0);
    T1 temp;

    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * X(real_ind[j]);
	Y(i) += alpha * temp;
      }

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += T1(rzero, imag_data[j]) * X(imag_ind[j]);
	Y(i) += alpha * temp;
      }
  }
  
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, ColComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i; long j;

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    typename ClassComplexType<T1>::Treal rzero(0);
    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, ColComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, ColComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < M.GetN(); i++)
      {
        for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
          Y(real_ind[j]) += alpha * real_data[j] * X(i);

	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  Y(imag_ind[j]) += alpha * T1(rzero, imag_data[j]) * X(i);
      }
  }

  
  /*** Symmetric complex sparse matrices ***/
  
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i; long j;
    typename ClassComplexType<T1>::Treal rzero(0);
    T1 zero(0, 0);
    T1 temp;

    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
            temp += real_data[j] * X(real_ind[j]);
            if (real_ind[j] != i)
              Y(real_ind[j]) += alpha * real_data[j] * X(i);
          }
        
        for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  {
            temp += T1(rzero, imag_data[j]) * X(imag_ind[j]);
            if (imag_ind[j] != i)
              Y(imag_ind[j]) += alpha * T1(rzero, imag_data[j]) * X(i);
          }
        
        Y(i) += alpha * temp;
      }
  }

  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i; long j;
    typename ClassComplexType<T1>::Treal rzero(0);
    T1 zero(0, 0);
    T1 temp;

    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i<ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
            temp += real_data[j] * X(real_ind[j]);
            if (real_ind[j] != i)
              Y(real_ind[j]) += alpha * real_data[j] * X(i);
          }

        for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  {
            temp += T1(rzero, imag_data[j]) * X(imag_ind[j]);
            if (imag_ind[j] != i)
              Y(imag_ind[j]) += alpha * T1(rzero, imag_data[j]) * X(i);
          }
        
        Y(i) += alpha * temp;
      }
  }

  
  /*** Complex sparse matrices, *Trans ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
	MltAddVector(alpha, M, X, beta, Y);
	return;
      }
    
    int i; long j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    typename ClassComplexType<T1>::Treal rzero(0);
    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    if (Trans.Trans())
      {
	for (i = 0; i < ma; i++)
	  for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	    Y(real_ind[j]) += alpha * real_data[j] * X(i);
	
	for (i = 0; i < ma; i++)
	  for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	    Y(imag_ind[j]) += alpha * T1(rzero, imag_data[j]) * X(i);
      }
    else
      {
	for (i = 0; i < ma; i++)
	  for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	    Y(real_ind[j]) += alpha * real_data[j] * X(i);
	
	for (i = 0; i < ma; i++)
	  for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	    Y(imag_ind[j]) += alpha * T1(rzero, - imag_data[j]) * X(i);
      }
  }
  
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, ColComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
	MltAddVector(alpha, M, X, beta, Y);
	return;
      }

    int i; long j;

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);
    
    T4 temp;
    typename ClassComplexType<T1>::Treal rzero(0);
    
    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, ColComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, ColComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();
    
    if (Trans.Trans())
      {
	for (i = 0; i < M.GetN(); i++)
	  {
	    SetComplexZero(temp);
	    for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	      temp += real_data[j] * X(real_ind[j]);
	    
	    for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	      temp += T1(rzero, imag_data[j]) * X(imag_ind[j]);
	    
	    Y(i) += alpha * temp;
	  }
      }
    else
      {
	for (i = 0; i < M.GetN(); i++)
	  {
	    SetComplexZero(temp);
	    for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	      temp += real_data[j] * X(real_ind[j]);
	    
	    for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	      temp += T1(rzero, -imag_data[j]) * X(imag_ind[j]);
	    
	    Y(i) += alpha * temp;
	  }
      }
  }
  
  
  /*** Symmetric complex sparse matrices, *Trans ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    if (!Trans.ConjTrans())
      {
	MltAddVector(alpha, M, X, beta, Y);
	return;
      }

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i; long j;
    T1 zero(0, 0);
    T1 temp;
    typename ClassComplexType<T1>::Treal rzero(0);

    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
            temp += real_data[j] * X(real_ind[j]);
            if (real_ind[j] != i)
              Y(real_ind[j]) += alpha * real_data[j] * X(i);
          }
	
        for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  {
            temp += T1(rzero, - imag_data[j]) * X(imag_ind[j]);
            if (imag_ind[j] != i)
              Y(imag_ind[j]) += alpha * T1(rzero, - imag_data[j]) * X(i);
          }
        
	Y(i) += alpha * temp;
      }
  }
  
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    if (!Trans.ConjTrans())
      {
	MltAddVector(alpha, M, X, beta, Y);
	return;
      }

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i; long j;
    T1 zero(0, 0);
    T1 temp;
    typename ClassComplexType<T1>::Treal rzero(0);
    
    long* real_ptr = M.GetRealPtr();
    long* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  {
            temp += real_data[j] * X(real_ind[j]);
            if (real_ind[j] != i)
              Y(real_ind[j]) += alpha * real_data[j] * X(i);
          }
	
        for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  {
            temp += T1(rzero, - imag_data[j]) * X(imag_ind[j]);
            if (imag_ind[j] != i)
              Y(imag_ind[j]) += alpha * T1(rzero, - imag_data[j]) * X(i);
          }
        
	Y(i) += alpha * temp;
      }
  }
  
  
  /*** ArrayRowSymComplexSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric,
		    ArrayRowSymComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    typename ClassComplexType<T1>::Treal val;
    T3 val_cplx;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		if (p == i)
		  C(i) += T1(0, val) * B(i);
		else
		  {
		    C(i) += T1(0, val) * B(p);
		    C(p) += T1(0, val) * B(i);
		  }
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * T1(0, A.ValueImag(i, k));
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		    ArrayRowSymComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    if (!Trans.ConjTrans())
      {
	MltAddVector(alpha, A, B, beta, C);
	return;
      }

    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);
    
    int m = A.GetM(),n,p;
    typename ClassComplexType<T1>::Treal val;
    T3 val_cplx;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		if (p == i)
		  C(i) -= T1(0, val) * B(i);
		else
		  {
		    C(i) -= T1(0, val) * B(p);
		    C(p) -= T1(0, val) * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * T1(0, A.ValueImag(i, k));
		if (p == i)
		  C(i) -= val_cplx * B(i);
		else
		  {
		    C(i) -= val_cplx * B(p);
		    C(p) -= val_cplx * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric,
		    ArrayColSymComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    typename ClassComplexType<T1>::Treal val;
    T3 val_cplx;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	    n = A.GetImagColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		if (p == i)
		  C(i) += T1(0, val) * B(i);
		else
		  {
		    C(i) += T1(0, val) * B(p);
		    C(p) += T1(0, val) * B(i);
		  }
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	    n = A.GetImagColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * T1(0, A.ValueImag(i, k));
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, Symmetric,
		    ArrayColSymComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    if (!Trans.ConjTrans())
      {
	MltAddVector(alpha, A, B, beta, C);
	return;
      }

    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);
    
    int m = A.GetM(),n,p;
    typename ClassComplexType<T1>::Treal val;
    T3 val_cplx;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	    n = A.GetImagColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		if (p == i)
		  C(i) -= T1(0, val) * B(i);
		else
		  {
		    C(i) -= T1(0, val) * B(p);
		    C(p) -= T1(0, val) * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	    n = A.GetImagColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * T1(0, A.ValueImag(i, k));
		if (p == i)
		  C(i) -= val_cplx * B(i);
		else
		  {
		    C(i) -= val_cplx * B(p);
		    C(p) -= val_cplx * B(i);
		  }
	      }
	  }
      }
  }


  /*** ArrayRowComplexSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, General,
		    ArrayRowComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);
    
    int m = A.GetM(), n, p;
    typename ClassComplexType<T1>::Treal val;
    T3 val_cplx;
    if (alpha == one)
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(i) += val * B(p);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(i) += T1(0, val) * B(p);
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		C(i) += val_cplx * B(p);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * T1(0, A.ValueImag(i, k));
		C(i) += val_cplx * B(p);
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, General,
		    ArrayRowComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    if (Trans.NoTrans())
      {
	MltAddVector(alpha, A, B, beta, C);
	return;
      }

    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);
    
    int m = A.GetM(),n,p;
    typename ClassComplexType<T1>::Treal val;
    T3 val_cplx;

    if (Trans.Trans())
      {
	if (alpha == one)
	  {
	    for (int i = 0 ; i < m ; i++)
	      {
		n = A.GetRealRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexReal(i, k);
		    val = A.ValueReal(i, k);
		    C(p) += val * B(i);
		  }
		n = A.GetImagRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexImag(i, k);
		    val = A.ValueImag(i, k);
		    C(p) += T1(0, val) * B(i);
		  }
	      }
	  }
	else
	  {
	    // alpha different from 1
	    for (int i = 0 ; i < m ; i++)
	      {
		n = A.GetRealRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexReal(i, k);
		    val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		    C(p) += val_cplx * B(i);
		  }
		
		n = A.GetImagRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexImag(i, k);
		    val_cplx = alpha * T1(0, A.ValueImag(i, k));
		    C(p) += val_cplx * B(i);
		  }
	      }
	  }
      }
    else
      {
	if (alpha == one)
	  {
	    for (int i = 0 ; i < m ; i++)
	      {
		n = A.GetRealRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexReal(i, k);
		    val = A.ValueReal(i, k);
		    C(p) += val * B(i);
		  }
		n = A.GetImagRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexImag(i, k);
		    val = A.ValueImag(i, k);
		    C(p) -= T1(0, val) * B(i);
		  }
	      }
	  }
	else
	  {
	    // alpha different from 1
	    for (int i = 0 ; i < m ; i++)
	      {
		n = A.GetRealRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexReal(i, k);
		    val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		    C(p) += val_cplx * B(i);
		  }
		n = A.GetImagRowSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexImag(i, k);
		    val_cplx = alpha * T1(0, -A.ValueImag(i, k));
		    C(p) += val_cplx * B(i);
		  }
	      }
	  }
      }
  }
  
  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, General,
		    ArrayColComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);
    
    int n, p;
    typename ClassComplexType<T1>::Treal val;
    T3 val_cplx;
    if (alpha == one)
      {
	for (int i = 0 ; i < A.GetN(); i++)
	  {
	    n = A.GetRealColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(p) += val * B(i);
	      }
	    n = A.GetImagColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(p) += T1(0, val) * B(i);
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0; i < A.GetN(); i++)
	  {
	    n = A.GetRealColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		C(p) += val_cplx * B(i);
	      }
	    n = A.GetImagColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * T1(0, A.ValueImag(i, k));
		C(p) += val_cplx * B(i);
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans, const Matrix<T1, General,
		    ArrayColComplexSparse, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& B,
		    const T4& beta,
		    Vector<T3, VectFull, Allocator3>& C)
  {
    if (Trans.NoTrans())
      {
	MltAddVector(alpha, A, B, beta, C);
	return;
      }

    T4 zero; T0 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    if (beta == zero)
      C.Fill(0);
    else
      Mlt(beta, C);
    
    int n, p;
    typename ClassComplexType<T1>::Treal val;
    T3 val_cplx;
    
    if (Trans.Trans())
      {
	if (alpha == one)
	  {
	    for (int i = 0; i < A.GetN(); i++)
	      {
		n = A.GetRealColumnSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexReal(i, k);
		    val = A.ValueReal(i, k);
		    C(i) += val * B(p);
		  }
		
		n = A.GetImagColumnSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexImag(i, k);
		    val = A.ValueImag(i, k);
		    C(i) += T1(0, val) * B(p);
		  }
	      }
	  }
	else
	  {
	    // alpha different from 1
	    for (int i = 0; i < A.GetN(); i++)
	      {
		n = A.GetRealColumnSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexReal(i, k);
		    val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		    C(i) += val_cplx * B(p);
		  }
		
		n = A.GetImagColumnSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexImag(i, k);
		    val_cplx = alpha * T1(0, A.ValueImag(i, k));
		    C(i) += val_cplx * B(p);
		  }
	      }
	  }
      }
    else
      {
	if (alpha == one)
	  {
	    for (int i = 0; i < A.GetN(); i++)
	      {
		n = A.GetRealColumnSize(i);
		for (int k = 0; k < n; k++)
		  {
		    p = A.IndexReal(i, k);
		    val = A.ValueReal(i, k);
		    C(i) += val * B(p);
		  }
		n = A.GetImagColumnSize(i);
		for (int k = 0; k < n; k++)
		  {
		    p = A.IndexImag(i, k);
		    val = A.ValueImag(i, k);
		    C(i) -= T1(0, val) * B(p);
		  }
	      }
	  }
	else
	  {
	    // alpha different from 1
	    for (int i = 0 ; i < A.GetN(); i++)
	      {
		n = A.GetRealColumnSize(i);
		for (int k = 0; k < n ; k++)
		  {
		    p = A.IndexReal(i, k);
		    val_cplx = alpha * T1(A.ValueReal(i, k), 0);
		    C(i) += val_cplx * B(p);
		  }
		n = A.GetImagColumnSize(i);
		for (int k = 0; k < n; k++)
		  {
		    p = A.IndexImag(i, k);
		    val_cplx = alpha * T1(0, -A.ValueImag(i, k));
		    C(i) += val_cplx * B(p);
		  }
	      }
	  }
      }
  }
  
  
  /////////
  // SOR //
  
  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    complex<T1> temp, zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    long* ptr_real = A.GetRealPtr();
    long* ptr_imag = A.GetImagPtr();
    int* ind_real = A.GetRealInd();
    int* ind_imag = A.GetImagInd();
    
    typename ClassComplexType<T0>::Treal* data_real = A.GetRealData();
    typename ClassComplexType<T0>::Treal* data_imag = A.GetImagData();
    
    T0 ajj;

    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = 0; j < ma; j++)
	  {
	    temp = zero;
            ajj = zero;
            for (long k = ptr_real[j]; k < ptr_real[j+1]; k++)
	      {
                if (ind_real[k] != j)
                  temp += data_real[k] * X(ind_real[k]);
                else
                  SetComplexReal(data_real[k], ajj);
              }

	    for (long k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
              {
                if (ind_imag[k] != j)
                  temp += complex<T1>(0, data_imag[k]) * X(ind_imag[k]);
                else
                  ajj += complex<T1>(0, data_imag[k]);
              }
            
#ifdef SELDON_CHECK_BOUNDS
            if (ajj == zero)
              throw WrongArgument("SOR", "Matrix must contain"
                                  " a non-null diagonal");
#endif
            
	    X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
	  }

    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = ma-1; j >= 0; j--)
	  {
	    temp = zero;
            ajj = zero;
            for (long k = ptr_real[j]; k < ptr_real[j+1]; k++)
	      {
                if (ind_real[k] != j)
                  temp += data_real[k] * X(ind_real[k]);
                else
                  SetComplexReal(data_real[k], ajj);
              }

	    for (long k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
              {
                if (ind_imag[k] != j)
                  temp += complex<T1>(0, data_imag[k]) * X(ind_imag[k]);
                else
                  ajj += complex<T1>(0, data_imag[k]);
              }
            
	    X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
	  }
  }


  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    complex<T1> temp, zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    T0 ajj;

    if (type_ssor%2 == 0)
      for (int i = 0; i < iter; i++)
        for (int j = 0; j < ma; j++)
          {
            temp = zero;
            ajj = zero;
            for (int k = 0; k < A.GetRealRowSize(j); k++)
              {
                if (A.IndexReal(j, k) != j)
                  temp += A.ValueReal(j,k) * X(A.IndexReal(j,k));
                else
                  ajj += A.ValueReal(j, k);
              }
            
            for (int k = 0; k < A.GetImagRowSize(j); k++)
              {
                if (A.IndexImag(j, k) != j)
                  temp += T0(0, A.ValueImag(j,k))
                    * X(A.IndexImag(j, k));
                else
                  ajj += T0(0, A.ValueImag(j,k));
              }
            
#ifdef SELDON_CHECK_BOUNDS
            if (ajj == zero)
              throw WrongArgument("SOR", "Matrix must contain"
                                  " a non-null diagonal");
#endif
            
            X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
          }

    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = ma-1; j >= 0; j--)
	  {
	      temp = zero;
	      ajj = zero;
	      for (int k = 0; k < A.GetRealRowSize(j); k++)
		{
		  if (A.IndexReal(j, k) != j)
                    temp += A.ValueReal(j,k) * X(A.IndexReal(j,k));
		  else
		    ajj += A.ValueReal(j, k);
		}
              
	      for (int k = 0; k < A.GetImagRowSize(j); k++)
		{
		  if (A.IndexImag(j, k) != j)
                    temp += T0(0, A.ValueImag(j,k))
                      * X(A.IndexImag(j, k));
		  else
		    ajj += T0(0, A.ValueImag(j,k));
		}
              
              X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
	  }
  }

  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor)
  {
    complex<T1> zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    long* ptr_real = A.GetRealPtr();
    long* ind_real = A.GetRealInd();
    typename ClassComplexType<T0>::Treal* data_real = A.GetRealData();

    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T0>::Treal* data_imag = A.GetImagData();

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A

    // Forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T0 ajj;
    T3 coef = (one - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (U + (1-omega)/omega D) X + B	      
	  for (int j = 0; j < ma; j++)
	    {
              long kr = ptr_real[j];
              while ( (kr < ptr_real[j+1]) && (ind_real[kr] < j))
                {
                  X(ind_real[kr]) -= data_real[kr]*X(j);
                  kr++;
                }
              
              if ( (kr < ptr_real[j+1]) && (ind_real[kr] == j))
                SetComplexReal(data_real[kr], ajj);

              long ki = ptr_imag[j];
              while ( (ki < ptr_imag[j+1]) && (ind_imag[ki] < j))
                {
                  X(ind_imag[ki]) -= T0(0, data_imag[ki])*X(j);
                  ki++;
                }

              if ( (ki < ptr_imag[j+1]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki]);
              
#ifdef SELDON_CHECK_BOUNDS
              if (ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
              X(j) = B(j) + coef * ajj * X(j);
	    }

          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              long kr = ptr_real[j];
              while ( (kr < ptr_real[j+1]) && (ind_real[kr] < j))
                kr++;
                             
              if ( (kr < ptr_real[j+1]) && (ind_real[kr] == j))
                {
                  SetComplexReal(data_real[kr], ajj);
                  kr++;
                }

              long ki = ptr_imag[j];
              while ( (ki < ptr_imag[j+1]) && (ind_imag[ki] < j))
                ki++;
              
              if ( (ki < ptr_imag[j+1]) && (ind_imag[ki] == j))
                {
                  ajj += T0(0, data_imag[ki]);
                  ki++;
                }
              
              X(j) *= omega/ajj;
              while (kr < ptr_real[j+1])
                {
                  X(ind_real[kr]) -= data_real[kr]*X(j);
                  kr++;
                }

              while (ki < ptr_imag[j+1])
                {
                  X(ind_imag[ki]) -= T0(0, data_imag[ki])*X(j);
                  ki++;
                }
            }
	}
    
    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (L + (1-omega)/omega D) X + B.
          for (int j = ma-1; j >= 0; j--)
            {
              long kr = ptr_real[j+1]-1;
              while ( (kr >= ptr_real[j]) && (ind_real[kr] > j) )
                {
                  X(ind_real[kr]) -= data_real[kr]*X(j);
                  kr--;
                }
              
              if ( (kr >= ptr_real[j]) && (ind_real[kr] == j))
                SetComplexReal(data_real[kr], ajj);
              
              long ki = ptr_imag[j+1]-1;
              while ( (ki >= ptr_imag[j]) && (ind_imag[ki] > j) )
                {
                  X(ind_imag[ki]) -= T0(0, data_imag[ki])*X(j);
                  ki--;
                }
              
              if ( (ki >= ptr_imag[j]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki]);
              
#ifdef SELDON_CHECK_BOUNDS
              if (ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
              X(j) = B(j) + coef * ajj * X(j);
            }
                    
          // Then we solve (D/omega - U) X = X.
          for (int j = ma-1; j >= 0; j--)
            {
              long kr = ptr_real[j+1]-1;
              while ( (kr >= ptr_real[j]) && (ind_real[kr] > j) )
                kr--;
              
              if ( (kr >= ptr_real[j]) && (ind_real[kr] == j))
                {
                  ajj = T0(data_real[kr], 0);
                  kr--;
                }
              
              long ki = ptr_imag[j+1]-1;
              while ( (ki >= ptr_imag[j]) && (ind_imag[ki] > j) )
                ki--;
              
              if ( (ki >= ptr_imag[j]) && (ind_imag[ki] == j))
                {
                  ajj += T0(0, data_imag[ki]);
                  ki--;
                }

              X(j) *= omega/ajj;
              while (kr >= ptr_real[j])
                {
                  X(ind_real[kr]) -= data_real[kr]*X(j);
                  kr--;
                }

              while (ki >= ptr_imag[j])
                {
                  X(ind_imag[ki]) -= T0(0, data_imag[ki])*X(j);
                  ki--;
                }
            }
        }
  }

  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    complex<T1> zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A

    // Forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T0 ajj;
    T3 coef = (one - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (U + (1-omega)/omega D) X + B	      
	  for (int j = 0; j < ma; j++)
	    {
              long kr = 0;
              while ( (kr < A.GetRealColumnSize(j)) && (A.IndexReal(j, kr) < j))
                {
                  X(A.IndexReal(j, kr)) -= A.ValueReal(j, kr) * X(j);
                  kr++;
                }
              
              if ( (kr < A.GetRealColumnSize(j)) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr), 0);
              
              long ki = 0;
              while ( (ki < A.GetImagColumnSize(j)) && (A.IndexImag(j, ki) < j))
                {
                  X(A.IndexImag(j, ki)) -= T0(0, A.ValueImag(j, ki))*X(j);
                  ki++;
                }

              if ( (ki < A.GetImagColumnSize(j)) && (A.IndexImag(j, ki) == j) )
                ajj += T0(0, A.ValueImag(j, ki));
              
#ifdef SELDON_CHECK_BOUNDS
              if (ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
              X(j) = B(j) + coef * ajj * X(j);
	    }

          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              long kr = 0;
              while ( (kr < A.GetRealColumnSize(j)) && (A.IndexReal(j, kr) < j))
                kr++;
                             
              if ( (kr < A.GetRealColumnSize(j)) && (A.IndexReal(j, kr) == j) )
                {
                  ajj = T0(A.ValueReal(j, kr), 0);
                  kr++;
                }

              long ki = 0;
              while ( (ki < A.GetImagColumnSize(j)) && (A.IndexImag(j, ki) < j))
                ki++;
              
              if ( (ki < A.GetImagColumnSize(j)) && (A.IndexImag(j, ki) == j))
                {
                  ajj += T0(0, A.ValueImag(j, ki));
                  ki++;
                }
              
              X(j) *= omega/ajj;
              while (kr < A.GetRealColumnSize(j))
                {
                  X(A.IndexReal(j, kr)) -= A.ValueReal(j, kr)*X(j);
                  kr++;
                }

              while (ki < A.GetImagColumnSize(j))
                {
                  X(A.IndexImag(j, ki)) -= T0(0, A.ValueImag(j, ki))*X(j);
                  ki++;
                }
            }
	}
    
    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (L + (1-omega)/omega D) X + B.
          for (int j = ma-1; j >= 0; j--)
            {
              long kr = A.GetRealColumnSize(j)-1;
              while ( (kr >= 0) && (A.IndexReal(j, kr) > j) )
                {
                  X(A.IndexReal(j, kr)) -= A.ValueReal(j, kr)*X(j);
                  kr--;
                }
              
              if ( (kr >= 0) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr), 0);
              
              long ki = A.GetImagColumnSize(j)-1;
              while ( (ki >= 0) && (A.IndexImag(j, ki) > j) )
                {
                  X(A.IndexImag(j, ki)) -= T0(0, A.ValueImag(j, ki))*X(j);
                  ki--;
                }
              
              if ( (ki >= 0) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki));
              
#ifdef SELDON_CHECK_BOUNDS
              if (ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
              X(j) = B(j) + coef * ajj * X(j);
            }
                    
          // Then we solve (D/omega - U) X = X.
          for (int j = ma-1; j >= 0; j--)
            {
              long kr = A.GetRealColumnSize(j)-1;
              while ( (kr >= 0) && (A.IndexReal(j, kr) > j) )
                kr--;
              
              if ( (kr >= 0) && (A.IndexReal(j, kr) == j))
                {
                  ajj = T0(A.ValueReal(j, kr), 0);
                  kr--;
                }
              
              long ki = A.GetImagColumnSize(j)-1;
              while ( (ki >= 0) && (A.IndexImag(j, ki) > j) )
                ki--;
              
              if ( (ki >= 0) && (A.IndexImag(j, ki) == j))
                {
                  ajj += T0(0, A.ValueImag(j, ki));
                  ki--;
                }

              X(j) *= omega/ajj;
              while (kr >= 0)
                {
                  X(A.IndexReal(j, kr)) -= A.ValueReal(j, kr)*X(j);
                  kr--;
                }

              while (ki >= 0)
                {
                  X(A.IndexImag(j, ki)) -= T0(0, A.ValueImag(j, ki))*X(j);
                  ki--;
                }
            }
        }
  }

  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    complex<T1> temp, zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T0>::Treal* data_real = A.GetRealData();

    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T0>::Treal* data_imag = A.GetImagData();

    T0 ajj;

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A, A is symmetric, so L = U^t

    // Forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T3 coef = (T3(1) - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we do X = (U + (1-omega)/omega D) X + B
	  for (int j = 0; j < ma; j++)
	    {
	      temp = zero;
              ajj = zero;
              long kr = ptr_real[j];
              if ((kr < ptr_real[j+1]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr++], 0);

              long ki = ptr_imag[j];
              if ((ki < ptr_imag[j+1]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki++]);
              
#ifdef SELDON_CHECK_BOUNDS
              if ( ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
	      for (long k = kr; k < ptr_real[j+1]; k++)
                temp += data_real[k] * X(ind_real[k]);

	      for (long k = ki; k < ptr_imag[j+1]; k++)
                temp += T0(0, data_imag[k]) * X(ind_imag[k]);
              
	      X(j) = coef * ajj * X(j) + B(j) - temp;
	    }

          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              ajj = zero;
              long kr = ptr_real[j];
              if ((kr < ptr_real[j+1]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr++], 0);

              long ki = ptr_imag[j];
              if ((ki < ptr_imag[j+1]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki++]);
              
              X(j) *= omega / ajj;
	      for (long k = kr; k < ptr_real[j+1]; k++)
                X(ind_real[k]) -= data_real[k]*X(j);

	      for (long k = ki; k < ptr_imag[j+1]; k++)
                X(ind_imag[k]) -= T0(0, data_imag[k])*X(j);
	    }
	}


    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (L + (1-omega)/omega D) X + B.
          for (int j = ma-1; j >= 0; j--)
	    {
              ajj = zero;
              long kr = ptr_real[j];
              if ((kr < ptr_real[j+1]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr++], 0);

              long ki = ptr_imag[j];
              if ((ki < ptr_imag[j+1]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki++]);
              
	      for (long k = kr; k < ptr_real[j+1]; k++)
		X(ind_real[k]) -= data_real[k]*X(j);

	      for (long k = ki; k < ptr_imag[j+1]; k++)
		X(ind_imag[k]) -= T0(0, data_imag[k])*X(j);
              
              X(j) = B(j) + coef * ajj * X(j);
	    }
          
          // Then we solve (D/omega - U) X = X.
	  for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              ajj = zero;
              long kr = ptr_real[j];
              if ((kr < ptr_real[j+1]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr++], 0);

              long ki = ptr_imag[j];
              if ((ki < ptr_imag[j+1]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki++]);
         
	      for (long k = kr; k < ptr_real[j+1]; k++)
                temp += data_real[k]*X(ind_real[k]);

	      for (long k = ki; k < ptr_imag[j+1]; k++)
                temp += T0(0, data_imag[k])*X(ind_imag[k]);
	      
              X(j) = (X(j) - temp) * omega / ajj;
	    }
	}
  }
  
  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    complex<T1> temp, zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    T0 ajj;

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A, A is symmetric, so L = U^t

    // Forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T3 coef = (one - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we do X = (U + (1-omega)/omega D) X + B
	  for (int j = 0; j < ma; j++)
	    {
	      temp = zero;
              ajj = zero;
              int kr = 0;
              if ((kr < A.GetRealRowSize(j)) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr++), 0);
              
              int ki = 0;
              if ((ki < A.GetImagRowSize(j)) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki++));
              
#ifdef SELDON_CHECK_BOUNDS
              if ( ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
	      for (int k = kr; k < A.GetRealRowSize(j); k++)
                temp += A.ValueReal(j, k) * X(A.IndexReal(j, k));

	      for (int k = ki; k < A.GetImagRowSize(j); k++)
                temp += T0(0, A.ValueImag(j, k)) * X(A.IndexImag(j, k));
              
	      X(j) = coef * ajj * X(j) + B(j) - temp;                            
	    }

          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              ajj = zero;
              int kr = 0;
              if ((kr < A.GetRealRowSize(j)) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr++), 0);

              int ki = 0;
              if ((ki < A.GetImagRowSize(j)) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki++));
              
              X(j) *= omega / ajj;
	      for (int k = kr; k < A.GetRealRowSize(j); k++)
                X(A.IndexReal(j, k)) -= A.ValueReal(j, k)*X(j);

	      for (int k = ki; k < A.GetImagRowSize(j); k++)
                X(A.IndexImag(j, k)) -= T0(0, A.ValueImag(j, k))*X(j);
	    }
	}


    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (L + (1-omega)/omega D) X + B.
          for (int j = ma-1; j >= 0; j--)
	    {
              ajj = zero;
              int kr = 0;
              if ((kr < A.GetRealRowSize(j)) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr++), 0);

              int ki = 0;
              if ((ki < A.GetImagRowSize(j)) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki++));
              
	      for (int k = kr; k < A.GetRealRowSize(j); k++)
		X(A.IndexReal(j, k)) -= A.ValueReal(j, k)*X(j);

	      for (int k = ki; k < A.GetImagRowSize(j); k++)
		X(A.IndexImag(j, k)) -= T0(0, A.ValueImag(j, k))*X(j);
              
              X(j) = B(j) + coef * ajj * X(j);
	    }
          
          // Then we solve (D/omega - U) X = X.
	  for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              ajj = zero;
              int kr = 0;
              if ((kr < A.GetRealRowSize(j)) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr++), 0);

              int ki = 0;
              if ((ki < A.GetImagRowSize(j)) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki++));
         
	      for (int k = kr; k < A.GetRealRowSize(j); k++)
                temp += A.ValueReal(j, k)*X(A.IndexReal(j, k));

	      for (int k = ki; k < A.GetImagRowSize(j); k++)
                temp += T0(0, A.ValueImag(j, k))*X(A.IndexImag(j, k));
	      
              X(j) = (X(j) - temp) * omega / ajj;
	    }
	}
  }

  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    complex<T1> temp, zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif
    
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T0>::Treal* data_real = A.GetRealData();
    
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T0>::Treal* data_imag = A.GetImagData();
    
    T0 ajj;

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A, A is symmetric, so L = U^t

    // Forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T3 coef = (one - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we do X = (U + (1-omega)/omega D) X + B
	  for (int j = 0; j < ma; j++)
	    {
              ajj = zero;
              long kr = ptr_real[j+1]-1;
              if ((kr >= ptr_real[j]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr--], 0);
              
              long ki = ptr_imag[j+1]-1;
              if ((ki >= ptr_imag[j]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki--]);
              
#ifdef SELDON_CHECK_BOUNDS
              if ( ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif

	      for (long k = ptr_real[j]; k <= kr; k++)
                X(ind_real[k]) -= data_real[k] * X(j);

	      for (long k = ptr_imag[j]; k <= ki; k++)
                X(ind_imag[k]) -= T0(0, data_imag[k]) * X(j);
              
	      X(j) = coef * ajj * X(j) + B(j);
	    }
          
          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              ajj = zero;
              long kr = ptr_real[j+1]-1;
              if ((kr >= ptr_real[j]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr--], 0);

              long ki = ptr_imag[j+1]-1;
              if ((ki >= ptr_imag[j]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki--]);
              
	      for (long k = ptr_real[j]; k <= kr; k++)
                X(j) -= data_real[k] * X(ind_real[k]);

	      for (long k = ptr_imag[j]; k <= ki; k++)
                X(j) -= T0(0, data_imag[k]) * X(ind_imag[k]);
              
              X(j) *= omega / ajj;
	    }
	}


    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (L + (1-omega)/omega D) X + B.
          for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              ajj = zero;
              long kr = ptr_real[j+1]-1;
              if ((kr >= ptr_real[j]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr--], 0);

              long ki = ptr_imag[j+1]-1;
              if ((ki >= ptr_imag[j]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki--]);
              
              for (long k = ptr_real[j]; k <= kr; k++)
		temp -= data_real[k] * X(ind_real[k]);

              for (long k = ptr_imag[j]; k <= ki; k++)
		temp -= T0(0, data_imag[k]) * X(ind_imag[k]);
              
              X(j) = B(j) + coef * ajj * X(j) + temp;
	    }
          
          // Then we solve (D/omega - U) X = X.
	  for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              ajj = zero;
              long kr = ptr_real[j+1]-1;
              if ((kr >= ptr_real[j]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr--], 0);

              long ki = ptr_imag[j+1]-1;
              if ((ki >= ptr_imag[j]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki--]);
              
              X(j) *= omega / ajj;
	      for (long k = ptr_real[j]; k <= kr; k++)
                X(ind_real[k]) -= data_real[k] * X(j);

	      for (long k = ptr_imag[j]; k <= ki; k++)
                X(ind_imag[k]) -= T0(0, data_imag[k]) * X(j);
	    }
	}
  }


  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    complex<T1> temp, zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    T0 ajj;

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A, A is symmetric, so L = U^t

    // Forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T3 coef = (T3(1) - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we do X = (U + (1-omega)/omega D) X + B
	  for (int j = 0; j < ma; j++)
	    {
              ajj = zero;
              int kr = A.GetRealColumnSize(j)-1;
              if ((kr >= 0) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr--), 0);

              int ki = A.GetImagColumnSize(j)-1;
              if ((ki >= 0) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki--));
              
#ifdef SELDON_CHECK_BOUNDS
              if ( ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif

	      for (int k = 0; k <= kr; k++)
                X(A.IndexReal(j, k)) -= A.ValueReal(j, k) * X(j);

	      for (int k = 0; k <= ki; k++)
                X(A.IndexImag(j, k)) -= T0(0, A.ValueImag(j, k)) * X(j);
              
	      X(j) = coef * ajj * X(j) + B(j);
	    }
          
          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              ajj = zero;
              int kr = A.GetRealColumnSize(j)-1;
              if ((kr >= 0) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr--), 0);

              int ki = A.GetImagColumnSize(j)-1;
              if ((ki >= 0) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki--));
              
	      for (int k = 0; k <= kr; k++)
                X(j) -= A.ValueReal(j, k) * X(A.IndexReal(j, k));

	      for (int k = 0; k <= ki; k++)
                X(j) -= T0(0, A.ValueImag(j, k)) * X(A.IndexImag(j, k));
              
              X(j) *= omega / ajj;
	    }
	}


    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (L + (1-omega)/omega D) X + B.
          for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              ajj = zero;
              int kr = A.GetRealColumnSize(j)-1;
              if ((kr >= 0) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr--), 0);

              int ki = A.GetImagColumnSize(j)-1;
              if ((ki >= 0) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki--));
              
              for (int k = 0; k <= kr; k++)
		temp -= A.ValueReal(j, k) * X(A.IndexReal(j, k));

              for (int k = 0; k <= ki; k++)
		temp -= T0(0, A.ValueImag(j, k)) * X(A.IndexImag(j, k));
              
              X(j) = B(j) + coef * ajj * X(j) + temp;
	    }
          
          // Then we solve (D/omega - U) X = X.
	  for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              ajj = zero;
              int kr = A.GetRealColumnSize(j)-1;
              if ((kr >= 0) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr--), 0);

              int ki = A.GetImagColumnSize(j)-1;
              if ((ki >= 0) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki--));
              
              X(j) *= omega / ajj;
	      for (int k = 0; k <= kr; k++)
                X(A.IndexReal(j, k)) -= A.ValueReal(j, k) * X(j);

	      for (int k = 0; k <= ki; k++)
                X(A.IndexImag(j, k)) -= T0(0, A.ValueImag(j, k)) * X(j);
	    }
	}
  }
  
  
  //! Successive overrelaxation.
  /*!
    Solving A^T X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);

    complex<T1> temp, zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    long* ptr_real = A.GetRealPtr();
    long* ptr_imag = A.GetImagPtr();
    int* ind_real = A.GetRealInd();
    int* ind_imag = A.GetImagInd();
    
    typename ClassComplexType<T0>::Treal* data_real = A.GetRealData();
    typename ClassComplexType<T0>::Treal* data_imag = A.GetImagData();
    
    T0 ajj;

    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = 0; j < ma; j++)
	  {
	    temp = zero;
            ajj = zero;
            for (long k = ptr_real[j]; k < ptr_real[j+1]; k++)
	      {
                if (ind_real[k] != j)
                  temp += data_real[k] * X(ind_real[k]);
                else
                  ajj = T0(data_real[k], 0);
              }

	    for (long k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
              {
                if (ind_imag[k] != j)
                  temp += complex<T1>(0, data_imag[k]) * X(ind_imag[k]);
                else
                  ajj += complex<T1>(0, data_imag[k]);
              }
            
#ifdef SELDON_CHECK_BOUNDS
            if (ajj == zero)
              throw WrongArgument("SOR", "Matrix must contain"
                                  " a non-null diagonal");
#endif
            
	    X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
	  }

    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = ma-1; j >= 0; j--)
	  {
	    temp = zero;
            ajj = zero;
            for (long k = ptr_real[j]; k < ptr_real[j+1]; k++)
	      {
                if (ind_real[k] != j)
                  temp += data_real[k] * X(ind_real[k]);
                else
                  ajj = T0(data_real[k], 0);
              }

	    for (long k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
              {
                if (ind_imag[k] != j)
                  temp += complex<T1>(0, data_imag[k]) * X(ind_imag[k]);
                else
                  ajj += complex<T1>(0, data_imag[k]);
              }
            
	    X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
	  }
  }


  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);

    complex<T1> temp, zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    T0 ajj;

    if (type_ssor%2 == 0)
      for (int i = 0; i < iter; i++)
        for (int j = 0; j < ma; j++)
          {
            temp = zero;
            ajj = zero;
            for (int k = 0; k < A.GetRealColumnSize(j); k++)
              {
                if (A.IndexReal(j, k) != j)
                  temp += A.ValueReal(j,k) * X(A.IndexReal(j,k));
                else
                  ajj += A.ValueReal(j, k);
              }
            
            for (int k = 0; k < A.GetImagColumnSize(j); k++)
              {
                if (A.IndexImag(j, k) != j)
                  temp += T0(0, A.ValueImag(j,k))
                    * X(A.IndexImag(j, k));
                else
                  ajj += T0(0, A.ValueImag(j,k));
              }
            
#ifdef SELDON_CHECK_BOUNDS
            if (ajj == zero)
              throw WrongArgument("SOR", "Matrix must contain"
                                  " a non-null diagonal");
#endif
            
            X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
          }

    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = ma-1; j >= 0; j--)
	  {
	      temp = zero;
	      ajj = zero;
	      for (int k = 0; k < A.GetRealColumnSize(j); k++)
		{
		  if (A.IndexReal(j, k) != j)
                    temp += A.ValueReal(j,k) * X(A.IndexReal(j,k));
		  else
		    ajj += A.ValueReal(j, k);
		}
              
	      for (int k = 0; k < A.GetImagColumnSize(j); k++)
		{
		  if (A.IndexImag(j, k) != j)
                    temp += T0(0, A.ValueImag(j,k))
                      * X(A.IndexImag(j, k));
		  else
		    ajj += T0(0, A.ValueImag(j,k));
		}
              
              X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
	  }
  }

  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);

    complex<T1> zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    typename ClassComplexType<T0>::Treal* data_real = A.GetRealData();

    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    typename ClassComplexType<T0>::Treal* data_imag = A.GetImagData();

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A

    // Forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T0 ajj;
    T3 coef = (one - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (U + (1-omega)/omega D) X + B	      
	  for (int j = 0; j < ma; j++)
	    {
              long kr = ptr_real[j];
              while ( (kr < ptr_real[j+1]) && (ind_real[kr] < j))
                {
                  X(ind_real[kr]) -= data_real[kr]*X(j);
                  kr++;
                }
              
              if ( (kr < ptr_real[j+1]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr], 0);

              long ki = ptr_imag[j];
              while ( (ki < ptr_imag[j+1]) && (ind_imag[ki] < j))
                {
                  X(ind_imag[ki]) -= T0(0, data_imag[ki])*X(j);
                  ki++;
                }

              if ( (ki < ptr_imag[j+1]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki]);
              
#ifdef SELDON_CHECK_BOUNDS
              if (ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
              X(j) = B(j) + coef * ajj * X(j);
	    }

          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              long kr = ptr_real[j];
              while ( (kr < ptr_real[j+1]) && (ind_real[kr] < j))
                kr++;
                             
              if ( (kr < ptr_real[j+1]) && (ind_real[kr] == j))
                {
                  ajj = T0(data_real[kr], 0);
                  kr++;
                }

              long ki = ptr_imag[j];
              while ( (ki < ptr_imag[j+1]) && (ind_imag[ki] < j))
                ki++;
              
              if ( (ki < ptr_imag[j+1]) && (ind_imag[ki] == j))
                {
                  ajj += T0(0, data_imag[ki]);
                  ki++;
                }
              
              X(j) *= omega/ajj;
              while (kr < ptr_real[j+1])
                {
                  X(ind_real[kr]) -= data_real[kr]*X(j);
                  kr++;
                }

              while (ki < ptr_imag[j+1])
                {
                  X(ind_imag[ki]) -= T0(0, data_imag[ki])*X(j);
                  ki++;
                }
            }
	}
    
    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (L + (1-omega)/omega D) X + B.
          for (int j = ma-1; j >= 0; j--)
            {
              long kr = ptr_real[j+1]-1;
              while ( (kr >= ptr_real[j]) && (ind_real[kr] > j) )
                {
                  X(ind_real[kr]) -= data_real[kr]*X(j);
                  kr--;
                }
              
              if ( (kr >= ptr_real[j]) && (ind_real[kr] == j))
                ajj = T0(data_real[kr], 0);
              
              long ki = ptr_imag[j+1]-1;
              while ( (ki >= ptr_imag[j]) && (ind_imag[ki] > j) )
                {
                  X(ind_imag[ki]) -= T0(0, data_imag[ki])*X(j);
                  ki--;
                }
              
              if ( (ki >= ptr_imag[j]) && (ind_imag[ki] == j))
                ajj += T0(0, data_imag[ki]);
              
#ifdef SELDON_CHECK_BOUNDS
              if (ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
              X(j) = B(j) + coef * ajj * X(j);
            }
                    
          // sThen we solve (D/omega - U) X = X.
          for (int j = ma-1; j >= 0; j--)
            {
              long kr = ptr_real[j+1]-1;
              while ( (kr >= ptr_real[j]) && (ind_real[kr] > j) )
                kr--;
              
              if ( (kr >= ptr_real[j]) && (ind_real[kr] == j))
                {
                  ajj = T0(data_real[kr], 0);
                  kr--;
                }
              
              long ki = ptr_imag[j+1]-1;
              while ( (ki >= ptr_imag[j]) && (ind_imag[ki] > j) )
                ki--;
              
              if ( (ki >= ptr_imag[j]) && (ind_imag[ki] == j))
                {
                  ajj += T0(0, data_imag[ki]);
                  ki--;
                }

              X(j) *= omega/ajj;
              while (kr >= ptr_real[j])
                {
                  X(ind_real[kr]) -= data_real[kr]*X(j);
                  kr--;
                }

              while (ki >= ptr_imag[j])
                {
                  X(ind_imag[ki]) -= T0(0, data_imag[ki])*X(j);
                  ki--;
                }
            }
        }
  }

  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega,int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);

    complex<T1> zero; T3 one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetM() || ma != B.GetM())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A

    // Forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T0 ajj;
    T3 coef = (one - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (U + (1-omega)/omega D) X + B	      
	  for (int j = 0; j < ma; j++)
	    {
              int kr = 0;
              while ( (kr < A.GetRealRowSize(j)) && (A.IndexReal(j, kr) < j))
                {
                  X(A.IndexReal(j, kr)) -= A.ValueReal(j, kr) * X(j);
                  kr++;
                }
              
              if ( (kr < A.GetRealRowSize(j)) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr), 0);
              
              int ki = 0;
              while ( (ki < A.GetImagRowSize(j)) && (A.IndexImag(j, ki) < j))
                {
                  X(A.IndexImag(j, ki)) -= T0(0, A.ValueImag(j, ki))*X(j);
                  ki++;
                }

              if ( (ki < A.GetImagRowSize(j)) && (A.IndexImag(j, ki) == j) )
                ajj += T0(0, A.ValueImag(j, ki));
              
#ifdef SELDON_CHECK_BOUNDS
              if (ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
              X(j) = B(j) + coef * ajj * X(j);
	    }

          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              int kr = 0;
              while ( (kr < A.GetRealRowSize(j)) && (A.IndexReal(j, kr) < j))
                kr++;
                             
              if ( (kr < A.GetRealRowSize(j)) && (A.IndexReal(j, kr) == j) )
                {
                  ajj = T0(A.ValueReal(j, kr), 0);
                  kr++;
                }

              int ki = 0;
              while ( (ki < A.GetImagRowSize(j)) && (A.IndexImag(j, ki) < j))
                ki++;
              
              if ( (ki < A.GetImagRowSize(j)) && (A.IndexImag(j, ki) == j))
                {
                  ajj += T0(0, A.ValueImag(j, ki));
                  ki++;
                }
              
              X(j) *= omega/ajj;
              while (kr < A.GetRealRowSize(j))
                {
                  X(A.IndexReal(j, kr)) -= A.ValueReal(j, kr)*X(j);
                  kr++;
                }

              while (ki < A.GetImagRowSize(j))
                {
                  X(A.IndexImag(j, ki)) -= T0(0, A.ValueImag(j, ki))*X(j);
                  ki++;
                }
            }
	}
    
    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we compute X = (L + (1-omega)/omega D) X + B.
          for (int j = ma-1; j >= 0; j--)
            {
              int kr = A.GetRealRowSize(j)-1;
              while ( (kr >= 0) && (A.IndexReal(j, kr) > j) )
                {
                  X(A.IndexReal(j, kr)) -= A.ValueReal(j, kr)*X(j);
                  kr--;
                }
              
              if ( (kr >= 0) && (A.IndexReal(j, kr) == j))
                ajj = T0(A.ValueReal(j, kr), 0);
              
              int ki = A.GetImagRowSize(j)-1;
              while ( (ki >= 0) && (A.IndexImag(j, ki) > j) )
                {
                  X(A.IndexImag(j, ki)) -= T0(0, A.ValueImag(j, ki))*X(j);
                  ki--;
                }
              
              if ( (ki >= 0) && (A.IndexImag(j, ki) == j))
                ajj += T0(0, A.ValueImag(j, ki));
              
#ifdef SELDON_CHECK_BOUNDS
              if (ajj == zero)
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              
              X(j) = B(j) + coef * ajj * X(j);
            }
                    
          // Then we solve (D/omega - U) X = X.
          for (int j = ma-1; j >= 0; j--)
            {
              int kr = A.GetRealRowSize(j)-1;
              while ( (kr >= 0) && (A.IndexReal(j, kr) > j) )
                kr--;
              
              if ( (kr >= 0) && (A.IndexReal(j, kr) == j))
                {
                  ajj = T0(A.ValueReal(j, kr), 0);
                  kr--;
                }
              
              int ki = A.GetImagRowSize(j)-1;
              while ( (ki >= 0) && (A.IndexImag(j, ki) > j) )
                ki--;
              
              if ( (ki >= 0) && (A.IndexImag(j, ki) == j))
                {
                  ajj += T0(0, A.ValueImag(j, ki));
                  ki--;
                }

              X(j) *= omega/ajj;
              while (kr >= 0)
                {
                  X(A.IndexReal(j, kr)) -= A.ValueReal(j, kr)*X(j);
                  kr--;
                }

              while (ki >= 0)
                {
                  X(A.IndexImag(j, ki)) -= T0(0, A.ValueImag(j, ki))*X(j);
                  ki--;
                }
            }
        }
  }

  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    SorVector(A, X, B, omega, iter, type_ssor);
  }
  
  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    SorVector(A, X, B, omega, iter, type_ssor);
  }

  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    SorVector(A, X, B, omega, iter, type_ssor);
  }


  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
		 Vector<complex<T2>, Storage2, Allocator2>& X,
		 const Vector<complex<T1>, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    SorVector(A, X, B, omega, iter, type_ssor);
  }

  
  // SOR //
  /////////
  
}

#define SELDON_FILE_FUNCTIONS_MATVECT_COMPLEX_CXX
#endif

