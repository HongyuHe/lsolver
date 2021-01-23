// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2001-2011 Marc Duruflé
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


#ifndef SELDON_FILE_FUNCTIONS_MATVECT_CXX
#define SELDON_FILE_FUNCTIONS_MATVECT_CXX


#include "Functions_MatVect.hxx"


/*
  Functions defined in this file:

  alpha M X + beta Y -> Y
  MltAdd(alpha, M, X, beta, Y)

  Gauss(M, X)

  GaussSeidel(M, Y, X, iter)

  SOR(M, Y, X, omega, iter)

  SolveLU(M, Y)

  GetAndSolveLU(M, Y)
*/


namespace Seldon
{

  /////////
  // MLT //

  
  // Y = M X for RowSparse matrices
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "Mlt(M, X, Y)");
#endif

    T4 zero, temp;
    SetComplexZero(zero);

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();
    
    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (long j = ptr[i]; j < ptr[i+1]; j++)
	  temp += data[j] * X(ind[j]);
	
	Y(i) = temp;
      }
  }
  

  // Y = M X for ColSparse matrices
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, ColSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "Mlt(M, X, Y)");
#endif

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();
    Y.Zero();

    for (int j = 0; j < M.GetN(); j++)
      {
	for (long k = ptr[j]; k < ptr[j+1]; k++)
	  Y(ind[k]) += data[k] * X(j);
      }
  }

  
  // Y = M X for RowSymSparse
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "Mlt(M, X, Y)");
#endif

    int i; long j;
    T4 zero, temp;
    SetComplexZero(zero);

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += data[j] * X(ind[j]);
	
	Y(i) = temp;
      }
    
    for (i = 0; i < ma-1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += data[j] * X(i);
  }
  

  // Y = M X for ColSymSparse
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, ColSymSparse, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "Mlt(M, X, Y)");
#endif

    int i; long j;
    T4 zero, temp;
    SetComplexZero(zero);

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  if (ind[j] != i)
            temp += data[j] * X(ind[j]);
        
	Y(i) = temp;
      }

    for (i = 0; i < ma; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	Y(ind[j]) += data[j] * X(i);
    
  }

  
  // Y = M X or M^T X or M^H X for RowSparse
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
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
    CheckDim(Trans, M, X, Y, "Mlt(SeldonTrans, M, X, Y)");
#endif

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();
    Y.Zero();
    
    if (Trans.Trans())
      {
	for (i = 0; i < ma; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    Y(ind[j]) += data[j] * X(i);
      }
    else
      {
	for (i = 0; i < ma; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    Y(ind[j]) += conjugate(data[j]) * X(i);
      }
  }
  
  
  // Y = M X or M^T X or M^H X for ColSparse
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, ColSparse, Allocator1>& M,
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
    CheckDim(Trans, M, X, Y, "Mlt(SeldonTrans, M, X, Y)");
#endif
    
    T4 temp, zero;
    SetComplexZero(zero);
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();
    
    if (Trans.Trans())
      {
	for (i = 0; i < M.GetN(); i++)
	  {
	    temp = zero;
	    for (j = ptr[i]; j < ptr[i + 1]; j++)
	      temp += data[j] * X(ind[j]);
	    
	    Y(i) = temp;
	  }
      }
    else
      {
	for (i = 0; i < M.GetN(); i++)
	  {
	    temp = zero;
	    for (j = ptr[i]; j < ptr[i + 1]; j++)
	      temp += conjugate(data[j]) * X(ind[j]);
	    
	    Y(i) = temp;
	  }
      }
  }


  // Y = M X or M^T X or M^H X for RowSymSparse
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
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
    CheckDim(Trans, M, X, Y, "Mlt(SeldonConjTrans, M, X, Y)");
#endif

    int i; long j;
    T4 zero, temp;
    SetComplexZero(zero);
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += conjugate(data[j]) * X(ind[j]);
	
        Y(i) = temp;
      }
    
    for (i = 0; i < ma - 1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += conjugate(data[j]) * X(i);
  }
  
    
  // Y = M X or M^T X or M^H X for ColSymSparse
  template <class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, ColSymSparse, Allocator1>& M,
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
    CheckDim(Trans, M, X, Y, "Mlt(SeldonConjTrans, M, X, Y)");
#endif

    int i; long j;
    T4 zero, temp;
    SetComplexZero(zero);

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += conjugate(data[j]) * X(ind[j]);
	
        Y(i) = temp;
      }
    
    for (i = 0; i < ma; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += conjugate(data[j]) * X(i);
  }

  
  // Y = M X for any matrix
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, Storage1, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();
    int na = M.GetN();
    
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "Mlt(M, X, Y)");
#endif

    // aborting computation if the matrix is sparse
    if (Storage1::Sparse)
      throw WrongArgument("Mlt", "This function is intended to dense"
                          " matrices only and not to sparse matrices");
    
    T4 zero, temp;
    SetComplexZero(zero);

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (int j = 0; j < na; j++)
	  temp += M(i, j) * X(j);
	
	Y(i) = temp;
      }
  }


  // Y = M X or M^T X or M^H X for any matrix
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T4, class Storage4, class Allocator4>
  void MltVector(const SeldonTranspose& Trans,
		 const Matrix<T1, Prop1, Storage1, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
	MltVector(M, X, Y);
	return;
      }

    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "Mlt(trans, M, X, Y)");
#endif
    
    if (Storage1::Sparse)
      throw WrongArgument("MltAdd", "This function is intended to dense"
                          " matrices only and not to sparse matrices");

    T4 zero, temp;
    SetComplexZero(zero);
    
    if (Trans.Trans())
      {
	for (int i = 0; i < na; i++)
	  {
	    temp = zero;
	    for (int j = 0; j < ma; j++)
	      temp += M(j, i) * X(j);
	    
	    Y(i) = temp;
	  }
      }
    else
      {
	for (int i = 0; i < na; i++)
	  {
	    temp = zero;
	    for (int j = 0; j < ma; j++)
	      temp += conjugate(M(j, i)) * X(j);
	    
	Y(i) = temp;
	  }
      }
  }
  
    
  // MLT //
  /////////


  ////////////
  // MLTADD //


  /*** PETSc matrices ***/


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
		    const Vector<T2, PETScSeq, Allocator2>& X,
		    const T3& beta, Vector<T4, PETScSeq, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif
    if (beta == T3(0))
      {
	if (alpha == T0(0))
	  {
	    Y.Zero();
	    return;
	  }
	else
	  {
	    MatMult(M.GetPetscMatrix(), X.GetPetscVector(), Y.GetPetscVector());
	    if (alpha != T0(1))
	      VecScale(Y.GetPetscVector(), alpha);
	    return;
	  }
      }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScSeq, Allocator2> tmp;
    tmp.Copy(Y);
    MatMult(M.GetPetscMatrix(), X.GetPetscVector(), tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
		    const Vector<T2, PETScPar, Allocator2>& X,
		    const T3& beta, Vector<T4, PETScPar, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif
    if (beta == T3(0))
      {
	if (alpha == T0(0))
	  {
	    Y.Zero();
	    return;
	  }
	else
	  {
	    MatMult(M.GetPetscMatrix(), X.GetPetscVector(), Y.GetPetscVector());
	    if (alpha != T0(1))
	      VecScale(Y.GetPetscVector(), alpha);
	    return;
	  }
      }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScPar, Allocator2> tmp;
    tmp.Copy(Y);
    MatMult(M.GetPetscMatrix(), X.GetPetscVector(), tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, PETScSeq, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Vector<T4, PETScSeq, Allocator4> X_Petsc;
    X_Petsc.Reallocate(X.GetM());
    for (int i = 0; i < X.GetM(); i++)
      X_Petsc.SetBuffer(i, X(i));
    X_Petsc.Flush();

    if (beta == T3(0))
      {
	if (alpha == T0(0))
	  {
	    Y.Zero();
	    return;
	  }
	else
	  {
	    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
		    Y.GetPetscVector());
	    if (alpha != T0(1))
	      VecScale(Y.GetPetscVector(), alpha);
	    return;
	  }
      }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScSeq, Allocator2> tmp;
    tmp.Copy(Y);
    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
            tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, PETScPar, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Vector<T4, PETScPar, Allocator4> X_Petsc;
    X_Petsc.SetCommunicator(M.GetCommunicator());
    X_Petsc.Reallocate(X.GetM());
    for (int i = 0; i < X.GetM(); i++)
      X_Petsc.SetBuffer(i, X(i));
    X_Petsc.Flush();

    if (beta == T3(0))
      {
	if (alpha == T0(0))
	  {
	    Y.Zero();
	    return;
	  }
	else
	  {
	    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
		    Y.GetPetscVector());
	    if (alpha != T0(1))
	      VecScale(Y.GetPetscVector(), alpha);
	    return;
	  }
      }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScPar, Allocator2> tmp;
    tmp.Copy(Y);
    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
            tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, PETScPar, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Vector<T4, PETScPar, Allocator4> X_Petsc;
     X_Petsc.SetCommunicator(M.GetCommunicator());
    X_Petsc.Reallocate(X.GetM());
    for (int i = 0; i < X.GetM(); i++)
      X_Petsc.SetBuffer(i, X(i));
    X_Petsc.Flush();

    if (beta == T3(0))
      {
	if (alpha == T0(0))
	  {
	    Y.Zero();
	    return;
	  }
	else
	  {
	    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
		    Y.GetPetscVector());
	    if (alpha != T0(1))
	      VecScale(Y.GetPetscVector(), alpha);
	    return;
	  }
      }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScPar, Allocator2> tmp;
    tmp.Copy(Y);
    tmp.SetCommunicator(M.GetCommunicator());
    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
            tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }



  /*** Sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    T4 zero, temp;
    SetComplexZero(zero);

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSparse, Allocator1>::pointer
      data = M.GetData();

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (long j = ptr[i]; j < ptr[i+1]; j++)
	  temp += data[j] * X(ind[j]);
	Y(i) += alpha * temp;
      }
  }


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Collection, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    typename T4::value_type zero;
    SetComplexZero(zero);
    typename T4::value_type temp;

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSparse, Allocator1>::pointer
      data = M.GetData();

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (long j = ptr[i]; j < ptr[i+1]; j++)
	  temp += data[j] * X(ind[j]);
	Y(i) += alpha * temp;
      }
  }

  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, ColSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, ColSparse, Allocator1>::pointer
      data = M.GetData();

    for (int j = 0; j < M.GetN(); j++)
      {
	for (long k = ptr[j]; k < ptr[j+1]; k++)
	  Y(ind[k]) += alpha * data[k] * X(j);
      }
  }

  
  /*** Symmetric sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i; long j;
    T4 zero;
    SetComplexZero(zero);
    T4 temp;

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSymSparse, Allocator1>::pointer
      data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += data[j] * X(ind[j]);
	Y(i) += alpha * temp;
      }
    for (i = 0; i < ma-1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += alpha * data[j] * X(i);
  }

  
  template <class T0,
            class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, ColSymSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i; long j;
    T4 zero;
    SetComplexZero(zero);
    T4 temp;

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, ColSymSparse, Allocator1>::pointer
      data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  if (ind[j] != i)
            temp += data[j] * X(ind[j]);
        
	Y(i) += alpha * temp;

        for (j = ptr[i]; j < ptr[i + 1]; j++)
	  Y(ind[j]) += alpha * data[j] * X(i);
      }
  }

  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
		    const Vector<T2, Collection, Allocator2>& X,
		    const T3& beta, Vector<T4, Collection, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i; long j;
    typename T4::value_type zero;
    SetComplexZero(zero);
    typename T4::value_type temp;

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSymSparse, Allocator1>::pointer
      data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += data[j] * X(ind[j]);
	Y(i) += alpha * temp;
      }
    for (i = 0; i < ma-1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += alpha * data[j] * X(i);
  }


  /*** Symmetric complex sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Collection, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i; long j;
    typename T4::value_type zero;
    SetComplexZero(zero);
    typename T4::value_type temp;

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSymSparse, Allocator1>::pointer
      data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += data[j] * X(ind[j]);
	Y(i) += alpha * temp;
      }
    for (i = 0; i < ma-1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += alpha * data[j] * X(i);
  }


  /*** Sparse matrices, *Trans ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
	MltAdd(alpha, M, X, beta, Y);
	return;
      }
    
    int i; long j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSparse, Allocator1>::pointer
      data = M.GetData();
    
    if (Trans.Trans())
      {
	for (i = 0; i < ma; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    Y(ind[j]) += alpha * data[j] * X(i);
      }
    else
      {
	for (i = 0; i < ma; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    Y(ind[j]) += alpha * conjugate(data[j]) * X(i);
      }
  }

  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, ColSparse, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
	MltAdd(alpha, M, X, beta, Y);
	return;
      }

    int i; long j;

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);
    
    T4 temp, zero;
    SetComplexZero(zero);
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, ColSparse, Allocator1>::pointer
      data = M.GetData();

    if (Trans.Trans())
      {
	for (i = 0; i < M.GetN(); i++)
	  {
	    temp = zero;
	    for (j = ptr[i]; j < ptr[i + 1]; j++)
	      temp += data[j] * X(ind[j]);
	    
	    Y(i) += alpha * temp;
	  }
      }
    else
      {
	for (i = 0; i < M.GetN(); i++)
	  {
	    temp = zero;
	    for (j = ptr[i]; j < ptr[i + 1]; j++)
	      temp += conjugate(data[j]) * X(ind[j]);
	    
	    Y(i) += alpha * temp;
	  }
      }
  }

  
  /*** Symmetric sparse matrices, *Trans ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
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
    T4 zero, temp;
    SetComplexZero(zero);
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += conjugate(data[j]) * X(ind[j]);
	
        Y(i) += alpha * temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += alpha * conjugate(data[j]) * X(i);
  }
  
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, ColSymSparse, Allocator1>& M,
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
    T4 zero, temp;
    SetComplexZero(zero);
    
    long* ptr = M.GetPtr();
    int* ind = M.GetInd();
    T1* data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += conjugate(data[j]) * X(ind[j]);
	
        Y(i) += alpha * temp;
        
        for (j = ptr[i]; j < ptr[i + 1]; j++)
          if (ind[j] != i)
            Y(ind[j]) += alpha * conjugate(data[j]) * X(i);
      }
  }

  
  // MLTADD //
  ////////////


  ////////////
  // MLTADD //


  /*! \brief Performs the multiplication of a matrix with a vector, and adds
    the result to another vector.
   */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, \f$ M \f$ is a \f$ m \times n
    \f$ matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$
    Y \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] M m by n matrix.
    \param[in] X vector of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector of length m, result of the product of \a M by \a
    X, times \a alpha, plus \a Y (on entry) times \a beta.
  */
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, Storage1, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    // aborting computation if the matrix is sparse
    if (Storage1::Sparse)
      throw WrongArgument("MltAdd", "This function is intended to dense"
                          " matrices only and not to sparse matrices");
    
    Mlt(beta, Y);
    
    T4 zero;
    SetComplexZero(zero);
    T4 temp;

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (int j = 0; j < na; j++)
	  temp += M(i, j) * X(j);
	Y(i) += alpha * temp;
      }
  }


  /*! \brief Performs the multiplication of a matrix with a vector, and adds
    the result to another vector.
   */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, \f$ M \f$ is a \f$ m \times n
    \f$ matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$
    Y \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] M m by n matrix.
    \param[in] X vector of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector of length m, result of the product of \a M by \a
    X, times \a alpha, plus \a Y (on entry) times \a beta.
  */
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, Storage1, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Collection, Allocator4>& Y)
  {
    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    if (Storage1::Sparse)
      throw WrongArgument("MltAdd", "This function is intended to dense"
                          " matrices only and not to sparse matrices");
    
    Mlt(beta, Y);

    typename T4::value_type zero;
    SetComplexZero(zero);
    typename T4::value_type temp;

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (int j = 0; j < na; j++)
	  temp += M(i, j) * X(j);
	Y(i) += alpha * temp;
      }
  }


  /*! \brief Performs the multiplication of a matrix collection with a vector
    collection, and adds the result to another vector.
  */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, \f$ M \f$ is a \f$ m \times n
    \f$ matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$
    Y \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] M m by n matrix colection.
    \param[in] X vector collection of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector collection of length m, result of the product of
    \a M by \a X, times \a alpha, plus \a Y (on entry) times \a beta.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, RowMajorCollection, Allocator1>& M,
		    const Vector<T2, Collection, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Collection, Allocator4>& Y)
  {
    int ma = M.GetMmatrix();
    int na = M.GetNmatrix();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif
    typedef typename T4::value_type value_type;

    Mlt(value_type(beta), Y);

    for (int i = 0; i < ma; i++)
      for (int j = 0; j < na; j++)
      	MltAdd(alpha, M.GetMatrix(i, j), X.GetVector(j), value_type(1.),
               Y.GetVector(i));
  }


   /*! \brief Performs the multiplication of a matrix collection with a vector
    collection, and adds the result to another vector.
  */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, \f$ M \f$ is a \f$ m \times n
    \f$ matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$
    Y \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] M m by n matrix colection.
    \param[in] X vector collection of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector collection of length m, result of the product of
    \a M by \a X, times \a alpha, plus \a Y (on entry) times \a beta.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Prop1, ColMajorCollection, Allocator1>& M,
		    const Vector<T2, Collection, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Collection, Allocator4>& Y)
  {
    int ma = M.GetMmatrix();
    int na = M.GetNmatrix();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif
    typedef typename T4::value_type value_type;

    Mlt(value_type(beta), Y);

    for (int i = 0; i < ma; i++)
      for (int j = 0; j < na; j++)
      	MltAdd(alpha, M.GetMatrix(i, j), X.GetVector(j), value_type(1.),
               Y.GetVector(i));
  }


  /*! \brief Performs the multiplication of a matrix (possibly transposed)
    with a vector, and adds the result to another vector.
   */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ or \f$ Y =
    \alpha M^T X + \beta Y \f$ where \f$ \alpha \f$ and \f$ \beta \f$ are
    scalars, \f$ M \f$ is a \f$ m \times n \f$ matrix or a \f$ n \times m \f$
    matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$ Y
    \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] Trans transposition status of \a M: it may be SeldonNoTrans,
    SeldonTrans or SeldonConjTrans
    \param[in] M m by n matrix, or n by m matrix if transposed.
    \param[in] X vector of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector of length m, result of the product of \a M by \a
    X, times \a alpha, plus \a Y (on entry) times \a beta.
  */
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const Matrix<T1, Prop1, Storage1, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
	MltAddVector(alpha, M, X, beta, Y);
	return;
      }

    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, trans, M, X, beta, Y)");
#endif
    
    if (Storage1::Sparse)
      throw WrongArgument("MltAdd", "This function is intended to dense"
                          " matrices only and not to sparse matrices");

    T3 zero3;
    SetComplexZero(zero3);
    T4 zero;
    SetComplexZero(zero);
    
    if (beta == zero3)
      Y.Zero();
    else
      Mlt(beta, Y);

    T4 temp;

    if (Trans.Trans())
      {
	for (int i = 0; i < na; i++)
	  {
	    temp = zero;
	    for (int j = 0; j < ma; j++)
	      temp += M(j, i) * X(j);
	    Y(i) += alpha * temp;
	  }
      }
    else
      {
	for (int i = 0; i < na; i++)
	  {
	    temp = zero;
	    for (int j = 0; j < ma; j++)
	      temp += conjugate(M(j, i)) * X(j);
	    Y(i) += alpha * temp;
	  }
      }
  }

  
  // MLTADD //
  ////////////


  ///////////
  // GAUSS //


  //! Solves M*Y = X with Gauss method.
  /*!
    Warning: M is modified. The results are stored in X.
    There is no partial pivoting performed here,
    the method will fail if a diagonal coefficient 
    generated during the factorisation is equal to 0
    For dense matrices, use rather GetLU/SolveLU
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Gauss(Matrix<T0, Prop0, Storage0, Allocator0>& M,
             Vector<T1, Storage1, Allocator1>& X)
  {
    int i, j, k;
    T1 r, S;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    int na = M.GetN();
    if (na != ma)
      throw WrongDim("Gauss(M, X)",
		     "The matrix must be squared.");

    CheckDim(M, X, "Gauss(M, X)");
#endif

    // aborting computation if the matrix is sparse
    if (Storage0::Sparse)
      throw WrongArgument("Gauss", "This function is intended to dense"
                          " matrices only and not to sparse matrices");

    for (k = 0; k < ma - 1; k++)
      for (i = k + 1; i < ma; i++)
	{
	  r = M(i, k) / M(k, k);
	  for (j = k + 1; j < ma; j++)
	    M(i, j) -= r * M(k, j);
	  X(i) -= r *= X(k);
	}

    X(ma - 1) = X(ma - 1) / M(ma - 1, ma - 1);
    for (k = ma - 2; k > -1; k--)
      {
	S = X(k);
	for (j = k + 1; j < ma; j++)
	  S -= M(k, j) * X(j);
	X(k) = S / M(k, k);
      }
  }


  // GAUSS //
  ///////////


  //////////////////
  // GAUSS-SEIDEL //


  //! Solve M*Y = X with Gauss-Seidel method.
  /*!
    Solving M Y = X by using Gauss-Seidel algorithm.
    iter is the number of iterations.
    type_algo = 2 forward sweep
    type_algo = 3 backward sweep
    type_algo = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void GaussSeidel(const Matrix<T0, Prop0, Storage0, Allocator0>& M,			  
                   Vector<T2, Storage2, Allocator2>& Y,
                   const Vector<T1, Storage1, Allocator1>& X,
                   int iter, int type_algo)
  {
    int i, j, k;
    T1 temp;

    int na = M.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    int ma = M.GetM();
    if (na != ma)
      throw WrongDim("GaussSeidel(M, X, Y, iter)",
		     "The matrix must be squared.");

    CheckDim(M, X, Y, "GaussSeidel(M, X, Y, iter)");
#endif

    // aborting computation if the matrix is sparse
    if (Storage0::Sparse)
      throw WrongArgument("GaussSeidel", "This function is intended to dense"
                          " matrices only and not to sparse matrices");
    
    T1 zero;
    SetComplexZero(zero);
    if (type_algo%2 == 0)
      for (i = 0; i < iter; i++)
        for (j = 0; j < na; j++)
          {
            temp = zero;
            for (k = 0; k < j; k++)
              temp -= M(j, k) * Y(k);
            for (k = j + 1; k < na; k++)
              temp -= M(j, k) * Y(k);
            Y(j) = (X(j) + temp) / M(j, j);
          }
    
    if (type_algo%3 == 0)
      for (i = 0; i < iter; i++)
        for (j = na-1; j >= 0; j--)
          {
            temp = zero;
            for (k = 0; k < j; k++)
              temp -= M(j, k) * Y(k);
            for (k = j + 1; k < na; k++)
              temp -= M(j, k) * Y(k);
            Y(j) = (X(j) + temp) / M(j, j);
          }    
  }


  // GAUSS-SEIDEL //
  //////////////////


  ///////////////////
  // S.O.R. METHOD //


  //! Solve M Y = X with S.O.R. method.
  /*!
    Solving M Y = X by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3>
  void SorVector(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		 Vector<T2, Storage2, Allocator2>& Y,
		 const Vector<T1, Storage1, Allocator1>& X,
		 const T3& omega, int iter, int type_ssor)
  {
    int i, j, k;
    T1 temp;
    T3 one;

    int na = M.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    int ma = M.GetM();
    if (na != ma)
      throw WrongDim("SOR(M, X, Y, omega, iter)",
		     "The matrix must be squared.");

    CheckDim(M, X, Y, "SOR(M, X, Y, omega, iter)");
#endif

    // aborting computation if the matrix is sparse
    if (Storage0::Sparse)
      throw WrongArgument("SOR", "This function is intended to dense"
                          " matrices only and not to sparse matrices");
    
    SetComplexOne(one);
    if (type_ssor%2 == 0)
      for (i = 0; i < iter; i++)
        for (j = 0; j < na; j++)
          {
            SetComplexZero(temp);
            for (k = 0; k < j; k++)
              temp -= M(j, k) * Y(k);
            for (k = j + 1; k < na; k++)
              temp -= M(j, k) * Y(k);
            Y(j) = (one - omega) * Y(j) + omega * (X(j) + temp) / M(j, j);
          }

    if (type_ssor%3 == 0)
      for (i = 0; i < iter; i++)
        for (j = na-1; j >= 0; j--)
          {
            SetComplexZero(temp);
            for (k = 0; k < j; k++)
              temp -= M(j, k) * Y(k);
            for (k = j + 1; k < na; k++)
              temp -= M(j, k) * Y(k);
            Y(j) = (one - omega) * Y(j) + omega * (X(j) + temp) / M(j, j);
          }
  }

  
  //! Solve M^T Y = X with S.O.R. method.
  /*!
    Solving M Y = X by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		 Vector<T2, Storage2, Allocator2>& Y,
		 const Vector<T1, Storage1, Allocator1>& X,
		 const T3& omega, int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(M, Y, X, omega, iter, type_ssor);
    
    int i, j, k;
    T1 temp;
    T3 one;

    int na = M.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    int ma = M.GetM();
    if (na != ma)
      throw WrongDim("SOR(M, X, Y, omega, iter)",
		     "The matrix must be squared.");

    CheckDim(M, X, Y, "SOR(M, X, Y, omega, iter)");
#endif

    // aborting computation if the matrix is sparse
    if (Storage0::Sparse)
      throw WrongArgument("SOR", "This function is intended to dense"
                          " matrices only and not to sparse matrices");
    
    SetComplexOne(one);
    if (type_ssor%2 == 0)
      for (i = 0; i < iter; i++)
        for (j = 0; j < na; j++)
          {
            SetComplexZero(temp);
            for (k = 0; k < j; k++)
              temp -= M(k, j) * Y(k);
            for (k = j + 1; k < na; k++)
              temp -= M(k, j) * Y(k);
            Y(j) = (one - omega) * Y(j) + omega * (X(j) + temp) / M(j, j);
          }

    if (type_ssor%3 == 0)
      for (i = 0; i < iter; i++)
        for (j = na-1; j >= 0; j--)
          {
            SetComplexZero(temp);
            for (k = 0; k < j; k++)
              temp -= M(k, j) * Y(k);
            for (k = j + 1; k < na; k++)
              temp -= M(k, j) * Y(k);
            Y(j) = (one - omega) * Y(j) + omega * (X(j) + temp) / M(j, j);
          }
  }

  
  // S.O.R. method //
  ///////////////////



  /////////////
  // SOLVELU //


  //! Solves a linear system whose matrix has been LU-factorized.
  /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
    \f$ X \f$ and \f$ Y \f$ are vectors. The matrix \a M cannot be provided as
    such to this function: it must already be factorized in LU form.
    \param[in] M the matrix of the linear system, already factorized in LU
    form. The lower part of \a M should be \a L, and the upper part should be
    \a U. The diagonal of \a M should be the diagonal of \a U. The diagonal
    elements of \a L are assumed to be ones.
    \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
    solution \f$ X \f$ of the system.
    \sa Seldon::GetLU(Matrix<T0, Prop0, Storage0, Allocator0>& A) to factorize
    a matrix before using this function.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SolveLuVector(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		     Vector<T1, Storage1, Allocator1>& Y)
  {
    int i, k;
    T1 temp;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    int na = M.GetN();
    if (na != ma)
      throw WrongDim("SolveLU(M, Y)",
		     "The matrix must be squared.");

    CheckDim(M, Y, "SolveLU(M, Y)");
#endif

    // aborting computation if the matrix is sparse
    if (Storage0::Sparse)
      throw WrongArgument("SolveLU", "This function is intended to dense"
                          " matrices only and not to sparse matrices");
    
    // Forward substitution.
    for (i = 0; i < ma; i++)
      {
	SetComplexZero(temp);
	for (k = 0; k < i; k++)
	  temp += M(i, k) * Y(k);
	Y(i) = (Y(i) - temp) / M(i, i);
      }
    // Back substitution.
    for (i = ma - 2; i > -1; i--)
      {
	SetComplexZero(temp);
	for (k = i + 1; k < ma; k++)
	  temp += M(i, k) * Y(k);
	Y(i) -= temp;
      }
  }


  // SOLVELU //
  /////////////

  
  ///////////////////
  // GetAndSolveLU //


  //! Solves a linear system using LU factorization.
  /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
    \f$ X \f$ and \f$ Y \f$ are vectors.
    \param[in] M the matrix of the linear system, to be factorized in LU
    form. On exit, \a M contains its LU factorization.
    \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
    solution \f$ X \f$ of the system.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetAndSolveLU(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                     Vector<T1, Storage1, Allocator1>& Y)
  {
#ifdef SELDON_WITH_LAPACK
    Vector<int> P;
    GetLU(M, P);
    SolveLU(M, P, Y);
#else
    GetLU(M);
    SolveLU(M, Y);
#endif
  }


  // GetAndSolveLU //
  ///////////////////

  
  ///////////
  // Solve //
  
  
  //! solves by triangular upper part or lower part of A
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void
  SolveTriangular(const SeldonUplo& Uplo,
		  const SeldonTranspose& TransA,
		  const SeldonDiag& DiagA,
		  const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
		  const Vector<T1, VectFull, Allocator1>& X,
		  Vector<T1, VectFull, Allocator1>& Y)
  {
    int ma = A.GetM();
    
#ifdef SELDON_CHECK_DIMENSIONS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("Solve(UpLo, TransA, DiagA, A, X, Y)",
		     "The matrix must be squared.");

    CheckDim(A, X, "Solve(UpLo, TransA, DiagA, A, X, Y)");
    CheckDim(A, Y, "Solve(UpLo, TransA, DiagA, A, X, Y)");
#endif
    
    Copy(X, Y);
    
    T0* data = A.GetData();
    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T0 zero; SetComplexZero(zero);
    T1 val;
    if (Uplo.Lower())
      {
        if (TransA.NoTrans())
          {
            if (DiagA.NonUnit())
              {
                for (int i = 0; i < ma; i++)
                  {
                    val = Y(i);
                    long j = ptr[i];
                    while (ind[j] < i)
                      {
                        val -= data[j]*Y(ind[j]);
                        j++;
                      }
                    
#ifdef SELDON_CHECK_BOUNDS
                    if ( (j >= ptr[i+1]) || (ind[j] != i) || (data[j] == zero))
                      throw WrongArgument("Solve", "Matrix must contain"
                                          " a non-null diagonal");
#endif
                    
                    Y(i) = val/data[j];
                  }
              }
            else
              {
                for (int i = 0; i < ma; i++)
                  {
                    val = Y(i);
                    long j = ptr[i];
                    while ((j < ptr[i+1]) && (ind[j] < i))
                      {
                        val -= data[j]*Y(ind[j]);
                        j++;
                      }
                    
                    Y(i) = val;
                  }
              }
          }
        else if (TransA.Trans())
          {
            if (DiagA.NonUnit())
              {
                for (int i = ma-1; i >= 0; i--)
                  {
                    long j = ptr[i+1]-1;
                    while (ind[j] > i)
                      j--;
                    
#ifdef SELDON_CHECK_BOUNDS
                    if ( (j < ptr[i]) || (ind[j] != i) || (data[j] == zero))
                      throw WrongArgument("Solve", "Matrix must contain"
                                          " a non-null diagonal");
#endif
                   
                    Y(i) /= data[j];
                    j--;
                    while (j >= ptr[i])
                      {
                        Y(ind[j]) -= data[j]*Y(i);
                        j--;
                      }
                  }
              }
            else
              {
                for (int i = ma-1; i >= 0; i--)
                  {
                    long j = ptr[i];
                    while ( (j < ptr[i+1]) && (ind[j] < i))
                      {
                        Y(ind[j]) -= data[j]*Y(i);
                        j++;
                      }
                  }
              }
          }
        else
          {
            if (DiagA.NonUnit())
              {
                for (int i = ma-1; i >= 0; i--)
                  {
                    long j = ptr[i+1]-1;
                    while (ind[j] > i)
                      j--;
                    
#ifdef SELDON_CHECK_BOUNDS
                    if ( (j < ptr[i]) || (ind[j] != i) || (data[j] == zero))
                      throw WrongArgument("Solve", "Matrix must contain"
                                          " a non-null diagonal");
#endif
                   
                    Y(i) /= conjugate(data[j]);
                    j--;
                    while (j >= ptr[i])
                      {
                        Y(ind[j]) -= conjugate(data[j])*Y(i);
                        j--;
                      }
                  }
              }
            else
              {
                for (int i = ma-1; i >= 0; i--)
                  {
                    long j = ptr[i];
                    while ((j < ptr[i+1]) && (ind[j] < i))
                      {
                        Y(ind[j]) -= conjugate(data[j])*Y(i);
                        j++;
                      }
                  }
              }
          }
      }
    else
      {
        if (TransA.NoTrans())
          {
            if (DiagA.NonUnit())
              {
                for (int i = ma-1; i >= 0; i--)
                  {
                    val = Y(i);
                    long j = ptr[i+1]-1;
                    while (ind[j] > i)
                      {
                        val -= data[j]*Y(ind[j]);
                        j--;
                      }
                    
#ifdef SELDON_CHECK_BOUNDS
                    if ( (j < ptr[i]) || (ind[j] != i) || (data[j] == zero))
                      throw WrongArgument("Solve", "Matrix must contain"
                                          " a non-null diagonal");
#endif
                   
                    Y(i) = val/data[j];
                  }
              }
            else
              {
                for (int i = ma-1; i >= 0; i--)
                  {
                    val = Y(i);
                    long j = ptr[i+1]-1;
                    while ( (j >= ptr[i]) && (ind[j] > i))
                      {
                        val -= data[j]*Y(ind[j]);
                        j--;
                      }
                    
                    Y(i) = val;
                  }
              }
          }
        else if (TransA.Trans())
          {
            if (DiagA.NonUnit())
              {
                for (int i = 0; i < ma; i++)
                  {
                    long j = ptr[i];
                    while (ind[j] < i)
                      j++;
                    
#ifdef SELDON_CHECK_BOUNDS
                    if ( (j >= ptr[i+1]) || (ind[j] != i) || (data[j] == zero) )
                      throw WrongArgument("Solve", "Matrix must contain"
                                          " a non-null diagonal");
#endif
                    
                    Y(i) /= data[j];
                    j++;
                    while (j < ptr[i+1])
                      {
                        Y(ind[j]) -= data[j]*Y(i);
                        j++;
                      }
                  }
              }
            else
              {
                for (int i = 0; i < ma; i++)
                  {
                    long j = ptr[i+1]-1;
                    while ( (j >= ptr[i]) && (ind[j] > i))
                      {
                        Y(ind[j]) -= data[j]*Y(i);
                        j--;
                      }         
                  }           
              }
          }
        else
          {
            if (DiagA.NonUnit())
              {
                for (int i = 0; i < ma; i++)
                  {
                    long j = ptr[i];
                    while (ind[j] < i)
                      j++;
                    
#ifdef SELDON_CHECK_BOUNDS
                    if ( (j >= ptr[i+1]) || (ind[j] != i) || (data[j] == zero) )
                      throw WrongArgument("Solve", "Matrix must contain"
                                          " a non-null diagonal");
#endif
                    
                    Y(i) /= conjugate(data[j]);
                    j++;
                    while (j < ptr[i+1])
                      {
                        Y(ind[j]) -= conjugate(data[j])*Y(i);
                        j++;
                      }
                  }
              }
            else
              {
                for (int i = 0; i < ma; i++)
                  {
                    long j = ptr[i+1]-1;
                    while ( (j >= ptr[i]) && (ind[j] > i))
                      {
                        Y(ind[j]) -= conjugate(data[j])*Y(i);
                        j--;
                      }         
                  }           
              }
          }
      }
  }
  
  
  // Solve //
  ///////////
  

  //////////////
  // CHECKDIM //


  //! Checks the compatibility of the dimensions.
  /*! Checks that M X + Y -> Y is possible according to the dimensions of
    the matrix M and the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param M matrix.
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function)
  {
    if (X.GetM() != M.GetN() || Y.GetM() != M.GetM())
      {
	
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
        string Mchar = to_str(&M), Xchar = to_str(&X), Ychar = to_str(&Y);
#else
        string Mchar("M"), Xchar("X"), Ychar("Y");
#endif
        
        throw WrongDim(function,
                       string("Operation M X + Y -> Y not permitted:")
                       + string("\n     M (") + Mchar + string(") is a ")
                       + to_str(M.GetM()) + string(" x ") + to_str(M.GetN())
                       + string(" matrix;\n     X (") + Xchar
                       + string(") is vector of length ")
                       + to_str(X.GetLength()) + string(";\n     Y (")
                       + Ychar + string(") is vector of length ")
                       + to_str(Y.GetLength()) + string("."));
      }
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that M X + Y -> Y is possible according to the dimensions of
    the matrix M and the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param M matrix.
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, RowMajorCollection, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Collection, Allocator2>& Y,
		string function)
  {
    if (X.GetNvector() != M.GetNmatrix() || Y.GetNvector() != M.GetMmatrix())
      {
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
	string Mchar = to_str(&M), Xchar = to_str(&X), Ychar = to_str(&Y);
#else
	string Mchar("M"), Xchar("X"), Ychar("Y");
#endif
	
	throw WrongDim(function,
		       string("Operation M X + Y -> Y not permitted:")
		       + string("\n     M (") + Mchar + string(") is a ")
		       + to_str(M.GetM()) + string(" x ") + to_str(M.GetN())
		       + string(" matrix;\n     X (") + Xchar
		       + string(") is vector of length ")
		       + to_str(X.GetNvector()) + string(";\n     Y (")
		       + Ychar + string(") is vector of length ")
		       + to_str(Y.GetNvector()) + string("."));
      }
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that M X + Y -> Y is possible according to the dimensions of
    the matrix M and the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param M matrix.
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, ColMajorCollection, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Collection, Allocator2>& Y,
		string function)
  {
    if (X.GetNvector() != M.GetNmatrix() || Y.GetNvector() != M.GetMmatrix())
      {
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
	string Mchar = to_str(&M), Xchar = to_str(&X), Ychar = to_str(&Y);
#else
	string Mchar("M"), Xchar("X"), Ychar("Y");
#endif
	
	throw WrongDim(function,
		       string("Operation M X + Y -> Y not permitted:")
		       + string("\n     M (") + Mchar + string(") is a ")
		       + to_str(M.GetM()) + string(" x ") + to_str(M.GetN())
		       + string(" matrix;\n     X (") + Xchar
		       + string(") is vector of length ")
		       + to_str(X.GetNvector()) + string(";\n     Y (")
		       + Ychar + string(") is vector of length ")
		       + to_str(Y.GetNvector()) + string("."));
      }
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that M X + Y -> Y is possible according to the dimensions of
    the matrix M and the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param M matrix.
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function)
  {
    if (X.GetM() != M.GetN() || Y.GetM() != M.GetM())
      {
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
	string Mchar = to_str(&M), Xchar = to_str(&X), Ychar = to_str(&Y);
#else
	string Mchar("M"), Xchar("X"), Ychar("Y");
#endif
	
	throw WrongDim(function,
		       string("Operation M X + Y -> Y not permitted:")
		       + string("\n     M (") + Mchar + string(") is a ")
		       + to_str(M.GetM()) + string(" x ") + to_str(M.GetN())
		       + string(" matrix;\n     X (") + Xchar
		       + string(") is vector of length ")
		       + to_str(X.GetLength()) + string(";\n     Y (")
		       + Ychar + string(") is vector of length ")
		       + to_str(Y.GetLength()) + string("."));
      }
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that M X + Y -> Y is possible according to the dimensions of
    the matrix M and the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param trans status of the matrix M, e.g. transposed.
    \param M matrix.
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "M X + Y -> Y".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const SeldonTranspose& trans,
		const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function, string op)
  {
    if (op == "M X + Y -> Y")
      if (trans.Trans())
	op = string("Operation M' X + Y -> Y not permitted:");
      else if (trans.ConjTrans())
	op = string("Operation M* X + Y -> Y not permitted:");
      else
	op = string("Operation M X + Y -> Y not permitted:");
    else
      op = string("Operation ") + op + string(" not permitted:");

    if (X.GetM() != M.GetN(trans) || Y.GetM() != M.GetM(trans))
      {    
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
	string Mchar = to_str(&M), Xchar = to_str(&X), Ychar = to_str(&Y);
#else
	string Mchar("M"), Xchar("X"), Ychar("Y");
#endif
	
	throw WrongDim(function, op + string("\n     M (") + Mchar
		       + string(") is a ") + to_str(M.GetM()) + string(" x ")
		     + to_str(M.GetN()) + string(" matrix;\n     X (")
		       + Xchar + string(") is vector of length ")
		       + to_str(X.GetLength()) + string(";\n     Y (")
		       + Ychar + string(") is vector of length ")
		       + to_str(Y.GetLength()) + string("."));
      }
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that M X is possible according to the dimensions of
    the matrix M and the vector X. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param M matrix.
    \param X vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed. Default: "M X".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		string function, string op)
  {
    if (X.GetM() != M.GetN())
      {
#ifndef SELDON_WITHOUT_TO_STR_CHECKDIM
	string Mchar = to_str(&M), Xchar = to_str(&X);
#else
	string Mchar("M"), Xchar("X");
#endif
	
	throw WrongDim(function, string("Operation ") + op + " not permitted:"
		       + string("\n     M (") + Mchar + string(") is a ")
		       + to_str(M.GetM()) + string(" x ") + to_str(M.GetN())
		       + string(" matrix;\n     X (") + Xchar
		       + string(") is vector of length ")
		       + to_str(X.GetLength()) + string("."));
      }
  }


  // CHECKDIM //
  //////////////

}  // namespace Seldon.


#endif
