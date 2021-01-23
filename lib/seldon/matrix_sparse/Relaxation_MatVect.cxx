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


#ifndef SELDON_FILE_RELAXATION_MATVECT_CXX

/*
  Functions defined in this file

  SOR(A, X, B, omega, iter, type_ssor)

*/

namespace Seldon
{

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
  void SorVector(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    T1 temp, zero; T3 one;
    T0 ajj;
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

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    typename Matrix<T0, Prop0, RowSparse, Allocator0>::pointer data
      = A.GetData();
    
    // Forward sweep.
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = 0; j < ma; j++)
	  {
	    temp = zero;
            long k = ptr[j];
            while (ind[k] < j)
              {
                temp += data[k] * X(ind[k]);
                k++;
              }
            
#ifdef SELDON_CHECK_BOUNDS
            if ( (k >= ptr[j+1]) || (ind[k] != j) || (data[k] == zero))
              throw WrongArgument("SOR", "Matrix must contain"
                                  " a non-null diagonal");
#endif
            
            ajj = data[k];
            k++;            
            while (k < ptr[j+1])
	      {
                temp += data[k] * X(ind[k]);
                k++;
              }

	    X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
	  }

    // Backward sweep.
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = ma-1 ; j >= 0; j--)
	  {
	    temp = zero;
            long k = ptr[j];
            while (ind[k] < j)
              {
                temp += data[k] * X(ind[k]);
                k++;
              }
            
            ajj = data[k];
            k++;            
            while (k < ptr[j+1])
	      {
                temp += data[k] * X(ind[k]);
                k++;
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
  void SorVector(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega,
		 int iter, int type_ssor)
  {
    T1 temp, zero; T3 one;
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

    // Forward sweep.
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = 0; j < ma; j++)
	  {
	    temp = zero;
            int k = 0;
            while (A.Index(j, k) < j)
              {
                temp += A.Value(j, k) * X(A.Index(j, k));
                k++;
              }
            
            ajj = A.Value(j, k);

#ifdef SELDON_CHECK_BOUNDS
	    if ((A.Index(j, k) != j) || (ajj == zero))
	      throw WrongArgument("SOR", "Matrix must contain"
				  "a non-null diagonal");
#endif
	    
            k++;
	    while (k < A.GetRowSize(j))
	      {
		temp += A.Value(j, k) * X(A.Index(j, k));
                k++;
	      }
            
	    X(j) = (one-omega) * X(j) + omega * (B(j) - temp) / ajj;
	  }
    
    // Backward sweep.
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = ma-1; j >= 0; j--)
	  {
	    temp = zero;
	    int k = 0;
	    while (A.Index(j, k) < j)
              {
                temp += A.Value(j, k) * X(A.Index(j, k));
                k++;
              }
            
            ajj = A.Value(j, k);
            k++;
	    while (k < A.GetRowSize(j))
	      {
		temp += A.Value(j, k) * X(A.Index(j, k));
                k++;
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
  void SorVector(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    T1 temp, zero; T3 one;
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

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T0* data = A.GetData();

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
              
#ifdef SELDON_CHECK_BOUNDS
              if ( (ptr[j] >= ptr[j+1]) || (ind[ptr[j]] != j) || (data[ptr[j]] == zero))
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              ajj = data[ptr[j]];

	      for (long k = ptr[j]+1; k < ptr[j+1]; k++)
                temp += data[k] * X(ind[k]);
              
	      X(j) = coef * ajj * X(j) + B(j) - temp;
	    }
          
          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              X(j) *= omega / data[ptr[j]];
	      for (long k = ptr[j]+1; k < ptr[j+1]; k++)
                X(ind[k]) -= data[k]*X(j);
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
	      for (long k = ptr[j]+1; k < ptr[j+1]; k++)
		X(ind[k]) -= data[k]*X(j);
              
              X(j) = B(j) + coef * data[ptr[j]] * X(j);
	    }

          // Then we solve (D/omega - U) X = X.
	  for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              ajj = data[ptr[j]];
	      for (long k = ptr[j]+1; k < ptr[j+1]; k++)
                temp += data[k]*X(ind[k]);
	      
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
  void SorVector(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    T1 temp, zero; T3 one;
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

    // Forward sweep.
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    T3 coef = (one - omega) / omega;
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
          // First we do X = (U + (1-omega)/omega D) X + B.
	  for (int j = 0; j < ma; j++)
	    {
	      temp = zero;
	      ajj = A.Value(j, 0);
	      
#ifdef SELDON_CHECK_BOUNDS
	      if ((A.Index(j, 0) != j) || (ajj == zero))
		throw WrongArgument("SOR", "Matrix must contain"
				    "a non-null diagonal");
#endif
	      
	      for (int k = 1; k < A.GetRowSize(j); k++)
                temp += A.Value(j, k) * X(A.Index(j, k));
              
	      X(j) = coef * ajj * X(j) + B(j) - temp;
	    }
          
          // Then we solve (D/omega - L) X = X.
	  for (int j = 0; j < ma; j++)
	    {
              X(j) *= omega / A.Value(j, 0);
	      for (int k = 1; k < A.GetRowSize(j); k++)
		X(A.Index(j, k)) -= A.Value(j, k)*X(j);
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
	      ajj = A.Value(j, 0);
	      for (int k = 1; k < A.GetRowSize(j); k++)
		X(A.Index(j, k)) -= A.Value(j, k) * X(j);
              
              X(j) = B(j) + coef * ajj * X(j);
	    }
          
          // Then we solve (D/omega - U) X = X.
	  for (int j = ma-1; j >= 0; j--)
	    {
	      temp = zero;
	      ajj = A.Value(j, 0);
	      for (int k = 1; k < A.GetRowSize(j); k++)
		temp += A.Value(j, k)*X(A.Index(j, k));
              
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
  void SorVector(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    T1 zero; T3 one;
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

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T0* data = A.GetData();

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
              long k = ptr[j];
              while (ind[k] < j)
                {
                  X(ind[k]) -= data[k]*X(j);
                  k++;
                }
              
#ifdef SELDON_CHECK_BOUNDS
              if ( (k >= ptr[j+1]) || (ind[k] != j) || (data[k] == zero))
                throw WrongArgument("SOR", "Matrix must contain"
                                    " a non-null diagonal");
#endif
              	      
	      ajj = data[k];
              X(j) = B(j) + coef * ajj * X(j);
	    }
          
          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              long k = ptr[j];
              while (ind[k] < j)
                k++;

              ajj = data[k]; k++;
              X(j) *= omega/ajj;
              while (k < ptr[j+1])
                {
                  X(ind[k]) -= data[k]*X(j);
                  k++;
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
              long k = ptr[j+1]-1;
              while (ind[k] > j)
                {
                  X(ind[k]) -= data[k]*X(j);
                  k--;
                }
              
#ifdef SELDON_CHECK_BOUNDS
              if ( (k < ptr[j]) || (ind[k] != j) || (data[k] == zero) )
                throw WrongArgument("SOR", "Matrix must contain"
                                    "a non-null diagonal");
#endif
              
              X(j) = B(j) + coef*data[k]*X(j);
            }
                    
          // Then we solve (D/omega - U) X = X.
          for (int j = ma-1; j >= 0; j--)
            {
              long k = ptr[j+1]-1;
              while (ind[k] > j)
                k--;
              
              X(j) *= omega/data[k];
              k--;
              while (k >= ptr[j])
                {
                  X(ind[k]) -= data[k]*X(j);
                  k--;
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
  void SorVector(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    T1 zero; T3 one;
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
              int k = 0;
              while (A.Index(j, k) < j)
                {
                  X(A.Index(j, k)) -= A.Value(j, k) * X(j);
                  k++;
                }
              
	      ajj = A.Value(j, k);

#ifdef SELDON_CHECK_BOUNDS
	      if ((A.Index(j, k) != j) || (ajj == zero))
		throw WrongArgument("SOR", "Matrix must contain"
				    "a non-null diagonal");
#endif
	      
              X(j) = B(j) + coef * ajj * X(j);
	    }
          
          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              int k = 0;
              while (A.Index(j, k) < j)
                k++;

              ajj = A.Value(j, k); k++;
              X(j) *= omega/ajj;
              while (k < A.GetColumnSize(j))
                {
                  X(A.Index(j, k)) -= A.Value(j, k) * X(j);
                  k++;
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
              int k = A.GetColumnSize(j) - 1;
              while (A.Index(j, k) > j)
                {
                  X(A.Index(j, k)) -= A.Value(j, k) * X(j);
                  k--;
                }
              
              X(j) = B(j) + coef*A.Value(j, k)*X(j);
            }
                    
          // Then we solve (D/omega - U) X = X.
          for (int j = ma-1; j >= 0; j--)
            {
              int k = A.GetColumnSize(j) - 1;
              while (A.Index(j, k) > j)
                k--;
              
              X(j) *= omega/A.Value(j, k);
              k--;
              while (k >= 0)
                {
                  X(A.Index(j, k)) -= A.Value(j, k) * X(j);
                  k--;
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
  void SorVector(const Matrix<T0, Prop0, ColSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    T1 temp, zero; T3 one;
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

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T0* data = A.GetData();
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
#ifdef SELDON_CHECK_BOUNDS
              if ( (ptr[j] >= ptr[j+1]) || (ind[ptr[j+1]-1] != j)
		   || (data[ptr[j+1]-1] == zero) )
                throw WrongArgument("SOR", "Matrix must contain"
                                    "a non-null diagonal");
#endif
              
              ajj = data[ptr[j+1]-1];

	      for (long k = ptr[j]; k < ptr[j+1]-1; k++)
                X(ind[k]) -= data[k] * X(j);
              
	      X(j) = coef * ajj * X(j) + B(j);
	    }

          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              ajj = data[ptr[j+1]-1];
	      for (long k = ptr[j]; k < ptr[j+1]-1; k++)
                X(j) -= data[k] * X(ind[k]);
              
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
              ajj = data[ptr[j+1]-1];
	      for (long k = ptr[j]; k < ptr[j+1]-1; k++)
		temp -= data[k] * X(ind[k]);
              
              X(j) = B(j) + coef * ajj * X(j) + temp;
	    }
          
          // Then we solve (D/omega - U) X = X.
	  for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              ajj = data[ptr[j+1]-1];
              X(j) *= omega / ajj;
	      for (long k = ptr[j]; k < ptr[j+1]-1; k++)
                X(ind[k]) -= data[k] * X(j);
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
  void SorVector(const Matrix<T0, Prop0, ArrayColSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    T1 temp, zero; T3 one;
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
              int kmax = A.GetColumnSize(j)-1;
              ajj = A.Value(j, kmax);

#ifdef SELDON_CHECK_BOUNDS
	      if ((A.Index(j, kmax) != j) || (ajj == zero))
		throw WrongArgument("SOR", "Matrix must contain"
				    "a non-null diagonal");
#endif
	      
	      for (int k = 0; k < kmax; k++)
                X(A.Index(j, k)) -= A.Value(j, k) * X(j);
              
	      X(j) = coef * ajj * X(j) + B(j);
	    }

          // Then we solve (D/omega - L) X = X
	  for (int j = 0; j < ma; j++)
	    {
              int kmax = A.GetColumnSize(j)-1;
              ajj = A.Value(j, kmax);
	      for (int k = 0; k < kmax; k++)
                X(j) -= A.Value(j, k) * X(A.Index(j, k));
              
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
              int kmax = A.GetColumnSize(j)-1;
              ajj = A.Value(j, kmax);
	      for (int k = 0; k < kmax; k++)
		temp -= A.Value(j, k) * X(A.Index(j, k));
              
              X(j) = B(j) + coef * ajj * X(j) + temp;
	    }
          
          // Then we solve (D/omega - U) X = X.
	  for (int j = ma-1; j >= 0; j--)
	    {
              temp = zero;
              int kmax = A.GetColumnSize(j)-1;
              ajj = A.Value(j, kmax);
              X(j) *= omega / ajj;
	      for (int k = 0; k < kmax; k++)
                X(A.Index(j, k)) -= A.Value(j, k) * X(j);
	    }
	}
  }


  /*********************************
   * S.O.R with transpose matrices *
   *********************************/
  
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);
    
    Matrix<T0, Prop0, ColSparse, Allocator0> Ac;
    Ac.SetData(A.GetN(), A.GetM(), A.GetDataSize(),
	       A.GetData(), A.GetPtr(), A.GetInd());
    
    SorVector(Ac, X, B, omega, iter, type_ssor);
    Ac.Nullify();
  }
  
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);

    Matrix<T0, Prop0, RowSparse, Allocator0> Ac;
    Ac.SetData(A.GetN(), A.GetM(), A.GetDataSize(),
	       A.GetData(), A.GetPtr(), A.GetInd());
    
    SorVector(Ac, X, B, omega, iter, type_ssor);
    Ac.Nullify();
  }

  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    SorVector(A, X, B, omega, iter, type_ssor);
  }


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ColSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    SorVector(A, X, B, omega, iter, type_ssor);
  }


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);

    Matrix<T0, Prop0, ArrayColSparse, Allocator0> Ac;
    Ac.SetData(A.GetN(), A.GetM(), A.GetData());
    
    SorVector(Ac, X, B, omega, iter, type_ssor);
    Ac.Nullify();
  }


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    if (transM.NoTrans())
      return SorVector(A, X, B, omega, iter, type_ssor);

    Matrix<T0, Prop0, ArrayRowSparse, Allocator0> Ac;
    Ac.SetData(A.GetN(), A.GetM(), A.GetData());
    
    SorVector(Ac, X, B, omega, iter, type_ssor);
    Ac.Nullify();
  }


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    SorVector(A, X, B, omega, iter, type_ssor);
  }


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const Matrix<T0, Prop0, ArrayColSymSparse, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor)
  {
    SorVector(A, X, B, omega, iter, type_ssor);
  }

} // end namespace

#define SELDON_FILE_RELAXATION_MATVECT_CXX
#endif
