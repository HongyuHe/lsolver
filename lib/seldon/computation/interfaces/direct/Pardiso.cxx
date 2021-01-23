// Copyright (C) 2003-20013 Marc Duruflé
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


#ifndef SELDON_FILE_PARDISO_CXX

#include "Pardiso.hxx"

namespace Seldon
{
  
  //! default constructor
  template<class T>
  MatrixPardiso<T>::MatrixPardiso()
  {
    size_matrix = 0;
    mtype = 0;
    maxfct = 1;
    mnum = 1;
    msglvl = 0;
    info_facto = 0;
    type_ordering = 2;
  }
    
  
  //! destructor
  template<class T>
  MatrixPardiso<T>::~MatrixPardiso()
  {
    Clear();
  }
  

  template<class T>
  bool MatrixPardiso<T>::UseInteger8() const  
  {
    if (sizeof(pardiso_int_t) == 8)
      return true;
    
    return false;
  }


  //! LU factorization is cleared
  template<class T>
  void MatrixPardiso<T>::Clear()
  {
    // releasing internal memory
    if (size_matrix > 0)
      {
        double ddum;
        pardiso_int_t nrhs = 0;
        pardiso_int_t phase = 0, error;
        // MKL version :
        call_pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size_matrix,
                     valA.GetData(), ptrA.GetData(), indA.GetData(),
                     perm.GetData(), &nrhs, iparm, &msglvl,
                     &ddum, &ddum, &error);
        
        // recent version :
        /* call_pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size_matrix,
                     valA.GetData(), ptrA.GetData(), indA.GetData(),
                     perm.GetData(), &nrhs, iparm, &msglvl,
                     &ddum, &ddum, &error, dparm); */
        
        size_matrix = 0;
        
        // and values of matrix
        ptrA.Clear();
        indA.Clear();
        valA.Clear();
      }    
  }
  

  //! sets ordering algorithm to use
  template<class T>
  void MatrixPardiso<T>::SelectOrdering(int type)
  {
    type_ordering = 2;
  }
  
  
  //! Sets ordering
  template<class T>
  void MatrixPardiso<T>::SetPermutation(const IVect& permut)
  {
    perm.Reallocate(permut.GetM());
    for (int i = 0; i < perm.GetM(); i++)
      perm(i) = permut(i) + 1;
  }
  
  
  //! Disables output messages
  template<class T>
  void MatrixPardiso<T>::HideMessages()
  {
    msglvl = 0;
  }
  
  
  //! Enables output messages
  template<class T>
  void MatrixPardiso<T>::ShowMessages()
  {
    msglvl = 1;
  }
  
  
  //! returns the error code obtained during numerical factorization
  template<class T>
  int MatrixPardiso<T>::GetInfoFactorization() const
  {
    return info_facto;
  }

  
  //! returns the size of memory used by numerical factorization in bytes
  template<class T>
  size_t MatrixPardiso<T>::GetMemorySize() const
  {
    size_t taille = sizeof(pardiso_int_t)*ptrA.GetM();
    if (size_matrix <= 0)
      return taille;
    
    taille += sizeof(pardiso_int_t)*indA.GetM();
    taille += sizeof(T)*valA.GetM();
    taille += sizeof(pardiso_int_t)*perm.GetM();
    taille += 1024*size_t(iparm[15]+iparm[16]);
    return taille;
  }
  
  
  //! performs analysis and factorization of matrix mat
  template<class T> template<class T0, class Prop, class Storage, class Allocator>
  void MatrixPardiso<T>::FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator>& mat,
                                         bool keep_matrix)
  {
    // previous factorization is cleared if present
    Clear();

    // conversion of sparse matrix into CSR form
    Prop prop;
    ConvertToCSR(mat, prop, ptrA, indA, valA);
    if (!keep_matrix)
      mat.Clear();
    
    FactorizeCSR(IsSymmetricMatrix(mat));
  }


  //! Performs analysis and factorisation of a matrix given in CSR form
  template<class T>
  void MatrixPardiso<T>::FactorizeCSR(Vector<pardiso_int_t>& Ptr, 
				      Vector<pardiso_int_t>& IndCol,
				      Vector<T>& Values, bool sym)
  {
    // previous factorization is cleared if present
    Clear();
    
    // structures are stored internally
    ptrA.SetData(Ptr); Ptr.Nullify();
    indA.SetData(IndCol); IndCol.Nullify();
    valA.SetData(Values); Values.Nullify();

    // then factorization is performed
    FactorizeCSR(sym);
  }


  template<class T>
  void MatrixPardiso<T>::FactorizeCSR(bool sym)
  {
    // checking that the matrix is non-empty
    size_matrix = ptrA.GetM()-1;
    if (size_matrix <= 0)
      return;
    
    // matrix type mtype
    double ddum;
    pardiso_int_t nrhs = 0;
    mtype = 0;
    if (IsComplexNumber(T()))
      {
        if (sym)
          mtype = 6;
        else
          mtype = 13;
      }
    else
      {
        if (sym)
          mtype = -2;
        else
          mtype = 11;
      }
    
    // initialization of solver
    
    // recent version of pardiso :
    /* pardiso_int_t solver = 0, error = 0;
    iparm[0] = 0; iparm[2] = 1;
    pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

    if (error != 0)
      {
        cout << "Error in the license file of Pardiso" << endl;
        if (error == -10)
          cout << "No license file pardiso.lic found" << endl;
        else if (error == -11)
          cout << "License is expired" << endl;
        else if (error == -12)
          cout << "Wrong username or hostname" << endl;
        
        abort();
        } */
    
    // MKL version of pardiso :
    pardiso_int_t error = 0; int mtype_ = mtype;
    int iparm_[64]; iparm_[0] = 0; iparm_[2] = 1;
    pardisoinit(pt, &mtype_, iparm_);
    
    for (int i = 0; i < 64; i++)
      iparm[i] = iparm_[i];
    
    for (int i = 0; i < ptrA.GetM(); i++)
      ptrA(i)++;
    
    for (int i = 0; i < indA.GetM(); i++)
      indA(i)++;
    
    iparm[4] = 0;
    if (perm.GetM() == size_matrix)
      iparm[4] = 1;
    else
      iparm[1] = type_ordering;
    
    // reordering and symbolic factorization
    pardiso_int_t phase = 11;
    // MKL version
    call_pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size_matrix,
                 valA.GetData(), ptrA.GetData(), indA.GetData(), 
                 perm.GetData(), &nrhs, iparm, &msglvl,
                 &ddum, &ddum, &error);
    
    // recent version
    /* call_pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size_matrix,
                 valA.GetData(), ptrA.GetData(), indA.GetData(), 
                 perm.GetData(), &nrhs, iparm, &msglvl,
                 &ddum, &ddum, &error, dparm); */
    
    if (error != 0)
      {
        info_facto = error;
        return;
      }

    if (msglvl >= 1)
      {
        cout << "Reordering completed" << endl;
        cout << "Number of non-zero factors  = " << iparm[17] << endl;
        cout << "Number of factorization MFLOPS = " << iparm[18] << endl;
      }
    
    // numerical factorization
    phase = 22;
    // MKL version
    call_pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size_matrix,
                 valA.GetData(), ptrA.GetData(), indA.GetData(), 
                 perm.GetData(), &nrhs, iparm, &msglvl,
                 &ddum, &ddum, &error);
    
    // recent version
    /*call_pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size_matrix,
                 valA.GetData(), ptrA.GetData(), indA.GetData(), 
                 perm.GetData(), &nrhs, iparm, &msglvl,
                 &ddum, &ddum, &error, dparm); */
    
    if (error != 0)
      {
        info_facto = error;
        cout << "Error during factorization " << error << endl;
        return;
      }
    
    if (msglvl >= 1)
      cout << "Factorization completed" << endl;
  }
  
    
  //! solves A x = b, x contains the source b on input, the solution x on output
  template<class T> template<class Allocator2>
  void MatrixPardiso<T>::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }
  
  
  //! solves A x = b or A^T x = b
  //! x contains the source b on input, the solution x on output
  template<class T> template<class Allocator2>
  void MatrixPardiso<T>::Solve(const SeldonTranspose& TransA,
                               Vector<T, VectFull, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), 1);
  }


  //! solves A x = b or A^T x = b
  //! x contains the source b on input, the solution x on output
  template<class T>
  void MatrixPardiso<T>::Solve(const SeldonTranspose& TransA,
                               T* x_ptr, int nrhs_)
  {
    if (size_matrix <= 0)
      return;
    
    Matrix<T, General, ColMajor> b(size_matrix, nrhs_);
    for (int i = 0; i < size_matrix; i++)
      for (int j = 0; j < nrhs_; j++)
        b(i, j) = x_ptr[i+j*size_matrix];
    
    pardiso_int_t nrhs = nrhs_;      
    
    pardiso_int_t phase = 33, error;
    iparm[5] = 0;
    if (TransA.Trans())
      iparm[11] = 2;
    else
      iparm[11] = 0;
    
    // MKL version
    call_pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size_matrix,
                 valA.GetData(), ptrA.GetData(), indA.GetData(), 
                 perm.GetData(), &nrhs, iparm, &msglvl,
                 b.GetData(), x_ptr, &error);
    
    // recent version
    /* call_pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size_matrix,
                 valA.GetData(), ptrA.GetData(), indA.GetData(), 
                 perm.GetData(), &nrhs, iparm, &msglvl,
                 b.GetData(), x.GetData(), &error, dparm); */
  }
  

  //! solves A x = b or A^T x = b (x, b are matrices)
  //! x contains the source b on input, the solution x on output  
  template<class T>
  template<class Allocator2, class Prop>
  void MatrixPardiso<T>::Solve(const SeldonTranspose& TransA,
                               Matrix<T, Prop, ColMajor, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), x.GetN());
  }


  //! Factorization of a matrix of same type T as for the Pardiso object 
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPardiso<T>& mat_lu, bool keep_matrix, T& x)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  

  //! Factorization of a complex matrix with a real Pardiso object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPardiso<T>& mat_lu,
             bool keep_matrix, complex<T>& x)
  {
    throw WrongArgument("GetLU(Matrix<complex<T> >& A, MatrixPardiso<T>& mat_lu, bool)",
			"The LU matrix must be complex");
  }


  //! Factorization of a real matrix with a complex Pardiso object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPardiso<complex<T> >& mat_lu,
             bool keep_matrix, T& x)
  {
    throw WrongArgument("GetLU(Matrix<T>& A, MatrixPardiso<complex<T> >& mat_lu, bool)",
			"The sparse matrix must be complex");
  }
  

  //! Factorization of a general matrix with Pardiso
  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixPardiso<T>& mat_lu,
	     bool keep_matrix)
  {
    // we check if the type of non-zero entries of matrix A
    // and of the Pardiso object (T) are different
    // we call one of the GetLUs written above
    // such a protection avoids to compile the factorisation of a complex
    // matrix with a real Pardiso object
    typename Matrix<T0, Prop, Storage, Allocator>::entry_type x;
    GetLU(A, mat_lu, keep_matrix, x);
  }


  //! LU resolution with a vector whose type is the same as for Pardiso object
  template<class T, class Allocator>
  void SolveLU(MatrixPardiso<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  //! LU resolution with a vector whose type is the same as for Pardiso object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPardiso<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! LU resolution with a matrix whose type is the same as for Pardiso object
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixPardiso<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(SeldonNoTrans, x);
  }


  //! LU resolution with a matrix whose type is the same as for Pardiso object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Allocator, class Prop>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPardiso<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }

  
  //! Solves A x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(MatrixPardiso<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x)
  {
    Matrix<double, General, ColMajor> y(x.GetM(), 2);
    
    for (int i = 0; i < x.GetM(); i++)
      {
	y(i, 0) = real(x(i));
	y(i, 1) = imag(x(i));
      }
    
    SolveLU(mat_lu, y);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<double>(y(i, 0), y(i, 1));
  }
  

  //! Solves A x = b or A^T x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPardiso<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x)
  {
    Matrix<double, General, ColMajor> y(x.GetM(), 2);
    
    for (int i = 0; i < x.GetM(); i++)
      {
	y(i, 0) = real(x(i));
	y(i, 1) = imag(x(i));
      }
    
    SolveLU(TransA, mat_lu, y);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<double>(y(i, 0), y(i, 1));    
  }
  

  //! Solves A x = b, where A is complex and x is real => Forbidden  
  template<class Allocator>
  void SolveLU(MatrixPardiso<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixPardiso<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }

  
  //! Solves A x = b or A^T x = b, where A is complex and x is real => Forbidden  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPardiso<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixPardiso<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }
  
}

#define SELDON_FILE_PARDISO_CXX
#endif
