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


#ifndef SELDON_FILE_PARDISO_HXX

// if PARDISO_INTSIZE64 is defined, pardiso_64 is used instead of pardiso
#ifdef PARDISO_INTSIZE64
#define pardiso_int_t int64_t
#define call_pardiso pardiso_64
#else
#define pardiso_int_t int
#define call_pardiso pardiso
#endif

// including Pardiso headers
extern "C"
{
  // MKL version :
  void pardisoinit(void*, int*, int*);
  void pardiso(void   *, int*, int*, int*,
               int*, int*,
               void *, int*, int*, int*,
               int*, int*, int*,
               void *, void *, int*);

  void pardiso_64(void   *, pardiso_int_t*, pardiso_int_t*, pardiso_int_t*,
                  pardiso_int_t*, pardiso_int_t*,
                  void *, pardiso_int_t*, pardiso_int_t*, pardiso_int_t*,
                  pardiso_int_t*, pardiso_int_t*, pardiso_int_t*,
                  void *, void *, pardiso_int_t*);
  
  // recent version :
  /* void pardisoinit (void   *, int *,   int *, int *, double *, int *);
  void pardiso     (void   *, int *,   int *, int *,    int *, int *, 
                    void *, int *,   int *, int *,   int *, int *,
                    int *, void *, void *, int *, double *); */
}

namespace Seldon
{
  template<class T>
  class MatrixPardiso : public VirtualSparseDirectSolver<T>
  {
  protected :
    void* pt[64]; //!< pointer to Pardiso object
    pardiso_int_t mtype; //!< matrix type 
    pardiso_int_t iparm[64]; //!< integer parameters
    double dparm[64];
    pardiso_int_t size_matrix;
    pardiso_int_t maxfct;
    //!< maximum number of factors with identical sparsity structure
    pardiso_int_t mnum; //!< factor number, 1 <= mnum <= maxfct
    pardiso_int_t msglvl; //!< verbosity level
    int info_facto;
    Vector<pardiso_int_t> ptrA, indA; //!< matrix in CSR form
    Vector<T> valA; //!< values of non-zero entries
    Vector<pardiso_int_t> perm; //!< permutation array
    pardiso_int_t type_ordering;
    
  public :
    MatrixPardiso();
    ~MatrixPardiso();

    bool UseInteger8() const;    
    void Clear();
    
    void SelectOrdering(int type);
    void SetPermutation(const IVect& permut);
    
    void HideMessages();
    void ShowMessages();
    
    int GetInfoFactorization() const;
    size_t GetMemorySize() const;
    
    template<class T0, class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator>& mat,
                         bool keep_matrix = false);

    void FactorizeCSR(Vector<pardiso_int_t>& Ptr,
		      Vector<pardiso_int_t>& IndCol,
		      Vector<T>& Values, bool sym);
    
    void FactorizeCSR(bool sym); 
   
    template<class Allocator2>
    void Solve(Vector<T, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
	       Vector<T, VectFull, Allocator2>& x);
    
    void Solve(const SeldonTranspose& TransA, T* x_ptr, int nrhs_);
    
    template<class Allocator2, class Prop>
    void Solve(const SeldonTranspose& TransA,
	       Matrix<T, Prop, ColMajor, Allocator2>& x);
    
  };


  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixPardiso<T>& mat_lu,
	     bool keep_matrix = false);

  template<class T, class Allocator>
  void SolveLU(MatrixPardiso<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPardiso<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixPardiso<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class T, class Allocator, class Prop>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPardiso<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x);
  
  template<class Allocator>
  void SolveLU(MatrixPardiso<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPardiso<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);
  
  template<class Allocator>
  void SolveLU(MatrixPardiso<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);
  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPardiso<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);
  
}

#define SELDON_FILE_PARDISO_HXX
#endif
