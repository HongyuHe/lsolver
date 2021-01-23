// Copyright (C) 2001-2010 Marc Duruflé
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

#ifndef SELDON_FILE_PASTIX_HXX

// including Pastix headers
extern "C"
{
#define _COMPLEX_H

#include "pastix.h"

#undef FLOAT
#undef INT
#undef UINT

}

namespace Seldon
{

  template<class T>
  class MatrixPastix : public VirtualSparseDirectSolver<T>
  {
  protected :
    //! pastix structure
    pastix_data_t* pastix_data;
    //! options (integers)
    pastix_int_t iparm[IPARM_SIZE];
    //! options (floats)
    double dparm[DPARM_SIZE];
    //! number of columns
    pastix_int_t n;
    //! permutation arrays
    Vector<pastix_int_t> perm, invp;
    //! local to global
    Vector<pastix_int_t> col_num;
    //! if true, resolution on several nodes
    bool distributed;
    //! if true, cholesky factorisation is computed
    bool cholesky;
    //! level of display
    int print_level;
    //! if true, solution is refined
    bool refine_solution;
    //! threshold for static pivoting
    double threshold_pivot;
    //! adjust threshold for static pivoting ?
    bool adjust_threshold_pivot;
    MPI_Comm comm_;
    
  public :

    MatrixPastix();
    ~MatrixPastix();

    bool UseInteger8() const;
    void Clear();

    void CallPastix(const MPI_Comm&, pastix_int_t* colptr, pastix_int_t* row,
                    T* val, T* b, pastix_int_t nrhs);

    void HideMessages();
    void ShowMessages();
    void ShowFullHistory();

    void SelectOrdering(int type);
    void SetPermutation(const IVect& permut);
    void SetCholeskyFacto(bool chol);
    
    void SetPivotThreshold(double);
    void RefineSolution();
    void DoNotRefineSolution();
    size_t GetMemorySize() const;
    int GetInfoFactorization() const;
    
    template<class T0, class Prop, class Storage, class Allocator, class Tint>
    void FindOrdering(Matrix<T0, Prop, Storage, Allocator> & mat,
		      Vector<Tint>& numbers, bool keep_matrix = false);

    template<class T0, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, General, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    void FactorizeCSC(Vector<pastix_int_t>& Ptr, Vector<pastix_int_t>& IndRow,
		      Vector<T>& Val, bool sym);

    template<class T0, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, Symmetric, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Allocator2>
    void Solve(Vector<T, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
	       Vector<T, VectFull, Allocator2>& x);

    void Solve(const SeldonTranspose&, T* x_ptr, int nrhs);
    
    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
	       Matrix<T, General, ColMajor, Allocator2>& x);

    template<class Allocator2>
    void Mlt(const SeldonTranspose& TransA,
             Vector<T, VectFull, Allocator2>& x);

    void Mlt(const SeldonTranspose&, T* x_ptr, int nrhs);
    
    void SetNumberOfThreadPerNode(int);
    
    void FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
 				    Vector<int>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    void FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
				    Vector<int64_t>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    template<class Tint>
    void FactorizeParallel(MPI_Comm& comm_facto,
                           Vector<Tint>&, Vector<Tint>&,
                           Vector<T>&,
                           const Vector<int>& glob_num,
                           bool sym, bool keep_matrix = false);

    void FactorizeParallel(MPI_Comm& comm_facto,
                           Vector<pastix_int_t>&,
                           Vector<pastix_int_t>&,
                           Vector<T>&, const Vector<int>& glob_num,
                           bool sym, bool keep_matrix = false);

    template<class Allocator2>
    void SolveDistributed(MPI_Comm& comm_facto,
                          const SeldonTranspose& TransA,
			  Vector<T, VectFull, Allocator2>& x,
                          const Vector<int>& glob_num);

    template<class Allocator2>
    void SolveDistributed(MPI_Comm& comm_facto,
                          const SeldonTranspose& TransA,
			  Matrix<T, General, ColMajor, Allocator2>& x,
                          const Vector<int>& glob_num);
    
    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA,
			  T* x_ptr, int nrhs,
			  const IVect& glob_num);
    
    template<class Allocator2>
    void MltDistributed(MPI_Comm& comm_facto,
			const SeldonTranspose& TransA,
			Vector<T, VectFull, Allocator2>& x,
			const Vector<int>& glob_num);

    template<class Allocator2>
    void MltDistributed(MPI_Comm& comm_facto,
			const SeldonTranspose& TransA,
			Matrix<T, General, ColMajor, Allocator2>& x,
			const Vector<int>& glob_num);
    
    void MltDistributed(MPI_Comm& comm_facto,
			const SeldonTranspose& TransA,
			T* x_ptr, int nrhs,
			const IVect& glob_num);

    
  };

  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixPastix<T>& mat_lu,
	     bool keep_matrix = false);

  template<class T, class Allocator>
  void SolveLU(MatrixPastix<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixPastix<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixPastix<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);
  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixPastix<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);

  template<class T, class Prop, class Storage, class Allocator>
  void GetCholesky(Matrix<T, Prop, Storage, Allocator>& A,
                   MatrixPastix<T>& mat_chol, bool keep_matrix = false);

  template<class T, class Allocator>
  void
  SolveCholesky(const SeldonTranspose& TransA,
                MatrixPastix<T>& mat_chol, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void
  MltCholesky(const SeldonTranspose& TransA,
              MatrixPastix<T>& mat_chol, Vector<T, VectFull, Allocator>& x);
  
}

#define SELDON_FILE_PASTIX_HXX
#endif

