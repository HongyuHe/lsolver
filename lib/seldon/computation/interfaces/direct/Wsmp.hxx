// Copyright (C) 2015-2015 Marc Duruflé
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

#ifndef SELDON_FILE_WSMP_HXX

extern "C"
{
  void wsetmaxthrds_(int*);
  void wsmp_initialize_();
  void wsmp_clear_();

#ifndef SELDON_WITH_MPI
  void wssmp_(int*, int*, int*, double*, double*, int*, int*,
              double*, int*, int*, double*, int*, int*, int*, double*);

  void zssmp_(int*, int*, int*, void*, void*, int*, int*,
              void*, int*, int*, void*, int*, int*, int*, double*);

  void wgsmp_(int*, int*, int*, double*, double*, int*, int*,
              double*, int*, double*);
  
  void zgsmp_(int*, int*, int*, void*, void*, int*, int*,
              void*, int*, void*);
#else
  void pwssmp_(int*, int*, int*, double*, double*, int*, int*,
               double*, int*, int*, double*, int*, int*, int*, double*);
  
  void pzssmp_(int*, int*, int*, void*, void*, int*, int*,
               void*, int*, int*, void*, int*, int*, int*, double*);

  void pwgsmp_(int*, int*, int*, double*, double*, int*, int*,
               double*, int*, double*);
  
  void pzgsmp_(int*, int*, int*, void*, void*, int*, int*,
               void*, int*, void*);
#endif
}

namespace Seldon
{

  template<class T>
  class MatrixWsmp : public VirtualSparseDirectSolver<T>
  {
  protected:    
    int n;
    IVect iparm; Vector<double> dparm;
    bool refine_solution, use_pivoting;
    bool cholesky, symmetric, distributed;
    IVect Ptr, IndRow; Vector<T> Val;
    IVect permut, inverse_permut;
    double threshold_pivot;

    void CallWssmp(int*, int*, int*, T*, T*, int*, int*,
                   T*, int*, int*, T*, int*, int*, int*, double*);

    void CallWgsmp(int*, int*, int*, T*, T*, int*, int*,
                   T*, int*, double*);
    
  public:
    MatrixWsmp();
    ~MatrixWsmp();

    bool UseInteger8() const;
    void Clear();
    
    void ShowMessages();
    void HideMessages();
    void ShowFullHistory();
    
    void SetPivotThreshold(double);
    void SetCholeskyFacto(bool chol);    
    int GetInfoFactorization() const;
    
    void RefineSolution();
    void DoNotRefineSolution();
    size_t GetMemorySize() const;
    
    void SetNumberOfThreadPerNode(int nb_threads);

    template<class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T, Symmetric, Storage, Allocator>& mat,
			 bool keep_matrix = false);

    template<class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T, General, Storage, Allocator>& mat,
			 bool keep_matrix = false);

    void FactorizeCSR(Vector<int>& Ptr, Vector<int>& IndRow,
		      Vector<T>& Val, bool sym);

    void FactorizeUnsymmetric();
    void FactorizeSymmetric();
    
    void Solve(Vector<T>& b);
    void Solve(const SeldonTranspose& trans, Vector<T>& b);
    void Solve(const SeldonTranspose&, T* x, int nrhs);
    void Solve(const SeldonTranspose& trans, Matrix<T, General, ColMajor>& b);

#ifdef SELDON_WITH_MPI
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
                           Vector<int>&, Vector<int>&,
                           Vector<T>&, const Vector<int>& glob_number,
                           bool sym, bool keep_matrix = false);

    void SolveDistributed(MPI_Comm& comm_facto,
                          const SeldonTranspose& TransA,
			  Vector<T>& x, const Vector<int>& glob_num);

    void SolveDistributed(MPI_Comm& comm_facto,
                          const SeldonTranspose& TransA,
			  Matrix<T, General, ColMajor>& x, const Vector<int>& glob_num);

    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA,
			  T* x_ptr, int nrhs, const Vector<int>& glob_num);
#endif
    
  };


  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T, Prop, Storage, Allocator>& A, MatrixWsmp<T>& mat_lu,
	     bool keep_matrix = false);

  template<class T, class Allocator>
  void SolveLU(MatrixWsmp<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixWsmp<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixWsmp<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixWsmp<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixWsmp<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixWsmp<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixWsmp<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixWsmp<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x);

}

#define SELDON_FILE_WSMP_HXX
#endif
