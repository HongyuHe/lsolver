// Copyright (C) 2003-2015 Marc Duruflé
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


#ifndef SELDON_FILE_MUMPS_HXX

// including Mumps headers
extern "C"
{
#include "dmumps_c.h"
#include "zmumps_c.h"
}

namespace Seldon
{
  template<class T>
  class TypeMumps
  {
  };


  //! class containing MUMPS data structure
  template<>
  class TypeMumps<double>
  {
  public :
    typedef DMUMPS_STRUC_C data;
    typedef double* pointer;
  };


  //! class containing MUMPS data structure
  template<>
  class TypeMumps<complex<double> >
  {
  public :
    typedef ZMUMPS_STRUC_C data;
    typedef mumps_double_complex* pointer;
  };


  //! object used to solve linear system by calling mumps subroutines
  template<class T>
  class MatrixMumps : public VirtualSparseDirectSolver<T>
  {
  protected :
    int type_ordering; //!< ordering scheme (AMD, Metis, etc)
    bool parallel_ordering;
    //! object containing Mumps data structure
    typename TypeMumps<T>::data struct_mumps;
    //! double* or complex<double>*
    typedef typename TypeMumps<T>::pointer pointer;
    int print_level;
    int info_facto;
    bool out_of_core;
    Vector<MUMPS_INT> num_row_glob, num_col_glob;
    Vector<MUMPS_INT> perm;
    double coef_overestimate;
    double coef_increase_memory;
    double coef_max_overestimate;
    double threshold_pivot;
    
    // internal methods
    void CallMumps();
    void IterateFacto();
    
    void InitMatrix(bool sym, bool dist = false);

  public :
    MatrixMumps();
    ~MatrixMumps();

    bool UseInteger8() const; 
    void Clear();

    void SelectOrdering(int num_ordering);
    void SelectParallelOrdering(int num_ordering);
    void SetPermutation(const IVect& permut);
    void SetPivotThreshold(double);
    
    void HideMessages();
    void ShowMessages();

    void EnableOutOfCore();
    void DisableOutOfCore();
    
    void SetCoefficientEstimationNeededMemory(double);
    void SetMaximumCoefficientEstimationNeededMemory(double);
    void SetIncreaseCoefficientEstimationNeededMemory(double);

    size_t GetMemorySize() const;
    int GetInfoFactorization() const;

    template<class T0, class Prop, class Storage, class Allocator>
    void FindOrdering(Matrix<T0, Prop, Storage, Allocator> & mat,
		      IVect& numbers, bool keep_matrix = false);

    template<class T0, class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    void FactorizeCoordinate1(int n, Vector<MUMPS_INT>& num_row,
			      Vector<MUMPS_INT>& num_col,
			      Vector<T>& values, bool sym);
    
    template<class Prop, class Storage, class Allocator>
    void PerformAnalysis(Matrix<T, Prop, Storage, Allocator> & mat);

    template<class Prop, class Allocator>
    void PerformFactorization(Matrix<T, Prop, RowSparse, Allocator> & mat);

    template<class Prop, class Allocator>
    void PerformFactorization(Matrix<T, Prop, RowSymSparse, Allocator> & mat);

    template<class Prop, class Storage, class Allocator>
    void PerformFactorization(Matrix<T, Prop, Storage, Allocator> & mat);

    template<class Prop1, class Storage1, class Allocator1,
	     class Prop2, class Storage2, class Allocator2>
    void GetSchurMatrix(Matrix<T, Prop1, Storage1, Allocator1>& mat,
			const IVect& num,
			Matrix<T, Prop2, Storage2, Allocator2> & mat_schur,
			bool keep_matrix = false);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
	       Vector<T, VectFull, Allocator2>& x);

    void Solve(const SeldonTranspose&, T* x_ptr, int nrhs);
    
    template<class Allocator2>
    void Solve(Vector<T, VectFull, Allocator2>& x);

    template<class Allocator2, class Prop>
    void Solve(const SeldonTranspose& TransA,
	       Matrix<T, Prop, ColMajor, Allocator2>& x);

#ifdef SELDON_WITH_MPI
    void FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
				    Vector<int>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    void FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
				    Vector<int64_t>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    void FactorizeParallel(MPI_Comm& comm_facto, Vector<long>& Ptr,
                           Vector<MUMPS_INT>& IndRow, Vector<T>& Val,
                           const Vector<int>& glob_number,
                           bool sym, bool keep_matrix = false);

    void PerformAnalysisDistributed(MPI_Comm& comm_facto, Vector<long>& Ptr,
				    Vector<int>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);
    
    void PerformAnalysisDistributed(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
				    Vector<int64_t>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);
    
    void PerformFactorizationDistributed(MPI_Comm& comm_facto, Vector<long>& Ptr,
					 Vector<int>& IndRow, Vector<T>& Val,
					 const Vector<int>& glob_number,
					 bool sym, bool keep_matrix = false);
    
    void PerformFactorizationDistributed(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
					 Vector<int64_t>& IndRow, Vector<T>& Val,
					 const Vector<int>& glob_number,
					 bool sym, bool keep_matrix = false);
    
    template<class Allocator2>
    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA,
			  Vector<T, VectFull, Allocator2>& x,
			  const IVect& glob_num);
    
    template<class Allocator2, class Prop>
    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA,
			  Matrix<T, Prop, ColMajor, Allocator2>& x,
			  const IVect& glob_num);

    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA,
			  T* x_ptr, int nrhs,
			  const IVect& glob_num);
#endif

  };

  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixMumps<T>& mat_lu,
	     bool keep_matrix = false);
  
  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T, Symmetric, Storage, Allocator>& A,
                      MatrixMumps<T>& mat_lu, const IVect& num,
                      MatrixFull& schur_matrix, bool keep_matrix = false);

  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T, General, Storage, Allocator>& A,
                      MatrixMumps<T>& mat_lu, const IVect& num,
                      MatrixFull& schur_matrix, bool keep_matrix = false);

  template<class T, class Allocator>
  void SolveLU(MatrixMumps<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixMumps<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixMumps<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixMumps<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class T, class Allocator, class Prop>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class Prop, class Allocator>
  void SolveLU(MatrixMumps<double>& mat_lu,
	       Matrix<complex<double>, Prop, ColMajor, Allocator>& x);

  template<class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<double>& mat_lu,
	       Matrix<complex<double>, Prop, ColMajor, Allocator>& x);

  template<class Prop, class Allocator>
  void SolveLU(MatrixMumps<complex<double> >& mat_lu,
	       Matrix<double, Prop, ColMajor, Allocator>& x);

  template<class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<complex<double> >& mat_lu,
	       Matrix<double, Prop, ColMajor, Allocator>& x);
  
}

#define SELDON_FILE_MUMPS_HXX
#endif



