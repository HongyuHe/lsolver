// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_SUPERLU_HXX

#ifdef SUPERLU_INTSIZE64
#define _LONGINT
#endif

#include "superlu_interface.h"

// here superlu_int_t is introduced in order to use
// int64_t instead of long int which are different types
// In SparseDistributedSolver int64_t is used, hence
// this need of int64_t
#ifdef SUPERLU_INTSIZE64
#define superlu_int_t int64_t
#undef _LONGINT
#else
#define superlu_int_t int
#endif

namespace Seldon
{

#ifndef SELDON_WITH_SUPERLU_DIST
  void SetComplexOne(superlu::doublecomplex& one);
#endif
  
  //! class interfacing SuperLU functions
  template<class T>
  class MatrixSuperLU_Base : public VirtualSparseDirectSolver<T>
  {
  protected :
#ifndef SELDON_WITH_SUPERLU_DIST
    //! objects of SuperLU
    superlu::SuperMatrix L, U, B;
    superlu::GlobalLU_t Glu; //!< object of SuperLU
    
#ifdef SELDON_WITH_SUPERLU_MT
    superlu::superlumt_options_t options; //!< options
    superlu::Gstat_t stat; //!< statistics
    int_t nprocs;
    
    double diag_pivot_thresh, drop_tol;
    superlu::yes_no_t usepr, refact;
    superlu::fact_t fact;
    
    superlu::SCPformat *Lstore;  //!< object of SuperLU
    superlu::NCPformat *Ustore;  //!< object of SuperLU

#else
    superlu::SCformat *Lstore;  //!< object of SuperLU
    superlu::NCformat *Ustore;  //!< object of SuperLU
    superlu::superlu_options_t options; //!< options
    superlu::SuperLUStat_t stat; //!< statistics
#endif

#else

    int nprow, npcol;
    superlu::gridinfo_t grid;
    superlu::superlu_options_t options;
    superlu::SuperLUStat_t stat;

    superlu::SuperMatrix A;
    superlu::ScalePermstruct_t ScalePermstruct;
    superlu::LUstruct_t LUstruct;
    superlu::SOLVEstruct_t SOLVEstruct;

#endif

    //! permutation array
    Vector<int_t> perm_r, perm_c;
    
    colperm_t permc_spec; //!< ordering scheme
    int_t n; //!< number of rows
    int_t nloc;
    bool display_info; //!< display information about factorization ?
    //! Error code returned by SuperLU.
    int info_facto;

  public :
    MatrixSuperLU_Base();
    ~MatrixSuperLU_Base();

    const Vector<int_t>& GetRowPermutation() const;
    const Vector<int_t>& GetColPermutation() const;

    void Init(int size, int_t& panel_size, int_t& relax);
    void SetNumberOfThreadPerNode(int p);
    
    void SelectOrdering(int type);
    void SetPermutation(const IVect&);

    bool UseInteger8() const;
    void Clear();
    void HideMessages();
    void ShowMessages();

    int GetInfoFactorization() const;
  };


  //! empty matrix
  template<class T>
  class MatrixSuperLU : public MatrixSuperLU_Base<T>
  {
  };


  //! class interfacing SuperLU functions in double precision
  template<>
  class MatrixSuperLU<double> : public MatrixSuperLU_Base<double>
  {
  public:
    MatrixSuperLU() : MatrixSuperLU_Base<double>() {}
    
    size_t GetMemorySize() const;

#ifndef SELDON_WITH_SUPERLU_DIST    
    template<class Prop, class Allocator>
    void GetLU(Matrix<double, Prop, ColSparse, Allocator>& Lmat,
               Matrix<double, Prop, ColSparse, Allocator>& Umat,
               bool permuted = true);

    template<class Prop, class Allocator>
    void GetLU(Matrix<double, Prop, RowSparse, Allocator>& Lmat,
               Matrix<double, Prop, RowSparse, Allocator>& Umat,
               bool permuted = true);
#endif

    template<class T0, class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    void FactorizeCSC(Vector<long>& Ptr, Vector<int>& IndRow,
		      Vector<double>& Val, bool sym);

    template<class Allocator2>
    void Solve(Vector<double, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
               Vector<double, VectFull, Allocator2>& x);

    void Solve(const SeldonTranspose& Trans,
	       double* x_ptr, int nrhs_);
    
    template<class Allocator2>
    void Solve(Matrix<double, General, ColMajor, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
               Matrix<double, General, ColMajor, Allocator2>& x);

#ifdef SELDON_WITH_SUPERLU_DIST
    void FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
				    Vector<int>& IndRow, Vector<double>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    void FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
				    Vector<int64_t>& IndRow, Vector<double>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    template<class Tint>
    void FactorizeParallel(MPI_Comm& comm_facto,
                           Vector<Tint>&, Vector<Tint>&,
                           Vector<double>&,
                           const Vector<int>& glob_num,
                           bool sym, bool keep_matrix = false);
    
    void FactorizeParallel(MPI_Comm& comm_facto,
                           Vector<superlu_int_t>&, Vector<superlu_int_t>&,
                           Vector<double>&,
                           const Vector<int>& glob_num,
                           bool sym, bool keep_matrix = false);
    
    template<class Allocator2>
    void SolveDistributed(MPI_Comm& comm_facto,
                          const SeldonTranspose& TransA,
			  Vector<double, VectFull, Allocator2>& x,
			  const Vector<int>& glob_num);

    template<class Allocator2>
    void SolveDistributed(MPI_Comm& comm_facto,
                          const SeldonTranspose& TransA,
			  Matrix<double, General, ColMajor, Allocator2>& x,
			  const Vector<int>& glob_num);

    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA,
			  double* x_ptr, int nrhs_,
			  const IVect& glob_num);
    
#endif

  };


  //! class interfacing SuperLU functions in complex double precision
  template<>
  class MatrixSuperLU<complex<double> >
    : public MatrixSuperLU_Base<complex<double> >
  {
  public:
    MatrixSuperLU() : MatrixSuperLU_Base<complex<double> >() {}

    size_t GetMemorySize() const;

#ifndef SELDON_WITH_SUPERLU_DIST        
    template<class Prop, class Allocator>
    void GetLU(Matrix<complex<double>, Prop, ColSparse, Allocator>& Lmat,
               Matrix<complex<double>, Prop, ColSparse, Allocator>& Umat,
               bool permuted = true);

    template<class Prop, class Allocator>
    void GetLU(Matrix<complex<double>, Prop, RowSparse, Allocator>& Lmat,
               Matrix<complex<double>, Prop, RowSparse, Allocator>& Umat,
               bool permuted = true);

#endif

    template<class T0, class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, Prop,
			 Storage, Allocator> & mat,
			 bool keep_matrix = false);

    void FactorizeCSC(Vector<long>& Ptr, Vector<int>& IndRow,
		      Vector<complex<double> >& Val, bool sym);
    
    template<class Allocator2>
    void Solve(Vector<complex<double>, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
               Vector<complex<double>, VectFull, Allocator2>& x);

    void Solve(const SeldonTranspose& Trans,
	       complex<double>* x_ptr, int nrhs_);

    template<class Allocator2>
    void Solve(Matrix<complex<double>, General, ColMajor, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
               Matrix<complex<double>, General, ColMajor, Allocator2>& x);

#ifdef SELDON_WITH_SUPERLU_DIST
    void FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
				    Vector<int>& IndRow, Vector<complex<double> >& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    void FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
				    Vector<int64_t>& IndRow, Vector<complex<double> >& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    template<class Tint>
    void FactorizeParallel(MPI_Comm& comm_facto,
                           Vector<Tint>&, Vector<Tint>&,
                           Vector<complex<double> >&,
                           const Vector<int>& glob_num,
                           bool sym, bool keep_matrix = false);
    
    void FactorizeParallel(MPI_Comm& comm_facto,
                           Vector<superlu_int_t>&, Vector<superlu_int_t>&,
                           Vector<complex<double> >&,
                           const Vector<int>& glob_num,
                           bool sym, bool keep_matrix = false);
    
    template<class Allocator2>
    void SolveDistributed(MPI_Comm& comm_facto,
                          const SeldonTranspose& TransA,
			  Vector<complex<double>, VectFull, Allocator2>& x,
			  const Vector<int>& glob_num);

    template<class Allocator2>
    void SolveDistributed(MPI_Comm& comm_facto,
                          const SeldonTranspose& TransA,
			  Matrix<complex<double>, General, ColMajor, Allocator2>& x,
			  const Vector<int>& glob_num);

    void SolveDistributed(MPI_Comm& comm_facto,
			  const SeldonTranspose& TransA,
			  complex<double>* x_ptr, int nrhs_,
			  const IVect& glob_num);
#endif

  };

  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixSuperLU<T>& mat_lu,
	     bool keep_matrix = false);

  template<class T, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x);
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixSuperLU<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixSuperLU<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);
  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);
  
}

#define SELDON_FILE_SUPERLU_HXX
#endif
