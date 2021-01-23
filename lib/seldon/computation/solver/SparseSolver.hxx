// Copyright (C) 2010 Vivien Mallet
// Copyright (C) 2010 Marc Durufle
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


#ifndef SELDON_FILE_COMPUTATION_SPARSESOLVER_HXX

#include "Ordering.hxx"

namespace Seldon
{

  //! Base class for an interface with a direct solver
  template<class T>
  class VirtualSparseDirectSolver
  {
  public:
    virtual ~VirtualSparseDirectSolver();
    
    virtual void SetPivotThreshold(double);
    virtual void RefineSolution();
    virtual void DoNotRefineSolution();
    
    virtual void SetCoefficientEstimationNeededMemory(double coef);
    virtual void SetMaximumCoefficientEstimationNeededMemory(double coef);
    virtual void SetIncreaseCoefficientEstimationNeededMemory(double coef);

    virtual void SelectOrdering(int);
    virtual void SelectParallelOrdering(int);
    virtual void SetPermutation(const Vector<int>&);
    virtual void SetNumberOfThreadPerNode(int n);
    
#ifdef SELDON_WITH_MPI
    virtual void 
    FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
			       Vector<int>& IndRow, Vector<T>& Val,
			       const Vector<int>& glob_number,
			       bool sym, bool keep_matrix = false);

    virtual void 
    FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
			       Vector<int64_t>& IndRow, Vector<T>& Val,
			       const Vector<int>& glob_number,
			       bool sym, bool keep_matrix = false);

    virtual void 
    PerformAnalysisDistributed(MPI_Comm& comm_facto, Vector<long>& Ptr,
			       Vector<int>& IndRow, Vector<T>& Val,
			       const Vector<int>& glob_number,
			       bool sym, bool keep_matrix = false);
    
    virtual void 
    PerformAnalysisDistributed(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
			       Vector<int64_t>& IndRow, Vector<T>& Val,
			       const Vector<int>& glob_number,
			       bool sym, bool keep_matrix = false);

    virtual void 
    PerformFactorizationDistributed(MPI_Comm& comm_facto, Vector<long>& Ptr,
				    Vector<int>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);

    virtual void 
    PerformFactorizationDistributed(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
				    Vector<int64_t>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix = false);
    
    virtual void SolveDistributed(MPI_Comm& comm_facto,
				  const SeldonTranspose& TransA,
				  T* x_ptr, int nrhs,
				  const IVect& glob_num);    
#endif


    // pure virtual methods (must be overloaded)

    // (returns true if using 64-bit integers)
    virtual bool UseInteger8() const = 0; 
    virtual size_t GetMemorySize() const = 0;
    virtual int GetInfoFactorization() const = 0;
    
    virtual void Clear() = 0;
    virtual void HideMessages() = 0;
    virtual void ShowMessages() = 0;
    
    virtual void Solve(const SeldonTranspose&, T* x_ptr, int nrhs) = 0;
    
  };
  
  
  //! Default solver in Seldon
  template<class T, class Allocator 
	   = typename SeldonDefaultAllocator<ArrayRowSparse, T>::allocator>
  class SparseSeldonSolver : public VirtualSparseDirectSolver<T>
  {
  protected :
    //! Verbosity level
    int print_level;
    //! Threshold for pivoting
    double permtol;
    //! Symmetric matrix.
    Matrix<T, Symmetric, ArrayRowSymSparse, Allocator> mat_sym;
    //! Unsymmetric matrix.
    Matrix<T, General, ArrayRowSparse, Allocator> mat_unsym;
    //! Permutation arrays.
    IVect permutation_row, permutation_col;
    //! if true the factorisation is contained in mat_sym
    bool symmetric_matrix;
    
  public :
    
    SparseSeldonSolver();
    
    bool UseInteger8() const;
    void Clear();
    
    void HideMessages();
    void ShowMessages();
    int GetPrintLevel() const;
    
    size_t GetMemorySize() const;
    int GetInfoFactorization() const;
    
    double GetPivotThreshold() const;
    void SetPivotThreshold(double);
    
    template<class T0, class Storage0, class Allocator0>
    void FactorizeMatrix(const IVect& perm,
                         Matrix<T0, General, Storage0, Allocator0>& mat,
                         bool keep_matrix = false);

    template<class T0, class Storage0, class Allocator0>
    void FactorizeMatrix(const IVect& perm,
                         Matrix<T0, Symmetric, Storage0, Allocator0>& mat,
                         bool keep_matrix = false);

    template<class T1>
    void Solve(Vector<T1>& z);

    template<class T1>
    void Solve(const SeldonTranspose& TransA, Vector<T1>& z);

    void Solve(const SeldonTranspose&, T* x_ptr, int nrhs);
    
  };

  template<class T, class Allocator>
  void GetLU(Matrix<T, General, ArrayRowSparse, Allocator>& A,
	     IVect& iperm, IVect& rperm, 
	     const double& permtol, int print_level);

  template<class T1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void SolveLuVector(const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		     Vector<T2, Storage2, Allocator2>& x);

  template<class T1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void SolveLuVector(const SeldonTranspose& transA,
		     const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
		     Vector<T2, Storage2, Allocator2>& x);
  
  template<class T, class Allocator>
  void GetLU(Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A, int print_level);

  template<class T1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void SolveLuVector(const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
		     Vector<T2, Storage2, Allocator2>& x);

  template<class T1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void SolveLuVector(const SeldonTranspose& transA,
		     const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
		     Vector<T2, Storage2, Allocator2>& x);

  template<class T0, class Prop, class Storage, class Allocator,
	   class T, class Alloc2>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A,
	     SparseSeldonSolver<T, Alloc2>& mat_lu,
	     bool keep_matrix = false);

  template<class T0, class Prop, class Storage, class Allocator,
	   class T, class Alloc2>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A,
	     SparseSeldonSolver<T, Alloc2>& mat_lu,
	     IVect& permut, bool keep_matrix = false);

  template<class T, class Alloc2, class T1, class Allocator>
  void SolveLU(SparseSeldonSolver<T, Alloc2>& mat_lu,
	       Vector<T1, VectFull, Allocator>& x);

  template<class T, class Alloc2, class T1, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       SparseSeldonSolver<T, Alloc2>& mat_lu,
	       Vector<T1, VectFull, Allocator>& x);
    
}  // namespace Seldon.


#define SELDON_FILE_COMPUTATION_SPARSESOLVER_HXX
#endif
