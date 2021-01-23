#ifndef SELDON_FILE_DISTRIBUTED_SOLVER_HXX

namespace Seldon
{

  //! general class for direct solver
  /*!
    This class will call one of the direct solver interfaced
    in Seldon (in sequential for Matrix and in distributed
    for DistributedMatrix). The user can also provide
    left/right scaling to apply to the matrix before
    factorization.
   */
  template<class T>
  class SparseDistributedSolver : public SparseDirectSolver<T>
  {
  protected :
    typedef typename ClassComplexType<T>::Treal Treal;
    
    bool diagonal_scaling_left; //!< left scaling ?
    bool diagonal_scaling_right; //!< right scaling ?
    Vector<Treal> diagonal_scale_left; //!< left scaling
    Vector<Treal> diagonal_scale_right; //!< right scaling

#ifdef SELDON_WITH_MPI
    // data associated with distributed matrix
    // see DistributedMatrix.hxx for a detailed description
    int nodl_scalar_, nb_unknowns_scal_;
    MPI_Comm comm_;
    IVect* ProcSharingRows_;
    Vector<IVect>* SharingRowNumbers_;
    IVect global_col_numbers, local_col_numbers;
    
    template<class T2>
    void AssembleVec(Vector<T2>& X) const;

    template<class T2>
    void AssembleVec(Matrix<T2, General, ColMajor>& A) const;
#endif

    template<class MatrixSparse>
    void ScaleMatrixRowCol(MatrixSparse& A);

  public :
    
    SparseDistributedSolver();
    
    void SetPrintLevel(int print);
    
    void Clear();

    template<class Prop0, class Storage0, class Allocator0>
    void Factorize(Matrix<T, Prop0, Storage0, Allocator0>& A,
		   bool keep_matrix = false, bool scale_matrix = false);

    template<class Prop0, class Storage0, class Allocator0>
    void PerformAnalysis(Matrix<T, Prop0, Storage0, Allocator0>& A);

    template<class Prop0, class Storage0, class Allocator0>
    void PerformFactorization(Matrix<T, Prop0, Storage0, Allocator0>& A,
			      bool scale_matrix = false);

#ifdef SELDON_WITH_MPI
    template<class Prop0, class Storage0, class Allocator0>
    void Factorize(DistributedMatrix<T, Prop0, Storage0, Allocator0>& A,
                   bool keep_matrix = false, bool scale_matrix = false);

    template<class Prop0, class Storage0, class Allocator0>
    void PerformAnalysis(DistributedMatrix<T, Prop0, Storage0, Allocator0>& A);
        
    template<class Prop0, class Storage0, class Allocator0>
    void PerformFactorization(DistributedMatrix<T, Prop0, Storage0, Allocator0>& A,
			      bool scale_matrix = false);
#endif
    
    template<class T1>
    void Solve(Vector<T1>& x_solution, const Vector<T1>& b_rhs);    

    template<class T1>
    void Solve(Vector<T1>& x_solution);

    template<class T1>
    void TransSolve(Vector<T1>& x_solution);

    template<class T1>
    void Solve(const SeldonTranspose&, Vector<T1>& x_solution, bool assemble = true);
    
    template<class T1>
    void Solve(Matrix<T1, General, ColMajor>& x_solution);    

    template<class T1>
    void TransSolve(Matrix<T1, General, ColMajor>& x_solution);    

    template<class T1>
    void Solve(const SeldonTranspose&,
               Matrix<T1, General, ColMajor>& x_solution);    
    
    template<class MatrixSparse, class MatrixFull>
    void GetSchurComplement(MatrixSparse& mat_direct,
			    const IVect& num, MatrixFull& mat_schur);
    
    size_t GetMemorySize() const;
        
  };
  
}

#define SELDON_FILE_DISTRIBUTED_SOLVER_HXX
#endif
  
