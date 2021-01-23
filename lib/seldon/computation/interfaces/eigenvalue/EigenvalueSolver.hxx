#ifndef SELDON_FILE_EIGENVALUE_SOLVER_HXX

namespace Seldon
{
  
  //! Base class to solve an eigenvalue problem
  /*!
    Resolution of a standard eigenvalue problem : K x = lambda x
    or a generalized eigenvalue problem K x = lambda M x
    M is called mass matrix, K stiffness matrix, lambda is the eigenvalue
    and x the eigenvector.
    
    This class should not be instantiated directly, but rather derived classes
    like DenseEigenProblem, SparseEigenProblem, MatrixFreeEigenProblem.
  */
  template<class T, class MatStiff = Matrix<T>,
           class MatMass = Matrix<double> >
  class EigenProblem_Base
  {
  public :
    //! several available modes to find eigenvalues (Arpack)
    /*!
      REGULAR_MODE : Regular mode
                 standard problem => no linear system to solve
                 generalized problem => M^-1 K x = lambda x (inverse of M required)
      SHIFTED_MODE : Shifted mode 
                 standard problem => (K - sigma I)^-1 X = lambda X
                 generalized problem => (K - sigma M)^-1 M X = lambda X
      BUCKLING_MODE : Buckling mode  (real symmetric problem)
                 generalized problem => (K - sigma M)^-1 K X = lambda X
      CAYLEY_MODE : Cayley mode (real symmetric problem)
                 generalized problem => (K - sigma M)^-1 (K + sigma M) X = lambda X
      INVERT_MODE : Shifted mode on matrix M^-1 K      
      IMAG_SHIFTED_MODE : using Imag( (K - sigma M)^-1 M)
                     instead of Real( (K - sigma M)^-1 M)
                         mode 4 in Arpack (dnaupd.f)
    */
    enum {REGULAR_MODE, SHIFTED_MODE, IMAG_SHIFTED_MODE, INVERT_MODE,
          BUCKLING_MODE, CAYLEY_MODE};
    
    //! parts of the spectrum (near from 0, at infinity or around a given value)
    /*!
      SMALL_EIGENVALUES : seeking eigenvalues near 0
      LARGE_EIGENVALUES : seeking largest eigenvalues
      CENTERED_EIGENVALUES : seeking eigenvalues near the shift sigma
    */
    enum {SMALL_EIGENVALUES, LARGE_EIGENVALUES, CENTERED_EIGENVALUES};

    //! different sorting strategies
    enum {SORTED_REAL, SORTED_IMAG, SORTED_MODULUS};
    
    //! different solvers (Anasazi)
    /*!
      SOLVER_LOBPCG : Locally Optimal Block Preconditioned Conjugate Gradient
      SOLVER_BKS : Block Krylov Schur
      SOLVER_BD : Block Davidson
    */
    enum {SOLVER_LOBPCG, SOLVER_BKS, SOLVER_BD};
    
    //! orthogonalization managers (Anasazi)
    enum {ORTHO_DGKS, ORTHO_SVQB};
    
    //! type for number stored in mass matrix
    typedef typename MatMass::value_type MassValue;
        
  protected :
    //! mode used to find eigenvalues (regular, shifted, Cayley, etc)
    int eigenvalue_computation_mode;
    
    //! number of eigenvalues to be computed
    int nb_eigenvalues_wanted; 
    
    //! additional number of eigenvalues
    /*! Sometimes Arpack finds more converged eigenvalues than asked
      it is needed to store these eigenvalues and eigenvalues
      to avoid segmentation fault */
    int nb_add_eigenvalues;

    //! which spectrum ? Near from Zero ? Near from Infinity ? or near from a value ?
    int type_spectrum_wanted;
    
    //! large eigenvalues because of their real part, imaginary part or magnitude ?
    int type_sort_eigenvalues;
    
    //! if true, the generalized problem is reduced to a standard problem
    /*!
      If matrix M is symmetric definite positive, one may compute Cholesky
      factorisation of M = L L^T, and find eigenvalues of the standard problem :
      L^-1 K L^-T x = lambda x
     */    
    bool use_cholesky;
    
    //! if true, the generalized problem is reduced to a standard problem
    /*!
      if M is diagonal, one can seek eigenvalues of the 
      standard problem M^-1/2  K  M^-1/2  x = lambda x
    */
    bool diagonal_mass;
      
    //! threshold for Arpack's iterative process
    double stopping_criterion;
    
    //! Maximal number of iterations
    int nb_maximum_iterations;
    
    //! number of matrix-vector products
    int nb_prod;
    
    //! size of the problem
    int n_;
    
    //! shift sigma (if type_spectrum = centered_eigenvalues)
    T shift, shift_imag;
    
    //! number of Arnoldi vectors
    int nb_arnoldi_vectors;
    
    //! if true nb_arnoldi_vectors is automatically computed
    bool automatic_selection_arnoldi_vectors;
    
    int print_level; //!< print level
    
    //! diagonal D^1/2 if the mass matrix is diagonal positive
    Vector<MassValue> sqrt_diagonal_mass;

    //! if true consider Real( (a M + bK)^-1) or Imag( (a M + b K)^-1 )
    //! or the whole system, a and/or b being complex
    bool complex_system;    
    
    //! mass matrix
    MatMass* Mh;

    //! stiffness matrix
    MatStiff* Kh;

#ifdef SELDON_WITH_MPI
    //! communicator used to compute eigenvalues
    MPI_Intracomm comm;
#endif

    //! which solver ?
    int type_solver;
    
    /************************** 
     * Parameters for Anasazi *
     **************************/
    
    //! orthogonalization manager
    int ortho_manager;
    
    //! number of blocks for blocked solvers
    int nb_blocks;
    
    //! restart parameter for blocked solvers
    int restart_number;
    
    /************************
     * Parameters for Feast *
     ************************/
    
    //! interval where eigenvalues are searched
    double emin_interval, emax_interval;
    
  public :

    EigenProblem_Base();
    
    // initialization
    void Init(int n);
    
    void InitMatrix(MatStiff&);
    void InitMatrix(MatStiff&, MatMass& );
    
    // basic functions
    int GetComputationalMode() const;
    void SetComputationalMode(int mode);
    
#ifdef SELDON_WITH_MPI
    void SetCommunicator(MPI_Comm& comm_);
    MPI_Intracomm& GetCommunicator();
#endif
    
    int GetNbAskedEigenvalues() const;
    void SetNbAskedEigenvalues(int n);

    int GetNbAdditionalEigenvalues() const;
    void SetNbAdditionalEigenvalues(int n);
    
    int GetNbBlocks() const;
    void SetNbBlocks(int);
    int GetNbMaximumRestarts() const;
    void SetNbMaximumRestarts(int);
    int GetOrthoManager() const;

    int GetEigensolverType() const;
    void SetEigensolverType(int type);
    
    int GetTypeSpectrum() const;
    int GetTypeSorting() const;
    T GetShiftValue() const;
    T GetImagShiftValue() const;
    void SetShiftValue(const T&);
    void SetImagShiftValue(const T&);

    void SetTypeSpectrum(int type, const T& val,
                         int type_sort = SORTED_MODULUS);
    
    void SetTypeSpectrum(int type, const complex<T>& val,
                         int type_sort = SORTED_MODULUS);

    double GetLowerBoundInterval() const;
    double GetUpperBoundInterval() const;

    void SetIntervalSpectrum(double, double);
            
    void SetCholeskyFactoForMass(bool chol = true);
    bool UseCholeskyFactoForMass() const;
    
    void SetDiagonalMass(bool diag = true);
    bool DiagonalMass() const;
    
    void SetStoppingCriterion(double eps);
    double GetStoppingCriterion() const;
    void SetNbMaximumIterations(int n);
    int GetNbMaximumIterations() const;
    
    int GetNbMatrixVectorProducts() const;
    
    int GetNbArnoldiVectors() const;
    void SetNbArnoldiVectors(int n);
    
    int GetM() const;
    int GetN() const;
    
    int GetPrintLevel() const;
    void SetPrintLevel(int lvl);
    
    void IncrementProdMatVect();
    
    void PrintErrorInit() const;
    bool IsSymmetricProblem() const;
    bool IsHermitianProblem() const;
    
    // mass matrix stuff
    void FactorizeDiagonalMass();
    template<class T0> void MltInvSqrtDiagonalMass(Vector<T0>& X);
    template<class T0> void MltSqrtDiagonalMass(Vector<T0>& X);
    
    void ComputeDiagonalMass();
    void ComputeMassForCholesky();
    
    void ComputeMassMatrix();
    void MltMass(const Vector<T>& X, Vector<T>& Y);
    
    // stiffness matrix stuff
    void ComputeStiffnessMatrix();    
    void ComputeStiffnessMatrix(const T& a, const T& b);
    
    void MltStiffness(const Vector<T>& X, Vector<T>& Y);
    void MltStiffness(const T& a, const T& b, const Vector<T>& X, Vector<T>& Y);
        
    // functions to overload (factorisation of mass and/or stiffness matrix)
    void ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b);
    void ComputeAndFactorizeStiffnessMatrix(const complex<T>& a,
                                            const complex<T>& b,
                                            bool real_p = true);
    
    template<class T0>
    void ComputeSolution(const Vector<T0>& X, Vector<T0>& Y);
    
    template<class TransA, class T0>
    void ComputeSolution(const TransA& transA,
                         const Vector<T0>& X, Vector<T0>& Y);
    
    void FactorizeCholeskyMass();
    
    template<class TransStatus>
    void MltCholeskyMass(const TransStatus& transA, Vector<T>& X);
    
    template<class TransStatus>
    void SolveCholeskyMass(const TransStatus& transA, Vector<T>& X);
    
    void Clear();
    
    // modification of eigenvectors to take into account 
    // the use of matrix D^-1/2 K D^-1/2 or L^-1 K L^-T instead of K
    // => eigenvectors are recovered by multiplying them by matrix D^1/2 or by L^T
    template<class EigenPb, class Vector1, class Matrix1, class T0>
    friend void ApplyScalingEigenvec(EigenPb& var,
                                     Vector1& eigen_values,
                                     Vector1& lambda_imag,
                                     Matrix1& eigen_vectors,
                                     const T0& shiftr, const T0& shifti);

    template<class EigenPb, class Vector1, class Matrix1, class T0>
    friend void ApplyScalingEigenvec(EigenPb& var,
                                     Vector1& eigen_values,
                                     Vector1& lambda_imag,
                                     Matrix1& eigen_vectors,
                                     const complex<T0>& shiftr,
                                     const complex<T0>& shifti);
    
  };    
  
  
  // sorting eigenvalues
  template<class T, class Storage1, class Storage2,
           class Alloc1, class Alloc2>
  void SortEigenvalues(Vector<T>& lambda_r, Vector<T>& lambda_i,
                       Matrix<T, General, Storage1, Alloc1>& eigen_old,
                       Matrix<T, General, Storage2, Alloc2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const T& shift_r, const T& shift_i);
  
  // sorting eigenvalues
  template<class T, class Storage1, class Storage2,
           class Alloc1, class Alloc2>
  void SortEigenvalues(Vector<complex<T> >& lambda_r,
                       Vector<complex<T> >& lambda_i,
                       Matrix<complex<T>, General, Storage1, Alloc1>& eigen_old,
                       Matrix<complex<T>, General, Storage2, Alloc2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const complex<T>& shift_r, const complex<T>& shift_i);
  
  
  //! computation of a few eigenvalues for dense matrices
  template<class T, class Prop, class Storage,
           class Tmass = double, class PropM = Symmetric,
           class StorageM = RowSymPacked>
  class DenseEigenProblem
    : public EigenProblem_Base<T, Matrix<T, Prop, Storage>,
                               Matrix<Tmass, PropM, StorageM> >
  {
  protected :
    //! LU factorisation of matrix
    Matrix<T, Prop, Storage> mat_lu;
    
    //! pivot used by the LU factorisation
    Vector<int> pivot;

    //! Cholesky factorisation of mass matrix
    Matrix<Tmass, PropM, StorageM> mat_chol;
    
    //! temporary vectors for Cholesky
    Vector<Tmass> Xchol_real, Xchol_imag;
    
  public :

    DenseEigenProblem();

    void FactorizeCholeskyMass();
    
    template<class TransStatus, class T0>
    void MltCholeskyMass(const TransStatus& transA, Vector<T0>& X);
    
    template<class TransStatus, class T0>
    void SolveCholeskyMass(const TransStatus& transA, Vector<T0>& X);

    template<class TransStatus>
    void MltCholeskyMass(const TransStatus& transA,
                         Vector<complex<double> >& X);
    
    template<class TransStatus>
    void SolveCholeskyMass(const TransStatus& transA,
                           Vector<complex<double> >& X);
    
    void ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b);
    void ComputeAndFactorizeStiffnessMatrix(const complex<T>& a,
                                            const complex<T>& b,
                                            bool real_p = true);
    
    template<class T0>
    void ComputeAndFactorizeComplexMatrix(const complex<T0>& a,
                                          const complex<T0>& b,
                                          bool real_p = true);

    void ComputeAndFactorizeComplexMatrix(const complex<double>& a,
                                          const complex<double>& b,
                                            bool real_p = true);

    template<class T0>
    void ComputeSolution(const Vector<T0>& X, Vector<T0>& Y);

    template<class TransA, class T0>
    void ComputeSolution(const TransA& transA,
                         const Vector<T0>& X, Vector<T0>& Y);
    
    void Clear();
    
  };
  
	
  //! computation of a few eigenvalues for sparse matrices
  template<class T, class MatStiff,
           class MatMass = Matrix<double, Symmetric, ArrayRowSymSparse> >
  class SparseEigenProblem
    : public EigenProblem_Base<T, MatStiff, MatMass>
  {
  public :
    typedef typename MatMass::value_type MassValue;
    
  protected :
    typedef typename ClassComplexType<T>::Tcplx Complexe;
    
    //! LU factorisation of sparse matrix
    SparseDirectSolver<T> mat_lu;
    
    //! factorisation of complex system
    SparseDirectSolver<Complexe> mat_lu_cplx;
    
    //! Cholesky factorisation of mass matrix if required
    SparseCholeskySolver<double> chol_facto_mass_matrix;
    
    //! temporary vectors for Cholesky
    Vector<double> Xchol_real, Xchol_imag;
    
    //! if true, we take imaginary part of (K - sigma M)^-1
    bool imag_solution;
    
  public :
    
    SparseEigenProblem();
    
    void SelectCholeskySolver(int type);
    void FactorizeCholeskyMass();
    
    template<class Storage, class Alloc>
    void FactorizeCholeskyMass(double&, Matrix<double, Symmetric, Storage, Alloc>& M);
    
    template<class T0, class T1, class Prop, class Storage, class Alloc>
    void FactorizeCholeskyMass(T0&, Matrix<T1, Prop, Storage, Alloc>& M);
    
    template<class TransStatus, class T0>
    void MltCholeskyMass(const TransStatus& transA, Vector<T0>& X);
    
    template<class TransStatus, class T0>
    void SolveCholeskyMass(const TransStatus& transA, Vector<T0>& X);

    template<class TransStatus>
    void MltCholeskyMass(const TransStatus& transA,
                         Vector<complex<double> >& X);

    template<class TransStatus>
    void MltCholeskyMass(const TransStatus& transA,
                         Vector<complex<double> >& X, double&);

    template<class TransStatus>
    void MltCholeskyMass(const TransStatus& transA,
                         Vector<complex<double> >& X, complex<double>&);

    template<class TransStatus>
    void SolveCholeskyMass(const TransStatus& transA,
			   Vector<complex<double> >& X);
    
    template<class TransStatus>
    void SolveCholeskyMass(const TransStatus& transA,
                           Vector<complex<double> >& X, double&);

    template<class TransStatus>
    void SolveCholeskyMass(const TransStatus& transA,
                           Vector<complex<double> >& X, complex<double>&);
    
    void ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b);
    void ComputeAndFactorizeStiffnessMatrix(const complex<T>& a,
                                            const complex<T>& b,
                                            bool real_p = true);
    
    template<class T0>
    void ComputeAndFactorizeComplexMatrix(const complex<T0>& a,
                                          const complex<T0>& b,
                                          bool real_p = true);

    void ComputeAndFactorizeComplexMatrix(const complex<double>& a,
                                          const complex<double>& b,
                                            bool real_p = true);

    template<class T0>
    void ComputeSolution(const Vector<T0>& X, Vector<T0>& Y);
    
    template<class TransA, class T0>
    void ComputeSolution(const TransA& transA,
                         const Vector<T0>& X, Vector<T0>& Y);
    
    template<class TransA>
    void ComputeComplexSolution(const TransA&,
                                const Vector<double>& X, Vector<double>& Y);
    
    template<class TransA>
    void ComputeComplexSolution(const TransA&,
                                const Vector<complex<double> >& X,
                                Vector<complex<double> >& Y);
    
    void Clear();
    
  };
 
  
  //! computation of a few eigenvalues with only matrix-vector product
  template<class T, class MatStiff>
  class MatrixFreeEigenProblem
    : public EigenProblem_Base<T, MatStiff,
                               Matrix<double, Symmetric, ArrayRowSymSparse> >
  {
  public :

    MatrixFreeEigenProblem();
    
  };
  
  
  //! list of availables eigenvalue solvers
  class TypeEigenvalueSolver
  {
  public :
    enum {DEFAULT, ARPACK, ANASAZI, FEAST};
    
    static int default_solver;
    
    static int GetDefaultSolver();
    
  };
  
  template<class EigenPb, class Vector1, class Vector2,
            class T, class Prop, class Storage, class Alloc3>
  void GetEigenvaluesEigenvectors(EigenPb& var_eig, Vector1& lambda,
				  Vector2& lambda_imag,
				  Matrix<T, Prop, Storage, Alloc3>& eigen_vec);
  
}

#define SELDON_FILE_EIGENVALUE_SOLVER_HXX
#endif

