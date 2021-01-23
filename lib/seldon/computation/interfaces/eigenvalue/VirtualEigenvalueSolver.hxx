#ifndef SELDON_FILE_VIRTUAL_EIGENVALUE_SOLVER_HXX

#ifdef SELDON_WITH_SLEPC    
#include <petscconf.h>

#if defined(PETSC_USE_COMPLEX)
#define Petsc_Scalar complex<double>
#else
#define Petsc_Scalar double
#endif

#define Petsc_Error_Code int
#define Petsc_Int int

#endif

namespace Seldon
{

  //! Parameters for Anasazi package
  class AnasaziParam
  {
  public:
    //! different solvers
    /*!
      SOLVER_LOBPCG : Locally Optimal Block Preconditioned Conjugate Gradient
      SOLVER_BKS : Block Krylov Schur
      SOLVER_BD : Block Davidson
    */
    enum {SOLVER_LOBPCG, SOLVER_BKS, SOLVER_BD};
    
    //! orthogonalization managers
    enum {ORTHO_DGKS, ORTHO_SVQB};

  protected:
    //! orthogonalization manager
    int ortho_manager;
    
    //! number of blocks for blocked solvers
    int nb_blocks;
    
    //! restart parameter for blocked solvers
    int restart_number;

    //! which solver ?
    int type_solver;
        
  public:
    AnasaziParam();
    
    int GetNbBlocks() const;
    void SetNbBlocks(int);
    int GetNbMaximumRestarts() const;
    void SetNbMaximumRestarts(int);
    int GetOrthoManager() const;

    int GetEigensolverType() const;
    void SetEigensolverType(int type);

  };


  //! Parameters for Slepc package
  class SlepcParam
  {
  protected:
    //! which solver ?
    int type_solver;
    int block_size;
    int nstep;
    int type_extraction;
    int quadrature_rule;
    int nstep_inner, nstep_outer;
    int npoints;
    int moment_size, block_max_size, npart;
    double delta_rank, delta_spur;
    int borth;
    int double_exp;
    int init_size;
    int krylov_restart, restart_number, restart_add;
    int nb_conv_vector, nb_conv_vector_proj;
    bool non_locking_variant;
    double restart_ratio;
    string method;
    int shift_type;
    
  public:
    enum {SHIFT_CONSTANT, SHIFT_RAYLEIGH, SHIFT_WILKINSON};
    
    enum {QUADRULE_TRAPEZE, QUADRULE_CHEBY};
    
    enum {EXTRACT_RITZ, EXTRACT_HANKEL};
    
    enum {POWER, SUBSPACE, ARNOLDI, LANCZOS, KRYLOVSCHUR, GD, JD,
	  RQCG, LOBPCG, CISS, LAPACK, ARPACK, BLZPACK, TRLAN, BLOPEX, PRIMME, FEAST};

    SlepcParam();
    
    const char* GetEigensolverChar() const;
    int GetEigensolverType() const;
    void SetEigensolverType(int type);

    void SetBlockSize(int n);
    int GetBlockSize() const;
    void SetMaximumBlockSize(int n);
    int GetMaximumBlockSize() const;

    void SetNumberOfSteps(int n);
    int GetNumberOfSteps() const;

    void SetExtractionType(int);
    int GetExtractionType() const;

    void SetQuadratureRuleType(int);
    int GetQuadratureRuleType() const;

    void SetInnerSteps(int);
    void SetOuterSteps(int);
    int GetInnerSteps() const;
    int GetOuterSteps() const;

    int GetNumberIntegrationPoints() const;
    void SetNumberIntegrationPoints(int);
    int GetMomentSize() const;
    void SetMomentSize(int);
    int GetNumberPartitions() const;
    void SetNumberPartitions(int);

    void SetThresholdRank(double);
    double GetThresholdRank() const;
    void SetThresholdSpurious(double);
    double GetThresholdSpurious() const;

    int GetBorthogonalization() const;
    void SetBorthogonalization(int);
    int GetDoubleExpansion() const;
    void SetDoubleExpansion(int);
    int GetInitialSize() const;
    void SetInitialSize(int);
    int GetKrylovRestart() const;
    void SetKrylovRestart(int);
    int GetRestartNumber() const;
    void SetRestartNumber(int);
    int GetRestartNumberAdd() const;
    void SetRestartNumberAdd(int);

    int GetNumberConvergedVectors() const;
    void SetNumberConvergedVectors(int);
    int GetNumberConvergedVectorsProjected() const;
    void SetNumberConvergedVectorsProjected(int);

    bool UseNonLockingVariant() const;
    void SetNonLockingVariant(bool n);
    double GetRestartRatio() const;
    void SetRestartRatio(double d);
    string GetMethod() const;
    void SetMethod(string s);
    int GetShiftType() const;
    void SetShiftType(int n);
    
  };


  //! parameters for Feast package
  class FeastParam
  {
  protected:
    bool evaluate_number_eigenval;
    int nb_points_quadrature;
    int type_integration;
    
  public:
    FeastParam();

    void EnableEstimateNumberEigenval(bool p = true);
    bool EstimateNumberEigenval() const;

    int GetNumOfQuadraturePoints() const;
    void SetNumOfQuadraturePoints(int);
    
    int GetTypeIntegration() const;
    void SetTypeIntegration(int);
    
  };

  template<class T0, class T1>
  void to_complex_eigen(const T0& x, T1& y);
  
  template<class T0, class T1>
  void to_complex_eigen(const complex<T0>& x, T1& y);

  template<class T0, class T1>
  void to_complex_eigen(const T0& x, complex<T1>& y);

  template<class T0, class T1>
  void to_complex_eigen(const complex<T0>& x, complex<T1>& y);


  //! base class for setting comparison between eigenvalues
  template<class T>
  class EigenvalueComparisonClass  
  {
  public:
    virtual int CompareEigenvalue(const T& Lr, const T& Li, const T& Lr2, const T& Li2) = 0;
    
  };
  
  //! Base class to solve an eigenvalue problem
  /*!
    Resolution of a standard eigenvalue problem : K x = lambda x
    or a generalized eigenvalue problem K x = lambda M x
    M is called mass matrix, K stiffness matrix, lambda is the eigenvalue
    and x the eigenvector.
    
    This class should not be instantiated directly, but rather derived classes
    like DenseEigenProblem, SparseEigenProblem, VirtualEigenProblem.
  */
  template<class T>
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
    enum {SORTED_REAL, SORTED_IMAG, SORTED_MODULUS, SORTED_USER};

    //! real part, imaginary part or complex solution
    enum { REAL_PART, IMAG_PART, COMPLEX_PART};
            
  protected :
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;
    
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
    Treal stopping_criterion;
    
    //! Maximal number of iterations
    int nb_maximum_iterations;
    
    //! number of matrix-vector products
    int nb_prod;
    
    //! size of the problem
    int n_, nglob;
    
    //! shift sigma (if type_spectrum = centered_eigenvalues)
    T shift, shift_imag; 

    //! number of Arnoldi vectors
    int nb_arnoldi_vectors;
    
    //! if true nb_arnoldi_vectors is automatically computed
    bool automatic_selection_arnoldi_vectors;
    
    int print_level; //!< print level
    int global_print_level; //!< global print level (for Feast)
    
    //! if true consider Real( (a M + bK)^-1) or Imag( (a M + b K)^-1 )
    //! or the whole system, a and/or b being complex
    bool complex_system; int selected_part;    
    
#ifdef SELDON_WITH_MPI
    //! communicator used to compute eigenvalues (linear system is split into processors)
    MPI_Comm comm;

    //! communicator for feast (each processor has a linear system)
    MPI_Comm comm_global;
#endif

    //! interval where eigenvalues are searched (real symmetric or hermitian)
    Treal emin_interval, emax_interval;

    //! disk where eigenvalues are searched (non-symmetric)
    Tcplx center_spectrum; Treal radius_spectrum;

    //! parameters for an ellipse
    Treal ratio_ellipse, angle_ellipse;
    
    //! additional parameters for Anasazi
    AnasaziParam anasazi_param;

    //! additional parameters for Slepc
    SlepcParam slepc_param;
    
    //! additional parameters for Feast
    FeastParam feast_param;
    
    //! class for comparing eigenvalues
    EigenvalueComparisonClass<T>* compar_eigenval;
    
  public :

    EigenProblem_Base();
    virtual ~EigenProblem_Base();
    
    // initialization
    void Init(int n);
    
    // basic functions
    int GetComputationalMode() const;
    void SetComputationalMode(int mode);
    
#ifdef SELDON_WITH_MPI
    void SetCommunicator(const MPI_Comm& comm_);
    MPI_Comm& GetCommunicator();

    void SetGlobalCommunicator(const MPI_Comm& comm_);
    MPI_Comm& GetGlobalCommunicator();
#endif

    int GetRankCommunicator() const;
    int GetGlobalRankCommunicator() const;
    
    int GetNbAskedEigenvalues() const;
    void SetNbAskedEigenvalues(int n);

    int GetNbAdditionalEigenvalues() const;
    void SetNbAdditionalEigenvalues(int n);

    AnasaziParam& GetAnasaziParameters();
    SlepcParam& GetSlepcParameters();
    FeastParam& GetFeastParameters();
    
    int GetTypeSpectrum() const;
    int GetTypeSorting() const;
    T GetShiftValue() const;
    T GetImagShiftValue() const;
    void SetShiftValue(const T&);
    void SetImagShiftValue(const T&);
    
    Treal GetCircleRadiusSpectrum() const;
    Tcplx GetCircleCenterSpectrum() const;
    Treal GetRatioEllipseSpectrum() const;
    Treal GetAngleEllipseSpectrum() const;
    
    void GetComplexShift(const Treal&, const Treal&, Tcplx&) const;
    void GetComplexShift(const Tcplx&, const Tcplx&, Tcplx&) const;

    void SetTypeSpectrum(int type, const T& val,
                         int type_sort = SORTED_MODULUS);
    
    void SetTypeSpectrum(int type, const complex<T>& val,
                         int type_sort = SORTED_MODULUS);
    
    void SetUserComparisonClass(EigenvalueComparisonClass<T>* ev);
    
    Treal GetLowerBoundInterval() const;
    Treal GetUpperBoundInterval() const;

    void SetIntervalSpectrum(Treal, Treal);
    void SetCircleSpectrum(const Tcplx& z, const Treal& r);
    void SetEllipseSpectrum(const Tcplx& z, const Treal& r,
                            const Treal& ratio, const Treal& teta);
    
    void SetCholeskyFactoForMass(bool chol = true);
    bool UseCholeskyFactoForMass() const;
    
    void SetDiagonalMass(bool diag = true);
    bool DiagonalMass() const;
    
    void SetStoppingCriterion(Treal eps);
    Treal GetStoppingCriterion() const;
    void SetNbMaximumIterations(int n);
    int GetNbMaximumIterations() const;
    
    int GetNbMatrixVectorProducts() const;
    
    int GetNbArnoldiVectors() const;
    void SetNbArnoldiVectors(int n);
    
    int GetM() const;
    int GetGlobalM() const;
    int GetN() const;
    
    int GetPrintLevel() const;
    int GetGlobalPrintLevel() const;
    void SetPrintLevel(int lvl);
    void SetGlobalPrintLevel(int lvl);
    
    void IncrementProdMatVect();
    
    void PrintErrorInit() const;
    
    virtual bool IsSymmetricProblem() const = 0;
    virtual bool IsHermitianProblem() const = 0;
    
    // mass matrix stuff
    virtual void ComputeDiagonalMass() = 0;
    virtual void FactorizeDiagonalMass() = 0;
    virtual void GetSqrtDiagonal(Vector<T>&) = 0;
    
    virtual void MltInvSqrtDiagonalMass(Vector<Treal>& X) = 0;
    virtual void MltSqrtDiagonalMass(Vector<Treal>& X) = 0;

    virtual void MltInvSqrtDiagonalMass(Vector<Tcplx>& X) = 0;
    virtual void MltSqrtDiagonalMass(Vector<Tcplx>& X) = 0;
    
    virtual void ComputeMassForCholesky();    
    virtual void ComputeMassMatrix();

    virtual void MltMass(const Vector<Treal>& X, Vector<Treal>& Y) = 0;
    virtual void MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y) = 0;
    virtual void MltMass(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y) = 0;
    virtual void MltMass(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y) = 0;
    
    // stiffness matrix stuff
    virtual void ComputeStiffnessMatrix();    
    virtual void ComputeStiffnessMatrix(const T& a, const T& b);
    
    virtual void MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y) = 0;
    virtual void MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y) = 0;

    virtual void MltStiffness(const T& a, const T& b,
			      const Vector<Treal>& X, Vector<Treal>& Y) = 0;

    virtual void MltStiffness(const T& a, const T& b,
			      const Vector<Tcplx>& X, Vector<Tcplx>& Y) = 0;
    
    virtual void MltStiffness(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y) = 0;
    virtual void MltStiffness(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y) = 0;
    
    // functions to overload (factorisation of mass and/or stiffness matrix)
    virtual void ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b,
						    int real_p = COMPLEX_PART);
    
    virtual void ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b,
						    int real_p = COMPLEX_PART);
    
    virtual void ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y);
    virtual void ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    virtual void ComputeSolution(const SeldonTranspose&,
				 const Vector<Treal>& X, Vector<Treal>& Y);

    virtual void ComputeSolution(const SeldonTranspose&,
				 const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    virtual void FactorizeCholeskyMass();
    
    virtual void MltCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    virtual void MltCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);
    
    virtual void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    virtual void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);

    virtual void DistributeEigenvectors(Matrix<T, General, ColMajor>& eigen_vec);
    
    virtual void Clear();

#ifdef SELDON_WITH_SLEPC    
    template<class T0>
    static Petsc_Error_Code
    GetComparisonEigenvalueSlepcGen(Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, T0,
                                    Petsc_Int*, void*);

    static Petsc_Error_Code
    GetComparisonEigenvalueSlepcGen(Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, Petsc_Scalar,
                                    Petsc_Int*, void*);
    
    static Petsc_Error_Code
    GetComparisonEigenvalueSlepc(Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, Petsc_Scalar,
                                 Petsc_Int*, void*);
#endif
    
    void FillComplexEigenvectors(int m, const Tcplx& Emid, Treal eps,
                                 const Vector<Tcplx>& lambda_cplx,
                                 const Matrix<Tcplx, General, ColMajor>& Ecplx,
                                 Vector<T>& Lr, Vector<T>& Li,
                                 Matrix<T, General, ColMajor>& E);

    // modification of eigenvectors to take into account 
    // the use of matrix D^-1/2 K D^-1/2 or L^-1 K L^-T instead of K
    // => eigenvectors are recovered by multiplying them by matrix D^1/2 or by L^T
    template<class T0, class Prop, class Storage>
    friend void ApplyScalingEigenvec(EigenProblem_Base<T0>& var,
                                     Vector<T0>& eigen_values,
                                     Vector<T0>& lambda_imag,
                                     Matrix<T0, Prop, Storage>& eigen_vectors,
                                     const T0& shiftr, const T0& shifti);

    template<class T0, class Prop, class Storage>
    friend void ApplyScalingEigenvec(EigenProblem_Base<complex<T0> >& var,
                                     Vector<complex<T0> >& eigen_values,
                                     Vector<complex<T0> >& lambda_imag,
                                     Matrix<complex<T0>, Prop, Storage>& eigen_vectors,
                                     const complex<T0>& shiftr,
                                     const complex<T0>& shifti);
        
  };    

  
  //! base class for eigenvalue problems
  template<class T, class StiffValue = T,
	   class MassValue = typename ClassComplexType<T>::Treal>  
  class VirtualEigenProblem : public EigenProblem_Base<T>
  {
  protected:
    //! mass matrix
    VirtualMatrix<MassValue>* Mh;

    //! stiffness matrix
    VirtualMatrix<StiffValue>* Kh;

    //! diagonal D^1/2 if the mass matrix is diagonal positive
    Vector<MassValue> sqrt_diagonal_mass;

    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;

  public:    
    VirtualEigenProblem();
    
    void InitMatrix(VirtualMatrix<StiffValue>& K, int n = -1);
    void InitMatrix(VirtualMatrix<StiffValue>& K, VirtualMatrix<MassValue>& M, int n = -1);

    void SetMatrix(VirtualMatrix<StiffValue>& K, VirtualMatrix<MassValue>& M);
    
    void SetTypeSpectrum(int type, const T& val,
                         int type_sort = EigenProblem_Base<T>::SORTED_MODULUS);

    void SetTypeSpectrum(int type, const complex<T>& val,
                         int type_sort = EigenProblem_Base<T>::SORTED_MODULUS);

    bool IsSymmetricProblem() const;
    bool IsHermitianProblem() const;

    void ComputeDiagonalMass();
    void FactorizeDiagonalMass();
    void GetSqrtDiagonal(Vector<T>&);
    
    void MltInvSqrtDiagonalMass(Vector<Treal>& X);
    void MltSqrtDiagonalMass(Vector<Treal>& X);

    void MltInvSqrtDiagonalMass(Vector<Tcplx>& X);
    void MltSqrtDiagonalMass(Vector<Tcplx>& X);

    void MltMass(const Vector<Treal>& X, Vector<Treal>& Y);    
    void MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y);    
    void MltMass(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y);    
    void MltMass(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y);    
    
    void MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y);
    void MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void MltStiffness(const T& coef_mass, const T& coef_stiff,
		      const Vector<Treal>& X, Vector<Treal>& Y);

    void MltStiffness(const T& coef_mass, const T& coef_stiff,
		      const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    void MltStiffness(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y);
    void MltStiffness(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    void Clear();
    
  };
  
  
  // sorting eigenvalues
  template<class T, class Storage1, class Storage2>
  void SortEigenvalues(Vector<T>& lambda_r, Vector<T>& lambda_i,
                       Matrix<T, General, Storage1>& eigen_old,
                       Matrix<T, General, Storage2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const T& shift_r, const T& shift_i);
  
  // sorting eigenvalues
  template<class T, class Storage1, class Storage2>
  void SortEigenvalues(Vector<complex<T> >& lambda_r,
                       Vector<complex<T> >& lambda_i,
                       Matrix<complex<T>, General, Storage1>& eigen_old,
                       Matrix<complex<T>, General, Storage2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const complex<T>& shift_r, const complex<T>& shift_i);
  
  
  //! computation of a few eigenvalues for dense matrices
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass = typename ClassComplexType<T>::Treal,
	   class PropM = Symmetric, class StorageM = RowSymPacked>
  class DenseEigenProblem : public VirtualEigenProblem<T, Tstiff, Tmass>
  {
  protected :
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;

    //! LU factorisation of a real matrix
    Matrix<Treal, Prop, Storage> mat_lu_real;

    //! LU factorisation of a complex matrix
    Matrix<Tcplx, Prop, Storage> mat_lu_cplx;

    //! mass matrix
    Matrix<Tmass, PropM, StorageM>* Mh;

    //! stiffness matrix
    Matrix<Tstiff, Prop, Storage>* Kh;
    
    //! pivot used by the LU factorisation
    Vector<int> pivot;
    
    //! Cholesky factorisation of mass matrix
    Matrix<Treal, PropM, StorageM> mat_chol;
    
    //! temporary vectors for Cholesky
    Vector<Treal> Xchol_real, Xchol_imag;
    
  public :

    DenseEigenProblem();

    void InitMatrix(Matrix<Tstiff, Prop, Storage>&);
    void InitMatrix(Matrix<Tstiff, Prop, Storage>&, Matrix<Tmass, PropM, StorageM>& );

    void ComputeDiagonalMass();
    void FactorizeCholeskyMass();
    
    void MltCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    void MltCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);
    
    void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);

    void MltMass(const Vector<Treal>& X, Vector<Treal>& Y);    
    void MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y);    
    void MltMass(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y);    
    void MltMass(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y);    
    
    void MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y);
    void MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    void MltStiffness(const T& coef_mass, const T& coef_stiff,
		      const Vector<Treal>& X, Vector<Treal>& Y);

    void MltStiffness(const T& coef_mass, const T& coef_stiff,
		      const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void MltStiffness(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y);
    void MltStiffness(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void ComputeAndFactoRealMatrix(const Treal&, const Treal& a,
				   const Treal& b, int which);
    
    void ComputeAndFactoRealMatrix(const Tcplx&, const Treal& a,
				   const Treal& b, int which);

    void ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b,
					    int which_part =
					    EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b,
                                            int which_part =
					    EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y);
    void ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void ComputeSolution(const SeldonTranspose& transA,
                         const Vector<Treal>& X, Vector<Treal>& Y);

    void ComputeSolution(const SeldonTranspose& transA,
                         const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    void Clear();
    
  };
  
	
  //! computation of a few eigenvalues for sparse matrices
  template<class T, class MatStiff,
           class MatMass = Matrix<typename ClassComplexType<T>::Treal,
				  Symmetric, ArrayRowSymSparse> >
  class SparseEigenProblem
    : public VirtualEigenProblem<T, typename MatStiff::entry_type,
				 typename MatMass::entry_type>
  {
  protected :
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;

#ifdef SELDON_WITH_MPI
    //! LU factorisation of sparse matrix
    SparseDistributedSolver<Treal> mat_lu_real;
    
    //! factorisation of complex system
    SparseDistributedSolver<Tcplx> mat_lu_cplx;

    //! Cholesky factorisation of mass matrix if required
    DistributedCholeskySolver<Treal> chol_facto_mass_matrix;
#else
    //! LU factorisation of sparse matrix
    SparseDirectSolver<Treal> mat_lu_real;
    
    //! factorisation of complex system
    SparseDirectSolver<Tcplx> mat_lu_cplx;

    //! Cholesky factorisation of mass matrix if required
    SparseCholeskySolver<Treal> chol_facto_mass_matrix;
#endif
    
    //! temporary vectors for Cholesky
    Vector<Treal> Xchol_real, Xchol_imag;

    MatMass* Mh;
    MatStiff* Kh;

    bool distributed; int nloc;
    Vector<int> local_col_numbers;
    Vector<int>* ProcSharingRows;
    Vector<Vector<int> >* SharingRowNumbers;
    int nodl_scalar_, nb_unknowns_scal_;
    
    int RetrieveLocalNumbers(MatStiff& K);
    
  public :
    
    SparseEigenProblem();
    
    void SelectCholeskySolver(int type);
    
    void InitMatrix(MatStiff&);
    void InitMatrix(MatStiff&, MatMass&);

    void ComputeDiagonalMass();
    void FactorizeCholeskyMass();
    
    template<class Storage>
    void FactorizeCholeskyMass(Matrix<Treal, Symmetric, Storage>& M);
    
    template<class T0, class Prop, class Storage>
    void FactorizeCholeskyMass(Matrix<T0, Prop, Storage>& M);

#ifdef SELDON_WITH_MPI
    template<class Storage>
    void FactorizeCholeskyMass(DistributedMatrix<Treal, Symmetric, Storage>& M);
    
    template<class T0, class Prop, class Storage>
    void FactorizeCholeskyMass(DistributedMatrix<T0, Prop, Storage>& M);
#endif

    void MltCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    void MltCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);
    
    void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);

    void MltMass(const Vector<Treal>& X, Vector<Treal>& Y);    
    void MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y);    
    void MltMass(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y);    
    void MltMass(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y);    
    
    void MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y);
    void MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    void MltStiffness(const T& coef_mass, const T& coef_stiff,
		      const Vector<Treal>& X, Vector<Treal>& Y);

    void MltStiffness(const T& coef_mass, const T& coef_stiff,
		      const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void MltStiffness(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y);
    void MltStiffness(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    void ComputeAndFactoRealMatrix(const Treal&, const Treal& a,
				   const Treal& b, int which);

    void ComputeAndFactoRealMatrix(const Tcplx&, const Treal& a,
				   const Treal& b, int which);
    
    void ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b,
					    int which =
					    EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b,
                                            int which =
					    EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y);
    void ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    void ComputeSolution(const SeldonTranspose& transA,
                         const Vector<Treal>& X, Vector<Treal>& Y);

    void ComputeSolution(const SeldonTranspose& transA,
                         const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void DistributeEigenvectors(Matrix<T, General, ColMajor>& eigen_vec);
    
    void Clear();
    
  };
 
    
  //! list of availables eigenvalue solvers
  class TypeEigenvalueSolver
  {
  public :
    enum {DEFAULT, ARPACK, ANASAZI, FEAST, SLEPC};
    
    static int GetDefaultSolver();
    
  };
  
  template<class T, class Prop, class Storage>
  void GetEigenvaluesEigenvectors(EigenProblem_Base<T>& var_eig,
				  Vector<T>& lambda, Vector<T>& lambda_imag,
				  Matrix<T, Prop, Storage>& eigen_vec,
				  int type_solver = TypeEigenvalueSolver::DEFAULT);
  
}

#define SELDON_FILE_VIRTUAL_EIGENVALUE_SOLVER_HXX
#endif

