#ifndef SELDON_FILE_NON_LINEAR_EIGENVALUE_SOLVER_HXX

namespace Seldon
{
  //! Parameters for Slepc package
  class SlepcParamNep
  {
  protected:
    //! which solver ?
    int type_solver;

  public :
    enum {RII, SLP, NARNOLDI, CISS, INTERPOL, NLEIGS};
    
    SlepcParamNep();
    
    int GetEigensolverType() const;
    void SetEigensolverType(int type);

  };

  //! Base class for non-linear eigenvalue solver
  template<class T>
  class NonLinearEigenProblem_Base
  {
  public:
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;

    //! parts of the spectrum (near from 0, at infinity or around a given value)
    /*!
      SMALL_EIGENVALUES : seeking eigenvalues near 0
      LARGE_EIGENVALUES : seeking largest eigenvalues
      CENTERED_EIGENVALUES : seeking eigenvalues near the shift sigma
    */
    enum {SMALL_EIGENVALUES, LARGE_EIGENVALUES, CENTERED_EIGENVALUES};

    //! different sorting strategies
    enum {SORTED_REAL, SORTED_IMAG, SORTED_MODULUS, SORTED_USER};

  protected:
    //! number of eigenvalues to be computed
    int nb_eigenvalues_wanted; 
    
    //! which spectrum ? Near from Zero ? Near from Infinity ? or near from a value ?
    int type_spectrum_wanted;
    T shift;
    bool use_spectral_transfo;
    
    //! large eigenvalues because of their real part, imaginary part or magnitude ?
    int type_sort_eigenvalues;
    
    //! size of the problem
    int n_, nglob;
    
#ifdef SELDON_WITH_MPI
    //! communicator used to compute eigenvalues
    MPI_Comm comm;
#endif

    //! class for comparing eigenvalues
    EigenvalueComparisonClass<T>* compar_eigenval;
    
    //! Slepc parameters
    SlepcParamNep slepc_param;

    //! threshold for iterative process
    Treal stopping_criterion;
    
    //! Maximal number of iterations
    int nb_maximum_iterations;
    
    //! number of matrix-vector products
    int nb_prod;
    
  public:
    NonLinearEigenProblem_Base();
    virtual ~NonLinearEigenProblem_Base();

    void Init(int n);

#ifdef SELDON_WITH_MPI
    void SetCommunicator(const MPI_Comm& comm_);
    MPI_Comm& GetCommunicator();
#endif

    int GetNbAskedEigenvalues() const;
    void SetNbAskedEigenvalues(int n);
    bool UseSpectralTransformation() const;
    void SetSpectralTransformation(bool t = true);
    
    int GetTypeSpectrum() const;
    T GetCenterSpectrum() const;
    int GetTypeSorting() const;

    void SetTypeSpectrum(int type, const T& val,
                         int type_sort = SORTED_MODULUS);

    SlepcParamNep& GetSlepcParameters();
    
    void SetUserComparisonClass(EigenvalueComparisonClass<T>* ev);

    void SetStoppingCriterion(Treal eps);
    Treal GetStoppingCriterion() const;
    void SetNbMaximumIterations(int n);
    int GetNbMaximumIterations() const;

    int GetM() const;
    int GetGlobalM() const;
    int GetN() const;
    
    int GetNbMatrixVectorProducts() const;
    void IncrementProdMatVect();

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

    virtual void ComputeOperator(const T& L);
    virtual void MltOperator(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
    virtual void ComputeJacobian(const T& L);
    virtual void MltJacobian(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
    virtual void ComputePreconditioning();
    virtual void ApplyPreconditioning(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
  };
  
}
 
#define SELDON_FILE_NON_LINEAR_EIGENVALUE_SOLVER_HXX
#endif
