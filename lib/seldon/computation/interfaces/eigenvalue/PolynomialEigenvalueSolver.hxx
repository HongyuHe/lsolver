#ifndef SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_HXX

namespace Seldon
{
  //! Parameters for Slepc package
  class SlepcParamPep
  {
  protected:
    //! which solver ?
    int type_solver;

  public :
    enum {TOAR, STOAR, QARNOLDI, LINEAR, JD};

    SlepcParamPep();
    
    int GetEigensolverType() const;
    void SetEigensolverType(int type);

  };

  //! Base class for polynomial eigenvalue solver
  template<class T>
  class PolynomialEigenProblem_Base
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
    
    //! polynomial degree
    int pol_degree;

#ifdef SELDON_WITH_MPI
    //! communicator used to compute eigenvalues (linear system is split into processors)
    MPI_Comm comm;
#endif

    //! class for comparing eigenvalues
    EigenvalueComparisonClass<T>* compar_eigenval;
    
    //! Slepc parameters
    SlepcParamPep slepc_param;

    //! threshold for Arpack's iterative process
    Treal stopping_criterion;
    
    //! Maximal number of iterations
    int nb_maximum_iterations;
    
    //! number of linear solves
    int nb_linear_solves; int display_every;
    int print_level;

    //! mass diagonal ?
    bool diagonal_mass;
    Vector<T> invDiag;
    
  public:
    PolynomialEigenProblem_Base();
    virtual ~PolynomialEigenProblem_Base();

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

    SlepcParamPep& GetSlepcParameters();
    
    void SetUserComparisonClass(EigenvalueComparisonClass<T>* ev);

    void SetStoppingCriterion(Treal eps);
    Treal GetStoppingCriterion() const;
    void SetNbMaximumIterations(int n);
    int GetNbMaximumIterations() const;

    int GetM() const;
    int GetGlobalM() const;
    int GetN() const;
    int GetPolynomialDegree() const;
    
    int GetNbLinearSolves() const;
    void IncrementLinearSolves();
    void SetPrintLevel(int);
    int GetPrintLevel() const;
    
    void SetDiagonalMass(bool diag = true);
    bool DiagonalMass();

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

    virtual void ComputeOperator(int num, const Vector<T>& coef);
    virtual void MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
    virtual void FactorizeMass();
    virtual void SolveMass(const Vector<T>& x, Vector<T>& y);
    
    virtual void FactorizeOperator(const Vector<T>& coef);
    virtual void SolveOperator(const Vector<T>& X, Vector<T>& Y);
    
  };

  //! matrix-free implementation of polynomial eigenvalue solver
  template<class T>
  class PolynomialEigenProblem : public PolynomialEigenProblem_Base<T>
  {
  protected:
    Vector<VirtualMatrix<T>* > list_op;
    Vector<Vector<T> > list_coef;
    
  public:
    void InitMatrix(const Vector<VirtualMatrix<T>* >& op, int n = -1);
    
    void ComputeOperator(int num, const Vector<T>& coef);
    void MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
  };

  //! implementation of polynomial eigenvalue solver for dense problem
  template<class T, class Prop, class Storage>
  class PolynomialDenseEigenProblem : public PolynomialEigenProblem<T>
  {
  protected:
    Vector<Matrix<T, Prop, Storage>* > list_mat;
    Matrix<T, Prop, Storage> mat_lu;
    Vector<int> pivot;
    
  public:
    void InitMatrix(const Vector<Matrix<T, Prop, Storage>* >& op);
    
    void FactorizeMass();
    void SolveMass(const Vector<T>& x, Vector<T>& y);

    void FactorizeOperator(const Vector<T>& coef);
    void SolveOperator(const Vector<T>& X, Vector<T>& Y);
    
  };

  //! implementation of sparse polynomial eigenvalue solver
  template<class T, class MatStiff,
           class MatMass = Matrix<T, Symmetric, ArrayRowSymSparse> >
  class PolynomialSparseEigenProblem : public PolynomialEigenProblem<T>
  {
  protected:

#ifdef SELDON_WITH_MPI
    //! LU factorization of sparse matrix
    SparseDistributedSolver<T> mat_lu;
#else
    //! LU factorization of sparse matrix
    SparseDirectSolver<T> mat_lu;
#endif
    
    // stiffness matrices and mass matrix
    Vector<MatStiff*> Kh;
    MatMass* Mh;
    
    // parallel stuff
    bool distributed; int nloc;
    Vector<int> local_col_numbers;
    Vector<int>* ProcSharingRows;
    Vector<Vector<int> >* SharingRowNumbers;
    int nodl_scalar_, nb_unknowns_scal_;

  public:
    PolynomialSparseEigenProblem();
    
    int RetrieveLocalNumbers(MatStiff& K);
    void InitMatrix(Vector<MatStiff*>&, MatMass&);
    
    void FactorizeMass();
    void SolveMass(const Vector<T>& x, Vector<T>& y);
    
    void FactorizeOperator(const Vector<T>& coef);
    void SolveOperator(const Vector<T>& X, Vector<T>& Y);

    void MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
  };

}
 
#define SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_HXX
#endif
