#ifndef SELDON_FILE_VIRTUAL_EIGENVALUE_SOLVER_CXX

#include "VirtualEigenvalueSolver.hxx"

#ifdef SELDON_WITH_SLEPC
#include <slepceps.h>
#endif

namespace Seldon
{

  /****************
   * AnasaziParam *
   ****************/


  //! default constructor
  AnasaziParam::AnasaziParam()
  {
    type_solver = SOLVER_LOBPCG;
    ortho_manager = ORTHO_DGKS;
    nb_blocks = 2;
    restart_number = 20;
  }


  //! returns the number of blocks used in blocked solvers
  int AnasaziParam::GetNbBlocks() const
  {
    return nb_blocks;
  }


  //! returns the number of blocks used in blocked solvers
  void AnasaziParam::SetNbBlocks(int n)
  {
    nb_blocks = n;
  }
  
  
  //! returns the restart parameter used in blocked solvers
  int AnasaziParam::GetNbMaximumRestarts() const
  {
    return restart_number;
  }
  

  //! sets the restart parameter used in blocked solvers
  void AnasaziParam::SetNbMaximumRestarts(int m)
  {
    restart_number = m;
  }

  
  //! returns orthogonalization manager set in Anasazi
  int AnasaziParam::GetOrthoManager() const
  {
    return ortho_manager;
  }
  
  
  //! returns the solver used in Anasazi
  int AnasaziParam::GetEigensolverType() const
  {
    return type_solver;
  }


  //! sets the solver used in Anasazi
  void AnasaziParam::SetEigensolverType(int type)
  {
    type_solver = type;
  }
  

  /****************
   * SlepcParam *
   ****************/


  //! Default constructor
  SlepcParam::SlepcParam()
  {
    type_solver = KRYLOVSCHUR;
    block_size = -1; nstep = -1;
    type_extraction = -1; quadrature_rule = -1;
    nstep_inner = -1; nstep_outer = -1;
    npoints = -1;
    moment_size = -1; block_max_size = -1; npart = 1;
    delta_rank = 0.0; delta_spur = 0.0;
    borth = -1; double_exp = -1;
    init_size = -1;
    krylov_restart = -1; restart_number = -1; restart_add = -1;
    nb_conv_vector = -1; nb_conv_vector_proj = -1;
    non_locking_variant = false;
    restart_ratio = -1.0;
  }
  
  
  // definition available in Slepc.cxx if Slepc is enabled
#ifndef SELDON_WITH_SLEPC
  const char* SlepcParam::GetEigensolverChar() const
  {
    return "";
  }
#endif

  
  /**************
   * FeastParam *
   **************/


  FeastParam::FeastParam()
  {
    evaluate_number_eigenval = false;
    nb_points_quadrature = 0;
    type_integration = -1;
  }  
  

  /******************
   * Initialization *
   ******************/
  
  
  //! default constructor
  template<class T>
  EigenProblem_Base<T>::EigenProblem_Base()
  {
    eigenvalue_computation_mode = 1;
    nb_eigenvalues_wanted = 0;
    nb_add_eigenvalues = 0;
    // default => we want largest eigenvalues by magnitude
    type_spectrum_wanted = LARGE_EIGENVALUES;
    type_sort_eigenvalues = SORTED_MODULUS;
    
    use_cholesky = false;   
    diagonal_mass = false;
    stopping_criterion = 1e-6;
    nb_maximum_iterations = 1000;
    nb_prod = 0;
    n_ = 0;

    shift = T(0);
    shift_imag = T(0);

    nb_arnoldi_vectors = 0;
    automatic_selection_arnoldi_vectors = true;
    
    print_level = 0;      
    
    complex_system = false;
    selected_part = COMPLEX_PART;

#ifdef SELDON_WITH_MPI
    // for parallel execution, default communicator : all the processors
    comm = MPI_COMM_WORLD;
    comm_global = MPI_COMM_WORLD;
#endif
    
    emin_interval = Treal(0);
    emax_interval = Treal(0);

    center_spectrum = Tcplx(0, 0);
    radius_spectrum = Treal(0);
  }
  

  //! Destructor
  template<class T>
  EigenProblem_Base<T>::~EigenProblem_Base()
  {
    compar_eigenval = NULL;
  }
  
  
  //! initialisation of the size of the eigenvalue problem
  template<class T>
  void EigenProblem_Base<T>::Init(int n)
  {
    n_ = n;
    nb_prod = 0;

    // counting the size of the global system for parallel computation
    nglob = n;
    
#ifdef SELDON_WITH_MPI
    MPI_Allreduce(&n, &nglob, 1, MPI_INTEGER, MPI_SUM, comm);    
#endif

    if (nb_eigenvalues_wanted >= (nglob - 2))
      {
        cout << "Too many wanted eigenvalues " << endl;
        cout << nb_eigenvalues_wanted <<
          " asked eigenvalues, but the rank of the matrix is lower than "
             << n_ << endl;
        
        abort();
      }
    
    if (automatic_selection_arnoldi_vectors)
      nb_arnoldi_vectors = min(nglob, 2*nb_eigenvalues_wanted+2);
    
    //cout << "n = " << n << endl;
    //cout << "nb_arnoldi_vectors = " << nb_arnoldi_vectors << endl;
  }
  
  
  /*******************
   * Basic functions *
   *******************/
  
  
  //! returns the spectral transformation used for evaluation of eigenvalues
  template<class T>
  int EigenProblem_Base<T>::GetComputationalMode() const
  {
    return eigenvalue_computation_mode;
  }
  
  
  //! sets the spectral transformation used for evaluation of eigenvalues
  template<class T>
  void EigenProblem_Base<T>::SetComputationalMode(int mode)
  {
    eigenvalue_computation_mode = mode;
  }
  
  
  //! returns the number of eigenvalues asked by the user
  template<class T>
  int EigenProblem_Base<T>::GetNbAskedEigenvalues() const
  {
    return nb_eigenvalues_wanted;
  }


  //! returns the additional number of eigenvalues
  template<class T>
  int EigenProblem_Base<T>::GetNbAdditionalEigenvalues() const
  {
    return nb_add_eigenvalues;
  }


#ifdef SELDON_WITH_MPI
  //! returns the MPI communicator shared by processors
  template<class T>
  MPI_Comm& EigenProblem_Base<T>::GetCommunicator()
  {
    return comm;
  }
  
  //! sets the MPI communicator shared by processors
  template<class T>
  void EigenProblem_Base<T>::SetCommunicator(const MPI_Comm& comm_)
  {
    comm = comm_;
  }


  //! returns the MPI communicator used by Feast
  template<class T>
  MPI_Comm& EigenProblem_Base<T>::GetGlobalCommunicator()
  {
    return comm_global;
  }
  
  //! sets the MPI communicator used by Feast
  template<class T>
  void EigenProblem_Base<T>::SetGlobalCommunicator(const MPI_Comm& comm_)
  {
    comm_global = comm_;
  }
#endif

  
  template<class T>
  int EigenProblem_Base<T>::GetRankCommunicator() const
  {
#ifdef SELDON_WITH_MPI
    int rank_proc;
    MPI_Comm_rank(comm, &rank_proc);
    return rank_proc;
#endif

    return 0;
  }


  template<class T>
  int EigenProblem_Base<T>::GetGlobalRankCommunicator() const
  {
#ifdef SELDON_WITH_MPI
    int rank_proc;
    MPI_Comm_rank(comm_global, &rank_proc);
    return rank_proc;
#endif
    
    return 0;
  }

  
  //! sets the number of eigenvalues to compute
  template<class T>
  void EigenProblem_Base<T>::SetNbAskedEigenvalues(int n)
  {
    nb_eigenvalues_wanted = n;
  }


  //! sets the number of additional eigenvalues
  template<class T>
  void EigenProblem_Base<T>::SetNbAdditionalEigenvalues(int n)
  {
    nb_add_eigenvalues = n;
  }
  

  //! returns parameters specific to Anasazi
  template<class T>
  AnasaziParam& EigenProblem_Base<T>::GetAnasaziParameters()
  {
    return anasazi_param;
  }

  
  //! returns parameters specific to Slepc
  template<class T>
  SlepcParam& EigenProblem_Base<T>::GetSlepcParameters()
  {
    return slepc_param;
  }
  

  //! returns parameters specific to Feast
  template<class T>
  FeastParam& EigenProblem_Base<T>::GetFeastParameters()
  {
    return feast_param;
  }
  

  //! returns the spectrum desired (large, small eigenvalues, etc)
  template<class T>
  int EigenProblem_Base<T>::GetTypeSpectrum() const
  {
    return type_spectrum_wanted;
  }

  
  //! returns how eigenvalues are sorted (real, imaginary part or modulus)
  template<class T>
  int EigenProblem_Base<T>::GetTypeSorting() const
  {
    return type_sort_eigenvalues;
  }

  
  //! returns the shift value used
  /*!
    If type_spectrum_wanted is set to CENTERED_EIGENVALUES,
    we search closest eigenvalues to the shift value.
    Matrix (A - (shift + i shift_imag)*I)^{-1} will be used instead of A
  */
  template<class T>
  T EigenProblem_Base<T>::GetShiftValue() const
  {
    return shift;
  }

  
  //! returns the imaginary part of shift value used
  /*!
    If type_spectrum_wanted is set to CENTERED_EIGENVALUES,
    we search closest eigenvalues to the shift value.
    Matrix (A - (shift + i shift_imag)*I)^{-1} will be used instead of A
    shift_imag is accessed only for real unsymmetric problems
  */
  template<class T>
  T EigenProblem_Base<T>::GetImagShiftValue() const
  {
    return shift_imag;
  }
  
  
  //! Sets the real part of shift value
  template<class T>
  void EigenProblem_Base<T>::SetShiftValue(const T& val)
  {
    shift = val;
  }

  
  //! Sets the imaginary part of shift value
  template<class T>
  void EigenProblem_Base<T>::SetImagShiftValue(const T& val)
  {
    shift_imag = val;
  }

  
  //! forms the complex shift from real and imaginary part
  template<class T>
  void EigenProblem_Base<T>
  ::GetComplexShift(const Treal& sr, const Treal& si, Tcplx& s) const
  {
    s = complex<Treal>(sr, si);
  }


  //! forms the complex shift from real and imaginary part
  template<class T>
  void EigenProblem_Base<T>
  ::GetComplexShift(const Tcplx& sr, const Tcplx& si, Tcplx& s) const
  {
    s = sr;
  }


  //! sets which eigenvalues are searched
  /*!
    You can ask small eigenvalues, large, or eigenvalues
    close to the shift.
  */
  template<class T>
  void EigenProblem_Base<T>
  ::SetTypeSpectrum(int type, const T& val, int type_sort)
  {
    type_spectrum_wanted = type;
    shift = val;
    type_sort_eigenvalues = type_sort;
  }

  
  //! sets which eigenvalues are searched
  /*!
    You can ask small eigenvalues, large, or eigenvalues
    close to the shift.
  */
  template<class T>
  void EigenProblem_Base<T>
  ::SetTypeSpectrum(int type, const complex<T>& val, int type_sort)
  {
    // for real unsymmetric eigenproblems, you can
    // specify a complex shift
    type_spectrum_wanted = type;
    shift = real(val);
    shift_imag = imag(val);
        
    type_sort_eigenvalues = type_sort;
  }
  

  //! sets the class where two eigenvalues can be compared
  template<class T>
  void EigenProblem_Base<T>::SetUserComparisonClass(EigenvalueComparisonClass<T>* ev)
  {
    compar_eigenval = ev;
  }
  
  
  //! returns lower bound of the interval where eigenvalues are searched
  template<class T>
  typename ClassComplexType<T>::Treal EigenProblem_Base<T>
  ::GetLowerBoundInterval() const
  {
    return emin_interval;
  }


  //! returns upper bound of the interval where eigenvalues are searched
  template<class T>
  typename ClassComplexType<T>::Treal EigenProblem_Base<T>
  ::GetUpperBoundInterval() const
  {
    return emax_interval;
  }
  
  
  //! sets the interval where eigenvalues are searched
  template<class T>
  void EigenProblem_Base<T>
  ::SetIntervalSpectrum(typename ClassComplexType<T>::Treal l0,
			typename ClassComplexType<T>::Treal l1)
  {
    emin_interval = l0;
    emax_interval = l1;
  }
    
  
  //! sets a circle where eigenvalues are searched
  template<class T>
  void EigenProblem_Base<T>::SetCircleSpectrum(const Tcplx& z, const Treal& r)
  {
    center_spectrum = z; radius_spectrum = r;
    ratio_ellipse = 100; angle_ellipse = 0;
  }


  //! sets an ellipse where eigenvalues are searched
  template<class T>
  void EigenProblem_Base<T>::SetEllipseSpectrum(const Tcplx& z, const Treal& r,
                                                const Treal& ratio, const Treal& teta)
  {
    center_spectrum = z; radius_spectrum = r;
    ratio_ellipse = ratio; angle_ellipse = teta;
  }
  

  //! indicates the use of Cholesky factorisation in order to 
  //! solve a standard eigenvalue problem instead of a generalized one
  template<class T>
  void EigenProblem_Base<T>
  ::SetCholeskyFactoForMass(bool chol)
  {
    use_cholesky = chol;
  }
  
  
  //! returns true if Cholesky factorisation has to be used for mass matrix
  template<class T>
  bool EigenProblem_Base<T>::UseCholeskyFactoForMass() const
  {
    return use_cholesky;
  }
    
  
  //! indicates that the mass matrix is diagonal
  template<class T>
  void EigenProblem_Base<T>::SetDiagonalMass(bool diag)
  {
    diagonal_mass = diag;
  }
  
  
  //! returns true if the mass matrix is diagonal
  template<class T>
  bool EigenProblem_Base<T>::DiagonalMass() const
  {
    return diagonal_mass;
  }
  
    
  //! modifies the stopping critertion
  template<class T>
  void EigenProblem_Base<T>
  ::SetStoppingCriterion(typename ClassComplexType<T>::Treal eps)
  {
    stopping_criterion = eps;
  }
  
  
  //! returns the stopping criterion
  template<class T>
  typename ClassComplexType<T>::Treal EigenProblem_Base<T>
  ::GetStoppingCriterion() const
  {
    return stopping_criterion;
  }
    
  
  //! sets the maximal number of iterations allowed for the iterative algorithm
  template<class T>
  void EigenProblem_Base<T>::SetNbMaximumIterations(int n)
  {
    nb_maximum_iterations = n;
  }
  
  
  //! returns the maximal number of iterations allowed for the iterative algorithm
  template<class T>
  int EigenProblem_Base<T>::GetNbMaximumIterations() const
  {
    return nb_maximum_iterations;
  }
  
  
  //! returns the number of matrix-vector products performed 
  //! since last call to Init
  template<class T>
  int EigenProblem_Base<T>::GetNbMatrixVectorProducts() const
  {
    return nb_prod;
  }
    
  
  //! returns the number of Arnoldi vectors to use
  template<class T>
  int EigenProblem_Base<T>::GetNbArnoldiVectors() const
  {
    return nb_arnoldi_vectors;
  }
  
  
  //! sets the number of Arnoldi vectors to use
  template<class T>
  void EigenProblem_Base<T>::SetNbArnoldiVectors(int n)
  {
    automatic_selection_arnoldi_vectors = false;
    nb_arnoldi_vectors = n;
  }
  
  
  //! returns number of rows
  template<class T>
  int EigenProblem_Base<T>::GetM() const
  {
    return n_;
  }


  //! returns global number of rows
  template<class T>
  int EigenProblem_Base<T>::GetGlobalM() const
  {
    return nglob;
  }

  
  //! returns number of columns
  template<class T>
  int EigenProblem_Base<T>::GetN() const
  {
    return n_;
  }
  
  
  //! returns level of verbosity
  template<class T>
  int EigenProblem_Base<T>::GetPrintLevel() const
  {
    return print_level;
  }


  //! returns level of verbosity
  template<class T>
  int EigenProblem_Base<T>::GetGlobalPrintLevel() const
  {
    return global_print_level;
  }
  
  
  //! sets the level of verbosity
  template<class T>
  void EigenProblem_Base<T>::SetPrintLevel(int lvl)
  {
    print_level = lvl;
  }


  //! sets the level of verbosity
  template<class T>
  void EigenProblem_Base<T>::SetGlobalPrintLevel(int lvl)
  {
    global_print_level = lvl;
  }
  
  
  //! increment of the number of matrix vector products
  template<class T>
  void EigenProblem_Base<T>::IncrementProdMatVect()
  {
    nb_prod++;
#ifdef SELDON_WITH_MPI
    int rank; MPI_Comm_rank(comm, & rank);
#else
    int rank(0);
#endif

    if (print_level >= 3)
      {
        if (nb_prod%10 == 0)
          if (rank == 0)
            cout<<" Iteration number " << nb_prod << endl;
      }
    else if (print_level >= 1)
      {
        if (nb_prod%100 == 0)
          if (rank == 0)
            cout<<" Iteration number " << nb_prod << endl;
      }			
  }
  

  //! prints error of initialization and aborts program
  template<class T>
  void EigenProblem_Base<T>::PrintErrorInit() const
  {
    cout << "InitMatrix has not been called" << endl;
    abort();
  }
  
  
  //! computation of mass matrix
  template<class T>
  void EigenProblem_Base<T>::ComputeMassForCholesky()
  {
    // nothing to do, we consider that mass matrix
    // is already computed
  }
  
  
  //! computation of mass matrix M
  template<class T>
  void EigenProblem_Base<T>::ComputeMassMatrix()
  {
    // mass matrix already computed in Mh
  }
  
  
  //! computation of stiffness matrix K
  template<class T>
  void EigenProblem_Base<T>::ComputeStiffnessMatrix()
  {
    // nothing to do, already computed in Kh
  }
  
  
  //! computation of matrix a M + b*K
  template<class T>
  void EigenProblem_Base<T>
  ::ComputeStiffnessMatrix(const T& a, const T& b)
  {
    // nothing to do, we use Kh and Mh for the matrix vector product
  }
  
  
  //! computation of matrix a M + b K and factorisation of this matrix
  /*!
    The factorisation process can be also the construction of preconditioning
    if an iterative solver is used to solve linear system (a M + b K) y = x 
  */
  template<class T>
  void EigenProblem_Base<T>
  ::ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b,
				       int which_part)
  {
    abort();
  }
  
  
  //! computation of matrix a M + b K and factorisation of this matrix
  /*!
    The factorisation process can be also the construction of preconditioning
    if an iterative solver is used to solve linear system (a M + b K) y = x 
  */
  template<class T>
  void EigenProblem_Base<T>
  ::ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b,
				       int which_part)
  {
    abort();
  }
  
  
  //! solving the linear system (a M + b K) Y = X
  template<class T>
  void EigenProblem_Base<T>
  ::ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    abort();
  }
  

  //! solving the linear system (a M + b K) Y = X
  template<class T>
  void EigenProblem_Base<T>
  ::ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    abort();
  }

  
  //! solving the linear system (a M + b K) Y = X
  template<class T>
  void EigenProblem_Base<T>
  ::ComputeSolution(const SeldonTranspose&,
		    const Vector<Treal>& X, Vector<Treal>& Y)
  {
    abort();
  }


  //! solving the linear system (a M + b K) Y = X
  template<class T>
  void EigenProblem_Base<T>
  ::ComputeSolution(const SeldonTranspose&,
		    const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    abort();
  }

  
  //! computation of Cholesky factorisation of M from matrix M
  template<class T>
  void EigenProblem_Base<T>::FactorizeCholeskyMass()
  {
    abort();
  }
  
  
  //! computation of L X or L^T x if M = L L^T
  template<class T>
  void EigenProblem_Base<T>
  ::MltCholeskyMass(const SeldonTranspose& TransA, Vector<Treal>& X)
  {
    abort();
  }


  //! computation of L X or L^T x if M = L L^T
  template<class T>
  void EigenProblem_Base<T>
  ::MltCholeskyMass(const SeldonTranspose& TransA, Vector<Tcplx>& X)
  {
    abort();
  }
  
  
  //! computation of L^-1 X or L^-T x if M = L L^T
  template<class T>
  void EigenProblem_Base<T>
  ::SolveCholeskyMass(const SeldonTranspose& TransA, Vector<Treal>& X)
  {
    abort();
  }


  //! computation of L^-1 X or L^-T x if M = L L^T
  template<class T>
  void EigenProblem_Base<T>
  ::SolveCholeskyMass(const SeldonTranspose& TransA, Vector<Tcplx>& X)
  {
    abort();
  }


  //! changes final eigenvectors if needed
  template<class T>
  void EigenProblem_Base<T>::DistributeEigenvectors(Matrix<T, General, ColMajor>& eigen_vec)
  {
  }

  
  //! memory release
  template<class T>
  void EigenProblem_Base<T>::Clear()
  {
  }

#ifdef SELDON_WITH_SLEPC
  template<class T> template<class T0>
  PetscErrorCode EigenProblem_Base<T>::
  GetComparisonEigenvalueSlepcGen(PetscScalar, PetscScalar, PetscScalar, PetscScalar, T0,
                                  PetscInt*, void*)
  {
    cout << "Incompatibles types between T and PetscScalar" << endl;
    abort();
    return 0;
  }

  
  template<class T>
  PetscErrorCode EigenProblem_Base<T>::
  GetComparisonEigenvalueSlepcGen(PetscScalar Lr, PetscScalar Li, PetscScalar Lr2, PetscScalar Li2, PetscScalar,
                                  PetscInt* res, void* ctx)
  {
    EigenProblem_Base<T>& var = *reinterpret_cast<EigenProblem_Base<T>* >(ctx);
    if (var.compar_eigenval == NULL)
      {
        cout << "Invalid pointer for compar_eigenval" << endl;
        abort();
      }
    
    T Lr_, Li_, Lr2_, Li2_;
    to_complex_eigen(Lr, Lr_); to_complex_eigen(Li, Li_);
    to_complex_eigen(Lr2, Lr2_); to_complex_eigen(Li2, Li2_);
    int y = var.compar_eigenval->CompareEigenvalue(Lr_, Li_, Lr2_, Li2_);
    *res = y;
    return 0;
  }

  template<class T>
  PetscErrorCode EigenProblem_Base<T>::
  GetComparisonEigenvalueSlepc(PetscScalar Lr, PetscScalar Li, PetscScalar Lr2, PetscScalar Li2,
                               PetscInt* res, void* ctx)
  {
    T z; SetComplexZero(z);
    GetComparisonEigenvalueSlepcGen(Lr, Li, Lr2, Li2, z, res, ctx);
    return 0;
  }
#endif


  template<>
  void EigenProblem_Base<double>
  ::FillComplexEigenvectors(int m, const complex<double>& Emid, double eps,
                            const Vector<complex<double> >& lambda_cplx,
                            const Matrix<complex<double>, General, ColMajor>& Ecplx,
                            Vector<double>& Lr, Vector<double>& Li,
                            Matrix<double, General, ColMajor>& E)
  {
    // if imag_pos is true, we will select only eigenvalues with Im(lambda) >= 0
    // if imag_pos is false, we will select only eigenvalues with Im(lambda) <= 0
    bool imag_pos = true;
    if (imag(Emid) < 0)
      imag_pos = false;
    
    // type_eigenval : 0 for real, 1 for complex, -1 to be dropped
    // some eigenvalues are dropped because by construction eigenvalues are complex conjugate
    Vector<int> type_eigenval(m);
    type_eigenval.Fill(-1);
    int nb_val = 0;
    for (int i = 0; i < m; i++)
      {
        if (abs(imag(lambda_cplx(i))) <= eps)
          {
            type_eigenval(i) = 0;
            nb_val++;
          }
        else
          {
            if (imag_pos)
              {
                if (imag(lambda_cplx(i)) > 0)
                  {
                    type_eigenval(i) = 1;
                    nb_val += 2;
                  }
              }
            else
              {
                if (imag(lambda_cplx(i)) < 0)
                  {
                    type_eigenval(i) = 1;
                    nb_val += 2;
                  }
              }
          }
      }
    
    // we fill output arrays Lr, Li and E
    int N = Ecplx.GetM();
    Lr.Reallocate(nb_val);
    Li.Reallocate(nb_val);
    E.Reallocate(N, nb_val);
    nb_val = 0;
    for (int i = 0; i < m; i++)
      {
        if (type_eigenval(i) == 0)
          {
            Lr(nb_val) = realpart(lambda_cplx(i));
            Li(nb_val) = 0.0;
            for (int j = 0; j < N; j++)
              E(j, nb_val) = realpart(Ecplx(j, i));
            
            nb_val++;
          }
        else if (type_eigenval(i) == 1)
          {
            Lr(nb_val) = real(lambda_cplx(i));
            Li(nb_val) = imag(lambda_cplx(i));
            Lr(nb_val+1) = Lr(nb_val);
            Li(nb_val+1) = -Li(nb_val);
            for (int j = 0; j < N; j++)
              {
                E(j, nb_val) = real(Ecplx(j, i));
                E(j, nb_val+1) = imag(Ecplx(j, i));
              }

            nb_val += 2;
          }
      }
  }


  template<>
  void EigenProblem_Base<complex<double> >
  ::FillComplexEigenvectors(int m, const complex<double>& Emid, double eps,
                            const Vector<complex<double> >& lambda_cplx,
                            const Matrix<complex<double>, General, ColMajor>& Ecplx,
                            Vector<complex<double> >& Lr, Vector<complex<double> >& Li,
                            Matrix<complex<double>, General, ColMajor>& E)
  {
    // we fill output arrays Lr and E
    int N = Ecplx.GetM();
    Lr.Reallocate(m);
    E.Reallocate(N, m);
    for (int i = 0; i < m; i++)
      {
        Lr(i) = lambda_cplx(i);
        for (int j = 0; j < N; j++)
          E(j, i) = Ecplx(j, i);
      }
  }
  

  /***********************
   * VirtualEigenProblem *
   ***********************/
  

  //! default constructor
  template<class T, class StiffValue, class MassValue>
  VirtualEigenProblem<T, StiffValue, MassValue>::VirtualEigenProblem()
  {
    Mh = NULL;
    Kh = NULL;
  }
  

  //! initialization of a standard eigenvalue problem
  /*!
    Stiffness matrix K is given in argument.
    we will search (lambda, x) such as K x = lambda x
    n is an optional argument for parallel computation (local number of rows)
  */
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::InitMatrix(VirtualMatrix<StiffValue>& K, int n)
  {
    Kh = &K;
    Mh = NULL;
    this->diagonal_mass = true;
    if ( (!K.IsSymmetric()) && (!K.IsComplex()) && (this->shift_imag != T(0)) )
      {
        // for real unsymmetric problems, if sigma is complex
        // we have to use mode 3 or 4 in Arpack => generalized problem
        this->diagonal_mass = false;
      }

    if (n == -1)
      this->Init(K.GetM());
    else
      this->Init(n);
  }
  
  
  //! initialization of a generalized eigenvalue problem
  /*!
    Mass matrix M and stiffness matrix K are given in argument
    we will search (lambda, x) such as K x = lambda M x
  */
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::InitMatrix(VirtualMatrix<StiffValue>& K, VirtualMatrix<MassValue>& M, int n)
  {
    Kh = &K;
    Mh = &M;
    this->diagonal_mass = false;

    if (n == -1)
      this->Init(K.GetM());
    else
      this->Init(n);
  }
  
  
  //! sets pointers to the stiffness and mass matrix
  /*!
    Pointers for mass and stiffness are set, but not the number of rows
    Init(n) must be called thereafter
   */
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::SetMatrix(VirtualMatrix<StiffValue>& K, VirtualMatrix<MassValue>& M)
  {
    Kh = &K;
    Mh = &M;
  }
  
		  
  //! sets which eigenvalues are searched
  /*!
    You can ask small eigenvalues, large, or eigenvalues
    close to the shift.
  */
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::SetTypeSpectrum(int type, const T& val, int type_sort)
  {
    EigenProblem_Base<T>::SetTypeSpectrum(type, val, type_sort);
  }
  
  
  //! sets which eigenvalues are searched
  /*!
    You can ask small eigenvalues, large, or eigenvalues
    close to the shift.
  */
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::SetTypeSpectrum(int type, const complex<T>& val, int type_sort)
  {
    EigenProblem_Base<T>::SetTypeSpectrum(type, val, type_sort);
    
    if (Kh != NULL)
      {
        if ( (!Kh->IsSymmetric())
             && (!Kh->IsComplex()) && (this->shift_imag != T(0)) )
          {
            // for real unsymmetric problems, if sigma is complex
            // we have to use mode 3 or 4 in Arpack => generalized problem
            this->diagonal_mass = false;
          }
      }
  }
  
  
  //! returns true if the matrix is symmetric
  template<class T, class StiffValue, class MassValue>
  bool VirtualEigenProblem<T, StiffValue, MassValue>::IsSymmetricProblem() const
  {
    if (Kh != NULL)
      {
        if (Kh->IsSymmetric())
          {
            if (Mh == NULL)
              return true;
            else
              return Mh->IsSymmetric();
          }
      }
    else
      this->PrintErrorInit();
    
    return false;
  }


  //! returns true if the matrix is hermitian
  template<class T, class StiffValue, class MassValue>
  bool VirtualEigenProblem<T, StiffValue, MassValue>::IsHermitianProblem() const
  {
    if (Kh != NULL)
      {
	if (Kh->IsComplex())
	  return false;
	else
	  {
            if (Kh->IsSymmetric())
              {
                if (Mh == NULL)
                  return true;
                else
                  {
		    if (Mh->IsComplex())
		      return false;

		    return Mh->IsSymmetric();
		  }
              }
          }
      }
    else
      this->PrintErrorInit();
    
    return false;
  }


  //! computation of diagonal of mass matrix
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>::ComputeDiagonalMass()
  {
    Vector<MassValue>& D = sqrt_diagonal_mass;
    if (Mh == NULL)
      {
        // M = identity
        D.Reallocate(this->n_);
        D.Fill(1.0);
      }
    else
      {
	cout << "not implemented for this kind of eigenproblem" << endl;
	abort();
      }
  }

  
  //! computation of D^1/2 from D
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>::FactorizeDiagonalMass()
  {
    Vector<MassValue>& D = sqrt_diagonal_mass;
    sqrt_diagonal_mass.Reallocate(D.GetM());
    for (int i = 0; i < D.GetM(); i++)
      sqrt_diagonal_mass(i) = sqrt(D(i));
  }

  
  //! fills D^1/2
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>::GetSqrtDiagonal(Vector<T>& D)
  {
    D.Reallocate(sqrt_diagonal_mass.GetM());
    for (int i = 0; i < D.GetM(); i++)
      D(i) = sqrt_diagonal_mass(i);
  }
  
  
  //! multiplication of X by D^-1/2
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltInvSqrtDiagonalMass(Vector<Treal>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) /= realpart(sqrt_diagonal_mass(i));
  }
  
  
  //! multiplication of X by D^1/2
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltSqrtDiagonalMass(Vector<Treal>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) *= realpart(sqrt_diagonal_mass(i));
  }



  //! multiplication of X by D^-1/2
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltInvSqrtDiagonalMass(Vector<Tcplx>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) /= sqrt_diagonal_mass(i);
  }
  
  
  //! multiplication of X by D^1/2
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltSqrtDiagonalMass(Vector<Tcplx>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) *= sqrt_diagonal_mass(i);
  }


  //! matrix vector product with mass matrix Y = M X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltMass(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mh->MltVector(X, Y);
  }


  //! matrix vector product with mass matrix Y = M X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mh->MltVector(X, Y);
  }


  //! matrix vector product with mass matrix Y = M^T X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltMass(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mh->MltVector(trans, X, Y);
  }


  //! matrix vector product with mass matrix Y = M^T X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltMass(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mh->MltVector(trans, X, Y);
  }
  

  //! matrix vector product with stifness matrix Y = K X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Kh->MltVector(X, Y);
  }


  //! matrix vector product with stifness matrix Y = K X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Kh->MltVector(X, Y);
  }
  
  
  //! matrix vector product with stifness and mass matrix Y = (a M + b K) X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltStiffness(const T& coef_massb, const T& coef_stiffb,
		 const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Kh->MltVector(X, Y);
    
    Treal coef_mass, coef_stiff;
    SetComplexReal(coef_massb, coef_mass);
    SetComplexReal(coef_stiffb, coef_stiff);
    if (coef_mass != T(0))
      {
        if (Mh == NULL)
          for (int i = 0; i < Y.GetM(); i++)
            Y(i) += coef_mass*X(i);
        else
          Mh->MltAddVector(coef_mass, X, coef_stiff, Y);
      }
    else
      {
        if (coef_stiff != T(1))
          Mlt(coef_stiff, Y);
      }
  }


  //! matrix vector product with stifness and mass matrix Y = (a M + b K) X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltStiffness(const T& coef_mass, const T& coef_stiff,
		 const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Kh->MltVector(X, Y);
    
    if (coef_mass != T(0))
      {
        if (Mh == NULL)
          for (int i = 0; i < Y.GetM(); i++)
            Y(i) += coef_mass*X(i);
        else
          Mh->MltAddVector(coef_mass, X, coef_stiff, Y);
      }
    else
      {
        if (coef_stiff != T(1))
          Mlt(coef_stiff, Y);
      }
  }
  

  //! matrix vector product with stifness matrix Y = K^T X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltStiffness(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Kh->MltVector(trans, X, Y);
  }


  //! matrix vector product with stifness matrix Y = K^T X
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>
  ::MltStiffness(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Kh->MltVector(trans, X, Y);
  }


  //! memory release
  template<class T, class StiffValue, class MassValue>
  void VirtualEigenProblem<T, StiffValue, MassValue>::Clear()
  {
    sqrt_diagonal_mass.Clear();
  }


  /********************************************
   * Modification of eigenvalues/eigenvectors *
   ********************************************/
  
  
  //! modification of eigenvectors to take into account scaling by mass matrix
  /*!
    One may desire to use matrix D^-1/2 K D^-1/2 or L^-1 K L^-T instead of K
    in order to solve a standard eigenvalue problem instead of a generalized one.
    => eigenvectors are recovered by multiplying them by matrix D^1/2 or by L^T
      with this function
  */
  template<class T0, class Prop, class Storage>
  void ApplyScalingEigenvec(EigenProblem_Base<T0>& var,
			    Vector<T0>& eigen_values, Vector<T0>& lambda_imag,
                            Matrix<T0, Prop, Storage>& eigen_vectors,
                            const T0& shiftr, const T0& shifti)
  {    
    if (var.DiagonalMass())
      {
	Vector<T0> sqrt_Dh;
	var.GetSqrtDiagonal(sqrt_Dh);
        // scaling to have true eigenvectors
        for (int i = 0; i < eigen_vectors.GetM(); i++)
          for (int j = 0; j < eigen_vectors.GetN(); j++)
            eigen_vectors(i, j) /= sqrt_Dh(i);
      }      
    else if (var.UseCholeskyFactoForMass())
      {
        Vector<T0> Xcol(eigen_vectors.GetM());
        for (int j = 0; j < eigen_vectors.GetN(); j++)
          {
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              Xcol(i) = eigen_vectors(i, j);
            
            var.SolveCholeskyMass(SeldonTrans, Xcol);
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              eigen_vectors(i, j) = Xcol(i);
          }
      }

    var.DistributeEigenvectors(eigen_vectors);
    
    if (var.GetComputationalMode() != var.REGULAR_MODE)
      {
        bool use_rayleigh_coef = false;
        if ( (var.eigenvalue_computation_mode == var.INVERT_MODE)
             && (var.GetTypeSpectrum() != var.CENTERED_EIGENVALUES))
          {
            // nothing to change
          }
        else if ((var.DiagonalMass())|| (var.UseCholeskyFactoForMass())
            || (var.eigenvalue_computation_mode == var.INVERT_MODE))
          {
            if ( (var.GetImagShiftValue() != T0(0)) && (var.eigenvalue_computation_mode == var.INVERT_MODE))
              use_rayleigh_coef = true;

            // shift-invert mode, we have to modify eigenvalues
            if (!use_rayleigh_coef)
              for (int i = 0; i < eigen_values.GetM(); i++)
                {
                  if ((eigen_values(i) == 0) && (lambda_imag(i) == 0))
                    {
                      eigen_values(i) = 0;
                      lambda_imag(i) = 0;
                    }
                  else
                    {
                      complex<T0> val = 1.0/complex<T0>(eigen_values(i), lambda_imag(i))
                        + complex<T0>(shiftr, shifti);
                      
                      eigen_values(i) = real(val);
                      lambda_imag(i) = imag(val);
                    }
                }            
          }
        else if (var.GetImagShiftValue() != T0(0))
          {
            use_rayleigh_coef = true;
          }

        if (use_rayleigh_coef)
          {
            int n = eigen_vectors.GetM();
            Vector<T0> X(n), Ax(n), Mx(n), Y(n);
            int j = 0;
            Ax.Fill(T0(0));
            Mx.Fill(T0(0));
            while (j < eigen_values.GetM())
              {
                if (lambda_imag(j) == T0(0))
                  {
                    // real eigenvalue
                    // lambda is retrieved by computing Rayleigh quotient
                    for (int i = 0; i < eigen_vectors.GetM(); i++)
                      X(i) = eigen_vectors(i,j);
                    
                    var.MltMass(X, Mx);
                    var.MltStiffness(X, Ax);
                    eigen_values(j) = DotProd(X, Ax)/DotProd(X, Mx);
                    
                    // next eigenvalue
                    j++;
                  }
                else
                  {
                    if (j == eigen_values.GetM() - 1)
                      {
                        eigen_values(j) = 0.0;
                        lambda_imag(j) = 0.0;
                          
                        break;
                      }
                    
                    // conjugate pair of eigenvalues
                    for (int i = 0; i < eigen_vectors.GetM(); i++)
                      {
                        X(i) = eigen_vectors(i, j);
                        Y(i) = eigen_vectors(i, j+1);
                      }
                    
                    // complex Rayleigh quotient
                    var.MltStiffness(X, Ax);
                    T0 numr = DotProd(X, Ax);
                    T0 numi = DotProd(Y, Ax);
                    
                    var.MltStiffness(Y, Ax);
                    numr += DotProd(Y, Ax);
                    numi -= DotProd(X, Ax);
                    
                    var.MltMass(X, Mx);
                    T0 denr = DotProd(X, Mx);
                    T0 deni = DotProd(Y, Mx);
                    
                    var.MltMass(Y, Mx);
                    denr += DotProd(Y, Mx);
                    deni -= DotProd(X, Mx);
                    
                    complex<T0> val = complex<T0>(numr, numi)/complex<T0>(denr, deni);
                    
                    eigen_values(j) = real(val);
                    eigen_values(j+1) = real(val);

                    lambda_imag(j) = -imag(val);
                    lambda_imag(j+1) = imag(val);
                    
                    // next eigenvalue
                    j += 2;
                  }
              }
          }
      }
  }


  //! modification of eigenvectors to take into account scaling by mass matrix
  /*!
    One may desire to use matrix D^-1/2 K D^-1/2 or L^-1 K L^-T instead of K
    in order to solve a standard eigenvalue problem instead of a generalized one.
    => eigenvectors are recovered by multiplying them by matrix D^1/2 or by L^T
      with this function
  */
  template<class T0, class Prop, class Storage>
  void ApplyScalingEigenvec(EigenProblem_Base<complex<T0> >& var,
			    Vector<complex<T0> >& eigen_values,
			    Vector<complex<T0> >& lambda_imag,
                            Matrix<complex<T0>, Prop, Storage>& eigen_vectors,
                            const complex<T0>& shiftr, const complex<T0>& shifti)
  {
    
    if (var.DiagonalMass())
      {
	Vector<complex<T0> > sqrt_Dh;
	var.GetSqrtDiagonal(sqrt_Dh);
        // scaling to have true eigenvectors
        for (int i = 0; i < eigen_vectors.GetM(); i++)
          for (int j = 0; j < eigen_vectors.GetN(); j++)
            eigen_vectors(i, j) /= sqrt_Dh(i);
      }      
    else if (var.UseCholeskyFactoForMass())
      {
        Vector<complex<T0> > Xcol(eigen_vectors.GetM());
        for (int j = 0; j < eigen_vectors.GetN(); j++)
          {
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              Xcol(i) = eigen_vectors(i, j);
            
            var.SolveCholeskyMass(SeldonTrans, Xcol);
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              eigen_vectors(i, j) = Xcol(i);
          }
      }

    if (var.GetComputationalMode() != var.REGULAR_MODE)
      {
        if ( (var.eigenvalue_computation_mode == var.INVERT_MODE)
             && (var.GetTypeSpectrum() != var.CENTERED_EIGENVALUES))
          {
            // nothing to change
	  }
        else if ((var.DiagonalMass())|| (var.UseCholeskyFactoForMass())
            || (var.eigenvalue_computation_mode == var.INVERT_MODE))
          {
            // shift-invert mode, we have to modify eigenvalues
            for (int i = 0; i < eigen_values.GetM(); i++)
              {
                complex<T0> val = 1.0/eigen_values(i) + shiftr;
                
                eigen_values(i) = val;
              }
            
          }
      }
  }
  
  
  /***********************
   * Sorting eigenvalues *
   ***********************/
  
  
  //! sorting eigenvalues
  template<class T, class Storage1, class Storage2>
  void SortEigenvalues(Vector<T>& lambda_r, Vector<T>& lambda_i,
                       Matrix<T, General, Storage1>& eigen_old,
                       Matrix<T, General, Storage2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const T& shift_r, const T& shift_i)
  {
    int n = min(lambda_r.GetM(), long(eigen_old.GetN()));
				 
    IVect permutation(n);
    permutation.Fill();
    eigen_new.Reallocate(eigen_old.GetM(), n);
    
    // creating a vector that can be sorted
    Vector<T>  L(n); L.Fill();
    if (type_spectrum == EigenProblem_Base<T>::CENTERED_EIGENVALUES)
      {
        // eigenvalues closest to shift are placed at beginning
        switch (type_sort)
          {
          case EigenProblem_Base<T>::SORTED_REAL :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(shift_r - lambda_r(i));
            }
            break;
          case EigenProblem_Base<T>::SORTED_IMAG :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(shift_i - lambda_i(i));
            }
            break;
          case EigenProblem_Base<T>::SORTED_MODULUS :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(complex<T>(shift_r - lambda_r(i), shift_i - lambda_i(i)));
            }
            break;
          }
      }
    else
      {
        // smallest eigenvalues are placed at beginning
        switch (type_sort)
          {
          case EigenProblem_Base<T>::SORTED_REAL :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(lambda_r(i));
            }
            break;
          case EigenProblem_Base<T>::SORTED_IMAG :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(lambda_i(i));
            }
            break;
          case EigenProblem_Base<T>::SORTED_MODULUS :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(complex<T>(lambda_r(i), lambda_i(i)));
            }
            break;            
          }
      }
    
    // sorting L, and retrieving permutation array
    Sort(L, permutation);
    
    // permuting eigenvalues and eigenvectors
    Vector<T> oldLambda_r = lambda_r, oldLambda_i = lambda_i;
    for (int i = 0; i < n; i++)
      {
        lambda_r(i) = oldLambda_r(permutation(i));
        lambda_i(i) = oldLambda_i(permutation(i));
        for (int j = 0; j < eigen_old.GetM(); j++)
          eigen_new(j, i) = eigen_old(j, permutation(i));
      }
    
  }


  //! sorting eigenvalues
  template<class T, class Storage1, class Storage2>
  void SortEigenvalues(Vector<complex<T> >& lambda_r, Vector<complex<T> >& lambda_i,
                       Matrix<complex<T>, General, Storage1>& eigen_old,
                       Matrix<complex<T>, General, Storage2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const complex<T>& shift_r, const complex<T>& shift_i)
  {
    // complex case, ignoring lambda_i and shift_i
    int n = min(lambda_r.GetM(), long(eigen_old.GetN()));
    
    IVect permutation(n);
    permutation.Fill();
    eigen_new.Reallocate(eigen_old.GetM(), n);
    
    // creating a vector that can be sorted
    Vector<T>  L(n);
    if (type_spectrum == EigenProblem_Base<T>::CENTERED_EIGENVALUES)
      {
        // eigenvalues closest to shift are placed at beginning
        switch (type_sort)
          {
          case EigenProblem_Base<T>::SORTED_REAL :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(real(shift_r - lambda_r(i)));
            }
            break;
          case EigenProblem_Base<T>::SORTED_IMAG :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(imag(shift_r - lambda_r(i)));
            }
            break;
          case EigenProblem_Base<T>::SORTED_MODULUS :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(shift_r - lambda_r(i));
            }
            break;
          }
      }
    else
      {
        // smallest eigenvalues are placed at beginning
        switch (type_sort)
          {
          case EigenProblem_Base<T>::SORTED_REAL :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(real(lambda_r(i)));
            }
            break;
          case EigenProblem_Base<T>::SORTED_IMAG :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(imag(lambda_r(i)));
            }
            break;
          case EigenProblem_Base<T>::SORTED_MODULUS :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(lambda_r(i));
            }
            break;            
          }
      }
    
    // sorting L, and retrieving permutation array
    Sort(L, permutation);
    
    // permuting eigenvalues and eigenvectors
    Vector<complex<T> > oldLambda_r = lambda_r, oldLambda_i = lambda_i;
    for (int i = 0; i < n; i++)
      {
        lambda_r(i) = oldLambda_r(permutation(i));
        lambda_i(i) = oldLambda_i(permutation(i));
        for (int j = 0; j < eigen_old.GetM(); j++)
          eigen_new(j, i) = eigen_old(j, permutation(i));
      }
    
  }
  
  
  /*********************
   * DenseEigenProblem *
   *********************/
  
  
  //! default constructor
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  DenseEigenProblem() : VirtualEigenProblem<T, Tstiff, Tmass>()
  {    
    Mh = NULL;
    Kh = NULL;
  }
  
  
  //! Sets stiffness matrix, eigenproblem is K x = lambda x
  template<class T, class Tstiff, class Prop, class Storage,
	   class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::InitMatrix(Matrix<Tstiff, Prop, Storage>& K)
  {
    VirtualEigenProblem<T, Tstiff, Tmass>::InitMatrix(K);
    Kh = &K;
    Mh = NULL;
  }
  
  
  //! Sets stiffness and mass matrix, eigenproblem is K x = lambda M x
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::InitMatrix(Matrix<Tstiff, Prop, Storage>& K,
	       Matrix<Tmass, PropM, StorageM>& M)
  {
    VirtualEigenProblem<T, Tstiff, Tmass>::InitMatrix(K, M);
    Kh = &K;
    Mh = &M;
  }
  

  //! computation of diagonal of mass matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::ComputeDiagonalMass()
  {
    Vector<Tmass>& D = this->sqrt_diagonal_mass;
    if (Mh == NULL)
      {
        // M = identity
        D.Reallocate(this->n_);
        D.Fill(1.0);
      }
    else
      {
        D.Reallocate(this->n_);
        for (int i = 0; i < this->n_; i++)
          D(i) = (*Mh)(i, i);
      }
  }

  
  //! Cholesky factorisation of mass matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  FactorizeCholeskyMass()
  {
    if (Mh == NULL)
      {
        mat_chol.Reallocate(this->n_, this->n_);
        mat_chol.SetIdentity();
      }
    else
      {
        mat_chol = *(Mh);
        GetCholesky(mat_chol);
        Xchol_real.Reallocate(mat_chol.GetM());
        Xchol_imag.Reallocate(mat_chol.GetM());
      }
  }
  
  
  //! computation of L X or L^T X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  MltCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X)
  {
    MltCholesky(transA, mat_chol, X);
  }


  //! computation of L X or L^T X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  MltCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    MltCholesky(transA, mat_chol, Xchol_real);
    MltCholesky(transA, mat_chol, Xchol_imag);

    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<Treal>(Xchol_real(i), Xchol_imag(i));    
  }
  
  
  //! computation of L^-1 X or L^-T X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  SolveCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X)
  {
    SolveCholesky(transA, mat_chol, X);
  }


  //! computation of L^-1 X or L^-T X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  SolveCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    SolveCholesky(transA, mat_chol, Xchol_real);
    SolveCholesky(transA, mat_chol, Xchol_imag);

    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<Treal>(Xchol_real(i), Xchol_imag(i));
  }
  

  //! matrix-vector product Y = M X where M is the mass matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltMass(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mlt(*Mh, X, Y);
  }


  //! matrix-vector product Y = M X where M is the mass matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mlt(*Mh, X, Y);
  }


  //! matrix-vector product Y = M^T X where M is the mass matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltMass(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mlt(trans, *Mh, X, Y);
  }


  //! matrix-vector product Y = M^T X where M is the mass matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltMass(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mlt(trans, *Mh, X, Y);
  }
    

  //! matrix-vector product Y = K X where K is the stiffness matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Mlt(*Kh, X, Y);
  }


  //! matrix-vector product Y = K X where K is the stiffness matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Mlt(*Kh, X, Y);
  }
  

  //! matrix-vector product Y = (coef_mas M + coef_stiff K) X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltStiffness(const T& coef_massb, const T& coef_stiffb,
		 const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Mlt(*Kh, X, Y);
    
    Treal coef_mass, coef_stiff;
    SetComplexReal(coef_massb, coef_mass);
    SetComplexReal(coef_stiffb, coef_stiff);
    if (coef_mass != T(0))
      {
        if (Mh == NULL)
          for (int i = 0; i < Y.GetM(); i++)
            Y(i) += coef_mass*X(i);
        else
          MltAdd(coef_mass, *Mh, X, coef_stiff, Y);
      }
    else
      {
        if (coef_stiff != T(1))
          Mlt(coef_stiff, Y);
      }
  }


  //! matrix-vector product Y = (coef_mas M + coef_stiff K) X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltStiffness(const T& coef_mass, const T& coef_stiff,
		 const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Mlt(*Kh, X, Y);
    
    if (coef_mass != T(0))
      {
        if (Mh == NULL)
          for (int i = 0; i < Y.GetM(); i++)
            Y(i) += coef_mass*X(i);
        else
          MltAdd(coef_mass, *Mh, X, coef_stiff, Y);
      }
    else
      {
        if (coef_stiff != T(1))
          Mlt(coef_stiff, Y);
      }
  }


  //! matrix-vector product Y = K^T X where K is the stiffness matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltStiffness(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Mlt(SeldonTrans, *Kh, X, Y);
  }


  //! matrix-vector product Y = K^T X where K is the stiffness matrix
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::MltStiffness(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      Mlt(SeldonTrans, *Kh, X, Y);
  }
  
  
  //! computation and factorisation of a M + b K
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b,
				       int which_part)
  {
    ComputeAndFactoRealMatrix(T(0), a, b, which_part);
  }


  //! computation and factorisation of a M + b K
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::ComputeAndFactoRealMatrix(const Tcplx&, const Treal& a, const Treal& b, int which)
  {
    cout << "Provide coefficients a and b of the same type as T" << endl;
    abort();
  }
  

  //! computation and factorisation of a M + b K
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>
  ::ComputeAndFactoRealMatrix(const Treal&, const Treal& a,
			      const Treal& b, int which_part)  
  {
    this->selected_part = which_part;
    if (Kh == NULL)
      this->PrintErrorInit();
    
    this->complex_system = false;
    // computation of mat_lu = a M + b K
    Copy(*Kh, mat_lu_real);
    Mlt(b, mat_lu_real);
    if (Mh == NULL)
      {
        for (int i = 0; i < this->n_; i++)
          mat_lu_real(i, i) += a;
      }
    else
      Add(a, *Mh, mat_lu_real);
    
    // factorisation
    GetLU(mat_lu_real, pivot);
  }
  
  
  //! computation and factorisation of a M + b K
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b,
                                     int which_part)
  {
    this->selected_part = which_part;
    if (Kh == NULL)
      this->PrintErrorInit();
 
    this->complex_system = true;
    // inverse of (a M + b K), then we take real_part or imaginary part
    Matrix<Tcplx, Prop, Storage> InvMat(this->n_, this->n_);
    for (int i = 0; i < this->n_; i++)
      for (int j = 0; j < this->n_; j++)
        InvMat(i, j) = (*Kh)(i, j);
    
    Mlt(b, InvMat);
    if (Mh == NULL)
      {
        for (int i = 0; i < this->n_; i++)
          InvMat(i, i) += a;
      }
    else
      {
        for (int i = 0; i < this->n_; i++)
          for (int j = 0; j < this->n_; j++)
            InvMat(i, j) += a * (*this->Mh)(i, j);
      }
    
    if (which_part == EigenProblem_Base<T>::COMPLEX_PART)
      {
	mat_lu_cplx = InvMat;
	InvMat.Clear();
	GetLU(mat_lu_cplx, pivot);
      }
    else
      {
	// inverse
	GetInverse(InvMat);
  
	// then extracting real or imaginary part
	mat_lu_real.Reallocate(this->n_, this->n_);
	if (which_part == EigenProblem_Base<T>::REAL_PART)
	  {
	    for (int i = 0; i < this->n_; i++)
	      for (int j = 0; j < this->n_; j++)
		mat_lu_real(i, j) = real(InvMat(i, j));
	  }
	else
	  {
	    for (int i = 0; i < this->n_; i++)
	      for (int j = 0; j < this->n_; j++)
		mat_lu_real(i, j) = imag(InvMat(i, j));
	  }
      }
  }
  

  //! solution of (a M + b K) Y = X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (this->complex_system)
      {
	if (this->selected_part == EigenProblem_Base<T>::COMPLEX_PART)
	  {
	    cout << "The result can not be a real vector" << endl;
	    abort();
	  }

	Mlt(mat_lu_real, X, Y);
      }
    else
      {
        Copy(X, Y);
        SolveLU(mat_lu_real, pivot, Y);
      }
  }


  //! solution of (a M + b K) Y = X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (this->complex_system)
      {
	if (this->selected_part == EigenProblem_Base<T>::COMPLEX_PART)
	  {
	    Copy(X, Y);
	    SolveLU(mat_lu_cplx, pivot, Y);
	  }
	else
	  Mlt(mat_lu_real, X, Y);
      }
    else
      {
        Copy(X, Y);
        SolveLU(mat_lu_real, pivot, Y);
      }
  }
   

  //! solution of (a M + b K) Y = X or transpose system
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeSolution(const SeldonTranspose& transA,
		  const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (this->complex_system)
      {
	if (this->selected_part == EigenProblem_Base<T>::COMPLEX_PART)
	  {
	    cout << "The result can not be a real vector" << endl;
	    abort();
	  }

	Mlt(transA, mat_lu_real, X, Y);
      }
    else
      {
        Copy(X, Y);
        SolveLU(transA, mat_lu_real, pivot, Y);
      }
  }
  

  //! solution of (a M + b K) Y = X
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeSolution(const SeldonTranspose& transA,
		  const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (this->complex_system)
      {
	if (this->selected_part == EigenProblem_Base<T>::COMPLEX_PART)
	  {
	    Copy(X, Y);
	    SolveLU(transA, mat_lu_cplx, pivot, Y);
	  }
	else
	  Mlt(transA, mat_lu_real, X, Y);
      }
    else
      {
        Copy(X, Y);
        SolveLU(transA, mat_lu_real, pivot, Y);
      }
  }

  
  //! clearing variables used for eigenvalue resolution
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Tstiff, Prop, Storage, Tmass, PropM, StorageM>::
  Clear()
  {
    EigenProblem_Base<T>::Clear();
    
    mat_lu_real.Clear();
    mat_lu_cplx.Clear();
    mat_chol.Clear();
  }
  
  
  /**********************
   * SparseEigenProblem *
   **********************/
  
  
  //! default constructor
  template<class T, class MatStiff, class MatMass>
  SparseEigenProblem<T, MatStiff, MatMass>::
  SparseEigenProblem()
    : VirtualEigenProblem<T, typename MatStiff::entry_type,
			  typename MatMass::entry_type>()
  {
    mat_lu_real.RefineSolution();
    mat_lu_cplx.RefineSolution();
    Mh = NULL;
    Kh = NULL;
    ProcSharingRows = NULL;
    SharingRowNumbers = NULL;
    nodl_scalar_ = nb_unknowns_scal_ = 0;
    nloc = 0;
  }
  
  
  //! sets Cholesky solver to use
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::SelectCholeskySolver(int type)
  {
    chol_facto_mass_matrix.SelectDirectSolver(type);
  }


  template<class T, class MatStiff, class MatMass>
  int SparseEigenProblem<T, MatStiff, MatMass>::RetrieveLocalNumbers(MatStiff& K)
  {
    nloc = K.GetM();
    
#ifdef SELDON_WITH_MPI
    try
      {
	DistributedMatrix_Base<typename MatStiff::entry_type>& A
	  = dynamic_cast<DistributedMatrix_Base<typename MatStiff::entry_type>& >(K);

	MPI_Comm comm = A.GetCommunicator();
	this->SetCommunicator(comm);
	int nb_proc;
	MPI_Comm_size(comm, &nb_proc);

	// only one processor => sequential case
	if (nb_proc <= 1)
	  return -1;

	// parallel case
	int m = A.GetLocalM();
	const IVect& OverlapRow = A.GetOverlapRowNumber();
	int noverlap = OverlapRow.GetM();
	int n = m - noverlap;
	local_col_numbers.Reallocate(n);
	Vector<bool> OverlappedRow(m); OverlappedRow.Fill(false);
	for (int i = 0; i < noverlap; i++)
	  OverlappedRow(OverlapRow(i)) = true;
	
	int ncol = 0;
	for (int i = 0; i < m; i++)
	  if (!OverlappedRow(i))
	    local_col_numbers(ncol++) = i;

	ProcSharingRows = &A.GetProcessorSharingRows();
	SharingRowNumbers = &A.GetSharingRowNumbers();
	nodl_scalar_ = A.GetNodlScalar();
	nb_unknowns_scal_ = A.GetNbScalarUnknowns();
	
	return n;
      }
    catch (const std::bad_cast&)
      {
	// a sequential matrix has been provided
	this->SetCommunicator(MPI_COMM_SELF);

	local_col_numbers.Clear();
      }
#endif

    return -1;
  }

  
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::InitMatrix(MatStiff& K)
  {
    int n = RetrieveLocalNumbers(K);
    distributed = false;
    if (n >= 0)
      distributed = true;
    
    VirtualEigenProblem<T, typename MatStiff::entry_type,
			typename MatMass::entry_type>::InitMatrix(K, n);
    Kh = &K;
    Mh = NULL;
  }
  
  
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::InitMatrix(MatStiff& K, MatMass& M)
  {
    int n = RetrieveLocalNumbers(K);
    distributed = false;
    if (n >= 0)
      distributed = true;
    
    VirtualEigenProblem<T, typename MatStiff::entry_type,
			typename MatMass::entry_type>::InitMatrix(K, M, n);
    
    Kh = &K;
    Mh = &M;
  }

  
  //! computation of diagonal of mass matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::ComputeDiagonalMass()
  {
    Vector<typename MatMass::entry_type>& D = this->sqrt_diagonal_mass;
    if (Mh == NULL)
      {
        // M = identity
        D.Reallocate(this->n_);
        D.Fill(1.0);
      }
    else
      {
        D.Reallocate(nloc);
        for (int i = 0; i < nloc; i++)
          D(i) = (*Mh)(i, i);
	
#ifdef SELDON_WITH_MPI
	// D is assembled for distributed matrices
	if (distributed)
	  {
	    Vector<typename MatMass::entry_type> M(D);
	    AssembleVector(M, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 15);

	    D.Reallocate(local_col_numbers.GetM());
	    for (int i = 0; i < local_col_numbers.GetM(); i++)
	      D(i) = M(local_col_numbers(i));
	  }
#endif
      }
  }


  //! computes Cholesky factorisation of M from matrix M
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::FactorizeCholeskyMass()
  {
    if (Mh == NULL)
      this->PrintErrorInit();
    
    if (this->print_level > 2)
      chol_facto_mass_matrix.ShowMessages();
    
    FactorizeCholeskyMass(*Mh);
    
    if (this->print_level < 4)
      chol_facto_mass_matrix.HideMessages();
    
    Xchol_real.Reallocate(nloc);
    Xchol_imag.Reallocate(nloc);
  }
  
  
  //! intermediary function
  template<class T, class MatStiff, class MatMass>
  template<class Storage>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  FactorizeCholeskyMass(Matrix<Treal, Symmetric, Storage>& M)
  {
    chol_facto_mass_matrix.Factorize(M, true);    
  }

  
  //! intermediary function
  template<class T, class MatStiff, class MatMass>
  template<class T0, class Prop, class Storage>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  FactorizeCholeskyMass(Matrix<T0, Prop, Storage>& M)
  {
    cout << "Cholesky factorisation has not been implemented "
	 << "for complex matrices" << endl;
    abort();    
  }
  

#ifdef SELDON_WITH_MPI
  //! intermediary function
  template<class T, class MatStiff, class MatMass>
  template<class Storage>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  FactorizeCholeskyMass(DistributedMatrix<Treal, Symmetric, Storage>& M)
  {
    chol_facto_mass_matrix.Factorize(M, true);
  }

  
  //! intermediary function
  template<class T, class MatStiff, class MatMass>
  template<class T0, class Prop, class Storage>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  FactorizeCholeskyMass(DistributedMatrix<T0, Prop, Storage>& M)
  {
    cout << "Cholesky factorisation has not been implemented "
	 << "for complex matrices" << endl;
    abort();    
  }
#endif
  
  
  //! computes L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  MltCholeskyMass(const SeldonTranspose& TransA, Vector<Treal>& X)
  {
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
	Xchol_real.Reallocate(nloc);
	Xchol_real.Zero();
	for (int i = 0; i < this->n_; i++)
	  Xchol_real(local_col_numbers(i)) = X(i);
	
	// actually assemble is not needed since the direct solver
	// should extract only local values
	//AssembleVector(Xchol_real, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
	// this->comm, nodl_scalar_, nb_unknowns_scal_, 16);	    
#endif
      }
    else
      Xchol_real = X;

    chol_facto_mass_matrix.Mlt(TransA, Xchol_real, false);

    if (distributed)
      {
#ifdef SELDON_WITH_MPI
	for (int i = 0; i < this->n_; i++)
	  X(i) = Xchol_real(local_col_numbers(i));
#endif
      }
    else	      
      Copy(Xchol_real, X);	        
  }
  
  
  //! computes L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  SolveCholeskyMass(const SeldonTranspose& TransA, Vector<Treal>& X)
  {
     if (distributed)
      {
#ifdef SELDON_WITH_MPI
	Xchol_real.Reallocate(nloc);
	Xchol_real.Zero();
	for (int i = 0; i < this->n_; i++)
	  Xchol_real(local_col_numbers(i)) = X(i);
	
	// actually assemble is not needed since the direct solver
	// should extract only local values
	//AssembleVector(Xchol_real, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
	// this->comm, nodl_scalar_, nb_unknowns_scal_, 16);	    
#endif
      }
     else
       Xchol_real = X;
    
    chol_facto_mass_matrix.Solve(TransA, Xchol_real, false);
    
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
	for (int i = 0; i < this->n_; i++)
	  X(i) = Xchol_real(local_col_numbers(i));
#endif
      }
    else	      
      Copy(Xchol_real, X);	        
  }

  
  //! computes L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  MltCholeskyMass(const SeldonTranspose& TransA, Vector<Tcplx>& X)  
  {
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
	Xchol_real.Reallocate(nloc);
	Xchol_imag.Reallocate(nloc);
	for (int i = 0; i < this->n_; i++)
	  {
	    Xchol_real(local_col_numbers(i)) = real(X(i));
	    Xchol_imag(local_col_numbers(i)) = imag(X(i));
	  }
	
	// actually assemble is not needed since the direct solver
	// should extract only local values
	//AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
	// this->comm, nodl_scalar_, nb_unknowns_scal_, 16);	    
#endif
      }
    else
      for (int i = 0; i < X.GetM(); i++)
	{
	  Xchol_real(i) = real(X(i));
	  Xchol_imag(i) = imag(X(i));
	}
    
    chol_facto_mass_matrix.Mlt(TransA, Xchol_real, false);
    chol_facto_mass_matrix.Mlt(TransA, Xchol_imag, false);

    if (distributed)
      {
#ifdef SELDON_WITH_MPI
	for (int i = 0; i < this->n_; i++)
	  X(i) = complex<Treal>(Xchol_real(local_col_numbers(i)), Xchol_imag(local_col_numbers(i)));
#endif
      }
    else	      
      for (int i = 0; i < X.GetM(); i++)
	X(i) = complex<Treal>(Xchol_real(i), Xchol_imag(i));    
  }
  
  
  //! computes L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  SolveCholeskyMass(const SeldonTranspose& TransA, Vector<Tcplx>& X)
  {
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
	Xchol_real.Reallocate(nloc);
	Xchol_imag.Reallocate(nloc);
	for (int i = 0; i < this->n_; i++)
	  {
	    Xchol_real(local_col_numbers(i)) = real(X(i));
	    Xchol_imag(local_col_numbers(i)) = imag(X(i));
	  }
	
	// actually assemble is not needed since the direct solver
	// should extract only local values
	//AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
	// this->comm, nodl_scalar_, nb_unknowns_scal_, 16);	    
#endif
      }
    else
      for (int i = 0; i < X.GetM(); i++)
	{
	  Xchol_real(i) = real(X(i));
	  Xchol_imag(i) = imag(X(i));
	}

    chol_facto_mass_matrix.Solve(TransA, Xchol_real, false);
    chol_facto_mass_matrix.Solve(TransA, Xchol_imag, false);

    if (distributed)
      {
#ifdef SELDON_WITH_MPI
	for (int i = 0; i < this->n_; i++)
	  X(i) = complex<Treal>(Xchol_real(local_col_numbers(i)), Xchol_imag(local_col_numbers(i)));
#endif
      }
    else	      
      for (int i = 0; i < X.GetM(); i++)
	X(i) = complex<Treal>(Xchol_real(i), Xchol_imag(i));
  }


  //! multiplication by mass matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::MltMass(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Mh == NULL)
      Seldon::Copy(X, Y);
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Treal> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Mh->MltVector(Xcplx, Ycplx);

	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  Mh->MltVector(X, Y);
      }
  }


  //! multiplication by mass matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Mh == NULL)
      Seldon::Copy(X, Y);
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Tcplx> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Mh->MltVector(Xcplx, Ycplx);

	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  Mh->MltVector(X, Y);
      }
  }


  //! multiplication by mass matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::MltMass(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Mh == NULL)
      Seldon::Copy(X, Y);
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Treal> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Mh->MltVector(trans, Xcplx, Ycplx);

	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  Mh->MltVector(trans, X, Y);
      }
  }


  //! multiplication by mass matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::MltMass(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Mh == NULL)
      Seldon::Copy(X, Y);
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Tcplx> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Mh->MltVector(trans, Xcplx, Ycplx);

	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  Mh->MltVector(trans, X, Y);
      }
  }


  //! multiplication by stiffness matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Treal> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Kh->MltVector(Xcplx, Ycplx);

	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  Kh->MltVector(X, Y);
      }
  }


  //! multiplication by stiffness matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Tcplx> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Kh->MltVector(Xcplx, Ycplx);

	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  Kh->MltVector(X, Y);
      }
  }


  //! multiplication by stiffness and mass matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::MltStiffness(const T& coef_massb, const T& coef_stiffb,
		 const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      {
	if (distributed)
	  {

#ifdef SELDON_WITH_MPI
            Treal coef_mass, coef_stiff;
            SetComplexReal(coef_massb, coef_mass);
            SetComplexReal(coef_stiffb, coef_stiff);

	    Vector<Treal> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Kh->MltVector(Xcplx, Ycplx);
	    if (Mh == NULL)
	      for (int i = 0; i < nloc; i++)
		Ycplx(i) += coef_mass*Xcplx(i);
	    else
	      Mh->MltAddVector(coef_mass, Xcplx, coef_stiff, Ycplx);
	    
	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  VirtualEigenProblem<T, typename MatStiff::entry_type, typename MatMass::entry_type>
	    ::MltStiffness(coef_massb, coef_stiffb, X, Y);
      }
  }


  //! multiplication by stiffness and mass matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::MltStiffness(const T& coef_mass, const T& coef_stiff,
		 const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Tcplx> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Kh->MltVector(Xcplx, Ycplx);
	    if (Mh == NULL)
	      for (int i = 0; i < nloc; i++)
		Ycplx(i) += coef_mass*Xcplx(i);
	    else
	      Mh->MltAddVector(coef_mass, Xcplx, coef_stiff, Ycplx);
	    
	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  VirtualEigenProblem<T, typename MatStiff::entry_type, typename MatMass::entry_type>
	    ::MltStiffness(coef_mass, coef_stiff, X, Y);
      }
  }


  //! multiplication by stiffness matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::MltStiffness(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Treal> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Kh->MltVector(trans, Xcplx, Ycplx);

	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  Kh->MltVector(trans, X, Y);
      }
  }


  //! multiplication by stiffness matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::MltStiffness(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (Kh == NULL)
      this->PrintErrorInit();
    else
      {
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<Tcplx> Xcplx(nloc), Ycplx(nloc);
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = X(i);
	    
	    AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    Kh->MltVector(trans, Xcplx, Ycplx);

	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Ycplx(local_col_numbers(i)); 
#endif
	  }
	else
	  Kh->MltVector(trans, X, Y);
      }
  }

  
  //! computes and factorizes a M + b K
  //! M is the mass matrix and K the stiffness matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b, int which)
  {
    ComputeAndFactoRealMatrix(T(0), a, b, which);
  }


  //! computes and factorizes a M + b K
  //! where M is the mass matrix and K the stiffness matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactoRealMatrix(const Tcplx&, const Treal& a, const Treal& b, int which)
  {
    cout << "Provide coefficients a and b of the same type as T" << endl;
    abort();
  }

  
  //! computes and factorizes a M + b K
  //! where M is the mass matrix and K the stiffness matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactoRealMatrix(const Treal&, const Treal& a, const Treal& b, int which)
  {
    this->complex_system = false;
    
    if (Kh == NULL)
      this->PrintErrorInit();
    
    if (this->print_level > 2)
      mat_lu_real.ShowMessages();
    
    // retrieving symmetry of mass matrix 
    bool sym_mh = true;
    if (Mh != NULL)
      {
	if (!IsSymmetricMatrix(*Mh))
	  sym_mh = false;
      }
    
    T zero; SetComplexZero(zero);
    
    if (b == zero)
      {
	// only mass matrix must be factorized
	if (sym_mh)
	  {
#ifdef SELDON_WITH_MPI
	    DistributedMatrix<Treal, Symmetric, ArrayRowSymSparse> A;
	    try
	      {
		DistributedMatrix_Base<typename MatStiff::entry_type>& K
		  = dynamic_cast<DistributedMatrix_Base<typename MatStiff::entry_type>& >(*Kh);

		A.Init(K);
	      }
	    catch (const std::bad_cast&)
	      {
	      }
#else
	    Matrix<Treal, Symmetric, ArrayRowSymSparse> A;
#endif
	    if (Mh == NULL)
	      {
		A.Reallocate(nloc, nloc);
		A.SetIdentity();
	      }
	    else
	      Copy(*Mh, A);
	    
	    Mlt(a, A);
	    mat_lu_real.Factorize(A);
	  }
	else
	  {
#ifdef SELDON_WITH_MPI
	    DistributedMatrix<Treal, General, ArrayRowSparse> A;
#else
	    Matrix<Treal, General, ArrayRowSparse> A;
#endif
	    Copy(*Mh, A);
	    Mlt(a, A);
	    mat_lu_real.Factorize(A);
	  }
      }
    else if (IsSymmetricMatrix(*Kh) && sym_mh)
      {
	// forming a M + b K and factorizing it when the result is symmetric
#ifdef SELDON_WITH_MPI
        DistributedMatrix<Treal, Symmetric, ArrayRowSymSparse> A;
#else
	Matrix<Treal, Symmetric, ArrayRowSymSparse> A;
#endif
	
        Copy(*Kh, A);
        Mlt(b, A);
	if (a != zero)
	  {
	    if (Mh == NULL)
	      {
#ifdef SELDON_WITH_MPI
		if (distributed)
		  {
		    for (int i = 0; i < this->local_col_numbers.GetM(); i++)
		      {
			int iloc = this->local_col_numbers(i);
			A.AddInteraction(iloc, iloc, a);
		      }
		  }
		else
		  for (int i = 0; i < nloc; i++)
		    A.AddInteraction(i, i, a);
#else
		for (int i = 0; i < nloc; i++)
		  A.AddInteraction(i, i, a);
#endif
	      }
	    else
	      Add(a, *Mh, A);
	  }
	
        mat_lu_real.Factorize(A);
      } 
    else
      {
	// forming a M + b K and factorizing it when the result is unsymmetric
#ifdef SELDON_WITH_MPI
        DistributedMatrix<Treal, General, ArrayRowSparse> A;
#else
	Matrix<Treal, General, ArrayRowSparse> A;
#endif
        Copy(*Kh, A);
        Mlt(b, A);
	if (a != zero)
	  {
	    if (Mh == NULL)
	      {
#ifdef SELDON_WITH_MPI
		if (distributed)
		  {
		    for (int i = 0; i < this->local_col_numbers.GetM(); i++)
		      {
			int iloc = this->local_col_numbers(i);
			A.AddInteraction(iloc, iloc, a);
		      }
		  }
		else
		  for (int i = 0; i < nloc; i++)
		    A.AddInteraction(i, i, a);
#else
		for (int i = 0; i < nloc; i++)
		  A.AddInteraction(i, i, a);
#endif
	      }
	    else
	      Add(a, *Mh, A);
	  }
	
        mat_lu_real.Factorize(A);
      }      

    if (this->print_level < 4)
      mat_lu_real.HideMessages();
  }
  
  
  //! computes and factorizes matrix (a M + b K) for complex values of a and b
  template<class T, class MatStiff, class MatMass> 
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b, int which)
  {
    this->complex_system = true;
    this->selected_part = which;
    
    if (Kh == NULL)
      this->PrintErrorInit();
    
    if (this->print_level > 2)
      mat_lu_cplx.ShowMessages();
    
    // retrieving symmetry of mass matrix 
    bool sym_mh = true;
    if (Mh != NULL)
      {
	if (!IsSymmetricMatrix(*Mh))
	  sym_mh = false;
      }

    Tcplx zero(0, 0);
    
    if (b == zero)
      {
	// only mass matrix must be factorized
	if (sym_mh)
	  {
#ifdef SELDON_WITH_MPI
	    DistributedMatrix<Tcplx, Symmetric, ArrayRowSymSparse> A;
	    try
	      {
		DistributedMatrix_Base<typename MatStiff::entry_type>& K
		  = dynamic_cast<DistributedMatrix_Base<typename MatStiff::entry_type>& >(*Kh);

		A.Init(K);
	      }
	    catch (const std::bad_cast&)
	      {
	      }
#else
	    Matrix<Tcplx, Symmetric, ArrayRowSymSparse> A;
#endif
	    if (Mh == NULL)
	      {
		A.Reallocate(nloc, nloc);
		A.SetIdentity();
	      }
	    else
	      Copy(*Mh, A);
	    
	    Mlt(a, A);
	    mat_lu_cplx.Factorize(A);
	  }
	else
	  {
#ifdef SELDON_WITH_MPI
	    DistributedMatrix<Tcplx, General, ArrayRowSparse> A;
#else
	    Matrix<Tcplx, General, ArrayRowSparse> A;
#endif
	    Copy(*Mh, A);
	    Mlt(a, A);
	    mat_lu_cplx.Factorize(A);
	  }
      }
    else if (IsSymmetricMatrix(*Kh) && sym_mh)
      {
	// forming a M + b K
#ifdef SELDON_WITH_MPI
	DistributedMatrix<Tcplx, Symmetric, ArrayRowSymSparse> A;
#else
	Matrix<Tcplx, Symmetric, ArrayRowSymSparse> A;
#endif
	Copy(*Kh, A);
	Mlt(b, A);
	if (a != zero)
	  {
	    if (Mh == NULL)
	      {
#ifdef SELDON_WITH_MPI
		if (distributed)
		  {
		    for (int i = 0; i < this->local_col_numbers.GetM(); i++)
		      {
			int iloc = this->local_col_numbers(i);
			A.AddInteraction(iloc, iloc, a);
		      }
		  }
		else
		  for (int i = 0; i < nloc; i++)
		    A.AddInteraction(i, i, a);
#else
		for (int i = 0; i < nloc; i++)
		  A.AddInteraction(i, i, a);
#endif
	      }
	    else
	      Add(a, *Mh, A);
	  }
	
	mat_lu_cplx.Factorize(A);
      }
    else
      {
	// forming a M + b K
#ifdef SELDON_WITH_MPI
	DistributedMatrix<Tcplx, General, ArrayRowSparse> A;
#else
	Matrix<Tcplx, General, ArrayRowSparse> A;
#endif
	Copy(*Kh, A);
	Mlt(b, A);
	if (a != zero)
	  {
	    if (Mh == NULL)
	      {
#ifdef SELDON_WITH_MPI
		if (distributed)
		  {
		    for (int i = 0; i < this->local_col_numbers.GetM(); i++)
		      {
			int iloc = this->local_col_numbers(i);
			A.AddInteraction(iloc, iloc, a);
		      }
		  }
		else
		  for (int i = 0; i < nloc; i++)
		    A.AddInteraction(i, i, a);
#else
		for (int i = 0; i < nloc; i++)
		  A.AddInteraction(i, i, a);
#endif
	      }
	    else
	      Add(a, *Mh, A);
	  }
	
	mat_lu_cplx.Factorize(A);
      }

    if (this->print_level < 4)
      mat_lu_cplx.HideMessages();
  }


  //! solves (a M + b K) Y = X with stored factorization 
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    ComputeSolution(SeldonNoTrans, X, Y);
  }


  //! solves (a M + b K) Y = X with stored factorization 
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    ComputeSolution(SeldonNoTrans, X, Y);
  }

  
  //! solves (a M + b K) Y = X with stored factorization 
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeSolution(const SeldonTranspose& transA, 
		  const Vector<Treal>& X, Vector<Treal>& Y)
  {
    if (this->complex_system)
      {
	if (this->selected_part == EigenProblem_Base<T>::COMPLEX_PART)
	  {
	    cout << "The result can not be a real vector" << endl;
	    abort();
	  }
	
	Vector<Tcplx> Xcplx(nloc);
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Xcplx.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xcplx(local_col_numbers(i)) = Tcplx(X(i), 0);

	    // actually assemble is not needed since the direct solver
	    // should extract only local values
	    //AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
	    // this->comm, nodl_scalar_, nb_unknowns_scal_, 16);	    
#endif
	  }
	else
	  for (int i = 0; i < nloc; i++)
	    Xcplx(i) = Tcplx(X(i), 0);
	
	mat_lu_cplx.Solve(transA, Xcplx, false);
	
	if (this->selected_part == EigenProblem_Base<T>::IMAG_PART)
	  {
	    if (distributed)
	      {
#ifdef SELDON_WITH_MPI
		for (int i = 0; i < this->n_; i++)
		  Y(i) = imag(Xcplx(local_col_numbers(i)));
#endif
	      }
	    else
	      for (int i = 0; i < nloc; i++)
		Y(i) = imag(Xcplx(i));
	  }
	else
	  {
	    if (distributed)
	      {
#ifdef SELDON_WITH_MPI
		for (int i = 0; i < this->n_; i++)
		  Y(i) = real(Xcplx(local_col_numbers(i)));
#endif
	      }
	    else	      
	      for (int i = 0; i < this->n_; i++)
		Y(i) = real(Xcplx(i));
	  }
      }
    else
      {
	Vector<Treal> Xreal(nloc);
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Xreal.Zero();
	    for (int i = 0; i < this->n_; i++)
	      Xreal(local_col_numbers(i)) = X(i);

	    // actually assemble is not needed since the direct solver
	    // should extract only local values
	    // AssembleVector(Xreal, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
	    // this->comm, nodl_scalar_, nb_unknowns_scal_, 16);	    
#endif
	  }
	else
	  Xreal = X;
	
        mat_lu_real.Solve(transA, Xreal, false);

	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    for (int i = 0; i < this->n_; i++)
	      Y(i) = Xreal(local_col_numbers(i));
#endif
	  }
	else	      
	  Copy(Xreal, Y);
      }
  }
  
  
  //! solves (a M + b K) Y = X or transpose system
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::ComputeSolution(const SeldonTranspose& transA,
		    const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    if (this->complex_system)
      {
	if (this->selected_part == EigenProblem_Base<T>::COMPLEX_PART)
	  {
	    Vector<Tcplx> Xcplx(nloc);
	    if (distributed)
	      {
#ifdef SELDON_WITH_MPI
		Xcplx.Zero();
		for (int i = 0; i < this->n_; i++)
		  Xcplx(local_col_numbers(i)) = X(i);

		// actually assemble is not needed since the direct solver
		// should extract only original values
		//AssembleVector(Xcplx, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
		//	    this->comm, nodl_scalar_, nb_unknowns_scal_, 16);	    
#endif
	      }
	    else
	      Xcplx = X;
	    
	    mat_lu_cplx.Solve(transA, Xcplx, false);

	    if (distributed)
	      {
#ifdef SELDON_WITH_MPI
		for (int i = 0; i < this->n_; i++)
		  Y(i) = Xcplx(local_col_numbers(i));
#endif
	      }
	    else	      
	      Copy(Xcplx, Y);	    
	  }
	else
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
      }
    else
      {
	cout << "not implemented" << endl;
	abort();
        //Copy(X, Y);
        //mat_lu_real.Solve(transA, Y);
      }
  }
  

  //! changes final eigenvectors if needed
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>
  ::DistributeEigenvectors(Matrix<T, General, ColMajor>& eigen_vec)
  {
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
	int nvec = eigen_vec.GetN();
	eigen_vec.Resize(nloc, nvec);

	Vector<T> X(nloc);
	for (int j = 0; j < nvec; j++)
	  {
	    X.Zero();
	    for (int i = 0; i < this->n_; i++)
	      X(local_col_numbers(i)) = eigen_vec(i, j);
	    
	    AssembleVector(X, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
	    
	    for (int i = 0; i < nloc; i++)
	      eigen_vec(i, j) = X(i);
	  }
#endif
      }
  }

  //! clears memory used by the object
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::Clear()
  {
    mat_lu_real.Clear();
    mat_lu_cplx.Clear();
    chol_facto_mass_matrix.Clear();
    Xchol_real.Clear();
    Xchol_imag.Clear();    
  }

  
  int TypeEigenvalueSolver::GetDefaultSolver()
  {
#ifdef SELDON_WITH_ARPACK
    return ARPACK;
#endif
    
#ifdef SELDON_WITH_ANASAZI
    return ANASAZI;
#endif
    
#ifdef SELDON_WITH_FEAST
    return FEAST;
#endif

#ifdef SELDON_WITH_SLEPC
    return SLEPC;
#endif

    return -1;
  }
  

  template<class T, class Prop, class Storage>
  void GetEigenvaluesEigenvectors(EigenProblem_Base<T>& var_eig,
				  Vector<T>& lambda, Vector<T>& lambda_imag,
				  Matrix<T, Prop, Storage>& eigen_vec,
				  int type_solver)
  {
    if (type_solver == TypeEigenvalueSolver::DEFAULT)
      type_solver = TypeEigenvalueSolver::GetDefaultSolver();

    if (type_solver == TypeEigenvalueSolver::ARPACK)
      {
#ifdef SELDON_WITH_ARPACK
        T zero; SetComplexZero(zero);
        Matrix<T, General, ColMajor> eigen_old;
        FindEigenvaluesArpack(var_eig, lambda, lambda_imag, eigen_old);
        
        // eigenvalues are sorted by ascending order
        SortEigenvalues(lambda, lambda_imag, eigen_old,
                        eigen_vec, var_eig.LARGE_EIGENVALUES,
                        var_eig.GetTypeSorting(), zero, zero);
#else
        cout << "Recompile with Arpack" << endl;
        abort();
#endif
      }
    else if (type_solver == TypeEigenvalueSolver::ANASAZI)
      {
#ifdef SELDON_WITH_ANASAZI
	T zero; SetComplexZero(zero);
        Matrix<T, General, ColMajor> eigen_old;
        FindEigenvaluesAnasazi(var_eig, lambda, lambda_imag, eigen_old);
        
        // eigenvalues are sorted by ascending order
        SortEigenvalues(lambda, lambda_imag, eigen_old,
                        eigen_vec, var_eig.LARGE_EIGENVALUES,
                        var_eig.GetTypeSorting(), zero, zero);
#else
        cout << "Recompile with Anasazi" << endl;
        abort();
#endif
      }
    else if (type_solver == TypeEigenvalueSolver::FEAST)
      {
#ifdef SELDON_WITH_FEAST
        T zero; SetComplexZero(zero);
        Matrix<T, General, ColMajor> eigen_old;
        FindEigenvaluesFeast(var_eig, lambda, lambda_imag, eigen_old);
        
        // eigenvalues are sorted by ascending order
        SortEigenvalues(lambda, lambda_imag, eigen_old,
                        eigen_vec, var_eig.LARGE_EIGENVALUES,
                        var_eig.GetTypeSorting(), zero, zero);
#else
        cout << "Recompile with MKL or Feast" << endl;
        abort();
#endif
      }
    else if (type_solver == TypeEigenvalueSolver::SLEPC)
      {
#ifdef SELDON_WITH_SLEPC
	T zero; SetComplexZero(zero);
        Matrix<T, General, ColMajor> eigen_old;
        FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_old);
        
        // eigenvalues are sorted by ascending order
        SortEigenvalues(lambda, lambda_imag, eigen_old,
                        eigen_vec, var_eig.LARGE_EIGENVALUES,
                        var_eig.GetTypeSorting(), zero, zero);
#else
        cout << "Recompile with Slepc" << endl;
        abort();
#endif
      }
    else
      {
        cout << "Recompile with eigenvalue solver" << endl;
        abort();
      }
    
  }
}

#define SELDON_FILE_VIRTUAL_EIGENVALUE_SOLVER_CXX
#endif
