#ifndef SELDON_FILE_EIGENVALUE_SOLVER_CXX

#include "EigenvalueSolver.hxx"

namespace Seldon
{
  
  /******************
   * Initialization *
   ******************/
  
  
  //! default constructor
  template<class T, class MatStiff, class MatMass>
  EigenProblem_Base<T, MatStiff, MatMass>::EigenProblem_Base()
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
    Mh = NULL;
    Kh = NULL;

#ifdef SELDON_WITH_MPI
    // for parallel execution, default communicator : all the processors
    comm = MPI_COMM_WORLD;
#endif
    
    type_solver = SOLVER_LOBPCG;
    ortho_manager = ORTHO_DGKS;
    nb_blocks = 2;
    restart_number = 20;
  }
  
  
  //! initialisation of the size of the eigenvalue problem
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::Init(int n)
  {
    n_ = n;
    nb_prod = 0;

    // counting the size of the global system for parallel computation
    int nglob = n;
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
  
  
  //! initialization of a standard eigenvalue problem
  /*!
    Stiffness matrix K is given in argument.
    we will search (lambda, x) such as K x = lambda x
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  InitMatrix(MatStiff& K)
  {
    Kh = &K;
    Mh = NULL;
    this->diagonal_mass = true;
    if ( (!IsSymmetricMatrix(K))
         && (!IsComplexMatrix(K)) && (shift_imag != T(0)) )
      {
        // for real unsymmetric problems, if sigma is complex
        // we have to use mode 3 or 4 in Arpack => generalized problem
        this->diagonal_mass = false;
      }
    
    this->Init(K.GetM());
  }
  
  
  //! initialization of a generalized eigenvalue problem
  /*!
    Mass matrix M and stiffness matrix K are given in argument
    we will search (lambda, x) such as K x = lambda M x
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  InitMatrix(MatStiff& K, MatMass& M)
  {
    Kh = &K;
    Mh = &M;
    this->diagonal_mass = false;
    this->Init(K.GetM());
  }
  
  
  /*******************
   * Basic functions *
   *******************/
  
  
  //! returns the spectral transformation used for evaluation of eigenvalues
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetComputationalMode() const
  {
    return eigenvalue_computation_mode;
  }
  
  
  //! sets the spectral transformation used for evaluation of eigenvalues
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetComputationalMode(int mode)
  {
    eigenvalue_computation_mode = mode;
  }
  
  
  //! returns the number of eigenvalues asked by the user
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetNbAskedEigenvalues() const
  {
    return nb_eigenvalues_wanted;
  }


  //! returns the additional number of eigenvalues
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetNbAdditionalEigenvalues() const
  {
    return nb_add_eigenvalues;
  }


  //! returns the number of blocks used in blocked solvers
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetNbBlocks() const
  {
    return nb_blocks;
  }


  //! returns the number of blocks used in blocked solvers
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbBlocks(int n)
  {
    nb_blocks = n;
  }
  
  
  //! returns the restart parameter used in blocked solvers
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetNbMaximumRestarts() const
  {
    return restart_number;
  }
  

  //! sets the restart parameter used in blocked solvers
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbMaximumRestarts(int m)
  {
    restart_number = m;
  }

  
  //! returns orthogonalization manager set in Anasazi
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetOrthoManager() const
  {
    return ortho_manager;
  }
  
  
  //! returns the solver used in Anasazi
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetEigensolverType() const
  {
    return type_solver;
  }


  //! sets the solver used in Anasazi
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetEigensolverType(int type)
  {
    type_solver = type;
  }
  

#ifdef SELDON_WITH_MPI
  //! returns the MPI communicator shared by processors
  template<class T, class MatStiff, class MatMass>
  MPI_Intracomm& EigenProblem_Base<T, MatStiff, MatMass>::GetCommunicator()
  {
    return comm;
  }
  
  //! sets the MPI communicator shared by processors
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetCommunicator(MPI_Comm& comm_)
  {
    comm = comm_;
  }
#endif
  
  //! sets the number of eigenvalues to compute
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbAskedEigenvalues(int n)
  {
    nb_eigenvalues_wanted = n;
  }


  //! sets the number of additional eigenvalues
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbAdditionalEigenvalues(int n)
  {
    nb_add_eigenvalues = n;
  }
  
  
  //! returns the spectrum desired (large, small eigenvalues, etc)
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetTypeSpectrum() const
  {
    return type_spectrum_wanted;
  }

  
  //! returns how eigenvalues are sorted (real, imaginary part or modulus)
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetTypeSorting() const
  {
    return type_sort_eigenvalues;
  }

  
  //! returns the shift value used
  /*!
    If type_spectrum_wanted is set to CENTERED_EIGENVALUES,
    we search closest eigenvalues to the shift value.
    Matrix (A - (shift + i shift_imag)*I)^{-1} will be used instead of A
  */
  template<class T, class MatStiff, class MatMass>
  T EigenProblem_Base<T, MatStiff, MatMass>::GetShiftValue() const
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
  template<class T, class MatStiff, class MatMass>
  T EigenProblem_Base<T, MatStiff, MatMass>::GetImagShiftValue() const
  {
    return shift_imag;
  }
  
  
  //! Sets the real part of shift value
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetShiftValue(const T& val)
  {
    shift = val;
  }

  
  //! Sets the imaginary part of shift value
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetImagShiftValue(const T& val)
  {
    shift_imag = val;
  }

  
  //! sets which eigenvalues are searched
  /*!
    You can ask small eigenvalues, large, or eigenvalues
    close to the shift.
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SetTypeSpectrum(int type, const T& val, int type_sort)
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
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SetTypeSpectrum(int type, const complex<T>& val, int type_sort)
  {
    // for real unsymmetric eigenproblems, you can
    // specify a complex shift
    type_spectrum_wanted = type;
    shift = real(val);
    shift_imag = imag(val);

    if (Kh != NULL)
      {
        if ( (!IsSymmetricMatrix(*Kh))
             && (!IsComplexMatrix(*Kh)) && (shift_imag != T(0)) )
          {
            // for real unsymmetric problems, if sigma is complex
            // we have to use mode 3 or 4 in Arpack => generalized problem
            this->diagonal_mass = false;
          }
      }
    
    type_sort_eigenvalues = type_sort;
  }
  
  
  //! returns lower bound of the interval where eigenvalues are searched
  template<class T, class MatStiff, class MatMass>
  double EigenProblem_Base<T, MatStiff, MatMass>
  ::GetLowerBoundInterval() const
  {
    return emin_interval;
  }


  //! returns upper bound of the interval where eigenvalues are searched
  template<class T, class MatStiff, class MatMass>
  double EigenProblem_Base<T, MatStiff, MatMass>
  ::GetUpperBoundInterval() const
  {
    return emax_interval;
  }
  
  
  //! sets the interval where eigenvalues are searched
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>
  ::SetIntervalSpectrum(double l0, double l1)
  {
    emin_interval = l0;
    emax_interval = l1;
  }
    

  //! indicates the use of Cholesky factorisation in order to 
  //! solve a standard eigenvalue problem instead of a generalized one
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SetCholeskyFactoForMass(bool chol)
  {
    use_cholesky = chol;
  }
  
  
  //! returns true if Cholesky factorisation has to be used for mass matrix
  template<class T, class MatStiff, class MatMass>
  bool EigenProblem_Base<T, MatStiff, MatMass>::
  UseCholeskyFactoForMass() const
  {
    return use_cholesky;
  }
    
  
  //! indicates that the mass matrix is diagonal
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetDiagonalMass(bool diag)
  {
    diagonal_mass = diag;
  }
  
  
  //! returns true if the mass matrix is diagonal
  template<class T, class MatStiff, class MatMass>
  bool EigenProblem_Base<T, MatStiff, MatMass>::DiagonalMass() const
  {
    return diagonal_mass;
  }
  
    
  //! modifies the stopping critertion
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SetStoppingCriterion(double eps)
  {
    stopping_criterion = eps;
  }
  
  
  //! returns the stopping criterion
  template<class T, class MatStiff, class MatMass>
  double EigenProblem_Base<T, MatStiff, MatMass>::
  GetStoppingCriterion() const
  {
    return stopping_criterion;
  }
    
  
  //! sets the maximal number of iterations allowed for the iterative algorithm
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbMaximumIterations(int n)
  {
    nb_maximum_iterations = n;
  }
  
  
  //! returns the maximal number of iterations allowed for the iterative algorithm
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::
  GetNbMaximumIterations() const
  {
    return nb_maximum_iterations;
  }
  
  
  //! returns the number of matrix-vector products performed 
  //! since last call to Init
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::
  GetNbMatrixVectorProducts() const
  {
    return nb_prod;
  }
    
  
  //! returns the number of Arnoldi vectors to use
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetNbArnoldiVectors() const
  {
    return nb_arnoldi_vectors;
  }
  
  
  //! sets the number of Arnoldi vectors to use
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbArnoldiVectors(int n)
  {
    automatic_selection_arnoldi_vectors = false;
    nb_arnoldi_vectors = n;
  }
  
  
  //! returns number of rows
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetM() const
  {
    return n_;
  }
    
  
  //! returns number of columns
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetN() const
  {
    return n_;
  }
  
  
  //! returns level of verbosity
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetPrintLevel() const
  {
    return print_level;
  }
  
  
  //! sets the level of verbosity
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetPrintLevel(int lvl)
  {
    print_level = lvl;
  }
  
  
  //! increment of the number of matrix vector products
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::IncrementProdMatVect()
  {
    nb_prod++;
    if (print_level >= 3)
      {
        if (nb_prod%10 == 0)
#ifdef SELDON_WITH_MPI
          if (comm.Get_rank() == 0)
#endif
            cout<<" Iteration number " << nb_prod << endl;
      }
    else if (print_level >= 1)
      {
        if (nb_prod%100 == 0)
#ifdef SELDON_WITH_MPI
          if (comm.Get_rank() == 0)
#endif
            cout<<" Iteration number " << nb_prod << endl;
      }			
  }
  

  //! prints error of initialization and aborts program
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::PrintErrorInit() const
  {
    cout << "InitMatrix has not been called" << endl;
    abort();
  }
  
  
  //! returns true if the matrix is symmetric
  template<class T, class MatStiff, class MatMass>
  bool EigenProblem_Base<T, MatStiff, MatMass>::IsSymmetricProblem() const
  {
    if (Kh != NULL)
      {
        if (IsSymmetricMatrix(*Kh))
          {
            if (Mh == NULL)
              return true;
            else
              return IsSymmetricMatrix(*Mh);
          }
      }
    else
      PrintErrorInit();
    
    return false;
  }


  //! returns true if the matrix is hermitian
  template<class T, class MatStiff, class MatMass>
  bool EigenProblem_Base<T, MatStiff, MatMass>::IsHermitianProblem() const
  {
    if (Kh != NULL)
      {
	if (IsComplexMatrix(*Kh))
	  return false;
	else
	  {
            if (IsSymmetricMatrix(*Kh))
              {
                if (Mh == NULL)
                  return true;
                else
                  {
		    if (IsComplexMatrix(*Mh))
		      return false;

		    return IsSymmetricMatrix(*Mh);
		  }
              }
          }
      }
    else
      PrintErrorInit();
    
    return false;
  }
  
    
  /*********************
   * Mass matrix stuff *
   *********************/
  
  
  //! computation of D^1/2 from D
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::FactorizeDiagonalMass()
  {
    for (int i = 0; i < sqrt_diagonal_mass.GetM(); i++)
      sqrt_diagonal_mass(i) = sqrt(sqrt_diagonal_mass(i));
  }
  
  
  //! multiplication of X by D^-1/2
  template<class T, class MatStiff, class MatMass> template<class T0>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltInvSqrtDiagonalMass(Vector<T0>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) /= sqrt_diagonal_mass(i);
  }
  
  
  //! multiplication of X by D^1/2
  template<class T, class MatStiff, class MatMass> template<class T0>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltSqrtDiagonalMass(Vector<T0>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) *= sqrt_diagonal_mass(i);
  }
    
  
  //! computation of diagonal of mass matrix
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeDiagonalMass()
  {
    if (Mh == NULL)
      {
        // M = identity
        sqrt_diagonal_mass.Reallocate(this->n_);
        sqrt_diagonal_mass.Fill(1.0);
      }
    else
      {
        sqrt_diagonal_mass.Reallocate(this->n_);
        for (int i = 0; i < this->n_; i++)
          sqrt_diagonal_mass(i) = (*Mh)(i, i);
      }
  }
  
  
  //! computation of mass matrix
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeMassForCholesky()
  {
    // nothing to do, we consider that mass matrix
    // is already computed
  }
  
  
  //! computation of mass matrix M
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::ComputeMassMatrix()
  {
    // mass matrix already computed in Mh
  }
  
  
  //! matrix vector product with mass matrix Y = M X
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltMass(const Vector<T>& X, Vector<T>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mlt(*Mh, X, Y);
  }
  
  
  /**************************
   * Stiffness matrix stuff *
   **************************/
  
  
  //! computation of stiffness matrix K
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::ComputeStiffnessMatrix()
  {
    // nothing to do, already computed in Kh
  }
  
  
  //! computation of matrix a M + b*K
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeStiffnessMatrix(const T& a, const T& b)
  {
    // nothing to do, we use Kh and Mh for the matrix vector product
  }
	
  
  //! matrix vector product with stifness matrix Y = K X
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltStiffness(const Vector<T>& X, Vector<T>& Y)
  {
    if (Kh == NULL)
      PrintErrorInit();
    else
      Mlt(*Kh, X, Y);
  }
  
  
  //! matrix vector product with stifness and mass matrix Y = (a M + b K) X
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltStiffness(const T& coef_mass, const T& coef_stiff,
               const Vector<T>& X, Vector<T>& Y)
  {
    if (Kh == NULL)
      PrintErrorInit();
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
    
  
  
  /*************************
   * Functions to overload *
   *************************/

  
  //! computation of matrix a M + b K and factorisation of this matrix
  /*!
    The factorisation process can be also the construction of preconditioning
    if an iterative solver is used to solve linear system (a M + b K) y = x 
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b)
  {
    abort();
  }
  
  
  //! computation of matrix a M + b K and factorisation of this matrix
  /*!
    The factorisation process can be also the construction of preconditioning
    if an iterative solver is used to solve linear system (a M + b K) y = x 
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const complex<T>& a, const complex<T>& b,
                                     bool real_part)
  {
    abort();
  }
  
  
  //! solving the linear system (a M + b K) Y = X
  template<class T, class MatStiff, class MatMass> template<class T0>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeSolution(const Vector<T0>& X, Vector<T0>& Y)
  {
    abort();
  }
  
  
  //! solving the linear system (a M + b K) Y = X
  template<class T, class MatStiff, class MatMass>
  template<class TransA, class T0>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeSolution(const TransA&, const Vector<T0>& X, Vector<T0>& Y)
  {
    abort();
  }

  
  //! computation of Cholesky factorisation of M from matrix M
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::FactorizeCholeskyMass()
  {
    abort();
  }
  
  
  //! computation of L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltCholeskyMass(const TransStatus& TransA, Vector<T>& X)
  {
    abort();
  }
  
  
  //! computation of L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SolveCholeskyMass(const TransStatus& TransA, Vector<T>& X)
  {
    abort();
  }
  
  
  //! memory release
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::Clear()
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
  template<class EigenPb, class Vector1, class Matrix1, class T0>
  void ApplyScalingEigenvec(EigenPb& var, Vector1& eigen_values, Vector1& lambda_imag,
                            Matrix1& eigen_vectors,
                            const T0& shiftr, const T0& shifti)
  {
    
    if (var.DiagonalMass())
      {
        // scaling to have true eigenvectors
        for (int i = 0; i < var.sqrt_diagonal_mass.GetM(); i++)
          for (int j = 0; j < eigen_vectors.GetN(); j++)
            eigen_vectors(i,j) /= var.sqrt_diagonal_mass(i);
      }      
    else if (var.UseCholeskyFactoForMass())
      {
        Vector<T0> Xcol(eigen_vectors.GetM());
        for (int j = 0; j < eigen_vectors.GetN(); j++)
          {
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              Xcol(i) = eigen_vectors(i,j);
            
            var.SolveCholeskyMass(SeldonTrans, Xcol);
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              eigen_vectors(i,j) = Xcol(i);
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
  template<class EigenPb, class Vector1, class Matrix1, class T0>
  void ApplyScalingEigenvec(EigenPb& var, Vector1& eigen_values, Vector1& lambda_imag,
                            Matrix1& eigen_vectors,
                            const complex<T0>& shiftr, const complex<T0>& shifti)
  {
    
    if (var.DiagonalMass())
      {
        // scaling to have true eigenvectors
        for (int i = 0; i < var.sqrt_diagonal_mass.GetM(); i++)
          for (int j = 0; j < eigen_vectors.GetN(); j++)
            eigen_vectors(i,j) /= var.sqrt_diagonal_mass(i);
      }      
    else if (var.UseCholeskyFactoForMass())
      {
        Vector<complex<T0> > Xcol(eigen_vectors.GetM());
        for (int j = 0; j < eigen_vectors.GetN(); j++)
          {
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              Xcol(i) = eigen_vectors(i,j);
            
            var.SolveCholeskyMass(SeldonTrans, Xcol);
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              eigen_vectors(i,j) = Xcol(i);
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
  template<class T, class Storage1, class Storage2,
           class Alloc1, class Alloc2>
  void SortEigenvalues(Vector<T>& lambda_r, Vector<T>& lambda_i,
                       Matrix<T, General, Storage1, Alloc1>& eigen_old,
                       Matrix<T, General, Storage2, Alloc2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const T& shift_r, const T& shift_i)
  {
    int n = min(lambda_r.GetM(), eigen_old.GetN());;
				 
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
  template<class T, class Storage1, class Storage2,
           class Alloc1, class Alloc2>
  void SortEigenvalues(Vector<complex<T> >& lambda_r, Vector<complex<T> >& lambda_i,
                       Matrix<complex<T>, General, Storage1, Alloc1>& eigen_old,
                       Matrix<complex<T>, General, Storage2, Alloc2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const complex<T>& shift_r, const complex<T>& shift_i)
  {
    // complex case, ignoring lambda_i and shift_i
    int n = min(lambda_r.GetM(), eigen_old.GetN());;
    
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
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  DenseEigenProblem() : EigenProblem_Base<T, Matrix<T, Prop, Storage>,
                                          Matrix<Tmass, PropM, StorageM> >()
  {    
  }
  
  
  //! Cholesky factorisation of mass matrix
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  FactorizeCholeskyMass()
  {
    if (this->Mh == NULL)
      {
        mat_chol.Reallocate(this->n_, this->n_);
        mat_chol.SetIdentity();
      }
    else
      {
        mat_chol = *(this->Mh);
        GetCholesky(mat_chol);
        Xchol_real.Reallocate(mat_chol.GetM());
        Xchol_imag.Reallocate(mat_chol.GetM());
      }
  }
  
  
  //! computation of L X or L^T X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransStatus, class T0>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  MltCholeskyMass(const TransStatus& transA, Vector<T0>& X)
  {
    MltCholesky(transA, mat_chol, X);
  }
  
  
  //! computation of L^-1 X or L^-T X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransStatus, class T0>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  SolveCholeskyMass(const TransStatus& transA, Vector<T0>& X)
  {
    SolveCholesky(transA, mat_chol, X);
  }
  
  
  //! computation of L X or L^T X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransStatus>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  MltCholeskyMass(const TransStatus& transA, Vector<complex<double> >& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    MltCholesky(transA, mat_chol, Xchol_real);
    MltCholesky(transA, mat_chol, Xchol_imag);

    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<double>(Xchol_real(i), Xchol_imag(i));    
  }
  
  
  //! computation of L^-1 X or L^-T X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransStatus>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  SolveCholeskyMass(const TransStatus& transA, Vector<complex<double> >& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    SolveCholesky(transA, mat_chol, Xchol_real);
    SolveCholesky(transA, mat_chol, Xchol_imag);

    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<double>(Xchol_real(i), Xchol_imag(i));
  }

  
  //! computation and factorisation of a M + b K
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b)
  {
    if (this->Kh == NULL)
      this->PrintErrorInit();
    
    this->complex_system = false;
    // computation of mat_lu = a M + b K
    mat_lu = *(this->Kh);
    Mlt(b, mat_lu);
    if (this->Mh == NULL)
      {
        for (int i = 0; i < this->n_; i++)
          mat_lu(i, i) += a;
      }
    else
      Add(a, *(this->Mh), mat_lu);
    
    // factorisation
    GetLU(mat_lu, pivot);
  }
  
  
  //! computation and factorisation of a M + b K
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeStiffnessMatrix(const complex<T>& a, const complex<T>& b,
                                     bool real_part)
  {
    ComputeAndFactorizeComplexMatrix(a, b);
  }
  
  
  //! computation and factorisation of a M + b K
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeComplexMatrix(const complex<double>& a, const complex<double>& b,
                                   bool real_part)
  {
    if (this->Kh == NULL)
      this->PrintErrorInit();
 
    this->complex_system = true;
    // inverse of (a M + b K), then we take real_part or imaginary part
    Matrix<complex<double>, Prop, Storage> InvMat(this->n_, this->n_);
    for (int i = 0; i < this->n_; i++)
      for (int j = 0; j < this->n_; j++)
        InvMat(i, j) = (*this->Kh)(i, j);
    
    Mlt(b, InvMat);
    if (this->Mh == NULL)
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
    
    // inverse
    GetInverse(InvMat);
  
    // then extracting real or imaginary part
    mat_lu.Reallocate(this->n_, this->n_);
    if (real_part)
      {
        for (int i = 0; i < this->n_; i++)
          for (int j = 0; j < this->n_; j++)
            mat_lu(i, j) = real(InvMat(i, j));
      }
    else
      {
        for (int i = 0; i < this->n_; i++)
          for (int j = 0; j < this->n_; j++)
            mat_lu(i, j) = imag(InvMat(i, j));
      }
  }
  

  //! computation and factorisation of a M + b K
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM> template<class T0>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeComplexMatrix(const complex<T0>& a,
                                   const complex<T0>& b, bool real_p)
  {
    // this case should not appear
    cout << "Case not handled" << endl;
    abort();
  }
  
  
  //! solution of (a M + b K) Y = X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM> template<class T0>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeSolution(const Vector<T0>& X, Vector<T0>& Y)
  {
    if (this->complex_system)
      Mlt(mat_lu, X, Y);
    else
      {
        Copy(X, Y);
        SolveLU(mat_lu, pivot, Y);
      }
  }
   

  //! solution of (a M + b K) Y = X or transpose system
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransA, class T0>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeSolution(const TransA& transA, const Vector<T0>& X, Vector<T0>& Y)
  {
    if (this->complex_system)
      Mlt(transA, mat_lu, X, Y);
    else
      {
        Copy(X, Y);
        SolveLU(transA, mat_lu, pivot, Y);
      }
  }
  
  
  //! clearing variables used for eigenvalue resolution
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  Clear()
  {
    EigenProblem_Base<T, Matrix<T, Prop, Storage>,
      Matrix<Tmass, PropM, StorageM> >::Clear();
    
    mat_lu.Clear();
    mat_chol.Clear();
  }
  
  
  /**********************
   * SparseEigenProblem *
   **********************/
  
  
  //! default constructor
  template<class T, class MatStiff, class MatMass>
  SparseEigenProblem<T, MatStiff, MatMass>::
  SparseEigenProblem()
    : EigenProblem_Base<T, MatStiff, MatMass>()
  {
    imag_solution = false;
    mat_lu.RefineSolution();
    mat_lu_cplx.RefineSolution();
  }
    
  
  //! sets Cholesky solver to use
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::SelectCholeskySolver(int type)
  {
    chol_facto_mass_matrix.SelectDirectSolver(type);
  }
  

  //! computes Cholesky factorisation of M from matrix M
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  FactorizeCholeskyMass()
  {
    if (this->Mh == NULL)
      this->PrintErrorInit();
    
    if (this->print_level > 0)
      chol_facto_mass_matrix.ShowMessages();
    
    MassValue x_test;
    FactorizeCholeskyMass(x_test, *this->Mh);
    
    if (this->print_level < 2)
      chol_facto_mass_matrix.HideMessages();
    
    Xchol_real.Reallocate(this->n_);
    Xchol_imag.Reallocate(this->n_);
  }
  
  
  //! intermediary function
  template<class T, class MatStiff, class MatMass>
  template<class Storage, class Alloc>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  FactorizeCholeskyMass(double&, Matrix<double, Symmetric, Storage, Alloc>& M)
  {
    chol_facto_mass_matrix.Factorize(M, true);    
  }


  //! intermediary function
  template<class T, class MatStiff, class MatMass>
  template<class T0, class T1, class Prop, class Storage, class Alloc>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  FactorizeCholeskyMass(T0&, Matrix<T1, Prop, Storage, Alloc>& M)
  {
    cout << "Cholesky factorisation has not been implemented "
	 << "for complex matrices" << endl;
    abort();    
  }
  
  
  //! computes L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus, class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  MltCholeskyMass(const TransStatus& TransA, Vector<T0>& X)
  {
    chol_facto_mass_matrix.Mlt(TransA, X);
  }
  
  
  //! computes L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus, class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  SolveCholeskyMass(const TransStatus& TransA, Vector<T0>& X)
  {
    chol_facto_mass_matrix.Solve(TransA, X);
  }


  //! computes L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  MltCholeskyMass(const TransStatus& TransA, Vector<complex<double> >& X)
  {
    MassValue x_test;
    MltCholeskyMass(TransA, X, x_test);
  }
  

  //! computes L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  MltCholeskyMass(const TransStatus& TransA, Vector<complex<double> >& X,
		  complex<double>&)
  {
    cout << "Cholesky factorisation has not been implemented "
	 << "for complex matrices" << endl;
    abort();
  }
  
  
  //! computes L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  MltCholeskyMass(const TransStatus& TransA, Vector<complex<double> >& X, double&)  
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    chol_facto_mass_matrix.Mlt(TransA, Xchol_real);
    chol_facto_mass_matrix.Mlt(TransA, Xchol_imag);
    
    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<double>(Xchol_real(i), Xchol_imag(i));
    
  }
  
  
  //! computes L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  SolveCholeskyMass(const TransStatus& TransA, Vector<complex<double> >& X)
  {
    MassValue x_test;
    SolveCholeskyMass(TransA, X, x_test);
  }
  

  //! computes L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  SolveCholeskyMass(const TransStatus& TransA, Vector<complex<double> >& X,
		    complex<double>&)
  {
    cout << "Cholesky factorisation has not been "
	 << "implemented for complex matrices" << endl;
    abort();
  }
  
  
  //! computes L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  SolveCholeskyMass(const TransStatus& TransA, Vector<complex<double> >& X, double&)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    chol_facto_mass_matrix.Solve(TransA, Xchol_real);
    chol_facto_mass_matrix.Solve(TransA, Xchol_imag);
    
    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<double>(Xchol_real(i), Xchol_imag(i));
  }
  
  
  //! computes and factorizes a M + b K
  //! where M is the mass matrix and K the stiffness matrix
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b)
  {
    this->complex_system = false;
    if (this->Kh == NULL)
      this->PrintErrorInit();
    
    if (this->print_level > 0)
      mat_lu.ShowMessages();
    
    // retrieving symmetry of mass matrix 
    bool sym_mh = true;
    if (this->Mh != NULL)
      {
	if (!IsSymmetricMatrix(*this->Mh))
	  sym_mh = false;
      }

    T zero; SetComplexZero(zero);
    
    if (b == zero)
      {
	// only mass matrix must be factorized
	if (sym_mh)
	  {
	    Matrix<T, Symmetric, ArrayRowSymSparse> A;
	    if (this->Mh == NULL)
	      {
		A.Reallocate(this->n_, this->n_);
		A.SetIdentity();
	      }
	    else
	      Copy(*(this->Mh), A);
	    
	    Mlt(a, A);
	    mat_lu.Factorize(A);
	  }
	else
	  {
	    Matrix<T, General, ArrayRowSparse> A;
	    Copy(*(this->Mh), A);
	    Mlt(a, A);
	    mat_lu.Factorize(A);
	  }
      }
    else if (IsSymmetricMatrix(*this->Kh) && sym_mh)
      {
	// forming a M + b K and factorizing it when the result is symmetric
        Matrix<T, Symmetric, ArrayRowSymSparse> A;
        Copy(*(this->Kh), A);
        Mlt(b, A);
	if (a != zero)
	  {
	    if (this->Mh == NULL)
	      {
		for (int i = 0; i < this->n_; i++)
		  A.AddInteraction(i, i, a);
	      }
	    else
	      Add(a, *(this->Mh), A);
	  }
	
        mat_lu.Factorize(A);
      } 
    else
      {
	// forming a M + b K and factorizing it when the result is unsymmetric
        Matrix<T, General, ArrayRowSparse> A;
        Copy(*(this->Kh), A);
        Mlt(b, A);
	if (a != zero)
	  {
	    if (this->Mh == NULL)
	      {
		for (int i = 0; i < this->n_; i++)
		  A.AddInteraction(i, i, a);
	      }
	    else
	      Add(a, *(this->Mh), A);
	  }
	
        mat_lu.Factorize(A);
      }      

    if (this->print_level < 2)
      mat_lu.HideMessages();
  }
  
  
  //! computes and factorizes matrix (a M + b K) for complex values of a and b
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const complex<T>& a,
                                     const complex<T>& b,
                                     bool real_p)
  {
    ComputeAndFactorizeComplexMatrix(a, b, real_p);
  }
  
  
  //! intermediary function
  template<class T, class MatStiff, class MatMass> template<class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeComplexMatrix(const complex<T0>& a,
                                   const complex<T0>& b,
                                   bool real_p)
  {
    // this case should not appear
    cout << "Case not handled" << endl;
    abort();
  }

  
  //! computes and factorizes matrix (a M + b K) for complex values of a and b
  template<class T, class MatStiff, class MatMass> 
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeComplexMatrix(const complex<double>& a,
                                   const complex<double>& b,
                                   bool real_p)
  {
    this->complex_system = true;
    if (this->Kh == NULL)
      this->PrintErrorInit();
    
    imag_solution = !real_p;

    if (this->print_level > 0)
      mat_lu_cplx.ShowMessages();
    
    // retrieving symmetry of mass matrix 
    bool sym_mh = true;
    if (this->Mh != NULL)
      {
	if (!IsSymmetricMatrix(*this->Mh))
	  sym_mh = false;
      }

    complex<double> zero(0, 0);
    
    if (b == zero)
      {
	// only mass matrix must be factorized
	if (sym_mh)
	  {
	    Matrix<Complexe, Symmetric, ArrayRowSymSparse> A;
	    if (this->Mh == NULL)
	      {
		A.Reallocate(this->n_, this->n_);
		A.SetIdentity();
	      }
	    else
	      Copy(*(this->Mh), A);
	    
	    Mlt(a, A);
	    mat_lu_cplx.Factorize(A);
	  }
	else
	  {
	    Matrix<Complexe, General, ArrayRowSparse> A;
	    Copy(*(this->Mh), A);
	    Mlt(a, A);
	    mat_lu_cplx.Factorize(A);
	  }
      }
    else if (IsSymmetricMatrix(*this->Kh) && sym_mh)
      {
	// forming a M + b K
	Matrix<Complexe, Symmetric, ArrayRowSymSparse> A;
	Copy(*(this->Kh), A);
	Mlt(b, A);
	if (a != zero)
	  {
	    if (this->Mh == NULL)
	      {
		for (int i = 0; i < this->n_; i++)
		  A.AddInteraction(i, i, a);
	      }
	    else
	      Add(a, *(this->Mh), A);
	  }
	
	mat_lu_cplx.Factorize(A);
      }
    else
      {
	// forming a M + b K
	Matrix<Complexe, General, ArrayRowSparse> A;
	Copy(*(this->Kh), A);
	Mlt(b, A);
	if (a != zero)
	  {
	    if (this->Mh == NULL)
	      {
		for (int i = 0; i < this->n_; i++)
		  A.AddInteraction(i, i, a);
	      }
	    else
	      Add(a, *(this->Mh), A);
	  }
	
	mat_lu_cplx.Factorize(A);
      }

    if (this->print_level < 2)
      mat_lu_cplx.HideMessages();
  }

  
  //! solves (a M + b K) Y = X with stored factorization 
  template<class T, class MatStiff, class MatMass> template<class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeSolution(const Vector<T0>& X, Vector<T0>& Y)
  {
    if (this->complex_system)
      {
        ComputeComplexSolution(SeldonNoTrans, X, Y);
      }
    else
      {
        Copy(X, Y);
        mat_lu.Solve(Y);
      }
  }
  
  
  //! solves (a M + b K) Y = X or transpose system
  template<class T, class MatStiff, class MatMass>
  template<class TransA, class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeSolution(const TransA& transA, const Vector<T0>& X, Vector<T0>& Y)
  {
    if (this->complex_system)
      {
        ComputeComplexSolution(transA, X, Y);
      }
    else
      {
        Copy(X, Y);
        mat_lu.Solve(transA, Y);
      }
  }
  
  
  //! solves (a M + b K) Y = X when a and b are complex
  template<class T, class MatStiff, class MatMass> template<class TransA>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeComplexSolution(const TransA& transA,
                         const Vector<double>& X, Vector<double>& Y)
  {
    Vector<complex<T> > Xcplx(this->n_);
    for (int i = 0; i < this->n_; i++)
      Xcplx(i) = X(i);
    
    mat_lu_cplx.Solve(Xcplx);
    
    if (imag_solution)
      for (int i = 0; i < this->n_; i++)
        Y(i) = imag(Xcplx(i));
    else
      for (int i = 0; i < this->n_; i++)
        Y(i) = real(Xcplx(i));
    
  }


  //! intermediary function
  template<class T, class MatStiff, class MatMass> template<class TransA>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeComplexSolution(const TransA& transA,
                         const Vector<complex<double> >& X,
                         Vector<complex<double> >& Y)
  {
    Copy(X, Y);
    mat_lu_cplx.Solve(transA, Y);
  }

  
  //! clears memory used by the object
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::Clear()
  {
    mat_lu.Clear();
    mat_lu_cplx.Clear();
    chol_facto_mass_matrix.Clear();
    Xchol_real.Clear();
    Xchol_imag.Clear();    
  }

  
  //! default constructor
  template<class T, class MatStiff>
  MatrixFreeEigenProblem<T, MatStiff>::
  MatrixFreeEigenProblem() 
    : EigenProblem_Base<T, MatStiff, Matrix<double, Symmetric, ArrayRowSymSparse> >()
  {
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  int TypeEigenvalueSolver::default_solver(0);  
#endif


  int TypeEigenvalueSolver::GetDefaultSolver()
  {
#ifdef SELDON_WITH_ANASAZI
    return ANASAZI;
#endif
    
#ifdef SELDON_WITH_FEAST
    return FEAST;
#endif
    
#ifdef SELDON_WITH_ARPACK
    return ARPACK;
#endif
    
    return -1;
  }
  

  template<class EigenPb, class Vector1, class Vector2,
            class T, class Prop, class Storage, class Alloc3>
  void GetEigenvaluesEigenvectors(EigenPb& var_eig, Vector1& lambda,
				  Vector2& lambda_imag,
				  Matrix<T, Prop, Storage, Alloc3>& eigen_vec)
  {
    int type_solver = TypeEigenvalueSolver::default_solver;
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
        cout << "Recompile with MKL" << endl;
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

#define SELDON_FILE_EIGENVALUE_SOLVER_CXX
#endif
