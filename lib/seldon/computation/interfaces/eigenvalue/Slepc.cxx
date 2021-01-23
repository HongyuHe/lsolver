#ifndef SELDON_FILE_SLEPC_CXX

#include <slepceps.h>
#include <slepcpep.h>
#include <slepcnep.h>
#include <slepc/private/stimpl.h>

namespace Seldon
{

  const char* SlepcParam::GetEigensolverChar() const
  {
    switch (type_solver)
      {
      case POWER : return EPSPOWER;
      case SUBSPACE : return EPSSUBSPACE;
      case ARNOLDI : return EPSARNOLDI;
      case LANCZOS : return EPSLANCZOS;
      case KRYLOVSCHUR : return EPSKRYLOVSCHUR;
      case GD : return EPSGD;
      case JD : return EPSJD;
      case RQCG : return EPSRQCG;
      case LOBPCG : return EPSLOBPCG;
      case CISS : return EPSCISS;
      case LAPACK : return EPSLAPACK;
      case ARPACK : return EPSARPACK;
      case BLZPACK : return EPSBLZPACK;
      case TRLAN : return EPSTRLAN;
      case BLOPEX : return EPSBLOPEX;
      case PRIMME : return EPSPRIMME;
        //case FEAST : return EPSFEAST;
      }
    
    return "";
  }


  void SetParametersSlepc(const SlepcParam& param, EPS& eps)
  {
    int ierr = EPSSetType(eps, param.GetEigensolverChar());
    if (ierr != 0)
      {
	cout << "Chosen type = " << param.GetEigensolverChar() << " Not enabled in Slepc ?" << endl;
	abort();
      }

    int solver = param.GetEigensolverType();
    if (solver == param.BLOPEX)
      {
#ifdef SLEPC_HAVE_BLOPEX
	PetscInt bs = param.GetBlockSize();
	if (bs > 0)
	  EPSBLOPEXSetBlockSize(eps, bs);
#else
	cout << "Slepc not compiled with BLOPEX" << endl;
	abort();
#endif
      }
    else if (solver == param.BLZPACK)
      {
#ifdef SLEPC_HAVE_BLZPACK
	PetscInt bs = param.GetBlockSize();
	if (bs > 0)
	  EPSBlzpackSetBlockSize(eps, bs);

	bs = param.GetNumberOfSteps();
	if (bs > 0)
	  EPSBlzpackSetNSteps(eps, bs);
#else
	cout << "Slepc not compiled with BLZPACK" << endl;
	abort();
#endif
      }
    else if (solver == param.CISS)
      {
	int type_extraction = param.GetExtractionType();
	if (type_extraction >= 0)
	  {
	    if (type_extraction == param.EXTRACT_RITZ)
	      EPSCISSSetExtraction(eps, EPS_CISS_EXTRACTION_RITZ);
	    else
	      EPSCISSSetExtraction(eps, EPS_CISS_EXTRACTION_HANKEL);
	  }
	
	int type_quad = param.GetQuadratureRuleType();
	if (type_quad >= 0)
	  {
	    if (type_quad == param.QUADRULE_TRAPEZE)
	      EPSCISSSetQuadRule(eps, EPS_CISS_QUADRULE_TRAPEZOIDAL);
	    else
	      EPSCISSSetQuadRule(eps, EPS_CISS_QUADRULE_CHEBYSHEV);
	  }
	
	PetscScalar a; PetscBool complex_number = PETSC_FALSE;
	if (IsComplexNumber(a))
	  complex_number = PETSC_TRUE;
	
	if (param.GetInnerSteps() > 0)
	  EPSCISSSetRefinement(eps, param.GetInnerSteps(), param.GetOuterSteps());
	
	if (param.GetNumberIntegrationPoints() > 0)
	  EPSCISSSetSizes(eps, param.GetNumberIntegrationPoints(), param.GetBlockSize(),
			  param.GetMomentSize(), param.GetNumberPartitions(),
			  param.GetMaximumBlockSize(), complex_number);
	
	if (param.GetThresholdRank() > 0)
	  EPSCISSSetThreshold(eps, param.GetThresholdRank(), param.GetThresholdSpurious());
      }
    else if (solver == param.FEAST)
      {
#ifdef SLEPC_HAVE_FEAST
	if (param.GetNumberIntegrationPoints() > 0)
	  EPSFEASTSetNumPoints(eps, param.GetNumberIntegrationPoints());
#else
	cout << "Slepc not compiled with FEAST" << endl;
	abort();
#endif
      }    
    else if (solver == param.GD)
      {
	if (param.GetBorthogonalization() >= 0)
	  {
	    PetscBool borth = PETSC_FALSE;
	    if (param.GetBorthogonalization() >= 1)
	      borth = PETSC_TRUE;
	    
	    EPSGDSetBOrth(eps, borth);
	  }

	PetscInt bs = param.GetBlockSize();
	if (bs > 0)
	  EPSGDSetBlockSize(eps, bs);

	if (param.GetDoubleExpansion() >= 0)
	  {
	    PetscBool exp = PETSC_FALSE;
	    if (param.GetDoubleExpansion() >= 1)
	      exp = PETSC_TRUE;
	    
	    EPSGDSetDoubleExpansion(eps, exp);
	  }
	
	if (param.GetInitialSize() > 0)
	  EPSGDSetInitialSize(eps, param.GetInitialSize());

	if (param.GetKrylovRestart() >= 0)
	  {
	    PetscBool restart = PETSC_FALSE;
	    if (param.GetKrylovRestart() >= 1)
	      restart = PETSC_TRUE;

	    EPSGDSetKrylovStart(eps, restart);
	  }

	if (param.GetRestartNumber() > 0)
	  EPSGDSetRestart(eps, param.GetRestartNumber(), param.GetRestartNumberAdd());
	
        // deprecated function : SetWindowSizes 
	//if (param.GetNumberConvergedVectors() > 0)
        // EPSGDSetWindowSizes(eps, param.GetNumberConvergedVectors(),
	//		      param.GetNumberConvergedVectorsProjected());
      }
    else if (solver == param.JD)
      {
	if (param.GetBorthogonalization() >= 0)
	  {
	    PetscBool borth = PETSC_FALSE;
	    if (param.GetBorthogonalization() >= 1)
	      borth = PETSC_TRUE;

	    EPSJDSetBOrth(eps, borth);
	  }

	PetscInt bs = param.GetBlockSize();
	if (bs > 0)
	  EPSJDSetBlockSize(eps, bs);

	if (param.GetInitialSize() > 0)
	  EPSJDSetInitialSize(eps, param.GetInitialSize());

	if (param.GetKrylovRestart() >= 0)
	  {
	    PetscBool restart = PETSC_FALSE;
	    if (param.GetKrylovRestart() >= 1)
	      restart = PETSC_TRUE;

	    EPSJDSetKrylovStart(eps, restart);
	  }

	if (param.GetRestartNumber() > 0)
	  EPSJDSetRestart(eps, param.GetRestartNumber(), param.GetRestartNumberAdd());

        // deprecated function : SetWindowSizes 	
	//if (param.GetNumberConvergedVectors() > 0)
	//  EPSJDSetWindowSizes(eps, param.GetNumberConvergedVectors(),
	//		      param.GetNumberConvergedVectorsProjected());
      }
    else if (solver == param.KRYLOVSCHUR)
      {
	if (param.UseNonLockingVariant())
	  EPSKrylovSchurSetLocking(eps, PETSC_FALSE);
	else
	  EPSKrylovSchurSetLocking(eps, PETSC_TRUE);

	if (param.GetRestartRatio() > 0)
	  EPSKrylovSchurSetRestart(eps, param.GetRestartRatio());	
      }
    else if (solver == param.LOBPCG)
      {
	if (param.GetBlockSize() > 0)
	  EPSLOBPCGSetBlockSize(eps, param.GetBlockSize());
	
	if (param.UseNonLockingVariant())
	  EPSLOBPCGSetLocking(eps, PETSC_FALSE);
	else
	  EPSLOBPCGSetLocking(eps, PETSC_TRUE);

	if (param.GetRestartRatio() > 0)
	  EPSLOBPCGSetRestart(eps, param.GetRestartRatio());	
      }
    else if (solver == param.PRIMME)
      {
#ifdef SLEPC_HAVE_PRIMME
	if (param.GetBlockSize() > 0)
	  EPSPRIMMESetBlockSize(eps, param.GetBlockSize());

	if (param.GetMethod().size() > 1)
	  {
	    if (param.GetMethod() == "DYNAMIC")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_DYNAMIC);
	    else if (param.GetMethod() == "DEFAULT_MIN_TIME")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_DEFAULT_MIN_TIME);
	    else if (param.GetMethod() == "DEFAULT_MIN_MATVECS")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_DEFAULT_MIN_MATVECS);
	    else if (param.GetMethod() == "ARNOLDI")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_ARNOLDI);
	    else if (param.GetMethod() == "GD")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_GD);
	    else if (param.GetMethod() == "GD_PLUSK")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_GD_PLUSK);
	    else if (param.GetMethod() == "GD_OLSEN_PLUSK")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_GD_OLSEN_PLUSK);
	    else if (param.GetMethod() == "JD_OLSEN_PLUSK")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_JD_OLSEN_PLUSK);
	    else if (param.GetMethod() == "RQI")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_RQI);
	    else if (param.GetMethod() == "JDQR")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_JDQR);
	    else if (param.GetMethod() == "JDQMR")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_JDQMR);
	    else if (param.GetMethod() == "JDQMR_ETOL")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_JDQMR_ETOL);
	    else if (param.GetMethod() == "SUBSPACE_ITERATION")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_SUBSPACE_ITERATION);
	    else if (param.GetMethod() == "LOBPCG_ORTHOBASIS")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_LOBPCG_ORTHOBASIS);
	    else if (param.GetMethod() == "LOBPCG_ORTHOBASISW")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_LOBPCG_ORTHOBASISW);
	  }
#else
	cout << "Slepc not compiled with PRIMME" << endl;
	abort();
#endif
      }
    else if (solver == param.POWER)
      {
	if (param.GetShiftType() >= 0)
	  {
	    if (param.GetShiftType() == param.SHIFT_CONSTANT)
	      EPSPowerSetShiftType(eps, EPS_POWER_SHIFT_CONSTANT);
	    else if (param.GetShiftType() == param.SHIFT_RAYLEIGH)
	      EPSPowerSetShiftType(eps, EPS_POWER_SHIFT_RAYLEIGH);
	    else if (param.GetShiftType() == param.SHIFT_WILKINSON)
	      EPSPowerSetShiftType(eps, EPS_POWER_SHIFT_WILKINSON);
	  }
      }
    else if (solver == param.RQCG)
      {
	if (param.GetNumberOfSteps() > 0)
	  EPSRQCGSetReset(eps, param.GetNumberOfSteps());	
      }
  }

  
  //! filling the vector y from pointer contained in x (Petsc vector)
  void CopyPointerPetsc(const Vec& x, Vector<PetscScalar>& y)
  {
    // it is assumed that the vector x is stored in a contiguous array
    PetscInt n;
    VecGetLocalSize(x, &n);
    
    PetscScalar* x_array;
    VecGetArrayRead(x, (const PetscScalar**)&x_array);
    
    y.SetData(n, x_array);
  }


  void AllocatePetscVector(const MPI_Comm& comm, Vec& Vr, int n, int nglob,
			   Vector<PetscScalar>& Vr_vec)
  {
    VecCreate(comm, &Vr);
    VecSetSizes(Vr, n, nglob);
    VecSetFromOptions(Vr);

    CopyPointerPetsc(Vr, Vr_vec);
  }

  
  //! matrix-vector product y = mat x (mat : stiffness operator)
  PetscErrorCode MatMult_Matrix(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec0, yvec;
    CopyPointerPetsc(x, xvec0);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    EigenProblem_Base<PetscScalar>& var_eig
      = *reinterpret_cast<EigenProblem_Base<PetscScalar>* >(ctx);

    var_eig.IncrementProdMatVect();
    Vector<PetscScalar> xvec(xvec0);
    if (var_eig.DiagonalMass() || (var_eig.UseCholeskyFactoForMass()))
      {
	// standard eigenvalue problem
	if (var_eig.GetComputationalMode() == var_eig.REGULAR_MODE)
	  {
	    if (var_eig.DiagonalMass())
	      var_eig.MltInvSqrtDiagonalMass(xvec);
	    else
	      var_eig.SolveCholeskyMass(Seldon::SeldonTrans, xvec);
	    
	    var_eig.MltStiffness(xvec, yvec);
	    
	    if (var_eig.DiagonalMass())
	      var_eig.MltInvSqrtDiagonalMass(yvec);
	    else
	      var_eig.SolveCholeskyMass(Seldon::SeldonNoTrans, yvec);
	  }
	else
	  {
	    if (var_eig.DiagonalMass())
	      var_eig.MltSqrtDiagonalMass(xvec);
	    else
	      var_eig.MltCholeskyMass(Seldon::SeldonNoTrans, xvec);
	    
	    var_eig.ComputeSolution(xvec, yvec);
	    
	    if (var_eig.DiagonalMass())
	      var_eig.MltSqrtDiagonalMass(yvec);
	    else
	      var_eig.MltCholeskyMass(Seldon::SeldonTrans, yvec);
	  }	
      }
    else
      {
	if (var_eig.GetComputationalMode() == var_eig.INVERT_MODE)
	  {
	    if (var_eig.GetTypeSpectrum() != var_eig.CENTERED_EIGENVALUES)
	      {
		var_eig.MltStiffness(xvec0, xvec);
		var_eig.ComputeSolution(xvec, yvec);
	      }
	    else
	      {
		var_eig.MltMass(xvec0, xvec);
		var_eig.ComputeSolution(xvec, yvec);
	      }
	  }
	else if (var_eig.GetComputationalMode() == var_eig.REGULAR_MODE)
	  {
	    var_eig.MltStiffness(xvec, yvec);
	  }
	else
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
      }

    xvec0.Nullify(); yvec.Nullify();
    return 0;
  }


  //! matrix-vector product y = mat x (mat : mass operator)
  PetscErrorCode MatMult_MassMatrix(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    EigenProblem_Base<PetscScalar>& var_eig
      = *reinterpret_cast<EigenProblem_Base<PetscScalar>* >(ctx);

    var_eig.MltMass(xvec, yvec);
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }


  //! matrix-vector product with transpose, y = mat^T x
  PetscErrorCode MatMultTranspose_Matrix(Mat A, Vec x, Vec y)
  {
    throw Undefined("Function MatMultTranspose_Matrix not implemented");
    return 0;
  }

  
  //! retrieves diagonal of A, d = diag(A)
  PetscErrorCode MatGetDiagonal_Matrix(Mat A, Vec d)
  {
    throw Undefined("Function MatGetDiagonal_Matrix not implemented");
    return 0;
  }


  template<class T>
  bool PutEigenpairLapackForm(int num, int nev, T& Lr, T& Li, Vector<T>& Vr, Vector<T>& Vi,
			      T& Lr_next, T& Li_next, Vector<T>& Vr_next, Vector<T>& Vi_next,
			      Vector<T>& eigen_values, Vector<T>& lambda_imag,
			      Matrix<T, General, ColMajor>& eigen_vectors)
  {
    bool eigen_pair = false;
    if ((Li != T(0)) && (num < nev-1))
      eigen_pair = true;

    int n = Vr.GetM();
    if (eigen_pair)
      {
	eigen_values(num) = Lr;
	lambda_imag(num) = Li;
	eigen_values(num+1) = Lr_next;
	lambda_imag(num+1) = Li_next;
	for (int j = 0; j < n; j++)	  
	  {
	    eigen_vectors(j, num) = Vr(j);
	    eigen_vectors(j, num+1) = Vi(j);
	  }
      }
    else
      {
	eigen_values(num) = Lr;
	lambda_imag(num) = Li;
	for (int j = 0; j < n; j++)
	  eigen_vectors(j, num) = Vr(j);
      }

    return eigen_pair;
  }


  template<class T>
  bool PutEigenpairLapackForm(int num, int nev, complex<T>& Lr, complex<T>& Li,
			      Vector<complex<T> >& Vr, Vector<complex<T> >& Vi,
			      complex<T>& Lr_next, complex<T>& Li_next,
			      Vector<complex<T> >& Vr_next, Vector<complex<T> >& Vi_next,
			      Vector<complex<T> >& eigen_values, Vector<complex<T> >& lambda_imag,
			      Matrix<complex<T>, General, ColMajor>& eigen_vectors)
  {
    int n = Vr.GetM();
    eigen_values(num) = Lr;
    lambda_imag(num) = Li;
    for (int j = 0; j < n; j++)
      eigen_vectors(j, num) = Vr(j);
    
    return false;
  }


  void FindEigenvaluesSlepc_(EigenProblem_Base<PetscScalar>& var,
			     Vector<PetscScalar>& eigen_values,
			     Vector<PetscScalar>& lambda_imag,
			     Matrix<PetscScalar, General, ColMajor>& eigen_vectors)
  {
    // initializing of computation
    PetscScalar shiftr = var.GetShiftValue(), shifti = var.GetImagShiftValue();    
    PetscScalar zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
        
    int print_level = var.GetPrintLevel();
    SlepcParam& param = var.GetSlepcParameters();
    int rank_proc; MPI_Comm_rank(var.GetCommunicator(), &rank_proc);

    Mat stiff, mass;
    int nev = var.GetNbAskedEigenvalues();
    int ncv = var.GetNbArnoldiVectors();
    int n = var.GetM();
    int nglob = var.GetGlobalM();

    // creation of a matrix-free Petsc structure
    MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
		   reinterpret_cast<void*>(&var), &stiff);

    MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
		   reinterpret_cast<void*>(&var), &mass);

    //MatSetFromOptions(stiff);     MatSetFromOptions(mass);
    MatShellSetOperation(stiff, MATOP_MULT, (void(*)())MatMult_Matrix);
    MatShellSetOperation(stiff, MATOP_MULT_TRANSPOSE, (void(*)())MatMultTranspose_Matrix);
    MatShellSetOperation(stiff, MATOP_GET_DIAGONAL, (void(*)())MatGetDiagonal_Matrix);

    MatShellSetOperation(mass, MATOP_MULT, (void(*)())MatMult_MassMatrix);
    
    // creation of the eigensolver
    EPS eps;
    EPSCreate(var.GetCommunicator(), &eps);
    SetParametersSlepc(param, eps);
    
    // type of eigenproblem (hermitian/generalized)
    bool generalized = true;
    bool isherm = var.IsHermitianProblem();
    if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
      generalized = false;
    else if (var.GetComputationalMode() == var.INVERT_MODE)
      {
	generalized = false;
	isherm = false;
      }
    
    if (generalized)
      EPSSetOperators(eps, stiff, mass);
    else
      EPSSetOperators(eps, stiff, NULL);

    if (isherm)
      {
	if (generalized)
	  EPSSetProblemType(eps, EPS_GHEP);
	else
	  EPSSetProblemType(eps, EPS_HEP);
      }
    else
      {
	if (generalized)
	  EPSSetProblemType(eps, EPS_PGNHEP);
	else
	  EPSSetProblemType(eps, EPS_NHEP);
      }

    EPSWhich which(EPS_LARGEST_MAGNITUDE);
    switch (var.GetTypeSorting())
      {
      case EigenProblem_Base<PetscScalar>::SORTED_REAL : which = EPS_LARGEST_REAL; break;
      case EigenProblem_Base<PetscScalar>::SORTED_IMAG : which = EPS_LARGEST_IMAGINARY; break;
      case EigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = EPS_LARGEST_MAGNITUDE; break;
      case EigenProblem_Base<PetscScalar>::SORTED_USER : which = EPS_WHICH_USER; break;
      }
    
    if (var.GetTypeSpectrum() == var.SMALL_EIGENVALUES)
      {
        switch (var.GetTypeSorting())
          {
          case EigenProblem_Base<PetscScalar>::SORTED_REAL : which = EPS_SMALLEST_REAL; break;
          case EigenProblem_Base<PetscScalar>::SORTED_IMAG : which = EPS_SMALLEST_IMAGINARY; break;
          case EigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = EPS_SMALLEST_MAGNITUDE; break;
          }
      }

    if ((var.GetComputationalMode() == var.REGULAR_MODE) &&
	(var.GetTypeSpectrum() == var.CENTERED_EIGENVALUES))
      {
	PetscScalar target = shiftr;
        switch (var.GetTypeSorting())
          {
          case EigenProblem_Base<PetscScalar>::SORTED_REAL : which = EPS_TARGET_REAL; break;
          case EigenProblem_Base<PetscScalar>::SORTED_IMAG : which = EPS_TARGET_IMAGINARY; target = shifti; break;
          case EigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = EPS_TARGET_MAGNITUDE; break;
          }
	
	EPSSetTarget(eps, target);
      }
    
    EPSSetWhichEigenpairs(eps, which);
    if (which == EPS_WHICH_USER)
      EPSSetEigenvalueComparison(eps, &EigenProblem_Base<PetscScalar>::GetComparisonEigenvalueSlepc, &var);
    
    double tol = var.GetStoppingCriterion();
    int nb_max_iter = var.GetNbMaximumIterations();
    EPSSetTolerances(eps, tol, nb_max_iter);
    EPSSetDimensions(eps, nev, ncv, PETSC_DEFAULT);
    EPSSetFromOptions(eps);
    
    // computing needed operators
    if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
      {
	if (var.DiagonalMass())
	  {
            // computation of M
            var.ComputeDiagonalMass();
	    
            // computation of M^{-1/2}
            var.FactorizeDiagonalMass();
          }
        else 
          {
            // computation of M for Cholesky factorisation
            var.ComputeMassForCholesky();
            
            // computation of Cholesky factorisation M = L L^T
            var.FactorizeCholeskyMass();
          }

	if (var.GetComputationalMode() != var.REGULAR_MODE)
	  var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);
      }
    else
      {
	if (var.GetComputationalMode() == var.INVERT_MODE)
          {
	    if (var.GetTypeSpectrum() != var.CENTERED_EIGENVALUES)
              {
                // computation and factorisation of mass matrix
                var.ComputeAndFactorizeStiffnessMatrix(one, zero);
                
                // computation of stiffness matrix
                var.ComputeStiffnessMatrix();
              }
            else
              {
                // computation and factorization of (K - sigma M)
                var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);

                // computation of M
                var.ComputeMassMatrix();
	      }
	  }
	else if (var.GetComputationalMode() == var.REGULAR_MODE)
          {
	    // factorization of the mass matrix
            var.ComputeAndFactorizeStiffnessMatrix(one, zero);
					 
            // computation of stiffness and mass matrix
            var.ComputeStiffnessMatrix();
            var.ComputeMassMatrix();

	    cout << "Case not implemented" << endl;
	    abort();
	  }
	else
	  {
	    cout << "Case not implemented" << endl;
	    abort();
	  }
      }

    // the eigenvalue problem is solved
    int ierr = EPSSolve(eps);
    if (ierr != 0)
      {
	cout << "Error during solution of eigensystem =  " << ierr << endl;
	abort();
      }

    if (print_level >= 4)
      {
	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
	EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD);
	PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
      }
    
    EPSConvergedReason reason;
    EPSGetConvergedReason(eps, &reason);
    if (reason < 0)
      {
	if (rank_proc == 0)
	  cout << "The solver did not converge" << endl;
        
        if (reason == EPS_DIVERGED_ITS)
          throw SolverMaximumIterationError("FindEigenvaluesSlepc", "Maximum number of iterations reached");
        else
          throw SolverDivergenceError("FindEigenvaluesSlepc", "The solver diverged");        
      }
    
    // eigenvalues and eigenvectors are extracted
    Vec Vr, Vi, Vr_next, Vi_next;
    Vector<PetscScalar> Vr_vec, Vi_vec, Vr_vec_next, Vi_vec_next;
    AllocatePetscVector(var.GetCommunicator(), Vr, n, nglob, Vr_vec);
    AllocatePetscVector(var.GetCommunicator(), Vi, n, nglob, Vi_vec);
    AllocatePetscVector(var.GetCommunicator(), Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(var.GetCommunicator(), Vi_next, n, nglob, Vi_vec_next);
    
    eigen_values.Reallocate(nev);
    lambda_imag.Reallocate(nev);    
    eigen_vectors.Reallocate(n, nev);
    int num = 0;
    PetscScalar Lr, Li, Lr_next, Li_next;
    bool eigen_pair = true;
    while (num < nev)
      {
	if (eigen_pair)
	  EPSGetEigenpair(eps, num, &Lr, &Li, Vr, Vi);

	if (num < nev-1)
	  EPSGetEigenpair(eps, num+1, &Lr_next, &Li_next, Vr_next, Vi_next);
	
	eigen_pair = PutEigenpairLapackForm(num, nev, Lr, Li, Vr_vec, Vi_vec,
					    Lr_next, Li_next, Vr_vec_next, Vi_vec_next,
					    eigen_values, lambda_imag, eigen_vectors);
	
	if (eigen_pair)
	  num += 2;
	else
	  {
	    Lr = Lr_next; Li = Li_next;
	    Vr_vec = Vr_vec_next; Vi_vec = Vi_vec_next;
	    num++;
	  }
      }

    Vr_vec.Nullify();
    Vi_vec.Nullify();
    Vr_vec_next.Nullify();
    Vi_vec_next.Nullify();

    // modifies eigenvalues and eigenvectors if needed
    ApplyScalingEigenvec(var, eigen_values, lambda_imag, eigen_vectors,
                         shiftr, shifti);
    
    // temporary objects are destroyed
    VecDestroy(&Vr);
    VecDestroy(&Vi);
    EPSDestroy(&eps);
    MatDestroy(&stiff);
    MatDestroy(&mass);

    // clears eigenproblem
    var.Clear();

  }
  
  /*****************************
   * Interface with PEP solver *
   *****************************/


  void ApplySpectralTransform(double shift, Vector<double>& eigen_values, Vector<double>& lambda_imag)
  {
    for (int i = 0; i < eigen_values.GetM(); i++)
      {
        if (lambda_imag(i) == 0.0)
          eigen_values(i) = shift + 1.0/eigen_values(i);
        else
          {
            complex<double> val(eigen_values(i), lambda_imag(i));
            complex<double> z = shift + 1.0/val;
            eigen_values(i) = real(z); lambda_imag(i) = imag(z);
          }
      }
  }

  void ApplySpectralTransform(complex<double> shift, Vector<complex<double> >& eigen_values, Vector<complex<double> >& lambda_imag)
  {
    for (int i = 0; i < eigen_values.GetM(); i++)
      eigen_values(i) = shift + 1.0/eigen_values(i);
  }

  
#ifdef SELDON_WITH_SLEPC_PEP
  void SetParametersSlepc(const SlepcParamPep& param, PEP& pep)
  {
    switch(param.GetEigensolverType())
      {
      case SlepcParamPep::TOAR :  PEPSetType(pep, PEPTOAR); break;
      case SlepcParamPep::STOAR :  PEPSetType(pep, PEPSTOAR); break;
      case SlepcParamPep::QARNOLDI :  PEPSetType(pep, PEPQARNOLDI); break;
      case SlepcParamPep::LINEAR :  PEPSetType(pep, PEPLINEAR); break;
      case SlepcParamPep::JD :  PEPSetType(pep, PEPJD); break;
      }        
  }
  
  struct MySlepcOperator
  {
    PolynomialEigenProblem_Base<PetscScalar>* var;
    int num_op;    
  };

  //! matrix-vector product y = op x
  PetscErrorCode MatMult_PepOperator(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    MySlepcOperator& op
      = *reinterpret_cast<MySlepcOperator*>(ctx);
    
    op.var->MltOperator(op.num_op, SeldonNoTrans, xvec, yvec);
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  //! matrix-vector product y = op^T x
  PetscErrorCode MatMult_PepOperatorTrans(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    MySlepcOperator& op
      = *reinterpret_cast<MySlepcOperator*>(ctx);
    
    op.var->MltOperator(op.num_op, SeldonTrans, xvec, yvec);
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  //! solving op y = x
  PetscErrorCode MatSolve_PepOperator(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    MySlepcOperator& op
      = *reinterpret_cast<MySlepcOperator*>(ctx);
    
    if (op.var->UseSpectralTransformation())
      op.var->SolveOperator(xvec, yvec);
    else
      op.var->SolveMass(xvec, yvec);
    
    op.var->IncrementLinearSolves();
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  //! solving op y = x
  PetscErrorCode MatSolveTrans_PepOperator(Mat mat, Vec x, Vec y)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  //! functions to compute eigenvalues with PEP solver of Slepc
  void FindEigenvaluesSlepc_(PolynomialEigenProblem_Base<PetscScalar>& var,
                             Vector<PetscScalar>& eigen_values,
                             Vector<PetscScalar>& lambda_imag,
                             Matrix<PetscScalar, General, ColMajor>& eigen_vectors)
  {
    // initializing of computation
    PetscScalar shift = var.GetCenterSpectrum();
    PetscScalar zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);

    // degree of polynom 
    int deg_pol = var.GetPolynomialDegree();
    int nev = var.GetNbAskedEigenvalues();
    int n = var.GetM();
    int nglob = var.GetGlobalM();
    
    // creation of shell matrices
    Vector<Mat> op(deg_pol+1);
    Vector<MySlepcOperator> my_op(deg_pol+1);
    for (int k = 0; k <= deg_pol; k++)
      {
        my_op(k).var = &var;
        my_op(k).num_op = k;
        MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
                       reinterpret_cast<void*>(&my_op(k)), &op(k));
        
        MatShellSetOperation(op(k), MATOP_MULT, (void(*)())MatMult_PepOperator);
        MatShellSetOperation(op(k), MATOP_MULT_TRANSPOSE, (void(*)())MatMult_PepOperatorTrans);        
      }

    // creation of the eigensolver
    PEP solver; ST st;
    PEPCreate(var.GetCommunicator(), &solver);
    
    PEPSetOperators(solver, deg_pol+1, op.GetData());
    PEPGetST(solver, &st);
    st->matsolve = (PetscErrorCode(*)(Mat, Vec, Vec))MatSolve_PepOperator;
    st->matsolve_trans = (PetscErrorCode(*)(Mat, Vec, Vec))MatSolveTrans_PepOperator;
    //st->D = NULL;
    STSetMatMode(st, ST_MATMODE_SHELL);
    //STSetType(st, STSINVERT);
    
    // sorting and selection of spectrum
    PEPWhich which(PEP_LARGEST_MAGNITUDE);
    switch (var.GetTypeSorting())
      {
      case PolynomialEigenProblem_Base<PetscScalar>::SORTED_REAL : which = PEP_LARGEST_REAL; break;
      case PolynomialEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = PEP_LARGEST_IMAGINARY; break;
      case PolynomialEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = PEP_LARGEST_MAGNITUDE; break;
      case PolynomialEigenProblem_Base<PetscScalar>::SORTED_USER : which = PEP_WHICH_USER; break;
      }
    
    if (var.GetTypeSpectrum() == var.SMALL_EIGENVALUES)
      {
        switch (var.GetTypeSorting())
          {
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_REAL : which = PEP_SMALLEST_REAL; break;
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = PEP_SMALLEST_IMAGINARY; break;
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = PEP_SMALLEST_MAGNITUDE; break;
          }
      }

    if ((!var.UseSpectralTransformation()) &&
        (var.GetTypeSpectrum() == var.CENTERED_EIGENVALUES))
      {
        switch (var.GetTypeSorting())
          {
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_REAL : which = PEP_TARGET_REAL; break;
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = PEP_TARGET_IMAGINARY; break;
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = PEP_TARGET_MAGNITUDE; break;
          }
	
	PEPSetTarget(solver, shift);
      }
    
    PEPSetWhichEigenpairs(solver, which);
    if (which == PEP_WHICH_USER)
      PEPSetEigenvalueComparison(solver, &PolynomialEigenProblem_Base<PetscScalar>::GetComparisonEigenvalueSlepc, &var);
    
    // tolerance and number of iterations
    double tol = var.GetStoppingCriterion();
    int nb_max_iter = var.GetNbMaximumIterations();
    PEPSetTolerances(solver, tol, nb_max_iter);    
    PEPSetDimensions(solver, nev, PETSC_DECIDE, PETSC_DECIDE);
    
    if (!var.UseSpectralTransformation())
      var.FactorizeMass();
    else
      {
        // binomial coefficients are computed
        Matrix<int> binom_coef(deg_pol+1, deg_pol+1);
        binom_coef.Zero();
        binom_coef(0, 0) = 1;
        binom_coef(1, 0) = 1; binom_coef(1, 1) = 1;
        for (int n = 2; n <= deg_pol; n++)
          {
            binom_coef(n, 0) = 1; binom_coef(n, n) = 1;
            for (int k = 1; k < n; k++)
              binom_coef(n, k) = binom_coef(n-1, k-1) + binom_coef(n-1, k);
          }
        
        // setting operators to compute
        Vector<PetscScalar> coef(deg_pol+1);
        for (int k = 0; k <= deg_pol; k++)
          {
            coef.Zero(); PetscScalar pow_shift = 1.0;
            for (int j = 0; j <= deg_pol-k; j++)
              {
                coef(j+k) = double(binom_coef(j+k, k))*pow_shift;
                pow_shift *= shift;
              }
            
            var.ComputeOperator(deg_pol-k, coef);
            if (k == 0)
              var.FactorizeOperator(coef);
          }
      }
    
    PEPSetScale(solver, PEP_SCALE_SCALAR, PETSC_DECIDE, PETSC_NULL, PETSC_NULL,
                PETSC_DECIDE, PETSC_DECIDE);
    
    // other parameters of PEP
    SetParametersSlepc(var.GetSlepcParameters(), solver);
    
    // the eigenvalue problem is solved
    int ierr = PEPSolve(solver);
    if (ierr != 0)
      {
	cout << "Error during solution of eigensystem =  " << ierr << endl;
	abort();
      }

    PEPConvergedReason reason;
    PEPGetConvergedReason(solver, &reason);
    if (reason < 0)
      {
        cout << "Failed to converged " << reason << endl;
        abort();
      }

    int print_level = var.GetPrintLevel();
    if (print_level >= 4)
      {
	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
	PEPErrorView(solver, PEP_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD);
	PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
      }

    int rank_proc = 0;
#ifdef SELDON_WITH_MPI
    MPI_Comm_rank(var.GetCommunicator(), &rank_proc);
#endif
    
    if ((print_level >= 1) && (rank_proc == 0))
      cout << "Number of linear solves = " << var.GetNbLinearSolves() << endl;
    
    // eigenvalues and eigenvectors are extracted
    Vec Vr, Vi, Vr_next, Vi_next;
    Vector<PetscScalar> Vr_vec, Vi_vec, Vr_vec_next, Vi_vec_next;
    AllocatePetscVector(var.GetCommunicator(), Vr, n, nglob, Vr_vec);
    AllocatePetscVector(var.GetCommunicator(), Vi, n, nglob, Vi_vec);
    AllocatePetscVector(var.GetCommunicator(), Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(var.GetCommunicator(), Vi_next, n, nglob, Vi_vec_next);
    
    eigen_values.Reallocate(nev);
    lambda_imag.Reallocate(nev);    
    eigen_vectors.Reallocate(n, nev);
    int num = 0;
    PetscScalar Lr, Li, Lr_next, Li_next;
    bool eigen_pair = true;
    while (num < nev)
      {
	if (eigen_pair)
	  PEPGetEigenpair(solver, num, &Lr, &Li, Vr, Vi);

	if (num < nev-1)
	  PEPGetEigenpair(solver, num+1, &Lr_next, &Li_next, Vr_next, Vi_next);
	
	eigen_pair = PutEigenpairLapackForm(num, nev, Lr, Li, Vr_vec, Vi_vec,
					    Lr_next, Li_next, Vr_vec_next, Vi_vec_next,
					    eigen_values, lambda_imag, eigen_vectors);
	
	if (eigen_pair)
	  num += 2;
	else
	  {
	    Lr = Lr_next; Li = Li_next;
	    Vr_vec = Vr_vec_next; Vi_vec = Vi_vec_next;
	    num++;
	  }
      }
    
    Vr_vec.Nullify();
    Vi_vec.Nullify();
    Vr_vec_next.Nullify();
    Vi_vec_next.Nullify();
    
    if (var.UseSpectralTransformation())
      ApplySpectralTransform(shift, eigen_values, lambda_imag);

    // temporary objects are destroyed
    VecDestroy(&Vr);
    VecDestroy(&Vi);
    PEPDestroy(&solver);
    for (int k = 0; k <= deg_pol; k++)
      MatDestroy(&op(k));
  }
#endif

  
  /*****************************
   * Interface with NEP solver *
   *****************************/

  
  void SetParametersSlepc(const SlepcParamNep& param, NEP& nep)
  {
    switch(param.GetEigensolverType())
      {
      case SlepcParamNep::RII :  NEPSetType(nep, NEPRII); break;
      case SlepcParamNep::SLP :  NEPSetType(nep, NEPSLP); break;
      case SlepcParamNep::NARNOLDI :  NEPSetType(nep, NEPNARNOLDI); break;
      case SlepcParamNep::CISS :  NEPSetType(nep, NEPCISS); break;
      case SlepcParamNep::INTERPOL :  NEPSetType(nep, NEPINTERPOL); break;
      case SlepcParamNep::NLEIGS :  NEPSetType(nep, NEPNLEIGS); break;
      }        

    //RG rg;
    //NEPGetRG(nep, &rg);
    //RGSetType(rg, RGINTERVAL);
    //RGIntervalSetEndpoints(rg, -10.0, -2.0, 0, 0);
  }
  
  struct NepSlepcOperator
  {
    bool jacobian;
    PetscScalar L;
    NonLinearEigenProblem_Base<PetscScalar>* var;    
  };

  PetscErrorCode FormFunctionNEP(NEP nep, PetscScalar lambda, Mat A, Mat B, void* ctx)
  {
    void* ctxF;
    MatShellGetContext(A, &ctxF);
    
    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctx);
    
    op.L = lambda;
    op.var->ComputeOperator(op.L);
    
    return 0;
  }

  PetscErrorCode FormJacobianNEP(NEP nep, PetscScalar lambda, Mat A, void* ctx)
  {
    void* ctxF;
    MatShellGetContext(A, &ctxF);
    
    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctx);
    
    op.L = lambda;
    op.var->ComputeJacobian(op.L);
    
    return 0;
  }

  //! matrix-vector product y = op x
  PetscErrorCode MatMult_NepOperator(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctx);
    
    if (op.jacobian)
      op.var->MltJacobian(op.L, SeldonNoTrans, xvec, yvec);
    else
      op.var->MltOperator(op.L, SeldonNoTrans, xvec, yvec);
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  //! matrix-vector product y = op^T x
  PetscErrorCode MatMultTrans_NepOperator(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctx);
    
    if (op.jacobian)
      op.var->MltJacobian(op.L, SeldonTrans, xvec, yvec);
    else
      op.var->MltOperator(op.L, SeldonTrans, xvec, yvec);
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  
  //! functions to compute eigenvalues with NEP solver of Slepc
  void FindEigenvaluesSlepc_(NonLinearEigenProblem_Base<PetscScalar>& var,
                             Vector<PetscScalar>& eigen_values,
                             Vector<PetscScalar>& lambda_imag,
                             Matrix<PetscScalar, General, ColMajor>& eigen_vectors)
  {
    // initializing of computation
    PetscScalar shift = var.GetCenterSpectrum();
    PetscScalar zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    // dimensions
    int nev = var.GetNbAskedEigenvalues();
    int n = var.GetM();
    int nglob = var.GetGlobalM();
    
    // creation of shell matrices
    Mat EvalF, EvalJacob;
    NepSlepcOperator operT;
    operT.jacobian = false; operT.var = &var;
    MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
                   reinterpret_cast<void*>(&operT), &EvalF);
    
    MatShellSetOperation(EvalF, MATOP_MULT, (void(*)())MatMult_NepOperator);
    MatShellSetOperation(EvalF, MATOP_MULT_TRANSPOSE, (void(*)())MatMultTrans_NepOperator); 

    NepSlepcOperator operTp;
    operTp.jacobian = true; operT.var = &var;
    MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
                   reinterpret_cast<void*>(&operTp), &EvalJacob);
    
    MatShellSetOperation(EvalJacob, MATOP_MULT, (void(*)())MatMult_NepOperator);
    MatShellSetOperation(EvalJacob, MATOP_MULT_TRANSPOSE, (void(*)())MatMultTrans_NepOperator);        
    
    // creation of the eigensolver
    NEP solver;
    NEPCreate(var.GetCommunicator(), &solver);
    
    NEPSetFunction(solver, EvalF, EvalF, FormFunctionNEP, NULL);
    NEPSetJacobian(solver, EvalJacob, FormJacobianNEP, NULL);
    
    // sorting and selection of spectrum
    NEPWhich which(NEP_LARGEST_MAGNITUDE);
    switch (var.GetTypeSorting())
      {
      case NonLinearEigenProblem_Base<PetscScalar>::SORTED_REAL : which = NEP_LARGEST_REAL; break;
      case NonLinearEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = NEP_LARGEST_IMAGINARY; break;
      case NonLinearEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = NEP_LARGEST_MAGNITUDE; break;
      case NonLinearEigenProblem_Base<PetscScalar>::SORTED_USER : which = NEP_WHICH_USER; break;
      }
    
    if (var.GetTypeSpectrum() == var.SMALL_EIGENVALUES)
      {
        switch (var.GetTypeSorting())
          {
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_REAL : which = NEP_SMALLEST_REAL; break;
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = NEP_SMALLEST_IMAGINARY; break;
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = NEP_SMALLEST_MAGNITUDE; break;
          }
      }

    if ((!var.UseSpectralTransformation()) &&
        (var.GetTypeSpectrum() == var.CENTERED_EIGENVALUES))
      {
        switch (var.GetTypeSorting())
          {
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_REAL : which = NEP_TARGET_REAL; break;
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = NEP_TARGET_IMAGINARY; break;
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = NEP_TARGET_MAGNITUDE; break;
          }
	
	NEPSetTarget(solver, shift);
      }
    
    NEPSetWhichEigenpairs(solver, which);
    if (which == NEP_WHICH_USER)
      NEPSetEigenvalueComparison(solver, &NonLinearEigenProblem_Base<PetscScalar>::GetComparisonEigenvalueSlepc, &var);
    
    // tolerance and number of iterations
    double tol = var.GetStoppingCriterion();
    int nb_max_iter = var.GetNbMaximumIterations();
    NEPSetTolerances(solver, tol, nb_max_iter);    
    NEPSetDimensions(solver, nev, PETSC_DEFAULT, PETSC_DEFAULT);
    
    // other parameters of NEP
    SetParametersSlepc(var.GetSlepcParameters(), solver);
    
    // the eigenvalue problem is solved
    int ierr = NEPSolve(solver);
    if (ierr != 0)
      {
	cout << "Error during solution of eigensystem =  " << ierr << endl;
	abort();
      }
    
    int print_level = 1;
    if (print_level >= 4)
      {
	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
	NEPErrorView(solver, NEP_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD);
	PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
      }
    
    // eigenvalues and eigenvectors are extracted
    Vec Vr, Vi, Vr_next, Vi_next;
    Vector<PetscScalar> Vr_vec, Vi_vec, Vr_vec_next, Vi_vec_next;
    AllocatePetscVector(var.GetCommunicator(), Vr, n, nglob, Vr_vec);
    AllocatePetscVector(var.GetCommunicator(), Vi, n, nglob, Vi_vec);
    AllocatePetscVector(var.GetCommunicator(), Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(var.GetCommunicator(), Vi_next, n, nglob, Vi_vec_next);
    
    eigen_values.Reallocate(nev);
    lambda_imag.Reallocate(nev);    
    eigen_vectors.Reallocate(n, nev);
    int num = 0;
    PetscScalar Lr, Li, Lr_next, Li_next;
    bool eigen_pair = true;
    while (num < nev)
      {
	if (eigen_pair)
	  NEPGetEigenpair(solver, num, &Lr, &Li, Vr, Vi);

	if (num < nev-1)
	  NEPGetEigenpair(solver, num+1, &Lr_next, &Li_next, Vr_next, Vi_next);
	
	eigen_pair = PutEigenpairLapackForm(num, nev, Lr, Li, Vr_vec, Vi_vec,
					    Lr_next, Li_next, Vr_vec_next, Vi_vec_next,
					    eigen_values, lambda_imag, eigen_vectors);
	
	if (eigen_pair)
	  num += 2;
	else
	  {
	    Lr = Lr_next; Li = Li_next;
	    Vr_vec = Vr_vec_next; Vi_vec = Vi_vec_next;
	    num++;
	  }
      }
    
    Vr_vec.Nullify();
    Vi_vec.Nullify();
    Vr_vec_next.Nullify();
    Vi_vec_next.Nullify();
    
    if (var.UseSpectralTransformation())
      ApplySpectralTransform(shift, eigen_values, lambda_imag);

    // temporary objects are destroyed
    VecDestroy(&Vr);
    VecDestroy(&Vi);
    NEPDestroy(&solver);
    MatDestroy(&EvalF);
    MatDestroy(&EvalJacob);
  }
  
#ifdef SELDON_WITH_VIRTUAL
  void FindEigenvaluesSlepc(EigenProblem_Base<PetscScalar>& var,
			    Vector<PetscScalar>& eigen_values,
			    Vector<PetscScalar>& lambda_imag,
			    Matrix<PetscScalar, General, ColMajor>& eigen_vectors)
  {
    FindEigenvaluesSlepc_(var, eigen_values, lambda_imag, eigen_vectors);
  }  

  void FindEigenvaluesSlepc(PolynomialEigenProblem_Base<Petsc_Scalar>& var,
                            Vector<Petsc_Scalar>& eigen_values,
                            Vector<Petsc_Scalar>& lambda_imag,
                            Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors)
  {
    FindEigenvaluesSlepc_(var, eigen_values, lambda_imag, eigen_vectors);
  }

  void FindEigenvaluesSlepc(NonLinearEigenProblem_Base<Petsc_Scalar>& var,
                            Vector<Petsc_Scalar>& eigen_values,
                            Vector<Petsc_Scalar>& lambda_imag,
                            Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors)
  {
    FindEigenvaluesSlepc_(var, eigen_values, lambda_imag, eigen_vectors);
  }
#else
  template<class EigenProblem, class T, class Allocator1,
           class Allocator2, class Allocator3>
  void FindEigenvaluesSlepc(EigenProblem& var,
			    Vector<T, VectFull, Allocator1>& eigen_values,
			    Vector<T, VectFull, Allocator2>& lambda_imag,
			    Matrix<T, General, ColMajor, Allocator3>& eigen_vectors)
  {
    cout << "Recompile with SELDON_WITH_VIRTUAL" << endl;
    abort();
  }
#endif
  
}

#define SELDON_FILE_SLEPC_HXX
#endif
