#ifndef SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_CXX

#include "PolynomialEigenvalueSolver.hxx"

namespace Seldon
{
  SlepcParamPep::SlepcParamPep()
  {
    type_solver = TOAR;
  }

  int SlepcParamPep::GetEigensolverType() const
  {
    return type_solver;
  }

  void SlepcParamPep::SetEigensolverType(int type)
  {
    type_solver = type;
  }
  
  /*******************************
   * PolynomialEigenProblem_Base *
   *******************************/

  //! default constructor
  template<class T>
  PolynomialEigenProblem_Base<T>::PolynomialEigenProblem_Base()
  {
    nb_eigenvalues_wanted = 0;
    type_spectrum_wanted = LARGE_EIGENVALUES;
    SetComplexZero(shift);
    use_spectral_transfo = false;
    type_sort_eigenvalues = SORTED_MODULUS;
    n_ = 0; nglob = 0;
    pol_degree = 0;
    
#ifdef SELDON_WITH_MPI
    // for parallel execution, default communicator : all the processors
    comm = MPI_COMM_WORLD;
#endif
    
    compar_eigenval = NULL;
    nb_linear_solves = 0;
    SetPrintLevel(0);
    diagonal_mass = false;

    stopping_criterion = 1e-6;
    nb_maximum_iterations = 1000;
  }

  //! destructor
  template<class T>
  PolynomialEigenProblem_Base<T>::~PolynomialEigenProblem_Base()
  {
    compar_eigenval = NULL;
  }

  //! initialisation of the size of the eigenvalue problem
  template<class T>
  void PolynomialEigenProblem_Base<T>::Init(int n)
  {
    n_ = n;
    nb_linear_solves = 0;
    
    // counting the size of the global system for parallel computation
    nglob = n;
    
#ifdef SELDON_WITH_MPI
    MPI_Allreduce(&n, &nglob, 1, MPI_INTEGER, MPI_SUM, comm);    
#endif

  }

#ifdef SELDON_WITH_MPI
  //! changes MPI communicator
  template<class T>
  void PolynomialEigenProblem_Base<T>::SetCommunicator(const MPI_Comm& comm_)
  {
    comm = comm_;
  }
  
  //! returns the MPI communicator
  template<class T>
  MPI_Comm& PolynomialEigenProblem_Base<T>::GetCommunicator()
  {
    return comm;
  }
#endif

  //! returns the number of asked eigenvalues
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetNbAskedEigenvalues() const
  {
    return nb_eigenvalues_wanted;
  }
  
  //! sets the number of eigenvalues to compute
  template<class T>
  void PolynomialEigenProblem_Base<T>::SetNbAskedEigenvalues(int n)
  {
    nb_eigenvalues_wanted = n;
  }

  //! returns true if a spectral transformation has to be used
  template<class T>
  bool PolynomialEigenProblem_Base<T>::UseSpectralTransformation() const
  {
    return use_spectral_transfo;
  }

  //! enables a spectral transformation
  template<class T>
  void PolynomialEigenProblem_Base<T>::SetSpectralTransformation(bool t)
  {
    use_spectral_transfo = t;
  }
  
  //! returns the type of spectrum desired
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetTypeSpectrum() const
  {
    return type_spectrum_wanted;
  }
  
  //! returns the center of the spectrum (if selected)
  template<class T>
  T PolynomialEigenProblem_Base<T>::GetCenterSpectrum() const
  {
    return shift;
  }

  //! returns how the eigenvalues are sorted
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetTypeSorting() const
  {
    return type_sort_eigenvalues;
  }
  
  //! sets the spectrum and sorting 
  template<class T>
  void PolynomialEigenProblem_Base<T>
  ::SetTypeSpectrum(int type, const T& val, int type_sort)
  {
    type_spectrum_wanted = type;
    shift = val;
    type_sort_eigenvalues = type_sort;
  }

  //! returns object storing slepc parameters
  template<class T>
  SlepcParamPep& PolynomialEigenProblem_Base<T>::GetSlepcParameters()
  {
    return slepc_param;
  }
  
  //! sets how the eigenvalues are sorted
  template<class T>
  void PolynomialEigenProblem_Base<T>::SetUserComparisonClass(EigenvalueComparisonClass<T>* ev)
  {
    compar_eigenval = ev;
  }
  
  //! sets the stopping criterion
  template<class T>
  void PolynomialEigenProblem_Base<T>::SetStoppingCriterion(Treal eps)
  {
    stopping_criterion = eps;
  }
  
  //! returns the stopping criterion
  template<class T>
  typename ClassComplexType<T>::Treal PolynomialEigenProblem_Base<T>::GetStoppingCriterion() const
  {
    return stopping_criterion;
  }
  
  //! sets the maximal number of iterations
  template<class T>
  void PolynomialEigenProblem_Base<T>::SetNbMaximumIterations(int n)
  {
    nb_maximum_iterations = n;
  }
  
  //! returns the maximal number of iterations
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetNbMaximumIterations() const
  {
    return nb_maximum_iterations;
  }

  //! returns the number of rows
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetM() const
  {
    return n_;
  }
  
  //! returns the global number of rows
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetGlobalM() const
  {
    return nglob;
  }

  //! returns the number of columns
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetN() const
  {
    return n_;
  }

  //! returns the polynomial degree
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetPolynomialDegree() const
  {
    return pol_degree;
  }
  
  //! returns the number of linear solves
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetNbLinearSolves() const
  {
    return nb_linear_solves;
  }

  //! increments the number of linear solves
  template<class T>    
  void PolynomialEigenProblem_Base<T>::IncrementLinearSolves()
  {
    nb_linear_solves++;

#ifdef SELDON_WITH_MPI
    int rank; MPI_Comm_rank(comm, & rank);
#else
    int rank(0);
#endif

    if (nb_linear_solves%display_every == 0)
      if (rank == 0)
        cout<<" Iteration number " << nb_linear_solves << endl;
  }

  
  template<class T>    
  void PolynomialEigenProblem_Base<T>::SetPrintLevel(int lvl)
  {
    print_level = lvl;
    display_every = 1000000;
    if (lvl == 1)
      display_every = 1000;
    else if (lvl == 2)
      display_every = 100;
    else if (lvl >= 3)
      display_every = 10;
  }
  

  template<class T>    
  int PolynomialEigenProblem_Base<T>::GetPrintLevel() const
  {
    return print_level;
  }

  
  //! sets a diagonal mass
  template<class T>    
  void PolynomialEigenProblem_Base<T>::SetDiagonalMass(bool diag)
  {
    diagonal_mass = diag;
  }

  //! returns true if the mass is diagonal
  template<class T>
  bool PolynomialEigenProblem_Base<T>::DiagonalMass()
  {
    return diagonal_mass;
  }

#ifdef SELDON_WITH_SLEPC
  //! compares two eigenvalues
  template<class T> template<class T0>
  Petsc_Error_Code PolynomialEigenProblem_Base<T>::
  GetComparisonEigenvalueSlepcGen(Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, T0,
                                  Petsc_Int*, void*)
  {
    cout << "Incompatibles types between T and PetscScalar" << endl;
    abort();
    return 0;
  }

  
  //! compares two eigenvalues
  template<class T>
  Petsc_Error_Code PolynomialEigenProblem_Base<T>::
  GetComparisonEigenvalueSlepcGen(Petsc_Scalar Lr, Petsc_Scalar Li, Petsc_Scalar Lr2, Petsc_Scalar Li2, Petsc_Scalar,
                                  Petsc_Int* res, void* ctx)
  {
    PolynomialEigenProblem_Base<T>& var = *reinterpret_cast<PolynomialEigenProblem_Base<T>* >(ctx);
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

  //! compares two eigenvalues
  template<class T>
  Petsc_Error_Code PolynomialEigenProblem_Base<T>::
  GetComparisonEigenvalueSlepc(Petsc_Scalar Lr, Petsc_Scalar Li, Petsc_Scalar Lr2, Petsc_Scalar Li2,
                               Petsc_Int* res, void* ctx)
  {
    T z; SetComplexZero(z);
    GetComparisonEigenvalueSlepcGen(Lr, Li, Lr2, Li2, z, res, ctx);
    return 0;
  }
#endif

  //! to overload
  template<class T>
  void PolynomialEigenProblem_Base<T>::ComputeOperator(int num, const Vector<T>& coef)
  {
    cout << "ComputeOperator not overloaded" << endl;
    abort();
  }

  //! to overload
  template<class T>
  void PolynomialEigenProblem_Base<T>::MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "MltOperator not overloaded" << endl;
    abort();
  }

  //! to overload
  template<class T>    
  void PolynomialEigenProblem_Base<T>::FactorizeMass()
  {
    cout << "FactorizeMass not overloaded" << endl;
    abort();
  }
  
  //! to overload for non-diagonal mass
  template<class T>    
  void PolynomialEigenProblem_Base<T>::SolveMass(const Vector<T>& x, Vector<T>& y)
  {
    if (DiagonalMass())
      for (int i = 0; i < x.GetM(); i++)
        y(i) = x(i)*this->invDiag(i);
    else
      {
        cout << "SolveMass not overloaded" << endl;
        abort();
      }
  }
  
  template<class T>
  void PolynomialEigenProblem_Base<T>::FactorizeOperator(const Vector<T>& coef)
  {
    cout << "FactorizeOperator not overloaded" << endl;
    abort();
  }

  template<class T>
  void PolynomialEigenProblem_Base<T>::SolveOperator(const Vector<T>& X, Vector<T>& Y)
  {
    cout << "SolveOperator not overloaded" << endl;
    abort(); 
  }

  /**************************
   * PolynomialEigenProblem *
   **************************/
  
  //! inits the operators of the polynomial
  template<class T>
  void PolynomialEigenProblem<T>::InitMatrix(const Vector<VirtualMatrix<T>* >& op, int n)
  {
    if (n == -1)
      this->Init(op(0)->GetM());
    else
      this->Init(n);
    
    list_op = op;
    this->pol_degree = list_op.GetM()-1;
  }
  
  //! computes the operator with coefficients stored in coef
  template<class T>
  void PolynomialEigenProblem<T>::ComputeOperator(int num, const Vector<T>& coef)
  {
    // default implementation : we store the coefficients, no matrix stored
    if (num >= list_coef.GetM())
      list_coef.Resize(num+1);
    
    list_coef(num) = coef;
  }
  
  //! Computes Y = A X where A is the operator num
  template<class T>
  void PolynomialEigenProblem<T>::MltOperator(int num, const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y)
  {
    if (list_op.GetM() <= 0)
      return;
    
    if (list_coef.GetM() == 0)
      {
        // no coeffients, the operator num is stored in list_op(num)
        list_op(num)->MltVector(trans, X, Y);
        return;
      }

    // coefficients, the operator num is a linear combination of stored operators
    T zero; SetComplexZero(zero);
    T one; SetComplexOne(one);
    if (list_coef(num)(0) != zero)
      {
        list_op(0)->MltVector(trans, X, Y);
        Mlt(list_coef(num)(0), Y);
      }
    else
      Y.Zero();
    
    for (int i = 1; i < list_op.GetM(); i++)
      if (list_coef(num)(i) != zero)
        list_op(i)->MltAddVector(list_coef(num)(i), trans, X, one, Y);
  }

  /*******************************
   * PolynomialDenseEigenProblem *
   *******************************/

  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>::InitMatrix(const Vector<Matrix<T, Prop, Storage>* >& op)
  {
    Vector<VirtualMatrix<T>* > op0(op.GetM());
    for (int i = 0; i < op0.GetM(); i++)
      op0(i) = op(i);
    
    PolynomialEigenProblem<T>::InitMatrix(op0);
    list_mat = op;
  }
  
  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>::FactorizeMass()
  {    
    if (this->DiagonalMass())
      {
        T one; SetComplexOne(one);
        this->invDiag.Reallocate(this->n_);
        Matrix<T, Prop, Storage>& M = *list_mat(this->pol_degree);
        for (int i = 0; i < this->n_; i++)
          this->invDiag(i) = one / M(i, i);
      }
    else
      {
        mat_lu = *list_mat(this->pol_degree);
        GetLU(mat_lu, pivot);
      }
  }

  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>::SolveMass(const Vector<T>& x, Vector<T>& y)
  {
    if (this->DiagonalMass())
      for (int i = 0; i < x.GetM(); i++)
        y(i) = x(i)*this->invDiag(i);
    else
      {
        y = x;
        SolveLU(mat_lu, pivot, y);
      }
  }

  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>::FactorizeOperator(const Vector<T>& coef)
  {
    T zero; SetComplexZero(zero);
    if (coef(0) == zero)
      {
        mat_lu.Reallocate(this->n_, this->n_);
        mat_lu.Zero();
      }
    else
      {
        mat_lu = *list_mat(0);
        Mlt(coef(0), mat_lu);
      }

    for (int k = 1; k < list_mat.GetM(); k++)
      if (coef(k) != zero)
        Add(coef(k), *list_mat(k), mat_lu);
    
    GetLU(mat_lu, pivot);
  }

  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>::SolveOperator(const Vector<T>& X, Vector<T>& Y)
  {
    Y = X;
    SolveLU(mat_lu, pivot, Y);
  }
  
  
  /*********************************
   *  PolynomialSparseEigenProblem *
   *********************************/

  template<class T, class MatStiff, class MatMass>
  PolynomialSparseEigenProblem<T, MatStiff, MatMass>::PolynomialSparseEigenProblem()
  {
    Mh = NULL;
    ProcSharingRows = NULL;
    SharingRowNumbers = NULL;
    nodl_scalar_ = nb_unknowns_scal_ = 0;
    nloc = 0;
  }
    
  template<class T, class MatStiff, class MatMass>
  int PolynomialSparseEigenProblem<T, MatStiff, MatMass>::RetrieveLocalNumbers(MatStiff& K)
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
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::InitMatrix(Vector<MatStiff*>& K, MatMass& M)
  {
    int n = RetrieveLocalNumbers(*K(0));
    distributed = false;
    if (n >= 0)
      distributed = true;
    
    Vector<VirtualMatrix<T>* > list_op(K.GetM() + 1);
    for (int i = 0; i < K.GetM(); i++)
      list_op(i) = K(i);
    
    list_op(K.GetM()) = &M;
    PolynomialEigenProblem<T>::InitMatrix(list_op, n);
    
    Kh.Reallocate(K.GetM());
    for (int i = 0; i < K.GetM(); i++)
      Kh(i) = K(i);
    
    Mh = &M;
  }
  
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::FactorizeMass()
  {
    if (this->DiagonalMass())
      {
        T one; SetComplexOne(one);
        this->invDiag.Reallocate(nloc);
        for (int i = 0; i < nloc; i++)
          this->invDiag(i) = (*Mh)(i, i);
        
	// D is assembled for distributed matrices
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<T> M(this->invDiag);
	    AssembleVector(M, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 15);
            
	    this->invDiag.Reallocate(local_col_numbers.GetM());
	    for (int i = 0; i < local_col_numbers.GetM(); i++)
	      this->invDiag(i) = one / M(local_col_numbers(i));
#endif
	  }
        else
          for (int i = 0; i < nloc; i++)
            this->invDiag(i) = one / this->invDiag(i);
      }
    else
      {
        mat_lu.Factorize(*Mh);
      }
  }
  
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::SolveMass(const Vector<T>& x, Vector<T>& y)
  {
    if (this->DiagonalMass())
      for (int i = 0; i < x.GetM(); i++)
        y(i) = x(i)*this->invDiag(i);
    else
      this->SolveOperator(x, y);

  }
  
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::FactorizeOperator(const Vector<T>& coef)
  {
    MatStiff A;
    T zero; SetComplexZero(zero);
    A = *Kh(0);
    Mlt(coef(0), A);
    for (int k = 1; k < Kh.GetM(); k++)
      if (coef(k) != zero)
        Add(coef(k), *Kh(k), A);
    
    Add(coef(Kh.GetM()), *Mh, A);
    
    mat_lu.Factorize(A);
  }

  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::SolveOperator(const Vector<T>& x, Vector<T>& y)
  {
    Vector<T> X(nloc);
    if (distributed)
      {
        X.Zero();
        for (int i = 0; i < this->n_; i++)
          X(local_col_numbers(i)) = x(i);                
      }
    else
      X = x;
    
    mat_lu.Solve(SeldonNoTrans, X, false);  
    if (distributed)
      {
        for (int i = 0; i < this->n_; i++)
          y(i) = X(local_col_numbers(i));
      }
    else
      y = X;
  }

  //! Computes Y = A X where A is the operator num
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::MltOperator(int num, const SeldonTranspose& trans, const Vector<T>& x, Vector<T>& y)
  {
    if (Kh.GetM() <= 0)
      return;
    
    Vector<T> X, Y;
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
        X.Reallocate(nloc); Y.Reallocate(nloc);
        X.Zero();
        for (int i = 0; i < this->n_; i++)
          X(local_col_numbers(i)) = x(i);
        
        AssembleVector(X, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
          this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
        
        if (this->list_coef.GetM() == 0)
          Kh(num)->MltVector(trans, X, Y);
        else
          {
            // coefficients, the operator num is a linear combination of stored operators
            T zero; SetComplexZero(zero);
            T one; SetComplexOne(one);
            if (this->list_coef(num)(0) != zero)
              {
                Kh(0)->MltVector(trans, X, Y);
                Mlt(this->list_coef(num)(0), Y);
              }
            else
              Y.Zero();
            
            for (int i = 1; i < this->list_op.GetM(); i++)
              if (this->list_coef(num)(i) != zero)
                this->list_op(i)->MltAddVector(this->list_coef(num)(i), trans, X, one, Y);
          }
        
        for (int i = 0; i < this->n_; i++)
          y(i) = Y(local_col_numbers(i));
#endif
      }
    else
      PolynomialEigenProblem<T>::MltOperator(num, trans, x, y);
  }
  
}
 
#define SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_CXX
#endif
