#ifndef SELDON_FILE_NON_LINEAR_EIGENVALUE_SOLVER_CXX

#include "NonLinearEigenvalueSolver.hxx"

namespace Seldon
{
  SlepcParamNep::SlepcParamNep()
  {
    type_solver = NLEIGS;
  }

  int SlepcParamNep::GetEigensolverType() const
  {
    return type_solver;
  }

  void SlepcParamNep::SetEigensolverType(int type)
  {
    type_solver = type;
  }
  
  /******************************
   * NonLinearEigenProblem_Base *
   ******************************/
  
  //! default constructor
  template<class T>
  NonLinearEigenProblem_Base<T>::NonLinearEigenProblem_Base()
  {
    nb_eigenvalues_wanted = 0;
    type_spectrum_wanted = LARGE_EIGENVALUES;
    SetComplexZero(shift);
    use_spectral_transfo = false;
    type_sort_eigenvalues = SORTED_MODULUS;
    n_ = 0; nglob = 0;
    
#ifdef SELDON_WITH_MPI
    // for parallel execution, default communicator : all the processors
    comm = MPI_COMM_WORLD;
#endif
    
    compar_eigenval = NULL;
    nb_prod = 0;
    
    stopping_criterion = 1e-6;
    nb_maximum_iterations = 1000;
  }

  //! destructor
  template<class T>
  NonLinearEigenProblem_Base<T>::~NonLinearEigenProblem_Base()
  {
    compar_eigenval = NULL;
  }

  //! initialisation of the size of the eigenvalue problem
  template<class T>
  void NonLinearEigenProblem_Base<T>::Init(int n)
  {
    n_ = n;
    nb_prod = 0;
    
    // counting the size of the global system for parallel computation
    nglob = n;
    
#ifdef SELDON_WITH_MPI
    MPI_Allreduce(&n, &nglob, 1, MPI_INTEGER, MPI_SUM, comm);    
#endif

  }

#ifdef SELDON_WITH_MPI
  //! changes MPI communicator
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetCommunicator(const MPI_Comm& comm_)
  {
    comm = comm_;
  }
  
  //! returns the MPI communicator
  template<class T>
  MPI_Comm& NonLinearEigenProblem_Base<T>::GetCommunicator()
  {
    return comm;
  }
#endif

  //! returns the number of asked eigenvalues
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetNbAskedEigenvalues() const
  {
    return nb_eigenvalues_wanted;
  }
  
  //! sets the number of eigenvalues to compute
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetNbAskedEigenvalues(int n)
  {
    nb_eigenvalues_wanted = n;
  }

  //! returns true if a spectral transformation has to be used
  template<class T>
  bool NonLinearEigenProblem_Base<T>::UseSpectralTransformation() const
  {
    return use_spectral_transfo;
  }

  //! enables a spectral transformation
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetSpectralTransformation(bool t)
  {
    use_spectral_transfo = t;
  }
  
  //! returns the type of spectrum desired
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetTypeSpectrum() const
  {
    return type_spectrum_wanted;
  }
  
  //! returns the center of the spectrum (if selected)
  template<class T>
  T NonLinearEigenProblem_Base<T>::GetCenterSpectrum() const
  {
    return shift;
  }

  //! returns how the eigenvalues are sorted
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetTypeSorting() const
  {
    return type_sort_eigenvalues;
  }
  
  //! sets the spectrum and sorting 
  template<class T>
  void NonLinearEigenProblem_Base<T>
  ::SetTypeSpectrum(int type, const T& val, int type_sort)
  {
    type_spectrum_wanted = type;
    shift = val;
    type_sort_eigenvalues = type_sort;
  }

  //! returns object storing slepc parameters
  template<class T>
  SlepcParamNep& NonLinearEigenProblem_Base<T>::GetSlepcParameters()
  {
    return slepc_param;
  }
  
  //! sets how the eigenvalues are sorted
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetUserComparisonClass(EigenvalueComparisonClass<T>* ev)
  {
    compar_eigenval = ev;
  }
  
  //! sets the stopping criterion
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetStoppingCriterion(Treal eps)
  {
    stopping_criterion = eps;
  }
  
  //! returns the stopping criterion
  template<class T>
  typename ClassComplexType<T>::Treal NonLinearEigenProblem_Base<T>::GetStoppingCriterion() const
  {
    return stopping_criterion;
  }
  
  //! sets the maximal number of iterations
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetNbMaximumIterations(int n)
  {
    nb_maximum_iterations = n;
  }
  
  //! returns the maximal number of iterations
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetNbMaximumIterations() const
  {
    return nb_maximum_iterations;
  }

  //! returns the number of rows
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetM() const
  {
    return n_;
  }
  
  //! returns the global number of rows
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetGlobalM() const
  {
    return nglob;
  }

  //! returns the number of columns
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetN() const
  {
    return n_;
  }

  //! returns the number of operator evaluation
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetNbMatrixVectorProducts() const
  {
    return nb_prod;
  }

  //! increments the number of operator evaluation
  template<class T>    
  void NonLinearEigenProblem_Base<T>::IncrementProdMatVect()
  {
    nb_prod++;
  }
  
#ifdef SELDON_WITH_SLEPC
  //! compares two eigenvalues
  template<class T> template<class T0>
  Petsc_Error_Code NonLinearEigenProblem_Base<T>::
  GetComparisonEigenvalueSlepcGen(Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, Petsc_Scalar, T0,
                                  Petsc_Int*, void*)
  {
    cout << "Incompatibles types between T and PetscScalar" << endl;
    abort();
    return 0;
  }

  
  //! compares two eigenvalues
  template<class T>
  Petsc_Error_Code NonLinearEigenProblem_Base<T>::
  GetComparisonEigenvalueSlepcGen(Petsc_Scalar Lr, Petsc_Scalar Li, Petsc_Scalar Lr2, Petsc_Scalar Li2, Petsc_Scalar,
                                  Petsc_Int* res, void* ctx)
  {
    NonLinearEigenProblem_Base<T>& var = *reinterpret_cast<NonLinearEigenProblem_Base<T>* >(ctx);
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
  Petsc_Error_Code NonLinearEigenProblem_Base<T>::
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
  void NonLinearEigenProblem_Base<T>::ComputeOperator(const T& L)
  {
    cout << "ComputeOperator not overloaded" << endl;
    abort();
  }
  
  //! to overload
  template<class T>
  void NonLinearEigenProblem_Base<T>::MltOperator(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "MltOperator not overloaded" << endl;
    abort();
  }
  
  //! to overload
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputeJacobian(const T& L)
  {
    cout << "ComputeJacobian not overloaded" << endl;
    abort();
  }
  
  //! to overload
  template<class T>
  void NonLinearEigenProblem_Base<T>::MltJacobian(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "MltJacobian not overloaded" << endl;
    abort();
  }
    
  //! to overload
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputePreconditioning()
  {
    cout << "ComputePreconditioning not overloaded" << endl;
    abort();
  }
  
  //! to overload
  template<class T>
  void NonLinearEigenProblem_Base<T>::ApplyPreconditioning(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "ApplyPreconditioning not overloaded" << endl;
    abort();
  }
  
}
 
#define SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_CXX
#endif
