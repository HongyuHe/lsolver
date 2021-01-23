#ifndef SELDON_FILE_SLEPC_HXX

// we include only basic defines of Petsc to know the type of PetscScalar
// with this trick, we can recompile easily Seldon when Slepc is enabled/disabled
#include <petscconf.h>

#if defined(PETSC_USE_COMPLEX)
#define Petsc_Scalar complex<double>
#else
#define Petsc_Scalar double
#endif

extern "C"
{
  int SlepcInitialize(int*, char***, const char[], const char[]);
  int SlepcFinalize();  
}

namespace Seldon
{

  void FindEigenvaluesSlepc_(EigenProblem_Base<Petsc_Scalar>& var,
			     Vector<Petsc_Scalar>& eigen_values,
			     Vector<Petsc_Scalar>& lambda_imag,
			     Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors);

  void FindEigenvaluesSlepc_(PolynomialEigenProblem_Base<Petsc_Scalar>& var,
                                Vector<Petsc_Scalar>& eigen_values,
                                Vector<Petsc_Scalar>& lambda_imag,
                                Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors);
  
#ifdef SELDON_WITH_VIRTUAL
  template<class T>
  inline void FindEigenvaluesSlepc(EigenProblem_Base<T>& var,
				   Vector<T>& eigen_values,
				   Vector<T>& lambda_imag,
				   Matrix<T, General, ColMajor>& eigen_vectors)
  {
    cout << "Link with the correct Slepc/Petsc library (real/complex) " << endl;
    cout << "sizeof(T) = " << sizeof(T) << endl;
    cout << "sizeof(Petsc_Scalar) = " << sizeof(Petsc_Scalar) << endl;
  }
  
  void FindEigenvaluesSlepc(EigenProblem_Base<Petsc_Scalar>& var,
			    Vector<Petsc_Scalar>& eigen_values,
			    Vector<Petsc_Scalar>& lambda_imag,
			    Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors);  

#ifdef SELDON_WITH_SLEPC_PEP
  template<class T>
  inline void FindEigenvaluesSlepc(PolynomialEigenProblem_Base<T>& var,
                                   Vector<T>& eigen_values,
                                   Vector<T>& lambda_imag,
                                   Matrix<T, General, ColMajor>& eigen_vectors)
  {
    cout << "Link with the correct Slepc/Petsc library (real/complex) " << endl;
    cout << "sizeof(T) = " << sizeof(T) << endl;
    cout << "sizeof(Petsc_Scalar) = " << sizeof(Petsc_Scalar) << endl;
  }
  
  void FindEigenvaluesSlepc(PolynomialEigenProblem_Base<Petsc_Scalar>& var,
                            Vector<Petsc_Scalar>& eigen_values,
                            Vector<Petsc_Scalar>& lambda_imag,
                            Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors);  
#endif
  
  template<class T>
  inline void FindEigenvaluesSlepc(NonLinearEigenProblem_Base<T>& var,
                                   Vector<T>& eigen_values,
                                   Vector<T>& lambda_imag,
                                   Matrix<T, General, ColMajor>& eigen_vectors)
  {
    cout << "Link with the correct Slepc/Petsc library (real/complex) " << endl;
    cout << "sizeof(T) = " << sizeof(T) << endl;
    cout << "sizeof(Petsc_Scalar) = " << sizeof(Petsc_Scalar) << endl;
  }
  
  void FindEigenvaluesSlepc(NonLinearEigenProblem_Base<Petsc_Scalar>& var,
                            Vector<Petsc_Scalar>& eigen_values,
                            Vector<Petsc_Scalar>& lambda_imag,
                            Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors);  
#else
  template<class EigenProblem, class T, class Allocator1,
           class Allocator2, class Allocator3>
  void FindEigenvaluesSlepc(EigenProblem& var,
			    Vector<T, VectFull, Allocator1>& eigen_values,
			    Vector<T, VectFull, Allocator2>& lambda_imag,
			    Matrix<T, General, ColMajor, Allocator3>& eigen_vectors);
#endif

}

#define SELDON_FILE_SLEPC_HXX
#endif
