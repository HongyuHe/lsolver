#ifndef SELDON_FILE_ANASAZI_HXX

namespace Seldon
{
  
  template<class T>
  void SetComplexEigenvalue(const T& x, const T& y, T& z, T& zimag);
  
  template<class T>
  void SetComplexEigenvalue(const T& x, const T& y, complex<T>& z, complex<T>& zimag);

#ifdef SELDON_WITH_VIRTUAL
  template<class T>
  void FindEigenvaluesAnasazi(EigenProblem_Base<T>& var,
			      Vector<T>& eigen_values,
			      Vector<T>& lambda_imag,
			      Matrix<T, General, ColMajor>& eigen_vectors);
#else
  template<class EigenProblem, class T, class Allocator1,
           class Allocator2, class Allocator3>
  void FindEigenvaluesAnasazi(EigenProblem& var,
			      Vector<T, VectFull, Allocator1>& eigen_values,
			      Vector<T, VectFull, Allocator2>& lambda_imag,
			      Matrix<T, General, ColMajor, Allocator3>& eigen_vectors);
#endif
  
}

#define SELDON_FILE_ANASAZI_HXX
#endif
