#ifndef SELDON_FILE_VIRTUAL_EIGENVALUE_SOLVER_INLINE_CXX

namespace Seldon
{

  template<class T0, class T1>
  inline void to_complex_eigen(const T0& x, T1& y) { y = T1(x); }
  
  template<class T0, class T1>
  inline void to_complex_eigen(const complex<T0>& x, T1& y) { y = realpart(x); }

  template<class T0, class T1>
  inline void to_complex_eigen(const T0& x, complex<T1>& y) { y = complex<T1>(T1(x), T1(0)); }

  template<class T0, class T1>
  inline void to_complex_eigen(const complex<T0>& x, complex<T1>& y) { y = x; }

  inline int SlepcParam::GetEigensolverType() const
  {
    return type_solver;
  }

  
  inline void SlepcParam::SetEigensolverType(int type)
  {
    type_solver = type;
  }


  inline void SlepcParam::SetBlockSize(int n) { block_size = n; }
  inline int SlepcParam::GetBlockSize() const { return block_size; }
  inline void SlepcParam::SetMaximumBlockSize(int n) { block_max_size = n; }
  inline int SlepcParam::GetMaximumBlockSize() const { return block_max_size; }

  inline void SlepcParam::SetNumberOfSteps(int n) { nstep = n; }
  inline int SlepcParam::GetNumberOfSteps() const { return nstep; }

  inline void SlepcParam::SetExtractionType(int n) { type_extraction = n; }
  inline int SlepcParam::GetExtractionType() const { return type_extraction; }

  inline void SlepcParam::SetQuadratureRuleType(int n) { quadrature_rule = n; }
  inline int SlepcParam::GetQuadratureRuleType() const { return quadrature_rule; }

  inline void SlepcParam::SetInnerSteps(int n) { nstep_inner = n; }
  inline void SlepcParam::SetOuterSteps(int n) { nstep_outer = n; }
  inline int SlepcParam::GetInnerSteps() const { return nstep_inner; }
  inline int SlepcParam::GetOuterSteps() const { return nstep_outer; }

  inline int SlepcParam::GetNumberIntegrationPoints() const { return npoints; }
  inline void SlepcParam::SetNumberIntegrationPoints(int n) { npoints = n; }
  inline int SlepcParam::GetMomentSize() const { return moment_size; }
  inline void SlepcParam::SetMomentSize(int n) { moment_size = n; }
  inline int SlepcParam::GetNumberPartitions() const { return npart; }
  inline void SlepcParam::SetNumberPartitions(int n) { npart = n; }

  inline void SlepcParam::SetThresholdRank(double d) { delta_rank = d; }
  inline double SlepcParam::GetThresholdRank() const { return delta_rank; }
  inline void SlepcParam::SetThresholdSpurious(double d) { delta_spur = d; }
  inline double SlepcParam::GetThresholdSpurious() const { return delta_spur; }
  
  inline int SlepcParam::GetBorthogonalization() const { return borth; }
  inline void SlepcParam::SetBorthogonalization(int n) { borth = n; }
  inline int SlepcParam::GetDoubleExpansion() const { return double_exp; }
  inline void SlepcParam::SetDoubleExpansion(int n) { double_exp = n; } 
  inline int SlepcParam::GetInitialSize() const { return init_size; }
  inline void SlepcParam::SetInitialSize(int n) { init_size = n; }
  inline int SlepcParam::GetKrylovRestart() const { return krylov_restart; }
  inline void SlepcParam::SetKrylovRestart(int n) { krylov_restart = n; }
  inline int SlepcParam::GetRestartNumber() const { return restart_number; }
  inline void SlepcParam::SetRestartNumber(int n) { restart_number = n; }
  inline int SlepcParam::GetRestartNumberAdd() const { return restart_add; }
  inline void SlepcParam::SetRestartNumberAdd(int n) { restart_add = n; }
  
  inline int SlepcParam::GetNumberConvergedVectors() const { return nb_conv_vector; }
  inline void SlepcParam::SetNumberConvergedVectors(int n) { nb_conv_vector = n; }
  inline int SlepcParam::GetNumberConvergedVectorsProjected() const { return nb_conv_vector_proj; }
  inline void SlepcParam::SetNumberConvergedVectorsProjected(int n) { nb_conv_vector_proj = n; }
  
  inline bool SlepcParam::UseNonLockingVariant() const { return non_locking_variant; }
  inline void SlepcParam::SetNonLockingVariant(bool n) { non_locking_variant = n; }

  inline double SlepcParam::GetRestartRatio() const { return restart_ratio; }
  inline void SlepcParam::SetRestartRatio(double d) { restart_ratio = d; }
  inline string SlepcParam::GetMethod() const { return method; }
  inline void SlepcParam::SetMethod(string s) { method = s; }
  inline int SlepcParam::GetShiftType() const { return shift_type; }
  inline void SlepcParam::SetShiftType(int n) { shift_type = n; }

  inline void FeastParam::EnableEstimateNumberEigenval(bool p) { evaluate_number_eigenval = p; }
  inline bool FeastParam::EstimateNumberEigenval() const { return evaluate_number_eigenval; }
  
  inline int FeastParam::GetNumOfQuadraturePoints() const { return nb_points_quadrature; }
  inline void FeastParam::SetNumOfQuadraturePoints(int n) { nb_points_quadrature = n; }

  inline int FeastParam::GetTypeIntegration() const { return type_integration; }
  inline void FeastParam::SetTypeIntegration(int t) { type_integration = t; }
  
  template<class T>
  inline typename ClassComplexType<T>::Tcplx EigenProblem_Base<T>::GetCircleCenterSpectrum() const
  {
    return center_spectrum;
  }

  template<class T>
  inline typename ClassComplexType<T>::Treal EigenProblem_Base<T>::GetCircleRadiusSpectrum() const
  {
    return radius_spectrum;
  }

  template<class T>
  inline typename ClassComplexType<T>::Treal EigenProblem_Base<T>::GetRatioEllipseSpectrum() const
  {
    return ratio_ellipse;
  }

  template<class T>
  inline typename ClassComplexType<T>::Treal EigenProblem_Base<T>::GetAngleEllipseSpectrum() const
  {
    return angle_ellipse;
  }
 
}

#define SELDON_FILE_VIRTUAL_EIGENVALUE_SOLVER_INLINE_CXX
#endif

