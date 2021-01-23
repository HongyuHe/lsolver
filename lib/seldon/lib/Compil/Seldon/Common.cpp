#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "share/MatrixFlag.cxx"
#include "share/Errors.cxx"
#include "share/Common.cxx"
#include "share/Allocator.cxx"
#endif

namespace Seldon
{
  
  SELDON_EXTERN template class NaNAlloc<char>;
  SELDON_EXTERN template class NaNAlloc<bool>;
  SELDON_EXTERN template class NaNAlloc<unsigned>;
  SELDON_EXTERN template class NaNAlloc<int>;
  SELDON_EXTERN template class NaNAlloc<long>;
  SELDON_EXTERN template class NaNAlloc<unsigned long>;
  SELDON_EXTERN template class NaNAlloc<unsigned long long>;
  SELDON_EXTERN template class NaNAlloc<long long>;
  SELDON_EXTERN template class NaNAlloc<float>;
  SELDON_EXTERN template class NaNAlloc<double>;
  SELDON_EXTERN template class NaNAlloc<long double>;
  SELDON_EXTERN template class NaNAlloc<complex<float> >;
  SELDON_EXTERN template class NaNAlloc<complex<double> >;
  SELDON_EXTERN template class NaNAlloc<complex<long double> >;
  
  SELDON_EXTERN template class MallocObject<Vector<Real_wp> >;
  SELDON_EXTERN template class MallocObject<Vector<Complex_wp> >;

  SELDON_EXTERN template class MallocObject<Vector<Vector<Real_wp>, VectFull, MallocObject<Vector<Real_wp> > > >;
  SELDON_EXTERN template class MallocObject<Vector<Vector<Complex_wp>, VectFull, MallocObject<Vector<Complex_wp> > > >;
  
}
