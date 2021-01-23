#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "array/Array3D.cxx"
#include "array/Array4D.cxx"
#endif

namespace Seldon
{
  
  SELDON_EXTERN template class Array3D<int>;
  SELDON_EXTERN template class Array3D<float>;
  SELDON_EXTERN template class Array3D<double>;
  SELDON_EXTERN template class Array3D<complex<float> >;
  SELDON_EXTERN template class Array3D<complex<double> >;

#ifdef SELDON_WITH_MULTIPLE
  SELDON_EXTERN template class Array3D<Real_wp>;
  SELDON_EXTERN template class Array3D<Complex_wp>;
#endif

  SELDON_EXTERN template void Array3D<int>::Fill(const int&);

  SELDON_EXTERN template void Array3D<Real_wp>::Fill(const int&);
  SELDON_EXTERN template void Array3D<Real_wp>::Fill(const Real_wp&);

  SELDON_EXTERN template void Array3D<Complex_wp>::Fill(const int&);
  SELDON_EXTERN template void Array3D<Complex_wp>::Fill(const Real_wp&);
  SELDON_EXTERN template void Array3D<Complex_wp>::Fill(const Complex_wp&);

  SELDON_EXTERN template void MltScalar(const Real_wp&, Array3D<Real_wp>&);

  SELDON_EXTERN template void MltScalar(const Real_wp&, Array3D<Complex_wp>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Array3D<Complex_wp>&);

  SELDON_EXTERN template class Array4D<Real_wp>;
  SELDON_EXTERN template class Array4D<Complex_wp>;
  
  SELDON_EXTERN template void Array4D<Real_wp>::Fill(const int&);
  SELDON_EXTERN template void Array4D<Real_wp>::Fill(const Real_wp&);

  SELDON_EXTERN template void Array4D<Complex_wp>::Fill(const int&);
  SELDON_EXTERN template void Array4D<Complex_wp>::Fill(const Real_wp&);
  SELDON_EXTERN template void Array4D<Complex_wp>::Fill(const Complex_wp&);

  SELDON_EXTERN template ostream& operator<<(ostream&, const Array3D<int>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Array3D<Real_wp>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const Array3D<Complex_wp>&);

}
