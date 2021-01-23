#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#include "vector/Vector2.hxx"
#include "vector/Vector3.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "vector/Vector2.cxx"
#include "vector/Vector3.cxx"
#endif

namespace Seldon
{
  
  // dense vectors

  SELDON_EXTERN template class Vector2<Real_wp>;
  SELDON_EXTERN template class Vector2<Complex_wp>;

  SELDON_EXTERN template void Vector2<Real_wp>::Flatten(Vector<Real_wp>&) const;
  SELDON_EXTERN template void Vector2<Real_wp>::Flatten(long, long, Vector<Real_wp>&) const;

  SELDON_EXTERN template void Vector2<Complex_wp>::Flatten(Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Vector2<Complex_wp>::Flatten(long, long, Vector<Complex_wp>&) const;

  SELDON_EXTERN template class Vector3<Real_wp>;
  SELDON_EXTERN template class Vector3<Complex_wp>;

  SELDON_EXTERN template void Vector3<Real_wp>::Flatten(Vector<Real_wp>&) const;
  SELDON_EXTERN template void Vector3<Real_wp>::Flatten(int, int, Vector<Real_wp>&) const;
  SELDON_EXTERN template void Vector3<Real_wp>::Flatten(int, int, int, int, Vector<Real_wp>&) const;

  SELDON_EXTERN template void Vector3<Complex_wp>::Flatten(Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Vector3<Complex_wp>::Flatten(int, int, Vector<Complex_wp>&) const;
  SELDON_EXTERN template void Vector3<Complex_wp>::Flatten(int, int, int, int, Vector<Complex_wp>&) const;
  
}
