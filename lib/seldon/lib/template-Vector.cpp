// Copyright (C) 2010 Vivien Mallet
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#include "SeldonHeader.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "share/Errors.hxx"
#include "share/Allocator.cxx"
#include "vector/VectorInline.cxx"
#include "vector/Vector.cxx"
#include "vector/Functions_Arrays.cxx"
#include "share/Common.cxx"
#endif




#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;


namespace Seldon
{

  SELDON_EXTERN template class MallocAlloc<@scalar>;
  SELDON_EXTERN template class Vector_Base<@scalar, MallocAlloc<@scalar> >;
  SELDON_EXTERN template class Vector<@scalar, VectFull, MallocAlloc<@scalar> >;

  // Function templates.
  SELDON_EXTERN template void Vector<@scalar, VectFull, MallocAlloc<@scalar> >::SetData(const Vector<@scalar, VectFull, MallocAlloc<@scalar> >&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, MallocAlloc<@scalar> >::PushBack(const @scalar&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, MallocAlloc<@scalar> >::PushBack(const Vector<@scalar, VectFull, MallocAlloc<@scalar> >&);
  SELDON_EXTERN template void Vector<@real_complex, VectFull, MallocAlloc<@real_complex> >::Fill(const @real_complex&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, MallocAlloc<@scalar> >::Fill(const int&);
#ifndef SWIG
  SELDON_EXTERN template Vector<@scalar, VectFull, MallocAlloc<@scalar> >& Vector<@scalar, VectFull, MallocAlloc<@scalar> >::operator= (const @scalar&);
#endif
  SELDON_EXTERN template Vector<@scalar, VectFull, MallocAlloc<@scalar> >& Vector<@scalar, VectFull, MallocAlloc<@scalar> >::operator*= (const @scalar&);

#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream& out, const Vector<@scalar, VectFull, MallocAlloc<@scalar> >& V);
#endif

  SELDON_EXTERN template class NewAlloc<@scalar>;
  SELDON_EXTERN template class Vector_Base<@scalar, NewAlloc<@scalar> >;
  SELDON_EXTERN template class Vector<@scalar, VectFull, NewAlloc<@scalar> >;

  // Function templates.
  SELDON_EXTERN template void Vector<@scalar, VectFull, NewAlloc<@scalar> >::SetData(const Vector<@scalar, VectFull, NewAlloc<@scalar> >&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, NewAlloc<@scalar> >::PushBack(const @scalar&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, NewAlloc<@scalar> >::PushBack(const Vector<@scalar, VectFull, NewAlloc<@scalar> >&);
  SELDON_EXTERN template void Vector<@real_complex, VectFull, NewAlloc<@real_complex> >::Fill(const @real_complex&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, NewAlloc<@scalar> >::Fill(const int&);
#ifndef SWIG
  SELDON_EXTERN template Vector<@scalar, VectFull, NewAlloc<@scalar> >& Vector<@scalar, VectFull, NewAlloc<@scalar> >::operator= (const @scalar&);
#endif
  SELDON_EXTERN template Vector<@scalar, VectFull, NewAlloc<@scalar> >& Vector<@scalar, VectFull, NewAlloc<@scalar> >::operator*= (const @scalar&);

#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream& out, const Vector<@scalar, VectFull, NewAlloc<@scalar> >& V);
#endif

  // Functions Sort, RemoveDuplicate, QuickSort, etc
  SELDON_EXTERN template void QuickSort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void QuickSort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void QuickSort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void MergeSort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void MergeSort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void MergeSort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(int, int, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(int, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(int, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(int, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Sort(Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Assemble(int&, Vector<int, VectFull, NewAlloc<int> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void Assemble(int&, Vector<int, VectFull, NewAlloc<int> >&);

  SELDON_EXTERN template void Assemble(Vector<int, VectFull, NewAlloc<int> >&);

  SELDON_EXTERN template void RemoveDuplicate(int&, Vector<int, VectFull, NewAlloc<int> >&, Vector<@real, VectFull, NewAlloc<@real> >&);

  SELDON_EXTERN template void RemoveDuplicate(int&, Vector<int, VectFull, NewAlloc<int> >&);

  SELDON_EXTERN template void RemoveDuplicate(Vector<int, VectFull, NewAlloc<int> >&);
}
