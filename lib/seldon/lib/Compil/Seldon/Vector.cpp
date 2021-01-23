#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "vector/Vector.cxx"
#include "vector/SparseVector.cxx"
#include "vector/Functions_Arrays.cxx"
#include "computation/basic_functions/Functions_Vector.cxx"
#endif

namespace Seldon
{
  
  // dense vectors

  SELDON_EXTERN template class Vector_Base<bool>;
  SELDON_EXTERN template class Vector_Base<char>;
  SELDON_EXTERN template class Vector_Base<int>;
  SELDON_EXTERN template class Vector_Base<long>;
  SELDON_EXTERN template class Vector_Base<long long>;
  SELDON_EXTERN template class Vector_Base<double>;
  SELDON_EXTERN template class Vector_Base<complex<double> >;
  SELDON_EXTERN template class Vector_Base<float>;
  SELDON_EXTERN template class Vector_Base<complex<float> >;

  SELDON_EXTERN template class Vector<bool>;
  SELDON_EXTERN template class Vector<char>;
  SELDON_EXTERN template class Vector<int>;
  SELDON_EXTERN template class Vector<long>;
  SELDON_EXTERN template class Vector<long long>;
  SELDON_EXTERN template class Vector<double>;
  SELDON_EXTERN template class Vector<complex<double> >;
  SELDON_EXTERN template class Vector<float>;
  SELDON_EXTERN template class Vector<complex<float> >;

#ifdef SELDON_WITH_COMPILATION32
  SELDON_EXTERN template void Vector<unsigned long>::Resize(size_t);
#endif

#ifdef SELDON_WITH_MULTIPLE
  SELDON_EXTERN template class Vector<Real_wp>;
  SELDON_EXTERN template class Vector<Complex_wp>;

  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<double>&);
#endif

  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<bool>&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<char>&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<float>&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<complex<float> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<int>&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<int64_t>&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Real_wp>&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Complex_wp>&);

  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Vector<int> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Vector<Real_wp> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Vector<Complex_wp> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Vector<complex<float> > >&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Vector<Matrix<Complex_wp, General, ColMajor> >&);
  
  // sparse vectors
  SELDON_EXTERN template class Vector<int, VectSparse>;
  SELDON_EXTERN template class Vector<Real_wp, VectSparse>;
  SELDON_EXTERN template class Vector<Complex_wp, VectSparse>;

  SELDON_EXTERN template void Vector<int, VectSparse>::SetData(const Vector<int, VectSparse>&);
  SELDON_EXTERN template void Vector<int, VectSparse>::SetData(Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void Vector<Real_wp, VectSparse>::SetData(const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void Vector<Complex_wp, VectSparse>::SetData(const Vector<Complex_wp, VectSparse>&);

  SELDON_EXTERN template Vector<Real_wp, VectSparse>& Vector<Real_wp, VectSparse>::operator=(const Real_wp&);
  
  SELDON_EXTERN template void Vector<Real_wp, VectSparse>::RemoveSmallEntry(const Real_wp&);
  SELDON_EXTERN template void Vector<Complex_wp, VectSparse>::RemoveSmallEntry(const Real_wp&);
  
  SELDON_EXTERN template void Vector<Real_wp, VectSparse>::Copy(const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void Vector<Complex_wp, VectSparse>::Copy(const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template void Vector<Complex_wp, VectSparse>::Copy(const Vector<Complex_wp, VectSparse>&);
  
  SELDON_EXTERN template void Vector<Real_wp, VectSparse>::AddInteractionRow(size_t, const Vector<int>&, const Vector<Real_wp>&, bool);
  SELDON_EXTERN template void Vector<Complex_wp, VectSparse>::AddInteractionRow(size_t, const Vector<int>&, const Vector<Complex_wp>&, bool);

  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<int, VectSparse>&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Real_wp, VectSparse>&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Complex_wp, VectSparse>&);
  
  // vector of  strings  
  SELDON_EXTERN template class Vector_Base<string>;
  SELDON_EXTERN template void Vector<string>::Resize(size_t);
  SELDON_EXTERN template void Vector<string>::ReadText(string);
  SELDON_EXTERN template void Vector<string>::WriteText(string) const;

  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<string>&);

  // vector of vectors
  SELDON_EXTERN template class Vector_Base<Vector<bool> >;
  SELDON_EXTERN template void Vector<Vector<bool> >::Resize(size_t);
  
  SELDON_EXTERN template class Vector_Base<Vector<int> >;
  SELDON_EXTERN template void Vector<Vector<int> >::Resize(size_t);

  SELDON_EXTERN template class Vector_Base<Vector<Real_wp> >;
  SELDON_EXTERN template class Vector_Base<Vector<Real_wp>, MallocObject<Vector<Real_wp> > >;
  SELDON_EXTERN template void Vector<Vector<Real_wp> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Real_wp>, VectFull, MallocObject<Vector<Real_wp> > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<Real_wp> > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<Real_wp>, VectFull, MallocObject<Vector<Real_wp> > >, VectFull, MallocObject<Vector<Vector<Real_wp>, VectFull, MallocObject<Vector<Real_wp> > > > >::Resize(size_t);
  
  SELDON_EXTERN template class Vector_Base<Vector<Complex_wp> >;
  SELDON_EXTERN template class Vector_Base<Vector<Complex_wp>, MallocObject<Vector<Real_wp> > >;
  SELDON_EXTERN template void Vector<Vector<Complex_wp> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Complex_wp>, VectFull, MallocObject<Vector<Complex_wp> > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<Complex_wp> > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<Complex_wp>, VectFull, MallocObject<Vector<Complex_wp> > >, VectFull, MallocObject<Vector<Vector<Complex_wp>, VectFull, MallocObject<Vector<Complex_wp> > > > >::Resize(size_t);
  
  // vector of sparse vectors  
  SELDON_EXTERN template void Vector<Vector<Real_wp, VectSparse>, VectSparse>::ReallocateVector(size_t);
  SELDON_EXTERN template void Vector<Vector<Complex_wp, VectSparse>, VectSparse>::ReallocateVector(size_t);
  SELDON_EXTERN template void Vector<Vector<Real_wp, VectSparse>, VectSparse>::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Complex_wp, VectSparse>, VectSparse>::Resize(size_t);

  SELDON_EXTERN template void Vector<Vector<Real_wp, VectSparse> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Complex_wp, VectSparse> >::Resize(size_t);

  SELDON_EXTERN template void Vector<Vector<Real_wp, VectSparse>, VectSparse>::Clear();
  SELDON_EXTERN template void Vector<Vector<Complex_wp, VectSparse>, VectSparse>::Clear();

  SELDON_EXTERN template void Vector<Real_wp, VectSparse>::SetData(Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void Vector<Complex_wp, VectSparse>::SetData(Vector<Complex_wp>&, Vector<int>&);
  
  SELDON_EXTERN template void Vector<Vector<string> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<string> > >::Resize(size_t);

  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<Vector<string> >&);

  // Functions Sort, RemoveDuplicate, QuickSort, etc
  SELDON_EXTERN template void QuickSort(long, long, Vector<Real_wp>&);  
  SELDON_EXTERN template void QuickSort(long, long, Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void QuickSort(long, long, Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void QuickSort(long, long, Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void QuickSort(long, long, Vector<Real_wp>&, Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void MergeSort(long, long, Vector<int>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<Real_wp>&);  
  SELDON_EXTERN template void MergeSort(long, long, Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<Real_wp>&, Vector<int>&, Vector<int>&);

  SELDON_EXTERN template void MergeSort(long, long, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<int>&, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<int>&, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MergeSort(long, long, Vector<int>&, Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void Sort(long, long, Vector<int>&);
  SELDON_EXTERN template void Sort(long, long, Vector<Real_wp>&);  
  SELDON_EXTERN template void Sort(long, long, Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(long, long, Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, long, Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(long, long, Vector<Real_wp>&, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<Complex_wp>&, Vector<int>&);

  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<long>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<int>&, Vector<long>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(long, long, Vector<int>&, Vector<int>&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void Sort(long, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<Real_wp>&);  
  SELDON_EXTERN template void Sort(long, Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(long, Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&);

  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<long>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<int>&, Vector<long>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<Complex_wp>&, Vector<int>&);
  
  SELDON_EXTERN template void Sort(Vector<int>&);
  SELDON_EXTERN template void Sort(Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void Sort(Vector<int>&, Vector<long>&);
  SELDON_EXTERN template void Sort(Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void Sort(Vector<int>&, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void Sort(Vector<int>&, Vector<int>&, Vector<long>&);
  SELDON_EXTERN template void Sort(Vector<int>&, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(Vector<int>&, Vector<int>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void Sort(Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(Vector<Real_wp>&, Vector<int>&);
  SELDON_EXTERN template void Sort(Vector<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Sort(Vector<Real_wp>&, Vector<int>&, Vector<int>&);
    
  SELDON_EXTERN template void Sort(Vector<Real_wp>&, Vector<Vector<Real_wp> >&);

#ifdef SELDON_WITH_MULTIPLE
  SELDON_EXTERN template void Sort(Vector<double>&, Vector<Vector<double> >&);
#endif

  SELDON_EXTERN template void Sort(Vector<Real_wp>&, Vector<Complex_wp>&, Vector<Vector<Real_wp> >&);
  SELDON_EXTERN template void Sort(Vector<Real_wp>&, Vector<Vector<Complex_wp> >&);
  SELDON_EXTERN template void Sort(Vector<Real_wp>&, Vector<Complex_wp>&, Vector<Vector<Complex_wp> >&);
  
  SELDON_EXTERN template void Assemble(long&, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Assemble(long&, Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Assemble(long&, Vector<int>&);
  SELDON_EXTERN template void Assemble(Vector<int>&);

  SELDON_EXTERN template void RemoveDuplicate(long&, Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void RemoveDuplicate(long&, Vector<int>&, Vector<Real_wp>&);
  SELDON_EXTERN template void RemoveDuplicate(long&, Vector<int>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void RemoveDuplicate(long&, Vector<int>&);
  SELDON_EXTERN template void RemoveDuplicate(Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void RemoveDuplicate(Vector<int>&);
  SELDON_EXTERN template void RemoveDuplicate(Vector<Real_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void RemoveDuplicate(Vector<Real_wp>&);
  
  SELDON_EXTERN template bool HasElement(const Vector<int>&, const int &);
  SELDON_EXTERN template bool HasElement(const Vector<Real_wp>&, const Real_wp &);
  SELDON_EXTERN template bool HasElement(const Vector<Complex_wp>&, const Complex_wp &);
  
  SELDON_EXTERN template void SwapPointer(Vector<int>&, Vector<int>&);
  SELDON_EXTERN template void SwapPointer(Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SwapPointer(Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void SwapPointer(Vector<Vector<int> >&, Vector<Vector<int> >&);
  SELDON_EXTERN template void SwapPointer(Vector<Vector<Real_wp, VectSparse> >&, Vector<Vector<Real_wp, VectSparse> >&);
  SELDON_EXTERN template void SwapPointer(Vector<Vector<Complex_wp, VectSparse> >&, Vector<Vector<Complex_wp, VectSparse> >&);

}
