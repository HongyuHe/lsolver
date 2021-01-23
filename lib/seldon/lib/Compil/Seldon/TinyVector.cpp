#include "SeldonFlag.hxx"

#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "vector/Vector.cxx"
#include "vector/TinyVector.cxx"
#include "matrix/TinyMatrix.cxx"
#endif

namespace Seldon
{
  
  /* TinyVector */
  
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<TinyVector<Real_wp, 1> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<TinyVector<Real_wp, 2> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<TinyVector<Real_wp, 3> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<TinyVector<Complex_wp, 1> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<TinyVector<Complex_wp, 2> >&);
  SELDON_EXTERN template ostream& operator <<(ostream&, const Vector<TinyVector<Complex_wp, 3> >&);

  SELDON_EXTERN template class Vector_Base<TinyVector<Real_wp, 1> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Real_wp, 2> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Real_wp, 3> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Real_wp, 4> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Real_wp, 5> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Real_wp, 6> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Real_wp, 7> >;

  SELDON_EXTERN template class Vector_Base<TinyVector<Complex_wp, 1> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Complex_wp, 2> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Complex_wp, 3> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Complex_wp, 4> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Complex_wp, 5> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Complex_wp, 6> >;
  SELDON_EXTERN template class Vector_Base<TinyVector<Complex_wp, 7> >;
  
  SELDON_EXTERN template Real_wp Norm2_Column(const Matrix<Real_wp>&, int, int);

  /* TinyMatrix */

  SELDON_EXTERN template void GetCholesky(TinyMatrix<Real_wp, Symmetric, 1, 1>&);
  SELDON_EXTERN template void SolveCholesky(const class_SeldonTrans&, const TinyMatrix<Real_wp, Symmetric, 1, 1>&, TinyVector<Real_wp, 1>&);
  SELDON_EXTERN template void SolveCholesky(const class_SeldonNoTrans&, const TinyMatrix<Real_wp, Symmetric, 1, 1>&, TinyVector<Real_wp, 1>&);

  // 2x2 matrices

  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Real_wp, General, 2, 2>&, int, int);
  
  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, General, 2, 2>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, General, 2, 2>&, TinyMatrix<Real_wp, General, 2, 2>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, General, 2, 2>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, General, 2, 2>&, TinyMatrix<Complex_wp, General, 2, 2>&);

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, Symmetric, 2, 2>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, Symmetric, 2, 2>&, TinyMatrix<Real_wp, Symmetric, 2, 2>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, Symmetric, 2, 2>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, Symmetric, 2, 2>&, TinyMatrix<Complex_wp, Symmetric, 2, 2>&);

  SELDON_EXTERN template void GetEigenvalues(TinyMatrix<Real_wp, General, 2, 2>&, TinyVector<Real_wp, 2>&, TinyVector<Real_wp, 2>&);
  SELDON_EXTERN template void GetEigenvalues(TinyMatrix<Complex_wp, General, 2, 2>&, TinyVector<Complex_wp, 2>&);

  SELDON_EXTERN template void GetEigenvalues(TinyMatrix<Real_wp, Symmetric, 2, 2>&, TinyVector<Real_wp, 2>&);
  SELDON_EXTERN template void GetEigenvalues(TinyMatrix<Complex_wp, Symmetric, 2, 2>&, TinyVector<Complex_wp, 2>&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(TinyMatrix<Real_wp, Symmetric, 2, 2>&, TinyVector<Real_wp, 2>&, TinyMatrix<Real_wp, General, 2, 2>&);

  SELDON_EXTERN template void GetSquareRoot(TinyMatrix<Real_wp, Symmetric, 2, 2>&);

  SELDON_EXTERN template void GetEigenvaluesEigenvectors(TinyMatrix<Complex_wp, Symmetric, 2, 2>&, TinyVector<Complex_wp, 2>&, TinyMatrix<Complex_wp, General, 2, 2>&);

  SELDON_EXTERN template void GetSquareRoot(TinyMatrix<Complex_wp, Symmetric, 2, 2>&);

  SELDON_EXTERN template void GetCholesky(TinyMatrix<Real_wp, Symmetric, 2, 2>&);
  SELDON_EXTERN template void SolveCholesky(const class_SeldonTrans&, const TinyMatrix<Real_wp, Symmetric, 2, 2>&, TinyVector<Real_wp, 2>&);
  SELDON_EXTERN template void SolveCholesky(const class_SeldonNoTrans&, const TinyMatrix<Real_wp, Symmetric, 2, 2>&, TinyVector<Real_wp, 2>&);
  SELDON_EXTERN template void MltCholesky(const class_SeldonTrans&, const TinyMatrix<Real_wp, Symmetric, 2, 2>&, TinyVector<Real_wp, 2>&);
  SELDON_EXTERN template void MltCholesky(const class_SeldonNoTrans&, const TinyMatrix<Real_wp, Symmetric, 2, 2>&, TinyVector<Real_wp, 2>&);

  SELDON_EXTERN template int IntersectionEdges(const TinyVector<Real_wp, 2>& pointA, const TinyVector<Real_wp, 2>& pointB,
                                               const TinyVector<Real_wp, 2>& pt1, const TinyVector<Real_wp, 2>& pt2,
                                               TinyVector<Real_wp, 2>& res, const Real_wp& threshold);

  SELDON_EXTERN template int IntersectionDroites(const TinyVector<Real_wp, 2>& pointA, const TinyVector<Real_wp, 2>& pointB,
                                                 const TinyVector<Real_wp, 2>& pt1, const TinyVector<Real_wp, 2>& pt2,
                                                 TinyVector<Real_wp, 2>& res, const Real_wp& threshold);
  
  // 3x3 matrices

  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Real_wp, General, 3, 3>&, int, int);

  SELDON_EXTERN template Real_wp Det(const TinyMatrix<Real_wp, General, 3, 3>&);
  SELDON_EXTERN template Real_wp Det(const TinyMatrix<Real_wp, Symmetric, 3, 3>&);

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, General, 3, 3>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, General, 3, 3>&, TinyMatrix<Real_wp, General, 3, 3>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, General, 3, 3>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, General, 3, 3>&, TinyMatrix<Complex_wp, General, 3, 3>&);

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, Symmetric, 3, 3>&, TinyMatrix<Real_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, Symmetric, 3, 3>&, TinyMatrix<Complex_wp, Symmetric, 3, 3>&);

  SELDON_EXTERN template void GetNormalProjector(const TinyVector<Real_wp, 3>&, TinyMatrix<Real_wp, General, 3, 3>&);
  SELDON_EXTERN template void GetNormalProjector(const TinyVector<Real_wp, 3>&, TinyMatrix<Real_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void GetNormalProjector(const TinyVector<Real_wp, 3>&, TinyMatrix<Complex_wp, General, 3, 3>&);
  SELDON_EXTERN template void GetNormalProjector(const TinyVector<Real_wp, 3>&, TinyMatrix<Complex_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void GetNormalProjector(const TinyVector<Complex_wp, 3>&, TinyMatrix<Complex_wp, General, 3, 3>&);
  SELDON_EXTERN template void GetNormalProjector(const TinyVector<Complex_wp, 3>&, TinyMatrix<Complex_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void GetTangentialProjector(const TinyVector<Real_wp, 3>&, TinyMatrix<Real_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void GetTangentialProjector(const TinyVector<Real_wp, 3>&, TinyMatrix<Complex_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void GetVectorPlane(const TinyVector<Real_wp, 3>&, TinyVector<Real_wp, 3>&, TinyVector<Real_wp, 3>&);
  
  SELDON_EXTERN template void GetEigenvalues(TinyMatrix<Real_wp, Symmetric, 3, 3>&, TinyVector<Real_wp, 3>&);
  SELDON_EXTERN template void GetEigenvalues(TinyMatrix<Complex_wp, Symmetric, 3, 3>&, TinyVector<Complex_wp, 3>&);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(TinyMatrix<Real_wp, Symmetric, 3, 3>&, TinyVector<Real_wp, 3>&, TinyMatrix<Real_wp, General, 3, 3>&);

  SELDON_EXTERN template void GetSquareRoot(TinyMatrix<Real_wp, Symmetric, 3, 3>&);

  SELDON_EXTERN template void GetEigenvaluesEigenvectors(TinyMatrix<Complex_wp, Symmetric, 3, 3>&, TinyVector<Complex_wp, 3>&, TinyMatrix<Complex_wp, General, 3, 3>&);

  SELDON_EXTERN template void GetSquareRoot(TinyMatrix<Complex_wp, Symmetric, 3, 3>&);

  SELDON_EXTERN template void GetCholesky(TinyMatrix<Real_wp, Symmetric, 3, 3>&);
  SELDON_EXTERN template void SolveCholesky(const class_SeldonTrans&, const TinyMatrix<Real_wp, Symmetric, 3, 3>&, TinyVector<Real_wp, 3>&);
  SELDON_EXTERN template void SolveCholesky(const class_SeldonNoTrans&, const TinyMatrix<Real_wp, Symmetric, 3, 3>&, TinyVector<Real_wp, 3>&);
  SELDON_EXTERN template void MltCholesky(const class_SeldonTrans&, const TinyMatrix<Real_wp, Symmetric, 3, 3>&, TinyVector<Real_wp, 3>&);
  SELDON_EXTERN template void MltCholesky(const class_SeldonNoTrans&, const TinyMatrix<Real_wp, Symmetric, 3, 3>&, TinyVector<Real_wp, 3>&);

  // 4x4 matrices

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, General, 4, 4>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, General, 4, 4>&, TinyMatrix<Real_wp, General, 4, 4>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, General, 4, 4>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, General, 4, 4>&, TinyMatrix<Complex_wp, General, 4, 4>&);

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, Symmetric, 4, 4>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, Symmetric, 4, 4>&, TinyMatrix<Real_wp, Symmetric, 4, 4>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, Symmetric, 4, 4>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, Symmetric, 4, 4>&, TinyMatrix<Complex_wp, Symmetric, 4, 4>&);

  // 5x5 matrices

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, General, 5, 5>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, General, 5, 5>&, TinyMatrix<Real_wp, General, 5, 5>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, General, 5, 5>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, General, 5, 5>&, TinyMatrix<Complex_wp, General, 5, 5>&);

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, Symmetric, 5, 5>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, Symmetric, 5, 5>&, TinyMatrix<Real_wp, Symmetric, 5, 5>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, Symmetric, 5, 5>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, Symmetric, 5, 5>&, TinyMatrix<Complex_wp, Symmetric, 5, 5>&);

  SELDON_EXTERN template void GetEigenvaluesEigenvectors(TinyMatrix<Real_wp, Symmetric, 5, 5>&, TinyVector<Real_wp, 5>&, TinyMatrix<Real_wp, General, 5, 5>&);

  SELDON_EXTERN template void GetSquareRoot(TinyMatrix<Real_wp, Symmetric, 5, 5>&);

  SELDON_EXTERN template void GetEigenvaluesEigenvectors(TinyMatrix<Complex_wp, Symmetric, 5, 5>&, TinyVector<Complex_wp, 5>&, TinyMatrix<Complex_wp, General, 5, 5>&);

  SELDON_EXTERN template void GetSquareRoot(TinyMatrix<Complex_wp, Symmetric, 5, 5>&);

  // 6x6 matrices

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, General, 6, 6>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, General, 6, 6>&, TinyMatrix<Real_wp, General, 6, 6>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, General, 6, 6>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, General, 6, 6>&, TinyMatrix<Complex_wp, General, 6, 6>&);

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, Symmetric, 6, 6>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, Symmetric, 6, 6>&, TinyMatrix<Real_wp, Symmetric, 6, 6>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, Symmetric, 6, 6>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, Symmetric, 6, 6>&, TinyMatrix<Complex_wp, Symmetric, 6, 6>&);
  
  // 7x7 matrices

  SELDON_EXTERN template void GetInverse(TinyMatrix<Real_wp, General, 7, 7>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Real_wp, General, 7, 7>&, TinyMatrix<Real_wp, General, 7, 7>&);
  SELDON_EXTERN template void GetInverse(TinyMatrix<Complex_wp, General, 7, 7>&);
  SELDON_EXTERN template void GetInverse(const TinyMatrix<Complex_wp, General, 7, 7>&, TinyMatrix<Complex_wp, General, 7, 7>&);

  // Norm2_Column
  
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Real_wp, General, 5, 4>&, int, int);
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Real_wp, General, 7, 7>&, int, int);
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Complex_wp, General, 5, 4>&, int, int);
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Complex_wp, General, 7, 7>&, int, int);
  
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Real_wp, Symmetric, 2, 2>&, int, int);
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Real_wp, Symmetric, 3, 3>&, int, int);
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Real_wp, Symmetric, 5, 5>&, int, int);
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Complex_wp, Symmetric, 2, 2>&, int, int);
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Complex_wp, Symmetric, 3, 3>&, int, int);
  SELDON_EXTERN template Real_wp Norm2_Column(const TinyMatrix<Complex_wp, Symmetric, 5, 5>&, int, int);
  
}
