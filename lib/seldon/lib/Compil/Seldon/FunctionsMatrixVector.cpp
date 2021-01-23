#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/basic_functions/Functions_MatVect.cxx"
#include "computation/basic_functions/Functions_Base.cxx"
#ifdef SELDON_WITH_MKL
#include "computation/interfaces/Mkl_Sparse.cxx"
#endif
#endif

namespace Seldon
{
  /* Matrix-vector products with sparse matrices */
  
  /* RowSparse */

  // Mkl like functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);  
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp , General, RowSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp , General, RowSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp , General, RowSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp , General, RowSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  /* ColSparse */
  
  // Mkl like functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ColSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);  
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, ColSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ColSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ColSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp , General, ColSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ColSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp , General, ColSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ColSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
      
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ColSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, ColSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp , General, ColSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ColSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ColSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp , General, ColSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);


  /* RowSymSparse */
  
  // Mkl like functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);  
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp , Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
    
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp , Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp , Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp , Symmetric, RowSymSparse>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  

  /* Matrix-vector products with dense matrices */


  /* RowMajor */

  // functions that should be blas
  
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  
  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
 
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  /* ColMajor */

  // functions that should be blas
  
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  
  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
 
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  /* RowSymPacked */

  // function that should be Blas
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);

  // for complex Real_wp (Blas functions not present)
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);

  
  /* RowSym */

  // functions that should be Blas
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSym>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  
  // for complex Real_wp
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, Symmetric, RowSym>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);


  SELDON_EXTERN template void MltVector(Matrix<Complex_wp, Symmetric, RowSym> const&, Vector<Complex_wp> const&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, Matrix<Complex_wp, Symmetric, RowSym> const&, Vector<Complex_wp> const&, Vector<Complex_wp>&);
  
  /* RowHermPacked */
    
  // functions that should be Blas
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);

  /* RowHerm */

  // functions that should be Blas
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  /* RowLoTriang */
  
  // functions that should be blas
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);  

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
   
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriang>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  /* RowLoTriangPacked */

  // functions that should be Blas
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowLoTriangPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  /* RowUpTriang */

  // functions that should be Blas
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
    
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriang>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(Complex_wp const&, Matrix<Real_wp, General, RowUpTriang> const&, Vector<Complex_wp> const&, Complex_wp const&, Vector<Complex_wp>&);
  
  /* RowUpTriangPacked */
  
  // functions that should be Blas
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);

  // other functions
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, RowUpTriangPacked>&, const Vector<Complex_wp >&, const Real_wp&, Vector<Complex_wp >&);

  /* Functions with dense matrices */

  
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp >&);
  SELDON_EXTERN template void SolveLU(const Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLU(const Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void Gauss(Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Gauss(Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp >&);
  SELDON_EXTERN template void Gauss(Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Gauss(Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, const Vector<Real_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, const Vector<Real_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&, const Vector<Real_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&, const Vector<Real_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Symmetric, RowSym>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);

  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, const Vector<Real_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowMajor>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, const Vector<Real_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, General, ColMajor>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&, const Vector<Real_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&, const Vector<Real_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSym>&, Vector<Complex_wp>&, const Vector<Complex_wp>&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, RowSym>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void SorVector(const SeldonTranspose&, const Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, Real_wp const&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Real_wp, General, RowMajor>&, Vector<Real_wp>&, const Vector<Real_wp>&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Complex_wp, General, RowMajor>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&, const Vector<Real_wp>&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Real_wp>&, const Vector<Real_wp>&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Complex_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Real_wp, Symmetric, RowSym>&, Vector<Real_wp>&, const Vector<Real_wp>&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Complex_wp, Symmetric, RowSym>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Complex_wp, Hermitian, RowHermPacked>&, Vector<Complex_wp >&, const Vector<Complex_wp >&, int, int);
  SELDON_EXTERN template void GaussSeidel(const Matrix<Complex_wp, Hermitian, RowHerm>&, Vector<Complex_wp >&, const Vector<Complex_wp>&, int, int);
    
}
