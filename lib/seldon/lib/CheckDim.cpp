#include "SeldonHeader.hxx"
#include "matrix_sparse/Matrix_ArraySparse.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "matrix_sparse/Matrix_ArraySparse.cxx"
#include "computation/basic_functions/Functions_Vector.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#endif

namespace Seldon
{
  /* Vectors */
  
  SELDON_EXTERN template void 
  CheckDim(const Vector<double, VectFull>&,
	   const Vector<double, VectFull>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Vector<complex<double>, VectFull>&,
	   const Vector<complex<double>, VectFull>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Vector<double, VectSparse>&,
	   const Vector<double, VectSparse>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Vector<complex<double>, VectSparse>&,
	   const Vector<complex<double>, VectSparse>&, string, string);
  
  /* ColMajor */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ColMajor>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ColMajor>&,
	   const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ColMajor>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ColMajor>&,
           const Matrix<double, General, ColMajor>&,
           const Matrix<double, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<double, General, ColMajor>&,
           const SeldonTranspose&,
           const Matrix<double, General, ColMajor>&,
           const Matrix<double, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<double, General, ColMajor>&,
           const SeldonTranspose&,
           const Matrix<double, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ColMajor>&,
           const Matrix<double, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColMajor>&,
           const Matrix<complex<double>, General, ColMajor>&,
           const Matrix<complex<double>, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<complex<double>, General, ColMajor>&,
           const SeldonTranspose&,
           const Matrix<complex<double>, General, ColMajor>&,
           const Matrix<complex<double>, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<complex<double>, General, ColMajor>&,
           const SeldonTranspose&,
           const Matrix<complex<double>, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColMajor>&,
           const Matrix<complex<double>, General, ColMajor>&, string);
  
  /* RowMajor */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, RowMajor>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, RowMajor>&,
	   const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, RowMajor>&,
           const Matrix<double, General, RowMajor>&,
           const Matrix<double, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<double, General, RowMajor>&,
           const SeldonTranspose&,
           const Matrix<double, General, RowMajor>&,
           const Matrix<double, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<double, General, RowMajor>&,
           const SeldonTranspose&,
           const Matrix<double, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, RowMajor>&,
           const Matrix<double, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowMajor>&,
           const Matrix<complex<double>, General, RowMajor>&,
           const Matrix<complex<double>, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<complex<double>, General, RowMajor>&,
           const SeldonTranspose&,
           const Matrix<complex<double>, General, RowMajor>&,
           const Matrix<complex<double>, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<complex<double>, General, RowMajor>&,
           const SeldonTranspose&,
           const Matrix<complex<double>, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowMajor>&,
           const Matrix<complex<double>, General, RowMajor>&, string);
  
  /* RowSymPacked */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, RowSymPacked>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSymPacked>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSymPacked>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, RowSymPacked>&,
	   const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSymPacked>&,
	   const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSymPacked>&,
	   const Vector<double>&, string, string);

  /* ColSymPacked */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ColSymPacked>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSymPacked>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSymPacked>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ColSymPacked>&,
	   const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSymPacked>&,
	   const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSymPacked>&,
	   const Vector<double>&, string, string);
  
  /* RowSym */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, RowSym>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSym>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSym>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, RowSym>&,
	   const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSym>&,
	   const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSym>&,
	   const Vector<double>&, string, string);
  
  /* ColSym */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ColSym>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSym>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSym>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ColSym>&,
	   const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSym>&,
	   const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSym>&,
	   const Vector<double>&, string, string);

  /* RowHerm */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Hermitian, RowHerm>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Hermitian, RowHerm>&,
	   const Vector<complex<double> >&, string, string);

  /* ColHerm */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Hermitian, ColHerm>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Hermitian, ColHerm>&,
	   const Vector<complex<double> >&, string, string);

  /* RowHermPacked */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Hermitian, RowHermPacked>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Hermitian, RowHermPacked>&,
	   const Vector<complex<double> >&, string, string);

  /* ColHermPacked */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Hermitian, ColHermPacked>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Hermitian, ColHermPacked>&,
	   const Vector<complex<double> >&, string, string);
  
  /* RowSparse */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, RowSparse>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ArrayRowSparse>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ArrayRowSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, RowSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ArrayRowSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ArrayRowSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ColSparse>&,
			 const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ArrayColSparse>&,
	   const Vector<double>&, const Vector<double>&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ArrayColSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ColSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ArrayColSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ArrayColSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, RowMajor>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, RowSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, RowSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ArrayRowSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ArrayRowSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, RowSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, Symmetric, RowSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, ArrayRowSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, Symmetric, ArrayRowSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ColSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ColSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ArrayColSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ArrayColSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&, const Matrix<double, Symmetric, ColSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, Symmetric, ColSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, ArrayColSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, Symmetric, ArrayColSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  /* Vectors */
  
  SELDON_EXTERN template void 
  CheckDim(const Vector<float, VectFull>&,
	   const Vector<float, VectFull>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Vector<complex<float>, VectFull>&,
	   const Vector<complex<float>, VectFull>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Vector<float, VectSparse>&,
	   const Vector<float, VectSparse>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Vector<complex<float>, VectSparse>&,
	   const Vector<complex<float>, VectSparse>&, string, string);
  
  /* ColMajor */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, General, ColMajor>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, ColMajor>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, ColMajor>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, General, ColMajor>&,
	   const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, ColMajor>&,
	   const Vector<complex<float> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, ColMajor>&,
	   const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<float, General, ColMajor>&,
	   const Vector<float>&, const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<float>, General, ColMajor>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<float>, General, ColMajor>&,
	   const Vector<float>&, const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, General, ColMajor>&,
           const Matrix<float, General, ColMajor>&,
           const Matrix<float, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<float, General, ColMajor>&,
           const SeldonTranspose&,
           const Matrix<float, General, ColMajor>&,
           const Matrix<float, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<float, General, ColMajor>&,
           const SeldonTranspose&,
           const Matrix<float, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, General, ColMajor>&,
           const Matrix<float, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, ColMajor>&,
           const Matrix<complex<float>, General, ColMajor>&,
           const Matrix<complex<float>, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<complex<float>, General, ColMajor>&,
           const SeldonTranspose&,
           const Matrix<complex<float>, General, ColMajor>&,
           const Matrix<complex<float>, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<complex<float>, General, ColMajor>&,
           const SeldonTranspose&,
           const Matrix<complex<float>, General, ColMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, ColMajor>&,
           const Matrix<complex<float>, General, ColMajor>&, string);
  
  /* RowMajor */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, General, RowMajor>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, RowMajor>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, RowMajor>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, General, RowMajor>&,
	   const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, RowMajor>&,
	   const Vector<complex<float> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, RowMajor>&,
	   const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, General, RowMajor>&,
           const Matrix<float, General, RowMajor>&,
           const Matrix<float, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<float, General, RowMajor>&,
           const SeldonTranspose&,
           const Matrix<float, General, RowMajor>&,
           const Matrix<float, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<float, General, RowMajor>&,
           const SeldonTranspose&,
           const Matrix<float, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, General, RowMajor>&,
           const Matrix<float, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, RowMajor>&,
           const Matrix<complex<float>, General, RowMajor>&,
           const Matrix<complex<float>, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<complex<float>, General, RowMajor>&,
           const SeldonTranspose&,
           const Matrix<complex<float>, General, RowMajor>&,
           const Matrix<complex<float>, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
           const Matrix<complex<float>, General, RowMajor>&,
           const SeldonTranspose&,
           const Matrix<complex<float>, General, RowMajor>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, General, RowMajor>&,
           const Matrix<complex<float>, General, RowMajor>&, string);
  
  /* RowSymPacked */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, Symmetric, RowSymPacked>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, RowSymPacked>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, RowSymPacked>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, Symmetric, RowSymPacked>&,
	   const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, RowSymPacked>&,
	   const Vector<complex<float> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, RowSymPacked>&,
	   const Vector<float>&, string, string);

  /* ColSymPacked */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, Symmetric, ColSymPacked>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, ColSymPacked>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, ColSymPacked>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, Symmetric, ColSymPacked>&,
	   const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, ColSymPacked>&,
	   const Vector<complex<float> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, ColSymPacked>&,
	   const Vector<float>&, string, string);
  
  /* RowSym */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, Symmetric, RowSym>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, RowSym>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, RowSym>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, Symmetric, RowSym>&,
	   const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, RowSym>&,
	   const Vector<complex<float> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, RowSym>&,
	   const Vector<float>&, string, string);
  
  /* ColSym */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, Symmetric, ColSym>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, ColSym>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, ColSym>&,
	   const Vector<float>&, const Vector<float>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<float, Symmetric, ColSym>&,
	   const Vector<float>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, ColSym>&,
	   const Vector<complex<float> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Symmetric, ColSym>&,
	   const Vector<float>&, string, string);

  /* RowHerm */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Hermitian, RowHerm>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Hermitian, RowHerm>&,
	   const Vector<complex<float> >&, string, string);

  /* ColHerm */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Hermitian, ColHerm>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Hermitian, ColHerm>&,
	   const Vector<complex<float> >&, string, string);

  /* RowHermPacked */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Hermitian, RowHermPacked>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Hermitian, RowHermPacked>&,
	   const Vector<complex<float> >&, string, string);

  /* ColHermPacked */
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Hermitian, ColHermPacked>&,
	   const Vector<complex<float> >&, const Vector<complex<float> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<float>, Hermitian, ColHermPacked>&,
	   const Vector<complex<float> >&, string, string);
  
}


