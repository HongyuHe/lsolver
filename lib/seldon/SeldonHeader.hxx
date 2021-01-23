// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_SELDONHEADER_HXX

#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <stdexcept>

// For the compiled library.
#ifdef SELDON_WITH_COMPILED_LIBRARY
#define SELDON_EXTERN extern
#else
#define SELDON_EXTERN
#endif

// For backward compatibility.
#ifdef SELDON_WITH_CBLAS
#define SELDON_WITH_BLAS
#endif

// CBLAS enums are defined here to avoid SELDON_WITH_BLAS in MatrixFlag.hxx
// as a result, compiled version of Seldon is quick to switch from a Blas to a non-Blas version
extern "C"
{
  enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
  enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
  enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
  enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
  enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};
}

#ifdef SELDON_WITH_BLAS
extern "C"
{
#include "computation/interfaces/cblas.h"
}
#endif

//////////////////
// DEBUG LEVELS //
//////////////////

#ifdef SELDON_DEBUG_LEVEL_4
#ifndef SELDON_DEBUG_LEVEL_3
#define SELDON_DEBUG_LEVEL_3
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_3
#ifndef SELDON_CHECK_BOUNDS
#define SELDON_CHECK_BOUNDS
#endif
#ifndef SELDON_DEBUG_LEVEL_2
#define SELDON_DEBUG_LEVEL_2
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_2
#ifndef SELDON_CHECK_DIMENSIONS
#define SELDON_CHECK_DIMENSIONS
#endif
#ifndef SELDON_DEBUG_LEVEL_1
#define SELDON_DEBUG_LEVEL_1
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_1
#ifndef SELDON_LAPACK_CHECK_INFO
#define SELDON_LAPACK_CHECK_INFO
#endif
#ifndef SELDON_CHECK_MEMORY
#define SELDON_CHECK_MEMORY
#endif
#ifndef SELDON_CHECK_IO
#define SELDON_CHECK_IO
#endif
#ifndef SELDON_DEBUG_LEVEL_0
#define SELDON_DEBUG_LEVEL_0
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_0
#ifndef SELDON_DEBUG_LEVEL_1
#ifndef SELDON_WITHOUT_THROW
#define SELDON_WITHOUT_THROW
#endif
#endif
#endif

// Convenient macros to catch exceptions.
#ifndef TRY
#define TRY try {
#endif
#ifndef END
#define END                                                             \
  }                                                                     \
    catch(Seldon::Error& Err)                                           \
      {                                                                 \
        Err.CoutWhat();                                                 \
        return 1;                                                       \
      }                                                                 \
    catch (std::exception& Err)                                         \
      {                                                                 \
        std::cout << "C++ exception: " << Err.what() << std::endl;      \
        return 1;                                                       \
      }                                                                 \
    catch (std::string& str)                                            \
      {                                                                 \
        std::cout << str << std::endl;                                  \
        return 1;                                                       \
      }                                                                 \
    catch (const char* str)                                             \
      {                                                                 \
        std::cout << str << std::endl;                                  \
        return 1;                                                       \
      }                                                                 \
    catch(...)                                                          \
      {                                                                 \
        std::cout << "Unknown exception..." << std::endl;               \
        return 1;                                                       \
      }
#endif

//! To display a message... call Hermes!
#ifndef ERR
#define ERR(x) std::cout << "Hermes - " #x << std::endl
#endif
//! To display a variable (with its name); same as DISPLAY.
#ifndef DISP
#define DISP(x) std::cout << #x ": " << x << std::endl
#endif
//! To display a variable (with its name); same as DISP.
#ifndef DISPLAY
#define DISPLAY(x) std::cout << #x ": " << x << std::endl
#endif

// For backward compatibility. These lines should be removed one day.
#define Vect_Full VectFull
#define Vect_Sparse VectSparse

//! Seldon namespace.
namespace Seldon
{
  using namespace std;
}

// Exceptions and useful functions.
#include "share/Errors.hxx"
#include "share/Common.hxx"

// Default allocator.
#ifndef SELDON_DEFAULT_ALLOCATOR
#define SELDON_DEFAULT_ALLOCATOR NewAlloc
#endif
// Memory management.
#include "share/Allocator.hxx"
#include "share/DefaultAllocator.hxx"

// Storage type.
#include "share/Storage.hxx"

#include "share/MatrixFlag.hxx"

// Properties.
#include "share/Properties.hxx"


namespace Seldon
{


  // Base structure for all vectors.
  template <class T, class Allocator>
  class Vector_Base;

  // Vector class - specialized for each used type.
  template <class T, class Storage = VectFull,
	    class Allocator = typename
	    SeldonDefaultAllocator<Storage, T>::allocator >
  class Vector;

  // Full vector.
  template <class T, class Allocator>
  class Vector<T, VectFull, Allocator>;

  // Sparse vector.
  template <class T, class Allocator>
  class Vector<T, VectSparse, Allocator>;

  // Vector collection.
  template <class T, class Allocator>
  class Vector<T, Collection, Allocator>;

  // Matrix class - specialized for each used type.
  template <class T, class Prop = General,
	    class Storage = RowMajor,
	    class Allocator = typename
	    SeldonDefaultAllocator<Storage, T>::allocator >
  class Matrix;

  // column-major matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColMajor, Allocator>;

  // row-major matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowMajor, Allocator>;

  // column-major symmetric packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymPacked, Allocator>;

  // row-major symmetric packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymPacked, Allocator>;

  // column-major upper-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColUpTriangPacked, Allocator>;

  // column-major lower-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColLoTriangPacked, Allocator>;

  // row-major upper-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowUpTriangPacked, Allocator>;

  // row-major lower-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowLoTriangPacked, Allocator>;

  // column-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSparse, Allocator>;

  // row-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSparse, Allocator>;

  // column-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymSparse, Allocator>;

  // row-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymSparse, Allocator>;

  // column-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSparse, Allocator>;

  // row-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSparse, Allocator>;

  // column-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSymSparse, Allocator>;

  // row-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymSparse, Allocator>;

  // 3D array.
  template <class T, int N, class Allocator
	    = typename SeldonDefaultAllocator<VectFull, T>::allocator>
  class Array;

  // 3D array.
  template <class T, class Allocator
	    = typename SeldonDefaultAllocator<VectFull, T>::allocator>
  class Array3D;

  // 4D array.
  template <class T, class Allocator
	    = typename SeldonDefaultAllocator<VectFull, T>::allocator>
  class Array4D;


} // namespace Seldon.


const int ARRAY_MINRANK = 3;
const int ARRAY_MAXRANK = 9;

#include "array/Array3D.hxx"
#include "array/Array4D.hxx"
#include "array/Array.hxx"
#include "matrix/Matrix_Base.hxx"
#include "matrix/Matrix_Pointers.hxx"
#include "matrix/Matrix_Triangular.hxx"
#include "matrix/Matrix_Symmetric.hxx"
#include "matrix/Matrix_Hermitian.hxx"
#include "matrix_sparse/Matrix_Sparse.hxx"
#include "matrix_sparse/Matrix_SymSparse.hxx"
#include "matrix/Matrix_SymPacked.hxx"
#include "matrix/Matrix_HermPacked.hxx"
#include "matrix/Matrix_TriangPacked.hxx"
#include "vector/Vector.hxx"
#include "vector/SparseVector.hxx"
#include "vector/Functions_Arrays.hxx"
#include "matrix/Functions.hxx"
#include "matrix_sparse/IOMatrixMarket.hxx"
#include "matrix_sparse/Matrix_Conversions.hxx"
#include "computation/basic_functions/Functions_Vector.hxx"
#include "computation/basic_functions/Functions_MatVect.hxx"
#include "computation/basic_functions/Functions_Matrix.hxx"
#include "computation/basic_functions/Functions_Base.hxx"

#include "matrix/SubMatrix_Base.hxx"
#include "matrix/SubMatrix.hxx"

// Blas interface.
#ifdef SELDON_WITH_BLAS

#include "computation/interfaces/Blas_1.hxx"
#include "computation/interfaces/Blas_2.hxx"
#include "computation/interfaces/Blas_3.hxx"

#endif

// Lapack interface.
#ifdef SELDON_WITH_LAPACK

#undef LAPACK_INTEGER
#define LAPACK_INTEGER int
#undef LAPACK_REAL
#define LAPACK_REAL float
#undef LAPACK_DOUBLEREAL
#define LAPACK_DOUBLEREAL double
#undef LAPACK_COMPLEX
#define LAPACK_COMPLEX void
#undef LAPACK_DOUBLECOMPLEX
#define LAPACK_DOUBLECOMPLEX void
#undef LAPACK_LOGICAL
#define LAPACK_LOGICAL int
#undef LAPACK_L_FP
#define LAPACK_L_FP int*
#undef LAPACK_FTNLEN
#define LAPACK_FTNLEN int*
extern "C"
{
#include "computation/interfaces/clapack.h"
}
#ifdef SELDON_LAPACK_CHECK_INFO
#ifndef SELDON_CHECK_INFO
#define SELDON_CHECK_INFO(f, lf) info.Check(f, lf)
#endif
#else
#ifndef SELDON_CHECK_INFO
#define SELDON_CHECK_INFO(f, lf)
#endif
#endif

#include "computation/interfaces/Lapack_LinearEquations.hxx"
#include "computation/interfaces/Lapack_LeastSquares.hxx"
#include "computation/interfaces/Lapack_Eigenvalues.hxx"

#endif // SELDON_WITH_LAPACK.

#ifdef SELDON_WITH_MKL
#include "computation/interfaces/Mkl_Sparse.hxx"
#endif

namespace Seldon
{

  typedef Array3D<int> IArr3D;
  typedef Vector<int, VectFull> IVect;
  typedef Vector<float, VectFull> SVect;
  typedef Vector<double, VectFull> DVect;
  typedef Vector<complex<float>, VectFull> CVect;
  typedef Vector<complex<double>, VectFull> ZVect;

  typedef Matrix<int, General, ColMajor> IGCMat;
  typedef Matrix<float, General, ColMajor> SGCMat;
  typedef Matrix<double, General, ColMajor> DGCMat;
  typedef Matrix<complex<float>, General, ColMajor> CGCMat;
  typedef Matrix<complex<double>, General, ColMajor> ZGCMat;

  typedef Matrix<int, General, RowMajor> IGRMat;
  typedef Matrix<float, General, RowMajor> SGRMat;
  typedef Matrix<double, General, RowMajor> DGRMat;
  typedef Matrix<complex<float>, General, RowMajor> CGRMat;
  typedef Matrix<complex<double>, General, RowMajor> ZGRMat;

  typedef Matrix<int, General, RowSparse> IGRSMat;
  typedef Matrix<float, General, RowSparse> SGRSMat;
  typedef Matrix<double, General, RowSparse> DGRSMat;
  typedef Matrix<complex<float>, General, RowSparse> CGRSMat;
  typedef Matrix<complex<double>, General, RowSparse> ZGRSMat;

  typedef Matrix<int, General, ColSparse> IGCSMat;
  typedef Matrix<float, General, ColSparse> SGCSMat;
  typedef Matrix<double, General, ColSparse> DGCSMat;
  typedef Matrix<complex<float>, General, ColSparse> CGCSMat;
  typedef Matrix<complex<double>, General, ColSparse> ZGCSMat;

  // Vector class - specialized for each used type.
  template <class T, class Storage, class Allocator>
  class Vector
  {
    // Nothing in it: no default vector is supplied so as to avoid suprises!
  };


  // Matrix class - specialized for each used type.
  template <class T, class Prop, class Storage, class Allocator>
  class Matrix
  {
    // Nothing in it: no default matrix is supplied so as to avoid suprises!
  };

} // namespace Seldon.

// small matrices and vectors
#include "vector/TinyVector.hxx"
#include "matrix/TinyMatrix.hxx"

#define SELDON_FILE_SELDONHEADER_HXX
#endif
