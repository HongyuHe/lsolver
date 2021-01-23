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


#ifndef SELDON_FILE_SELDON_HXX


#include "SeldonHeader.hxx"
#include "SeldonInline.hxx"

#include "share/Allocator.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "share/Common.cxx"
#include "share/MatrixFlag.cxx"
#include "share/Errors.cxx"
#endif

#include "array/Array3D.cxx"
#include "array/Array4D.cxx"
#include "array/Array.cxx"
#include "matrix/Matrix_Base.cxx"
#include "matrix/Matrix_Pointers.cxx"
#include "matrix/Matrix_Triangular.cxx"
#include "matrix/Matrix_Symmetric.cxx"
#include "matrix/Matrix_Hermitian.cxx"
#include "matrix_sparse/Matrix_Sparse.cxx"
#include "matrix_sparse/Matrix_SymSparse.cxx"
#include "matrix/Matrix_SymPacked.cxx"
#include "matrix/Matrix_HermPacked.cxx"
#include "matrix/Matrix_TriangPacked.cxx"
#include "vector/Vector.cxx"
#include "vector/VectorCollection.cxx"
#include "vector/Functions_Arrays.cxx"
#include "vector/SparseVector.cxx"
#include "matrix/Functions.cxx"
#include "matrix_sparse/IOMatrixMarket.cxx"
#include "matrix_sparse/Matrix_Conversions.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#include "computation/basic_functions/Functions_Vector.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"
#include "computation/basic_functions/Functions_Base.cxx"

#include "vector/TinyVector.cxx"
#include "matrix/TinyMatrix.cxx"

// Blas interface.
#ifdef SELDON_WITH_BLAS
#include "computation/interfaces/Blas_1.cxx"
#include "computation/interfaces/Blas_2.cxx"
#include "computation/interfaces/Blas_3.cxx"
#endif

// Lapack interface.
#ifdef SELDON_WITH_LAPACK
#include "computation/interfaces/Lapack_LinearEquations.cxx"
#include "computation/interfaces/Lapack_LeastSquares.cxx"
#include "computation/interfaces/Lapack_Eigenvalues.cxx"
#endif // SELDON_WITH_LAPACK.

// MKL additional functions
#ifdef SELDON_WITH_MKL
#include "computation/interfaces/Mkl_Sparse.cxx"
#endif

#define SELDON_FILE_SELDON_HXX
#endif
