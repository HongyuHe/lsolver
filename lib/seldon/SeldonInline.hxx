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


#ifndef SELDON_FILE_SELDON_INLINE_HXX


#include "share/CommonInline.cxx"

// Memory management.
#include "share/AllocatorInline.cxx"

// Storage type.
#include "share/StorageInline.cxx"
#include "share/MatrixFlagInline.cxx"

#include "array/Array3D_Inline.cxx"
#include "array/Array4D_Inline.cxx"
#include "array/ArrayInline.cxx"
#include "matrix/Matrix_BaseInline.cxx"
#include "matrix/Matrix_PointersInline.cxx"
#include "matrix/Matrix_TriangularInline.cxx"
#include "matrix/Matrix_SymmetricInline.cxx"
#include "matrix/Matrix_HermitianInline.cxx"
#include "matrix_sparse/Matrix_SparseInline.cxx"
#include "matrix_sparse/Matrix_SymSparseInline.cxx"
#include "matrix/Matrix_SymPackedInline.cxx"
#include "matrix/Matrix_HermPackedInline.cxx"
#include "matrix/Matrix_TriangPackedInline.cxx"
#include "vector/VectorInline.cxx"
#include "vector/VectorCollectionInline.cxx"
#include "vector/SparseVectorInline.cxx"

#include "matrix/SubMatrix_BaseInline.cxx"
#include "matrix/SubMatrixInline.cxx"

#include "vector/TinyVectorInline.cxx"
#include "matrix/TinyMatrixInline.cxx"

#include "computation/basic_functions/Functions_BaseInline.cxx"

#define SELDON_FILE_SELDON_INLINE_HXX
#endif
