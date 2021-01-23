#ifndef __SELDONLIB
#define __SELDONLIB

// common flags are included
#include "SeldonFlag.hxx"

#ifdef SELDON_WITH_COMPILED_LIBRARY

// compilation FAST or TEST

// only header files are included, we assume here that all the needed functions and classes
// have been compiled during the compilation of cpp files in lib/Compil/Seldon/
#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "SeldonComplexMatrixHeader.hxx"
#include "SeldonComplexMatrixInline.hxx"

// lonely files
#include "vector/Vector2.hxx"
#include "vector/Vector3.hxx"
#include "matrix_sparse/IOMatrixMarket.hxx"

#include "matrix_sparse/BandMatrix.hxx"
#include "matrix_sparse/BandMatrixInline.cxx"

#else

// usual compilation (without library)
#include "SeldonSolver.hxx"
#include "SeldonComplexMatrix.hxx"

// lonely files
#include "vector/Vector2.cxx"
#include "vector/Vector3.cxx"
#include "matrix_sparse/IOMatrixMarket.cxx"
#include "matrix_sparse/BandMatrix.cxx"

#endif

#endif // __SELDONLIB

