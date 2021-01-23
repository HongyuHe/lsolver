#ifndef __SELDONFLAG
#define __SELDONFLAG

#include <list>
#include <complex>

// you choose here the default allocator for non-standard types
// for basic types (such as int, float, double, etc), it is MallocAlloc
// defined by the class SeldonDefaultAllocator
#define SELDON_DEFAULT_ALLOCATOR NewAlloc

// you choose here if you want to stop the program when an exception is raised
// seldon will call abort() when encountering an exception if this flag is set
// if you want to catch exceptions, it is better to unset this flag
// this flag can also be adjusted in the Makefile
#define SELDON_WITH_ABORT

// if this flag is set, any call to FillRand will not call the function srand()
// it is assumed that the user has called srand at the beginning of its program for instance
// and he does not want to call again this function for each vector or matrix filled randomly
#define SELDON_WITHOUT_REINIT_RANDOM

// if this flag is set, some virtual functions are present in Seldon
// it is usually better because iterative solvers or eigenvalue solvers can be compiled once
// and new structures of matrices will not require a new instantiation of these solvers
#define SELDON_WITH_VIRTUAL

// if this flag is set, preconditionings present in Seldon are included
#define SELDON_WITH_PRECONDITIONING


// this flag is useful to avoid too many instantiations in Seldon when the functions Mlt/MltAdd/Add are called
// it is highly recommended to set this flag
#define SELDON_WITH_REDUCED_TEMPLATE

// float number and complex number in working precision
namespace Seldon
{
  typedef double Real_wp;
  typedef std::complex<double> Complex_wp;

  // uncomment next line if the working precision is multiple precision (MPFR for example)
  // #define SELDON_WITH_MULTIPLE
}

// this flag should be used if you cannot call zdotcsub, cdotcsub directly
// if this flag is defined the functions cblas_zdotc_sub, etc are called in DotProd (Blas_1.cxx)
// #define SELDON_WITHOUT_CBLAS_LIB

// four levels are present (1, 2, 3, 4)
// with the following checks (each next level includes checks of all previous levels)
// no flag set : no check is performed
// level 1 : Lapack errors, allocations, input/output files present
// level 2 : dimensions of vectors/matrices compatible
// level 3 : bounds of arrays (when accessing to an element of an array)
// level 4 : no additional check
// usually this flag is set in the compilation command or in the Makefile
// #define SELDON_DEBUG_LEVEL_1

// this flag is usually defined in the compilation line (-DSELDON_WITH_MPI) or in the Makefile
// #define SELDON_WITH_MPI

// Mpi is included if the flag SELDON_WITH_MPI has been set
#ifdef SELDON_WITH_MPI
#include "mpi.h"
#endif

// this flag has to be set if you want to write dense vectors/matrices in HDF5 format
// #define SELDON_WITH_HDF5

// this flag can be modified if you want to use another allocator for Vector2
// (for the second level SeldonDefaultAllocator is used)
// #define SELDON_VECTOR2_DEFAULT_ALLOCATOR_1 NewAlloc

// this flag can be modified if you want to use another allocator for Vector3 (for the two first levels)
// (for the third level SeldonDefaultAllocator is used)
// #define SELDON_VECTOR3_DEFAULT_ALLOCATOR_1 NewAlloc
// #define SELDON_VECTOR3_DEFAULT_ALLOCATOR_2 NewAlloc


// usually these two defines (Seldon will be interfaced with Blas/Lapack)
// are set in the compilation line or in the Makefile
// #define SELDON_WITH_BLAS
// #define SELDON_WITH_LAPACK


// this flag will activate the interface between Seldon and MKL
// mainly sparse Blas functions will be used for matrix-vector and matrix-matrix operations
// involving sparse matrices
// this flag can also be set in the Makefile
// #define SELDON_WITH_MKL


// this flag allows the user to use the parenthesis operator () to modify
// sparse matrices such as ArrayRowSparse and ArrayRowSymSparse
// A(i, j) = val; will be allowed for these matrices
// it is advised to not set this flag and use
// A.Get(i, j) = val
// such that the access operator () does not modify the matrix in any case
// If you have set this flag and you write if (A(i, j) == 3.0)
// the non-zero entry (i, j) will be created if it does not exist 
// #define SELDON_WITH_MODIFIABLE_PARENTHESIS_OPERATOR


// these flags can be set in the Makefile or in the compilation command
// they activate the interface between Seldon and the direct solver chosen
// #define SELDON_WITH_MUMPS
// #define SELDON_WITH_UMFPACK
// #define SELDON_WITH_SUPERLU
// #define SELDON_WITH_SUPERLU_MT 
// #define SELDON_WITH_SUPERLU_DIST 
// #define SELDON_WITH_PASTIX
// #define SELDON_WITH_PARDISO
// #define SELDON_WITH_WSMP
// #define SELDON_WITH_CHOLMOD

// for solvers compiled in 64-bits, you have to define the following flags
// #define PARDISO_INTSIZE64
// #define SUPERLU_INTSIZE64
// #define UMFPACK_INTSIZE64


// these flags can be set in the Makefile or in the compilation command
// they activate the interface between Seldon and the eigenvalue solver chosen
// #define SELDON_WITH_ARPACK
// #define SELDON_WITH_ANASAZI
// #define SELDON_WITH_FEAST

// this flag can be set in the Makefile or in the compilation command
// it activates the interface with pep solver in Slepc
// the path apply_patch_slepc.py must be executed in Slepc directory before the compilation
//#define SELDON_WITH_SLEPC_PEP

#endif // __SELDONFLAG
