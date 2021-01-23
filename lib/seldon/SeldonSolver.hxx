// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_SELDON_SOLVER_HXX
#define SELDON_FILE_SELDON_SOLVER_HXX

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#include "Seldon.hxx"

// additional classes and functions for sparse matrices
#include "matrix_sparse/Matrix_Conversions.cxx"
#include "matrix_sparse/Matrix_ArraySparse.cxx"
#include "matrix_sparse/Permutation_ScalingMatrix.cxx"
#include "matrix_sparse/Relaxation_MatVect.cxx"
#include "matrix_sparse/Functions_MatrixArray.cxx"


// interfaces with direct solvers
#ifndef SELDON_WITH_COMPILED_LIBRARY
#ifdef SELDON_WITH_MUMPS
#include "computation/interfaces/direct/Mumps.cxx"
#endif

#ifdef SELDON_WITH_UMFPACK
#include "computation/interfaces/direct/UmfPack.cxx"
#endif

#ifdef SELDON_WITH_SUPERLU
#include "computation/interfaces/direct/SuperLU.cxx"
#endif

#ifdef SELDON_WITH_PASTIX
#include "computation/interfaces/direct/Pastix.cxx"
#endif

#ifdef SELDON_WITH_PARDISO
#include "computation/interfaces/direct/Pardiso.cxx"
#endif
#endif

#ifdef SELDON_WITH_WSMP
#include "computation/interfaces/direct/Wsmp.cxx"
#endif

#ifdef SELDON_WITH_PRECONDITIONING
#include "SeldonPreconditioner.hxx"
#endif

#include "computation/solver/Ordering.cxx"
#include "computation/solver/SparseSolver.cxx"
#include "computation/interfaces/direct/SparseDirectSolver.cxx"

// iterative solvers and preconditioning
#include "computation/solver/iterative/Iterative.cxx"

// Cholesky Solver
#ifndef SELDON_WITH_COMPILED_LIBRARY
#ifdef SELDON_WITH_CHOLMOD
#include "computation/interfaces/direct/Cholmod.cxx"
#endif
#endif

#include "computation/solver/SparseCholeskyFactorisation.cxx"

#ifdef SELDON_WITH_MPI
#include "SeldonDistributed.hxx"
#endif

#ifdef SELDON_WITH_HYPRE
#include "computation/solver/preconditioner/Hypre.cxx"
#endif

// eigenvalue stuff
#ifdef SELDON_WITH_ARPACK
#include "computation/interfaces/eigenvalue/Arpack.cxx"
#include "computation/interfaces/eigenvalue/ArpackSolver.cxx"
#endif

#ifdef SELDON_WITH_ANASAZI
#include "computation/interfaces/eigenvalue/MyMultiVec.cpp"
#include "computation/interfaces/eigenvalue/Anasazi.cxx"
#endif

#ifdef SELDON_WITH_FEAST
#include "computation/interfaces/eigenvalue/Feast.cxx"
#endif

#ifdef SELDON_WITH_SLEPC
#include "computation/interfaces/eigenvalue/Slepc.cxx"
#endif

#ifdef SELDON_WITH_VIRTUAL
#include "computation/interfaces/eigenvalue/VirtualEigenvalueSolver.cxx"
#include "computation/interfaces/eigenvalue/PolynomialEigenvalueSolver.cxx"
#include "computation/interfaces/eigenvalue/NonLinearEigenvalueSolver.cxx"
#else
#include "computation/interfaces/eigenvalue/EigenvalueSolver.cxx"
#endif

#endif
