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


#ifndef SELDON_FILE_SELDON_SOLVER_HEADER_HXX
#define SELDON_FILE_SELDON_SOLVER_HEADER_HXX

#include "SeldonHeader.hxx"

// additional classes and functions for sparse matrices
#include "matrix_sparse/Matrix_Conversions.hxx"
#include "matrix_sparse/Matrix_ArraySparse.hxx"
#include "matrix_sparse/Permutation_ScalingMatrix.hxx"
#include "matrix_sparse/Relaxation_MatVect.hxx"
#include "matrix_sparse/Functions_MatrixArray.hxx"

// iterative solvers and preconditioning
#include "computation/solver/iterative/Iterative.hxx"

// interfaces with direct solvers
#include "computation/solver/SparseSolver.hxx"

#ifdef SELDON_WITH_MUMPS
#include "computation/interfaces/direct/Mumps.hxx"
#endif

#ifdef SELDON_WITH_UMFPACK
#include "computation/interfaces/direct/UmfPack.hxx"
#endif

#ifdef SELDON_WITH_SUPERLU
#include "computation/interfaces/direct/SuperLU.hxx"
#endif

#ifdef SELDON_WITH_PASTIX
#include "computation/interfaces/direct/Pastix.hxx"
#endif

#ifdef SELDON_WITH_PARDISO
#include "computation/interfaces/direct/Pardiso.hxx"
#endif

#ifdef SELDON_WITH_WSMP
#include "computation/interfaces/direct/Wsmp.hxx"
#endif

#ifdef SELDON_WITH_PRECONDITIONING
#include "SeldonPreconditionerHeader.hxx"
#endif

#include "computation/interfaces/direct/SparseDirectSolver.hxx"

// Cholesky Solver
#ifdef SELDON_WITH_CHOLMOD
#include "computation/interfaces/direct/Cholmod.hxx"
#endif

#include "computation/solver/SparseCholeskyFactorisation.hxx"

#ifdef SELDON_WITH_MPI
#include "SeldonDistributedHeader.hxx"
#endif

#ifdef SELDON_WITH_HYPRE
#include "computation/solver/preconditioner/Hypre.hxx"
#endif

// eigenvalue stuff
#ifdef SELDON_WITH_VIRTUAL
#include "computation/interfaces/eigenvalue/VirtualEigenvalueSolver.hxx"
#include "computation/interfaces/eigenvalue/PolynomialEigenvalueSolver.hxx"
#include "computation/interfaces/eigenvalue/NonLinearEigenvalueSolver.hxx"
#else
#include "computation/interfaces/eigenvalue/EigenvalueSolver.hxx"
#endif

#ifdef SELDON_WITH_ARPACK
#include "computation/interfaces/eigenvalue/Arpack.hxx"
#include "computation/interfaces/eigenvalue/ArpackSolver.hxx"
#endif

#ifdef SELDON_WITH_ANASAZI
#include "computation/interfaces/eigenvalue/Anasazi.hxx"
#endif

#ifdef SELDON_WITH_FEAST
#include "computation/interfaces/eigenvalue/Feast.hxx"
#endif

#ifdef SELDON_WITH_SLEPC
#include "computation/interfaces/eigenvalue/Slepc.hxx"
#endif

#endif
