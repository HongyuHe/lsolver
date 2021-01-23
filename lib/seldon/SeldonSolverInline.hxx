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


#ifndef SELDON_FILE_SELDON_SOLVER_INLINE_HXX

#include "SeldonInline.hxx"

#include "matrix_sparse/Matrix_ArraySparseInline.cxx"

#ifdef SELDON_WITH_PRECONDITIONING
#include "SeldonPreconditionerInline.hxx"
#endif

#include "computation/solver/SparseSolverInline.cxx"
#include "computation/interfaces/direct/SparseDirectSolverInline.cxx"

// iterative solvers and preconditioning
#include "computation/solver/iterative/IterativeInline.cxx"

#ifdef SELDON_WITH_VIRTUAL
#include "computation/interfaces/eigenvalue/VirtualEigenvalueSolverInline.cxx"
#endif

#ifdef SELDON_WITH_MPI
#include "SeldonDistributedInline.hxx"
#endif

#ifdef SELDON_WITH_HYPRE
#include "computation/solver/preconditioner/HypreInline.cxx"
#endif

#ifdef SELDON_WITH_ARPACK
#include "computation/interfaces/eigenvalue/ArpackInline.cxx"
#endif

#ifdef SELDON_WITH_ANASAZI
#include "computation/interfaces/eigenvalue/MyMultiVecInline.cpp"
#endif

#define SELDON_FILE_SELDON_SOLVER_INLINE_HXX
#endif
