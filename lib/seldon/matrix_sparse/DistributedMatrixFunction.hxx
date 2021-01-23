// Copyright (C) 2013-2015 INRIA
// Author(s): Marc Duruflé
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

#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_FUNCTION_HXX

namespace Seldon
{

  template<class T, class Prop>
  void DistributeCentralizedMatrix(Matrix<T, Prop, RowSparse>& A,
				   const MPI_Comm& comm, int root_proc,
				   const Vector<int>& num_proc,
				   DistributedMatrixIntegerArray& Binfo,
				   DistributedMatrix<T, Prop, RowSparse>& B,
				   bool keep_matrix = false);

  template<class T, class Prop>
  void DistributeCentralizedMatrix(Matrix<T, Prop, RowSymSparse>& A,
				   const MPI_Comm& comm, int root_proc,
				   const Vector<int>& num_proc,
				   DistributedMatrixIntegerArray& Binfo,
				   DistributedMatrix<T, Prop, RowSymSparse>& B,
				   bool keep_matrix = false);

  template<class T, class Prop>
  void DistributeCentralizedMatrix(Matrix<T, Prop, ArrayRowSparse>& A,
				   const MPI_Comm& comm, int root_proc,
				   const Vector<int>& num_proc,
				   DistributedMatrixIntegerArray& Binfo,
				   DistributedMatrix<T, Prop, ArrayRowSparse>& B,
				   bool keep_matrix = false);

  template<class T, class Prop>
  void DistributeCentralizedMatrix(Matrix<T, Prop, ArrayRowSymSparse>& A,
				   const MPI_Comm& comm, int root_proc,
				   const Vector<int>& num_proc,
				   DistributedMatrixIntegerArray& Binfo,
				   DistributedMatrix<T, Prop, ArrayRowSymSparse>& B,
				   bool keep_matrix = false);

}

#define SELDON_FILE_DISTRIBUTED_MATRIX_FUNCTION_HXX
#endif

