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

#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_FUNCTION_CXX

namespace Seldon
{

  template<class T, class Prop>
  void DistributeCentralizedMatrix(Matrix<T, Prop, RowSparse>& A,
				   const MPI_Comm& comm, int root_proc,
				   const Vector<int>& num_proc,
				   DistributedMatrixIntegerArray& Binfo,
				   DistributedMatrix<T, Prop, RowSparse>& B,
				   bool keep_matrix)
  {
    Binfo.Clear();
    B.Clear();
    
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();

    Vector<int> numB(3), IndColB, NumProcB;
    Vector<T> ValB; Vector<long> PtrB;
    Vector<int64_t> xtmp;
    if (rank_proc == root_proc)
      {
	Vector<int> glob_to_local(A.GetM());
	glob_to_local.Fill(-1);
	
	// we count the number of rows for each proc
	// and construct glob_to_local
	Vector<int> nb_rows_proc(nb_proc);
	nb_rows_proc.Zero();
	for (int i = 0; i < num_proc.GetM(); i++)
	  glob_to_local(i) = nb_rows_proc(num_proc(i))++;

	// storing row numbers for all the processors
	Vector<Vector<int> > all_rows(nb_proc);
	for (int i = 0; i < nb_proc; i++)
	  all_rows(i).Reallocate(nb_rows_proc(i));

	nb_rows_proc.Zero();
	for (int i = 0; i < num_proc.GetM(); i++)
	  {
	    all_rows(num_proc(i))(nb_rows_proc(num_proc(i))) = i;
	    nb_rows_proc(num_proc(i))++;
	  }
	
	Vector<long> Ptr; Vector<int> IndCol, NumProc; Vector<T> Val;
	for (int p = 0; p < nb_proc; p++)
	  {
	    // we count the number of values to send
	    int nb_rows = nb_rows_proc(p);
	    Ptr.Reallocate(nb_rows+1);
	    long nb_val = 0;
	    Ptr(0) = 0;
	    for (int j = 0; j < all_rows(p).GetM(); j++)
	      {
		int i = all_rows(p)(j);
		nb_val += ptr[i+1] - ptr[i];
		Ptr(j+1) = nb_val;
	      }
	    
	    // then we store the values to send
	    NumProc.Reallocate(nb_val);
	    IndCol.Reallocate(nb_val);
	    Val.Reallocate(nb_val);
	    nb_val = 0;
	    for (int j = 0; j < all_rows(p).GetM(); j++)
	      {
		int i = all_rows(p)(j);
		for (long k = ptr[i]; k < ptr[i+1]; k++)
		  {
		    NumProc(nb_val) = num_proc(ind[k]);
		    if (NumProc(nb_val) == p)
		      IndCol(nb_val) = glob_to_local(ind[k]);
		    else
		      IndCol(nb_val) = ind[k];
		    
		    Val(nb_val) = data[k];
		    nb_val++;
		  }
	      }

	    Vector<int> num(3);
	    num(0) = nb_rows; num(1) = nb_val; num(2) = A.GetM();

	    // the values are sent to the appropriate processor
	    if (p == root_proc)
	      {
		numB = num; Binfo.GlobalRowNumbers = all_rows(p);
		PtrB = Ptr; NumProcB = NumProc;
		IndColB = IndCol;
		ValB = Val; 
	      }
	    else
	      {
		MpiSsend(comm, num, xtmp, 3, p, 19);
		MpiSsend(comm, all_rows(p), xtmp, nb_rows, p, 20);
		MpiSsend(comm, Ptr, xtmp, Ptr.GetM(), p, 21);
		MpiSsend(comm, NumProc, xtmp, NumProc.GetM(), p, 22);
		MpiSsend(comm, IndCol, xtmp, IndCol.GetM(), p, 23);
		MpiSsend(comm, Val, xtmp, Val.GetM(), p, 24);
	      }
	  }
      }
    else
      {
	MPI_Status status;
	MpiRecv(comm, numB, xtmp, 3, root_proc, 19, status);

	Binfo.GlobalRowNumbers.Reallocate(numB(0));
	PtrB.Reallocate(numB(0)+1);
	NumProcB.Reallocate(numB(1));
	IndColB.Reallocate(numB(1));
	ValB.Reallocate(numB(1));
	MpiRecv(comm, Binfo.GlobalRowNumbers, xtmp, numB(0), root_proc, 20, status);
	MpiRecv(comm, PtrB, xtmp, PtrB.GetM(), root_proc, 21, status);
	MpiRecv(comm, NumProcB, xtmp, NumProcB.GetM(), root_proc, 22, status);
	MpiRecv(comm, IndColB, xtmp, IndColB.GetM(), root_proc, 23, status);
	MpiRecv(comm, ValB, xtmp, ValB.GetM(), root_proc, 24, status);
      }

    // centralized matrix is destructed if required
    if (!keep_matrix)
      A.Clear();

    // filling Binfo
    Binfo.nodl_scalar = numB(0); Binfo.nb_unknowns_scal = 1;
    Binfo.nglob = numB(2); Binfo.comm = comm;
    
    // then we construct the distributed matrix from arrays
    // PtrB, NumProcB, IndColB and ValB
    
    // we initialize the distributed matrix B with arrays contained in Binfo
    B.Init(Binfo);
    B.Reallocate(numB(0), numB(0));
    
    // we extract the local part of the matrix
    long nb_local_values = 0;
    for (int i = 0; i < NumProcB.GetM(); i++)
      if (NumProcB(i) == rank_proc)
	nb_local_values++;

    Vector<long> PtrLocB(numB(0)+1);
    Vector<int> IndLocB(nb_local_values);
    Vector<T> ValLocB(nb_local_values);
    PtrLocB.Zero();
    
    nb_local_values = 0; PtrLocB(0) = 0;
    for (int i = 0; i < numB(0); i++)
      {
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) == rank_proc)
	    {
	      IndLocB(nb_local_values) = IndColB(j);
	      ValLocB(nb_local_values) = ValB(j);
	      nb_local_values++;
	    }

	PtrLocB(i+1) = nb_local_values;
      }

    B.SetData(numB(0), numB(0), ValLocB, PtrLocB, IndLocB);
    B.ReallocateDist(numB(0), numB(0));
    
    // then we add the distant part
    for (int i = 0; i < numB(0); i++)
      {
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) != rank_proc)
	    B.AddDistantInteraction(i, IndColB(j), NumProcB(j), ValB(j));
      }    
  }

  
  template<class T, class Prop>
  void DistributeCentralizedMatrix(Matrix<T, Prop, RowSymSparse>& A,
				   const MPI_Comm& comm, int root_proc,
				   const Vector<int>& num_proc,
				   DistributedMatrixIntegerArray& Binfo,
				   DistributedMatrix<T, Prop, RowSymSparse>& B,
				   bool keep_matrix)
  {
    Binfo.Clear();
    B.Clear();
    
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);

    long* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();

    Vector<int> numB(3), IndColB, NumProcB;
    Vector<T> ValB; Vector<long> PtrB;
    Vector<int64_t> xtmp;
    if (rank_proc == root_proc)
      {
	Vector<int> glob_to_local(A.GetM());
	glob_to_local.Fill(-1);
	
	// we count the number of rows for each proc
	// and construct glob_to_local
	Vector<int> nb_rows_proc(nb_proc);
	nb_rows_proc.Zero();
	for (int i = 0; i < num_proc.GetM(); i++)
	  glob_to_local(i) = nb_rows_proc(num_proc(i))++;

	// storing row numbers for all the processors
	Vector<Vector<int> > all_rows(nb_proc);
	for (int i = 0; i < nb_proc; i++)
	  all_rows(i).Reallocate(nb_rows_proc(i));

	nb_rows_proc.Zero();
	for (int i = 0; i < num_proc.GetM(); i++)
	  {
	    all_rows(num_proc(i))(nb_rows_proc(num_proc(i))) = i;
	    nb_rows_proc(num_proc(i))++;
	  }
	
	Vector<long> Ptr; Vector<int> IndCol, NumProc; Vector<T> Val;
	for (int p = 0; p < nb_proc; p++)
	  {
	    // we count the number of values to send
	    int nb_rows = nb_rows_proc(p);
	    Ptr.Reallocate(nb_rows+1);
	    long nb_val = 0;
	    Ptr(0) = 0;
	    for (int j = 0; j < all_rows(p).GetM(); j++)
	      {
		int i = all_rows(p)(j);
		nb_val += ptr[i+1] - ptr[i];
		Ptr(j+1) = nb_val;
	      }
	    
	    // then we store the values to send
	    NumProc.Reallocate(nb_val);
	    IndCol.Reallocate(nb_val);
	    Val.Reallocate(nb_val);
	    nb_val = 0;
	    for (int j = 0; j < all_rows(p).GetM(); j++)
	      {
		int i = all_rows(p)(j);
		for (long k = ptr[i]; k < ptr[i+1]; k++)
		  {
		    NumProc(nb_val) = num_proc(ind[k]);
		    if (NumProc(nb_val) == p)
		      IndCol(nb_val) = glob_to_local(ind[k]);
		    else
		      IndCol(nb_val) = ind[k];
		    
		    Val(nb_val) = data[k];
		    nb_val++;
		  }
	      }

	    Vector<int> num(3);
	    num(0) = nb_rows; num(1) = nb_val; num(2) = A.GetM();

	    // the values are sent to the appropriate processor
	    if (p == root_proc)
	      {
		numB = num; Binfo.GlobalRowNumbers = all_rows(p);
		PtrB = Ptr; NumProcB = NumProc;
		IndColB = IndCol;
		ValB = Val; 
	      }
	    else
	      {
		MpiSsend(comm, num, xtmp, 3, p, 19);
		MpiSsend(comm, all_rows(p), xtmp, nb_rows, p, 20);
		MpiSsend(comm, Ptr, xtmp, Ptr.GetM(), p, 21);
		MpiSsend(comm, NumProc, xtmp, NumProc.GetM(), p, 22);
		MpiSsend(comm, IndCol, xtmp, IndCol.GetM(), p, 23);
		MpiSsend(comm, Val, xtmp, Val.GetM(), p, 24);
	      }
	  }
      }
    else
      {
	MPI_Status status;
	MpiRecv(comm, numB, xtmp, 3, root_proc, 19, status);

	Binfo.GlobalRowNumbers.Reallocate(numB(0));
	PtrB.Reallocate(numB(0)+1);
	NumProcB.Reallocate(numB(1));
	IndColB.Reallocate(numB(1));
	ValB.Reallocate(numB(1));
	MpiRecv(comm, Binfo.GlobalRowNumbers, xtmp, numB(0), root_proc, 20, status);
	MpiRecv(comm, PtrB, xtmp, PtrB.GetM(), root_proc, 21, status);
	MpiRecv(comm, NumProcB, xtmp, NumProcB.GetM(), root_proc, 22, status);
	MpiRecv(comm, IndColB, xtmp, IndColB.GetM(), root_proc, 23, status);
	MpiRecv(comm, ValB, xtmp, ValB.GetM(), root_proc, 24, status);
      }

    // centralized matrix is destructed if required
    if (!keep_matrix)
      A.Clear();

    // filling Binfo
    Binfo.nodl_scalar = numB(0); Binfo.nb_unknowns_scal = 1;
    Binfo.nglob = numB(2); Binfo.comm = comm;
    
    // then we construct the distributed matrix from arrays
    // PtrB, NumProcB, IndColB and ValB
    
    // we initialize the distributed matrix B with arrays contained in Binfo
    B.Init(Binfo);
    B.Reallocate(numB(0), numB(0));
    
    // we extract the local part of the matrix
    long nb_local_values = 0;
    for (int i = 0; i < NumProcB.GetM(); i++)
      if (NumProcB(i) == rank_proc)
	nb_local_values++;

    Vector<long> PtrLocB(numB(0)+1); Vector<int> IndLocB(nb_local_values);
    Vector<T> ValLocB(nb_local_values);
    PtrLocB.Zero();
    
    nb_local_values = 0; PtrLocB(0) = 0;
    for (int i = 0; i < numB(0); i++)
      {
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) == rank_proc)
	    {
	      IndLocB(nb_local_values) = IndColB(j);
	      ValLocB(nb_local_values) = ValB(j);
	      nb_local_values++;
	    }

	PtrLocB(i+1) = nb_local_values;
      }

    B.SetData(numB(0), numB(0), ValLocB, PtrLocB, IndLocB);
    B.ReallocateDist(numB(0), numB(0));
    
    // then we add the distant part
    for (int i = 0; i < numB(0); i++)
      {
	for (int j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) != rank_proc)
	    B.AddDistantInteraction(i, IndColB(j), NumProcB(j), ValB(j));
      }    
  }


  template<class T, class Prop>
  void DistributeCentralizedMatrix(Matrix<T, Prop, ArrayRowSparse>& A,
				   const MPI_Comm& comm, int root_proc,
				   const Vector<int>& num_proc,
				   DistributedMatrixIntegerArray& Binfo,
				   DistributedMatrix<T, Prop, ArrayRowSparse>& B,
				   bool keep_matrix)
  {
    Binfo.Clear();
    B.Clear();
    
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);

    Vector<int> numB(3), IndColB, NumProcB;
    Vector<T> ValB; Vector<long> PtrB;
    Vector<int64_t> xtmp;
    if (rank_proc == root_proc)
      {
	Vector<int> glob_to_local(A.GetM());
	glob_to_local.Fill(-1);
	
	// we count the number of rows for each proc
	// and construct glob_to_local
	Vector<int> nb_rows_proc(nb_proc);
	nb_rows_proc.Zero();
	for (int i = 0; i < num_proc.GetM(); i++)
	  glob_to_local(i) = nb_rows_proc(num_proc(i))++;

	// storing row numbers for all the processors
	Vector<Vector<int> > all_rows(nb_proc);
	for (int i = 0; i < nb_proc; i++)
	  all_rows(i).Reallocate(nb_rows_proc(i));

	nb_rows_proc.Zero();
	for (int i = 0; i < num_proc.GetM(); i++)
	  {
	    all_rows(num_proc(i))(nb_rows_proc(num_proc(i))) = i;
	    nb_rows_proc(num_proc(i))++;
	  }
	
	Vector<long> Ptr; Vector<int> IndCol, NumProc; Vector<T> Val;
	for (int p = 0; p < nb_proc; p++)
	  {
	    // we count the number of values to send
	    int nb_rows = nb_rows_proc(p);
	    Ptr.Reallocate(nb_rows+1);
	    long nb_val = 0;
	    Ptr(0) = 0;
	    for (int j = 0; j < all_rows(p).GetM(); j++)
	      {
		int i = all_rows(p)(j);
		nb_val += A.GetRowSize(i);
		Ptr(j+1) = nb_val;
	      }
	    
	    // then we store the values to send
	    NumProc.Reallocate(nb_val);
	    IndCol.Reallocate(nb_val);
	    Val.Reallocate(nb_val);
	    nb_val = 0;
	    for (int j = 0; j < all_rows(p).GetM(); j++)
	      {
		int i = all_rows(p)(j);
		for (int k = 0; k < A.GetRowSize(i); k++)
		  {
		    NumProc(nb_val) = num_proc(A.Index(i, k));
		    if (NumProc(nb_val) == p)
		      IndCol(nb_val) = glob_to_local(A.Index(i, k));
		    else
		      IndCol(nb_val) = A.Index(i, k);
		    
		    Val(nb_val) = A.Value(i, k);
		    nb_val++;
		  }
	      }
	    
	    Vector<int> num(3);
	    num(0) = nb_rows; num(1) = nb_val; num(2) = A.GetM();

	    // the values are sent to the appropriate processor
	    if (p == root_proc)
	      {
		numB = num; Binfo.GlobalRowNumbers = all_rows(p);
		PtrB = Ptr; NumProcB = NumProc;
		IndColB = IndCol;
		ValB = Val; 
	      }
	    else
	      {
		MpiSsend(comm, num, xtmp, 3, p, 19);
		MpiSsend(comm, all_rows(p), xtmp, nb_rows, p, 20);
		MpiSsend(comm, Ptr, xtmp, Ptr.GetM(), p, 21);
		MpiSsend(comm, NumProc, xtmp, NumProc.GetM(), p, 22);
		MpiSsend(comm, IndCol, xtmp, IndCol.GetM(), p, 23);
		MpiSsend(comm, Val, xtmp, Val.GetM(), p, 24);
	      }
	  }
      }
    else
      {
	MPI_Status status;
	MpiRecv(comm, numB, xtmp, 3, root_proc, 19, status);

	Binfo.GlobalRowNumbers.Reallocate(numB(0));
	PtrB.Reallocate(numB(0)+1);
	NumProcB.Reallocate(numB(1));
	IndColB.Reallocate(numB(1));
	ValB.Reallocate(numB(1));
	MpiRecv(comm, Binfo.GlobalRowNumbers, xtmp, numB(0), root_proc, 20, status);
	MpiRecv(comm, PtrB, xtmp, PtrB.GetM(), root_proc, 21, status);
	MpiRecv(comm, NumProcB, xtmp, NumProcB.GetM(), root_proc, 22, status);
	MpiRecv(comm, IndColB, xtmp, IndColB.GetM(), root_proc, 23, status);
	MpiRecv(comm, ValB, xtmp, ValB.GetM(), root_proc, 24, status);
      }

    // centralized matrix is destructed if required
    if (!keep_matrix)
      A.Clear();

    // filling Binfo
    Binfo.nodl_scalar = numB(0); Binfo.nb_unknowns_scal = 1;
    Binfo.nglob = numB(2); Binfo.comm = comm;
    
    // then we construct the distributed matrix from arrays
    // PtrB, NumProcB, IndColB and ValB
    
    // we initialize the distributed matrix B with arrays contained in Binfo
    B.Init(Binfo);
    B.Reallocate(numB(0), numB(0));
    
    // we extract the local part of the matrix
    for (int i = 0; i < numB(0); i++)
      {
	int size_loc = 0;
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) == rank_proc)
	    size_loc++;

	B.ReallocateRow(i, size_loc);
	size_loc = 0;
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) == rank_proc)
	    {
	      B.Index(i, size_loc) = IndColB(j);
	      B.Value(i, size_loc) = ValB(j);
	      size_loc++;
	    }
      }
    
    // then we add the distant part
    for (int i = 0; i < numB(0); i++)
      {
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) != rank_proc)
	    B.AddDistantInteraction(i, IndColB(j), NumProcB(j), ValB(j));
      }    
  }


  template<class T, class Prop>
  void DistributeCentralizedMatrix(Matrix<T, Prop, ArrayRowSymSparse>& A,
				   const MPI_Comm& comm, int root_proc,
				   const Vector<int>& num_proc,
				   DistributedMatrixIntegerArray& Binfo,
				   DistributedMatrix<T, Prop, ArrayRowSymSparse>& B,
				   bool keep_matrix)
  {
    Binfo.Clear();
    B.Clear();
    
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);

    Vector<int> numB(3), IndColB, NumProcB;
    Vector<T> ValB; Vector<long> PtrB;
    Vector<int64_t> xtmp;
    if (rank_proc == root_proc)
      {
	Vector<int> glob_to_local(A.GetM());
	glob_to_local.Fill(-1);
	
	// we count the number of rows for each proc
	// and construct glob_to_local
	Vector<int> nb_rows_proc(nb_proc);
	nb_rows_proc.Zero();
	for (int i = 0; i < num_proc.GetM(); i++)
	  glob_to_local(i) = nb_rows_proc(num_proc(i))++;

	// storing row numbers for all the processors
	Vector<Vector<int> > all_rows(nb_proc);
	for (int i = 0; i < nb_proc; i++)
	  all_rows(i).Reallocate(nb_rows_proc(i));

	nb_rows_proc.Zero();
	for (int i = 0; i < num_proc.GetM(); i++)
	  {
	    all_rows(num_proc(i))(nb_rows_proc(num_proc(i))) = i;
	    nb_rows_proc(num_proc(i))++;
	  }
	
	Vector<long> Ptr; Vector<int> IndCol, NumProc; Vector<T> Val;
	for (int p = 0; p < nb_proc; p++)
	  {
	    // we count the number of values to send
	    int nb_rows = nb_rows_proc(p);
	    Ptr.Reallocate(nb_rows+1);
	    long nb_val = 0;
	    Ptr(0) = 0;
	    for (int j = 0; j < all_rows(p).GetM(); j++)
	      {
		int i = all_rows(p)(j);
		nb_val += A.GetRowSize(i);
		Ptr(j+1) = nb_val;
	      }
	    
	    // then we store the values to send
	    NumProc.Reallocate(nb_val);
	    IndCol.Reallocate(nb_val);
	    Val.Reallocate(nb_val);
	    nb_val = 0;
	    for (int j = 0; j < all_rows(p).GetM(); j++)
	      {
		int i = all_rows(p)(j);
		for (int k = 0; k < A.GetRowSize(i); k++)
		  {
		    NumProc(nb_val) = num_proc(A.Index(i, k));
		    if (NumProc(nb_val) == p)
		      IndCol(nb_val) = glob_to_local(A.Index(i, k));
		    else
		      IndCol(nb_val) = A.Index(i, k);
		    
		    Val(nb_val) = A.Value(i, k);
		    nb_val++;
		  }
	      }
	    
	    Vector<int> num(3);
	    num(0) = nb_rows; num(1) = nb_val; num(2) = A.GetM();

	    // the values are sent to the appropriate processor
	    if (p == root_proc)
	      {
		numB = num; Binfo.GlobalRowNumbers = all_rows(p);
		PtrB = Ptr; NumProcB = NumProc;
		IndColB = IndCol;
		ValB = Val; 
	      }
	    else
	      {
		MpiSsend(comm, num, xtmp, 3, p, 19);
		MpiSsend(comm, all_rows(p), xtmp, nb_rows, p, 20);
		MpiSsend(comm, Ptr, xtmp, Ptr.GetM(), p, 21);
		MpiSsend(comm, NumProc, xtmp, NumProc.GetM(), p, 22);
		MpiSsend(comm, IndCol, xtmp, IndCol.GetM(), p, 23);
		MpiSsend(comm, Val, xtmp, Val.GetM(), p, 24);
	      }
	  }
      }
    else
      {
	MPI_Status status;
	MpiRecv(comm, numB, xtmp, 3, root_proc, 19, status);

	Binfo.GlobalRowNumbers.Reallocate(numB(0));
	PtrB.Reallocate(numB(0)+1);
	NumProcB.Reallocate(numB(1));
	IndColB.Reallocate(numB(1));
	ValB.Reallocate(numB(1));
	MpiRecv(comm, Binfo.GlobalRowNumbers, xtmp, numB(0), root_proc, 20, status);
	MpiRecv(comm, PtrB, xtmp, PtrB.GetM(), root_proc, 21, status);
	MpiRecv(comm, NumProcB, xtmp, NumProcB.GetM(), root_proc, 22, status);
	MpiRecv(comm, IndColB, xtmp, IndColB.GetM(), root_proc, 23, status);
	MpiRecv(comm, ValB, xtmp, ValB.GetM(), root_proc, 24, status);
      }

    // centralized matrix is destructed if required
    if (!keep_matrix)
      A.Clear();

    // filling Binfo
    Binfo.nodl_scalar = numB(0); Binfo.nb_unknowns_scal = 1;
    Binfo.nglob = numB(2); Binfo.comm = comm;
    
    // then we construct the distributed matrix from arrays
    // PtrB, NumProcB, IndColB and ValB
    
    // we initialize the distributed matrix B with arrays contained in Binfo
    B.Init(Binfo);
    B.Reallocate(numB(0), numB(0));
    
    // we extract the local part of the matrix
    for (int i = 0; i < numB(0); i++)
      {
	int size_loc = 0;
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) == rank_proc)
	    size_loc++;

	B.ReallocateRow(i, size_loc);
	size_loc = 0;
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) == rank_proc)
	    {
	      B.Index(i, size_loc) = IndColB(j);
	      B.Value(i, size_loc) = ValB(j);
	      size_loc++;
	    }
      }
    
    // then we add the distant part
    for (int i = 0; i < numB(0); i++)
      {
	for (long j = PtrB(i); j < PtrB(i+1); j++)
	  if (NumProcB(j) != rank_proc)
	    B.AddDistantInteraction(i, IndColB(j), NumProcB(j), ValB(j));
      }    
  }
  
}

#define SELDON_FILE_DISTRIBUTED_MATRIX_FUNCTION_CXX
#endif

