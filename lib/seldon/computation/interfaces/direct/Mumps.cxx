// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_MUMPS_CXX

#include "Mumps.hxx"

extern "C"
{
  // including mpi from sequential version of Mumps if the
  // compilation is not made on a parallel machine
#ifndef SELDON_WITH_MPI
#include "mpi.h"
#endif
}

namespace Seldon
{

  //! Mumps is called in double precision
  template<>
  void MatrixMumps<double>::CallMumps()
  {
    dmumps_c(&struct_mumps);
  }


  //! Mumps is called in complex double precision
  template<>
  void MatrixMumps<complex<double> >::CallMumps()
  {
    zmumps_c(&struct_mumps);
  }


  //! Default constructor
  template<class T>
  MatrixMumps<T>::MatrixMumps()
  {
    struct_mumps.comm_fortran = -987654;

    // parameters for mumps
    struct_mumps.job = -1;
    struct_mumps.par = 1;
    struct_mumps.sym = 0; // 0 -> unsymmetric matrix

    struct_mumps.info[8] = 0;
    struct_mumps.info[9] = 0;

    // other parameters
    struct_mumps.n = 0;
    type_ordering = 7; // default : we let Mumps choose the ordering
    print_level = -1;
    info_facto = 0;
    out_of_core = false;
    parallel_ordering = false;
    coef_overestimate = 1.3;
    coef_increase_memory = 1.5;
    coef_max_overestimate = 50.0;
    threshold_pivot = 0.01;
  }


  template<class T>
  bool MatrixMumps<T>::UseInteger8() const
  {
    // Blas is interfaced with 32 bits integers
    // so we consider that Mumps is compiled with usual integers
    return false;
  }


  //! Function used to force factorisation when estimated space was too small
  template<class T>
  void MatrixMumps<T>::IterateFacto()
  {
    /*DISP(struct_mumps.info[0]);
    DISP(struct_mumps.icntl[23]);
    DISP(struct_mumps.cntl[2]);
    if (struct_mumps.info[0] == -10)
      {
        int npiv = struct_mumps.infog[27];
        cout << "npiv =" << npiv << endl;
        for (int i = 0; i < npiv; i++)
          DISP(struct_mumps.pivnul_list[i]);
          } */

    // if error -9 occurs, retrying with larger size
    int init_percentage = struct_mumps.icntl[13];
    int new_percentage = init_percentage;
    while (((struct_mumps.info[0] == -9) || (struct_mumps.info[0] == -8)
            || (struct_mumps.info[0] == -17) || (struct_mumps.info[0] == -20))
           && (new_percentage < coef_max_overestimate*100.0))
      {
        new_percentage = int(new_percentage*coef_increase_memory);
        struct_mumps.icntl[13] = new_percentage;
        struct_mumps.job = 2;
        CallMumps();
      }
    
    struct_mumps.icntl[13] = init_percentage;
    info_facto = struct_mumps.info[0];
  }
  
  
  //! Calls initialization routine provided by Mumps
  template<class T>
  void MatrixMumps<T>
  ::InitMatrix(bool sym, bool distributed)
  {
    // we clear previous factorization
    Clear();

#ifdef SELDON_WITH_MPI
    // MPI initialization for parallel version
    if (distributed)
      {
        // for distributed matrices
        // we consider that the communicator has been set previously
      }
    else
      {
	// centralized matrix => a linear system per processor
        struct_mumps.comm_fortran = MPI_Comm_c2f(MPI_COMM_SELF);
      }
#endif

    // symmetry is specified during the initialization stage
    struct_mumps.job = -1;
    if (sym)
      {
        struct_mumps.sym = 2; // general symmetric matrix
        threshold_pivot = 0.01;
      }
    else
      struct_mumps.sym = 0; // unsymmetric matrix

    // mumps is called
    CallMumps();
    
    // pivot threshold
    //struct_mumps.icntl[23] = 1;
    struct_mumps.cntl[0] = threshold_pivot;
    struct_mumps.icntl[13] = int(100.0*(coef_overestimate-1.0));
    //struct_mumps.cntl[2] = 1e-30;
    if (parallel_ordering)
      {
        struct_mumps.icntl[6] = 0;
        struct_mumps.icntl[27] = 2;
        struct_mumps.icntl[28] = type_ordering;
      }
    else
      {
        struct_mumps.icntl[6] = type_ordering;
        if (type_ordering == 1)
          struct_mumps.perm_in = perm.GetData();
      }
    
    // setting out of core parameters
    if (out_of_core)
      struct_mumps.icntl[21] = 1;
    else
      struct_mumps.icntl[21] = 0;

    struct_mumps.icntl[17] = 0;

    // the print level is set in mumps
    if (print_level >= 0)
      {
	struct_mumps.icntl[0] = 6;
	struct_mumps.icntl[1] = 0;
	struct_mumps.icntl[2] = 6;
	struct_mumps.icntl[3] = 2;
      }
    else
      {
	struct_mumps.icntl[0] = -1;
	struct_mumps.icntl[1] = -1;
	struct_mumps.icntl[2] = -1;
	struct_mumps.icntl[3] = 0;
      }
  }


  //! Selects another ordering scheme.
  template<class T>
  void MatrixMumps<T>::SelectOrdering(int num_ordering)
  {
    type_ordering = num_ordering;
  }


  //! Selects another ordering scheme (for distributed matrices)
  template<class T>
  void MatrixMumps<T>::SelectParallelOrdering(int num_ordering)
  {
    parallel_ordering = true;
    type_ordering = num_ordering;
  }


  //! Provides the permutation array
  template<class T>
  void MatrixMumps<T>::SetPermutation(const IVect& permut)
  {
    type_ordering = 1;
    perm.Reallocate(permut.GetM());
    for (int i = 0; i < perm.GetM(); i++)
      perm(i) = permut(i) + 1;
  }

  
  template<class T>
  void MatrixMumps<T>::SetPivotThreshold(double eps)
  {
    threshold_pivot = eps;
  }
  
  
  //! Destructor
  template<class T>
  MatrixMumps<T>::~MatrixMumps()
  {
    Clear();
  }


  //! Clears factorization
  template<class T>
  void MatrixMumps<T>::Clear()
  {
    if (struct_mumps.n > 0)
      {
	struct_mumps.job = -2;
	// Mumps variables are deleted.
        CallMumps();

	num_row_glob.Clear();
	num_col_glob.Clear();
	struct_mumps.n = 0;
      }
  }


  //! Informs Mumps that no message should be displayed
  template<class T>
  void MatrixMumps<T>::HideMessages()
  {
    print_level = -1;

    struct_mumps.icntl[0] = -1;
    struct_mumps.icntl[1] = -1;
    struct_mumps.icntl[2] = -1;
    struct_mumps.icntl[3] = 0;

  }


  //! Informs Mumps to display standard output
  template<class T>
  void MatrixMumps<T>::ShowMessages()
  {
    print_level = 0;

    struct_mumps.icntl[0] = 6;
    struct_mumps.icntl[1] = 0;
    struct_mumps.icntl[2] = 6;
    struct_mumps.icntl[3] = 2;

  }


  //! Enables writing on the disk (out of core).
  template<class T>
  void MatrixMumps<T>::EnableOutOfCore()
  {
    out_of_core = true;
  }


  //! Disables writing on the disk (incore).
  template<class T>
  void MatrixMumps<T>::DisableOutOfCore()
  {
    out_of_core = false;
  }

  
  //! Sets the coefficient used to overestimate the needed memory
  template<class T>
  void MatrixMumps<T>::SetCoefficientEstimationNeededMemory(double coef)
  {
    coef_overestimate = coef;
  }
  
  
  //! Sets the maximal allowed coefficient for the memory space multiplication
  template<class T>
  void MatrixMumps<T>::SetMaximumCoefficientEstimationNeededMemory(double coef)
  {
    coef_max_overestimate = coef;
  }
  
  
  //! Sets multiplication factor for each try to factorize the matrix
  template<class T>
  void MatrixMumps<T>::SetIncreaseCoefficientEstimationNeededMemory(double coef)
  {
    coef_increase_memory = coef;
  }
  

  //! Computes an ordering for matrix renumbering.
  /*!
    \param[in,out] mat matrix whose we want to find the ordering
    \param[out] numbers new row numbers
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class T0, class Prop,class Storage,class Allocator>
  void MatrixMumps<T>::FindOrdering(Matrix<T0, Prop, Storage, Allocator> & mat,
				    IVect& numbers, bool keep_matrix)
  {
    InitMatrix(IsSymmetricMatrix(mat));

    int n = mat.GetM();
    // conversion in coordinate format
    Vector<MUMPS_INT> num_row, num_col; Vector<T> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);
    
    long nnz = num_col.GetSize();
    // no values needed to renumber
    values.Clear();
    if (!keep_matrix)
      mat.Clear();

    struct_mumps.n = n; struct_mumps.nz = nnz;
    struct_mumps.nnz = nnz;
    struct_mumps.irn = num_row.GetData();
    struct_mumps.jcn = num_col.GetData();

    // Call the MUMPS package.
    struct_mumps.job = 1; // we analyse the system
    CallMumps();

    info_facto = struct_mumps.info[0];

    numbers.Reallocate(n);
    for (int i = 0; i < n; i++)
      numbers(i) = struct_mumps.sym_perm[i]-1;
  }

  
  //! Factorizes a given matrix
  /*!
    \param[in,out] mat matrix to factorize
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class T0, class Prop, class Storage, class Allocator>
  void MatrixMumps<T>::FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
				       bool keep_matrix)
  {
    int n = mat.GetM();
    bool sym = IsSymmetricMatrix(mat);
    // conversion in coordinate format with fortran convention (1-index)
    Vector<MUMPS_INT> num_row, num_col; Vector<T> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);
    if (!keep_matrix)
      mat.Clear();
    
    FactorizeCoordinate1(n, num_row, num_col, values, sym);
  }


  //! Factorizes a coordinate matrix with 1-index numbering
  template<class T>
  void MatrixMumps<T>::FactorizeCoordinate1(int n, Vector<MUMPS_INT>& num_row,
					    Vector<MUMPS_INT>& num_col,
					    Vector<T>& values, bool sym)
  {
    InitMatrix(sym);
    
    long nnz = values.GetM();
    //if (nnz > INT_MAX)
    //cout << "Factorisation d'une grosse matrice" << endl;

    struct_mumps.n = n; struct_mumps.nz = nnz;
    struct_mumps.nnz = nnz;
    struct_mumps.irn = num_row.GetData();
    struct_mumps.jcn = num_col.GetData();
    struct_mumps.a = reinterpret_cast<pointer>(values.GetData());

    // Call the MUMPS package.
    struct_mumps.job = 4; // we analyse and factorize the system
    CallMumps();
    
    IterateFacto();
  }


  //! Symbolic factorization
  template<class T> template<class Prop, class Storage, class Allocator>
  void MatrixMumps<T>
  ::PerformAnalysis(Matrix<T, Prop, Storage, Allocator> & mat)
  {
    InitMatrix(IsSymmetricMatrix(mat));

    int n = mat.GetM();
    long nnz = mat.GetNonZeros();
    // conversion in coordinate format with fortran convention (1-index)
    Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row_glob, num_col_glob, values, 1);

    struct_mumps.n = n; struct_mumps.nz = nnz;
    struct_mumps.nnz = nnz;
    struct_mumps.irn = num_row_glob.GetData();
    struct_mumps.jcn = num_col_glob.GetData();
    struct_mumps.a = reinterpret_cast<pointer>(values.GetData());

    // Call the MUMPS package.
    struct_mumps.job = 1; // we analyse the system
    CallMumps();

    info_facto = struct_mumps.info[0];
  }


  //! Numerical factorization
  /*!
    Be careful, because no conversion is performed in the method,
    so you have to choose RowSparse/ColSparse for unsymmetric matrices
    and RowSymSparse/ColSymSparse for symmetric matrices.
    The other formats should not work
  */
  template<class T> template<class Prop, class Allocator>
  void MatrixMumps<T>
  ::PerformFactorization(Matrix<T, Prop, RowSparse, Allocator> & mat)
  {
    // we consider that the values are corresponding
    // to the row/column numbers given for the analysis
    struct_mumps.a = reinterpret_cast<pointer>(mat.GetData());
    
    // Call the MUMPS package.
    struct_mumps.job = 2; // we factorize the system
    CallMumps();
    
    IterateFacto();
  }


    //! Numerical factorization
  /*!
    Be careful, because no conversion is performed in the method,
    so you have to choose RowSparse/ColSparse for unsymmetric matrices
    and RowSymSparse/ColSymSparse for symmetric matrices.
    The other formats should not work
  */
  template<class T> template<class Prop, class Allocator>
  void MatrixMumps<T>
  ::PerformFactorization(Matrix<T, Prop, RowSymSparse, Allocator> & mat)
  {
    // we consider that the values are corresponding
    // to the row/column numbers given for the analysis
    struct_mumps.a = reinterpret_cast<pointer>(mat.GetData());
    
    // Call the MUMPS package.
    struct_mumps.job = 2; // we factorize the system
    CallMumps();
    
    IterateFacto();
  }


  //! Numerical factorization
  /*!
    Be careful, because no conversion is performed in the method,
    so you have to choose RowSparse/ColSparse for unsymmetric matrices
    and RowSymSparse/ColSymSparse for symmetric matrices.
    The other formats should not work
  */
  template<class T> template<class Prop, class Storage, class Allocator>
  void MatrixMumps<T>
  ::PerformFactorization(Matrix<T, Prop, Storage, Allocator> & mat)
  {
    Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row_glob, num_col_glob, values, 1);
    
    // we consider that the values are corresponding
    // to the row/column numbers given for the analysis
    struct_mumps.a = reinterpret_cast<pointer>(values.GetData());
    
    // Call the MUMPS package.
    struct_mumps.job = 2; // we factorize the system
    CallMumps();
    
    IterateFacto();
  }


  //! Returns memory used by the factorisation in bytes
  template<class T>
  size_t MatrixMumps<T>::GetMemorySize() const
  {
    size_t taille = sizeof(int)*(num_row_glob.GetM()
				 + num_col_glob.GetM() + perm.GetM());
   
    if (struct_mumps.n <= 0)
      return taille;
    
    size_t nnz = struct_mumps.info[8];
    if (struct_mumps.info[8] < 0)
      nnz = abs(struct_mumps.info[8])*size_t(1024*1024);

    taille += sizeof(T)*nnz ;
    nnz = struct_mumps.info[9];
    if (struct_mumps.info[9] < 0)
      nnz = abs(struct_mumps.info[9])*size_t(1024*1024);
    
    taille += sizeof(int)*nnz;
    return taille;
  }
  

  //! Returns information about factorization performed
  template<class T>
  int MatrixMumps<T>::GetInfoFactorization() const
  {
    return info_facto;
  }


  //! Computation of Schur complement.
  /*!
    \param[in,out] mat initial matrix.
    \param[in] num numbers to keep in Schur complement.
    \param[out] mat_schur Schur matrix.
    \param[in] keep_matrix if false, \a mat is cleared.
  */
  template<class T> template<class Prop1, class Storage1, class Allocator,
			     class Prop2, class Storage2, class Allocator2>
  void MatrixMumps<T>::
  GetSchurMatrix(Matrix<T, Prop1, Storage1, Allocator>& mat, const IVect& num,
		 Matrix<T, Prop2, Storage2, Allocator2> & mat_schur,
		 bool keep_matrix)
  {
    InitMatrix(IsSymmetricMatrix(mat));
    
    int n_schur = num.GetM();
    // Subscripts are changed to respect fortran convention
    Vector<MUMPS_INT> index_schur(n_schur);
    for (int i = 0; i < n_schur; i++)
      index_schur(i) = num(i)+1;

    // array that will contain values of Schur matrix
    Vector<T, VectFull, Allocator2> vec_schur(n_schur*n_schur);
    vec_schur.Fill(0);
    
    struct_mumps.icntl[18] = 1;
    struct_mumps.size_schur = n_schur;
    struct_mumps.listvar_schur = index_schur.GetData();
    struct_mumps.schur = reinterpret_cast<pointer>(vec_schur.GetData());

    int n = mat.GetM(); long nnz = mat.GetNonZeros();
    // conversion in coordinate format with fortran convention (1-index)
    Vector<MUMPS_INT> num_row, num_col;
    Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);
    if (!keep_matrix)
      mat.Clear();

    struct_mumps.n = n;
    struct_mumps.nz = nnz;     struct_mumps.nnz = nnz;
    struct_mumps.irn = num_row.GetData();
    struct_mumps.jcn = num_col.GetData();
    struct_mumps.a = reinterpret_cast<pointer>(values.GetData());

    // Call the MUMPS package.
    struct_mumps.job = 4; // we analyse and factorize the system
    CallMumps();
    
    IterateFacto();
    
    // resetting parameters related to Schur complement
    struct_mumps.icntl[18] = 0;
    struct_mumps.size_schur = 0;
    struct_mumps.listvar_schur = NULL;
    struct_mumps.schur = NULL;

    // schur complement stored by rows
    int nb = 0;
    mat_schur.Reallocate(n_schur, n_schur);
    if (IsSymmetricMatrix(mat))
      for (int i = 0; i < n_schur; i++)
	for (int j = 0; j <= i; j++)
	  {
	    mat_schur(i, j) = vec_schur(i*n_schur + j);
	    mat_schur(j, i) = vec_schur(i*n_schur + j);
	  }
    else
      for (int i = 0; i < n_schur; i++)
	for (int j = 0; j < n_schur;j++)
	  mat_schur(i, j) = vec_schur(nb++);
    
    vec_schur.Clear(); index_schur.Clear();
  }


  //! Solves a linear system using the computed factorization
  /*!
    \param[in] TransA solves A x = b or A^T x = b
    \param[in,out] x right-hand-side on input, solution on output
    It is assumed that a call to FactorizeMatrix has been done before
  */
  template<class T> template<class Allocator2>
  void MatrixMumps<T>::Solve(const SeldonTranspose& TransA,
			     Vector<T, VectFull, Allocator2>& x)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (x.GetM() != struct_mumps.n)
      throw WrongDim("Mumps::Solve(TransA, c)",
                     string("The length of x is equal to ")
                     + to_str(x.GetM())
                     + " while the size of the matrix is equal to "
                     + to_str(struct_mumps.n) + ".");
#endif
    
    Solve(TransA, x.GetData(), 1);
  }


  //! Solves linear system with multiple right hand sides
  template<class T>
  void MatrixMumps<T>::Solve(const SeldonTranspose& TransA, T* x_ptr, int nrhs)
  {
    if (TransA.Trans())
      struct_mumps.icntl[8] = 0;
    else
      struct_mumps.icntl[8] = 1;
    
    struct_mumps.nrhs = nrhs;
    struct_mumps.lrhs = struct_mumps.n;
    struct_mumps.rhs = reinterpret_cast<pointer>(x_ptr);
    struct_mumps.job = 3; // we solve system
    CallMumps();
  }


  //! Solves a linear system using the factorization computed
  template<class T> template<class Allocator2>
  void MatrixMumps<T>::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves a linear system using a computed factorization.
  /*!    
    \param[in,out] x on entry, the right-hand-side; on exit, the solution.
    It is assumed that 'FactorizeMatrix' has already been called.
  */
  template<class T>
  template<class Allocator2, class Prop>
  void MatrixMumps<T>::Solve(const SeldonTranspose& TransA,
			     Matrix<T, Prop, ColMajor, Allocator2>& x)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (x.GetM() != struct_mumps.n)
      throw WrongDim("Mumps::Solve", string("Row size of x is equal to ")
                     + to_str(x.GetM()) + " while size of matrix is equal to "
                     + to_str(struct_mumps.n));
#endif
    
    Solve(TransA, x.GetData(), x.GetN());
  }


#ifdef SELDON_WITH_MPI
  template<class T>
  void MatrixMumps<T>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
                             Vector<int>& IndRow, Vector<T>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    FactorizeParallel(comm_facto, Ptr, IndRow, Val, glob_number,
                      sym, keep_matrix);
  }
  
  
  template<class T>
  void MatrixMumps<T>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
                             Vector<int64_t>& IndRow, Vector<T>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    cout << "Mumps with 64-bits integers does not work in parallel" << endl;
    abort();
  }


  //! factorization of a given matrix in distributed form (parallel execution)
  /*!
    \param[in,out] comm_facto MPI communicator
    \param[in,out] Ptr start indices for each column
    \param[in,out] IndRow row indices
    \param[in,out] Val data
    \param[in] sym if true, the matrix is assumed to be symmetric (upper part
    is provided)
    \param[in] glob_number row numbers (in the global matrix)
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T>
  void MatrixMumps<T>::
  FactorizeParallel(MPI_Comm& comm_facto, Vector<long>& Ptr,
                    Vector<MUMPS_INT>& IndRow, Vector<T>& Val,
                    const Vector<int>& glob_number,
                    bool sym, bool keep_matrix)
  {
    // for distributed matrix, we convert C++ communicator
    struct_mumps.comm_fortran = MPI_Comm_c2f(comm_facto);
    
    // Initialization depending on symmetry of the matrix.
    InitMatrix(sym, true);
    
    // Fortran communicator
    struct_mumps.comm_fortran = MPI_Comm_c2f(comm_facto);

    // distributed matrix
    struct_mumps.icntl[17] = 3;
    
    // finding the size of the overall system
    int nmax = 0, N = 0;
    for (int i = 0; i < glob_number.GetM(); i++)
      nmax = max(glob_number(i)+1, nmax);

    MPI_Allreduce(&nmax, &N, 1, MPI_INTEGER, MPI_MAX, comm_facto);

    // number of non-zero entries on this processor
    long nnz = IndRow.GetSize();

    // adding 1 to have 1-based indexes
    Vector<MUMPS_INT> IndCol(nnz);
    for (long i = 0; i < IndRow.GetM(); i++)
      IndRow(i)++;

    for (int i = 0; i < Ptr.GetM()-1; i++)
      for (long j = Ptr(i); j < Ptr(i+1); j++)
        IndCol(j) = glob_number(i) + 1;
    
    if (!keep_matrix)
      Ptr.Clear();
    
    // Define the problem on the host
    struct_mumps.n = N;
    struct_mumps.nz_loc = nnz; struct_mumps.nnz_loc = nnz;
    struct_mumps.irn_loc = IndRow.GetData();
    struct_mumps.jcn_loc = IndCol.GetData();
    struct_mumps.a_loc = reinterpret_cast<pointer>(Val.GetData());
    
    // Call the MUMPS package.
    struct_mumps.job = 1; // we analyse the system
    CallMumps();
    
    // overestimating size in order to avoid error -9
    double coef = coef_overestimate;
    struct_mumps.icntl[22] = int(coef*struct_mumps.infog[25]);
    struct_mumps.job = 2; // we factorize the system
    CallMumps();
    
    int info = 0;
    MPI_Allreduce(&struct_mumps.info[0], &info, 1, MPI_INTEGER, MPI_MIN, comm_facto);
    
    bool test_loop = true;
    int init_percentage = struct_mumps.icntl[13];
    int new_percentage = init_percentage;
    // loop in order to complete factorisation
    while (test_loop)
      {
        if (((info == -9)||(info==-8)||(info==-17)||(info==-19)|| (info==-20))
             && (coef < coef_max_overestimate))
          {
            coef *= coef_increase_memory;
            // increasing icntl(23) if error -9 occured
            struct_mumps.icntl[22] = int(coef*struct_mumps.infog[25]);
            new_percentage = int(new_percentage*coef_increase_memory);
            struct_mumps.icntl[13] = new_percentage;        
          }
        else
          test_loop = false;
        
        if (test_loop)
          {
            CallMumps();
            
            MPI_Allreduce(&struct_mumps.info[0],
                          &info, 1, MPI_INTEGER, MPI_MIN, comm_facto);
          }
      }

    struct_mumps.icntl[13] = init_percentage;
    info_facto = info;
    
    int rank_proc; MPI_Comm_rank(comm_facto, &rank_proc);
    if ((rank_proc == 0) && (print_level >= 0))
      cout<<"Factorization completed"<<endl;

  }


  template<class T>
  void MatrixMumps<T>
  ::PerformAnalysisDistributed(MPI_Comm& comm_facto, Vector<long>& Ptr,
			       Vector<int>& IndRow, Vector<T>& Val,
			       const Vector<int>& glob_number,
			       bool sym, bool keep_matrix)
  {
    // for distributed matrix, we convert C++ communicator
    struct_mumps.comm_fortran = MPI_Comm_c2f(comm_facto);

    // Initialization depending on symmetry of the matrix.
    InitMatrix(sym, true);
    
    // Fortran communicator
    struct_mumps.comm_fortran = MPI_Comm_c2f(comm_facto);

    // distributed matrix
    struct_mumps.icntl[17] = 3;
    
    // finding the size of the overall system
    int nmax = 0, N = 0;
    for (int i = 0; i < glob_number.GetM(); i++)
      nmax = max(glob_number(i)+1, nmax);

    MPI_Allreduce(&nmax, &N, 1, MPI_INTEGER, MPI_MAX, comm_facto);

    // number of non-zero entries on this processor
    long nnz = IndRow.GetSize();

    // adding 1 to have 1-based indexes
    Vector<int> IndCol(nnz);
    for (long i = 0; i < IndRow.GetM(); i++)
      IndRow(i)++;

    for (int i = 0; i < Ptr.GetM()-1; i++)
      for (long j = Ptr(i); j < Ptr(i+1); j++)
        IndCol(j) = glob_number(i) + 1;
    
    if (!keep_matrix)
      Ptr.Clear();

    // Define the problem on the host
    struct_mumps.n = N;
    struct_mumps.nz_loc = nnz;
    struct_mumps.nnz_loc = nnz;
    struct_mumps.irn_loc = IndRow.GetData();
    struct_mumps.jcn_loc = IndCol.GetData();
    struct_mumps.a_loc = reinterpret_cast<pointer>(Val.GetData());
    
    // Call the MUMPS package.
    struct_mumps.job = 1; // we analyse the system
    CallMumps();
  }
  
  
  template<class T>
  void MatrixMumps<T>
  ::PerformAnalysisDistributed(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
			       Vector<int64_t>& IndRow, Vector<T>& Val,
			       const Vector<int>& glob_number,
			       bool sym, bool keep_matrix)
  {
    cout << "Mumps with 64-bits integers does not work in parallel" << endl;
    abort();
  }
  
  template<class T>
  void MatrixMumps<T>
  ::PerformFactorizationDistributed(MPI_Comm& comm_facto, Vector<long>& Ptr,
				    Vector<int>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix)
  {
    // number of non-zero entries on this processor
    long nnz = IndRow.GetSize();

    // adding 1 to have 1-based indexes
    Vector<int> IndCol(nnz);
    for (long i = 0; i < IndRow.GetM(); i++)
      IndRow(i)++;

    for (int i = 0; i < Ptr.GetM()-1; i++)
      for (long j = Ptr(i); j < Ptr(i+1); j++)
        IndCol(j) = glob_number(i) + 1;
    
    if (!keep_matrix)
      Ptr.Clear();

    struct_mumps.nz_loc = nnz;
    struct_mumps.nnz_loc = nnz;
    struct_mumps.irn_loc = IndRow.GetData();
    struct_mumps.jcn_loc = IndCol.GetData();
    struct_mumps.a_loc = reinterpret_cast<pointer>(Val.GetData());

    // overestimating size in order to avoid error -9
    double coef = coef_overestimate;
    struct_mumps.icntl[22] = int(coef*struct_mumps.infog[25]);
    struct_mumps.job = 2; // we factorize the system

    CallMumps();
    
    int info = 0;
    MPI_Allreduce(&struct_mumps.info[0], &info, 1, MPI_INTEGER, MPI_MIN, comm_facto);
    
    bool test_loop = true;
    int init_percentage = struct_mumps.icntl[13];
    int new_percentage = init_percentage;
    // loop in order to complete factorisation
    while (test_loop)
      {
        if (((info == -9)||(info==-8)||(info==-17)||(info==-20))
             && (coef < coef_max_overestimate))
          {
            coef *= coef_increase_memory;
            // increasing icntl(23) if error -9 occured
            struct_mumps.icntl[22] = int(coef*struct_mumps.infog[25]);
            new_percentage = int(new_percentage*coef_increase_memory);
            struct_mumps.icntl[13] = new_percentage;        
          }
        else
          test_loop = false;
        
        if (test_loop)
          {
            CallMumps();
            
            MPI_Allreduce(&struct_mumps.info[0],
                          &info, 1, MPI_INTEGER, MPI_MIN, comm_facto);
          }
      }

    struct_mumps.icntl[13] = init_percentage;
    info_facto = info;
    
    int rank_proc; MPI_Comm_rank(comm_facto, &rank_proc);
    if ((rank_proc == 0) && (print_level >= 0))
      cout<<"Factorization completed"<<endl;
  }
  
  
  template<class T>
  void MatrixMumps<T>
  ::PerformFactorizationDistributed(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
				    Vector<int64_t>& IndRow, Vector<T>& Val,
				    const Vector<int>& glob_number,
				    bool sym, bool keep_matrix)
  {
    cout << "Mumps with 64-bits integers does not work in parallel" << endl;
    abort();
  }
  
  
  //! solves linear system with parallel execution
  /*!
    \param[in] TransA we solve A x = b or A^T x = b
    \param[in,out] x right-hand-side then solution
    \param[in,out] glob_num global row numbers
  */
  template<class T> template<class Allocator2>
  void MatrixMumps<T>::SolveDistributed(MPI_Comm& comm_facto,
                                        const SeldonTranspose& TransA,
					Vector<T, VectFull, Allocator2>& x,
					const IVect& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), 1, glob_num);
  }

  
  //! solves linear system with parallel execution
  /*!
    \param[in] TransA we solve A x = b or A^T x = b
    \param[in,out] x right-hand-side then solution
    \param[in,out] glob_num global row numbers
  */
  template<class T> template<class Allocator2, class Prop>
  void MatrixMumps<T>::SolveDistributed(MPI_Comm& comm_facto,
                                        const SeldonTranspose& TransA,
					Matrix<T, Prop, ColMajor, Allocator2>& x,
					const IVect& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), x.GetN(), glob_num);
  }

  
  //! solves linear system with parallel execution
  template<class T>
  void MatrixMumps<T>::SolveDistributed(MPI_Comm& comm_facto,
                                        const SeldonTranspose& TransA,
					T* x_ptr, int nrhs,
					const IVect& glob_num)
  {
    // storing the global right hand side
    Matrix<T, General, ColMajor> rhs;
    int cplx = sizeof(T)/8;
    int nblock = 8; // 8 right hand sides at the same time
    // allocating the global right hand side
    int rank_proc; MPI_Comm_rank(comm_facto, &rank_proc);
    if (rank_proc == 0)
      {
        rhs.Reallocate(struct_mumps.n, min(nrhs, nblock));
        rhs.Zero();
      }
    
    Matrix<T, General, ColMajor> xp;
    int nb_procs; MPI_Comm_size(comm_facto, &nb_procs);
    MPI_Status status;
 
    // retrieving dofs of all processors
    Vector<IVect> nump(nb_procs);
    if (nb_procs > 1)
      {
        if (rank_proc == 0)
          {
            for (int i = 0; i < nb_procs; i++)
              {
                if (i != 0)
                  {
                    int nb_dof;
                    MPI_Recv(&nb_dof, 1, MPI_INTEGER, i, 33, comm_facto, &status);
                    
                    nump(i).Reallocate(nb_dof);
                    MPI_Recv(nump(i).GetData(), nb_dof, MPI_INTEGER,
                             i, 38, comm_facto, &status);
                  }
                else
                  nump(i) = glob_num;
                
              }
          }
        else
          {
            int nb = glob_num.GetM();
            MPI_Send(&nb, 1, MPI_INTEGER, 0, 33, comm_facto);
            MPI_Send(glob_num.GetData(), nb, MPI_INTEGER, 0, 38, comm_facto);
          }
      }
    
    // then loop over blocks of right hand sides
    int num_rhs = 0; long ldb = glob_num.GetM();
    int lvl = print_level; HideMessages();
    while (num_rhs < nrhs)
      {
        int nrhs_p = min(nrhs-num_rhs, nblock);
        if (rank_proc == 0)
          {
            // on the host, we retrieve datas of all the other processors
            if (nb_procs > 1)
              {
                // assembling the right hand side
                for (int i = 0; i < nb_procs; i++)
                  {
                    
                    if (i != 0)
                      {
                        // On the host processor receiving components of right
                        // hand side.
                        int nb_dof = nump(i).GetM();
                        xp.Reallocate(nb_dof, nrhs_p);
                        
                        if (nb_dof > 0)
                          MPI_Recv(xp.GetDataVoid(), cplx*nb_dof*nrhs_p,
                                   MPI_DOUBLE, i, 37, comm_facto, &status);
                      }
                    else
                      {
                        xp.Reallocate(glob_num.GetM(), nrhs_p);
                        for (int k = 0; k < nrhs_p; k++)
                          for (int j = 0; j < glob_num.GetM(); j++)
                            xp(j, k) = x_ptr[j + (num_rhs+k)*ldb];
                      }
                    
                    for (int k = 0; k < nrhs_p; k++)
                      for (int j = 0; j < nump(i).GetM(); j++)
                        rhs(nump(i)(j), k) = xp(j, k);
                  }
              }
            else
              {
                for (int k = 0; k < nrhs_p; k++)
                  for (int j = 0; j < glob_num.GetM(); j++)
                    rhs(j, k) = x_ptr[j + (num_rhs+k)*ldb];
              }
          
            struct_mumps.rhs = reinterpret_cast<pointer>(rhs.GetData());
          }
        else
          {
            // On other processors, we send right hand side.
            int nb_dof = glob_num.GetM();
            if (nb_dof > 0)              
              MPI_Send(&x_ptr[num_rhs*nb_dof], cplx*nb_dof*nrhs_p, MPI_DOUBLE, 0, 37, comm_facto);
          }
	
        // we solve system
        if (TransA.Trans())
          struct_mumps.icntl[8] = 0;
        else
          struct_mumps.icntl[8] = 1;
        
        struct_mumps.nrhs = nrhs_p;
        struct_mumps.lrhs = struct_mumps.n;
        struct_mumps.job = 3;
        CallMumps();
        
        // we distribute solution on all the processors
        if (nb_procs > 1)
          {
            if (rank_proc == 0)
              {
                for (int i = 0; i < nb_procs; i++)
                  {
                    if (i != 0)
                      {
                        int nb_dof = nump(i).GetM(); 
                        xp.Reallocate(nb_dof, nrhs_p);
                        for (int j = 0; j < nb_dof; j++)
                          for (int k = 0; k < nrhs_p; k++)
                            xp(j, k) = rhs(nump(i)(j), k);
                        
                        MPI_Send(xp.GetDataVoid(), cplx*nb_dof*nrhs_p,
                                  MPI_DOUBLE, i, 40, comm_facto);
                      }
                    else
                      {
                        for (int k = 0; k < nrhs_p; k++)
                          for (int j = 0; j < glob_num.GetM(); j++)
                            x_ptr[j + ldb*(num_rhs+k)] = rhs(glob_num(j), k);
                      }
                  }
              }
            else
              {
                int nb_dof = glob_num.GetM();
                xp.Reallocate(nb_dof, nrhs_p);
                MPI_Recv(xp.GetDataVoid(), cplx*nb_dof*nrhs_p,
                         MPI_DOUBLE, 0, 40, comm_facto, &status);
                
                for (int k = 0; k < nrhs_p; k++)
                  for (int j = 0; j < glob_num.GetM(); j++)
                    x_ptr[j + (num_rhs+k)*ldb] = xp(j, k);
              }
          }
        else
          {
            for (int k = 0; k < nrhs_p; k++)
              for (int i = 0; i < glob_num.GetM(); i++)
                x_ptr[i + ldb*(num_rhs+k)] = rhs(i, k);
          }
        
        num_rhs += nrhs_p;
      }
    
    if (lvl >= 0)
      ShowMessages();
    
    print_level = lvl;
  }
#endif


  //! Factorization of a matrix of same type T as for the Mumps object   
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixMumps<T>& mat_lu, bool keep_matrix, T& x)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  

  //! Factorization of a complex matrix with a real Mumps object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixMumps<T>& mat_lu, bool keep_matrix, complex<T>& x)
  {
    throw WrongArgument("GetLU(Matrix<complex<T> >& A, MatrixMumps<T>& mat_lu, bool)",
			"The LU matrix must be complex");
  }


  //! Factorization of a real matrix with a complex Mumps object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixMumps<complex<T> >& mat_lu, bool keep_matrix, T& x)
  {
    throw WrongArgument("GetLU(Matrix<T>& A, MatrixMumps<complex<T> >& mat_lu, bool)",
			"The sparse matrix must be complex");
  }
  

  //! Factorization of a general matrix with Mumps
  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixMumps<T>& mat_lu,
	     bool keep_matrix)
  {
    // we check if the type of non-zero entries of matrix A
    // and of the Mumps object (T) are different
    // we call one of the GetLUs written above
    // such a protection avoids to compile the factorisation of a complex
    // matrix with a real Mumps object
    typename Matrix<T0, Prop, Storage, Allocator>::entry_type x;
    GetLU(A, mat_lu, keep_matrix, x);
  }
  
  
  //! computes Schur complement of a symmetric matrix
  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T, Symmetric, Storage, Allocator>& A,
                      MatrixMumps<T>& mat_lu, const IVect& num,
                      MatrixFull& schur_matrix, bool keep_matrix)
  {
    mat_lu.GetSchurMatrix(A, num, schur_matrix, keep_matrix);
  }


  //! computes Schur complement of an unsymmetric matrix
  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T, General, Storage, Allocator>& A,
                      MatrixMumps<T>& mat_lu, const IVect& num,
                      MatrixFull& schur_matrix, bool keep_matrix)
  {
    mat_lu.GetSchurMatrix(A, num, schur_matrix, keep_matrix);
  }


  //! LU resolution with a vector whose type is the same as for Mumps object
  template<class T, class Allocator>
  void SolveLU(MatrixMumps<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  //! LU resolution with a vector whose type is the same as for Mumps object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! Solves A x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(MatrixMumps<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x)
  {
    Matrix<double, General, ColMajor> y(x.GetM(), 2);
    
    for (int i = 0; i < x.GetM(); i++)
      {
	y(i, 0) = real(x(i));
	y(i, 1) = imag(x(i));
      }
    
    SolveLU(mat_lu, y);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<double>(y(i, 0), y(i, 1));
  }
  

  //! Solves A x = b or A^T x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x)
  {
    Matrix<double, General, ColMajor> y(x.GetM(), 2);
    
    for (int i = 0; i < x.GetM(); i++)
      {
	y(i, 0) = real(x(i));
	y(i, 1) = imag(x(i));
      }
    
    SolveLU(TransA, mat_lu, y);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<double>(y(i, 0), y(i, 1));

  }


  //! Solves A x = b, where A is complex and x is real => Forbidden
  template<class Allocator>
  void SolveLU(MatrixMumps<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixMumps<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }

  
  //! Solves A x = b or A^T x = b, where A is complex and x is real => Forbidden  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixMumps<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }


  //! LU resolution with a matrix whose type is the same as for Mumps object
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixMumps<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(SeldonNoTrans, x);
  }


  //! LU resolution with a matrix whose type is the same as for Mumps object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Allocator, class Prop>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  template<class Prop, class Allocator>
  void SolveLU(MatrixMumps<double>& mat_lu,
	       Matrix<complex<double>, Prop, ColMajor, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixMumps<double>, Matrix<complex<double> >)", 
			"The result should be a real vector");
  }

  
  template<class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<double>& mat_lu,
	       Matrix<complex<double>, Prop, ColMajor, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixMumps<double>, Matrix<complex<double> >)", 
			"The result should be a real vector");
  }
  

  template<class Prop, class Allocator>
  void SolveLU(MatrixMumps<complex<double> >& mat_lu,
	       Matrix<double, Prop, ColMajor, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixMumps<complex<double> >, Matrix<double>)", 
			"The result should be a complex matrix");
  }

  
  template<class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixMumps<complex<double> >& mat_lu,
	       Matrix<double, Prop, ColMajor, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixMumps<complex<double> >, Matrix<double>)", 
			"The result should be a complex matrix");
  }

}

#define SELDON_FILE_MUMPS_CXX
#endif
