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


#ifndef SELDON_FILE_SUPERLU_CXX

#include "SuperLU.hxx"

namespace Seldon
{
#ifndef SELDON_WITH_SUPERLU_DIST      
  void SetComplexOne(superlu::doublecomplex& one)
  {
    one.r = 1.0;
    one.i = 0.0;
  }
  

  // The function comes from the Matlab interface to SuperLU. It is part of
  // SuperLU package. Its copyright is held by University of California
  // Berkeley, Xerox Palo Alto Research Center and Lawrence Berkeley National
  // Lab. It is released under a license compatible with the GNU LGPL.
  template<class T>
  void LUextract(superlu::SuperMatrix *L, superlu::SuperMatrix *U,
                 T *Lval, int *Lrow, long *Lcol,
                 T *Uval, int *Urow, long *Ucol, int *snnzL,
		 int *snnzU)
  {
    int         i, j, k;
    int         upper;
    int         fsupc, istart, nsupr;
    long         lastl = 0, lastu = 0;
    T      *SNptr;
    T one;
    
    SetComplexOne(one);

#ifdef SELDON_WITH_SUPERLU_MT
    superlu::SCPformat    *Lstore;
    superlu::NCPformat    *Ustore;
    Lstore = static_cast<superlu::SCPformat*>(L->Store);
    Ustore = static_cast<superlu::NCPformat*>(U->Store);
#else
    superlu::SCformat    *Lstore;
    superlu::NCformat    *Ustore;
    Lstore = static_cast<superlu::SCformat*>(L->Store);
    Ustore = static_cast<superlu::NCformat*>(U->Store);
#endif

    Lcol[0] = 0;
    Ucol[0] = 0;
    
    /* for each supernode */
    for (k = 0; k <= Lstore->nsuper; ++k) {
      
      fsupc = L_FST_SUPC(k);
      istart = L_SUB_START(fsupc);
      nsupr = L_SUB_START(fsupc+1) - istart;
      upper = 1;
      
      /* for each column in the supernode */
      for (j = fsupc; j < L_FST_SUPC(k+1); ++j) {
	SNptr = &(static_cast<T*>(Lstore->nzval))[L_NZ_START(j)];
	
	/* Extract U */
	for (i = U_NZ_START(j); i < U_NZ_START(j+1); ++i) {
	  Uval[lastu] = (static_cast<T*>(Ustore->nzval))[i];
	  Urow[lastu++] = U_SUB(i);
	}
	for (i = 0; i < upper; ++i) { /* upper triangle in the supernode */
	  Uval[lastu] = SNptr[i];
	  Urow[lastu++] = L_SUB(istart+i);
	}
	Ucol[j+1] = lastu;
	
	/* Extract L */
	Lval[lastl] = one; /* unit diagonal */
	Lrow[lastl++] = L_SUB(istart + upper - 1);
	for (i = upper; i < nsupr; ++i) {
	  Lval[lastl] = SNptr[i];
	  Lrow[lastl++] = L_SUB(istart+i);
	}
	Lcol[j+1] = lastl;
	
	++upper;
	
      } /* for j ... */
      
    } /* for k ... */
    
    *snnzL = lastl;
    *snnzU = lastu;
  }
#endif


  /**********************
   * MatrixSuperLU_Base *
   **********************/
  
  
  //! default constructor
  template<class T>
  MatrixSuperLU_Base<T>::MatrixSuperLU_Base()
  {
    n = 0;
    
#ifndef SELDON_WITH_SUPERLU_DIST    
    
    permc_spec = superlu::COLAMD;
    Lstore = NULL;
    Ustore = NULL;
    
#ifdef SELDON_WITH_SUPERLU_MT
    nprocs = 1;
    diag_pivot_thresh  = 0.01;
    usepr              = superlu::NO;
    drop_tol           = 0.0;
    
#else
    superlu::set_default_options(&options);
#endif

#else

    // distributed version
    superlu::set_default_options_dist(&options);
    
    options.ParSymbFact = superlu::YES;
    options.ColPerm = superlu::PARMETIS; 
    options.IterRefine = superlu::NOREFINE;
    
#endif

    ShowMessages();
    display_info = false;
    info_facto = 0;
  }


  //! destructor
  template<class T>
  MatrixSuperLU_Base<T>::~MatrixSuperLU_Base()
  {
    Clear();
  }


  //! same effect as a call to the destructor
  template<class T>
  void MatrixSuperLU_Base<T>::Clear()
  {
    if (n > 0)
      {
#ifndef SELDON_WITH_SUPERLU_DIST    
	// SuperLU objects are cleared
	superlu::Destroy_SuperNode_Matrix(&L);
	superlu::Destroy_CompCol_Matrix(&U);
	if (permc_spec != superlu::MY_PERMC)
	  {
	    perm_r.Clear();
	    perm_c.Clear();
	  }
        
        superlu::StatFree(&stat);

#else

        superlu::PStatFree(&stat);
        
        // permuted matrix A is cleared
        superlu::Destroy_CompRowLoc_Matrix_dist(&A);    

        superlu::ScalePermstructFree(&ScalePermstruct);
        superlu::Destroy_LU(n, &grid, &LUstruct);
        superlu::LUstructFree(&LUstruct);

        superlu::superlu_gridexit(&grid);        
#endif
        
	n = 0;
      }
  }


  //! inits SuperLU computation
  template<class T>
  void MatrixSuperLU_Base<T>::Init(int size, int_t& panel_size, int_t& relax)
  {
    n = size;

#ifndef SELDON_WITH_SUPERLU_DIST
    panel_size = superlu::sp_ienv(1);
    relax = superlu::sp_ienv(2);
    
#ifdef SELDON_WITH_SUPERLU_MT
    fact               = superlu::EQUILIBRATE;
    refact             = superlu::NO;    

    superlu::StatAlloc(n, nprocs, panel_size, relax, &stat);
    superlu::StatInit(n, nprocs, &stat);
#else
    superlu::StatInit(&stat);
#endif

#else
    // Initialize ScalePermstruct and LUstruct.
    superlu::ScalePermstructInit(n, n, &ScalePermstruct);
    superlu::LUstructInit(n, &LUstruct);
    
    // Initialize the statistics variables.
    superlu::PStatInit(&stat);

#endif

  }


  template<class T>
  void MatrixSuperLU_Base<T>::SetNumberOfThreadPerNode(int p)
  {
#ifdef SELDON_WITH_SUPERLU_MT
    nprocs = p;
#endif
  }


  //! no message from SuperLU
  template<class T>
  void MatrixSuperLU_Base<T>::HideMessages()
  {
    display_info = false;
  }


  //! allows messages from SuperLU
  template<class T>
  void MatrixSuperLU_Base<T>::ShowMessages()
  {
    display_info = true;
  }


  template<class T>
  bool MatrixSuperLU_Base<T>::UseInteger8() const  
  {
    if (sizeof(int_t) == 8)
      return true;
    
    return false;
  }

  
  //! returns status of factorisation
  template<class T>
  int MatrixSuperLU_Base<T>::GetInfoFactorization() const
  {
    return info_facto;
  }


  //! Returns the permutation of rows.
  /*! In order to retain the sparsity as much as possible, SuperLU permutes
    rows and columns before the factorization. This method returns the row
    permutation that was employed in the factorization. This method is
    obviously to be called after the factorization has been performed.
    \return The permutation of the rows.
  */
  template<class T>
  const Vector<int_t>& MatrixSuperLU_Base<T>::GetRowPermutation() const
  {
    return perm_r;
  }


  //! Returns the permutation of columns.
  /*! In order to retain the sparsity as much as possible, SuperLU permutes
    rows and columns before the factorization. This method returns the column
    permutation that was employed in the factorization. This method is
    obviously to be called after the factorization has been performed.
    \return The permutation of the columns.
  */
  template<class T>
  const Vector<int_t>& MatrixSuperLU_Base<T>::GetColPermutation() const
  {
    return perm_c;
  }


  template<class T>
  void MatrixSuperLU_Base<T>::SelectOrdering(int type)
  {
    permc_spec = (colperm_t) type;
  }


  template<class T>
  void MatrixSuperLU_Base<T>::SetPermutation(const IVect& permut)
  {
    permc_spec = superlu::MY_PERMC;
    perm_c.Reallocate(permut.GetM());
    for (int i = 0; i < permut.GetM(); i++)
      perm_c(i) = permut(i);
    
    perm_r.Reallocate(perm_c.GetM());
    perm_r.Fill();
  }
  
  
  /*************************
   * MatrixSuperLU<double> *
   *************************/


#ifndef SELDON_WITH_SUPERLU_DIST
  //! Returns the LU factorization.
  /*!
    \param[out] Lmat matrix L in the LU factorization.
    \param[out] Umat matrix U in the LU factorization.
    \param[in] permuted should the permuted matrices be provided? SuperLU
    permutes the rows and columns of the factorized matrix. If \a permuted is
    set to true, L and U are returned as SuperLU computed them, hence with
    permuted rows and columns. If \a permuted is set to false, the matrices L
    and U are "unpermuted" so that L times U is equal to the initial matrix.
  */
  template<class Prop, class Allocator>
  void MatrixSuperLU<double>
  ::GetLU(Matrix<double, Prop, ColSparse, Allocator>& Lmat,
          Matrix<double, Prop, ColSparse, Allocator>& Umat,
          bool permuted)
  {
#ifdef SELDON_WITH_SUPERLU_MT
    Lstore = static_cast<superlu::SCPformat*>(L.Store);
    Ustore = static_cast<superlu::NCPformat*>(U.Store);
#else
    Lstore = static_cast<superlu::SCformat*>(L.Store);
    Ustore = static_cast<superlu::NCformat*>(U.Store);
#endif

    int Lnnz = Lstore->nnz;
    int Unnz = Ustore->nnz;

    int m = U.nrow;
    int n = U.ncol;

    Vector<double, VectFull, Allocator> Lval(Lnnz);
    Vector<int> Lrow(Lnnz);
    Vector<long> Lcol(n + 1);

    Vector<double, VectFull, Allocator> Uval(Unnz);
    Vector<int> Urow(Unnz);
    Vector<long> Ucol(n + 1);

    int Lsnnz;
    int Usnnz;
    LUextract(&L, &U, Lval.GetData(), Lrow.GetData(), Lcol.GetData(),
              Uval.GetData(), Urow.GetData(), Ucol.GetData(), &Lsnnz, &Usnnz);

    Lmat.SetData(m, n, Lval, Lcol, Lrow);
    Umat.SetData(m, n, Uval, Ucol, Urow);

    if (!permuted)
      {
        Vector<int> row_perm_orig = perm_r;
        Vector<int> col_perm_orig = perm_c;

        Vector<int> row_perm(n);
        Vector<int> col_perm(n);
        row_perm.Fill();
        col_perm.Fill();

        Sort(row_perm_orig, row_perm);
        Sort(col_perm_orig, col_perm);

        ApplyInversePermutation(Lmat, row_perm, col_perm);
        ApplyInversePermutation(Umat, row_perm, col_perm);
      }
  }


  //! Returns the LU factorization.
  /*!
    \param[out] Lmat matrix L in the LU factorization.
    \param[out] Umat matrix U in the LU factorization.
    \param[in] permuted should the permuted matrices be provided? SuperLU
    permutes the rows and columns of the factorized matrix. If \a permuted is
    set to true, L and U are returned as SuperLU computed them, hence with
    permuted rows and columns. If \a permuted is set to false, the matrices L
    and U are "unpermuted" so that L times U is equal to the initial matrix.
    \note This method will first retrieve the L and U matrices in 'ColSparse'
    format and then convert them into 'RowSparse'.
  */
  template<class Prop, class Allocator>
  void MatrixSuperLU<double>
  ::GetLU(Matrix<double, Prop, RowSparse, Allocator>& Lmat,
          Matrix<double, Prop, RowSparse, Allocator>& Umat,
          bool permuted)
  {
    Lmat.Clear();
    Umat.Clear();

    Matrix<double, Prop, ColSparse, Allocator> Lmat_col;
    Matrix<double, Prop, ColSparse, Allocator> Umat_col;
    GetLU(Lmat_col, Umat_col, permuted);

    Copy(Lmat_col, Lmat);
    Lmat_col.Clear();
    Copy(Umat_col, Umat);
    Umat_col.Clear();
  }


  //! Returns the size of memory used by the current object
  size_t MatrixSuperLU<double>::GetMemorySize() const
  {
    size_t taille = sizeof(int)*(perm_r.GetM()+perm_c.GetM());
    if (this->n > 0)
      {
#ifdef SELDON_WITH_SUPERLU_MT
        superlu::superlu_memusage_t mem_usage;
        int_t panel_size = superlu::sp_ienv(1);
        superlu::superlu_dQuerySpace(nprocs, const_cast<superlu::SuperMatrix*>(&L),
                                     const_cast<superlu::SuperMatrix*>(&U),
                                     panel_size, &mem_usage);
#else
        superlu::mem_usage_t mem_usage;
        superlu::dQuerySpace(const_cast<superlu::SuperMatrix*>(&L),
                             const_cast<superlu::SuperMatrix*>(&U), &mem_usage);
#endif
        
        taille += mem_usage.total_needed;
      }
    
    return taille;
  }
  
  
  //! factorization of matrix in double precision using SuperLU
  template<class T0, class Prop, class Storage, class Allocator>
  void MatrixSuperLU<double>::
  FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
		  bool keep_matrix)
  {
    // clearing previous factorization
    Clear();

    // conversion in CSC format
    General prop;
    Vector<superlu_int_t> Ptr, IndRow;
    Vector<double> Val;

    ConvertToCSC(mat, prop, Ptr, IndRow, Val, false);
    if (!keep_matrix)
      mat.Clear();

    FactorizeCSC(Ptr, IndRow, Val, false);
  }

  
  //! factorization of matrix given in CSC form
  void MatrixSuperLU<double>
  ::FactorizeCSC(Vector<superlu_int_t>& Ptr, Vector<superlu_int_t>& IndRow,
		 Vector<double>& Val, bool sym)
  {  
    // initializing parameters
    int_t panel_size, relax;
    int_t lwork = 0;
    Init(Ptr.GetM()-1, panel_size, relax);

    superlu::SuperMatrix AA;
    int_t nnz = IndRow.GetM();
    superlu::dCreate_CompCol_Matrix(&AA, n, n, nnz, Val.GetData(),
                                    reinterpret_cast<int_t*>(IndRow.GetData()),
                                    reinterpret_cast<int_t*>(Ptr.GetData()),
                                    superlu::SLU_NC,
                                    superlu::SLU_D, superlu::SLU_GE);

    // we get renumbering vectors perm_r and perm_c
    options.ColPerm = permc_spec;    
    if (permc_spec != superlu::MY_PERMC)
      {
        perm_r.Reallocate(n);
        perm_c.Reallocate(n);
        perm_r.Fill();
        perm_c.Fill();
        
        superlu::get_perm_c(permc_spec, &AA, perm_c.GetData());        
      }    
    
#ifdef SELDON_WITH_SUPERLU_MT
    superlu::SuperMatrix AC;
    superlu::trans_t  trans = superlu::NOTRANS;
    
    superlu::pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
                          diag_pivot_thresh, usepr, drop_tol,
                          perm_c.GetData(), perm_r.GetData(),
                          NULL, lwork, &AA, &AC, &options, &stat);
    
    int_t info;
    superlu::pdgstrf(&options, &AC, perm_r.GetData(), &L, &U, &stat, &info);
    
    info_facto = info;
    superlu::pxgstrf_finalize(&options, &AC);
    
    if (info_facto == 0 && display_info)
      {
        superlu::PrintStat(&stat);
      }
    
#else
    superlu::SuperMatrix A;
    // original matrix AA is permuted to obtain matrix A
    Vector<int> etree(n);
    sp_preorder(&options, &AA, perm_c.GetData(), etree.GetData(), &A);
    
    // then calling factorisation on permuted matrix
    superlu::dgstrf(&options, &A, relax, panel_size, etree.GetData(),
                    NULL, lwork, perm_c.GetData(), perm_r.GetData(), &L, &U,
                    &Glu, &stat, &info_facto);
    
    if (info_facto == 0 && display_info)
      {
	superlu::mem_usage_t mem_usage;
	Lstore = (superlu::SCformat *) L.Store;
	Ustore = (superlu::NCformat *) U.Store;
	cout << "Number of nonzeros in factor L = " << Lstore->nnz << endl;
	cout << "Number of nonzeros in factor U = " << Ustore->nnz << endl;
	cout << "Number of nonzeros in L+U     = "
             << Lstore->nnz + Ustore->nnz << endl;
	superlu::dQuerySpace(&L, &U, &mem_usage);
	cout << "Memory used for factorization in MB: "
             << mem_usage.total_needed / (1024. * 1024.) << endl;
      }

    superlu::Destroy_CompCol_Permuted(&A);
#endif
    
    // clearing matrices
    superlu::Destroy_CompCol_Matrix(&AA);
    
    Ptr.Nullify(); IndRow.Nullify(); Val.Nullify();
  }


  //! Solves linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(Vector<double, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system A x = b or A^T x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(const SeldonTranspose& TransA,
                                    Vector<double, VectFull, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), 1);
  }

  
  //! Solves linear system A x = b or A^T x = b
  void MatrixSuperLU<double>::Solve(const SeldonTranspose& TransA,
				    double* x_ptr, int nrhs_)
  {
    superlu::trans_t trans;
    if (TransA.NoTrans())
      trans = superlu::NOTRANS;
    else
      trans = superlu::TRANS;
    
    int_t nb_rhs = nrhs_, info;
    // Putting right hand side on SuperLU structure.
    superlu::dCreate_Dense_Matrix(&B, n, nb_rhs, x_ptr, n,
                                  superlu::SLU_DN, superlu::SLU_D, superlu::SLU_GE);

#ifdef SELDON_WITH_SUPERLU_MT
    superlu::dgstrs(trans, &L, &U, perm_r.GetData(),
                    perm_c.GetData(), &B, &stat, &info);
#else
    superlu::dgstrs(trans, &L, &U, perm_c.GetData(),
                    perm_r.GetData(), &B, &stat, &info);
#endif
    
    superlu::Destroy_SuperMatrix_Store(&B);
  }


  //! Solves linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(Matrix<double, General, ColMajor, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system A x = b or A^T x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(const SeldonTranspose& TransA,
                                    Matrix<double, General, ColMajor, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), x.GetN());
  }
  

#else
  
  /**********************************
   * Distributed version for double *
   **********************************/


  //! Returns the size of memory used by the current object
  size_t MatrixSuperLU<double>::GetMemorySize() const
  {
    size_t size = 0;
    if (n > 0)
      {
        superlu::mem_usage_t mem_usage;
        superlu::dQuerySpace_dist(n, const_cast<superlu::LUstruct_t*>(&LUstruct),
                                  const_cast<superlu::gridinfo_t*>(&grid),
                                  const_cast<superlu::SuperLUStat_t*>(&stat), &mem_usage);
        
        size += mem_usage.total;
      }
    
    return size;
  }


  //! factorization of matrix in double precision using SuperLU
  template<class T0, class Prop, class Storage, class Allocator>
  void MatrixSuperLU<double>::
  FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
		  bool keep_matrix)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<double> Val; General prop;
    ConvertToCSC(mat, prop, Ptr, IndRow, Val, false);
    if (!keep_matrix)
      mat.Clear();
    
    FactorizeCSC(Ptr, IndRow, Val, IsSymmetricMatrix(mat));
  }

  
  void MatrixSuperLU<double>
  ::FactorizeCSC(Vector<long>& Ptr, Vector<int>& IndRow,
		 Vector<double>& Val, bool sym)
  {  
    Vector<int> glob_num(1);
    glob_num(0) = 0;
    MPI_Comm comm = MPI_COMM_SELF;
    FactorizeDistributedMatrix(comm, Ptr, IndRow, Val,
			       glob_num, sym, false);
  }


  //! Solves linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(Vector<double, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system A x = b or A^T x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(const SeldonTranspose& TransA,
                                    Vector<double, VectFull, Allocator2>& x)
  {
    Vector<int> glob_num(1);
    glob_num(0) = 0;
    MPI_Comm comm = MPI_COMM_SELF;
    SolveDistributed(comm, TransA, x, glob_num);
  }


  void MatrixSuperLU<double>::Solve(const SeldonTranspose& TransA,
				    double* x_ptr, int nrhs_)
  {
    Vector<int> glob_num(1);
    glob_num(0) = 0;
    MPI_Comm comm = MPI_COMM_SELF;
    SolveDistributed(comm, TransA, x_ptr, nrhs_, glob_num);
  }
  

  //! Solves linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(Matrix<double, General, ColMajor, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system A x = b or A^T x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(const SeldonTranspose& TransA,
                                    Matrix<double, General, ColMajor, Allocator2>& x)
  {
    Vector<int> glob_num(1);
    glob_num(0) = 0;
    MPI_Comm comm = MPI_COMM_SELF;
    SolveDistributed(comm, TransA, x, glob_num);
  }


  void MatrixSuperLU<double>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
                             Vector<int>& IndRow, Vector<double>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    Vector<int> PtrInt(Ptr.GetM());
    for (int i = 0; i < Ptr.GetM(); i++)
      PtrInt(i) = Ptr(i);

    FactorizeParallel(comm_facto, PtrInt, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }

  
  void MatrixSuperLU<double>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
                             Vector<int64_t>& IndRow, Vector<double>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    FactorizeParallel(comm_facto, Ptr, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }
  
  
  template<class Tint>
  void MatrixSuperLU<double>::
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<Tint>&, Vector<Tint>&,
                    Vector<double>&,
                    const Vector<int>& glob_num,
                    bool sym, bool keep_matrix)
  {
    cout << "Not available for this type of integers " << endl;
    cout << "Size of Tint = " << sizeof(Tint) << endl;
    abort();
  }


  void MatrixSuperLU<double>::
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<superlu_int_t>& Ptr, Vector<superlu_int_t>& Row,
                    Vector<double>& Val, const Vector<int>& glob_num,
                    bool sym, bool keep_matrix)
  {
#ifdef SELDON_WITH_SUPERLU_DOUBLE
    // previous factorization is cleared if present
    Clear();
    
    if (!sym)
      {
        cout << "Problem with option TRANS of SuperLU" << endl;
        abort();
      }
    
    // m_loc : local number of rows
    // N : global number of rows
    int m_loc_ = Ptr.GetM()-1;    
    int N = m_loc_; nloc = m_loc_;
    MPI_Allreduce(&m_loc_, &N, 1, MPI_INTEGER, MPI_SUM, comm_facto);
    int_t m_loc = m_loc_;
    
    // structures are initialized with N
    int_t panel_size, relax;
    Init(N, panel_size, relax);

    // 1-D grid
    MPI_Comm_size(comm_facto, &nprow);
    npcol = 1;
    
    // initialize the superlu process grid
    superlu::superlu_gridinit(comm_facto, nprow, npcol, &grid);
    
    // fills the superlu structure
    // global numbers are assumed to be consecutive, we provide the first row number
    int_t fst_row = glob_num(0);
    int_t nnz_loc = Row.GetM();
    //superlu::SuperMatrix A;
    superlu::
      dCreate_CompRowLoc_Matrix_dist(&A, n, n, nnz_loc, m_loc, fst_row,
                                     Val.GetData(), Row.GetData(), Ptr.GetData(), 
                                     superlu::SLU_NR_loc, superlu::SLU_D,
                                     superlu::SLU_GE);
    
    // completes factorization
    int_t nrhs = 0;
    options.Trans = superlu::NOTRANS;
    options.Fact = superlu::DOFACT;
    superlu::
      pdgssvx(&options, &A, &ScalePermstruct,
              NULL, m_loc, nrhs, &grid,
              &LUstruct, &SOLVEstruct, NULL, &stat, &info_facto);
    
    Ptr.Nullify(); Row.Nullify(); Val.Nullify();
#else
    cout << "Recompile Seldon with SELDON_WITH_SUPERLU_DOUBLE" << endl;
    cout << "Double and complex<double> cannot be handled simultaneously" << endl;
    abort();
#endif
  }

  
  template<class Allocator2>
  void MatrixSuperLU<double>::
  SolveDistributed(MPI_Comm& comm_facto,
                   const SeldonTranspose& TransA,
                   Vector<double, VectFull, Allocator2>& x,
		   const IVect& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), 1, glob_num);
  }

  
  template<class Allocator2>
  void MatrixSuperLU<double>::
  SolveDistributed(MPI_Comm& comm_facto,
                   const SeldonTranspose& TransA,
                   Matrix<double, General, ColMajor, Allocator2>& x,
		   const IVect& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), x.GetN(), glob_num);
  }


  void MatrixSuperLU<double>::SolveDistributed(MPI_Comm& comm_facto,
					       const SeldonTranspose& TransA,
					       double* x_ptr, int nrhs_,
					       const IVect& glob_num)
  {
#ifdef SELDON_WITH_SUPERLU_DOUBLE
    options.Fact = superlu::FACTORED;
    // inverting transpose because we have provided columns instead of rows
    if (TransA.NoTrans())
      options.Trans = superlu::TRANS;
    else
      options.Trans = superlu::NOTRANS;
       
    options.Trans = superlu::NOTRANS;
    Vector<double> berr(nrhs_);
    for (int k = 0; k < nrhs_; k++)
      {
	// workaround, loop over right hand sides
	int_t nrhs = 1; int_t info;
	// int_t nrhs = nrhs_, info;
	superlu::
	  pdgssvx(&options, &A, &ScalePermstruct, x_ptr,
		  n, nrhs, &grid, &LUstruct, &SOLVEstruct,
		  berr.GetData(), &stat, &info);

	x_ptr += nloc;
      }
#endif
  }
#endif


  /********************************
   * MatrixSuperLU<complexdouble> *
   ********************************/
  

#ifndef SELDON_WITH_SUPERLU_DIST    
  //! Returns the LU factorization.
  /*!
    \param[out] Lmat matrix L in the LU factorization.
    \param[out] Umat matrix U in the LU factorization.
    \param[in] permuted should the permuted matrices be provided? SuperLU
    permutes the rows and columns of the factorized matrix. If \a permuted is
    set to true, L and U are returned as SuperLU computed them, hence with
    permuted rows and columns. If \a permuted is set to false, the matrices L
    and U are "unpermuted" so that L times U is equal to the initial matrix.
  */
  template<class Prop, class Allocator>
  void MatrixSuperLU<complex<double> >
  ::GetLU(Matrix<complex<double>, Prop, ColSparse, Allocator>& Lmat,
	  Matrix<complex<double>, Prop, ColSparse, Allocator>& Umat,
	  bool permuted)
  {
#ifdef SELDON_WITH_SUPERLU_MT
    Lstore = static_cast<superlu::SCPformat*>(L.Store);
    Ustore = static_cast<superlu::NCPformat*>(U.Store);
#else
    Lstore = static_cast<superlu::SCformat*>(L.Store);
    Ustore = static_cast<superlu::NCformat*>(U.Store);
#endif

    int Lnnz = Lstore->nnz;
    int Unnz = Ustore->nnz;

    int m = U.nrow;
    int n = U.ncol;

    Vector<complex<double>, VectFull, Allocator> Lval(Lnnz);
    Vector<int> Lrow(Lnnz);
    Vector<long> Lcol(n + 1);

    Vector<complex<double>, VectFull, Allocator> Uval(Unnz);
    Vector<int> Urow(Unnz);
    Vector<long> Ucol(n + 1);

    int Lsnnz;
    int Usnnz;
    LUextract(&L, &U, reinterpret_cast<superlu::doublecomplex*>(Lval.GetData()),
	      Lrow.GetData(), Lcol.GetData(),
              reinterpret_cast<superlu::doublecomplex*>(Uval.GetData()),
	      Urow.GetData(), Ucol.GetData(), &Lsnnz, &Usnnz);

    Lmat.SetData(m, n, Lval, Lcol, Lrow);
    Umat.SetData(m, n, Uval, Ucol, Urow);

    if (!permuted)
      {
        Vector<int> row_perm_orig = perm_r;
        Vector<int> col_perm_orig = perm_c;

        Vector<int> row_perm(n);
        Vector<int> col_perm(n);
        row_perm.Fill();
        col_perm.Fill();

        Sort(row_perm_orig, row_perm);
        Sort(col_perm_orig, col_perm);

        ApplyInversePermutation(Lmat, row_perm, col_perm);
        ApplyInversePermutation(Umat, row_perm, col_perm);
      }

  }
  
  
  //! Returns the LU factorization.
  /*!
    \param[out] Lmat matrix L in the LU factorization.
    \param[out] Umat matrix U in the LU factorization.
    \param[in] permuted should the permuted matrices be provided? SuperLU
    permutes the rows and columns of the factorized matrix. If \a permuted is
    set to true, L and U are returned as SuperLU computed them, hence with
    permuted rows and columns. If \a permuted is set to false, the matrices L
    and U are "unpermuted" so that L times U is equal to the initial matrix.
    \note This method will first retrieve the L and U matrices in 'ColSparse'
    format and then convert them into 'RowSparse'.
  */
  template<class Prop, class Allocator>
  void MatrixSuperLU<complex<double> >
  ::GetLU(Matrix<complex<double>, Prop, RowSparse, Allocator>& Lmat,
	  Matrix<complex<double>, Prop, RowSparse, Allocator>& Umat,
	  bool permuted)
  {
    Lmat.Clear();
    Umat.Clear();

    Matrix<complex<double>, Prop, ColSparse, Allocator> Lmat_col;
    Matrix<complex<double>, Prop, ColSparse, Allocator> Umat_col;
    GetLU(Lmat_col, Umat_col, permuted);

    Copy(Lmat_col, Lmat);
    Lmat_col.Clear();
    Copy(Umat_col, Umat);
    Umat_col.Clear();
  }


  //! Returns the size of memory used by the current object
  size_t MatrixSuperLU<complex<double> >::GetMemorySize() const
  {
    size_t taille = sizeof(int)*(perm_r.GetM()+perm_c.GetM());
    if (this->n > 0)
      {
#ifdef SELDON_WITH_SUPERLU_MT
        superlu::superlu_memusage_t mem_usage;
        int_t panel_size = superlu::sp_ienv(1);
        superlu::superlu_zQuerySpace(nprocs, const_cast<superlu::SuperMatrix*>(&L),
                                     const_cast<superlu::SuperMatrix*>(&U),
                                     panel_size, &mem_usage);
#else
        superlu::mem_usage_t mem_usage;
        superlu::zQuerySpace(const_cast<superlu::SuperMatrix*>(&L),
                             const_cast<superlu::SuperMatrix*>(&U), &mem_usage);
#endif

        taille += mem_usage.total_needed;
      }
    
    return taille;
  }

  
  //! factorization of matrix in complex double precision using SuperLU
  template<class T0, class Prop, class Storage, class Allocator>
  void MatrixSuperLU<complex<double> >::
  FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
		  bool keep_matrix)
  {
    // clearing previous factorization
    Clear();

    // conversion in CSC format
    General prop;
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<complex<double> > Val;

    ConvertToCSC(mat, prop, Ptr, IndRow, Val, false);
    if (!keep_matrix)
      mat.Clear();

    FactorizeCSC(Ptr, IndRow, Val, false);
  }

  
  void MatrixSuperLU<complex<double> >
  ::FactorizeCSC(Vector<long>& Ptr, Vector<int>& IndRow,
		 Vector<complex<double> >& Val, bool sym)
  {
    // initializing parameters
    int_t panel_size, relax;
    int_t lwork = 0;
    Init(Ptr.GetM()-1, panel_size, relax);

    superlu::SuperMatrix AA;
    int_t nnz = IndRow.GetM();
    superlu::
      zCreate_CompCol_Matrix(&AA, n, n, nnz,
                             reinterpret_cast<superlu::doublecomplex*>(Val.GetData()),
                             reinterpret_cast<int_t*>(IndRow.GetData()),
                             reinterpret_cast<int_t*>(Ptr.GetData()),
                             superlu::SLU_NC, superlu::SLU_Z, superlu::SLU_GE);

    // We get renumbering vectors perm_r and perm_c.
    options.ColPerm = permc_spec;
    if (permc_spec != superlu::MY_PERMC)
      {
        perm_r.Reallocate(n);
        perm_c.Reallocate(n);
        perm_r.Fill();
        perm_c.Fill();
        
        superlu::get_perm_c(permc_spec, &AA, perm_c.GetData());        
      }
    
#ifdef SELDON_WITH_SUPERLU_MT
    superlu::SuperMatrix AC;
    superlu::trans_t  trans = superlu::NOTRANS;
    perm_r.Reallocate(n);
    perm_c.Reallocate(n);
    
    superlu::
      pzgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
                   diag_pivot_thresh, usepr, drop_tol,
                   perm_c.GetData(), perm_r.GetData(),
                   NULL, lwork, &AA, &AC, &options, &stat);
    
    int_t info;
    superlu::
      pzgstrf(&options, &AC, perm_r.GetData(), &L, &U, &stat, &info);
    
    info_facto = info;
    superlu::pxgstrf_finalize(&options, &AC);

    if (info_facto == 0 && display_info)
      {
	cout << "Memory used for factorization in MiB: "
	     << this->GetMemorySize() / (1024. * 1024.) << endl;
      }
    
#else
    superlu::SuperMatrix A;
    // permuting matrix 
    Vector<int> etree(n);
    superlu::sp_preorder(&options, &AA, perm_c.GetData(), etree.GetData(), &A);

    // factorisation
    superlu::zgstrf(&options, &A, relax, panel_size, etree.GetData(),
                    NULL, lwork, perm_c.GetData(), perm_r.GetData(), &L, &U,
                    &Glu, &stat, &info_facto);

    if (info_facto == 0 && display_info)
      {
	superlu::mem_usage_t mem_usage;
	Lstore = (superlu::SCformat *) L.Store;
	Ustore = (superlu::NCformat *) U.Store;
	cout << "Number of nonzeros in factor L = " << Lstore->nnz<<endl;
	cout << "Number of nonzeros in factor U = " << Ustore->nnz<<endl;
	cout << "Number of nonzeros in L+U     = "
             << Lstore->nnz + Ustore->nnz<<endl;
	superlu::zQuerySpace(&L, &U, &mem_usage);
	cout << "Memory used for factorization in MiB: "
	     << mem_usage.total_needed / (1024. * 1024.) << endl;
      }
        
    superlu::Destroy_CompCol_Permuted(&A);    
    
#endif
            
    // clearing matrices
    superlu::Destroy_CompCol_Matrix(&AA);

    Ptr.Nullify(); IndRow.Nullify(); Val.Nullify();    
  }


  //! Solves linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(Vector<complex<double>, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system A x = b or A^T x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(const SeldonTranspose& TransA,
        Vector<complex<double>, VectFull, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), 1);
  }


  //! Solves linear system A x = b or A^T x = b
  void MatrixSuperLU<complex<double> >
  ::Solve(const SeldonTranspose& TransA, complex<double>* x_ptr, int nrhs_)
  {
    superlu::trans_t trans = superlu::NOTRANS;
    if (TransA.Trans())
      trans = superlu::TRANS;
    
    int_t nb_rhs = nrhs_, info;
    superlu::
      zCreate_Dense_Matrix(&B, n, nb_rhs,
                           reinterpret_cast<superlu::doublecomplex*>(x_ptr),
                           n, superlu::SLU_DN, superlu::SLU_Z, superlu::SLU_GE);
    
#ifdef SELDON_WITH_SUPERLU_MT
    superlu::zgstrs(trans, &L, &U, perm_r.GetData(),
                    perm_c.GetData(), &B, &stat, &info);
#else
    superlu::zgstrs(trans, &L, &U, perm_c.GetData(),
                    perm_r.GetData(), &B, &stat, &info);
#endif

    superlu::Destroy_SuperMatrix_Store(&B);
  }


  //! Solves linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(Matrix<complex<double>, General, ColMajor, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system A x = b or A^T x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(const SeldonTranspose& TransA,
        Matrix<complex<double>, General, ColMajor, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), x.GetN());
  }

#else


  /*****************************************
   * Distributed version for complexdouble *
   *****************************************/


  //! Returns the size of memory used by the current object
  size_t MatrixSuperLU<complex<double> >::GetMemorySize() const
  {
    size_t size = 0;
#ifndef SELDON_WITH_SUPERLU_DOUBLE
    if (n > 0)
      {
        superlu::mem_usage_t mem_usage;
        superlu::zQuerySpace_dist(n, const_cast<superlu::LUstruct_t*>(&LUstruct),
                                  const_cast<superlu::gridinfo_t*>(&grid),
                                  const_cast<superlu::SuperLUStat_t*>(&stat), &mem_usage);
        
        size += mem_usage.total;
      }
#endif

    return size;
  }


  //! factorization of matrix in complex double precision using SuperLU
  template<class T0, class Prop, class Storage, class Allocator>
  void MatrixSuperLU<complex<double> >::
  FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
		  bool keep_matrix)
  {
    Vector<long> Ptr; Vector<int> IndRow;
    Vector<complex<double> > Val; General prop;
    ConvertToCSC(mat, prop, Ptr, IndRow, Val, false);
    if (!keep_matrix)
      mat.Clear();
    
    FactorizeCSC(Ptr, IndRow, Val, IsSymmetricMatrix(mat));
  }


  void MatrixSuperLU<complex<double> >
  ::FactorizeCSC(Vector<long>& Ptr, Vector<int>& IndRow,
		 Vector<complex<double> >& Val, bool sym)
  {  
    Vector<int> glob_num(1);
    glob_num(0) = 0;
    MPI_Comm comm = MPI_COMM_SELF;
    FactorizeDistributedMatrix(comm, Ptr, IndRow, Val,
                               glob_num, sym, false);
  }


  //! Solves linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(Vector<complex<double>, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system A x = b or A^T x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(const SeldonTranspose& TransA,
        Vector<complex<double>, VectFull, Allocator2>& x)
  {
    Vector<int> glob_num(1);
    glob_num(0) = 0;
    MPI_Comm comm = MPI_COMM_SELF;
    SolveDistributed(comm, TransA, x, glob_num);
  }


  void MatrixSuperLU<complex<double> >
  ::Solve(const SeldonTranspose& TransA,
	  complex<double>* x_ptr, int nrhs_)
  {
    Vector<int> glob_num(1);
    glob_num(0) = 0;
    MPI_Comm comm = MPI_COMM_SELF;
    SolveDistributed(comm, TransA, x_ptr, nrhs_, glob_num);
  }
  

  //! Solves linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(Matrix<complex<double>, General, ColMajor, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Solves linear system A x = b or A^T x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(const SeldonTranspose& TransA,
        Matrix<complex<double>, General, ColMajor, Allocator2>& x)
  {
    Vector<int> glob_num(1);
    glob_num(0) = 0;
    MPI_Comm comm = MPI_COMM_SELF;
    SolveDistributed(comm, TransA, x, glob_num);
  }


  void MatrixSuperLU<complex<double> >::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
                             Vector<int>& IndRow, Vector<complex<double> >& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    Vector<int> PtrInt(Ptr.GetM());
    for (int i = 0; i < Ptr.GetM(); i++)
      PtrInt(i) = Ptr(i);

    FactorizeParallel(comm_facto, PtrInt, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }

  
  void MatrixSuperLU<complex<double> >::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
                             Vector<int64_t>& IndRow, Vector<complex<double> >& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    FactorizeParallel(comm_facto, Ptr, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }
  
  
  template<class Tint>
  void MatrixSuperLU<complex<double> >::
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<Tint>&, Vector<Tint>&,
                    Vector<complex<double> >&,
                    const Vector<int>& glob_num,
                    bool sym, bool keep_matrix)
  {
    cout << "Not available for this type of integers " << endl;
    cout << "Size of Tint = " << sizeof(Tint) << endl;
    abort();
  }


  void MatrixSuperLU<complex<double> >::
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<superlu_int_t>& Ptr, Vector<superlu_int_t>& Row,
                    Vector<complex<double> >& Val,
                    const Vector<int>& glob_num,
                    bool sym, bool keep_matrix)
  {
#ifndef SELDON_WITH_SUPERLU_DOUBLE
    // previous factorization is cleared if present
    Clear();
    
    if (!sym)
      {
        cout << "Problem with option TRANS of SuperLU" << endl;
        abort();
      }
    
    // m_loc : local number of rows
    // N : global number of rows
    int m_loc_ = Ptr.GetM()-1;    
    int N = m_loc_; nloc = m_loc_;
    MPI_Allreduce(&m_loc_, &N, 1, MPI_INTEGER, MPI_SUM, comm_facto);
    int_t m_loc = m_loc_;

    // structures are initialized with N
    int_t panel_size, relax;
    Init(N, panel_size, relax);

    // 1-D grid
    MPI_Comm_size(comm_facto, &nprow);
    npcol = 1;
    
    // initialize the superlu process grid
    superlu::superlu_gridinit(comm_facto, nprow, npcol, &grid);
    
    // fills the superlu structure
    // global numbers are assumed to be consecutive, we provide the first row number
    int_t fst_row = glob_num(0);
    int_t nnz_loc = Row.GetM();
    //superlu::SuperMatrix A;
    superlu::
      zCreate_CompRowLoc_Matrix_dist(&A, n, n, nnz_loc, m_loc, fst_row,
                                     reinterpret_cast<superlu::doublecomplex*>
                                     (Val.GetData()),
                                     reinterpret_cast<int_t*>(Row.GetData()),
                                     reinterpret_cast<int_t*>(Ptr.GetData()), 
                                     superlu::SLU_NR_loc, superlu::SLU_Z,
                                     superlu::SLU_GE);
    
    // completes factorization
    int_t nrhs = 0;
    options.Trans = superlu::NOTRANS;
    options.Fact = superlu::DOFACT;
    superlu::
      pzgssvx(&options, &A, &ScalePermstruct,
              NULL, m_loc, nrhs, &grid,
              &LUstruct, &SOLVEstruct, NULL, &stat, &info_facto);
    
    Ptr.Nullify(); Row.Nullify(); Val.Nullify();
#else
    cout << "Recompile Seldon without SELDON_WITH_SUPERLU_DOUBLE" << endl;
    cout << "Double and complex<double> cannot be handled simultaneously" << endl;
    abort();
#endif
  }

  
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  SolveDistributed(MPI_Comm& comm_facto,
                   const SeldonTranspose& TransA,
                   Vector<complex<double>, VectFull, Allocator2>& x,
		   const Vector<int>& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), 1, glob_num);
  }

  
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  SolveDistributed(MPI_Comm& comm_facto,
                   const SeldonTranspose& TransA,
                   Matrix<complex<double>, General, ColMajor, Allocator2>& x,
		   const Vector<int>& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), x.GetN(), glob_num);
  }


  void MatrixSuperLU<complex<double> >::
  SolveDistributed(MPI_Comm& comm_facto,
		   const SeldonTranspose& TransA,
		   complex<double>* x_ptr, int nrhs_,
		   const IVect& glob_num)
  {
#ifndef SELDON_WITH_SUPERLU_DOUBLE
    options.Fact = superlu::FACTORED;
    // inverting transpose because we have provided columns instead of rows
    if (TransA.NoTrans())
      options.Trans = superlu::TRANS;
    else
      options.Trans = superlu::NOTRANS;
       
    // TRANS does not work, we put NOTRANS
    options.Trans = superlu::NOTRANS;
    Vector<double> berr(nrhs_);
    for (int k = 0; k < nrhs_; k++)
      {
	// workaround, loop over right hand sides
	int_t nrhs = 1; int_t info;
	//int_t nrhs = nrhs_; int_t info;
	superlu::
	  pzgssvx(&options, &A, &ScalePermstruct,
		  reinterpret_cast<superlu::doublecomplex*>(x_ptr),
		  n, nrhs, &grid, &LUstruct, &SOLVEstruct,
		  berr.GetData(), &stat, &info);

	x_ptr += nloc;
      }
#endif
  }
#endif

  
  /***************************
   * GetLU/SolveLU functions *
   ***************************/
  
  
  //! Factorization of a matrix of same type T as for the SuperLU object 
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixSuperLU<T>& mat_lu, bool keep_matrix, T& x)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  

  //! Factorization of a complex matrix with a real SuperLU object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixSuperLU<T>& mat_lu,
             bool keep_matrix, complex<T>& x)
  {
    throw WrongArgument("GetLU(Matrix<complex<T> >& A, MatrixSuperLU<T>& mat_lu, bool)",
			"The LU matrix must be complex");
  }


  //! Factorization of a real matrix with a complex SuperLU object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixSuperLU<complex<T> >& mat_lu, bool keep_matrix, T& x)
  {
    throw WrongArgument("GetLU(Matrix<T>& A, MatrixSuperLU<complex<T> >& mat_lu, bool)",
			"The sparse matrix must be complex");
  }
  

  //! Factorization of a general matrix with SuperLU
  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixSuperLU<T>& mat_lu,
	     bool keep_matrix)
  {
    // we check if the type of non-zero entries of matrix A
    // and of the SuperLU object (T) are different
    // we call one of the GetLUs written above
    // such a protection avoids to compile the factorisation of a complex
    // matrix with a real SuperLU object
    typename Matrix<T0, Prop, Storage, Allocator>::entry_type x;
    GetLU(A, mat_lu, keep_matrix, x);
  }
  

  //! LU resolution with a vector whose type is the same as for SuperLU object
  template<class T, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  //! LU resolution with a vector whose type is the same as for SuperLU object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! LU resolution with a matrix whose type is the same as for SuperLU object
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(SeldonNoTrans, x);
  }


  //! LU resolution with a matrix whose type is the same as for SuperLU object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! Solves A x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(MatrixSuperLU<double>& mat_lu,
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
	       MatrixSuperLU<double>& mat_lu,
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
  void SolveLU(MatrixSuperLU<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixSuperLU<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }


  //! Solves A x = b or A^T x = b, where A is complex and x is real => Forbidden  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixSuperLU<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }
  
}

#define SELDON_FILE_SUPERLU_CXX
#endif
