// Copyright (C) 2010 Lin Wu
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


#ifndef SELDON_FILE_ARPACKSOLVER_CXX

#include "ArpackSolver.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
ArpackLog debug_;
#endif

namespace Seldon
{

  //! Constructor.
  template<class T, class Y>
  ArpackSolver<T, Y>::ArpackSolver()
  {
  }


  //! Destructor.
  template<class T, class Y>
  ArpackSolver<T, Y>::~ArpackSolver()
  {
    Clear();
  }


  //! Initializations.
  /*!
    \param n dimension of the problem.
    \param nev number of eigenvalues to be computed.
    \param ncv number of Arnoldi vectors generated at each iteration.
    \param maxit maximum number of Arnoldi update iterations.
    \param tol stopping criterion (relative accuracy of Ritz values).
    \param solver_type type of solver. Set to "symmetric", "non-symmetric",
    "complex-single", or "complex-double"
    \param mode indicates the type of the eigenproblem (regular, shift and
    invert, etc).
    \param which specify which of the Ritz values of OP to compute.
    \param bmat standard ('I') or generalized ('G") eigenproblem?
    \param HowMny whether eigenvectors ('A') or schur vectors ('P') to be
    computed. Works with 'rvec' set to false.
    \param with_arpack_verbose whether ARPACK info is displayed.
   */
  template<class T, class Y>
  void ArpackSolver<T, Y>
  ::Init(int n, int nev, int ncv, int maxit, T tol,
	 string solver_type, int mode, string which,
	 char bmat, char HowMny, bool with_arpack_verbose)
  {
    n_ = n;
    nev_ = nev;
    ncv_ = ncv;
    maxit_ = maxit;
    tol_ = tol;
    solver_type_ = solver_type;
    mode_ = mode;
    which_ = which;
    bmat_ = bmat;
    HowMny_ = HowMny;

    // Initializations.
    info_ = 0;
    ido_ = 0;

    // Algorithm mode.
    ishfts_ = 1;

    iparam_[0] = ishfts_;
    iparam_[2] = maxit_;
    iparam_[6] = mode_;

    // Dimensions.
    ldv_ = n_;
    lworkl_ = ncv_ * (ncv_ + 8);

    // Check parameters.
    CheckParameter();

    // Allocations.
    Allocate();

    // Verbose.
    if (with_arpack_verbose)
      SetArpackVerbose();
  }


  //! Check parameters.
  template<class T, class Y>
  void ArpackSolver<T, Y>::CheckParameter()
  {
    if (solver_type_ != "symmetric" &&
	solver_type_ != "non-symmetric" &&
	solver_type_ != "complex-single" &&
	solver_type_ != "complex-double")
      throw Error("ArpackSolver::Continue",
		  "Unsupported solver type \"" + solver_type_ + "\".");
  }


  //! Clear memories.
  template<class T, class Y>
  void ArpackSolver<T, Y>::Clear()
  {
    n_ = 0;
    lworkl_ = 0;
    ldv_ = 0;

    Deallocate();
  }


  //! Allocates arrays.
  template<class T, class Y>
  void ArpackSolver<T, Y>::Allocate()
  {
    v_ = new Y[ldv_ * ncv_];
    workl_ = new Y[ncv_ * (ncv_ + 8)];
    workd_ = new Y[3 * n_];

    eig_val_ = new Y[ncv_];
    resid_ = new Y[n_];

    pselect_ = new int[ncv_];
  }


  //! Deallocates arrays.
  template<class T, class Y>
  void ArpackSolver<T, Y>::Deallocate()
  {
    if (v_) delete[] v_;
    if (workl_) delete[] workl_;
    if (workd_) delete[] workd_;
    
    if (eig_val_) delete[] eig_val_;
    if (resid_) delete[] resid_;

    if (pselect_) delete[] pselect_;
  }


  //! tells to Arpack to display debug information
  template<class T, class Y>
  void ArpackSolver<T, Y>::SetArpackVerbose()
  {
    SetArpackVerbose(solver_type_);
  }

  
  //! tells to Arpack to display debug information
  template<class T, class Y>
  void ArpackSolver<T, Y>::SetArpackVerbose(string solver_type)
  {  
    if (solver_type == "symmetric")
      {
	debug_.logfil = 6;
	debug_.ndigit = -3;
	debug_.msgets = 0;
	debug_.msaitr = 0; 
	debug_.msapps = 0;
	debug_.msaupd = 1;
	debug_.msaup2 = 0;
	debug_.mseigt = 0;
	debug_.mseupd = 0;
      }
    else if (solver_type == "non-symmetric")
      {
	debug_.logfil = 6;
	debug_.ndigit = -3;
	debug_.mgetv0 = 0;
	debug_.mnaupd = 1;
	debug_.mnaup2 = 0;
	debug_.mnaitr = 0;
	debug_.mneigt = 0;
	debug_.mnapps = 0;
	debug_.mngets = 0;
	debug_.mneupd = 0;
      }
    else if (solver_type == "complex-single" ||
	     solver_type == "complex-double")
      {
	debug_.logfil = 6;
	debug_.ndigit = -3;
	debug_.mgetv0 = 0;
	debug_.mcaupd = 1;
	debug_.mcaup2 = 0;
	debug_.mcaitr = 0;
	debug_.mceigt = 0;
	debug_.mcapps = 0;
	debug_.mcgets = 0;
	debug_.mceupd = 0;
      }
  }


  //! No ARPACK info is to be displayed.
  template<class T, class Y>
  void ArpackSolver<T, Y>::ClearArpackVerbose()
  {
    debug_.logfil = 6;
    debug_.ndigit = 0;
    debug_.mgetv0 = 0;
    debug_.msaupd = 0;
    debug_.msaup2 = 0;
    debug_.msaitr = 0;
    debug_.mseigt = 0;
    debug_.msapps = 0;
    debug_.msgets = 0;
    debug_.mseupd = 0;
    debug_.mnaupd = 0;
    debug_.mnaup2 = 0;
    debug_.mnaitr = 0;
    debug_.mneigt = 0;
    debug_.mnapps = 0;
    debug_.mngets = 0;
    debug_.mneupd = 0;
    debug_.mcaupd = 0;
    debug_.mcaup2 = 0;
    debug_.mcaitr = 0;
    debug_.mceigt = 0;
    debug_.mcapps = 0;
    debug_.mcgets = 0;
    debug_.mceupd = 0;
  }


  //! Gets the address of the first vector in the ARPACK working array.
  template<class T, class Y>
  Y* ArpackSolver<T, Y>::GetFirstWorkVector()
  {
    return workd_ + ipntr_[0] - 1;
  }


  //! Gets the address of the second vector in the ARPACK working array.
  template<class T, class Y>
  Y* ArpackSolver<T, Y>::GetSecondWorkVector()
  {
    return workd_ + ipntr_[1] - 1;
  }


  //! Gets one eigenvector.
  /*!
    \param index index of the eigenvector.
   */
  template<class T, class Y>
  Y* ArpackSolver<T, Y>::GetEigenVector(int index)
  {
    return v_ + index * ldv_;
  }


  //! Gets one eigenvalue.
  /*!
    \param index index of the eigenvalue.
   */
  template<class T, class Y>
  Y ArpackSolver<T, Y>::GetEigenValue(int index)
  {
    return eig_val_[index];
  }


  //! Gets the reverse communication flag.
  template<class T, class Y>
  int ArpackSolver<T, Y>::GetReverseCommunicationFlag()
  {
    return ido_;
  }


  //! Sets the reverse communication flag.
  /*!
    \param ido the reverse communication flag to be set.
  */
  template<class T, class Y>
  void ArpackSolver<T, Y>::SetReverseCommunicationFlag(int ido)
  {
    ido_ = ido;
  }

  
  //! Gets the info flag.
  template<class T, class Y>
  int ArpackSolver<T, Y>::GetInfoFlag()
  {
    return info_;
  }


  //! Sets the info flag.
  /*!
    \param info the info flag to set.
  */
  template<class T, class Y>
  void ArpackSolver<T, Y>::SetInfoFlag(int info)
  {
    info_ = info;
  }


  //! Gets the number of "converged" Ritz values.
  template<class T, class Y>
  int ArpackSolver<T, Y>::GetConvergedNumber()
  {
    return nconv_;
  }


  //! Calls ARPACK computation routine.
  template<class T, class Y>
  bool ArpackSolver<T, Y>::Continue()
  {
#ifdef SELDON_WITH_MPI
    int comm(0);
    if (solver_type_ == "symmetric")
      saupd(comm, ido_, bmat_, n_, (char *)which_.c_str(), nev_, tol_, resid_, ncv_,
	    v_, ldv_, iparam_, ipntr_, workd_, workl_, lworkl_, info_);
#else
    if (solver_type_ == "symmetric")
      saupd(ido_, bmat_, n_, (char *)which_.c_str(), nev_, tol_, resid_, ncv_,
	    v_, ldv_, iparam_, ipntr_, workd_, workl_, lworkl_, info_);
#endif
    
    return (ido_ != 99);
  }


  //! Post-processing.
  template<class T, class Y>
  bool ArpackSolver<T, Y>::Finish()
  {
    bool i_success = (info_ >= 0);

    if (i_success)
      {
	rvec_ = true;

#ifdef SELDON_WITH_MPI
    int comm(0);
    seupd(comm, int(rvec_), 'A', pselect_, eig_val_, v_, ldv_,
	  sigma_, bmat_, n_, (char *) which_.c_str(), nev_, tol_,
	  resid_, ncv_, v_, ldv_, iparam_, ipntr_, workd_, workl_,
	  lworkl_, ierr_);
#else
    seupd(int(rvec_), 'A', pselect_, eig_val_, v_, ldv_,
	  sigma_, bmat_, n_, (char *) which_.c_str(), nev_, tol_,
	  resid_, ncv_, v_, ldv_, iparam_, ipntr_, workd_, workl_,
	  lworkl_, ierr_);
#endif    
      }
    
    bool p_success = (ierr_ == 0);
    
    if (p_success)
      nconv_ = iparam_[4];

    return (i_success && p_success);
  }
  
} // namespace Seldon.


#define SELDON_FILE_ARPACKSOLVER_CXX
#endif
