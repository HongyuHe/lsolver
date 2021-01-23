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


#ifndef SELDON_FILE_ARPACKSOLVER_HXX

namespace Seldon
{

  
  //! Driver for ARPARK solvers. 
  template<class T, class Y>
  class ArpackSolver
  {
  protected:

    /*** Parameters ***/

    //! Dimension of the problem.
    int n_;
    //! Number of eigenvalues to be computed.
    int nev_;
    //! Number of Arnoldi vectors generated at each iteration.
    int ncv_;
    //! Maximum number of Arnoldi update iterations.
    int maxit_;
    //! Stopping criterion (relative accuracy of Ritz values).
    T tol_;
    //! Type of solver; set to "symmetric", "non-symmetric", "complex-single",
    //! or "complex-double".
    string solver_type_;
    //! Indicates the type of the eigenproblem (regular, shift and invert, etc).
    int mode_;
    //! Specify which of the Ritz values of OP to compute.      
    string which_;
    //! Standard ('I') or generalized ('G") eigenproblem?
    char bmat_;
    /*! Whether eigenvectors ('A') or schur vectors ('P') to be computed.
      Works with 'rvec' set to false. */
    char HowMny_;
    //! Shift.
    Y sigma_;
    //! Residual vector.
    Y *resid_;
    //! With shift?
    int ishfts_;
    //! Leading dimension of the Arnoldi basis / Schur vectors (v_).
    int ldv_;

    /*** Variables ***/

    /*! Logical array with dimension equal to the number of Arnoldi vectors
      (ncv_). */
    int *pselect_;
    //! Error flag on output.
    int ierr_;
    //! Eigenvectors/schur vectors to be computed?
    bool rvec_;
    //! ARPACK reverse communication flag.
    int ido_;
    //! ARPACK error flag.
    int info_;
    //! Dimension of array workl.
    int lworkl_;
    //! Dimension of array workv.
    int lworkv_;
    //! Dimension of array rwork.
    int lrwork_;
    //! Vector that handles original ARPACK parameters.
    int iparam_[11];
    //! Vector that handles original ARPACK pointers.
    int ipntr_[14];
    //! ARPACK internal vector.
    T *rwork_;
    //! ARPACK internal vector.
    Y *workl_;
    //! ARPACK internal vector.
    Y *workd_;
    //! ARPACK internal vector.
    Y *workev_;
    //! Arnoldi basis / Schur vectors.
    Y *v_;

    /*** Output variables ***/

    //! Number of "converged" Ritz values.
    int nconv_;
    //! Eigenvalues.
    Y *eig_val_;
    //! Eigenvectors.
    Y *eig_vec_;

  public:

    /*** Constructor and destructor ***/

    ArpackSolver();
    ~ArpackSolver();

    /*** Initialization ***/

    void Init(int n, int nev, int ncv, int maxit, T tol, string solver_type,
	      int mode, string which, char bmat, char HowMny,
	      bool with_arpack_verbose = false);
    void CheckParameter();
    void Clear();
    void Allocate();
    void Deallocate();
    void SetArpackVerbose();
    static void SetArpackVerbose(string s);
    void ClearArpackVerbose();

    /*** Access methods ***/

    Y* GetFirstWorkVector();
    Y* GetSecondWorkVector();
    Y* GetEigenVector(int index);
    Y GetEigenValue(int index);
    int GetReverseCommunicationFlag();
    void SetReverseCommunicationFlag(int ido);
    int GetInfoFlag();
    void SetInfoFlag(int info);
    int GetConvergedNumber();

    /*** Methods **/

    bool Continue();
    bool Finish();
  };


} // namespace Seldon.


#define SELDON_FILE_ARPACKSOLVER_HXX
#endif
