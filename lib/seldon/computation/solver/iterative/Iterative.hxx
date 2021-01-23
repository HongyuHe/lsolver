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


#ifndef SELDON_FILE_ITERATIVE_HXX

namespace Seldon
{
  //! Base class for preconditioners
  template<class T>
  class Preconditioner_Base
  {
  public :

    Preconditioner_Base();    
    
#ifdef SELDON_WITH_VIRTUAL
    virtual ~Preconditioner_Base();
    
    virtual void Solve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>& z);
    virtual void TransSolve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>& z);
    
    virtual void SetInputPreconditioning(const string&, const Vector<string>&);
    virtual void CopyParameter(const Preconditioner_Base<T>&);
#else
    // solving M z = r
    template<class Matrix1, class Vector1>
    void Solve(const Matrix1& A, const Vector1 & r, Vector1 & z);

    // solving M^t z = r
    template<class Matrix1, class Vector1>
    void TransSolve(const Matrix1& A, const Vector1& r, Vector1 & z);
#endif

  };


  //! Class containing parameters for an iterative resolution
  /*!
    Titer is the precision (float or double), the solved
    linear system can be real or complex
  */
  template<class Titer>
  class Iteration
  {
  protected :
    Titer tolerance; //!< stopping criterion
    Titer facteur_reste; //!< inverse of norm of first residual
    int max_iter; //!< maximum number of iterations
    int nb_iter; //!< number of iterations
    int error_code; //!< error code returned by iterative solver
    bool fail_convergence; //!< true if the iterative solver has converged
    //! print level
    /*!
      0 -> no display
      1 -> displays residual after each 100 iterations
      6 -> displays residual after each iteration
    */
    int print_level;
    bool init_guess_null; //!< true if initial guess is null
    int type_solver; //!< iterative solver used
    int parameter_restart; //!< restart parameter (for Gmres and Gcr)
    int type_preconditioning; //!< preconditioner used
    string file_name_history;
    
  public :

    Iteration();
    Iteration(int max_iteration, const Titer& tol);
    Iteration(const Iteration<Titer>& outer);

    int GetTypeSolver() const;
    int GetRestart() const;
    Titer GetFactor() const;
    Titer GetTolerance() const;
    int GetNumberIteration() const;

    void SetSolver(int type_resolution, int param_restart, int type_prec);
    void SetRestart(int m);
    void SetTolerance(Titer stopping_criterion);
    void SetMaxNumberIteration(int max_iteration);
    void SetNumberIteration(int nb);

    void ShowMessages();
    void ShowFullHistory();
    void SaveFullHistory(const string& file);
    void HideMessages();

    template<class Vector1>
    int Init(const Vector1& r);
    bool First() const;

    bool IsInitGuess_Null() const;
    void SetInitGuess(bool type);

    template<class Vector1>
    bool Finished(const Vector1& r) const;
    bool Finished(const Titer& r) const;

    void Fail(int i, const string& s);

    Iteration& operator ++ (void);

    int ErrorCode();

    template<class T1, class Matrix1, class Vector1>
    void MltAdd(const T1& alpha, const Matrix1& A, const Vector1& x,
                const T1& beta, Vector1& y);
    
    template<class Matrix1, class Vector1>
    void Mlt(const Matrix1& A, const Vector1& x, Vector1& y);
    
    template<class Matrix1, class Vector1>
    void Mlt(const class_SeldonTrans& trans,
	     const Matrix1& A, const Vector1& x, Vector1& y);

  };

  // declarations of all iterative solvers
#ifdef SELDON_WITH_VIRTUAL
  template<class T, class Vector1>
  int BiCg(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	   Preconditioner_Base<T>& M,
	   Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int BiCgStab(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	       Preconditioner_Base<T>& M,
	       Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int BiCgStabl(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
		Preconditioner_Base<T>& M,
		Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int BiCgcr(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	     Preconditioner_Base<T>& M,
	     Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int Cg(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	 Preconditioner_Base<T>& M,
	 Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int Cgne(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	   Preconditioner_Base<T>& M,
	   Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int Cgs(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	  Preconditioner_Base<T>& M,
	  Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int CoCg(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	   Preconditioner_Base<T>& M,
	   Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int Gcr(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	  Preconditioner_Base<T>& M,
	  Iteration<typename ClassComplexType<T>::Treal>& outer);

  template<class T, class Vector1>
  int Gmres(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	    Preconditioner_Base<T>& M,
	    Iteration<typename ClassComplexType<T>::Treal>& outer);

  template<class T, class Vector1>
  int Lsqr(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	   Preconditioner_Base<T>& M,
	   Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int MinRes(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	     Preconditioner_Base<T>& M,
	     Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int QCgs(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	   Preconditioner_Base<T>& M,
	   Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int Qmr(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	  Preconditioner_Base<T>& M,
	  Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int QmrSym(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	     Preconditioner_Base<T>& M,
	     Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int Symmlq(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	     Preconditioner_Base<T>& M,
	     Iteration<typename ClassComplexType<T>::Treal>& iter);

  template<class T, class Vector1>
  int TfQmr(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	    Preconditioner_Base<T>& M,
	    Iteration<typename ClassComplexType<T>::Treal>& iter);

#else
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int BiCg(const Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int BiCgStab(const Matrix1& A, Vector1& x, const Vector1& b,
	       Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int BiCgStabl(const Matrix1& A, Vector1& x, const Vector1& b,
		Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int BiCgcr(const Matrix1& A, Vector1& x, const Vector1& b,
	     Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Cg(const Matrix1& A, Vector1& x, const Vector1& b,
	 Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Cgne(const Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Cgs(const Matrix1& A, Vector1& x, const Vector1& b,
	  Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int CoCg(const Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Gcr(const Matrix1& A, Vector1& x, const Vector1& b,
	  Preconditioner& M, Iteration<Titer> & outer);

  template <class Titer, class MatrixSparse, class Vector1, class Preconditioner>
  int Gmres(const MatrixSparse& A, Vector1& x, const Vector1& b,
	    Preconditioner& M, Iteration<Titer> & outer);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Lsqr(const Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int MinRes(const Matrix1& A, Vector1& x, const Vector1& b,
	     Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int QCgs(const Matrix1& A, Vector1& x, const Vector1& b,
	   Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Qmr(const Matrix1& A, Vector1& x, const Vector1& b,
	  Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int QmrSym(const Matrix1& A, Vector1& x, const Vector1& b,
	     Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Symmlq(const Matrix1& A, Vector1& x, const Vector1& b,
	     Preconditioner& M, Iteration<Titer> & iter);

  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int TfQmr(const Matrix1& A, Vector1& x, const Vector1& b,
	    Preconditioner& M, Iteration<Titer> & iter);

#endif

} // end namespace

#define SELDON_FILE_ITERATIVE_HXX
#endif
