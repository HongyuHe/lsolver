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


#ifndef SELDON_FILE_ITERATIVE_SYMMLQ_CXX

namespace Seldon
{

  //! Solves a linear system by using Symmetric LQ (SymmLQ)
  /*!
    Solves the real symmetric linear system Ax = b with
    restarted Preconditioned Symmetric LQ Algorithm.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See C. PAIGE AND M. SAUNDERS,
    Solution of sparse indefinite systems of linear equations,
    SIAM J. Numer. Anal., 12 (1975), pp. 617-629.

    \param[in] A  Real Symmetric Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
#ifdef SELDON_WITH_VIRTUAL
  template<class T, class Vector1>
  int Symmlq(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	     Preconditioner_Base<T>& M,
	     Iteration<typename ClassComplexType<T>::Treal>& iter)
#else
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Symmlq(const Matrix1& A, Vector1& x, const Vector1& b,
	     Preconditioner& M, Iteration<Titer> & iter)
#endif
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    Complexe alpha, beta, ibeta, beta_old, beta1,
      ceta, ceta_oold, ceta_old, ceta_bar;
    Complexe c, cold, s, sold, coold, soold, rho0, rho1, rho2, rho3, dp;
    Complexe zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    ceta = zero;
    
    Vector1 r(b), z(b), u(b), v(b), w(b), u_old(b), v_old(b), w_bar(b);

    typedef typename ClassComplexType<Complexe>::Treal Treal;
    Treal np, s_prod;
    u_old.Fill(zero); v_old.Fill(zero); w.Fill(zero); w_bar.Fill(zero);

    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();

    Copy(b, r);
    // r = b - A x
    if (!iter.IsInitGuess_Null())
      iter.MltAdd(-one, A, x, one, r);
    else
      x.Fill(zero);

    ceta_oold = zero; ceta_old = zero;
    c = one; cold = one; s = zero; sold = zero;
    M.Solve(A, r, z);

    dp = DotProd(r, z);
    dp = sqrt(dp); beta = dp; beta1 = beta;
    s_prod = abs(beta1);

    Copy(r, v); Copy(z, u);
    ibeta = one/beta;
    Mlt(ibeta, v); Mlt(ibeta, u);
    Copy(u, w_bar);
    np = Norm2(b);

    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (!iter.Finished(np))
      {
	// update
	if (!iter.First())
	  {
	    Copy(v, v_old); Copy(u, u_old);
	    ibeta = one/beta;
	    // v = ibeta r
	    // u = ibeta z
	    Copy(r, v); Mlt(ibeta, v);
	    Copy(z, u); Mlt(ibeta, u);
	    // w = c*w_bar + s*u
	    Copy(w_bar, w); Mlt(c, w); Add(s, u, w);
	    // w_bar = -s*w_bar + c*u
	    Mlt(-s,w_bar); Add(c, u, w_bar);
	    // x = x+ceta*w
	    Add(ceta, w, x);

	    ceta_oold = ceta_old;
	    ceta_old = ceta;
	  }

	// product matrix vector r = A u
	iter.Mlt(A, u, r);
	alpha = DotProd(u, r);
	// preconditioning
	M.Solve(A, r, z);

	// r = r - alpha*v
	// z = z - alpha*u
	Add(-alpha, v, r);
	Add(-alpha, u, z);

	// r = r - beta*v_old
	// z = z - beta*u_old
	Add(-beta, v_old, r);
	Add(-beta, u_old, z);

	beta_old = beta;
	dp = DotProd(r, z);
	beta = sqrt(dp);

	// QR factorization
	coold = cold; cold = c; soold = sold; sold = s;
	rho0 = cold * alpha - coold * sold * beta_old;    // gamma_bar
	rho1 = sqrt(rho0*rho0 + beta*beta);               // gamma
	rho2 = sold * alpha + coold * cold * beta_old;    // delta
	rho3 = soold * beta_old;                          // epsilon

	// Givens rotation
	c = rho0 / rho1; s = beta / rho1;

	if (iter.First())
	  ceta = beta1/rho1;
	else
	  ceta = -(rho2*ceta_old + rho3*ceta_oold)/rho1;

	s_prod *= abs(s);
	if (c == zero)
	  np = s_prod*1e16;
	else
	  np = s_prod/abs(c);

	++iter;
      }
    if (c == zero)
      ceta_bar = ceta*1e15;
    else
      ceta_bar = ceta/c;

    Add(ceta_bar, w_bar, x);

    return iter.ErrorCode();
  }

} // end namespace

#define SELDON_FILE_ITERATIVE_SYMMLQ_CXX
#endif
