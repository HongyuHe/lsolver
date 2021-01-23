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


#ifndef SELDON_FILE_ITERATIVE_GCR_CXX

namespace Seldon
{

  //! Solves a linear system by using Generalized Conjugate Residual (GCR)
  /*!
    Solves the linear system Ax = b with restarted Preconditioned
    Generalized Conjugate Residual Algorithm.

    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.

    See: Y. Saad, Iterative Methods for Sparse Linear System, PWS Publishing
    Company, 1996

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] outer Iteration parameters
  */
#ifdef SELDON_WITH_VIRTUAL
  template<class T, class Vector1>
  int Gcr(const VirtualMatrix<T>& A, Vector1& x, const Vector1& b,
	  Preconditioner_Base<T>& M,
	  Iteration<typename ClassComplexType<T>::Treal>& outer)
#else
  template <class Titer, class Matrix1, class Vector1, class Preconditioner>
  int Gcr(const Matrix1& A, Vector1& x, const Vector1& b,
	  Preconditioner& M, Iteration<Titer> & outer)
#endif
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;

    typedef typename Vector1::value_type Complexe;
    int m = outer.GetRestart();
    // we initialize outer
    int success_init = outer.Init(b);
    if (success_init != 0)
      return outer.ErrorCode();

    std::vector<Vector1> p(m+1, b), w(m+1, b);

    Vector<Complexe> beta(m+1);

    Vector1 r(b), q(b), u(b);
    Complexe zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);

    for (int i = 0; i < (m+1); i++)
      {
	p[i].Fill(zero);
	w[i].Fill(zero);
      }

    // we compute initial residual
    Copy(b,u);
    if (!outer.IsInitGuess_Null())
      outer.MltAdd(-one, A, x, one, u);
    else
      x.Fill(zero);

    M.Solve(A, u, r);

    Complexe alpha, delta;
    
    typedef typename ClassComplexType<Complexe>::Treal Treal;
    Treal normr = Norm2(r);
    outer.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (! outer.Finished(r))
      {
	// m is the maximum number of inner iterations
	Iteration<Treal> inner(outer);
	inner.SetNumberIteration(outer.GetNumberIteration());
	inner.SetMaxNumberIteration(outer.GetNumberIteration()+m);
	Copy(r, p[0]);
	Mlt(Treal(1)/normr, p[0]);

	int j = 0;

	while (! inner.Finished(r) )
	  {
	    // product matrix vector u=A*p(j)
	    outer.Mlt(A, p[j], u);

	    // preconditioning
	    M.Solve(A, u, w[j]);

	    beta(j) = DotProdConj(w[j], w[j]);
	    if (beta(j) == zero)
	      {
		outer.Fail(1, "Gcr breakdown #1");
		break;
	      }

	    // new iterate x = x + alpha*p(j) new residual r = r - alpha*w(j)
	    // where alpha = (conj(r_j),A*p_j)/(A*p_j,A*p_j)
	    alpha = DotProdConj(w[j], r) / beta(j);
	    Add(alpha, p[j], x);
	    Add(-alpha, w[j], r);

	    ++inner;
	    ++outer;
	    // product Matrix vector u = A*r
	    outer.Mlt(A, r, u);
	    M.Solve(A, u, q);

	    Copy(r, p[j+1]);
	    // we compute direction p(j+1) = r(j+1) +
	    // \sum_{i=0..j} ( -(A*r_j+1,A*p_i)/(A*p_i,A*p_i) p(i))
	    for (int i = 0; i <= j; i++)
	      {
		delta = -DotProdConj(w[i], q)/beta(i);
		Add(delta, p[i], p[j+1]);
	      }

	    ++inner;
	    ++outer;
	    ++j;
	  }
	normr = Norm2(r);
      }

    return outer.ErrorCode();
  }

} // end namespace

#define SELDON_FILE_ITERATIVE_GCR_CXX
#endif
