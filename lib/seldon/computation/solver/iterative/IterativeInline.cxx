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


#ifndef SELDON_FILE_ITERATIVE_INLINE_CXX

namespace Seldon
{

  //! Default constructor
  template<class T>
  inline Preconditioner_Base<T>::Preconditioner_Base()
  {
  }


#ifdef SELDON_WITH_VIRTUAL
  //! destructor
  template<class T>
  inline Preconditioner_Base<T>::~Preconditioner_Base()
  {
  }
  
  
  //! Solves M z = r
  /*!
    Identity preconditioner M = I
  */
  template<class T>
  inline void Preconditioner_Base<T>
  ::Solve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>& z)
  {
    Copy(r, z);
  }

  //! Solves M^t z = r
  /*!
    Identity preconditioner M = I
  */  
  template<class T>
  inline void Preconditioner_Base<T>
  ::TransSolve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>& z)
  {
    Copy(r, z);
  }


  //! sets parameters of preconditioning
  template<class T>
  inline void Preconditioner_Base<T>
  ::SetInputPreconditioning(const string& keyword, const Vector<string>& param)
  {
  }


  //! copies parameters of another preconditioning
  template<class T>
  inline void Preconditioner_Base<T>::CopyParameter(const Preconditioner_Base<T>&)
  {
  }

  
#else

  
  //! Solves M z = r
  /*!
    Identity preconditioner M = I
  */
  template<class T> template<class Matrix1, class Vector1>
  inline void Preconditioner_Base<T>
  ::Solve(const Matrix1& A, const Vector1& r, Vector1& z)
  {
    Copy(r, z);
  }


  //! Solves M^t z = r
  /*!
    Identity preconditioner M = I
  */
  template<class T> template<class Matrix1, class Vector1>
  inline void Preconditioner_Base<T>
  ::TransSolve(const Matrix1& A, const Vector1 & r, Vector1 & z)
  {
    Solve(A, r, z);
  }
#endif


  //! Computes y = beta y + alpha A x
  template<class T> template<class T1, class Matrix1, class Vector1>
  inline void Iteration<T>
  ::MltAdd(const T1& alpha, const Matrix1& A, const Vector1& x,
	   const T1& beta, Vector1& y)
  {
#ifdef SELDON_WITH_VIRTUAL
    A.MltAddVector(alpha, x, beta, y);
#else
    Seldon::MltAdd(alpha, A, x, beta, y);
#endif
  }


  //! Computes y = A x
  template<class T> template<class Matrix1, class Vector1>
  inline void Iteration<T>
  ::Mlt(const Matrix1& A, const Vector1& x, Vector1& y)
  {
#ifdef SELDON_WITH_VIRTUAL
    A.MltVector(x, y);
#else
    Seldon::Mlt(A, x, y);
#endif
  }


  //! Computes y = A x or y = A^T x
  template<class T> template<class Matrix1, class Vector1>
  inline void Iteration<T>
  ::Mlt(const class_SeldonTrans& trans,
	const Matrix1& A, const Vector1& x, Vector1& y)
  {
#ifdef SELDON_WITH_VIRTUAL
    A.MltVector(trans, x, y);
#else
    Seldon::Mlt(trans, A, x, y);
#endif
  }


} // end namespace

#define SELDON_FILE_ITERATIVE_INLINE_CXX
#endif

