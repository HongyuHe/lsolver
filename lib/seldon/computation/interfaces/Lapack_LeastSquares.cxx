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


#ifndef SELDON_FILE_LAPACK_LEAST_SQUARES_CXX


#include "Lapack_LeastSquares.hxx"


namespace Seldon
{
  
    
  ///////////
  // GETQR //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<float, Prop0, ColMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<float, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    sgeqrf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<double, Prop0, ColMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    dgeqrf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<float>, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    cgeqrf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    zgeqrf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }

  
  /*** RowMajor ***/
  

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<float, Prop0, RowMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<float, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    sgelqf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<double, Prop0, RowMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    dgelqf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<float>, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    cgelqf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    zgelqf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  // GETQR //
  ///////////


  /////////////////
  // GETQR_PIVOT //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR_Pivot(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   Vector<int>& ipivot, LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 4 * max(m, n);
    ipivot.Fill(0);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    dgeqp3_(&m, &n, A.GetData(), &m, ipivot.GetData(), tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  // GETQR_PIVOT //
  /////////////////


  /////////////////
  // GETQ_FROMQR //


  /*** ColMajor ***/

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<float, Prop0, ColMajor, Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    int k = min(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < k)
      throw WrongDim("GetQ_FromQR", "tau not large enough");
#endif

    Vector<float, VectFull, Allocator1> work(lwork);
    sorgqr_(&m, &k, &k, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
    
    if (m < n)
      A.Resize(m, m);
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    int k = min(m, n);
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < k)
      throw WrongDim("GetQ_FromQR", "tau not large enough");
#endif

    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    cungqr_(&m, &k, &k, A.GetDataVoid(), &m, tau.GetData(),
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
    
    if (m < n)
      A.Resize(m, m);
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    int k = min(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < k)
      throw WrongDim("GetQ_FromQR", "tau not large enough");
#endif

    Vector<double, VectFull, Allocator1> work(lwork);
    dorgqr_(&m, &k, &k, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
    
    if (m < n)
      A.Resize(m, m);
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    int k = min(m, n);
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < k)
      throw WrongDim("GetQ_FromQR", "tau not large enough");
#endif

    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    zungqr_(&m, &k, &k, A.GetDataVoid(), &m, tau.GetData(),
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
    
    if (m < n)
      A.Resize(m, m);
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<float, Prop0, ColMajor, Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("GetQ_FromLQ", "tau not large enough");
#endif

    int k = min(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    sorglq_(&k, &n, &k, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());

    if (m > n)
      A.Resize(n, n);
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("GetQ_FromLQ", "tau not large enough");
#endif

    int k = min(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    cunglq_(&k, &n, &k, A.GetDataVoid(), &m, tau.GetDataVoid(),
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());

    if (m > n)
      A.Resize(n, n);
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("GetQ_FromLQ", "tau not large enough");
#endif

    int k = min(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    dorglq_(&k, &n, &k, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());

    if (m > n)
      A.Resize(n, n);
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("GetQ_FromLQ", "tau not large enough");
#endif

    int k = min(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    zunglq_(&k, &n, &k, A.GetDataVoid(), &m, tau.GetDataVoid(),
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());

    if (m > n)
      A.Resize(n, n);
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<float, Prop0, RowMajor, Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    int k = min(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < k)
      throw WrongDim("GetQ_FromQR", "tau not large enough");
#endif

    Vector<float, VectFull, Allocator1> work(lwork);
    sorglq_(&k, &m, &k, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
    
    if (n > m)
      A.Resize(m, m);
    
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    int k = min(m, n);
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < k)
      throw WrongDim("GetQ_FromQR", "tau not large enough");
#endif

    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    cunglq_(&k, &m, &k, A.GetDataVoid(), &n, tau.GetData(),
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());

    if (n > m)
      A.Resize(m, m);
  }
  
  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<double, Prop0, RowMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    int k = min(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < k)
      throw WrongDim("GetQ_FromQR", "tau not large enough");
#endif

    Vector<double, VectFull, Allocator1> work(lwork);
    dorglq_(&k, &m, &k, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
    
    if (n > m)
      A.Resize(m, m);
    
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    int k = min(m, n);
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < k)
      throw WrongDim("GetQ_FromQR", "tau not large enough");
#endif

    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    zunglq_(&k, &m, &k, A.GetDataVoid(), &n, tau.GetData(),
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());

    if (n > m)
      A.Resize(m, m);
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<float, Prop0, RowMajor, Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("GetQ_FromLQ", "tau not large enough");
#endif
    
    int k = min(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    sorgqr_(&n, &k, &k, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());

    if (m > n)
      A.Resize(n, n);    
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("GetQ_FromLQ", "tau not large enough");
#endif

    int k = min(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    cungqr_(&n, &k, &k, A.GetDataVoid(), &n, tau.GetDataVoid(),
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());

    if (m > n)
      A.Resize(n, n);
  }
  
  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<double, Prop0, RowMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("GetQ_FromLQ", "tau not large enough");
#endif
    
    int k = min(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    dorgqr_(&n, &k, &k, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());

    if (m > n)
      A.Resize(n, n);    
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("GetQ_FromLQ", "tau not large enough");
#endif

    int k = min(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    zungqr_(&n, &k, &k, A.GetDataVoid(), &n, tau.GetDataVoid(),
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());

    if (m > n)
      A.Resize(n, n);
  }

  
  // GETQ_FROMQR //
  /////////////////


  ///////////
  // GETLQ //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<float, Prop0, ColMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<float, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    sgelqf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<double, Prop0, ColMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    dgelqf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<float>, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    cgelqf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    zgelqf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  /*** RowMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<float, Prop0, RowMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<float, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization QR of A^t.
    sgeqrf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<double, Prop0, RowMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    dgeqrf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<float>, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    cgeqrf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    zgeqrf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  // GETLQ //
  ///////////


  /////////////////
  // MLTQ_FROMQR //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<float, Prop0, ColMajor, Allocator0>& A,		   
		   const Vector<float, VectFull, Allocator1>& tau,
		   Vector<float, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < lda)
      throw WrongDim("MltQ_FromQR", "b not large enough");
#endif
    
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.Char();
    sormqr_(&side, &trans_, &m, &n, &k, A.GetData(), &lda, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < lda)
      throw WrongDim("MltQ_FromQR", "b not large enough");
#endif
    
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.Char();
    cunmqr_(&side, &trans_, &m, &n, &k, A.GetDataVoid(), &lda, tau.GetDataVoid(),
	    b.GetDataVoid(), &m, work.GetDataVoid(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<double, Prop0, ColMajor, Allocator0>& A,		   
		   const Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < lda)
      throw WrongDim("MltQ_FromQR", "b not large enough");
#endif
    
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.Char();
    dormqr_(&side, &trans_, &m, &n, &k, A.GetData(), &lda, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < lda)
      throw WrongDim("MltQ_FromQR", "b not large enough");
#endif
    
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.Char();
    zunmqr_(&side, &trans_, &m, &n, &k, A.GetDataVoid(), &lda, tau.GetDataVoid(),
	    b.GetDataVoid(), &m, work.GetDataVoid(), &lwork,
	    &info.GetInfoRef());
  }

  
    template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<float, Prop0, ColMajor, Allocator0>& A,
		   const Vector<float, VectFull, Allocator1>& tau,
		   Matrix<float, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetM())
	 || (side.Right() && C.GetN() < A.GetM()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<float, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int lda = A.GetM();
    int k = min(A.GetM(), A.GetN());
    
    sormqr_(&side_, &trans_, &m, &n, &k, A.GetData(), &lda,
	    tau.GetData(), C.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Matrix<complex<float>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetM())
	 || (side.Right() && C.GetN() < A.GetM()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int lda = A.GetM();
    int k = min(A.GetM(), A.GetN());
    
    cunmqr_(&side_, &trans_, &m, &n, &k, A.GetDataVoid(), &lda,
	    tau.GetDataVoid(), C.GetDataVoid(), &m,
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   const Vector<double, VectFull, Allocator1>& tau,
		   Matrix<double, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetM())
	 || (side.Right() && C.GetN() < A.GetM()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<double, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int lda = A.GetM();
    int k = min(A.GetM(), A.GetN());
    
    dormqr_(&side_, &trans_, &m, &n, &k, A.GetData(), &lda,
	    tau.GetData(), C.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetM())
	 || (side.Right() && C.GetN() < A.GetM()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int lda = A.GetM();
    int k = min(A.GetM(), A.GetN());
    
    zunmqr_(&side_, &trans_, &m, &n, &k, A.GetDataVoid(), &lda,
	    tau.GetDataVoid(), C.GetDataVoid(), &m,
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<float, Prop0, ColMajor, Allocator0>& A,		   
		   const Vector<float, VectFull, Allocator1>& tau,
		   Vector<float, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetN())
      throw WrongDim("MltQ_FromLQ", "b not large enough");
#endif

    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.Char();
    sormlq_(&side, &trans_, &m, &n, &k, A.GetData(), &lda, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetN())
      throw WrongDim("MltQ_FromLQ", "b not large enough");
#endif

    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.Char();
    cunmlq_(&side, &trans_, &m, &n, &k, A.GetDataVoid(), &lda, tau.GetDataVoid(),
	    b.GetDataVoid(), &m, work.GetDataVoid(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<double, Prop0, ColMajor, Allocator0>& A,		   
		   const Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetN())
      throw WrongDim("MltQ_FromLQ", "b not large enough");
#endif

    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.Char();
    dormlq_(&side, &trans_, &m, &n, &k, A.GetData(), &lda, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetN())
      throw WrongDim("MltQ_FromLQ", "b not large enough");
#endif

    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.Char();
    zunmlq_(&side, &trans_, &m, &n, &k, A.GetDataVoid(), &lda, tau.GetDataVoid(),
	    b.GetDataVoid(), &m, work.GetDataVoid(), &lwork,
	    &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<float, Prop0, ColMajor, Allocator0>& A,
		   const Vector<float, VectFull, Allocator1>& tau,
		   Matrix<float, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetN())
	 || (side.Right() && C.GetN() < A.GetN()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<float, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int lda = A.GetM();
    int k = min(A.GetM(), A.GetN());
    
    sormlq_(&side_, &trans_, &m, &n, &k, A.GetData(), &lda,
	    tau.GetData(), C.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Matrix<complex<float>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetN())
	 || (side.Right() && C.GetN() < A.GetN()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int lda = A.GetM();
    int k = min(A.GetM(), A.GetN());
    
    cunmlq_(&side_, &trans_, &m, &n, &k, A.GetDataVoid(), &lda,
	    tau.GetDataVoid(), C.GetDataVoid(), &m,
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   const Vector<double, VectFull, Allocator1>& tau,
		   Matrix<double, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetN())
	 || (side.Right() && C.GetN() < A.GetN()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<double, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int lda = A.GetM();
    int k = min(A.GetM(), A.GetN());
    
    dormlq_(&side_, &trans_, &m, &n, &k, A.GetData(), &lda,
	    tau.GetData(), C.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetN())
	 || (side.Right() && C.GetN() < A.GetN()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int lda = A.GetM();
    int k = min(A.GetM(), A.GetN());
    
    zunmlq_(&side_, &trans_, &m, &n, &k, A.GetDataVoid(), &lda,
	    tau.GetDataVoid(), C.GetDataVoid(), &m,
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
  }
  
  
  /*** RowMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<float, Prop0, RowMajor, Allocator0>& A,		   
		   const Vector<float, VectFull, Allocator1>& tau,
		   Vector<float, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetM()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetM())
      throw WrongDim("MltQ_FromQR", "b not large enough");
#endif
    
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.RevChar();
    sormlq_(&side, &trans_, &m, &n, &k, A.GetData(), &lda, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetM()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetM())
      throw WrongDim("MltQ_FromQR", "b not large enough");
#endif
    
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_('C');
    if (trans.ConjTrans())
      trans_ = 'N';
    
    Conjugate(b);
    cunmlq_(&side, &trans_, &m, &n, &k, A.GetDataVoid(), &lda, tau.GetDataVoid(),
	    b.GetDataVoid(), &m, work.GetDataVoid(), &lwork,
	    &info.GetInfoRef());

    Conjugate(b);
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<double, Prop0, RowMajor, Allocator0>& A,		   
		   const Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetM()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetM())
      throw WrongDim("MltQ_FromQR", "b not large enough");
#endif
    
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.RevChar();
    dormlq_(&side, &trans_, &m, &n, &k, A.GetData(), &lda, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetM()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetM())
      throw WrongDim("MltQ_FromQR", "b not large enough");
#endif
    
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_('C');
    if (trans.ConjTrans())
      trans_ = 'N';
    
    Conjugate(b);
    zunmlq_(&side, &trans_, &m, &n, &k, A.GetDataVoid(), &lda, tau.GetDataVoid(),
	    b.GetDataVoid(), &m, work.GetDataVoid(), &lwork,
	    &info.GetInfoRef());

    Conjugate(b);
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<float, Prop0, RowMajor, Allocator0>& A,
		   const Vector<float, VectFull, Allocator1>& tau,
		   Matrix<float, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info)
  {
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetM())
	 || (side.Right() && C.GetN() < A.GetM()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<float, VectFull, Allocator1> work(lwork);
    char side_ = side.RevChar();
    char trans_ = trans.Char();
    int lda = A.GetN();
    int k = min(A.GetM(), A.GetN());
    
    sormlq_(&side_, &trans_, &n, &m, &k, A.GetData(), &lda,
	    tau.GetData(), C.GetData(), &n, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }
  

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Matrix<complex<float>, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info)
  {
        
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetM())
	 || (side.Right() && C.GetN() < A.GetM()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    char side_ = side.RevChar();
    char trans_ = trans.Char();
    int lda = A.GetN();
    int k = min(A.GetM(), A.GetN());
    
    cunmlq_(&side_, &trans_, &n, &m, &k, A.GetDataVoid(), &lda,
	    tau.GetDataVoid(), C.GetDataVoid(), &n,
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<double, Prop0, RowMajor, Allocator0>& A,
		   const Vector<double, VectFull, Allocator1>& tau,
		   Matrix<double, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info)
  {
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetM())
	 || (side.Right() && C.GetN() < A.GetM()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<double, VectFull, Allocator1> work(lwork);
    char side_ = side.RevChar();
    char trans_ = trans.Char();
    int lda = A.GetN();
    int k = min(A.GetM(), A.GetN());
    
    dormlq_(&side_, &trans_, &n, &m, &k, A.GetData(), &lda,
	    tau.GetData(), C.GetData(), &n, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }
  

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info)
  {
        
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetM())
	 || (side.Right() && C.GetN() < A.GetM()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    char side_ = side.RevChar();
    char trans_ = trans.Char();
    int lda = A.GetN();
    int k = min(A.GetM(), A.GetN());
    
    zunmlq_(&side_, &trans_, &n, &m, &k, A.GetDataVoid(), &lda,
	    tau.GetDataVoid(), C.GetDataVoid(), &n,
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
  }

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<float, Prop0, RowMajor, Allocator0>& A,		   
		   const Vector<float, VectFull, Allocator1>& tau,
		   Vector<float, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetM()))
      throw WrongDim("MltQ_FromLQ", "tau not large enough");

    if (b.GetM() < A.GetN())
      throw WrongDim("MltQ_FromLQ", "b not large enough");
#endif

    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.RevChar();
    sormqr_(&side, &trans_, &m, &n, &k, A.GetData(), &lda, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetM()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetN())
      throw WrongDim("MltQ_FromLQ", "b not large enough");
#endif

    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_('C');
    if (trans.ConjTrans())
      trans_ = 'N';
    
    Conjugate(b);
    cunmqr_(&side, &trans_, &m, &n, &k, A.GetDataVoid(), &lda, tau.GetDataVoid(),
	    b.GetDataVoid(), &m, work.GetDataVoid(), &lwork,
	    &info.GetInfoRef());
    
    Conjugate(b);
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<double, Prop0, RowMajor, Allocator0>& A,		   
		   const Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetM()))
      throw WrongDim("MltQ_FromLQ", "tau not large enough");

    if (b.GetM() < A.GetN())
      throw WrongDim("MltQ_FromLQ", "b not large enough");
#endif

    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_ = trans.RevChar();
    dormqr_(&side, &trans_, &m, &n, &k, A.GetData(), &lda, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   LapackInfo& info)
  {
    int lda = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(lda, A.GetM()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if (b.GetM() < A.GetN())
      throw WrongDim("MltQ_FromLQ", "b not large enough");
#endif

    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_('C');
    if (trans.ConjTrans())
      trans_ = 'N';
    
    Conjugate(b);
    zunmqr_(&side, &trans_, &m, &n, &k, A.GetDataVoid(), &lda, tau.GetDataVoid(),
	    b.GetDataVoid(), &m, work.GetDataVoid(), &lwork,
	    &info.GetInfoRef());
    
    Conjugate(b);
  }
  
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<float, Prop0, RowMajor, Allocator0>& A,
		   const Vector<float, VectFull, Allocator1>& tau,
		   Matrix<float, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info)
  {
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetN())
	 || (side.Right() && C.GetN() < A.GetN()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<float, VectFull, Allocator1> work(lwork);
    char side_ = side.RevChar();
    char trans_ = trans.Char();
    int lda = A.GetN();
    int k = min(A.GetM(), A.GetN());
    
    sormqr_(&side_, &trans_, &n, &m, &k, A.GetData(), &lda,
	    tau.GetData(), C.GetData(), &n, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }
  
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Matrix<complex<float>, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetN())
	 || (side.Right() && C.GetN() < A.GetN()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    char side_ = side.RevChar();
    char trans_ = trans.Char();
    int lda = A.GetN();
    int k = min(A.GetM(), A.GetN());
    
    cunmqr_(&side_, &trans_, &n, &m, &k, A.GetDataVoid(), &lda,
	    tau.GetDataVoid(), C.GetDataVoid(), &n,
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<double, Prop0, RowMajor, Allocator0>& A,
		   const Vector<double, VectFull, Allocator1>& tau,
		   Matrix<double, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info)
  {
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetN())
	 || (side.Right() && C.GetN() < A.GetN()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<double, VectFull, Allocator1> work(lwork);
    char side_ = side.RevChar();
    char trans_ = trans.Char();
    int lda = A.GetN();
    int k = min(A.GetM(), A.GetN());
    
    dormqr_(&side_, &trans_, &n, &m, &k, A.GetData(), &lda,
	    tau.GetData(), C.GetData(), &n, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }
  
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(A.GetM(), A.GetN()))
      throw WrongDim("MltQ_FromQR", "tau not large enough");

    if ( (side.Left() && C.GetM() < A.GetN())
	 || (side.Right() && C.GetN() < A.GetN()))
      throw WrongDim("MltQ_FromQR", "C not large enough");
#endif

    int m = C.GetM();
    int n = C.GetN();
    int lwork = max(A.GetM(), A.GetN());
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    char side_ = side.RevChar();
    char trans_ = trans.Char();
    int lda = A.GetN();
    int k = min(A.GetM(), A.GetN());
    
    zunmqr_(&side_, &trans_, &n, &m, &k, A.GetDataVoid(), &lda,
	    tau.GetDataVoid(), C.GetDataVoid(), &n,
	    work.GetDataVoid(), &lwork, &info.GetInfoRef());
  }
  
  
  // MLTQ_FROMQR //
  /////////////////


  /////////////
  // SOLVEQR //


  /*** ColMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();

    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('T');
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    // Computes Q^t b.
    sormqr_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Then solves R x = Q^t b.
    float alpha(1);
    cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
		CblasNonUnit, b.GetM(), nrhs,
		alpha, A.GetData(), A.GetM(), b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveQR", "tau not large enough");

    if (b.GetM() < m)
      throw WrongDim("SolveQR", "b not large enough");
#endif

    int nrhs = 1;
    char side('L');
    char trans('T');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    // Computes Q^t b.
    dormqr_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &lwork, work.GetData(), &lwork, &info.GetInfoRef());
    
    if (m > n)
      b.Resize(n);
    
    // Then solves R x = Q^t b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
		CblasNonUnit, b.GetM(), nrhs,
		alpha, A.GetData(), A.GetM(), b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<float>, VectFull, Allocator1>& tau,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveQR", "tau not large enough");

    if (b.GetM() < m)
      throw WrongDim("SolveQR", "b not large enough");
#endif

    int nrhs = 1;
    char side('L');
    char trans('C');
    int lwork = max(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    // Computes Q^t b.
    cunmqr_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    if (m > n)
      b.Resize(n);
    
    // Then solves R x = Q^t b.
    complex<float> alpha(1);
    cblas_ctrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
		CblasNonUnit, b.GetM(), nrhs,
		&alpha, A.GetData(), A.GetM(), b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveQR", "tau not large enough");

    if (b.GetM() < m)
      throw WrongDim("SolveQR", "b not large enough");
#endif

    int nrhs = 1;
    char side('L');
    char trans('C');
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    // Computes Q^t b.
    zunmqr_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    if (m > n)
      b.Resize(n);
    
    // Then solves R x = Q^t b.
    complex<double> alpha(1);
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
		CblasNonUnit, b.GetM(), nrhs,
		&alpha, A.GetData(), A.GetM(), b.GetData(), b.GetM());
  }


  /*** RowMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveQR", "tau not large enough");

    if (b.GetM() < m)
      throw WrongDim("SolveQR", "b not large enough");
#endif

    int k = tau.GetM();
    int nrhs = 1;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    // Computes Q b.
    sormlq_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    for (int i = m; i < n; i++)
      b(i) = 0;

    // Solves L^t y = b.
    float alpha(1);
    cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans,
		CblasNonUnit, n, nrhs,
		alpha, A.GetData(), n, b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveQR", "tau not large enough");

    if (b.GetM() < m)
      throw WrongDim("SolveQR", "b not large enough");
#endif
    
    int k = min(m, n);
    int nrhs = 1;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    // Computes Q b.
    dormlq_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());
    
    if (m > n)
      b.Resize(n);
    
    // Solves L^t y = b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans,
		CblasNonUnit, k, nrhs,
		alpha, A.GetData(), n, b.GetData(), b.GetM());
    
    
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveQR", "tau not large enough");

    if (b.GetM() < m)
      throw WrongDim("SolveQR", "b not large enough");
#endif
    
    int k = min(m, n);
    int nrhs = 2;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    Vector<double> rhs(2*m);
    for (int i = 0; i < m; i++)
      {
        rhs(i) = realpart(b(i));
        rhs(i+m) = imagpart(b(i));
      }
    
    // Computes Q b
    dormlq_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), rhs.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    for (int i = 0; i < n; i++)
      rhs(i+n) = rhs(i+m);

    rhs.Resize(2*n);
    
    // Solves L^t y = b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans,
		CblasNonUnit, k, nrhs,
		alpha, A.GetData(), n, rhs.GetData(), n);
    
    b.Reallocate(n);
    for (int i = 0; i < n; i++)
      b(i) = complex<double>(rhs(i), rhs(i+n));    
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<float>, VectFull, Allocator1>& tau,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveQR", "tau not large enough");

    if (b.GetM() < m)
      throw WrongDim("SolveQR", "b not large enough");
#endif

    int nrhs = 1;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    Conjugate(b);
    // Computes Q b.
    cunmlq_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    Conjugate(b);
    if (m > n)
      b.Resize(n);

    // Solves L^t y = b.
    complex<float> alpha(1);
    cblas_ctrsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans,
		CblasNonUnit, k, nrhs,
		&alpha, A.GetData(), n, b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveQR", "tau not large enough");

    if (b.GetM() < m)
      throw WrongDim("SolveQR", "b not large enough");
#endif

    int nrhs = 1;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    Conjugate(b);
    // Computes Q b.
    zunmlq_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    Conjugate(b);
    if (m > n)
      b.Resize(n);

    // Solves L^t y = b.
    complex<double> alpha(1);
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans,
		CblasNonUnit, k, nrhs,
		&alpha, A.GetData(), n, b.GetData(), b.GetM());
  }


  // SOLVEQR //
  /////////////


  /////////////
  // SOLVELQ //


  /*** ColMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {    
    int m = A.GetM();
    int n = A.GetN();
    if (m > n)
      throw Undefined("SolveLQ", "Function implemented only for n >= m");
    
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('T');
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    // Solves L y = b.
    float alpha(1);
    cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasNonUnit, m, nrhs,
		alpha, A.GetData(), m, b.GetData(), b.GetM());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Computes Q^t b.
    sormlq_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (m > n)
      throw Undefined("SolveLQ", "Function implemented only for n >= m");

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveLQ", "tau not large enough");

    if (b.GetM() < n)
      throw WrongDim("SolveLQ", "b not large enough");
#endif
 
    int k = tau.GetM();
    int nrhs = 1;
    char side('L');
    char trans('T');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    // Solves L y = b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasNonUnit, m, nrhs,
		alpha, A.GetData(), m, b.GetData(), b.GetM());
    
    // setting extra-components to 0
    for (int i = m; i < n; i++)
      b(i) = 0;

    // Applies Q^t b.
    dormlq_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<float>, VectFull, Allocator1>& tau,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (m > n)
      throw Undefined("SolveLQ", "Function implemented only for n >= m");
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveLQ", "tau not large enough");

    if (b.GetM() < n)
      throw WrongDim("SolveLQ", "b not large enough");
#endif

    int k = tau.GetM();
    int nrhs = 1;
    char side('L');
    char trans('C');
    int lwork = max(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    // Solve L y = b.
    complex<float> alpha(1);
    cblas_ctrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasNonUnit, m, nrhs,
		&alpha, A.GetData(), m, b.GetData(), b.GetM());

    // setting extra-components to 0
    for (int i = m; i < n; i++)
      b(i) = 0;

    // Applies Q^t.
    cunmlq_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (m > n)
      throw Undefined("SolveLQ", "Function implemented only for n >= m");
    
#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveLQ", "tau not large enough");

    if (b.GetM() < n)
      throw WrongDim("SolveLQ", "b not large enough");
#endif

    int k = tau.GetM();
    int nrhs = 1;
    char side('L');
    char trans('C');
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    // Solve L y = b.
    complex<double> alpha(1);
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasNonUnit, m, nrhs,
		&alpha, A.GetData(), m, b.GetData(), b.GetM());

    // setting extra-components to 0
    for (int i = m; i < n; i++)
      b(i) = 0;

    // Applies Q^t.
    zunmlq_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  /*** RowMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (m > n)
      throw Undefined("SolveLQ", "Function implemented only for n >= m");

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveLQ", "tau not large enough");

    if (b.GetM() < n)
      throw WrongDim("SolveLQ", "b not large enough");
#endif
    
    int k = tau.GetM();
    int nrhs = 1;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    // Solves R^t x = b.
    float alpha(1);
    cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans,
		CblasNonUnit, k, nrhs,
		alpha, A.GetData(), A.GetN(), b.GetData(), b.GetM());
    
    for (int i = m; i < n; i++)
      b(i) = 0;

    // Multiplies by Q.
    sormqr_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (m > n)
      throw Undefined("SolveLQ", "Function implemented only for n >= m");

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveLQ", "tau not large enough");

    if (b.GetM() < n)
      throw WrongDim("SolveLQ", "b not large enough");
#endif
    
    int k = tau.GetM();
    int nrhs = 1;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    // Solves R^t x = b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans,
		CblasNonUnit, k, nrhs,
		alpha, A.GetData(), A.GetN(), b.GetData(), b.GetM());

    for (int i = m; i < n; i++)
      b(i) = 0;

    // Multiplies by Q.
    dormqr_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<float>, VectFull, Allocator1>& tau,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (m > n)
      throw Undefined("SolveLQ", "Function implemented only for n >= m");

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveLQ", "tau not large enough");

    if (b.GetM() < n)
      throw WrongDim("SolveLQ", "b not large enough");
#endif
    
    int k = tau.GetM();
    int nrhs = 1;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<complex<float>, VectFull, Allocator1> work(lwork);
    // Solves R^t x = b.
    complex<float> alpha(1);
    cblas_ctrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans,
		CblasNonUnit,k, nrhs,
		&alpha, A.GetData(), A.GetN(), b.GetData(), b.GetM());

    for (int i = m; i < n; i++)
      b(i) = 0;
    
    Conjugate(b);
    // Computes Q b.
    cunmqr_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
    
    Conjugate(b);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (m > n)
      throw Undefined("SolveLQ", "Function implemented only for n >= m");

#ifdef SELDON_CHECK_DIMENSIONS
    if (tau.GetM() < min(m, n))
      throw WrongDim("SolveLQ", "tau not large enough");

    if (b.GetM() < n)
      throw WrongDim("SolveLQ", "b not large enough");
#endif
    
    int k = tau.GetM();
    int nrhs = 1;
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    // Solves R^t x = b.
    complex<double> alpha(1);
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans,
		CblasNonUnit,k, nrhs,
		&alpha, A.GetData(), A.GetN(), b.GetData(), b.GetM());

    for (int i = m; i < n; i++)
      b(i) = 0;
    
    Conjugate(b);
    // Computes Q b.
    zunmqr_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
    
    Conjugate(b);
  }


  // SOLVELQ //
  /////////////

  
  //! generic function to extract lower triangular part of a matrix
  /*!
    Equivalent Matlab function : B = tril(A)
  */
  template<class T, class Prop, class Storage, class Allocator,
	   class T2, class Prop2, class Storage2, class Allocator2>
  void GetLowerTriangular(const Matrix<T, Prop, Storage, Allocator>& A,
			  Matrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();

    T zero;
    SetComplexZero(zero);
    
    B.Clear();
    int k = min(m, n);
    B.Reallocate(m, k);
    B.Fill(zero);
    
    for (int i = 0; i < m; i++)
      for (int j = 0; j <= min(i, k-1); j++)
	B.Set(i, j, A(i, j));
    
  }


  //! generic function to extract upper triangular part of a matrix
  /*!
    Equivalent Matlab function : B = triu(A)
  */  
  template<class T, class Prop, class Storage, class Allocator,
	   class T2, class Prop2, class Storage2, class Allocator2>
  void GetUpperTriangular(const Matrix<T, Prop, Storage, Allocator>& A,
			  Matrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();

    T zero;
    SetComplexZero(zero);
    
    B.Clear();
    B.Reallocate(min(m, n), n);
    B.Fill(zero);
    
    for (int i = 0; i < min(m, n); i++)
      for (int j = i; j < n; j++)
	B.Set(i, j, A(i, j));
    
  }

  
  //! Generic function to explicit Q and R from QR factorisation of a matrix
  /*!
    Equivalent Matlab function : [Q, R] = qr(A)
   */
  template<class T, class Prop, class Storage, class Allocator,
	   class T2, class Prop2, class Storage2, class Allocator2,
	   class T3, class Prop3, class Storage3, class Allocator3>
  void GetQR(const Matrix<T, Prop, Storage, Allocator>& A,
	     Matrix<T2, Prop2, Storage2, Allocator2>& Q,
	     Matrix<T3, Prop3, Storage3, Allocator3>& R)
  {
    // QR-factorisation of A
    Vector<T> tau;
    Q = A;
    GetQR(Q, tau);
    
    // R is extracted
    GetUpperTriangular(Q, R);

    // then Q is explicited
    GetQ_FromQR(Q, tau);    
  }
  

  //! Generic function to explicit L, and Q from LQ factorisation of a matrix
  /*!
    Equivalent Matlab function : [L, Q] = lq(A)
   */
  template<class T, class Prop, class Storage, class Allocator,
	   class T2, class Prop2, class Storage2, class Allocator2,
	   class T3, class Prop3, class Storage3, class Allocator3>
  void GetLQ(const Matrix<T, Prop, Storage, Allocator>& A,
	     Matrix<T2, Prop2, Storage2, Allocator2>& L,
	     Matrix<T3, Prop3, Storage3, Allocator3>& Q)
  {
    // QR-factorisation of A
    Vector<T> tau;
    Q = A;
    GetLQ(Q, tau);
    
    // L is extracted
    GetLowerTriangular(Q, L);

    // then Q is explicited
    GetQ_FromLQ(Q, tau);
  }
  
} // end namespace

#define SELDON_FILE_LAPACK_LEAST_SQUARES_CXX
#endif
