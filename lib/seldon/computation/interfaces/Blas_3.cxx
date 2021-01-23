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


#ifndef SELDON_FILE_BLAS_3_CXX


#include "Blas_3.hxx"


namespace Seldon
{


  ////////////
  // MltAdd //


  /*** ColMajor and NoTrans ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const float& alpha,
		    const Matrix<float, Prop0, ColMajor, Allocator0>& A,
		    const Matrix<float, Prop1, ColMajor, Allocator1>& B,
		    const float& beta,
		    Matrix<float, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		C.GetM(), C.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetLD(), B.GetData(), B.GetLD(),
		beta, C.GetData(), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const double& alpha,
		    const Matrix<double, Prop0, ColMajor, Allocator0>& A,
		    const Matrix<double, Prop1, ColMajor, Allocator1>& B,
		    const double& beta,
		    Matrix<double, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		C.GetM(), C.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetLD(), B.GetData(), B.GetLD(),
		beta, C.GetData(), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const complex<float>& alpha,
		    const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		    const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
		    const complex<float>& beta,
		    Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		C.GetM(), C.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(B.GetData()), B.GetLD(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const complex<double>& alpha,
		    const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		    const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
		    const complex<double>& beta,
		    Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		C.GetM(), C.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(B.GetData()), B.GetLD(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetLD());
  }


  /*** ColMajor and TransA, TransB ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const float& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<float, Prop0, ColMajor, Allocator0>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<float, Prop1, ColMajor, Allocator1>& B,
		    const float& beta,
		    Matrix<float, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

    cblas_sgemm(CblasColMajor, TransA.Cblas(), TransB.Cblas(), C.GetM(), C.GetN(),
		A.GetN(TransA), alpha, A.GetData(), A.GetLD(),
		B.GetData(), B.GetLD(), beta, C.GetData(), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const double& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<double, Prop0, ColMajor, Allocator0>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<double, Prop1, ColMajor, Allocator1>& B,
		    const double& beta,
		    Matrix<double, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

    cblas_dgemm(CblasColMajor, TransA.Cblas(), TransB.Cblas(), C.GetM(), C.GetN(),
		A.GetN(TransA), alpha, A.GetData(), A.GetLD(),
		B.GetData(), B.GetLD(), beta, C.GetData(), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const complex<float>& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
		    const complex<float>& beta,
		    Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

    cblas_cgemm(CblasColMajor, TransA.Cblas(), TransB.Cblas(), C.GetM(), C.GetN(),
		A.GetN(TransA), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(B.GetData()), B.GetLD(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const complex<double>& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
		    const complex<double>& beta,
		    Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

    cblas_zgemm(CblasColMajor, TransA.Cblas(), TransB.Cblas(), C.GetM(), C.GetN(),
		A.GetN(TransA), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(B.GetData()), B.GetLD(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetLD());
  }


  /*** RowMajor and NoTrans ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const float& alpha,
		    const Matrix<float, Prop0, RowMajor, Allocator0>& A,
		    const Matrix<float, Prop1, RowMajor, Allocator1>& B,
		    const float& beta,
		    Matrix<float, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		C.GetM(), C.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetLD(), B.GetData(), B.GetLD(),
		beta, C.GetData(), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const double& alpha,
		    const Matrix<double, Prop0, RowMajor, Allocator0>& A,
		    const Matrix<double, Prop1, RowMajor, Allocator1>& B,
		    const double& beta,
		    Matrix<double, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		C.GetM(), C.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetLD(), B.GetData(), B.GetLD(),
		beta, C.GetData(), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const complex<float>& alpha,
		    const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		    const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
		    const complex<float>& beta,
		    Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		C.GetM(), C.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(B.GetData()), B.GetLD(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const complex<double>& alpha,
		    const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		    const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
		    const complex<double>& beta,
		    Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		C.GetM(), C.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(B.GetData()), B.GetLD(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetLD());
  }


  /*** RowMajor and TransA, TransB ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const float& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<float, Prop0, RowMajor, Allocator0>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<float, Prop1, RowMajor, Allocator1>& B,
		    const float& beta,
		    Matrix<float, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

    cblas_sgemm(CblasRowMajor, TransA.Cblas(), TransB.Cblas(), C.GetM(), C.GetN(),
		A.GetN(TransA), alpha, A.GetData(), A.GetLD(),
		B.GetData(), B.GetLD(), beta, C.GetData(), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const double& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<double, Prop0, RowMajor, Allocator0>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<double, Prop1, RowMajor, Allocator1>& B,
		    const double& beta,
		    Matrix<double, Prop2, RowMajor, Allocator2>& C)
  {
    
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

    cblas_dgemm(CblasRowMajor, TransA.Cblas(), TransB.Cblas(), C.GetM(), C.GetN(),
		A.GetN(TransA), alpha, A.GetData(), A.GetLD(),
		B.GetData(), B.GetLD(), beta, C.GetData(), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const complex<float>& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
		    const complex<float>& beta,
		    Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

    cblas_cgemm(CblasRowMajor, TransA.Cblas(), TransB.Cblas(), C.GetM(), C.GetN(),
		A.GetN(TransA), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(B.GetData()), B.GetLD(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAddMatrix(const complex<double>& alpha,
		    const SeldonTranspose& TransA,
		    const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		    const SeldonTranspose& TransB,
		    const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
		    const complex<double>& beta,
		    Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, TransB, B, C,
	     "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

    cblas_zgemm(CblasRowMajor, TransA.Cblas(), TransB.Cblas(), C.GetM(), C.GetN(),
		A.GetN(TransA), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(B.GetData()), B.GetLD(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetLD());
  }


  // MltAdd //
  ////////////



  ////////////
  // MltAdd //


  /*** ColSym and Upper ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const float& alpha,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_ssymm(CblasColMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const double& alpha,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_dsymm(CblasColMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_csymm(CblasColMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_zsymm(CblasColMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** ColSym and UpLo ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const float& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_ssymm(CblasColMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const double& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_dsymm(CblasColMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_csymm(CblasColMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_zsymm(CblasColMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** RowSym and Upper ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const float& alpha,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_ssymm(CblasRowMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetN(),
		beta, C.GetData(), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const double& alpha,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_dsymm(CblasRowMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetN(),
		beta, C.GetData(), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_csymm(CblasRowMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_zsymm(CblasRowMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  /*** RowSym and UpLo ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const float& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_ssymm(CblasRowMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetN(),
		beta, C.GetData(), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const double& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_dsymm(CblasRowMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetN(),
		beta, C.GetData(), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_csymm(CblasRowMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_zsymm(CblasRowMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  // MltAdd //
  ////////////



  ////////////
  // MltAdd //


  /*** ColHerm and Upper ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_chemm(CblasColMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_zhemm(CblasColMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** ColHerm and UpLo ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_chemm(CblasColMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_zhemm(CblasColMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** RowHerm and Upper ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_chemm(CblasRowMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, A, B, beta, C)");
#endif

    cblas_zhemm(CblasRowMajor, Side.Cblas(), CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  /*** RowHerm and UpLo ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_chemm(CblasRowMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, C, "MltAdd(side, alpha, uplo, A, B, beta, C)");
#endif

    cblas_zhemm(CblasRowMajor, Side.Cblas(), Uplo.Cblas(),
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  // MltAdd //
  ////////////



  /////////
  // Mlt //


  /*** ColUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_strmm(CblasColMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_dtrmm(CblasColMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_ctrmm(CblasColMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_ztrmm(CblasColMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_strmm(CblasColMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_dtrmm(CblasColMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_ctrmm(CblasColMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_ztrmm(CblasColMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_strmm(CblasColMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_dtrmm(CblasColMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_ctrmm(CblasColMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_ztrmm(CblasColMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_strmm(CblasColMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_dtrmm(CblasColMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_ctrmm(CblasColMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_ztrmm(CblasColMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** RowUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_strmm(CblasRowMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_dtrmm(CblasRowMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_ctrmm(CblasRowMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_ztrmm(CblasRowMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_strmm(CblasRowMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_dtrmm(CblasRowMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_ctrmm(CblasRowMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_ztrmm(CblasRowMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_strmm(CblasRowMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_dtrmm(CblasRowMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_ctrmm(CblasRowMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, A, B)");
#endif

    cblas_ztrmm(CblasRowMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_strmm(CblasRowMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_dtrmm(CblasRowMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_ctrmm(CblasRowMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Mlt(side, alpha, status, diag, A, B)");
#endif

    cblas_ztrmm(CblasRowMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  // Mlt //
  /////////



  ///////////
  // Solve //


  /*** ColUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_strsm(CblasColMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_dtrsm(CblasColMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_ctrsm(CblasColMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_ztrsm(CblasColMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_strsm(CblasColMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_dtrsm(CblasColMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_ctrsm(CblasColMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_ztrsm(CblasColMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_strsm(CblasColMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_dtrsm(CblasColMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_ctrsm(CblasColMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_ztrsm(CblasColMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_strsm(CblasColMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_dtrsm(CblasColMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_ctrsm(CblasColMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_ztrsm(CblasColMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** RowUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_strsm(CblasRowMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_dtrsm(CblasRowMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_ctrsm(CblasRowMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_ztrsm(CblasRowMajor, Side.Cblas(), CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_strsm(CblasRowMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_dtrsm(CblasRowMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_ctrsm(CblasRowMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_ztrsm(CblasRowMajor, Side.Cblas(), CblasUpper, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_strsm(CblasRowMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_dtrsm(CblasRowMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_ctrsm(CblasRowMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, A, B)");
#endif

    cblas_ztrsm(CblasRowMajor, Side.Cblas(), CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_strsm(CblasRowMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_dtrsm(CblasRowMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_ctrsm(CblasRowMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Side, A, B, "Solve(side, alpha, status, diag, A, B)");
#endif

    cblas_ztrsm(CblasRowMajor, Side.Cblas(), CblasLower, TransA.Cblas(), DiagA.Cblas(),
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  // Solve //
  ///////////


} // namespace Seldon.

#define SELDON_FILE_BLAS_3_CXX
#endif
