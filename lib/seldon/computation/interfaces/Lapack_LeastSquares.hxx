// Copyright (C) 2012 Vivien Mallet
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


#ifndef SELDON_FILE_LAPACK_LEAST_SQUARES_HXX

/*
  Functions included in this file:

  xGEQRF   (GetQR, GetLQ)
  xGELQF   (GetQR, GetLQ)
  xGEQP3   (GetQR_Pivot)
  xORGQR   (GetQ_FromQR)
  xUNGQR   (GetQ_FromQR)
  xUNMQR   (MltQ_FromQR)
  xORMQR   (MltQ_FromQR)
  xORMQR + xTRSM   (SolveQR)
  ZUNMQR + ZTRSM   (SolveQR)
  xORMLQ + xTRSM   (SolveQR)
  ZUNMLQ + ZTRSM   (SolveQR)
  xTRSM + xORMLQ   (SolveLQ)
  ZTRSM + ZUNMLQ   (SolveLQ)
  xTRSM + xORMQR   (SolveLQ)
  ZTRSM + ZUNMQR   (SolveLQ)

  GetLowerTriangular
  GetUpperTriangular
*/

namespace Seldon
{


  ///////////
  // GETQR //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<float, Prop0, ColMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);
  
  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<double, Prop0, ColMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<float>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);


  /*** RowMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<float, Prop0, RowMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<double, Prop0, RowMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<float>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);
  

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);

  
  // GETQR //
  ///////////


  /////////////////
  // GETQR_PIVOT //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR_Pivot(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   Vector<int>& ipivot, LapackInfo& info = lapack_info);


  // GETQR_PIVOT //
  /////////////////
  

  /////////////////
  // GETQ_FROMQR //


  /*** ColMajor ***/
  

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<float, Prop0, ColMajor, Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);
  
  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<float, Prop0, ColMajor, Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<float, Prop0, RowMajor, Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<double, Prop0, RowMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<float, Prop0, RowMajor, Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<double, Prop0, RowMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromLQ(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info);
  
  
  // GETQ_FROMQR //
  /////////////////


  ///////////
  // GETLQ //


  /*** ColMajor ***/

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<float, Prop0, ColMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<double, Prop0, ColMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);

  
  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<float>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);
  

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);


  /*** RowMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<float, Prop0, RowMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<double, Prop0, RowMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<float>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);
  

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info);
  
  
  // GETQ_LQ //
  ////////////


  /////////////////
  // MLTQ_FROMQR //


  /*** ColMajor ***/

  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<float, Prop0, ColMajor, Allocator0>& A,		   
		   const Vector<float, VectFull, Allocator1>& tau,
		   Vector<float, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<double, Prop0, ColMajor, Allocator0>& A,		   
		   const Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<float, Prop0, ColMajor, Allocator0>& A,
		   const Vector<float, VectFull, Allocator1>& tau,
		   Matrix<float, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Matrix<complex<float>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   const Vector<double, VectFull, Allocator1>& tau,
		   Matrix<double, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<float, Prop0, ColMajor, Allocator0>& A,		   
		   const Vector<float, VectFull, Allocator1>& tau,
		   Vector<float, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<double, Prop0, ColMajor, Allocator0>& A,		   
		   const Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<float, Prop0, ColMajor, Allocator0>& A,
		   const Vector<float, VectFull, Allocator1>& tau,
		   Matrix<float, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Matrix<complex<float>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   const Vector<double, VectFull, Allocator1>& tau,
		   Matrix<double, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<float, Prop0, RowMajor, Allocator0>& A,		   
		   const Vector<float, VectFull, Allocator1>& tau,
		   Vector<float, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<double, Prop0, RowMajor, Allocator0>& A,		   
		   const Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<float, Prop0, RowMajor, Allocator0>& A,
		   const Vector<float, VectFull, Allocator1>& tau,
		   Matrix<float, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Matrix<complex<float>, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<double, Prop0, RowMajor, Allocator0>& A,
		   const Vector<double, VectFull, Allocator1>& tau,
		   Matrix<double, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromQR(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<float, Prop0, RowMajor, Allocator0>& A,		   
		   const Vector<float, VectFull, Allocator1>& tau,
		   Vector<float, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<double, Prop0, RowMajor, Allocator0>& A,		   
		   const Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<float, Prop0, RowMajor, Allocator0>& A,
		   const Vector<float, VectFull, Allocator1>& tau,
		   Matrix<float, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<float>, VectFull, Allocator1>& tau,
		   Matrix<complex<float>, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<double, Prop0, RowMajor, Allocator0>& A,
		   const Vector<double, VectFull, Allocator1>& tau,
		   Matrix<double, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void MltQ_FromLQ(const SeldonSide& side, const SeldonTranspose& trans,
		   const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		   const Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, RowMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info);


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
	       LapackInfo& info = lapack_info);


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  
  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<float>, VectFull, Allocator1>& tau,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  
  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);


  /*** RowMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  
  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<float>, VectFull, Allocator1>& tau,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);
  

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);


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
	       LapackInfo& info = lapack_info);


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  
  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<float>, VectFull, Allocator1>& tau,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);
  
  
  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);


  /*** RowMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  
  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<float>, VectFull, Allocator1>& tau,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);
  
  
  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);
  
  
  // SOLVELQ //
  /////////////

  
  template<class T, class Prop, class Storage, class Allocator,
	   class T2, class Prop2, class Storage2, class Allocator2>
  void GetLowerTriangular(const Matrix<T, Prop, Storage, Allocator>& A,
			  Matrix<T2, Prop2, Storage2, Allocator2>& B);

  template<class T, class Prop, class Storage, class Allocator,
	   class T2, class Prop2, class Storage2, class Allocator2>
  void GetUpperTriangular(const Matrix<T, Prop, Storage, Allocator>& A,
			  Matrix<T2, Prop2, Storage2, Allocator2>& B);
  
  template<class T, class Prop, class Storage, class Allocator,
	   class T2, class Prop2, class Storage2, class Allocator2,
	   class T3, class Prop3, class Storage3, class Allocator3>
  void GetQR(const Matrix<T, Prop, Storage, Allocator>& A,
	     Matrix<T2, Prop2, Storage2, Allocator2>& Q,
	     Matrix<T3, Prop3, Storage3, Allocator3>& R);

  template<class T, class Prop, class Storage, class Allocator,
	   class T2, class Prop2, class Storage2, class Allocator2,
	   class T3, class Prop3, class Storage3, class Allocator3>
  void GetLQ(const Matrix<T, Prop, Storage, Allocator>& A,
	     Matrix<T2, Prop2, Storage2, Allocator2>& L,
	     Matrix<T3, Prop3, Storage3, Allocator3>& Q);
  
  
} // end namespace

#define SELDON_FILE_LAPACK_LEAST_SQUARES_HXX
#endif

