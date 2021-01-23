// Copyright (C) 2001-2010 Vivien Mallet
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


#ifndef SELDON_FILE_SHARE_MATRIXFLAG_INLINE_CXX


#include "MatrixFlag.hxx"


namespace Seldon
{


  /////////////////////
  // SELDONTRANSPOSE //
  /////////////////////


  inline SeldonTranspose::SeldonTranspose(int status)
  {
    status_ = status;
    if (status_ == 0)
      cblas_status_ = CblasTrans;
    else if (status_ == 1)
      cblas_status_ = CblasNoTrans;
    else
      cblas_status_ = CblasConjTrans;
  }


  inline SeldonTranspose::SeldonTranspose(const enum CBLAS_TRANSPOSE status):
    cblas_status_(status)
  {
    if (cblas_status_ == CblasTrans)
      status_ = 0;
    else if (cblas_status_ == CblasNoTrans)
      status_ = 1;
    else
      status_ = 2;
  }


  
  inline CBLAS_TRANSPOSE SeldonTranspose::Cblas() const
  {
    return cblas_status_;
  }


  inline char SeldonTranspose::Char() const
  {
    if (status_ == 0)
      return 'T';
    else if (status_ == 1)
      return 'N';
    else
      return 'C';
  }


  inline char SeldonTranspose::RevChar() const
  {
    if (status_ == 0)
      return 'N';
    else if (status_ == 1)
      return 'T';
    else
      return 'N';
  }


  inline bool SeldonTranspose::Trans() const
  {
    return (status_ == 0);
  }


  inline bool SeldonTranspose::NoTrans() const
  {
    return (status_ == 1);
  }


  inline bool SeldonTranspose::ConjTrans() const
  {
    return (status_ == 2);
  }


  inline class_SeldonTrans::class_SeldonTrans(): SeldonTranspose(0)
  {
  }


  inline class_SeldonNoTrans::class_SeldonNoTrans(): SeldonTranspose(1)
  {
  }


  inline class_SeldonConjTrans::class_SeldonConjTrans(): SeldonTranspose(2)
  {
  }


  ////////////////
  // SELDONDIAG //
  ////////////////


  inline SeldonDiag::SeldonDiag(int status)
  {
    status_ = status;
    if (status_ == 0)
      cblas_status_ = CblasNonUnit;
    else
      cblas_status_ = CblasUnit;
  }


  inline CBLAS_DIAG SeldonDiag::Cblas() const
  {
    return cblas_status_;
  }


  inline char SeldonDiag::Char() const
  {
    return (status_ == 0) ? 'N' : 'U';
  }


  inline bool SeldonDiag::NonUnit() const
  {
    return (status_ == 0);
  }


  inline bool SeldonDiag::Unit() const
  {
    return (status_ == 1);
  }


  inline class_SeldonNonUnit::class_SeldonNonUnit(): SeldonDiag(0)
  {
  }


  inline class_SeldonUnit::class_SeldonUnit(): SeldonDiag(1)
  {
  }


  ////////////////
  // SELDONUPLO //
  ////////////////


  inline SeldonUplo::SeldonUplo(int status)
  {
    status_ = status;
    if (status_ == 0)
      cblas_status_ = CblasUpper;
    else
      cblas_status_ = CblasLower;
  }


  inline CBLAS_UPLO SeldonUplo::Cblas() const
  {
    return cblas_status_;
  }


  inline bool SeldonUplo::Upper() const
  {
    return (status_ == 0);
  }


  inline bool SeldonUplo::Lower() const
  {
    return (status_ == 1);
  }


  inline char SeldonUplo::Char() const
  {
    return (status_ == 0) ? 'U' : 'L';
  }


  inline char SeldonUplo::RevChar() const
  {
    return (status_ == 0) ? 'L' : 'U';
  }


  ////////////////
  // SELDONNORM //
  ////////////////


  inline SeldonNorm::SeldonNorm(int status)
  {
    status_ = status;
  }


  inline char SeldonNorm::Char() const
  {
    return (status_ == 0) ? 'I' : '1';
  }


  inline char SeldonNorm::RevChar() const
  {
    return (status_ == 0) ? '1' : 'I';
  }


  /////////////////////
  // SELDONCONJUGATE //
  /////////////////////


  inline SeldonConjugate::SeldonConjugate(bool status)
  {
    status_ = status;
  }


  inline bool SeldonConjugate::Conj() const
  {
    return status_;
  }


  ////////////////
  // SELDONSIDE //
  ////////////////


  inline SeldonSide::SeldonSide(int status)
  {
    status_ = status;
    if (status_ == 0)
      cblas_status_ = CblasLeft;
    else
      cblas_status_ = CblasRight;
  }


  inline SeldonSide::SeldonSide(const enum CBLAS_SIDE status):
    cblas_status_(status)
  {
    if (cblas_status_ == CblasLeft)
      status_ = 0;
    else
      status_ = 1;
  }


  inline CBLAS_SIDE SeldonSide::Cblas() const
  {
    return cblas_status_;
  }


  inline char SeldonSide::Char() const
  {
    if (status_ == 0)
      return 'L';
    else
      return 'R';
  }
  

  inline char SeldonSide::RevChar() const
  {
    if (status_ == 0)
      return 'R';
    else
      return 'L';
  }

  
  inline bool SeldonSide::Left() const
  {
    return (status_ == 0);
  }


  inline bool SeldonSide::Right() const
  {
    return (status_ == 1);
  }


  inline class_SeldonLeft::class_SeldonLeft(): SeldonSide(0)
  {
  }


  inline class_SeldonRight::class_SeldonRight(): SeldonSide(1)
  {
  }

}


#define SELDON_FILE_SHARE_MATRIXFLAG_INLINE_CXX
#endif
