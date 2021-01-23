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


#ifndef SELDON_FILE_SHARE_MATRIXFLAG_CXX


#include "MatrixFlag.hxx"
#include "MatrixFlagInline.cxx"


namespace Seldon
{

#ifndef SELDON_WITH_COMPILED_LIBRARY
  class_SeldonTrans SeldonTrans;
  class_SeldonNoTrans SeldonNoTrans;
  class_SeldonConjTrans SeldonConjTrans;

  class_SeldonNonUnit SeldonNonUnit;
  class_SeldonUnit SeldonUnit;

  SeldonUplo SeldonUpper(0);
  SeldonUplo SeldonLower(1);

  SeldonNorm SeldonNormInf(0);
  SeldonNorm SeldonNorm1(1);

  SeldonConjugate SeldonUnconj(false);
  SeldonConjugate SeldonConj(true);

  class_SeldonLeft SeldonLeft;
  class_SeldonRight SeldonRight;
#endif


}


#define SELDON_FILE_SHARE_MATRIXFLAG_CXX
#endif
