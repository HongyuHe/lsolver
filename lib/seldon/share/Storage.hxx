// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2003-2011 Marc Duruflé
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


#ifndef SELDON_FILE_STORAGE_HXX

namespace Seldon
{


  //////////////////////
  // GENERAL MATRICES //
  //////////////////////


#ifndef SWIG
  class ColMajor
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };
#endif


  class RowMajor
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };



  /////////////
  // VECTORS //
  /////////////


  class VectFull;
  class VectSparse;
  class Collection;
  class DenseSparseCollection;
  class PETScSeq;
  class PETScPar;
  class PETScSeqDense;
  class PETScMPIDense;
  class PETScMPIAIJ;


  ////////////
  // SPARSE //
  ////////////


#ifndef SWIG
  class ColSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = true;
  };
#endif


  class RowSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = true;
  };


#ifndef SWIG
  class ColSymSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = true;
  };


  class RowSymSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = true;
  };

  class ArrayRowSparse : public RowSparse
  {
  };

  class ArrayColSparse : public ColSparse
  {
  };

  class ArrayRowSymSparse : public RowSymSparse
  {
  };

  class ArrayColSymSparse : public ColSymSparse
  {
  };

  
  ///////////////
  // SYMMETRIC //
  ///////////////


  class ColSymPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };


  class RowSymPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };


  class ColSym
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };


  class RowSym
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };



  ///////////////
  // HERMITIAN //
  ///////////////


  class ColHerm
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };


  class RowHerm
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };


  class ColHermPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };


  class RowHermPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static const bool Sparse = false;
  };



  ////////////////
  // TRIANGULAR //
  ////////////////


  class ColUpTriang
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static bool UpLo();
    static const bool Sparse = false;
  };


  class ColLoTriang
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static bool UpLo();
    static const bool Sparse = false;
  };


  class RowUpTriang
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static bool UpLo();
    static const bool Sparse = false;
  };


  class RowLoTriang
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static bool UpLo();
    static const bool Sparse = false;
  };


  class ColUpTriangPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static bool UpLo();
    static const bool Sparse = false;
  };


  class ColLoTriangPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static bool UpLo();
    static const bool Sparse = false;
  };


  class RowUpTriangPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static bool UpLo();
    static const bool Sparse = false;
  };


  class RowLoTriangPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
    static bool UpLo();
    static const bool Sparse = false;
  };
#endif


  //////////////////
  // SUB-MATRICES //
  //////////////////


  template <class M>
  class SubStorage
  {
  public:
    static const bool Sparse = false;
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static int GetEndLoop(int m, int n, int i);
  };


  ///////////
  // TYPES //
  ///////////


  class FloatDouble
  {
  };


  ////////////////
  // COLLECTION //
  ////////////////


  class ColMajorCollection
  {
  };


  class RowMajorCollection
  {
  };


  class ColSymPackedCollection
  {
  };


  class RowSymPackedCollection
  {
  };


  class ColUpTriangPackedCollection
  {
  };


  class RowUpTriangPackedCollection
  {
  };

  // complex matrices (in matrix_sparse/complex)

  class ColComplexSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static const bool Sparse = true;
  };


  class RowComplexSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static const bool Sparse = true;
  };


  class ColSymComplexSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static const bool Sparse = true;
  };


  class RowSymComplexSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static int GetBeginLoop(int i);
    static const bool Sparse = true;
  };


  class ArrayRowComplexSparse : public RowComplexSparse
  {
  };

  class ArrayRowSymComplexSparse : public RowSymComplexSparse
  {
  };

  class ArrayColComplexSparse : public ColComplexSparse
  {
  };

  class ArrayColSymComplexSparse : public ColSymComplexSparse
  {
  };



} // namespace Seldon.

#define SELDON_FILE_STORAGE_HXX
#endif
