// Copyright (C) 2010 Marc Duruflé
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


#ifndef SELDON_FILE_CHOLMOD_HXX

extern "C"
{
#include "cholmod.h"
}


namespace Seldon
{

  //! Object containing Cholesky factorization.
  class MatrixCholmod : public VirtualSparseDirectSolver<double>
  {
  protected :
    cholmod_common param_chol;
    cholmod_factor* L;
    int n;

  public :
    MatrixCholmod();
    ~MatrixCholmod();

    bool UseInteger8() const;
    void Clear();

    void HideMessages();
    void ShowMessages();
    void ShowFullHistory();
    
    size_t GetMemorySize() const;
    int GetInfoFactorization() const;

    template<class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<double, Prop, Storage, Allocator> & mat,
                         bool keep_matrix = false);

    void FactorizeCSR(Matrix<double, Symmetric, RowSymSparse> & Acsc);

    template<class Allocator>
    void Solve(const SeldonTranspose& TransA,
               Vector<double, VectFull, Allocator>& x);

    void Solve(const SeldonTranspose&, double* x_ptr, int nrhs);
    
    template<class Allocator>
    void Mlt(const SeldonTranspose& TransA,
             Vector<double, VectFull, Allocator>& x);
    
  };

  template<class T, class Prop, class Storage, class Allocator>
  void GetCholesky(Matrix<T, Prop, Storage, Allocator>& A,
                   MatrixCholmod& mat_chol, bool keep_matrix = false);

  template<class T, class Allocator>
  void
  SolveCholesky(const SeldonTranspose& TransA,
                MatrixCholmod& mat_chol, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void
  MltCholesky(const SeldonTranspose& TransA,
	      MatrixCholmod& mat_chol, Vector<T, VectFull, Allocator>& x);
  
}

#define SELDON_FILE_CHOLMOD_HXX
#endif

