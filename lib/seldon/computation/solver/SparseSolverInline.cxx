// Copyright (C) 2010 Vivien Mallet
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


#ifndef SELDON_FILE_COMPUTATION_SPARSESOLVER_INLINE_CXX


#include "SparseSolver.hxx"

namespace Seldon
{

  /*************************
   * Default Seldon solver *
   *************************/
  
  
  template<class T, class Allocator>
  inline SparseSeldonSolver<T, Allocator>::SparseSeldonSolver()
  {
    print_level = -1;
    symmetric_matrix = false;
    permtol = 0.1;
  }
  
    
  template<class T, class Allocator>
  inline bool SparseSeldonSolver<T, Allocator>::UseInteger8() const
  {
    return false;
  }


  template<class T, class Allocator>
  inline void SparseSeldonSolver<T, Allocator>::Clear()
  {
    mat_sym.Clear();
    mat_unsym.Clear();
  }
    
  
  template<class T, class Allocator>
  inline void SparseSeldonSolver<T, Allocator>::HideMessages()
  {
    print_level = -1;
  }
  
  
  template<class T, class Allocator>
  inline void SparseSeldonSolver<T, Allocator>::ShowMessages()
  {
    print_level = 1;
  }
  

  template<class T, class Allocator>
  inline int SparseSeldonSolver<T, Allocator>::GetPrintLevel() const
  {
    return print_level;
  }

  
  template<class T, class Allocator>
  inline double SparseSeldonSolver<T, Allocator>::GetPivotThreshold() const
  {
    return permtol;
  }
  
  
  template<class T, class Allocator>
  inline void SparseSeldonSolver<T, Allocator>::SetPivotThreshold(double a)
  {
    permtol = a;
  }
  

  template<class T, class Allocator>
  inline int SparseSeldonSolver<T, Allocator>::GetInfoFactorization() const
  {
    return 0;
  }

}  // namespace Seldon.


#define SELDON_FILE_COMPUTATION_SPARSESOLVER_INLINE_CXX
#endif
