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


#ifndef SELDON_FILE_DISTRIBUTED_CHOLESKY_SOLVER_HXX


namespace Seldon
{

  template<class T>
  class DistributedCholeskySolver : public SparseCholeskySolver<T>
  {
  protected :
    typedef typename ClassComplexType<T>::Treal Treal;

#ifdef SELDON_WITH_MPI
    // data associated with distributed matrix
    // see DistributedMatrix.hxx for a detailed description
    int nodl_scalar_, nb_unknowns_scal_;
    MPI_Comm comm_;
    IVect* ProcSharingRows_;
    Vector<IVect>* SharingRowNumbers_;
    IVect global_col_numbers, local_col_numbers;
    
    template<class T2>
    void AssembleVec(Vector<T2>& X) const;

    template<class T2>
    void AssembleVec(Matrix<T2, General, ColMajor>& A) const;
#endif

  public :
    
    DistributedCholeskySolver();

    template<class Prop0, class Storage0, class Allocator0>
    void Factorize(Matrix<T, Prop0, Storage0, Allocator0>& A,
		   bool keep_matrix = false);

#ifdef SELDON_WITH_MPI
    template<class Prop0, class Storage0, class Allocator0>
    void Factorize(DistributedMatrix<T, Prop0, Storage0, Allocator0>& A,
                   bool keep_matrix = false);

#endif
    
    template<class T1>
    void Solve(const SeldonTranspose&, Vector<T1>& x_solution, bool assemble = true);
    
    template<class T1>
    void Mlt(const SeldonTranspose&, Vector<T1>& x_solution, bool assemble = true);
        
  };
  
} // namespace Seldon.


#define SELDON_FILE_DISTRIBUTED_CHOLESKY_SOLVER_HXX
#endif
