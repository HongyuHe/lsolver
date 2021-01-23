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


#ifndef SELDON_FILE_UMFPACK_HXX

extern "C"
{
#include "umfpack.h"
#include "colamd.h"
#include "camd.h"
}

#ifdef UMFPACK_INTSIZE64
#define umfpack_int_t SuiteSparse_long
#else
#define umfpack_int_t int
#endif


namespace Seldon
{
  //!< base class to solve linear system by using UmfPack
  template<class T>
  class MatrixUmfPack_Base : public VirtualSparseDirectSolver<T>
  {
  public :
    Vector<double> Control, Info; //!< parameters for UmfPack
    void *Symbolic, *Numeric ; //!< pointers of UmfPack objects
    int n; //!< number of rows in the matrix
    int print_level;
    bool transpose; //! transpose system to solve ?
    int status_facto;

  public :
    MatrixUmfPack_Base();

    void HideMessages();
    void ShowMessages();
    void ShowFullHistory();
    bool UseInteger8() const;

    int GetInfoFactorization() const;
    size_t GetMemorySize() const;

    void SelectOrdering(int type);
    void SetPermutation(const IVect&);
  };

  //! empty class
  template<class T>
  class MatrixUmfPack : public MatrixUmfPack_Base<T>
  {
  };

  //! class to solve linear system in double precision with UmfPack
  template<>
  class MatrixUmfPack<double> : public MatrixUmfPack_Base<double>
  {

  protected :
    //! arrays containing matrix pattern in csc format
    umfpack_int_t* ind_, *ptr_;
    //! non-zero values
    double* data_;

  public :

    MatrixUmfPack();
    ~MatrixUmfPack();

    void Clear();

    template<class T0, class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    void FactorizeCSC(Vector<umfpack_int_t>& Ptr, Vector<umfpack_int_t>& IndRow,
		      Vector<double>& Val, bool sym);
    
    template<class Prop, class Allocator>
    void PerformAnalysis(Matrix<double, Prop, RowSparse, Allocator> & mat);

    template<class Prop, class Allocator>
    void
    PerformFactorization(Matrix<double, Prop, RowSparse, Allocator> & mat);

    template<class Allocator2>
    void Solve(Vector<double, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose&, Vector<double, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
	       Matrix<double, General, ColMajor, Allocator2>& x);

    void Solve(const SeldonTranspose& TransA, double* x_ptr, int nrhs);
    
  };


  //! class to solve linear system in complex double precision with UmfPack
  template<>
  class MatrixUmfPack<complex<double> >
    : public MatrixUmfPack_Base<complex<double> >
  {

  protected:
    //! arrays containing matrix pattern in csc format
    umfpack_int_t* ptr_, *ind_;
    //! non-zero values
    double* data_real_, *data_imag_;

  public :

    MatrixUmfPack();
    ~MatrixUmfPack();

    void Clear();

    template<class T0, class Prop, class Storage, class Allocator>
    void
    FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
                    bool keep_matrix = false);

    void FactorizeCSC(Vector<umfpack_int_t>& Ptr, Vector<umfpack_int_t>& IndRow,
		      Vector<complex<double> >& Val, bool sym);

    template<class Allocator2>
    void Solve(Vector<complex<double>, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose&, Vector<complex<double>, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
	       Matrix<complex<double>, General, ColMajor, Allocator2>& x);

    void Solve(const SeldonTranspose& TransA, complex<double>* x_ptr, int nrhs);
    
  };

  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixUmfPack<T>& mat_lu,
	     bool keep_matrix = false);
  
  template<class T, class Allocator>
  void SolveLU(MatrixUmfPack<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
               MatrixUmfPack<T>& mat_lu, Vector<T, VectFull, Allocator>& x);
  
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixUmfPack<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixUmfPack<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x);
  
  template<class Allocator>
  void SolveLU(MatrixUmfPack<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixUmfPack<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixUmfPack<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixUmfPack<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x);

}

#define SELDON_FILE_UMFPACK_HXX
#endif
