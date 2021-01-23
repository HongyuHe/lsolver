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


#ifndef SELDON_FILE_ILUT_PRECONDITIONING_HXX

namespace Seldon
{

  template<class T, class Allocator
	   = typename SeldonDefaultAllocator<ArrayRowSparse, T>::allocator>
  class IlutPreconditioning : public Preconditioner_Base<T>,
			      public VirtualSparseDirectSolver<T>
  {
  protected :
    //! Verbosity level.
    int print_level;
    //! True if symmetric matrix is constructed.
    bool symmetric_algorithm;
    //! Type of incomplete factorization.
    int type_ilu;
    /*! \brief Maximum number of elements on a row of L or U.
      For Ilu(k), fill_level = k
    */
    int fill_level;
    /*! \brief Additional number of elements for each row.
      This number is only used for ILUT(k)
    */
    int additional_fill;
    //! Size of block where the pivot is searched.
    int mbloc;
    //! Diagonal compensation parameter (alpha = 0 -> ILU, alpha = 1 -> MILU).
    typename ClassComplexType<T>::Treal alpha;
    //! Threshold used for dropping small terms.
    typename ClassComplexType<T>::Treal droptol;
    //! Threshold for pivoting.
    typename ClassComplexType<T>::Treal permtol;
    //! Permutation arrays.
    IVect permutation_row, permutation_col;
    //! Symmetric matrix.
    Matrix<T, Symmetric, ArrayRowSymSparse, Allocator> mat_sym;
    //! Unsymmetric matrix.
    Matrix<T, General, ArrayRowSparse, Allocator> mat_unsym;

  public :

    //! Available types of incomplete factorization.
    enum {ILUT, ILU_D, ILUT_K, ILU_0, MILU_0, ILU_K};

    IlutPreconditioning();

    void Clear();
    
    bool UseInteger8() const;
    void HideMessages();
    void ShowMessages();

    int GetFactorisationType() const;
    int GetFillLevel() const;
    int GetAdditionalFillNumber() const;
    int GetPrintLevel() const;
    int GetPivotBlockInteger() const;
    size_t GetMemorySize() const;
    int GetInfoFactorization() const;

    void SetFactorisationType(int);
    void SetFillLevel(int);
    void SetAdditionalFillNumber(int);
    void SetPrintLevel(int);
    void SetPivotBlockInteger(int);
    void SetSymmetricAlgorithm();
    void SetUnsymmetricAlgorithm();

    typename ClassComplexType<T>::Treal GetDroppingThreshold() const;
    typename ClassComplexType<T>::Treal GetDiagonalCoefficient() const;
    typename ClassComplexType<T>::Treal GetPivotThreshold() const;

    void SetDroppingThreshold(typename ClassComplexType<T>::Treal);
    void SetDiagonalCoefficient(typename ClassComplexType<T>::Treal);
    void SetPivotThreshold(double);

    template<class MatrixSparse>
    void FactorizeSymMatrix(const IVect& perm,
                            MatrixSparse& mat, bool keep_matrix = false);

    template<class MatrixSparse>
    void FactorizeUnsymMatrix(const IVect& perm,
                              MatrixSparse& mat, bool keep_matrix = false);

    template<class T0, class Storage0, class Allocator0>
    void FactorizeMatrix(const IVect& perm,
                         Matrix<T0, General, Storage0, Allocator0>& mat,
                         bool keep_matrix = false);

    template<class T0, class Storage0, class Allocator0>
    void FactorizeMatrix(const IVect& perm,
                         Matrix<T0, Symmetric, Storage0, Allocator0>& mat,
                         bool keep_matrix = false);

#ifdef SELDON_WITH_VIRTUAL
    void Solve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>& z);
    void TransSolve(const VirtualMatrix<T>&, const Vector<T>& r, Vector<T>&);
#else
    template<class Matrix1, class Vector1>
    void TransSolve(const Matrix1& A, const Vector1& r, Vector1& z);

    template<class Matrix1, class Vector1>
    void Solve(const Matrix1& A, const Vector1& r, Vector1& z);
#endif

    template<class Vector1>
    void TransSolve(Vector1& z);

    template<class Vector1>
    void Solve(Vector1& z);

    template<class Vector1>
    void Solve(const SeldonTranspose&, Vector1& z);

    void Solve(const SeldonTranspose&, T* x_ptr, int nrhs);
    
  };

  template<class real, class cplx, class Storage, class Allocator>
  void qsplit_ilut(Vector<cplx, Storage, Allocator>& a, IVect& ind, int first,
                   int n, int ncut, const real& abs_ncut);

  template<class cplx, class Allocator1, class Allocator2>
  void GetIlut(const IlutPreconditioning<cplx, Allocator1>& param,
               Matrix<cplx, General, ArrayRowSparse, Allocator2>& A,
               IVect& iperm, IVect& rperm);

  template<class cplx, class Allocator>
  void GetIluk(int lfil, Matrix<cplx, General, ArrayRowSparse, Allocator>& A);

  template<class cplx, class Allocator>
  void GetIlu0(Matrix<cplx, General, ArrayRowSparse, Allocator>& A);

  template<class cplx, class Allocator>
  void GetMilu0(Matrix<cplx, General, ArrayRowSparse, Allocator>& A);

  template<class cplx, class Allocator1, class Allocator2>
  void GetIlut(const IlutPreconditioning<cplx, Allocator1>& param,
               Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator2>& A);
  
  template<class cplx, class Allocator>
  void GetIluk(int lfil, int print_level,
               Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A);
  
  template<class cplx, class Allocator>
  void GetIlu0(Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A);
  
  template<class cplx, class Allocator>
  void GetMilu0(Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A);
  
  template<class MatrixSparse, class T, class Alloc2>
  void GetLU(MatrixSparse& A, IlutPreconditioning<T, Alloc2>& mat_lu,
	     IVect& permut, bool keep_matrix, T& x);
  
  template<class MatrixSparse, class T, class Alloc2>
  void GetLU(MatrixSparse& A, IlutPreconditioning<T, Alloc2>& mat_lu,
	     IVect& permut, bool keep_matrix, complex<T>& x);
  
  template<class MatrixSparse, class T, class Alloc2>
  void GetLU(MatrixSparse& A, IlutPreconditioning<complex<T>, Alloc2>& mat_lu,
	     IVect& permut, bool keep_matrix, T& x);
  
  template<class T0, class Prop, class Storage, class Allocator,
	   class T, class Alloc2>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A,
	     IlutPreconditioning<T, Alloc2>& mat_lu,
	     IVect& permut, bool keep_matrix = false);

}

#define SELDON_FILE_ILUT_PRECONDITIONING_HXX
#endif
