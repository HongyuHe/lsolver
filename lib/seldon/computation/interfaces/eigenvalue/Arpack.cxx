// Copyright (C) 2011 Marc Duruflé
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

#ifndef SELDON_FILE_ARPACK_CXX

#include "Arpack.hxx"

namespace Seldon
{
  
  //! computation of eigenvalues and related eigenvectors
  /*!
    \param[in,out] var eigenproblem to solve
    \param[out] eigen_values eigenvalue
    \param[out] eigen_vectors eigenvectors
  */
#ifdef SELDON_WITH_VIRTUAL
  template<class T>
  void FindEigenvaluesArpack(EigenProblem_Base<T>& var,
                             Vector<T>& eigen_values,
                             Vector<T>& lambda_imag,
                             Matrix<T, General, ColMajor>& eigen_vectors)
#else
  template<class EigenProblem, class T, class Allocator1,
           class Allocator2, class Allocator3>
  void FindEigenvaluesArpack(EigenProblem& var,
                             Vector<T, VectFull, Allocator1>& eigen_values,
                             Vector<T, VectFull, Allocator2>& lambda_imag,
                             Matrix<T, General, ColMajor, Allocator3>& eigen_vectors)
#endif
  {
    // declaration (and initialization) of all the variables needed by CallArpack
    int ido = 0, n = var.GetM();
    char bmat;
    string which("LM");
    switch (var.GetTypeSorting())
      {
      case EigenProblem_Base<T>::SORTED_REAL : which = "LR"; break;
      case EigenProblem_Base<T>::SORTED_IMAG : which = "LI"; break;
      case EigenProblem_Base<T>::SORTED_MODULUS : which = "LM"; break;
      }
    
    if (var.GetTypeSpectrum() == var.SMALL_EIGENVALUES)
      {
        switch (var.GetTypeSorting())
          {
          case EigenProblem_Base<T>::SORTED_REAL : which = "SR"; break;
          case EigenProblem_Base<T>::SORTED_IMAG : which = "SI"; break;
          case EigenProblem_Base<T>::SORTED_MODULUS : which = "SM"; break;
          }
      }
    
    // initializing of computation
    int nev = var.GetNbAskedEigenvalues();
    double tol = var.GetStoppingCriterion();
    
    Vector<T> resid;
    int ncv = var.GetNbArnoldiVectors();
    Matrix<T, General, ColMajor> v;
    int ldv = n;
		
    IVect iparam(11), ipntr(14);
    iparam.Fill(0); ipntr.Fill(0);
    int ishift = 1, nb_max_iter = var.GetNbMaximumIterations();
    int sym = 0;
    if (var.IsSymmetricProblem())
      sym = 1;
    
    int non_sym = 0;
    int sym_mode = sym;
    iparam(0) = ishift;
    iparam(2) = nb_max_iter;
    iparam(3) = 1;
    
    int lworkl = 3*ncv*ncv + 6*ncv, info(0);
    Vector<T> workd, workl, Xh, Yh, Zh;
    Vector<double> rwork;

    T shiftr = var.GetShiftValue(), shifti = var.GetImagShiftValue();
    
#ifdef SELDON_WITH_VIRTUAL
    typename ClassComplexType<T>::Tcplx shift_complex, cone;
    SetComplexOne(cone);
    var.GetComplexShift(shiftr, shifti, shift_complex);
#endif    

    T zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
        
    // mass matrix for diagonal case, and Cholesky case
    int print_level = var.GetPrintLevel();

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var.GetCommunicator();
    int comm_f = MPI_Comm_c2f(comm);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int comm_f(0);
#endif

    // checking if computation mode is compatible with the spectrum
    // the used wants to find
    if (var.GetComputationalMode() == var.REGULAR_MODE)
      {
        if (var.GetTypeSpectrum() == var.CENTERED_EIGENVALUES)
          {
            cout << "You can not use regular mode to find "
                 << "eigenvalues closest to a given value" << endl;
            cout << "Try to use shifted mode for example "<<endl;
            abort();
          }
      }
    else
      {
        if ((var.GetTypeSpectrum() != var.CENTERED_EIGENVALUES) &&
	    (var.GetComputationalMode() != var.INVERT_MODE))
          {
            cout << "To find large or small eigenvalues, use a regular mode" << endl;
            abort();
          }
        
        if (var.GetComputationalMode() == var.IMAG_SHIFTED_MODE)
          {
            cout << "Currently not working" << endl;
            abort();
          }
	
        if ( (var.GetComputationalMode() == var.BUCKLING_MODE)
             || (var.GetComputationalMode() == var.CAYLEY_MODE) )
          {
            if ( !var.IsSymmetricProblem() || IsComplexNumber(zero))
              {
                cout << "Cayley or Bucking mode are reserved for real symmetric "
                     << "generalized eigenproblems " << endl;                
                abort();
              }
          }
        else          
          {
            if (shifti != zero)
              {
                if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
                  {
                    cout << "Complex shifts for unsymmetric real problems are "
                         << "not relevant for diagonal or Cholesky-factorized"
                         << " mass " << endl;                   
                    cout << "Specify a general mass matrix" << endl;
                    cout << "Diagonal = " << var.DiagonalMass() << 
                      ", Cholesky : " << var.UseCholeskyFactoForMass() << endl;
                    
                    abort();
                  }
              }
          }
      }
					 

    // we want to find eigenvalues of the generalized problem
    // K U  =  lambda M U
    
    // A first solution is to compute a cholesky factorization of M = L L^t
    // the system is transformed into
    // L^-1 K L^-T V = lambda V
    // with V = L^T U
    // we can use then regular or shifted mode for this standard problem
    // Similarly, if matrix M is diagonal, we have the standard problem 
    // D^-1/2 K D^-1/2 V = lambda V
    // with V = D^1/2 U
    
    // a second solution is to consider the generalized problem
    if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
      {
        // solving a standard eigenvalue problem
        if (print_level >= 2)
#ifdef SELDON_WITH_MPI
          if (rank_proc == 0)
#endif
            cout << "We solve a standard eigenvalue problem " << endl;
        
        // we can reduce to a standard problem
        // we consider S  = M^{-1/2}  K  M^{-1/2}
        // if matrix M is diagonal
        // and S = L^-1 K L^-T    if Cholesky factorisation of matrix
        // M = L L^T is used        
        bmat = 'I';
        if (var.DiagonalMass())
          {
            // computation of M
            var.ComputeDiagonalMass();
	    
            // computation of M^{-1/2}
            var.FactorizeDiagonalMass();
          }
        else 
          {
            // computation of M for Cholesky factorisation
            var.ComputeMassForCholesky();
            
            // computation of Cholesky factorisation M = L L^T
            var.FactorizeCholeskyMass();
          }

        if (var.GetComputationalMode() == var.REGULAR_MODE)
          {
            shiftr = zero;
            // regular mode -> mode 1 for arpack
            iparam(6) = 1;

            // computation of the stiffness matrix
            var.ComputeStiffnessMatrix();
					 
            // Initialization of variables
            v.Reallocate(n, ncv);
            workd.Reallocate(3*n);
            workl.Reallocate(lworkl);
            Xh.Reallocate(n);
            Yh.Reallocate(n);
            rwork.Reallocate(ncv);
            resid.Reallocate(n);
            resid.Fill(zero);
            v.Fill(zero);
            workd.Fill(zero);
            workl.Fill(zero);
            rwork.Fill(0.0);
            Xh.Fill(zero);
            Yh.Fill(zero);
            
            if (print_level >= 2)
#ifdef SELDON_WITH_MPI
              if (rank_proc == 0)
#endif
                cout << "Starting Arpack iterations..." << endl;
            
            bool test_loop = true;
            while (test_loop)
              {
                CallArpack(comm_f, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                           iparam, ipntr, sym, workd, workl, lworkl, rwork, info);

                if ((ido == -1)||(ido == 1))
                  {
                    // matrix vector product with M^-1/2  K  M^-1/2
                    //  or L^-1  K L^-T 
                    for (int i = 0; i < n; i++)
                      Xh(i) = workd(ipntr(0)+i-1);

                    if (var.DiagonalMass())
                      var.MltInvSqrtDiagonalMass(Xh);
                    else
                      var.SolveCholeskyMass(SeldonTrans, Xh);

		    var.MltStiffness(Xh, Yh);
                    var.IncrementProdMatVect();

		    if (var.DiagonalMass())
                      var.MltInvSqrtDiagonalMass(Yh);
                    else
                      var.SolveCholeskyMass(SeldonNoTrans, Yh);

		    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                    
                  }
                else
                  test_loop = false;

              }
          }
        else
          {
            // shifted mode, we compute a factorization of (K - \sigma M)
            // where \sigma is the shift
            // shifted mode for regular problem -> mode 1 for arpack
            iparam(6) = 1;
            
            // computation and factorization of K - sigma M
            var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);
            
            // Initialization of variables
            v.Reallocate(n, ncv);
            workd.Reallocate(3*n);
            workl.Reallocate(lworkl);
            Xh.Reallocate(n);
            Yh.Reallocate(n);
            rwork.Reallocate(ncv);
            resid.Reallocate(n);
            resid.Fill(zero);
            v.Fill(zero);
            workd.Fill(zero);
            workl.Fill(zero);
            rwork.Fill(0.0);
            Xh.Fill(zero);
            Yh.Fill(zero);
		
            if (print_level >= 2)
#ifdef SELDON_WITH_MPI
              if (rank_proc == 0)
#endif
                cout << "Starting Arpack iterations..." << endl;
            
            bool test_loop = true;			 
            while (test_loop)
              {
                CallArpack(comm_f, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                           iparam, ipntr, sym, workd, workl, lworkl, rwork, info);

                if ((ido == -1)||(ido == 1))
                  {
                    // matrix vector product M^1/2 (K - sigma M)^-1 M^1/2  X
                    // or L^T (K - sigma M)^-1 L  X
                    for (int i = 0; i < n; i++)
                      Xh(i) = workd(ipntr(0)+i-1);

                    if (var.DiagonalMass())
                      var.MltSqrtDiagonalMass(Xh);
                    else
                      var.MltCholeskyMass(SeldonNoTrans, Xh);

                    var.ComputeSolution(Xh, Yh);
                    var.IncrementProdMatVect();
                    
                    if (var.DiagonalMass())
                      var.MltSqrtDiagonalMass(Yh);
                    else
                      var.MltCholeskyMass(SeldonTrans, Yh);

                    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                  }
                else
                  test_loop = false;
                
              }
          }
      }
    else
      {
        // generalized problem
        bmat = 'G';
        // different computational modes available in Arpack :
        if (var.GetComputationalMode() == var.INVERT_MODE)
          {
            // we consider standard problem M^-1 K U = lambda U
            // drawback : the matrix is non-symmetric
	    //           even if K and M are symmetric            
            bmat = 'I';
            
            if (var.GetTypeSpectrum() != var.CENTERED_EIGENVALUES)
              {
                // => regular mode for matrix M^-1 K
                iparam(6) = 1;
                
                // computation and factorisation of mass matrix
                var.ComputeAndFactorizeStiffnessMatrix(one, zero);
                
                // computation of stiffness matrix
                var.ComputeStiffnessMatrix();
              }
            else
              {
                if (shifti != zero)
                  {
                    // shifted mode for generalized eigenvalue problem -> mode 3 for arpack
                    iparam(6) = 3;
                    
                    // computation and factorization of (K - sigma M)^-1
#ifdef SELDON_WITH_VIRTUAL
                    var.ComputeAndFactorizeStiffnessMatrix(-shift_complex, cone,
                                                           EigenProblem_Base<T>::REAL_PART);
#else
                    var.ComputeAndFactorizeStiffnessMatrix(-complex<T>(shiftr, shifti),
                                                           complex<T>(one, zero), true);
#endif
                    
                    var.ComputeStiffnessMatrix();
                  }
                else
                  {
                    // => regular mode for matrix (M^-1 K - sigma I)^-1
                    // = (K - sigma M)^-1 M
                    iparam(6) = 1;
                    
                    // computation and factorization of (K - sigma M)
                    var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);
                    
                  }
                
                // computation of M
                var.ComputeMassMatrix();
              }
					 
            // Initialization of variables
            v.Reallocate(n, ncv);
            workd.Reallocate(3*n);
            workl.Reallocate(lworkl);
            Xh.Reallocate(n);
            Yh.Reallocate(n);
            Zh.Reallocate(n);
            rwork.Reallocate(ncv);
            resid.Reallocate(n);
            resid.Fill(zero);
            v.Fill(zero);
            workd.Fill(zero);
            workl.Fill(zero);
            rwork.Fill(0.0);
            Xh.Fill(zero);
            Yh.Fill(zero);
            Zh.Fill(zero);
            
            if (print_level >= 2)
#ifdef SELDON_WITH_MPI
              if (rank_proc == 0)
#endif
                cout << "Starting Arpack iterations..." << endl;

            bool test_loop = true;		
            sym_mode = 0;
            while (test_loop)
              {
                CallArpack(comm_f, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                           iparam, ipntr, non_sym, workd, workl, lworkl, rwork, info);
                
                if ((ido == -1)||(ido == 1))
                  {
                    // computation of (K - sigma M)^-1 M X
                    if ((iparam(6) == 3) && (ido == 1))
                      for (int i = 0; i < n; i++)
                        Xh(i) = workd(ipntr(2)+i-1);
                    else
                      for (int i = 0; i < n; i++)
                        Xh(i) = workd(ipntr(0)+i-1);
                    
                    if (var.GetTypeSpectrum() != var.CENTERED_EIGENVALUES)
                      {
			var.MltStiffness(Xh, Zh);
                        var.ComputeSolution(Zh, Yh);
                        var.IncrementProdMatVect();
                      }
                    else
                      {
                        var.MltMass(Xh, Zh);
                        var.ComputeSolution(Zh, Yh);
                        var.IncrementProdMatVect();
                      }
			
                    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                  }
                else if (ido == 2)
                  {
                    // computation of M X (Identity)
                    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = workd(ipntr(0)+i-1);
                  }
                else
                  test_loop = false;		
              }	    
          }
        else if (var.GetComputationalMode() == var.REGULAR_MODE)
          {
            // regular mode for generalized eigenvalue problem -> mode 2 for arpack
            iparam(6) = 2;
            bmat = 'G';
            
            // factorization of the mass matrix
            var.ComputeAndFactorizeStiffnessMatrix(one, zero);
					 
            // computation of stiffness and mass matrix
            var.ComputeStiffnessMatrix();
            var.ComputeMassMatrix();
            					 
            // Initialization of variables
            v.Reallocate(n, ncv);
            workd.Reallocate(3*n);
            workl.Reallocate(lworkl);
            Xh.Reallocate(n);
            Yh.Reallocate(n);
            Zh.Reallocate(n);
            rwork.Reallocate(ncv);
            resid.Reallocate(n);
            resid.Fill(zero);
            v.Fill(zero);
            workd.Fill(zero);
            workl.Fill(zero);
            rwork.Fill(0.0);
            Xh.Fill(zero);
            Yh.Fill(zero);
            Zh.Fill(zero);

            if (print_level >= 2)
#ifdef SELDON_WITH_MPI
              if (rank_proc == 0)
#endif
                cout << "Starting Arpack iterations..." << endl;

            bool test_loop = true;
            while (test_loop)
              {
                CallArpack(comm_f, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                           iparam, ipntr, sym, workd, workl, lworkl, rwork, info);
                
                if ((ido == -1)||(ido == 1))
                  {
                    // computation of M^-1 K X and K X
                    for (int i = 0; i < n; i++)
                      Xh(i) = workd(ipntr(0)+i-1);
                    
                    var.MltStiffness(Xh, Zh);
                    var.ComputeSolution(Zh, Yh);
                    var.IncrementProdMatVect();
							 
                    for (int i = 0; i < n; i++)
                      {
                        workd(ipntr(0)+i-1) = Zh(i);
                        workd(ipntr(1)+i-1) = Yh(i);
                      }
                  }
                else if (ido == 2)
                  {
                    // computation of M X
                    for (int i = 0; i < n; i++)
                      Xh(i) = workd(ipntr(0)+i-1);
                    
                    var.MltMass(Xh, Yh);
							 
                    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                  }
                else
                  test_loop = false;
						 
              }
          }
        else if ((var.GetComputationalMode() == var.SHIFTED_MODE)
                 || (var.GetComputationalMode() == var.IMAG_SHIFTED_MODE) )
          {
            // shifted mode for generalized eigenvalue problem -> mode 3 for arpack
            iparam(6) = 3;
            
            // computation and factorization of (K - sigma M)^-1
            if (shifti != zero)
              {
#ifdef SELDON_WITH_VIRTUAL
                if (var.GetComputationalMode() == var.SHIFTED_MODE)
                  var.ComputeAndFactorizeStiffnessMatrix(-shift_complex, cone,
							 EigenProblem_Base<T>::REAL_PART);
                else
                  {
                    iparam(6) = 4;
                    var.ComputeAndFactorizeStiffnessMatrix(-shift_complex, cone,
							   EigenProblem_Base<T>::IMAG_PART);
                  }
#else
                if (var.GetComputationalMode() == var.SHIFTED_MODE)
                  var.ComputeAndFactorizeStiffnessMatrix(-complex<T>(shiftr, shifti),
                                                         complex<T>(one, zero), true);
                else
                  {
                    iparam(6) = 4;
                    var.ComputeAndFactorizeStiffnessMatrix(-complex<T>(shiftr, shifti),
                                                           complex<T>(one, zero), false);
                  }
#endif
                
                var.ComputeStiffnessMatrix();
              }
            else
              var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);
            
            // computation of mass matrix
            var.ComputeMassMatrix();
            
            // Initialization of variables
            v.Reallocate(n, ncv);
            workd.Reallocate(3*n);
            workl.Reallocate(lworkl);
            Xh.Reallocate(n);
            Yh.Reallocate(n);
            Zh.Reallocate(n);
            rwork.Reallocate(ncv);
            resid.Reallocate(n);
            resid.Fill(zero);
            v.Fill(zero);
            workd.Fill(zero);
            workl.Fill(zero);
            rwork.Fill(0.0);
            Xh.Fill(zero);
            Yh.Fill(zero);
            Zh.Fill(zero);
		
            if (print_level >= 2)
#ifdef SELDON_WITH_MPI
              if (rank_proc == 0)
#endif
                cout << "Starting Arpack iterations..." << endl;

            bool test_loop = true;
            while (test_loop)
              {
                CallArpack(comm_f, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                           iparam, ipntr, sym, workd, workl, lworkl, rwork, info);
               
                if ((ido == -1)||(ido == 1))
                  {
                    // computation of Real( (K - sigma M)^-1  M) X
                    // or Imag( (K - sigma M)^-1  M) X
                    // if ido == 1, M X is already computed and available               
                    if (ido == -1)
                      {
                        for (int i = 0; i < n; i++)
                          Zh(i) = workd(ipntr(0)+i-1);
			
			var.MltMass(Zh, Xh);
		      }
                    else
                      for (int i = 0; i < n; i++)
                        Xh(i) = workd(ipntr(2)+i-1);
                    
                    var.ComputeSolution(Xh, Yh);
                    var.IncrementProdMatVect();
		    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                    
                  }
                else if (ido == 2)
                  {
                    // computation of M X
                    for (int i = 0; i < n; i++)
                      Xh(i) = workd(ipntr(0)+i-1);
		    
		    var.MltMass(Xh, Yh);
		    
                    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
		    
		  }
                else
                  test_loop = false;
		
              }					 
          }
        else if (var.GetComputationalMode() == var.BUCKLING_MODE)
          {
            // buckling mode for generalized eigenvalue problem -> mode 4 for arpack
            iparam(6) = 4;
            
            // computation and factorisation of (K - sigma M)
            var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);
            
            // computation of matrix K
            var.ComputeStiffnessMatrix();
            					 
            // Initialization of variables
            v.Reallocate(n, ncv);
            workd.Reallocate(3*n);
            workl.Reallocate(lworkl);
            Xh.Reallocate(n);
            Yh.Reallocate(n);
            Zh.Reallocate(n);
            rwork.Reallocate(ncv);
            resid.Reallocate(n);
            resid.Fill(zero);
            v.Fill(zero);
            workd.Fill(zero);
            workl.Fill(zero);
            rwork.Fill(0.0);
            Xh.Fill(zero);
            Yh.Fill(zero);
            Zh.Fill(zero);
            
            if (print_level >= 2)
#ifdef SELDON_WITH_MPI
              if (rank_proc == 0)
#endif
                cout << "Starting Arpack iterations..." << endl;

            bool test_loop = true;			 
            while (test_loop)
              {
                CallArpack(comm_f, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                           iparam, ipntr, sym, workd, workl, lworkl, rwork, info);
		
                if ((ido == -1)||(ido == 1))
                  {
                    // computation of (K - sigma M)^-1 K X or (K - sigma M)^-1 X
                    if (ido == -1)
                      {
                        for (int i = 0; i < n; i++)
                          Zh(i) = workd(ipntr(0)+i-1);
								 
                        var.MltStiffness(Zh, Xh);
                      }
                    else
                      for (int i = 0; i < n; i++)
                        Xh(i) = workd(ipntr(2)+i-1);
							 
                    var.ComputeSolution(Xh, Yh);
                    var.IncrementProdMatVect();
		    
		    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                  }
                else if (ido == 2)
                  {
                    // computation of K X
                    for (int i = 0; i < n; i++)
                      Xh(i) = workd(ipntr(0)+i-1);
		    
                    var.MltStiffness(Xh, Yh);
                    var.IncrementProdMatVect();
		    
		    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                  }
                else
                  test_loop = false;
		
              }
          }
        else if (var.GetComputationalMode() == var.CAYLEY_MODE)
          {
            // shifted mode for generalized eigenvalue problem -> mode 5 for arpack
            iparam(6) = 5;
            
            // computation and factorization of (K - sigma M)
            var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);
            
            // computation of M and (K + sigma M)
            var.ComputeMassMatrix();
            var.ComputeStiffnessMatrix(shiftr, one);
            					 
            // Initialization of variables
            v.Reallocate(n, ncv);
            workd.Reallocate(3*n);
            workl.Reallocate(lworkl);
            Xh.Reallocate(n);
            Yh.Reallocate(n);
            Zh.Reallocate(n);
            rwork.Reallocate(ncv);
            resid.Reallocate(n);
            resid.Fill(zero);
            v.Fill(zero);
            workd.Fill(zero);
            workl.Fill(zero);
            rwork.Fill(0.0);
            Xh.Fill(zero);
            Yh.Fill(zero);
            Zh.Fill(zero);
            
            if (print_level >= 2)
#ifdef SELDON_WITH_MPI
              if (rank_proc == 0)
#endif
                cout << "Starting Arpack iterations..." << endl;

            bool test_loop = true;			 
            while (test_loop)
              {
                CallArpack(comm_f, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                           iparam, ipntr, sym, workd, workl, lworkl, rwork, info);
                
                if ((ido == -1)||(ido == 1))
                  {
                    // computation of (K - sigma M)^-1 (K + sigma M) X
                    for (int i = 0; i < n; i++)
                      Zh(i) = workd(ipntr(0)+i-1);
							 
                    var.MltStiffness(shiftr, one, Zh, Xh);
                    var.ComputeSolution(Xh, Yh);
                    var.IncrementProdMatVect();
							 
                    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                  }
                else if (ido == 2)
                  {
                    // computation of M X
                    for (int i = 0; i < n; i++)
                      Xh(i) = workd(ipntr(0)+i-1);
							 
                    var.MltMass(Xh, Yh);
							 
                    for (int i = 0; i < n; i++)
                      workd(ipntr(1)+i-1) = Yh(i);
                  }
                else
                  test_loop = false;
						 
              }
          }	
      }    

    if (info != 0)
      {
        if (info == 1)
          {
            cout << "Maximum number of iterations reached" << endl;
            cout << "Try again with a larger number of iterations" 
                 << " or with a less restrictive stopping criterion" << endl;
            
            throw SolverMaximumIterationError("FindEigenvaluesArpack", "Maximum number of iterations reached");
          }
        else
          {
            cout << "Error during the computation of eigenvalues, info = "
		 << info << endl;
            cout << "Look at documentation of ARPACK " << endl;
            
            throw SolverDivergenceError("FindEigenvaluesArpack", "The solver diverged");        
          }
      }
    
    Xh.Clear();
    Yh.Clear();
    Zh.Clear();
    
    if ((Norm2(resid) == 0) && (n > ncv))
      {
	cout << "This eigenvalue problem rises an unknown bug in Arpack " << endl;
	cout << "Is the mass matrix or stiffness matrix null ?" << endl;
	cout << "Null eigenvalues are found." << endl;
        throw SolverDivergenceError("FindEigenvaluesArpack", "The solver diverged");        
      }
    
    // we recover eigenvalues and eigenvectors
    int nconv = nev+1+var.GetNbAdditionalEigenvalues();
    eigen_values.Reallocate(nconv);
    eigen_vectors.Reallocate(n, nconv);
    eigen_values.Fill(zero);
    eigen_vectors.Fill(zero);
    lambda_imag.Reallocate(nconv);
    lambda_imag.Fill(zero);
    
    char howmny = 'A';
    IVect selec(ncv); selec.Zero();
    int rvec = 1;
    int ldz = n;
    if (print_level >= 2)
#ifdef SELDON_WITH_MPI
      if (rank_proc == 0)
#endif
        cout << "recovering eigenvalues ..." << endl;

    CallArpack(comm_f, rvec, howmny, selec, eigen_values, lambda_imag, eigen_vectors,
               ldz, shiftr, shifti, bmat, n, which, nev, tol, resid, ncv,
               v, ldv, iparam, ipntr, sym_mode, workd, workl, lworkl, rwork, info);

    if (print_level >= 2)
#ifdef SELDON_WITH_MPI
      if (rank_proc == 0)
#endif
        cout << "end recover " << endl;

    eigen_values.Resize(nev);
    lambda_imag.Resize(nev);
    eigen_vectors.Resize(n, nev);
    ApplyScalingEigenvec(var, eigen_values, lambda_imag, eigen_vectors,
                         shiftr, shifti);

    var.Clear();
  }


  /*************************************************************************
   * Overload of CallArpack function to call the correct Arpack subroutine *
   *************************************************************************/
  
  
  //! calling arpack routine
  template<class T, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Alloc5>
  void CallArpack(int comm, int& ido, char& bmat, int& n, string& which, int& nev,
		  T& tol, Vector<T, VectFull, Allocator1>& resid,
		  int& ncv, Matrix<T, General, ColMajor, Allocator2>& v,
		  int& ldv, Vector<int>& iparam, Vector<int>& ipntr, int sym,
		  Vector<T, VectFull, Allocator3>& workd,
		  Vector<T, VectFull, Allocator4>& workl,
		  int& lworkl, Vector<T, VectFull, Alloc5>& rwork, int& info)
  {
    if (sym)
      {
        if ((which == "SR") || (which == "SI"))
          which = "SA";
        
        if ((which == "LR") || (which == "LI"))
          which = "LA";
        
        // real symmetric
        // call of dsaupd/ssaupd
#ifdef SELDON_WITH_MPI
        saupd(comm, ido, bmat, n, &which[0], nev, tol, resid.GetData(),
	      ncv, v.GetData(), ldv, iparam.GetData(), ipntr.GetData(), 
	      workd.GetData(), workl.GetData(), lworkl, info);
#else
        saupd(ido, bmat, n, &which[0], nev, tol, resid.GetData(),
	      ncv, v.GetData(), ldv, iparam.GetData(), ipntr.GetData(), 
	      workd.GetData(), workl.GetData(), lworkl, info);
#endif
      }
    else
      {
        // real non-symmetric
        // call of dnaupd/snaupd
#ifdef SELDON_WITH_MPI
        naupd(comm, ido, bmat, n, &which[0], nev, tol, resid.GetData(),
	      ncv, v.GetData(), ldv, iparam.GetData(), ipntr.GetData(), 
	      workd.GetData(), workl.GetData(), lworkl, info);
#else        
        naupd(ido, bmat, n, &which[0], nev, tol, resid.GetData(),
	      ncv, v.GetData(), ldv, iparam.GetData(), ipntr.GetData(), 
	      workd.GetData(), workl.GetData(), lworkl, info);
#endif
      }
  }
  
  
  //! calling arpack routine
  template<class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Alloc5>
  void CallArpack(int comm, int& ido, char& bmat, int& n,
		  string& which, int& nev, double& tol,
		  Vector<complex<double>, VectFull, Allocator1>& resid, int& ncv,
		  Matrix<complex<double>, General, ColMajor, Allocator2>& v,
		  int& ldv, Vector<int>& iparam, Vector<int>& ipntr, int sym,
		  Vector<complex<double>, VectFull, Allocator3>& workd,
		  Vector<complex<double>, VectFull, Allocator4>& workl,
		  int& lworkl, Vector<double, VectFull, Alloc5>& rwork, int& info)
  {
    // complex
    // call of znaupd
#ifdef SELDON_WITH_MPI
    pznaupd_(&comm, &ido, &bmat, &n, &which[0], &nev, &tol, resid.GetDataVoid(),
            &ncv, v.GetDataVoid(), &ldv, iparam.GetData(), ipntr.GetData(), 
            workd.GetDataVoid(), workl.GetDataVoid(), &lworkl, rwork.GetData(), &info);
#else
    znaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid.GetDataVoid(),
            &ncv, v.GetDataVoid(), &ldv, iparam.GetData(), ipntr.GetData(), 
            workd.GetDataVoid(), workl.GetDataVoid(), &lworkl, rwork.GetData(), &info);
#endif
  }
  
  
  //! calling Arpack routine
  template<class T, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Allocator5,
	   class Allocator6, class Allocator7, class Alloc8>
  void CallArpack(int comm, int& rvec, char& howmny, Vector<int>& selec,
		  Vector<T, VectFull, Allocator1>& lambda,
		  Vector<T, VectFull, Allocator2>& lambda_i,
		  Matrix<T, General, ColMajor, Allocator3>& eigen_vec,
		  int& ldz, T& shiftr, T& shifti, char& bmat, int& n,
		  string& which, int& nev, T& tol,
		  Vector<T, VectFull, Allocator4>& resid,
		  int& ncv, Matrix<T, General, ColMajor, Allocator5>&v,
		  int& ldv, Vector<int>& iparam, Vector<int>& ipntr,
		  int sym, Vector<T, VectFull, Allocator6>& workd,
		  Vector<T, VectFull, Allocator7>& workl,
		  int& lworkl, Vector<T, VectFull, Alloc8>& rwork, int& info)
  {
    if (sym)
      {
        // real symmetric
        // call of dseupd
#ifdef SELDON_WITH_MPI
        seupd(comm, rvec, howmny, selec.GetData(), lambda.GetData(), eigen_vec.GetData(),
	      ldz, shiftr, bmat, n, &which[0], nev, tol, resid.GetData(), ncv,
	      v.GetData(), ldv, iparam.GetData(), ipntr.GetData(),
	      workd.GetData(), workl.GetData(), lworkl, info);
#else
        seupd(rvec, howmny, selec.GetData(), lambda.GetData(), eigen_vec.GetData(),
	      ldz, shiftr, bmat, n, &which[0], nev, tol, resid.GetData(), ncv,
	      v.GetData(), ldv, iparam.GetData(), ipntr.GetData(),
	      workd.GetData(), workl.GetData(), lworkl, info);
#endif
      }
    else
      {
        Vector<double> workev(3*ncv);
        // real non-symmetric
        // call of dneupd
#ifdef SELDON_WITH_MPI
        neupd(comm, rvec, howmny, selec.GetData(), lambda.GetData(),
	      lambda_i.GetData(), eigen_vec.GetData(), ldz, shiftr, shifti,
	      workev.GetData(), bmat, n, &which[0], nev, tol,
	      resid.GetData(), ncv, v.GetData(), ldv, iparam.GetData(),
	      ipntr.GetData(), workd.GetData(), workl.GetData(), lworkl, info);
#else
        neupd(rvec, howmny, selec.GetData(), lambda.GetData(),
              lambda_i.GetData(), eigen_vec.GetData(), ldz, shiftr, shifti,
	      workev.GetData(), bmat, n, &which[0], nev, tol,
	      resid.GetData(), ncv, v.GetData(), ldv, iparam.GetData(),
	      ipntr.GetData(), workd.GetData(), workl.GetData(), lworkl, info);
#endif
      }
  }
  
  
  //! calling arpack routine
  template<class Allocator1, class Allocator2, class Allocator3,
	   class Allocator4, class Allocator5, class Allocator6,
	   class Allocator7, class Alloc8>
  void CallArpack(int comm, int& rvec, char& howmny, Vector<int>& selec,
		  Vector<complex<double>, VectFull, Allocator1>& lambda,
		  Vector<complex<double>, VectFull, Allocator2>& lambda_i,
		  Matrix<complex<double>, General, ColMajor, Allocator3>& eigen_vectors,
		  int& ldz, complex<double>& shiftr, complex<double>& shifti,
		  char& bmat, int& n, string& which, int& nev, double& tol,
		  Vector<complex<double>, VectFull, Allocator4>& resid,
		  int& ncv, Matrix<complex<double>, General, ColMajor, Allocator5>& v,
		  int& ldv, Vector<int>& iparam, Vector<int>& ipntr, int sym,
		  Vector<complex<double>, VectFull, Allocator6>& workd,
		  Vector<complex<double>, VectFull, Allocator7>& workl,
		  int& lworkl, Vector<double, VectFull, Alloc8>& rwork, int& info)
  {
    // complex
    // call of zneupd
    Vector<complex<double> > workev(2*ncv);
    workev.Zero();
#ifdef SELDON_WITH_MPI
    pzneupd_(&comm, &rvec, &howmny, selec.GetData(), lambda.GetDataVoid(),
             eigen_vectors.GetDataVoid(), &ldz, &shiftr, workev.GetDataVoid(),
             &bmat, &n, &which[0], &nev, &tol, resid.GetDataVoid(),
             &ncv, v.GetDataVoid(), &ldv, iparam.GetData(),
             ipntr.GetData(), workd.GetDataVoid(), workl.GetDataVoid(),
             &lworkl, rwork.GetData(), &info);
#else
    zneupd_(&rvec, &howmny, selec.GetData(), lambda.GetDataVoid(),
	    eigen_vectors.GetDataVoid(), &ldz, &shiftr, workev.GetDataVoid(),
	    &bmat, &n, &which[0], &nev, &tol, resid.GetDataVoid(),
	    &ncv, v.GetDataVoid(), &ldv, iparam.GetData(),
	    ipntr.GetData(), workd.GetDataVoid(), workl.GetDataVoid(),
	    &lworkl, rwork.GetData(), &info);
#endif
  }
  
}
			 
#define SELDON_FILE_ARPACK_CXX
#endif
