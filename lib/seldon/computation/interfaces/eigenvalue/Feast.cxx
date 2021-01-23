#ifndef SELDON_FILE_FEAST_CXX
#define SELDON_FILE_FEAST_CXX

namespace Seldon
{
  
  // main function to find eigenvalues and eigenvectors with Feast (MKL implementation)
#ifdef SELDON_WITH_VIRTUAL
  template<class T>
  void FindEigenvaluesFeast(EigenProblem_Base<T>& var,
                            Vector<T>& eigen_values,
                            Vector<T>& lambda_imag,
                            Matrix<T, General, ColMajor>& eigen_vectors)
#else
  template<class EigenProblem, class T, class Allocator1,
	   class Allocator2, class Allocator3>
  void FindEigenvaluesFeast(EigenProblem& var,
                            Vector<T, VectFull, Allocator1>& eigen_values,
                            Vector<T, VectFull, Allocator2>& lambda_imag,
                            Matrix<T, General, ColMajor, Allocator3>& eigen_vectors)
#endif
  {    
    // initialization of feast parameters
    IVect fpm(128);
    fpm.Fill(0);
    feastinit_(fpm.GetData());
    
    fpm(0) = 0;
    if (var.GetGlobalPrintLevel() > 0)
      fpm(0) = 1;

    double tol = var.GetStoppingCriterion();    
    fpm(2) = -log10(tol);
    
    for (int i = 64; i < fpm.GetM(); i++)
      fpm(i) = 0;
    
    fpm(3) = var.GetNbMaximumIterations();

    FeastParam& param = var.GetFeastParameters();
    if (param.EstimateNumberEigenval())
      fpm(13) = 2;

    // regular mode only
    // no shift is applied, the eigenvalues are directly searched into a circle
    var.SetComputationalMode(var.REGULAR_MODE);

#ifdef SELDON_WITH_MPI
    fpm(8) = MPI_Comm_c2f(var.GetGlobalCommunicator());
#endif
    
    
    // initialization of computation
    int n = var.GetM();

    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;
        
    int m0 = var.GetNbAskedEigenvalues()+1;

    Treal emin = var.GetLowerBoundInterval();
    Treal emax = var.GetUpperBoundInterval();
    Tcplx Emid = var.GetCircleCenterSpectrum();
    Treal r = var.GetCircleRadiusSpectrum();
    fpm(17) = var.GetRatioEllipseSpectrum();
    fpm(18) = var.GetAngleEllipseSpectrum();
    
    T zero; Tcplx one;
    SetComplexZero(zero);
    SetComplexOne(one);

    bool herm = var.IsHermitianProblem();
    bool sym = var.IsSymmetricProblem();
    if (!herm)
      {
        if (var.UseCholeskyFactoForMass())
          {
            cout << "Cholesky case not implemented in Feast interface" << endl;
            abort();
          }
      }

    if (herm)
      {
        if (param.GetTypeIntegration() >= 0)
          fpm(15) = param.GetTypeIntegration();
        
        if (param.GetNumOfQuadraturePoints() > 0)
          fpm(1) = param.GetNumOfQuadraturePoints();
      }
    else
      {
        if (param.GetTypeIntegration() >= 0)
          fpm(16) = param.GetTypeIntegration();

        if (param.GetNumOfQuadraturePoints() > 0)
          fpm(7) = param.GetNumOfQuadraturePoints();
      }

    // work arrays
    Matrix<T, General, ColMajor> work;
    Matrix<Tcplx, General, ColMajor> zwork;
    Matrix<Tcplx, General, ColMajor> workc;
    Vector<Tcplx> xc(n), yc(n);
    Vector<T> aq, bq;
    Vector<Tcplx> zaq, zbq;
    Vector<T> x(n), y(n);
    
    workc.Reallocate(n, m0);
    if (herm)
      {
        work.Reallocate(n, m0);        
        aq.Reallocate(m0*m0);
        bq.Reallocate(m0*m0);
      }
    else if (sym)
      {
        zwork.Reallocate(n, m0);        
        zaq.Reallocate(m0*m0);
        zbq.Reallocate(m0*m0);
      }
    else
      {
        zwork.Reallocate(n, 2*m0);
        zaq.Reallocate(m0*m0);
        zbq.Reallocate(m0*m0);
      }
    
    work.Zero(); zwork.Zero();
    workc.Zero(); xc.Zero(); yc.Zero();
    x.Zero(); y.Zero();
    aq.Zero(); bq.Zero(); zaq.Zero(); zbq.Zero();
    
    // output arrays
    Treal epsout(0); Tcplx ze;
    Vector<Treal> res(m0);

    int loop = 0, m = 0;
    Matrix<Tcplx, General, ColMajor> eigen_vec_cplx;
    Vector<Treal> lambda; Vector<Tcplx> lambda_cplx;

    if (herm)
      {
        eigen_vectors.Reallocate(n, m0);
        lambda.Reallocate(m0);
      }
    else if (sym)
      {
        eigen_vec_cplx.Reallocate(n, m0);
        lambda_cplx.Reallocate(m0);
      }
    else
      {
        eigen_vec_cplx.Reallocate(n, 2*m0);
        lambda_cplx.Reallocate(m0);
        res.Reallocate(2*m0);
      }

    res.Zero();
    eigen_vectors.Zero(); lambda.Zero();
    lambda_cplx.Zero(); eigen_vec_cplx.Zero();
    
    if (var.DiagonalMass())
      {
        // if the mass matrix is diagonal :
        // the diagonal is computed
        var.ComputeDiagonalMass();
        
        // and M^{-1/2} is computed
        var.FactorizeDiagonalMass();
      }
    else if (var.UseCholeskyFactoForMass())
      {
        // if the user wants to use the Cholesky factorization of the mass matrix:
        
        // computing mass matrix in a convenient form
        var.ComputeMassForCholesky();
        
        // performing Cholesky factorization
        var.FactorizeCholeskyMass();
      }
    else
      var.ComputeMassMatrix();
    
    // evaluating stiffness matrix K
    var.ComputeStiffnessMatrix();

    // main loop (reverse communication interface)
    int ijob = -1, info = 0;
    while (ijob != 0)
      { 
        // feast is called
        CallFeast(ijob, n, herm, sym, ze, work, workc, aq, bq, zwork, zaq, zbq,
                  fpm, epsout, loop, emin, emax, Emid, r, m0, lambda, eigen_vectors,
                  lambda_cplx, eigen_vec_cplx, m, res, info);
        
        if (ijob == 10)
          {
            // we have to factorize ze M - K
            var.ComputeAndFactorizeStiffnessMatrix(ze, -one);
          }
        else if (ijob == 11)
          {
            // solves (ze K - M) y = workc(n, m0)
            for (int k = 0; k < fpm(22); k++)
              {
                for (int i = 0; i < n; i++)
                  xc(i) = workc(i, k);
                
                if (var.DiagonalMass())
                  var.MltSqrtDiagonalMass(xc);
                else if (var.UseCholeskyFactoForMass())
                  var.MltCholeskyMass(SeldonNoTrans, xc);
                
                var.ComputeSolution(xc, yc);
                var.IncrementProdMatVect();
                
                if (var.DiagonalMass())
                  var.MltSqrtDiagonalMass(yc);
                else if (var.UseCholeskyFactoForMass())
                  var.MltCholeskyMass(SeldonTrans, yc);
                
                for (int i = 0; i < n; i++)
                  workc(i, k) = yc(i);
              }
          }
        else if (ijob == 20)
          {
            // factorize (ze K - M)^H
            // already done in ijob = 10
          }
        else if (ijob == 21)
          {
            // solves (ze K - M)^H y = workc(n, m0)
            for (int k = 0; k < fpm(22); k++)
              {
                for (int i = 0; i < n; i++)
                  xc(i) = workc(i, k);
                
                if (var.DiagonalMass())
                  var.MltSqrtDiagonalMass(xc);
                else if (var.UseCholeskyFactoForMass())
                  var.MltCholeskyMass(SeldonNoTrans, xc);
                
                Conjugate(xc);
                var.ComputeSolution(SeldonTrans, xc, yc);
                Conjugate(yc);
                
                if (var.DiagonalMass())
                  var.MltSqrtDiagonalMass(yc);
                else if (var.UseCholeskyFactoForMass())
                  var.MltCholeskyMass(SeldonTrans, yc);
                
                for (int i = 0; i < n; i++)
                  workc(i, k) = yc(i);
              }
          }
        else if (ijob == 30)
          {
            // multiplication by matrix A
            int i = fpm(23), j = fpm(23) + fpm(24)-1;
            if (herm)
              for (int k = i; k <= j; k++)
                {
                  GetCol(eigen_vectors, k-1, x);
                  
                  if (var.DiagonalMass())
                    var.MltInvSqrtDiagonalMass(x);
                  else if (var.UseCholeskyFactoForMass())
                    var.SolveCholeskyMass(SeldonTrans, x);
                  
                  var.MltStiffness(x, y);
                  
                  if (var.DiagonalMass())
                    var.MltInvSqrtDiagonalMass(y);
                  else if (var.UseCholeskyFactoForMass())
                    var.SolveCholeskyMass(SeldonNoTrans, y);
                  
                  for (int p = 0; p < n; p++)
                    work(p, k-1) = y(p);
                }
            else
              for (int k = i; k <= j; k++)
                {
                  GetCol(eigen_vec_cplx, k-1, xc);
                  
                  if (var.DiagonalMass())
                    var.MltInvSqrtDiagonalMass(xc);

                  var.MltStiffness(xc, yc);
                  
                  if (var.DiagonalMass())
                    var.MltInvSqrtDiagonalMass(yc);
                  
                  SetCol(yc, k-1, zwork);
                }
          }
        else if (ijob == 31)
          {
            // multiplication by matrix A^H
            int i = fpm(33), j = fpm(33) + fpm(34)-1;
            if (herm)
              {
                cout << "impossible"  << endl;
                abort();
              }
            else
              for (int k = i; k <= j; k++)
                {
                  GetCol(eigen_vec_cplx, k-1, xc);
                  
                  if (var.DiagonalMass())
                    var.MltInvSqrtDiagonalMass(xc);
                  
                  var.MltStiffness(SeldonConjTrans, xc, yc);
                  
                  if (var.DiagonalMass())
                    var.MltInvSqrtDiagonalMass(yc);
                  
                  SetCol(yc, k-1, zwork);
                }
          }
        else if (ijob == 40)
          {
            // multiplication by matrix B
            int i = fpm(23), j = fpm(23) + fpm(24)-1;
            if (herm)
              for (int k = i; k <= j; k++)
                {
                  for (int p = 0; p < n; p++)
                    x(p) = eigen_vectors(p, k-1);
                  
                  if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
                    Copy(x, y);
                  else
                    var.MltMass(x, y);
                  
                  for (int p = 0; p < n; p++)
                    work(p, k-1) = y(p);
                }
            else
              for (int k = i; k <= j; k++)
                {
                  GetCol(eigen_vec_cplx, k-1, xc);
                  
                  if (var.DiagonalMass())
                    Copy(xc, yc);
                  else
                    var.MltMass(xc, yc);

                  SetCol(yc, k-1, zwork);
                }
          }
        else if (ijob == 41)
          {
            // multiplication by matrix B^H
            int i = fpm(33), j = fpm(33) + fpm(34)-1;
            if (herm)
              for (int k = i; k <= j; k++)
                {
                  for (int p = 0; p < n; p++)
                    x(p) = eigen_vectors(p, k-1);
                  
                  if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
                    Copy(x, y);
                  else
                    var.MltMass(SeldonConjTrans, x, y);
                  
                  for (int p = 0; p < n; p++)
                    work(p, k-1) = y(p);
                }
            else
              for (int k = i; k <= j; k++)
                {
                  GetCol(eigen_vec_cplx, k-1, xc);
                  
                  if (var.DiagonalMass())
                    Copy(xc, yc);
                  else
                    var.MltMass(SeldonConjTrans, xc, yc);
                  
                  SetCol(yc, k-1, zwork);
                }
          }
      }
    
    if (info != 0)
      {
        if (info == 3)
          {
            cout << "Feast returns error code 3 " << endl;
            cout << "It usually means that there are more than "
                 << m0-1 << " eigenvalues in the ";

            if (sym)
              cout << "interval [" << emin << ", " << emax << "] " << endl;
            else
              cout << "disk of center " << Emid << " and radius" << r << endl;
            
            cout << "You can try to change the ";
            if (sym)
              cout << "interval ";
            else
              cout << "disk ";
            
            cout << "or increase the number of asked eigenvalues" << endl;
            
            abort();
          }
        else if (info == 1)
          {
            cout << "No eigenvalues have been found in the ";
            if (sym)
              cout << "interval [" << emin << ", " << emax << "] " << endl;
            else
              cout << "disk of center " << Emid << " and radius" << r << endl;
            
            cout << "You can try to change the ";
            if (sym)
              cout << "interval. " << endl;
            else
              cout << "disk. " << endl;

            abort();
          }
        else if (info == 2)
          {
            cout << "Feast did not converge, loop = " << loop << endl;
            cout << "Maximum number of iterations = " << var.GetNbMaximumIterations() << endl;
            cout << "Residual = " << res << endl;
          }
        else if (info == 5)
          {
            if (var.GetGlobalRankCommunicator() == 0)
              cout << "Estimated number of eigenvalues = " << m << endl;
            
            return;
          }
        else if (info == 6)
          {
            if (var.GetGlobalRankCommunicator() == 0)
              cout << "Warning : subspace is not biorthogonal" << endl;
            
            /*if (sym)
              for (int i = 0; i < m; i++)
                cout << "Errors on lambda " << i << " = " << res(i) << endl;
            else
              for (int i = 0; i < m; i++)
                cout << "Errors on lambda " << i << " = " << res(i) << " " << res(i+m0) << endl;
            */
          }
        else if (info > 6)
          {
            cout << "An error occurred during feast eigenvalue solving " << endl;
            cout << "info = " << info << endl;
            abort();
          }
      }
    
    if ((var.GetGlobalPrintLevel() > 0) && (var.GetGlobalRankCommunicator() == 0))
      cout << "Feast converged in " << loop << " cycles" << endl;
    
    work.Clear(); workc.Clear(); zwork.Clear();
    aq.Clear(); bq.Clear(); zaq.Clear(); zbq.Clear();
    
    eigen_values.Reallocate(m);   
    lambda_imag.Reallocate(m);
    lambda_imag.Zero(); 
    if (herm)
      {
        for (int i = 0; i < m; i++)
          eigen_values(i) = lambda(i);
        
        eigen_vectors.Resize(n, m);
      }
    else
      var.FillComplexEigenvectors(m, Emid, tol, lambda_cplx, eigen_vec_cplx,
                                  eigen_values, lambda_imag, eigen_vectors);
    
    lambda.Clear();
    lambda_cplx.Clear(); eigen_vec_cplx.Clear();
    
    T shiftr = var.GetShiftValue(), shifti = var.GetImagShiftValue();
    ApplyScalingEigenvec(var, eigen_values, lambda_imag, eigen_vectors,
                         shiftr, shifti);

  }
  

  void CallFeast(int& ijob, int n, bool herm, bool sym, complex<double>& ze,
		 Matrix<double, General, ColMajor>& work,
		 Matrix<complex<double>, General, ColMajor>& workc,
		 Vector<double>& aq, Vector<double>& bq,
		 Matrix<complex<double>, General, ColMajor>& zwork,
                 Vector<complex<double> >& zaq, Vector<complex<double> >& zbq,
		 IVect& fpm, double& epsout, int& loop,
		 double emin, double emax, complex<double> Emid, double r,
                 int m0, Vector<double>& lambda,
		 Matrix<double, General, ColMajor>& eigen_vectors,
                 Vector<complex<double> >& lambda_cplx,
                 Matrix<complex<double>, General, ColMajor>& eigen_vec_cplx,
		 int& m, Vector<double>& res, int& info)
  {
    if (sym)
      dfeast_srci_(&ijob, &n, reinterpret_cast<void*>(&ze), work.GetData(),
                   workc.GetDataVoid(), aq.GetData(), bq.GetData(),
                   fpm.GetData(), &epsout, &loop, &emin, &emax,
                   &m0, lambda.GetData(), eigen_vectors.GetData(), 
                   &m, res.GetData(), &info);
    else
      dfeast_grci_(&ijob, &n, reinterpret_cast<void*>(&ze), zwork.GetDataVoid(),
                   workc.GetDataVoid(), zaq.GetDataVoid(), zbq.GetDataVoid(),
                   fpm.GetData(), &epsout, &loop, &Emid, &r,
                   &m0, lambda_cplx.GetDataVoid(), eigen_vec_cplx.GetDataVoid(),
                   &m, res.GetData(), &info);
  }
  
  void CallFeast(int& ijob, int n, bool herm, bool sym, complex<double>& ze,
		 Matrix<complex<double>, General, ColMajor>& work,
		 Matrix<complex<double>, General, ColMajor>& workc,
		 Vector<complex<double> >& aq, Vector<complex<double> >& bq,
                 Matrix<complex<double>, General, ColMajor>& zwork,
                 Vector<complex<double> >& zaq, Vector<complex<double> >& zbq,
		 IVect& fpm, double& epsout, int& loop,
		 double emin, double emax, complex<double> Emid, double r,
                 int m0, Vector<double>& lambda,
		 Matrix<complex<double>, General, ColMajor>& eigen_vectors,
                 Vector<complex<double> >& lambda_cplx,
                 Matrix<complex<double>, General, ColMajor>& eigen_vec_cplx,
		 int& m, Vector<double>& res, int& info)
  {
    if (herm)
      zfeast_hrci_(&ijob, &n, reinterpret_cast<void*>(&ze), work.GetDataVoid(),
                   workc.GetDataVoid(), aq.GetDataVoid(), bq.GetDataVoid(),
                   fpm.GetData(), &epsout, &loop, &emin, &emax, &m0,
                   lambda.GetData(), eigen_vectors.GetDataVoid(), 
                   &m, res.GetData(), &info);
    else if (sym)
      zfeast_srci_(&ijob, &n, reinterpret_cast<void*>(&ze), zwork.GetDataVoid(),
                   workc.GetDataVoid(), zaq.GetDataVoid(), zbq.GetDataVoid(),
                   fpm.GetData(), &epsout, &loop, &Emid, &r, &m0,
                   lambda_cplx.GetData(), eigen_vec_cplx.GetDataVoid(), 
                   &m, res.GetData(), &info);
    else
      zfeast_grci_(&ijob, &n, reinterpret_cast<void*>(&ze), zwork.GetDataVoid(),
                   workc.GetDataVoid(), zaq.GetDataVoid(), zbq.GetDataVoid(),
                   fpm.GetData(), &epsout, &loop, &Emid, &r, &m0,
                   lambda_cplx.GetData(), eigen_vec_cplx.GetDataVoid(), 
                   &m, res.GetData(), &info);
  }
  
}

#endif
