
#ifndef SELDON_FILE_ARPACK_INLINE_CXX
#define SELDON_FILE_ARPACK_INLINE_CXX

namespace Seldon
{
  
#ifdef SELDON_WITH_MPI
  //! Symmetric problems in double precision.
  inline void saupd(ARPACK_INTEGER comm, ARPACK_INTEGER& ido, char bmat,
		    ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_DOUBLEREAL& tol, ARPACK_DOUBLEREAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
		    ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    pdsaupd_(&comm, &ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
	     V, &ldv, iparam, ipntr, workd, workl,
	     &lworkl, &info);
  }
  
  //! Symmetric problems in single precision.
  inline void saupd(ARPACK_INTEGER comm, ARPACK_INTEGER& ido, char bmat,
		    ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_REAL& tol, ARPACK_REAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
		    ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    pssaupd_(&comm, &ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
	     V, &ldv, iparam, ipntr, workd, workl,
	     &lworkl, &info);
  }
  
  //! Postprocessing for symmetric problems in double precision.
  inline void seupd(ARPACK_INTEGER comm,
		    ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
		    ARPACK_DOUBLEREAL* d, ARPACK_DOUBLEREAL* Z,
		    ARPACK_INTEGER ldz, ARPACK_DOUBLEREAL sigma, char bmat,
		    ARPACK_INTEGER n, char* which, ARPACK_INTEGER nev,
		    ARPACK_DOUBLEREAL tol, ARPACK_DOUBLEREAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
		    ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    pdseupd_(&comm, &rvec, &HowMny, select, d, Z, &ldz, &sigma, &bmat,
	     &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam,
	     ipntr, workd, workl, &lworkl, &info );
  }
  
  //! Postprocessing for symmetric problems in single precision.
  inline void seupd(ARPACK_INTEGER comm, ARPACK_LOGICAL rvec,
		    char HowMny, ARPACK_LOGICAL *select,
		    ARPACK_REAL* d, ARPACK_REAL* Z, ARPACK_INTEGER ldz,
		    ARPACK_REAL sigma, char bmat, ARPACK_INTEGER n,
		    char* which, ARPACK_INTEGER nev, ARPACK_REAL tol,
		    ARPACK_REAL* resid, ARPACK_INTEGER ncv, ARPACK_REAL* V,
		    ARPACK_INTEGER ldv, ARPACK_INTEGER* iparam,
		    ARPACK_INTEGER* ipntr, ARPACK_REAL* workd, ARPACK_REAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    psseupd_(&comm, &rvec, &HowMny, select, d, Z, &ldz, &sigma, &bmat,
	     &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam,
	     ipntr, workd, workl, &lworkl, &info);
  }
  
  //! Non-symmetric problems in double precision.
  inline void naupd(ARPACK_INTEGER comm,
		    ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_DOUBLEREAL& tol,
		    ARPACK_DOUBLEREAL* resid, ARPACK_INTEGER ncv,
		    ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv, ARPACK_INTEGER* iparam,
		    ARPACK_INTEGER* ipntr, ARPACK_DOUBLEREAL* workd,
		    ARPACK_DOUBLEREAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    pdnaupd_(&comm, &ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
	     V, &ldv, iparam, ipntr, workd, workl,
	     &lworkl, &info);
  }
  
  //! Non-symmetric problems in single precision.
  inline void naupd(ARPACK_INTEGER comm,
		    ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_REAL& tol, ARPACK_REAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
		    ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    
    psnaupd_(&comm, &ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
	     V, &ldv, iparam, ipntr, workd, workl,
	     &lworkl, &info);
    
  }
  
  //! Postprocessing for non-symmetric problems in double precision.
  inline void neupd(ARPACK_INTEGER comm,
		    ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
		    ARPACK_DOUBLEREAL* dr, ARPACK_DOUBLEREAL* di,
		    ARPACK_DOUBLEREAL* Z, ARPACK_INTEGER ldz, ARPACK_DOUBLEREAL sigmar,
		    ARPACK_DOUBLEREAL sigmai, ARPACK_DOUBLEREAL* workv, char bmat,
		    ARPACK_INTEGER n, char* which, ARPACK_INTEGER nev,
		    ARPACK_DOUBLEREAL tol, ARPACK_DOUBLEREAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
		    ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    pdneupd_(&comm, &rvec, &HowMny, select, dr, di, Z, &ldz, &sigmar,
	     &sigmai, workv, &bmat, &n, which, &nev, &tol,
	     resid, &ncv, V, &ldv, iparam, ipntr,
	     workd, workl, &lworkl, &info);
  }
  
  //! Postprocessing for non-symmetric problems in single precision.
  inline void neupd(ARPACK_INTEGER comm,
		    ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
		    ARPACK_REAL* dr, ARPACK_REAL* di, ARPACK_REAL* Z,
		    ARPACK_INTEGER ldz, ARPACK_REAL sigmar, ARPACK_REAL sigmai,
		    ARPACK_REAL* workv, char bmat, ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_REAL tol, ARPACK_REAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
		    ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    psneupd_(&comm, &rvec, &HowMny, select, dr, di, Z, &ldz, &sigmar,
	     &sigmai, workv, &bmat, &n, which, &nev, &tol,
	     resid, &ncv, V, &ldv, iparam, ipntr,
	     workd, workl, &lworkl, &info);
  }
#else
  //! Symmetric problems in double precision.
  inline void saupd(ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_DOUBLEREAL& tol, ARPACK_DOUBLEREAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
		    ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
	    V, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);
  }
  
  //! Symmetric problems in single precision.
  inline void saupd(ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_REAL& tol, ARPACK_REAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
		    ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    ssaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
	    V, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);
  }
  
  //! Postprocessing for symmetric problems in double precision.
  inline void seupd(ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
		    ARPACK_DOUBLEREAL* d, ARPACK_DOUBLEREAL* Z,
		    ARPACK_INTEGER ldz, ARPACK_DOUBLEREAL sigma, char bmat,
		    ARPACK_INTEGER n, char* which, ARPACK_INTEGER nev,
		    ARPACK_DOUBLEREAL tol, ARPACK_DOUBLEREAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
		    ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    dseupd_(&rvec, &HowMny, select, d, Z, &ldz, &sigma, &bmat,
	    &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam,
	    ipntr, workd, workl, &lworkl, &info );
  }
  
  //! Postprocessing for symmetric problems in single precision.
  inline void seupd(ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
		    ARPACK_REAL* d, ARPACK_REAL* Z, ARPACK_INTEGER ldz,
		    ARPACK_REAL sigma, char bmat, ARPACK_INTEGER n,
		    char* which, ARPACK_INTEGER nev, ARPACK_REAL tol,
		    ARPACK_REAL* resid, ARPACK_INTEGER ncv, ARPACK_REAL* V,
		    ARPACK_INTEGER ldv, ARPACK_INTEGER* iparam,
		    ARPACK_INTEGER* ipntr, ARPACK_REAL* workd, ARPACK_REAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    sseupd_(&rvec, &HowMny, select, d, Z, &ldz, &sigma, &bmat,
	    &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam,
	    ipntr, workd, workl, &lworkl, &info);
  }
  
  //! Non-symmetric problems in double precision.
  inline void naupd(ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_DOUBLEREAL& tol,
		    ARPACK_DOUBLEREAL* resid, ARPACK_INTEGER ncv,
		    ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv, ARPACK_INTEGER* iparam,
		    ARPACK_INTEGER* ipntr, ARPACK_DOUBLEREAL* workd,
		    ARPACK_DOUBLEREAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    dnaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
	    V, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);
    
  }
  
  //! Non-symmetric problems in single precision.
  inline void naupd(ARPACK_INTEGER& ido, char bmat, ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_REAL& tol, ARPACK_REAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
		    ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    
    snaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
	    V, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);
    
  }
  
  //! Postprocessing for non-symmetric problems in double precision.
  inline void neupd(ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
		    ARPACK_DOUBLEREAL* dr, ARPACK_DOUBLEREAL* di,
		    ARPACK_DOUBLEREAL* Z, ARPACK_INTEGER ldz, ARPACK_DOUBLEREAL sigmar,
		    ARPACK_DOUBLEREAL sigmai, ARPACK_DOUBLEREAL* workv, char bmat,
		    ARPACK_INTEGER n, char* which, ARPACK_INTEGER nev,
		    ARPACK_DOUBLEREAL tol, ARPACK_DOUBLEREAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_DOUBLEREAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr,
		    ARPACK_DOUBLEREAL* workd, ARPACK_DOUBLEREAL* workl,
		    ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    dneupd_(&rvec, &HowMny, select, dr, di, Z, &ldz, &sigmar,
	    &sigmai, workv, &bmat, &n, which, &nev, &tol,
	    resid, &ncv, V, &ldv, iparam, ipntr,
	    workd, workl, &lworkl, &info);
  }
  
  //! Postprocessing for non-symmetric problems in single precision.
  inline void neupd(ARPACK_LOGICAL rvec, char HowMny, ARPACK_LOGICAL *select,
		    ARPACK_REAL* dr, ARPACK_REAL* di, ARPACK_REAL* Z,
		    ARPACK_INTEGER ldz, ARPACK_REAL sigmar, ARPACK_REAL sigmai,
		    ARPACK_REAL* workv, char bmat, ARPACK_INTEGER n, char* which,
		    ARPACK_INTEGER nev, ARPACK_REAL tol, ARPACK_REAL* resid,
		    ARPACK_INTEGER ncv, ARPACK_REAL* V, ARPACK_INTEGER ldv,
		    ARPACK_INTEGER* iparam, ARPACK_INTEGER* ipntr, ARPACK_REAL* workd,
		    ARPACK_REAL* workl, ARPACK_INTEGER lworkl, ARPACK_INTEGER& info)
  {
    sneupd_(&rvec, &HowMny, select, dr, di, Z, &ldz, &sigmar,
	    &sigmai, workv, &bmat, &n, which, &nev, &tol,
	    resid, &ncv, V, &ldv, iparam, ipntr,
	    workd, workl, &lworkl, &info);
  }
#endif

} // end namespace Seldon


#endif
