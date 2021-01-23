// Headers from ARPACK.

// Modifications (by Lin Wu):
// Replacements:
//    integer       --> ARPACK_INTEGER
//    real          --> ARPACK_REAL
//    doublereal    --> ARPACK_DOUBLEREAL
//    complex       --> ARPACK_COMPLEX
//    doublecomplex --> ARPACK_DOUBLECOMPLEX
//    logical       --> ARPACK_LOGICAL

#ifndef __CARPACK_H
#define __CARPACK_H
 

struct ArpackLog
{
  ARPACK_INTEGER logfil, ndigit, mgetv0;
  ARPACK_INTEGER msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
  ARPACK_INTEGER mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
  ARPACK_INTEGER mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
};

extern ArpackLog debug_;

#ifdef SELDON_WITH_MPI
// double precision symmetric routines.
void pdsaupd_(ARPACK_INTEGER *comm, ARPACK_INTEGER *ido, char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_DOUBLEREAL *tol, ARPACK_DOUBLEREAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_DOUBLEREAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr, ARPACK_DOUBLEREAL *workd,
	     ARPACK_DOUBLEREAL *workl, ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);

void pdseupd_(ARPACK_INTEGER *comm, ARPACK_LOGICAL *rvec, char *HowMny, ARPACK_LOGICAL *select,
	     ARPACK_DOUBLEREAL *d, ARPACK_DOUBLEREAL *Z, ARPACK_INTEGER *ldz,
	     ARPACK_DOUBLEREAL *sigma, char *bmat, ARPACK_INTEGER *n,
	     char *which, ARPACK_INTEGER *nev, ARPACK_DOUBLEREAL *tol,
	     ARPACK_DOUBLEREAL *resid, ARPACK_INTEGER *ncv, ARPACK_DOUBLEREAL *V,
	     ARPACK_INTEGER *ldv, ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr,
	     ARPACK_DOUBLEREAL *workd, ARPACK_DOUBLEREAL *workl,
	     ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);


// double precision nonsymmetric routines.
void pdnaupd_(ARPACK_INTEGER *comm, ARPACK_INTEGER *ido, char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_DOUBLEREAL *tol, ARPACK_DOUBLEREAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_DOUBLEREAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr, ARPACK_DOUBLEREAL *workd,
	     ARPACK_DOUBLEREAL *workl, ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);
void pdneupd_(ARPACK_INTEGER *comm, ARPACK_LOGICAL *rvec, char *HowMny, ARPACK_LOGICAL *select,
	     ARPACK_DOUBLEREAL *dr, ARPACK_DOUBLEREAL *di, ARPACK_DOUBLEREAL *Z,
	     ARPACK_INTEGER *ldz, ARPACK_DOUBLEREAL *sigmar,
	     ARPACK_DOUBLEREAL *sigmai, ARPACK_DOUBLEREAL *workev,
	     char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_DOUBLEREAL *tol, ARPACK_DOUBLEREAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_DOUBLEREAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr,
	     ARPACK_DOUBLEREAL *workd, ARPACK_DOUBLEREAL *workl,
	     ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);


// single precision symmetric routines.
void pssaupd_(ARPACK_INTEGER *comm, ARPACK_INTEGER *ido, char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_REAL *tol, ARPACK_REAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_REAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr, ARPACK_REAL *workd,
	     ARPACK_REAL *workl, ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);
void psseupd_(ARPACK_INTEGER *comm, ARPACK_LOGICAL *rvec, char *HowMny, ARPACK_LOGICAL *select,
	     ARPACK_REAL *d, ARPACK_REAL *Z, ARPACK_INTEGER *ldz,
	     ARPACK_REAL *sigma, char *bmat, ARPACK_INTEGER *n,
	     char *which, ARPACK_INTEGER *nev, ARPACK_REAL *tol,
	     ARPACK_REAL *resid, ARPACK_INTEGER *ncv, ARPACK_REAL *V,
	     ARPACK_INTEGER *ldv, ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr,
	     ARPACK_REAL *workd, ARPACK_REAL *workl,
	     ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);


// single precision nonsymmetric routines.
void psnaupd_(ARPACK_INTEGER *comm, ARPACK_INTEGER *ido, char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_REAL *tol, ARPACK_REAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_REAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr, ARPACK_REAL *workd,
	     ARPACK_REAL *workl, ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);
void psneupd_(ARPACK_INTEGER *comm, ARPACK_LOGICAL *rvec, char *HowMny, ARPACK_LOGICAL *select,
	     ARPACK_REAL *dr, ARPACK_REAL *di, ARPACK_REAL *Z,
	     ARPACK_INTEGER *ldz, ARPACK_REAL *sigmar,
	     ARPACK_REAL *sigmai, ARPACK_REAL *workev, char *bmat,
	     ARPACK_INTEGER *n, char *which, ARPACK_INTEGER *nev,
	     ARPACK_REAL *tol, ARPACK_REAL *resid, ARPACK_INTEGER *ncv,
	     ARPACK_REAL *V, ARPACK_INTEGER *ldv, ARPACK_INTEGER *iparam,
	     ARPACK_INTEGER *ipntr, ARPACK_REAL *workd, ARPACK_REAL *workl,
	     ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);

// complex routines
void pznaupd_(int*, int*, char*, int*, char*, int*, double*, void*,
	     int*, void*, int*, int*, int*, void*, void*, int*, double*, int*);

void pzneupd_(int*, int*, char*, int*, void*, void*, int*, void*, void*, char*, int*,
	     char*, int*, double*, void*, int*, void*, int*, int*, int*, void*,
	     void*, int*, double*, int*);
#endif

// double precision symmetric routines.
void dsaupd_(ARPACK_INTEGER *ido, char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_DOUBLEREAL *tol, ARPACK_DOUBLEREAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_DOUBLEREAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr, ARPACK_DOUBLEREAL *workd,
	     ARPACK_DOUBLEREAL *workl, ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);

void dseupd_(ARPACK_LOGICAL *rvec, char *HowMny, ARPACK_LOGICAL *select,
	     ARPACK_DOUBLEREAL *d, ARPACK_DOUBLEREAL *Z, ARPACK_INTEGER *ldz,
	     ARPACK_DOUBLEREAL *sigma, char *bmat, ARPACK_INTEGER *n,
	     char *which, ARPACK_INTEGER *nev, ARPACK_DOUBLEREAL *tol,
	     ARPACK_DOUBLEREAL *resid, ARPACK_INTEGER *ncv, ARPACK_DOUBLEREAL *V,
	     ARPACK_INTEGER *ldv, ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr,
	     ARPACK_DOUBLEREAL *workd, ARPACK_DOUBLEREAL *workl,
	     ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);


// double precision nonsymmetric routines.
void dnaupd_(ARPACK_INTEGER *ido, char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_DOUBLEREAL *tol, ARPACK_DOUBLEREAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_DOUBLEREAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr, ARPACK_DOUBLEREAL *workd,
	     ARPACK_DOUBLEREAL *workl, ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);
void dneupd_(ARPACK_LOGICAL *rvec, char *HowMny, ARPACK_LOGICAL *select,
	     ARPACK_DOUBLEREAL *dr, ARPACK_DOUBLEREAL *di, ARPACK_DOUBLEREAL *Z,
	     ARPACK_INTEGER *ldz, ARPACK_DOUBLEREAL *sigmar,
	     ARPACK_DOUBLEREAL *sigmai, ARPACK_DOUBLEREAL *workev,
	     char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_DOUBLEREAL *tol, ARPACK_DOUBLEREAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_DOUBLEREAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr,
	     ARPACK_DOUBLEREAL *workd, ARPACK_DOUBLEREAL *workl,
	     ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);


// single precision symmetric routines.
void ssaupd_(ARPACK_INTEGER *ido, char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_REAL *tol, ARPACK_REAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_REAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr, ARPACK_REAL *workd,
	     ARPACK_REAL *workl, ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);
void sseupd_(ARPACK_LOGICAL *rvec, char *HowMny, ARPACK_LOGICAL *select,
	     ARPACK_REAL *d, ARPACK_REAL *Z, ARPACK_INTEGER *ldz,
	     ARPACK_REAL *sigma, char *bmat, ARPACK_INTEGER *n,
	     char *which, ARPACK_INTEGER *nev, ARPACK_REAL *tol,
	     ARPACK_REAL *resid, ARPACK_INTEGER *ncv, ARPACK_REAL *V,
	     ARPACK_INTEGER *ldv, ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr,
	     ARPACK_REAL *workd, ARPACK_REAL *workl,
	     ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);


// single precision nonsymmetric routines.
void snaupd_(ARPACK_INTEGER *ido, char *bmat, ARPACK_INTEGER *n, char *which,
	     ARPACK_INTEGER *nev, ARPACK_REAL *tol, ARPACK_REAL *resid,
	     ARPACK_INTEGER *ncv, ARPACK_REAL *V, ARPACK_INTEGER *ldv,
	     ARPACK_INTEGER *iparam, ARPACK_INTEGER *ipntr, ARPACK_REAL *workd,
	     ARPACK_REAL *workl, ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);
void sneupd_(ARPACK_LOGICAL *rvec, char *HowMny, ARPACK_LOGICAL *select,
	     ARPACK_REAL *dr, ARPACK_REAL *di, ARPACK_REAL *Z,
	     ARPACK_INTEGER *ldz, ARPACK_REAL *sigmar,
	     ARPACK_REAL *sigmai, ARPACK_REAL *workev, char *bmat,
	     ARPACK_INTEGER *n, char *which, ARPACK_INTEGER *nev,
	     ARPACK_REAL *tol, ARPACK_REAL *resid, ARPACK_INTEGER *ncv,
	     ARPACK_REAL *V, ARPACK_INTEGER *ldv, ARPACK_INTEGER *iparam,
	     ARPACK_INTEGER *ipntr, ARPACK_REAL *workd, ARPACK_REAL *workl,
	     ARPACK_INTEGER *lworkl, ARPACK_INTEGER *info);

// complex routines
void znaupd_(int*, char*, int*, char*, int*, double*, void*,
	     int*, void*, int*, int*, int*, void*, void*, int*, double*, int*);

void zneupd_(int*, char*, int*, void*, void*, int*, void*, void*, char*, int*,
	     char*, int*, double*, void*, int*, void*, int*, int*, int*, void*,
	     void*, int*, double*, int*);


#endif /* __CARPACK_H */
