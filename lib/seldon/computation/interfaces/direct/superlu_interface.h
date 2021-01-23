#ifndef SELDON_FILE_SUPERLU_INTERFACE_H /* allow multiple inclusions */
#define SELDON_FILE_SUPERLU_INTERFACE_H

#ifdef SELDON_WITH_SUPERLU_DIST

namespace superlu
{

#ifdef SELDON_WITH_SUPERLU_DOUBLE

#include "superlu_ddefs.h"
#undef Reduce

#else

#include "superlu_zdefs.h"
#undef Reduce

#endif
}

#else

namespace superlu
{
#ifdef SELDON_WITH_SUPERLU_MT
#include "slu_mt_zdefs.h"
#else
#include "slu_zdefs.h"
#endif
}

#endif

using superlu::int_t;
using superlu::colperm_t;

// since slu_zdefs.h and slu_ddefs.h can not be included
// at the same time, needed functions of slu_ddefs are copied here
// version : SuperLU 5.0, SuperLU_MT 3.0, SuperLU_DIST 4.1
extern "C"
{

namespace superlu
{  

#ifdef SELDON_WITH_SUPERLU_DIST
  int_t dQuerySpace_dist(int_t, LUstruct_t *, gridinfo_t *,
                         SuperLUStat_t *, mem_usage_t *);
  
  void dCreate_CompRowLoc_Matrix_dist(SuperMatrix *, int_t, int_t, int_t, int_t,
                                      int_t, double *, int_t *, int_t *,
                                      Stype_t, Dtype_t, Mtype_t);
  
  void  pdgssvx(superlu_options_t *, SuperMatrix *, 
                ScalePermstruct_t *, double *,
                int, int, gridinfo_t *, LUstruct_t *,
                SOLVEstruct_t *, double *, SuperLUStat_t *, int *);

 
#else
  
  void
  dCreate_CompCol_Matrix(SuperMatrix *, int_t, int_t, int_t, double *,
                         int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
  
  void
  dCreate_Dense_Matrix(SuperMatrix *, int_t, int_t, double *, int_t,
                       Stype_t, Dtype_t, Mtype_t);
    
#ifdef SELDON_WITH_SUPERLU_MT
  int_t  superlu_dQuerySpace (int_t, SuperMatrix *, SuperMatrix *, int_t, 
                              superlu_memusage_t *);
  
  void pdgstrf (superlumt_options_t *, SuperMatrix *, int_t *, 
                SuperMatrix *, SuperMatrix *, Gstat_t *, int_t *);

  void pdgstrf_init (int_t, fact_t, trans_t, yes_no_t, int_t, int_t, double, yes_no_t, double,
                     int_t *, int_t *, void *, int_t, SuperMatrix *,
                     SuperMatrix *, superlumt_options_t *, Gstat_t *);
  
  void dgstrs (trans_t, SuperMatrix *, SuperMatrix*, 
               int_t*, int_t*, SuperMatrix*, Gstat_t *, int_t *);  
#else
  void    dgstrf (superlu_options_t*, SuperMatrix*,
                  int, int, int*, void *, int, int *, int *, 
                  SuperMatrix *, SuperMatrix *, GlobalLU_t *,
                  SuperLUStat_t*, int *);
  
  void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                  SuperMatrix *, SuperLUStat_t*, int *);
  
  int     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
#endif

#endif
}
  
}

#endif /* __SUPERLU_INTERFACE */
