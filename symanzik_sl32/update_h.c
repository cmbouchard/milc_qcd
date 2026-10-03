/****** update_h.c  -- update the momentum matrices ******************/
/* MIMD version 7 */

/* Modifications:
   2/17/98  ANSI prototyping,
            improved gauge updating taken from "ks_imp_dyn" U.M.H.
   */

#include "symanzik_sl32_includes.h"

#if defined(USE_GF_GPU) && !defined(LAYOUT_HYPER_PRIME)
#error "QUDA gauge force needs LAYOUT=layout_hyper_prime.o ADDDEFINES=-DLAYOUT_HYPER_PRIME"
#endif

void update_h(Real eps) {
    /* gauge field force */
#ifdef USE_GF_GPU
    /* QUDA force; always acts on the site links and site mom */
    imp_gauge_force_gpu(eps,F_OFFSET(mom));
#else
    imp_gauge_force(eps,F_OFFSET(mom));
#endif
} /* update_h */
