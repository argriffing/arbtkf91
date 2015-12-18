#ifndef TKF91_RGENERATORS_H
#define TKF91_RGENERATORS_H

#include "expressions.h"
#include "rgenerators.h"


/* These are indices of registered generators. */
typedef struct
{
    slong m1_00;
    slong m0_10;
    slong m0_i0_incr[4];
    slong m2_01;
    slong m2_0j_incr[4];
    slong c0_incr[4];
    slong c1_incr[16];
    slong c2_incr[4];
} tkf91_rgenerators_struct;
typedef tkf91_rgenerators_struct tkf91_rgenerators_t[1];
typedef tkf91_rgenerators_struct * tkf91_rgenerators_ptr;



#ifdef __cplusplus
extern "C" {
#endif

void tkf91_rgenerators_init(
        tkf91_rgenerators_t x,
        rgen_reg_ptr g,
        tkf91_rationals_t r,
        tkf91_expressions_t p,
        slong *A, slong Alen,
        slong *B, slong Blen);
void tkf91_rgenerators_clear(tkf91_rgenerators_t x);

#ifdef __cplusplus
}
#endif

#endif
