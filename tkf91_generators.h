#ifndef TKF91_GENERATORS_H
#define TKF91_GENERATORS_H

#include "expressions.h"
#include "generators.h"


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
} tkf91_generators_struct;
typedef tkf91_generators_struct tkf91_generators_t[1];
typedef tkf91_generators_struct * tkf91_generators_ptr;



#ifdef __cplusplus
extern "C" {
#endif

void tkf91_generators_init(
        tkf91_generators_t x,
        generator_reg_t g,
        tkf91_expressions_t p,
        slong *A, slong Alen,
        slong *B, slong Blen);
void tkf91_generators_clear(tkf91_generators_t x);

#ifdef __cplusplus
}
#endif

#endif
