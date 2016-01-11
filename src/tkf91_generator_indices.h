#ifndef TKF91_GENERATOR_INDICES_H
#define TKF91_GENERATOR_INDICES_H

#include "flint/flint.h"


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
} tkf91_generator_indices_struct;
typedef tkf91_generator_indices_struct tkf91_generator_indices_t[1];
typedef tkf91_generator_indices_struct * tkf91_generator_indices_ptr;

#endif
