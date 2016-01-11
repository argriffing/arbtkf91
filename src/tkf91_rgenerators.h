#ifndef TKF91_RGENERATORS_H
#define TKF91_RGENERATORS_H

#include "expressions.h"
#include "generators.h"
#include "rgenerators.h"
#include "tkf91_generator_indices.h"


#ifdef __cplusplus
extern "C" {
#endif

void tkf91_rgenerators_init(
        tkf91_generator_indices_t x,
        rgen_reg_ptr g,
        tkf91_rationals_t r,
        tkf91_expressions_t p,
        const slong *A, slong Alen,
        const slong *B, slong Blen);
void tkf91_rgenerators_clear(tkf91_generator_indices_t x);

#ifdef __cplusplus
}
#endif

#endif
