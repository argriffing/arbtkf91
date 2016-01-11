#ifndef TKF91_GENERATORS_H
#define TKF91_GENERATORS_H

#include "expressions.h"
#include "generators.h"
#include "tkf91_generator_indices.h"




#ifdef __cplusplus
extern "C" {
#endif

void tkf91_generators_init(
        tkf91_generator_indices_t x,
        generator_reg_t g,
        tkf91_expressions_t p,
        slong *A, slong Alen,
        slong *B, slong Blen);
void tkf91_generators_clear(tkf91_generator_indices_t x);

#ifdef __cplusplus
}
#endif

#endif
