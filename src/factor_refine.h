#ifndef FACTOR_REFINEMENT_H
#define FACTOR_REFINEMENT_H

#include "flint/flint.h"
#include "flint/fmpz.h"

#ifdef __cplusplus
extern "C" {
#endif

void fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f);

#ifdef __cplusplus
}
#endif

#endif
