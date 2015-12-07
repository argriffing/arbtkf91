/* factor_refinement -- integer factor refinement
 */

#ifndef FACTOR_REFINEMENT_H
#define FACTOR_REFINEMENT_H

#include "flint/flint.h"
#include "flint/fmpz.h"

#ifdef __cplusplus
extern "C" {
#endif

void factor_refinement(fmpz **ybase, fmpz **yexp, slong *ylen,
        const fmpz *x, const slong xlen);

#ifdef __cplusplus
}
#endif

#endif
