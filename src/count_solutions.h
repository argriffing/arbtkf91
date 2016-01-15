#ifndef COUNT_SOLUTIONS_H
#define COUNT_SOLUTIONS_H

#include "flint/flint.h"
#include "flint/fmpz.h"

#include "dp.h"


#ifdef __cplusplus
extern "C" {
#endif


void count_solutions(fmpz_t res, dp_mat_t mat);


#ifdef __cplusplus
}
#endif

#endif

