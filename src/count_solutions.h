#ifndef COUNT_SOLUTIONS_H
#define COUNT_SOLUTIONS_H

#include "flint/flint.h"
#include "flint/fmpz.h"

#include "breadcrumbs.h"


#ifdef __cplusplus
extern "C" {
#endif


void count_solutions(fmpz_t res, const breadcrumb_mat_t mat);


#ifdef __cplusplus
}
#endif

#endif

