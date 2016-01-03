#ifndef BOUND_MAT_H
#define BOUND_MAT_H

#include "flint/flint.h"
#include "flint/fmpz_mat.h"

#include "tkf91_generator_indices.h"
#include "breadcrumbs.h"



#ifdef __cplusplus
extern "C" {
#endif


void
tkf91_dp_verify_symbolically(
        fmpz_mat_t mat,
        const tkf91_generator_indices_t g,
        breadcrumb_mat_t mask,
        const slong *A,
        const slong *B);


#ifdef __cplusplus
}
#endif

#endif
