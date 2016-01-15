#ifndef BOUND_MAT_H
#define BOUND_MAT_H

#include "flint/flint.h"
#include "flint/fmpz_mat.h"

#include "tkf91_generator_indices.h"
#include "dp.h"



#ifdef __cplusplus
extern "C" {
#endif


void
tkf91_dp_verify_symbolically(
        int *verified,
        fmpz_mat_t mat,
        const tkf91_generator_indices_t g,
        dp_mat_t tableau,
        expr_ptr * expressions_table,
        const slong *A,
        const slong *B);


#ifdef __cplusplus
}
#endif

#endif
