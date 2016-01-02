#ifndef TKF91_DP_F_H
#define TKF91_DP_F_H

#include "flint/flint.h"
#include "flint/fmpz_mat.h"

#include "femtocas.h"
#include "tkf91_generator_indices.h"
#include "tkf91_dp.h"


#ifdef __cplusplus
extern "C" {
#endif


void tkf91_dp_f(
        solution_t, const request_t,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g,
        slong *A, size_t szA,
        slong *B, size_t szB);


#ifdef __cplusplus
}
#endif

#endif