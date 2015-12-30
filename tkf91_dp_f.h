#ifndef TKF91_F_H
#define TKF91_F_H

#include "flint.h"
#include "arb_mat.h"

#include "femtocas.h"
#include "tkf91_generator_indices.h"


#ifdef __cplusplus
extern "C" {
#endif


void tkf91_dp_f(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g,
        int trace_flag,
        int png_flag,
        slong *A, size_t szA,
        slong *B, size_t szB);


#ifdef __cplusplus
}
#endif

#endif
