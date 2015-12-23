#ifndef TKF91_R_H
#define TKF91_R_H

#include "flint.h"
#include "arb_mat.h"

#include "femtocas.h"
#include "tkf91_generator_indices.h"


#ifdef __cplusplus
extern "C" {
#endif


void
tkf91_dp_r(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g,
        int trace_flag,
        slong *A, size_t szA,
        slong *B, size_t szB);


#ifdef __cplusplus
}
#endif

#endif
