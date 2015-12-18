#ifndef RGENERATORS_H
#define RGENERATORS_H

#include "flint/flint.h"
#include "flint/fmpz.h"

#include "expressions.h"


struct rgen_reg_struct_tag;
typedef struct rgen_reg_struct_tag * rgen_reg_ptr;


#ifdef __cplusplus
extern "C" {
#endif

void rgen_reg_init(rgen_reg_t g);
void rgen_reg_open(rgen_reg_t g, slong *pidx);
void rgen_reg_add_fmpq(rgen_reg_t g, fmpq_t value, slong count);
void rgen_reg_add_expr(rgen_reg_t g, expr_ptr p, slong count);
void rgen_reg_close(rgen_reg_t g);
void rgen_reg_finalize(rgen_reg_t g, reg_ptr reg);
slong rgen_reg_nrows(rgen_reg_t g);
slong rgen_reg_ncols(rgen_reg_t g);
void rgen_reg_get_matrix(fmpz_mat_t mat, rgen_reg_t g);
void rgen_reg_clear(rgen_reg_t g);

#ifdef __cplusplus
}
#endif

#endif
