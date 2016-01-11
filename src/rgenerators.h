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


rgen_reg_ptr rgen_reg_new();

void rgen_open(rgen_reg_ptr g, slong *pidx);
void rgen_add_fmpq(rgen_reg_ptr g, fmpq_t value, slong count);
void rgen_add_expr(rgen_reg_ptr g, expr_ptr p, slong count);
void rgen_close(rgen_reg_ptr g);

void rgen_reg_finalize(rgen_reg_ptr g, reg_ptr reg);
slong rgen_reg_nrows(rgen_reg_ptr g);
slong rgen_reg_ncols(rgen_reg_ptr g);
void rgen_reg_get_matrix(fmpz_mat_t mat, rgen_reg_ptr g);
void rgen_reg_clear(rgen_reg_ptr g);


#ifdef __cplusplus
}
#endif

#endif
