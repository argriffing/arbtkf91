#ifndef RGENERATORS_H
#define RGENERATORS_H

#include "flint/flint.h"
#include "flint/fmpz.h"

#include "expressions.h"


#define RGEN_STATUS_CLOSED 0
#define RGEN_STATUS_OPEN 1
#define RGEN_STATUS_FINALIZED 2

/*
 * Rational factor refining generators.
 */

typedef struct rgen_fmpq_node_struct
{
    fmpq_t value;
    fmpz_t count;
    struct rgen_fmpq_node_struct * next;
} rgen_fmpq_node_struct;
typedef rgen_fmpq_node_struct * rgen_fmpq_node_ptr;
typedef rgen_fmpq_node_struct rgen_fmpq_node_t[1];

typedef struct rgen_expr_node_struct
{
    expr_ptr p;
    fmpz_t count;
    struct rgen_expr_node_struct * next;
} rgen_expr_node_struct;
typedef rgen_expr_node_struct * rgen_expr_node_ptr;
typedef rgen_expr_node_struct rgen_expr_node_t[1];

typedef struct rgen_reg_node_struct
{
    rgen_fmpq_node_ptr p_fmpq_head, p_fmpq_tail;
    rgen_expr_node_ptr p_expr_head, p_expr_tail;
    struct rgen_reg_node_struct * next;
} rgen_reg_node_struct;
typedef rgen_reg_node_struct * rgen_reg_node_ptr;
typedef rgen_reg_node_struct rgen_reg_node_t[1];

/* TODO move this internal structure to the c code... but forward declare it */
/* one of these is created when the registry is 'finalized' */
typedef struct
{
    slong nbases;
    fmpz * base_integers;
    expr_ptr base_expressions;
} rgen_reg_refinement_struct;
typedef rgen_reg_refinement_struct rgen_reg_refinement_t[1];
typedef rgen_reg_refinement_struct * rgen_reg_refinement_ptr;

typedef struct
{
    int status;
    rgen_reg_node_ptr head;
    rgen_reg_node_ptr tail;
    slong size;
    rgen_refinement_ptr refinement;
    reg_ptr reg;
} rgen_reg_struct;
typedef rgen_reg_struct rgen_reg_t[1];
typedef rgen_reg_struct * rgen_reg_ptr;



#ifdef __cplusplus
extern "C" {
#endif

void rgen_fmpq_node_init(rgen_fmpq_node_t p);
void rgen_fmpq_node_clear(rgen_fmpq_node_t p);

void rgen_expr_node_init(rgen_expr_node_t p);
void rgen_expr_node_clear(rgen_expr_node_t p);

void rgen_reg_refinement_init(rgen_reg_refinement_t p);
void rgen_reg_refinement_clear(rgen_reg_refinement_t p);

void rgen_reg_init(rgen_reg_t g);
void rgen_reg_open(rgen_reg_t g, slong *pidx);
void rgen_reg_add_expr(rgen_reg_t g, expr_ptr p, slong count);
void rgen_reg_add_fmpq(rgen_reg_t g, fmpq_t value, slong count);
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
