#ifndef EXPRESSIONS_H
#define EXPRESSIONS_H

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "femtocas.h"


/*
 * This is a node in a registry of mathematical expressions.
 * It is not very clever, and it doesn't allocate or free memory by itself.
 */
typedef struct reg_node_struct
{
    slong index;
    expr_ptr p;
    struct reg_node_struct * next;
} reg_node_struct;
typedef reg_node_struct * reg_node_ptr;
typedef reg_node_struct reg_node_t[1];


/*
 * This structure manages the memory of all of its registry nodes
 * and all of the expression objects in those nodes.
 */
typedef struct
{
    reg_node_ptr head;
    reg_node_ptr tail;
    slong size;
} reg_struct;
typedef reg_struct reg_t[1];
typedef reg_struct * reg_ptr;


/*
 * This is just an aggregate of named expression pointers.
 * The structure does not manage the memory used by the expressions.
 */
typedef struct
{
    /* factors related to sequence length equilibrium frequency */
    expr_ptr one_minus_lambda_div_mu;
    expr_ptr lambda_div_mu;

    /* factors related to sequence composition */
    expr_ptr pi[4];

    /* factors related to the indel process involving beta */
    expr_ptr exp_neg_mu_tau;
    expr_ptr one_minus_lambda_beta;
    expr_ptr lambda_beta;
    expr_ptr the_long_beta_expression; /* 1 - exp(-mu*t) - mu*beta */
    expr_ptr mu_beta;

    /* factors related to point substitutions */
    expr_ptr match[4];
    expr_ptr mismatch[4];

} tkf91_expressions_struct;
typedef tkf91_expressions_struct * tkf91_expressions_ptr;
typedef tkf91_expressions_struct tkf91_expressions_t[1];



#ifdef __cplusplus
extern "C" {
#endif


void reg_node_init(reg_node_t x, slong index, expr_ptr p);
void reg_node_clear(reg_node_t x);


void reg_init(reg_ptr x);
void reg_clear(reg_ptr x);
expr_ptr reg_new(reg_ptr x);
expr_ptr * reg_vec(reg_ptr x);


void tkf91_expressions_init(
        tkf91_expressions_ptr p,
        reg_t reg,
        const fmpq * lambda,
        const fmpq * mu,
        const fmpq * tau,
        const fmpq * pi);

void tkf91_expressions_clear(tkf91_expressions_t);


#ifdef __cplusplus
}
#endif

#endif
