#include "flint/flint.h"
#include "flint/fmpz.h"

#include "expressions.h"
#include "rgenerators.h"


/* TODO remove declaration in this C file*/
typedef struct rgen_fmpq_node_struct
{
    fmpq_t value;
    fmpz_t count;
    struct rgen_fmpq_node_struct * next;
} rgen_fmpq_node_struct;
typedef rgen_fmpq_node_struct * rgen_fmpq_node_ptr;
typedef rgen_fmpq_node_struct rgen_fmpq_node_t[1];

void
rgen_fmpq_node_init(rgen_fmpq_node_t p)
{
    fmpq_init(p->value);
    fmpz_init(p->count);
    p->next = NULL;
}

void rgen_fmpq_node_clear(rgen_fmpq_node_t p)
{
    fmpq_clear(p->value);
    fmpz_clear(p->count);
    p->next = NULL;
}



/* TODO remove declaration in this C file*/
typedef struct rgen_expr_node_struct
{
    expr_ptr p;
    fmpz_t count;
    struct rgen_expr_node_struct * next;
} rgen_expr_node_struct;
typedef rgen_expr_node_struct * rgen_expr_node_ptr;
typedef rgen_expr_node_struct rgen_expr_node_t[1];

void
rgen_expr_node_init(rgen_expr_node_t p)
{
    fmpz_init(p->count);
    p->next = NULL;
}

void rgen_expr_node_clear(rgen_expr_node_t p)
{
    fmpz_clear(p->count);
    p->next = NULL;
}





/* TODO remove declaration in this C file*/
typedef struct rgen_reg_node_struct
{
    rgen_fmpq_node_ptr p_fmpq_head, p_fmpq_tail;
    rgen_expr_node_ptr p_expr_head, p_expr_tail;
    struct rgen_reg_node_struct * next;
} rgen_reg_node_struct;
typedef rgen_reg_node_struct * rgen_reg_node_ptr;
typedef rgen_reg_node_struct rgen_reg_node_t[1];


void
rgen_expr_node_init(rgen_expr_node_t p)
{
    p->p_fmpq_head = NULL;
    p->p_fmpq_tail = NULL;
    p->p_expr_head = NULL;
    p->p_expr_tail = NULL;
    p->next = NULL;
}

void
rgen_expr_node_clear(rgen_expr_node_t p)
{
    p->p_fmpq_head = NULL;
    p->p_fmpq_tail = NULL;
    p->p_expr_head = NULL;
    p->p_expr_tail = NULL;
    p->next = NULL;
}



/* TODO struct eventually should be defined in .c but declared in .h */
typedef struct
{
    slong nbases;
    fmpz * base_integers;
    expr_ptr base_expressions;
} rgen_reg_refinement_struct;
typedef rgen_reg_refinement_struct rgen_reg_refinement_t[1];
typedef rgen_reg_refinement_struct * rgen_reg_refinement_ptr;

void
rgen_reg_refinement_init(rgen_reg_refinement_t p)
{
    p->nbases = 0;
    p->base_integers = NULL;
    p->base_expressions = NULL;
}

void
rgen_reg_refinement_clear(rgen_reg_refinement_t p)
{
    _fmpz_vec_clear(p->base_integers, p->nbases);
    free(p->base_expressions);
    p->base_integers = NULL;
    p->base_expressions = NULL;
    p->nbases = 0;
}



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



void rgen_reg_init(rgen_reg_t g);
void rgen_reg_open(rgen_reg_t g, slong *pidx);
void rgen_reg_add_expr(rgen_reg_t g, expr_ptr p, slong count);
void rgen_reg_add_fmpq(rgen_reg_t g, fmpq_t value, slong count);
void rgen_reg_close(rgen_reg_t g);

void
rgen_reg_finalize(rgen_reg_t g, reg_ptr reg)
{
    if (g->status != RGEN_STATUS_CLOSED)
    {
        flint_printf("inappropriate status for this call\n");
        abort();
    }

    g->reg = reg;
    g->refinement = flint_malloc(sizeof(rgen_reg_struct));
    rgen_reg_refinement_init(g->refinement);


    /* TODO */

    g->status != RGEN_STATUS_FINALIZED;
}

slong
rgen_reg_nrows(rgen_reg_t g)
{
    if (g->status != RGEN_STATUS_FINALIZED)
    {
        flint_printf("inappropriate status for this call\n");
        abort();
    }
    return g->size;
}

slong
rgen_reg_ncols(rgen_reg_t g)
{
    if (g->status != RGEN_STATUS_FINALIZED)
    {
        flint_printf("inappropriate status for this call\n");
        abort();
    }
    return g->reg->size;
}


void
rgen_reg_get_matrix(fmpz_mat_t mat, rgen_reg_t g)
{
    if (g->status != RGEN_STATUS_FINALIZED)
    {
        flint_printf("inappropriate status for this call\n");
        abort();
    }
    if (fmpz_mat_nrows(mat) != rgen_reg_nrows(g) ||
        fmpz_mat_ncols(mat) != rgen_reg_ncols(g))
    {
        flint_printf("incompatible matrix dimensions\n");
        abort();
    }

}

void
rgen_reg_clear(rgen_reg_t g)
{
    if (g->status != RGEN_STATUS_FINALIZED)
    {
        flint_printf("inappropriate status for this call\n");
        abort();
    }

    /* clear and delete all of the nodes in the registry */
    rgen_reg_node_ptr node, next;
    node = g->head;
    while(node)
    {
        next = node->next;
        rgen_reg_node_clear(node);
        flint_free(node);
        node = next;
    }

    /* clear the refinement */
    rgen_reg_refinement_clear(g->refinement);
}
