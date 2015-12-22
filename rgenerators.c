#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "flint/fmpz_mat.h"

#include "factor_refinement.h"
#include "expressions.h"
#include "rgenerators.h"


#define RGEN_STATUS_CLOSED 0
#define RGEN_STATUS_OPEN 1
#define RGEN_STATUS_FINALIZED 2



typedef struct rgen_fmpq_node_struct
{
    fmpq_t value;
    fmpz_t count;
    struct rgen_fmpq_node_struct * next;
} rgen_fmpq_node_struct;
typedef rgen_fmpq_node_struct * rgen_fmpq_node_ptr;
typedef rgen_fmpq_node_struct rgen_fmpq_node_t[1];

void rgen_fmpq_node_init(rgen_fmpq_node_t p);
void rgen_fmpq_node_clear(rgen_fmpq_node_t p);


typedef struct rgen_expr_node_struct
{
    expr_ptr p;
    fmpz_t count;
    struct rgen_expr_node_struct * next;
} rgen_expr_node_struct;
typedef rgen_expr_node_struct * rgen_expr_node_ptr;
typedef rgen_expr_node_struct rgen_expr_node_t[1];

void rgen_expr_node_init(rgen_expr_node_t p);
void rgen_expr_node_clear(rgen_expr_node_t p);


typedef struct rgen_reg_node_struct
{
    rgen_fmpq_node_ptr p_fmpq_head, p_fmpq_tail;
    rgen_expr_node_ptr p_expr_head, p_expr_tail;
    struct rgen_reg_node_struct * next;
} rgen_reg_node_struct;
typedef rgen_reg_node_struct * rgen_reg_node_ptr;
typedef rgen_reg_node_struct rgen_reg_node_t[1];

void rgen_reg_node_init(rgen_reg_node_t p);
void rgen_reg_node_clear(rgen_reg_node_t p);
void rgen_reg_node_add_fmpq(rgen_reg_node_ptr p, fmpq_t value, slong count);
void rgen_reg_node_add_expr(rgen_reg_node_ptr p, expr_ptr expr, slong count);


typedef struct
{
    fmpz_factor_t fac;
    expr_ptr * base_expressions;
} rgen_reg_refinement_struct;
typedef rgen_reg_refinement_struct rgen_reg_refinement_t[1];
typedef rgen_reg_refinement_struct * rgen_reg_refinement_ptr;

void rgen_reg_refinement_init(rgen_reg_refinement_t p);
void rgen_reg_refinement_clear(rgen_reg_refinement_t p);


typedef struct rgen_reg_struct_tag
{
    int status;
    rgen_reg_node_ptr head;
    rgen_reg_node_ptr tail;
    slong size;
    rgen_reg_refinement_ptr refinement;
    reg_ptr reg;
} rgen_reg_struct;

void rgen_reg_assert_status(rgen_reg_ptr g, int status);






void
rgen_fmpq_node_init(rgen_fmpq_node_t p)
{
    fmpq_init(p->value);
    fmpz_init(p->count);
    p->next = NULL;
}

void
rgen_fmpq_node_clear(rgen_fmpq_node_t p)
{
    fmpq_clear(p->value);
    fmpz_clear(p->count);
    p->next = NULL;
}

void
rgen_expr_node_init(rgen_expr_node_t p)
{
    p->p = NULL;
    fmpz_init(p->count);
    p->next = NULL;
}

void
rgen_expr_node_clear(rgen_expr_node_t p)
{
    p->p = NULL;
    fmpz_clear(p->count);
    p->next = NULL;
}




void
rgen_reg_node_add_fmpq(rgen_reg_node_ptr p, fmpq_t value, slong count)
{
    rgen_fmpq_node_ptr node;
    node = flint_malloc(sizeof(rgen_fmpq_node_struct));
    rgen_fmpq_node_init(node);
    fmpq_set(node->value, value);
    fmpz_set_si(node->count, count);

    if (p->p_fmpq_head)
    {
        p->p_fmpq_tail->next = node;
        p->p_fmpq_tail = node;
    }
    else
    {
        p->p_fmpq_head = node;
        p->p_fmpq_tail = node;
    }
}

void
rgen_reg_node_add_expr(rgen_reg_node_ptr p, expr_ptr expr, slong count)
{
    rgen_expr_node_ptr node;
    node = flint_malloc(sizeof(rgen_expr_node_struct));
    rgen_expr_node_init(node);
    node->p = expr;
    fmpz_set_si(node->count, count);

    if (p->p_expr_head)
    {
        p->p_expr_tail->next = node;
        p->p_expr_tail = node;
    }
    else
    {
        p->p_expr_head = node;
        p->p_expr_tail = node;
    }
}

void
rgen_reg_node_init(rgen_reg_node_t p)
{
    p->p_fmpq_head = NULL;
    p->p_fmpq_tail = NULL;
    p->p_expr_head = NULL;
    p->p_expr_tail = NULL;
    p->next = NULL;
}

void
rgen_reg_node_clear(rgen_reg_node_t p)
{
    /* clear and delete all of the fmpq nodes */
    {
        rgen_fmpq_node_ptr node, next;
        node = p->p_fmpq_head;
        while(node)
        {
            next = node->next;
            rgen_fmpq_node_clear(node);
            flint_free(node);
            node = next;
        }
    }

    /* clear and delete all of the expr nodes */
    {
        rgen_expr_node_ptr node, next;
        node = p->p_expr_head;
        while(node)
        {
            next = node->next;
            rgen_expr_node_clear(node);
            flint_free(node);
            node = next;
        }
    }

    p->p_fmpq_head = NULL;
    p->p_fmpq_tail = NULL;
    p->p_expr_head = NULL;
    p->p_expr_tail = NULL;
    p->next = NULL;
}



void
rgen_reg_refinement_init(rgen_reg_refinement_t p)
{
    fmpz_factor_init(p->fac);
    p->base_expressions = NULL;
}

void
rgen_reg_refinement_clear(rgen_reg_refinement_t p)
{
    fmpz_factor_clear(p->fac);
    free(p->base_expressions);
    p->base_expressions = NULL;
}



void rgen_reg_assert_status(rgen_reg_ptr g, int status)
{
    if (g->status != status)
    {
        flint_printf("inappropriate status for this call\n");
        abort();
    }
}

rgen_reg_ptr
rgen_reg_new()
{
    rgen_reg_ptr g;
    g = flint_malloc(sizeof(rgen_reg_struct));
    g->status = RGEN_STATUS_CLOSED;
    g->head = NULL;
    g->tail = NULL;
    g->size = 0;
    g->refinement = NULL;
    g->reg = NULL;
    return g;
}


void rgen_open(rgen_reg_ptr g, slong *pidx)
{
    rgen_reg_assert_status(g, RGEN_STATUS_CLOSED);

    *pidx = g->size;

    rgen_reg_node_ptr node;
    node = flint_malloc(sizeof(rgen_reg_node_struct));
    rgen_reg_node_init(node);

    if (g->head)
    {
        g->tail->next = node;
        g->tail = node;
    }
    else
    {
        g->head = node;
        g->tail = node;
    }

    g->size += 1;
    g->status = RGEN_STATUS_OPEN;
}

void rgen_add_expr(rgen_reg_ptr g, expr_ptr expr, slong count)
{
    rgen_reg_assert_status(g, RGEN_STATUS_OPEN);
    rgen_reg_node_add_expr(g->tail, expr, count);
}

void rgen_add_fmpq(rgen_reg_ptr g, fmpq_t value, slong count)
{
    rgen_reg_assert_status(g, RGEN_STATUS_OPEN);
    rgen_reg_node_add_fmpq(g->tail, value, count);
}

void
rgen_close(rgen_reg_ptr g)
{
    rgen_reg_assert_status(g, RGEN_STATUS_OPEN);
    g->status = RGEN_STATUS_CLOSED;
}


void
rgen_reg_finalize(rgen_reg_ptr g, reg_ptr reg)
{
    rgen_reg_assert_status(g, RGEN_STATUS_CLOSED);

    g->reg = reg;
    g->refinement = flint_malloc(sizeof(rgen_reg_struct));
    rgen_reg_refinement_init(g->refinement);

    int i;
    rgen_reg_node_ptr gnode;
    rgen_fmpq_node_ptr fnode;

    /*
     * Allocate the fmpz vector as input for factor refinement.
     * This will have twice the total number of fmpq nodes;
     * it has an entry for each numerator and each denominator.
     * Semantically, the input for factor refinement will be treated
     * as a set, so there is no need to preserve any kind of ordering.
     */
    fmpz_factor_t unrefined;
    fmpz_factor_init(unrefined);
    unrefined->sign = 1;

    /*
     * Fill the fmpz vector by digging through the layers of node links
     * and scavenging the numerators and denominators.
     */
    for (gnode = g->head; gnode; gnode = gnode->next)
    {
        for (fnode = gnode->p_fmpq_head; fnode; fnode = fnode->next)
        {
            _fmpz_factor_append(unrefined, fmpq_numref(fnode->value), 1);
            _fmpz_factor_append(unrefined, fmpq_denref(fnode->value), 1);
        }
    }

    flint_printf("before refinement:\n");
    fmpz_factor_print(unrefined);
    flint_printf("\n");

    fmpz_factor_refine(g->refinement->fac, unrefined);
    fmpz_factor_clear(unrefined);

    flint_printf("after refinement:\n");
    fmpz_factor_print(g->refinement->fac);
    flint_printf("\n");

    /*
     * Create an array of uninitialized expression pointers
     * which will eventually hold the registered refined factor expressions.
     */
    g->refinement->base_expressions = flint_malloc(
            g->refinement->fac->num * sizeof(expr_ptr));

    /* Register the refined factors as new expressions. */
    expr_ptr expr;
    for (i = 0; i < g->refinement->fac->num; i++)
    {
        expr = reg_new(reg);
        expr_fmpz(expr, g->refinement->fac->p+i);
        g->refinement->base_expressions[i] = expr;
    }

    g->status = RGEN_STATUS_FINALIZED;
}

slong
rgen_reg_nrows(rgen_reg_ptr g)
{
    rgen_reg_assert_status(g, RGEN_STATUS_FINALIZED);
    return g->size;
}

slong
rgen_reg_ncols(rgen_reg_ptr g)
{
    rgen_reg_assert_status(g, RGEN_STATUS_FINALIZED);
    return g->reg->size;
}


void
rgen_reg_get_matrix(fmpz_mat_t mat, rgen_reg_ptr g)
{
    rgen_reg_assert_status(g, RGEN_STATUS_FINALIZED);

    if (fmpz_mat_nrows(mat) != rgen_reg_nrows(g) ||
        fmpz_mat_ncols(mat) != rgen_reg_ncols(g))
    {
        flint_printf("incompatible matrix dimensions\n");
        abort();
    }

    fmpz_mat_zero(mat);

    slong row, col;
    rgen_reg_node_ptr gnode;
    fmpz * entry;
    reg_node_ptr r;
    
    /*
     * Add entries that correspond to expressions
     * that are not associated with factor refinements of rational numbers.
     */
    rgen_expr_node_ptr enode;
    row = 0;
    for (gnode = g->head; gnode; gnode = gnode->next)
    {
        for (enode = gnode->p_expr_head; enode; enode = enode->next)
        {
            r = enode->p->userdata;
            col = r->index;
            entry = fmpz_mat_entry(mat, row, col);
            fmpz_add(entry, entry, enode->count);
        }
        row++;
    }

    /*
     * Add entries using the generator registry,
     * the expressions registry, the vector of refined factors,
     * and the vector giving the expressions corresponding to those factors.
     */
    rgen_fmpq_node_ptr fnode;
    int i;
    fmpz_t x;
    fmpz_init(x);
    row = 0;
    for (gnode = g->head; gnode; gnode = gnode->next)
    {
        for (fnode = gnode->p_fmpq_head; fnode; fnode = fnode->next)
        {
            /* numerator */
            fmpz_set(x, fmpq_numref(fnode->value));
            for (i = 0; i < g->refinement->fac->num; i++)
            {
                while (fmpz_divisible(x, g->refinement->fac->p+i))
                {
                    r = g->refinement->base_expressions[i]->userdata;
                    col = r->index;
                    entry = fmpz_mat_entry(mat, row, col);
                    fmpz_add(entry, entry, fnode->count);
                    fmpz_divexact(x, x, g->refinement->fac->p+i);
                }
            }
            if (!fmpz_is_one(x))
            {
                flint_printf("numerator could not be reduced\n");
                abort();
            }

            /* denominator */
            fmpz_set(x, fmpq_denref(fnode->value));
            for (i = 0; i < g->refinement->fac->num; i++)
            {
                while (fmpz_divisible(x, g->refinement->fac->p+i))
                {
                    r = g->refinement->base_expressions[i]->userdata;
                    col = r->index;
                    entry = fmpz_mat_entry(mat, row, col);
                    fmpz_sub(entry, entry, fnode->count);
                    fmpz_divexact(x, x, g->refinement->fac->p+i);
                }
            }
            if (!fmpz_is_one(x))
            {
                flint_printf("denominator could not be reduced\n");
                abort();
            }
        }
        row++;
    }
    fmpz_clear(x);
}

void
rgen_reg_clear(rgen_reg_ptr g)
{
    rgen_reg_assert_status(g, RGEN_STATUS_FINALIZED);

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
    flint_free(g->refinement);

    flint_free(g);
}
