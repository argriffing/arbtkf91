#include "flint/fmpz_vec.h"
#include "flint/fmpz_mat.h"
#include "expressions.h"
#include "generators.h"


/*
 * These functions are used only internally to this module.
 */
void generator_reg_node_init(generator_reg_node_t x, slong index,
        slong expressions_len);
void generator_reg_node_clear(generator_reg_node_t x,
        slong expressions_len);


void
generator_reg_node_init(generator_reg_node_t x, slong index,
        slong expressions_len)
{
    x->index = index;
    x->p = _fmpz_vec_init(expressions_len);
    x->next = NULL;
}

void
generator_reg_node_clear(generator_reg_node_t x,
        slong expressions_len)
{
    x->index = -1;
    _fmpz_vec_clear(x->p, expressions_len);
    x->p = NULL;
    x->next = NULL;
}


void
generator_reg_init(generator_reg_ptr x, slong expressions_len)
{
    x->head = NULL;
    x->tail = NULL;
    x->size = 0;
    x->expressions_len = expressions_len;
    x->open_for_adding = 0;
}

slong
generator_reg_expressions_len(generator_reg_ptr x)
{
    return x->expressions_len;
}

slong
generator_reg_generators_len(generator_reg_ptr x)
{
    if (x->open_for_adding)
    {
        flint_printf("the final number of generators is not available ");
        flint_printf("while the generator registry is open for adding\n");
        abort();
    }

    return x->size;
}

void
generator_reg_get_matrix(fmpz_mat_t mat, generator_reg_ptr x)
{
    if (x->open_for_adding)
    {
        flint_printf("the matrix is not available ");
        flint_printf("while the generator registry is open for adding\n");
        abort();
    }

    if (fmpz_mat_nrows(mat) != generator_reg_generators_len(x) ||
        fmpz_mat_ncols(mat) != generator_reg_expressions_len(x))
    {
        flint_printf("the dimensions of the generator matrix are wrong\n");
        abort();
    }

    slong i, j;
    fmpz_mat_zero(mat);
    generator_reg_node_ptr node;
    i = 0;
    for (node = x->head; node; node = node->next)
    {
        for (j = 0; j < fmpz_mat_ncols(mat); j++)
        {
            fmpz_set(fmpz_mat_entry(mat, i, j), node->p+j);
        }
        i += 1;
    }
}

void
generator_reg_clear(generator_reg_ptr x)
{
    if (x->open_for_adding)
    {
        flint_printf("the generator cannot be cleared ");
        flint_printf("while it is open for adding\n");
        abort();
    }

    /* clear and delete all of the nodes in the registry */
    generator_reg_node_ptr node, next;
    node = x->head;
    while(node)
    {
        next = node->next;
        generator_reg_node_clear(node, generator_reg_expressions_len(x));
        flint_free(node);
        node = next;
    }

    x->head = NULL;
    x->tail = NULL;
    x->size = 0;
}



void
gen_open(generator_reg_t g, slong *pidx)
{
    if (g->open_for_adding)
    {
        flint_printf("the generator is already open for adding\n");
        abort();
    }

    (*pidx) = g->size;
    size_t node_size = sizeof(generator_reg_node_struct);
    generator_reg_node_ptr node = flint_malloc(node_size);
    generator_reg_node_init(node, g->size, g->expressions_len);
    if (g->size)
    {
        g->tail->next = node;
    }
    else
    {
        g->head = node;
    }
    g->tail = node;
    g->size += 1;
    g->open_for_adding = 1;
}

void
gen_add(generator_reg_t g, expr_ptr expression, slong exponent)
{
    if (!g->open_for_adding)
    {
        flint_printf("the generator is not open for adding\n");
        abort();
    }

    /*
     * Trace back from the expression to the corresponding node
     * in the registry of expressions.
     * Then get the expression index associated with that node.
     */
    reg_node_ptr expression_registry_node = expression->userdata;
    slong idx = expression_registry_node->index;

    /*
     * Find the open generator in the generator registry.
     * The open generator has a vector of expression multiplicities.
     * Update that vector by adding the provided exponent to the
     * entry of the vector corresponding to the multiplicty of the
     * provided expression.
     */
    generator_reg_node_ptr active_node = g->tail;
    fmpz * multiplicity = active_node->p+idx;
    fmpz_add_si(multiplicity, multiplicity, exponent);
}

void
gen_close(generator_reg_t g)
{
    if (!g->open_for_adding)
    {
        flint_printf("attempting to close ");
        flint_printf("while the generator is not open for adding\n");
        abort();
    }
    g->open_for_adding = 0;
}
