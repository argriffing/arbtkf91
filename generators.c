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


/* These are convenience functions. */
void gen_open(generator_reg_t g, slong *pidx);
void gen_add(generator_reg_t g, expr_ptr expression, slong exponent);
void gen_close(generator_reg_t g);


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


/*
 * Convenience functions for building generators.
 * These are to be used only internally to this module.
 */

void gen_add_p0_bar(generator_reg_t g, tkf91_expressions_ptr p, slong k);
void gen_add_gamma_0(generator_reg_t g, tkf91_expressions_ptr p, slong k);
void gen_add_gamma_1(generator_reg_t g, tkf91_expressions_ptr p, slong k);
void gen_add_zeta_1(generator_reg_t g, tkf91_expressions_ptr p, slong k);
void gen_add_zeta_2(generator_reg_t g, tkf91_expressions_ptr p, slong k);
void gen_add_p1(generator_reg_t g, tkf91_expressions_ptr p, slong k);
void gen_add_p1_bar(generator_reg_t g, tkf91_expressions_ptr p, slong k);

void gen_add_p0_bar(generator_reg_t g, tkf91_expressions_ptr p, slong k) {
    gen_add(g, p->mu_beta, k);
}

void gen_add_gamma_0(generator_reg_t g, tkf91_expressions_ptr p, slong k) {
    gen_add(g, p->one_minus_lambda_div_mu, k);
}

void gen_add_gamma_1(generator_reg_t g, tkf91_expressions_ptr p, slong k) {
    gen_add(g, p->one_minus_lambda_div_mu, k);
    gen_add(g, p->lambda_div_mu, k);
}

void gen_add_zeta_1(generator_reg_t g, tkf91_expressions_ptr p, slong k) {
    gen_add(g, p->one_minus_lambda_beta, k);
}

void gen_add_zeta_2(generator_reg_t g, tkf91_expressions_ptr p, slong k) {
    gen_add(g, p->one_minus_lambda_beta, k);
    gen_add(g, p->lambda_beta, k);
}

void gen_add_p1(generator_reg_t g, tkf91_expressions_ptr p, slong k) {
    gen_add(g, p->exp_neg_mu_tau, k);
    gen_add(g, p->one_minus_lambda_beta, k);
}

void gen_add_p1_bar(generator_reg_t g, tkf91_expressions_ptr p, slong k) {
    gen_add(g, p->the_long_beta_expression, k);
    gen_add(g, p->one_minus_lambda_beta, k);
}





/* This convenience function is used only internally to this module. */
int tenacious_strict_gt(expr_t a, expr_t b);

int
tenacious_strict_gt(expr_t a, expr_t b)
{
    arb_t x, y;
    int disjoint;
    int result;
    slong level;

    /* check identity of pointers */
    if (a == b)
    {
        return 0;
    }

    /* evaluate at increasing levels of precision until there is no overlap */
    disjoint = 0;
    result = 0;
    for (level = 0; level < EXPR_CACHE_CAP && !disjoint; level++)
    {
        arb_init(x);
        arb_init(y);
        expr_eval(x, a, level);
        expr_eval(y, b, level);
        if (!arb_overlaps(x, y))
        {
            disjoint = 1;
            result = arb_gt(x, y);
        }
        arb_clear(x);
        arb_clear(y);
    }
    if (!disjoint)
    {
        flint_printf("tenacious strict comparison failed\n");
        abort();
    }
    /*
    else
    {
        flint_printf("detected inequality at level %wd\n", level);
    }
    */

    return result;
}


void
named_generators_init(
        named_generators_t x,
        generator_reg_t g,
        tkf91_expressions_t p,
        slong *A, slong Alen,
        slong *B, slong Blen)
{
    /*
     * The linear integer combinations defining the generators themselves
     * depend on the inputs to some extent.
     * For example the initialization depends on initial sequence characters,
     * and part of the recursion depends on a comparison
     * between an expression related to substitution transition probabilities
     * and an expression related to the birth/death indel parameters.
     */
    slong i, j;

    if (!Alen || !Blen)
    {
        flint_printf("expected both sequences to have length at least 1\n");
        abort();
    }

    /* M1(0, 0) = gamma_0 * zeta_1 */
    gen_open(g, &(x->m1_00));
    gen_add_gamma_0(g, p, 1);
    gen_add_zeta_1(g, p, 1);
    gen_close(g);

    /* M0(1, 0) = gamma_1 * zeta_1 * pi_{A_1} * \bar{p0} */
    /* Note that this depends on the first character of the first sequence. */
    gen_open(g, &(x->m0_10));
    gen_add_gamma_1(g, p, 1);
    gen_add_zeta_1(g, p, 1);
    gen_add(g, p->pi[A[0]], 1);
    gen_add_p0_bar(g, p, 1);
    gen_close(g);

    /* M0(i>1, 0) = pi_{A_i} * (...) * M0(i-1, 0) */
    /* There are four of these generators, one for each nucleotide. */
    for (i = 0; i < 4; i++)
    {
        gen_open(g, x->m0_i0_incr+i);
        gen_add(g, p->lambda_div_mu, 1); /* contribution from gamma */
        gen_add(g, p->lambda_beta, 1); /* contribution from zeta */
        gen_add(g, p->pi[i], 1);
        gen_add_p0_bar(g, p, 1);
        gen_close(g);
    }

    /* M2(0, 1) = gamma_0 * zeta_2 * pi_{B_1} */
    /* Note that this depends on the first character of the second sequence. */
    gen_open(g, &(x->m2_01));
    gen_add_gamma_0(g, p, 1);
    gen_add_zeta_2(g, p, 1);
    gen_add(g, p->pi[B[0]], 1);
    gen_close(g);

    /* M2(0, j>1) = pi_{B_i} * (...) * M2(0, i-1) */
    /* There are four of these generators, one for each nucleotide. */
    for (j = 0; j < 4; j++)
    {
        gen_open(g, x->m2_0j_incr+j);
        gen_add(g, p->lambda_beta, 1); /* contribution from zeta */
        gen_add(g, p->pi[j], 1);
        gen_close(g);
    }

    /* C0 multiplier for the M0(i, j) recursion */
    for (i = 0; i < 4; i++)
    {
        gen_open(g, x->c0_incr+i);
        gen_add(g, p->lambda_div_mu, 1);
        gen_add(g, p->pi[i], 1);
        gen_add_p0_bar(g, p, 1);
        gen_close(g);
    }

    /* C1 multiplier for the M1(i, j) recursion */
    /* Note that this depends on a precomputable argmax. */
    expr_t tmp_lhs, tmp_rhs;
    expr_t tmp_p1, tmp_p1_bar;

    expr_mul(tmp_p1, p->exp_neg_mu_tau, p->one_minus_lambda_beta);
    expr_mul(tmp_p1_bar, p->the_long_beta_expression, p->one_minus_lambda_beta);

    for (j = 0; j < 4; j++)
    {
        expr_mul(tmp_rhs, p->pi[j], tmp_p1_bar);

        for (i = 0; i < 4; i++)
        {
            gen_open(g, x->c1_incr+i*4+j);

            gen_add(g, p->lambda_div_mu, 1);
            gen_add(g, p->pi[i], 1);

            /*
             * compare
             * P_{ai->bj} * p1
             * vs
             * pi_{bj} * \bar{p1}
             */

            if (i == j)
            {
                /* generator for matching states */
                expr_mul(tmp_lhs, p->match[j], tmp_p1);
                if (tenacious_strict_gt(tmp_lhs, tmp_rhs))
                {
                    gen_add(g, p->match[j], 1);
                    gen_add_p1(g, p, 1);
                }
                else
                {
                    gen_add(g, p->pi[j], 1);
                    gen_add_p1_bar(g, p, 1);
                }
            }
            else
            {
                /* generator for mismatched states */
                expr_mul(tmp_lhs, p->mismatch[j], tmp_p1);
                if (tenacious_strict_gt(tmp_lhs, tmp_rhs))
                {
                    /* gen_add(g, p->mismatch[j], 1); */
                    gen_add(g, p->pi[j], 1);
                    gen_add(g, p->one_minus_exp_negdt, 1);
                    gen_add_p1(g, p, 1);
                }
                else
                {
                    gen_add(g, p->pi[j], 1);
                    gen_add_p1_bar(g, p, 1);
                }
            }
            expr_clear(tmp_lhs);
            gen_close(g);
        }
        expr_clear(tmp_rhs);
    }

    expr_clear(tmp_p1);
    expr_clear(tmp_p1_bar);

    /* C2 multiplier for the M2(i, j) recursion */
    for (i = 0; i < 4; i++)
    {
        gen_open(g, x->c2_incr+i);
        gen_add(g, p->pi[i], 1);
        gen_add(g, p->lambda_beta, 1);
        gen_close(g);
    }
}

void
named_generators_clear(named_generators_t x)
{
    UNUSED(x);
}
