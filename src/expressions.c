#include "expressions.h"
#include "unused.h"


void
reg_node_init(reg_node_t x, slong index, expr_ptr p)
{
    x->index = index;
    x->p = p;
    x->next = NULL;
}

void
reg_node_clear(reg_node_t x)
{
    x->index = -1;
    x->p = NULL;
    x->next = NULL;
}


void
reg_init(reg_ptr x)
{
    x->head = NULL;
    x->tail = NULL;
    x->size = 0;
}

void
reg_clear(reg_ptr x)
{
    slong i;
    reg_node_ptr node, next;
    i = 0;
    node = x->head;
    while (node)
    {
        next = node->next;
        if (node->index != i)
        {
            flint_printf("found a node out of order ");
            flint_printf("when clearing a list of expressions");
            abort();
        }
        if (next == NULL && node != x->tail)
        {
            flint_printf("reached the end of a list ");
            flint_printf("before reaching the tail");
            abort();
        }
        expr_clear(node->p);
        flint_free(node->p);
        reg_node_clear(node);
        flint_free(node);
        node = next;
        i++;
    }
    x->head = NULL;
    x->tail = NULL;
    x->size = 0;
}

expr_ptr
reg_new(reg_ptr x)
{
    expr_ptr expr = flint_malloc(sizeof(expr_struct));
    reg_node_ptr node = flint_malloc(sizeof(reg_node_struct));
    expr->userdata = node;
    reg_node_init(node, x->size, expr);
    if (x->size)
    {
        x->tail->next = node;
    }
    else
    {
        x->head = node;
    }
    x->tail = node;
    x->size += 1;
    return expr;
}

expr_ptr *
reg_vec(reg_ptr x)
{
    slong i;
    expr_ptr * vec = flint_malloc(sizeof(expr_ptr) * x->size);
    reg_node_ptr node;
    i = 0;
    for (node = x->head; node; node = node->next)
    {
        vec[i] = node->p;
        i++;
    }
    return vec;
}


void tkf91_expressions_init(
        tkf91_expressions_ptr p,
        reg_t reg,
        const tkf91_rationals_t r)
{
    /*
     * Create a bunch of static single assignment expressions.
     * Track all of them in an expressions registry;
     * this will let the expressions be cleared later,
     * and it will associate each expression with an array index.
     * Give some of the expressions special names to be used
     * later when the generators are defined.
     */
    slong i, j;

    /* factors related to sequence length equilibrium frequency */
    {
        p->lambda_div_mu = reg_new(reg);
        expr_fmpq(p->lambda_div_mu, r->lambda_div_mu);

        p->one_minus_lambda_div_mu = reg_new(reg);
        expr_fmpq(p->one_minus_lambda_div_mu, r->one_minus_lambda_div_mu);
    }

    /* factors related to sequence composition */
    {
        expr_ptr alias;
        for (i = 0; i < 4; i++)
        {
            alias = NULL;
            for (j = 0; j < i; j++)
            {
                if (fmpq_equal(r->pi+i, r->pi+j))
                {
                    alias = p->pi[j];
                    break;
                }
            }
            if (alias)
            {
                p->pi[i] = alias;
            }
            else
            {
                p->pi[i] = reg_new(reg);
                expr_fmpq(p->pi[i], r->pi+i);
            }
        }
    }

    /* factors related to the indel process involving beta */
    {
        expr_ptr x_mu = reg_new(reg);
        expr_fmpq(x_mu, r->mu);

        expr_ptr x_lambda = reg_new(reg);
        expr_fmpq(x_lambda, r->lambda);

        expr_ptr a = reg_new(reg);
        expr_exp_fmpq(a, r->beta_exponent);

        expr_ptr num = reg_new(reg);
        expr_complement(num, a);

        expr_ptr b = reg_new(reg);
        expr_mul(b, x_lambda, a);

        expr_ptr den = reg_new(reg);
        expr_sub(den, x_mu, b);

        p->exp_neg_mu_tau = reg_new(reg);
        expr_exp_fmpq(p->exp_neg_mu_tau, r->neg_mu_tau);

        p->beta = reg_new(reg);
        expr_div(p->beta, num, den);

        p->lambda_beta = reg_new(reg);
        expr_mul(p->lambda_beta, x_lambda, p->beta);

        p->one_minus_lambda_beta = reg_new(reg);
        expr_complement(p->one_minus_lambda_beta, p->lambda_beta);

        p->mu_beta = reg_new(reg);
        expr_mul(p->mu_beta, x_mu, p->beta);

        expr_ptr c = reg_new(reg);
        expr_add(c, p->exp_neg_mu_tau, p->mu_beta);

        p->the_long_beta_expression = reg_new(reg);
        expr_complement(p->the_long_beta_expression, c);
    }

    /* factors related to point substitutions */
    {
        expr_ptr exp_negdt = reg_new(reg);
        expr_exp_fmpq(exp_negdt, r->negdt);

        p->one_minus_exp_negdt = reg_new(reg);
        expr_complement(p->one_minus_exp_negdt, exp_negdt);

        int alias_found;
        expr_ptr match_alias, mismatch_alias;
        for (i = 0; i < 4; i++)
        {
            alias_found = 0;
            match_alias = NULL;
            mismatch_alias = NULL;
            for (j = 0; j < i; j++)
            {
                if (fmpq_equal(r->pi+i, r->pi+j))
                {
                    match_alias = p->match[j];
                    mismatch_alias = p->mismatch[j];
                    alias_found = 1;
                    break;
                }
            }
            if (alias_found)
            {
                p->match[i] = match_alias;
                p->mismatch[i] = mismatch_alias;
            }
            else
            {
                expr_ptr a = reg_new(reg);
                expr_mul(a, p->pi[i], p->one_minus_exp_negdt);

                p->match[i] = reg_new(reg);
                expr_add(p->match[i], exp_negdt, a);

                p->mismatch[i] = a;
            }
        }
    }
}

void
tkf91_expressions_clear(tkf91_expressions_t p)
{
    /*
     * The names in this structure are only aliases
     * for entries in a registry of expressions.
     */
    UNUSED(p);
}




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
