#include "expressions.h"


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


/* the rational intermediates struct should be used only within the module. */
typedef struct
{
    fmpq qi[4];
    fmpq_t negdt;
    fmpq_t lambda_div_mu;
    fmpq_t one_minus_lambda_div_mu;
    fmpq_t beta_exponent; /* (lambda - mu) * tau */
    fmpq_t neg_mu_tau;
} rational_intermediates_struct;

typedef rational_intermediates_struct rational_intermediates_t[1];

void rational_intermediates_clear(rational_intermediates_t q);
void rational_intermediates_init(rational_intermediates_t q,
        const fmpq *lambda, const fmpq *mu, const fmpq *tau, const fmpq *pi);


void
rational_intermediates_init(rational_intermediates_t q,
        const fmpq *lambda, const fmpq *mu, const fmpq *tau, const fmpq *pi)
{
    slong i;

    fmpq_t one;

    fmpq_init(one);
    fmpq_one(one);

    /* initialize qi */
    {
        for (i = 0; i < 4; i++)
        {
            fmpq_init(q->qi+i);
            fmpq_sub(q->qi+i, one, pi+i);
        }
    }

    /* initialize negdt */
    {
        fmpq_t dt;
        fmpq_init(dt);
        fmpq_init(q->negdt);
        fmpq_one(dt);
        for (i = 0; i < 4; i++)
        {
            fmpq_submul(dt, pi+i, pi+i);
        }
        fmpq_div(dt, tau, dt);
        fmpq_neg(q->negdt, dt);
        fmpq_clear(dt);
    }

    /* initialize rational values related to tkf91 gamma */
    {
        fmpq_init(q->lambda_div_mu);
        fmpq_init(q->one_minus_lambda_div_mu);
        fmpq_div(q->lambda_div_mu, lambda, mu);
        fmpq_sub(q->one_minus_lambda_div_mu, one, q->lambda_div_mu);
    }

    /* initialize rational values related to tkf91 beta */
    {
        fmpq_t lambda_minus_mu;
        fmpq_init(lambda_minus_mu);
        fmpq_sub(lambda_minus_mu, lambda, mu);
        fmpq_init(q->beta_exponent);
        fmpq_mul(q->beta_exponent, lambda_minus_mu, tau);
        fmpq_clear(lambda_minus_mu);
    }

    /* -mu*tau */
    {
        fmpq_t mu_tau;
        fmpq_init(mu_tau);
        fmpq_mul(mu_tau, mu, tau);
        fmpq_init(q->neg_mu_tau);
        fmpq_neg(q->neg_mu_tau, mu_tau);
        fmpq_clear(mu_tau);
    }

    fmpq_clear(one);
}

void
rational_intermediates_clear(rational_intermediates_t q)
{
    slong i;
    for (i = 0; i < 4; i++)
    {
        fmpq_clear(q->qi+i);
    }
    fmpq_clear(q->negdt);
    fmpq_clear(q->lambda_div_mu);
    fmpq_clear(q->one_minus_lambda_div_mu);
    fmpq_clear(q->beta_exponent);
    fmpq_clear(q->neg_mu_tau);
}



void
tkf91_expressions_init(
        tkf91_expressions_t p,
        reg_t reg,
        const fmpq * lambda,
        const fmpq * mu,
        const fmpq * tau,
        const fmpq * pi)
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

    rational_intermediates_t q;
    rational_intermediates_init(q, lambda, mu, tau, pi);

    /* factors related to sequence length equilibrium frequency */
    {
        p->lambda_div_mu = reg_new(reg);
        expr_fmpq(p->lambda_div_mu, q->lambda_div_mu);

        p->one_minus_lambda_div_mu = reg_new(reg);
        expr_fmpq(p->one_minus_lambda_div_mu, q->one_minus_lambda_div_mu);
    }

    /* factors related to sequence composition */
    {
        expr_ptr alias;
        for (i = 0; i < 4; i++)
        {
            alias = NULL;
            for (j = 0; j < i; j++)
            {
                if (fmpq_equal(pi+i, pi+j))
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
                expr_fmpq(p->pi[i], pi+i);
            }
        }
    }

    /* factors related to the indel process involving beta */
    {
        expr_ptr x_mu = reg_new(reg);
        expr_fmpq(x_mu, mu);

        expr_ptr x_lambda = reg_new(reg);
        expr_fmpq(x_lambda, lambda);

        expr_ptr a = reg_new(reg);
        expr_exp_fmpq(a, q->beta_exponent);

        expr_ptr num = reg_new(reg);
        expr_complement(num, a);

        expr_ptr b = reg_new(reg);
        expr_mul(b, x_lambda, a);

        expr_ptr den = reg_new(reg);
        expr_sub(den, x_mu, b);

        p->exp_neg_mu_tau = reg_new(reg);
        expr_exp_fmpq(p->exp_neg_mu_tau, q->neg_mu_tau);

        expr_ptr beta = reg_new(reg);
        expr_div(beta, num, den);

        p->lambda_beta = reg_new(reg);
        expr_mul(p->lambda_beta, x_lambda, beta);

        p->one_minus_lambda_beta = reg_new(reg);
        expr_complement(p->one_minus_lambda_beta, p->lambda_beta);

        p->mu_beta = reg_new(reg);
        expr_mul(p->mu_beta, x_mu, beta);

        expr_ptr c = reg_new(reg);
        expr_add(c, p->exp_neg_mu_tau, p->mu_beta);

        p->the_long_beta_expression = reg_new(reg);
        expr_complement(p->the_long_beta_expression, c);
    }

    /* factors related to point substitutions */
    {
        expr_ptr exp_negdt = reg_new(reg);
        expr_exp_fmpq(exp_negdt, q->negdt);

        expr_ptr one_minus_exp_negdt = reg_new(reg);
        expr_complement(one_minus_exp_negdt, exp_negdt);

        int alias_found;
        expr_ptr match_alias, mismatch_alias;
        for (i = 0; i < 4; i++)
        {
            alias_found = 0;
            match_alias = NULL;
            mismatch_alias = NULL;
            for (j = 0; j < i; j++)
            {
                if (fmpq_equal(pi+i, pi+j))
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
                expr_mul(a, p->pi[i], one_minus_exp_negdt);

                p->match[i] = reg_new(reg);
                expr_add(p->match[i], exp_negdt, a);

                p->mismatch[i] = a;
            }
        }
    }

    rational_intermediates_clear(q);
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
