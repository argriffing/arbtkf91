/*
 */

typedef struct
{
    /* factors related to sequence length equilibrium frequency */
    expr_ptr one_minus_lambda_div_mu;
    expr_ptr lambda_div_mu;

    /* factors related to the indel process involving beta */
    expr_ptr exp_neg_mu_tau;
    expr_ptr one_minus_lambda_beta;
    expr_ptr lambda_beta;
    expr_ptr the_long_beta_expression; /* 1 - exp(-mu*t) - mu*beta */
    expr_ptr mu_beta;

    /* factors related to sequence composition */
    expr_ptr pi[4];

    /* factors related to point substitutions */
    expr_ptr match[4];
    expr_ptr mismatch[4];

} basis_expressions_struct;

typedef basis_expressions_struct * basis_expressions_ptr;



/* convenience functions for building generators */

int gen_add_p0_bar(generator_reg_t g, basis_expressions_ptr p, slong k) {
    gen_add(g, p->mu_beta, k);
}

int gen_add_gamma_0(generator_reg_t g, basis_expressions_ptr p, slong k) {
    gen_add(g, p->one_minus_lambda_div_mu, k);
}

int gen_add_gamma_1(generator_reg_t g, basis_expressions_ptr p, slong k) {
    gen_add(g, p->one_minus_lambda_div_mu, k);
    gen_add(g, p->lambda_div_mu, k);
}

int gen_add_zeta_1(generator_reg_t g, basis_expressions_ptr p, slong k) {
    gen_add(g, p->one_minus_lambda_beta, k);
}

int gen_add_zeta_2(generator_reg_t g, basis_expressions_ptr p, slong k) {
    gen_add(g, p->one_minus_lambda_beta, k);
    gen_add(g, p->lambda_beta, k);
}

int gen_add_p1(generator_reg_t g, basis_expressions_ptr p, slong k) {
    gen_add(g, p->exp_neg_mu_tau, k);
    gen_add(g, p->one_minus_lambda_beta, k);
}

int gen_add_p1_bar(generator_reg_t g, basis_expressions_ptr p, slong k) {
    gen_add(g, p->the_long_beta_expression, k);
    gen_add(g, p->one_minus_lambda_beta, k);
}


typedef struct
{
    int m1_00;
    int m0_10;
    int m0_i0_incr[4];
    int m2_01;
    int m2_0j_incr[4];
    int c0_incr[4];
    int c1_match_incr[4];
    int c1_mismatch_incr[4];
    int c2_incr[4];
} named_generators_struct;

typedef named_generators_struct * named_generators_ptr;

int
tenacious_strict_gt(expr_t a, expr_t b)
{
    arb_t x, y;
    int disjoint;
    int result;

    /* check identity of pointers */
    if (a == b)
    {
        return 0;
    }

    /* evaluate at increasing levels of precision until there is no overlap */
    disjoint = 0;
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

    return result;
}

int create_generators(
        named_generators_ptr x,
        generator_reg_t g,
        basis_expressions_ptr p,
        slong *A, nrows)
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
    slong match_flag;

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
    gen_add(g, p->pi+A[0], 1);
    gen_add_p0_bar(g, p, 1);
    gen_close(g);

    /* M0(i>1, 0) = pi_{A_i} * \bar{x} M0(i-1, 0) */
    /* There are four of these generators, one for each nucleotide. */
    for (i = 0; i < 4; i++)
    {
        gen_open(g, x->m0_i0_incr+i);
        gen_add(g, p->pi+i, 1);
        gen_add_p0_bar(g, p, 1);
        gen_close(g);
    }

    /* M2(0, 1) = gamma_0 * zeta_2 * pi_{B_1} */
    /* Note that this depends on the first character of the second sequence. */
    gen_open(g, &(x->m2_01));
    gen_add_gamma_0(g, p, 1);
    gen_add_zeta_2(g, p, 1);
    gen_add(g, p->pi+B[0], 1);
    gen_close(g);

    /* M2(0, j>1) = pi_{B_i} * \bar{x} M2(0, i-1) */
    /* There are four of these generators, one for each nucleotide. */
    for (j = 0; j < 4; j++)
    {
        gen_open(g, x->m2_0j_incr+j);
        gen_add(g, p->pi+j, 1);
        gen_close(g);
    }

    /* C0 multiplier for the M0(i, j) recursion */
    for (i = 0; i < 4; i++)
    {
        gen_open(g, x->c0_incr+i);
        gen_add(g, p->lamda_div_mu, 1);
        gen_add(g, p->pi+i, 1);
        gen_add_p0_bar(g, p, 1);
        gen_close(g);
    }

    /* C1 multiplier for the M1(i, j) recursion */
    /* Note that this depends on a precomputable argmax. */
    expr_t tmp_lhs, tmp_rhs;
    expr_t tmp_p1, tmp_p1_bar;

    expr_mul(tmp_p1, p->exp_neg_mu_tau, p->one_minus_lambda_beta);
    expr_mul(tmp_p1_bar, p->exp_neg_mu_tau, p->one_minus_lambda_beta);

    for (j = 0; j < 4; j++)
    {
        /*
         * compare
         * P_{ai->bj} * p1
         * vs
         * pi_{bj} * \bar{p1}
         */
        expr_mul(tmp_rhs, p->pi+j, tmp_p1_bar);

        /* generator for matching states */
        gen_open(g, x->c1_match_incr+j);
        expr_mul(tmp_lhs, p->match+j, tmp_p1);
        if (tenacious_strict_gt(lhs, rhs))
        {
            gen_add(g, p->match+j, 1);
            gen_add_p1(g, p, 1);
        }
        else
        {
            gen_add(g, p->pi+j, 1);
            gen_add_p1_bar(g, p, 1);
        }
        expr_clear(tmp_lhs);
        gen_close(g);

        /* generator for mismatched states */
        gen_open(g, x->c1_mismatch_incr+j);
        expr_mul(tmp_lhs, p->mismatch+j, tmp_p1);
        if (tenacious_strict_gt(lhs, rhs))
        {
            gen_add(g, p->mismatch+j, 1);
            gen_add_p1(g, p, 1);
        }
        else
        {
            gen_add(g, p->pi+j, 1);
            gen_add_p1_bar(g, p, 1);
        }
        expr_clear(tmp_lhs);
        gen_close(g);

        expr_clear(tmp_rhs);
    }

    expr_clear(tmp_p1);
    expr_clear(tmp_p1_bar);

    /* C2 multiplier for the M2(i, j) recursion */
    for (i = 0; i < 4; i++)
    {
        gen_open(g, x->c2_incr+i);
        gen_add(g, p->pi+i, 1);
        gen_add(g, p->lamda_beta, 1);
        gen_close(g);
    }
}


void
gen_open(generator_reg_t g, slong *pidx)
{
    ;
}

void
gen_add(generator_reg_t g, expr_ptr expression, slong exponent)
{
    ;
}

void
gen_close(generator_reg_t g)
{
    ;
}
