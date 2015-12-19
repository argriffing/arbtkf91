#include "flint/fmpz_vec.h"
#include "flint/fmpz_mat.h"

#include "expressions.h"
#include "generators.h"
#include "rgenerators.h"
#include "tkf91_generators.h"
#include "tkf91_rgenerators.h"
#include "tkf91_generator_indices.h"


void rgen_add_p0_bar(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k);
void rgen_add_gamma_0(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k);
void rgen_add_gamma_1(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k);
void rgen_add_zeta_1(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k);
void rgen_add_zeta_2(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k);
void rgen_add_p1(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k);
void rgen_add_p1_bar(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k);



void rgen_add_p0_bar(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k) {
    rgen_add_fmpq(g, r->mu, k);
    rgen_add_expr(g, p->beta, k);
}

void rgen_add_gamma_0(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k) {
    UNUSED(p);
    rgen_add_fmpq(g, r->one_minus_lambda_div_mu, k);
}

void rgen_add_gamma_1(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k) {
    UNUSED(p);
    rgen_add_fmpq(g, r->one_minus_lambda_div_mu, k);
    rgen_add_fmpq(g, r->lambda_div_mu, k);
}

void rgen_add_zeta_1(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k) {
    UNUSED(r);
    rgen_add_expr(g, p->one_minus_lambda_beta, k);
}

void rgen_add_zeta_2(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k) {
    rgen_add_expr(g, p->one_minus_lambda_beta, k);
    rgen_add_fmpq(g, r->lambda, k);
    rgen_add_expr(g, p->beta, k);
}

void rgen_add_p1(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k) {
    UNUSED(r);
    rgen_add_expr(g, p->exp_neg_mu_tau, k);
    rgen_add_expr(g, p->one_minus_lambda_beta, k);
}

void rgen_add_p1_bar(rgen_reg_ptr g,
        tkf91_rationals_ptr r, tkf91_expressions_ptr p, slong k) {
    UNUSED(r);
    rgen_add_expr(g, p->the_long_beta_expression, k);
    rgen_add_expr(g, p->one_minus_lambda_beta, k);
}



void
tkf91_rgenerators_init(
        tkf91_generator_indices_t x,
        rgen_reg_ptr g,
        tkf91_rationals_t r,
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
    rgen_open(g, &(x->m1_00));
    rgen_add_gamma_0(g, r, p, 1);
    rgen_add_zeta_1(g, r, p, 1);
    rgen_close(g);

    /* M0(1, 0) = gamma_1 * zeta_1 * pi_{A_1} * \bar{p0} */
    /* Note that this depends on the first character of the first sequence. */
    rgen_open(g, &(x->m0_10));
    rgen_add_gamma_1(g, r, p, 1);
    rgen_add_zeta_1(g, r, p, 1);
    rgen_add_fmpq(g, r->pi+A[0], 1);
    rgen_add_p0_bar(g, r, p, 1);
    rgen_close(g);

    /* M0(i>1, 0) = pi_{A_i} * (...) * M0(i-1, 0) */
    /* There are four of these generators, one for each nucleotide. */
    for (i = 0; i < 4; i++)
    {
        rgen_open(g, x->m0_i0_incr+i);
        rgen_add_fmpq(g, r->lambda_div_mu, 1); /* contribution from gamma */
        rgen_add_fmpq(g, r->lambda, 1); /* contribution from zeta */
        rgen_add_expr(g, p->beta, 1); /* contribution from zeta */
        rgen_add_fmpq(g, r->pi+i, 1);
        rgen_add_p0_bar(g, r, p, 1);
        rgen_close(g);
    }

    /* M2(0, 1) = gamma_0 * zeta_2 * pi_{B_1} */
    /* Note that this depends on the first character of the second sequence. */
    rgen_open(g, &(x->m2_01));
    rgen_add_gamma_0(g, r, p, 1);
    rgen_add_zeta_2(g, r, p, 1);
    rgen_add_fmpq(g, r->pi+B[0], 1);
    rgen_close(g);

    /* M2(0, j>1) = pi_{B_i} * (...) * M2(0, i-1) */
    /* There are four of these generators, one for each nucleotide. */
    for (j = 0; j < 4; j++)
    {
        rgen_open(g, x->m2_0j_incr+j);
        rgen_add_fmpq(g, r->lambda, 1); /* contribution from zeta */
        rgen_add_expr(g, p->beta, 1); /* contribution from zeta */
        rgen_add_fmpq(g, r->pi+j, 1);
        rgen_close(g);
    }

    /* C0 multiplier for the M0(i, j) recursion */
    for (i = 0; i < 4; i++)
    {
        rgen_open(g, x->c0_incr+i);
        rgen_add_fmpq(g, r->lambda_div_mu, 1);
        rgen_add_fmpq(g, r->pi+i, 1);
        rgen_add_p0_bar(g, r, p, 1);
        rgen_close(g);
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
            rgen_open(g, x->c1_incr+i*4+j);

            rgen_add_fmpq(g, r->lambda_div_mu, 1);
            rgen_add_fmpq(g, r->pi+i, 1);

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
                    rgen_add_expr(g, p->match[j], 1);
                    rgen_add_p1(g, r, p, 1);
                }
                else
                {
                    rgen_add_fmpq(g, r->pi+j, 1);
                    rgen_add_p1_bar(g, r, p, 1);
                }
            }
            else
            {
                /* generator for mismatched states */
                expr_mul(tmp_lhs, p->mismatch[j], tmp_p1);
                if (tenacious_strict_gt(tmp_lhs, tmp_rhs))
                {
                    /* gen_add(g, p->mismatch[j], 1); */
                    rgen_add_fmpq(g, r->pi+j, 1);
                    rgen_add_expr(g, p->one_minus_exp_negdt, 1);
                    rgen_add_p1(g, r, p, 1);
                }
                else
                {
                    rgen_add_fmpq(g, r->pi+j, 1);
                    rgen_add_p1_bar(g, r, p, 1);
                }
            }
            expr_clear(tmp_lhs);
            rgen_close(g);
        }
        expr_clear(tmp_rhs);
    }

    expr_clear(tmp_p1);
    expr_clear(tmp_p1_bar);

    /* C2 multiplier for the M2(i, j) recursion */
    for (i = 0; i < 4; i++)
    {
        rgen_open(g, x->c2_incr+i);
        rgen_add_fmpq(g, r->pi+i, 1);
        rgen_add_fmpq(g, r->lambda, 1);
        rgen_add_expr(g, p->beta, 1);
        rgen_close(g);
    }
}

void
tkf91_rgenerators_clear(tkf91_generator_indices_t x)
{
    UNUSED(x);
}
