#include "flint/flint.h"
#include "flint/fmpq.h"

#include "tkf91_rationals.h"


void
tkf91_rationals_init(tkf91_rationals_t r,
        const fmpq_t lambda, const fmpq_t mu, const fmpq_t tau,
        const fmpq * pi)
{
    slong i;

    fmpq_t one;

    fmpq_init(one);
    fmpq_one(one);

    fmpq_init(r->lambda);
    fmpq_set(r->lambda, lambda);

    fmpq_init(r->mu);
    fmpq_set(r->mu, mu);

    fmpq_init(r->tau);
    fmpq_set(r->tau, tau);

    /* initialize pi */
    for (i = 0; i < 4; i++)
    {
        fmpq_init(r->pi+i);
        fmpq_set(r->pi+i, pi+i);
    }

    /* initialize qi */
    for (i = 0; i < 4; i++)
    {
        fmpq_init(r->qi+i);
        fmpq_sub(r->qi+i, one, pi+i);
    }

    /* initialize negdt */
    {
        fmpq_t dt;
        fmpq_init(dt);
        fmpq_init(r->negdt);
        fmpq_one(dt);
        for (i = 0; i < 4; i++)
        {
            fmpq_submul(dt, pi+i, pi+i);
        }
        fmpq_div(dt, tau, dt);
        fmpq_neg(r->negdt, dt);
        fmpq_clear(dt);
    }

    /* initialize rational values related to tkf91 gamma */
    {
        fmpq_init(r->lambda_div_mu);
        fmpq_init(r->one_minus_lambda_div_mu);
        fmpq_div(r->lambda_div_mu, lambda, mu);
        fmpq_sub(r->one_minus_lambda_div_mu, one, r->lambda_div_mu);
    }

    /* initialize rational values related to tkf91 beta */
    {
        fmpq_t lambda_minus_mu;
        fmpq_init(lambda_minus_mu);
        fmpq_sub(lambda_minus_mu, lambda, mu);
        fmpq_init(r->beta_exponent);
        fmpq_mul(r->beta_exponent, lambda_minus_mu, tau);
        fmpq_clear(lambda_minus_mu);
    }

    /* -mu*tau */
    {
        fmpq_t mu_tau;
        fmpq_init(mu_tau);
        fmpq_mul(mu_tau, mu, tau);
        fmpq_init(r->neg_mu_tau);
        fmpq_neg(r->neg_mu_tau, mu_tau);
        fmpq_clear(mu_tau);
    }

    fmpq_clear(one);
}


void
tkf91_rationals_clear(tkf91_rationals_t r)
{
    slong i;
    for (i = 0; i < 4; i++)
    {
        fmpq_clear(r->pi+i);
        fmpq_clear(r->qi+i);
    }
    fmpq_clear(r->lambda);
    fmpq_clear(r->mu);
    fmpq_clear(r->tau);
    fmpq_clear(r->negdt);
    fmpq_clear(r->lambda_div_mu);
    fmpq_clear(r->one_minus_lambda_div_mu);
    fmpq_clear(r->beta_exponent);
    fmpq_clear(r->neg_mu_tau);
}
