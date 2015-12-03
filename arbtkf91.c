/*
 */

typedef struct
{
    fmpq lambda[1];
    fmpq mu[1];
    fmpq t[1];
    fmpq pi[4];
}
user_params_struct;

typedef user_params_struct user_params_t[1];

void
user_params_init(user_params_t p)
{
    slong i;
    fmpq_init(p->lambda);
    fmpq_init(p->mu);
    fmpq_init(p->t);
    for (i = 0; i < 4; i++)
    {
        fmpq_init(p->pi+i);
    }
}

void
user_params_clear(user_params_t p)
{
    slong i;
    fmpq_clear(p->lambda);
    fmpq_clear(p->mu);
    fmpq_clear(p->t);
    for (i = 0; i < 4; i++)
    {
        fmpq_clear(p->pi+i);
    }
}



void
_arb_log1mexp(arb_t z, const arb_t x, slong prec)
{
    /* log(1 - exp(-x)) : see log1mexp-note.pdf in CRAN */
    if (arf_cmp_2exp_si(arb_midref(x), 0) < 0)
    {
        /* z = log(-expm1(-x)) : x < 1 */
        arb_neg(z, x);
        arb_expm1(z, z, prec);
        arb_neg(z, z);
        arb_log(z, z, prec);
    }
    else
    {
        /* z = log1p(-exp(-x)) : 1 <= x */
        arb_neg(z, x);
        arb_exp(z, z, prec);
        arb_neg(z, z);
        arb_log1p(z, z, prec);
    }
}

void
_arb_max(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    if (arb_overlaps(x, y))
    {
        arf_max(arb_midref(z), arb_midref(x), arb_midref(y));
        mag_max(arb_radref(z), arb_radref(x), arb_radref(y));
    }
    else
    {
        if (arf_cmp(arb_midref(x), arb_midref(y)) < 0)
        {
            arb_set(z, y);
        }
        else
        {
            arb_set(z, x);
        }
    }
}

void
_arb_mul_fmpq(arb_t z, const fmpq_t x, slong prec)
{
    arb_mul_fmpz(z, fmpq_numref(x), prec);
    arb_div_fmpz(z, fmpq_denref(x), prec);
}

void
_arb_addmul_fmpq(arb_t z, const arb_t x, const fmpq_t y, slong prec)
{
    arb_t a;
    arb_init(a);
    arb_init_fmpq(a, y, prec);
    arb_addmul(z, x, a, prec);
    arb_clear(a);
}

void
_arb_submul_fmpq(arb_t z, const arb_t x, const fmpq_t y, slong prec)
{
    arb_t a;
    arb_init(a);
    arb_init_fmpq(a, y, prec);
    arb_submul(z, x, a, prec);
    arb_clear(a);
}



void
_get_trans_element(arb_t elem, ulong i, ulong j, const arb_mat_t logP,
        const arb_t beta, const user_params_t p, slong prec)
{
    arb_t a, b, mut_x, logpj;
    fmpq mut;

    /* mut = mu * t */
    fmpq_init(mut);
    arb_init(mut_x);
    fmpq_mul(mut, p->mu, p->t);
    arb_set_fmpq(mut_x, mut);

    /* a = log(P_ij) - mu*t */
    arb_init(a);
    arb_sub(a, arb_mat_entry(logP, i, j), mut_x, prec);

    /* logpj = log(pi[j]) */
    arb_init(logpj);
    arb_set_fmpq(logpj, p->pi+j, prec);

    /* b = log(pi_j) + log1p(-(exp(-mu*t) + mu*beta)) */
    arb_init(b);
    arb_neg(b, mut_x);
    arb_exp(b, b, prec);
    arb_addmul_fmpq(b, beta, p->mu, prec);
    arb_neg(b, b);
    arb_log1p(b, b, prec);
    arb_add(b, logpj);

    /* elem = max(a, b) */
    _arb_max(elem, a, b);

    arb_clear(a);
    arb_clear(b);
    arb_clear(mut_x);
    arb_clear(logpj);
    fmpq_clear(mut);
}

void
_get_beta(arb_t beta, const user_params p, slong prec)
{
    arb_t a, b, e, ex;
    fmpq_t qe;

    /* e = (lambda - mu)*t */
    fmpq_init(qe);
    arb_init(e);
    fmpq_sub(qe, p->lambda, p->mu);
    fmpq_mul(qe, qe, p->t);
    arb_set_fmpq(e, qe, prec);

    /* a = 1 - exp((lambda - mu)*t) */
    arb_init(a);
    arb_expm1(a, e, prec);
    arb_neg(a);

    /* b = mu - lambda * exp((lambda - mu)*t) */
    arb_init(ex);
    arb_exp(ex, e, prec);
    arb_init(b);
    arb_set_fmpq(b, p->mu, prec);
    arb_submul_fmpq(b, ex, p->lambda, prec);

    /* beta = (1 - exp(...)) / (mu - lambda * exp(...)) */
    arb_div(beta, a, b, prec);

    arb_clean(a);
    arb_clean(b);
    arb_clean(e);
    arb_clean(ex);
    fmpq_clean(qe);
}

void
_get_logbeta(arb_t logbeta, const user_params_t p, slong prec)
{
    arb_t a, b, e, ex;
    fmpq_t qe;

    /* e = (lambda - mu)*t */
    fmpq_init(qe);
    arb_init(e);
    fmpq_sub(qe, p->lambda, p->mu);
    fmpq_mul(qe, qe, p->t);
    arb_set_fmpq(e, qe, prec);

    /* a = log(1 - exp((lambda - mu)*t)) */
    arb_init(a);
    arb_neg(a, e);
    _arb_log1mexp(a, a, prec);

    /* b = log(mu - lambda * exp((lambda - mu)*t)) */
    arb_init(ex);
    arb_exp(ex, e, prec);
    arb_init(b);
    arb_set_fmpq(b, p->mu, prec);
    arb_submul_fmpq(b, ex, lambda, prec);
    arb_log(b);

    /* logbeta = log(1 - exp(...)) - log(mu - lambda * exp(...)) */
    arb_sub(logbeta, a, b, prec);

    arb_clean(a);
    arb_clean(b);
    arb_clean(e);
    arb_clean(ex);
    fmpq_clean(qe);
}

void
_get_dt(fmpq_t dt, const fmpq * pi, const fmpq_t t)
{
    ulong i, n;

    /* number of states */
    n = 4;

    /* dt = t / (1 - pi[0]^2 + pi[1]^2 + ... + pi[n-1]^2) */
    fmpq_one(dt);
    for (i = 0; i < n; i++)
    {
        fmpq_submul(delta, pi[i], pi[i]);
    }
    fmpq_div(dt, t, dt);
}

void
_get_logP(arb_mat_t logP, const fmpq *pi, const fmpq *qi,
          const arb_t dt, slong prec)
{
    ulong i, j, n;
    arb_t a, b, p, q;
    arb_t match, mismatch;

    /* number of states */
    n = 4;

    arb_init(p);
    arb_init(q);
    arb_init(match);
    arb_init(mismatch);

    /* a = exp(-dt) */
    arb_init(a);
    arb_neg(a, dt);
    arb_exp(a, prec);

    /* b = log(1 - exp(-dt)) */
    arb_init(b);
    _arb_log1mexp(b, dt, prec);

    for (j = 0; j < n; j++)
    {
        arb_set_fmpq(p, pi+j, prec);
        arb_set_fmpq(q, qi+j, prec);

        /* match = log(pi[j] + (1-pi[j])*exp(-d*t)) */
        arb_set(match, p);
        arb_addmul(match, q, a, prec);
        arb_log(match, match, prec);

        /* mismatch = log(pi[j]*(1 - exp(-d*t))) */
        arb_log(mismatch, p, prec);
        arb_add(mismatch, mismatch, b, prec);

        for (i = 0; i < n; i++)
        {
            if (i == j)
            {
                arb_set(arb_mat_entry(logP, i, j), match);
            }
            else
            {
                arb_set(arb_mat_entry(logP, i, j), mismatch);
            }
        }
    }
}

void
_init_M0(arb_mat_t M0, ulong n, ulong m, const fmpq * pi)
{
    ulong i;

    /* allocate the matrix and fill with zeros */
    arb_mat_init(M0, n+1, m+1);
    arb_mat_zero(M0);

    /* fill the first column */
    for (i = 0; i < n+1; i++)
    {
        arb_mat_entry();
    }
}

void
_compute_elements(arb_ptr v, const user_params_t p, slong prec)
{
    slong i, h;
    fmpq_t dt;
    arb_mat_t logP;
    arb_t beta, logbeta;
    fmpq_t lambda_div_mu;
    fmpq qi[4];

    /* lambda_div_mu = lambda / mu */
    fmpq_init(lambda_div_mu);
    fmpq_div(lambda_div_mu, p->lambda, p->mu);

    /* define the nucleotide probability complements */
    for (i = 0; i < 4; i++)
    {
        fmpq_init(qi+i);
        fmpq_one(qi+i);
        fmpq_sub(qi+i, qi+i, pi+i);
    }

    /* compute the transition probability matrix entries */
    arb_mat_init(logP);
    fmpq_init(dt);
    _get_dt(dt, pi, t);
    _get_logP(logP, pi, qi, dt, prec);

    /* compute beta and logbeta */
    arb_init(beta);
    arb_init(logbeta);
    _get_beta(beta, lambda, mu, t, slong prec);
    _get_logbeta(logbeta, lambda, mu, t, slong prec);

    /* prepare to serialize the elements to a flat vector */

    h = 0;

    /* serialize a function of each transition matrix entry */
    /* N = 16 */
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            _get_trans_element(v+h, i, j, logP, beta, p, prec);
            h++;
        }
    }

    /* serialize the logarithms of nucleotide probabilities */
    /* N = 4 */
    for (i = 0; i < 4; i++)
    {
        arb_set_fmpq(v+h, pi+i, prec);
        arb_log(v+h, v+h, prec);
        h++;
    }

    /*
     * serialize a few miscellaneous elements:
     * log(lambda)
     * log(mu)
     * log(beta)
     * log1p(-lambda/mu)
     * log1p(-lambda*beta)
     */
    /* N = 5 */

    arb_set_fmpq(v+h, p->lambda, prec);
    arb_log(v+h, v+h, prec);
    h++;

    arb_set_fmpq(v+h, p->mu, prec);
    arb_log(v+h, v+h, prec);
    h++;

    arb_set(v+h, beta);
    h++;

    arb_set_fmpq(v+h, lambda_div_mu, prec);
    arb_neg(v+h, v+h);
    arb_log1p(v+h, v+h, prec);
    h++;

    arb_mul_fmpq(v+h, beta, lambda, prec);
    arb_neg(v+h, v+h);
    arb_log1p(v+h, v+h);
    h++;

    /* clean up */
    fmpq_clean(dt);
    arb_mat_clean(logP);
    arb_clean(beta);
    arb_clean(logbeta);
    fmpq_clean(lambda_div_mu);
    for (i = 0; i < 4; i++)
    {
        fmpq_clean(qi+i);
    }
}

int
main(int argc, char *argc[])
{
    /* use hardcoded parameters but command-line sequences */

    user_params_t p;

    user_params_init(p);

    /* lambda, mu, t = 1, 1/2, 1/10 */
    fmpq_set_si(p->lambda, 1, 1);
    fmpq_set_si(p->mu, 1, 2);
    fmpq_set_si(p->t, 1, 10);

    /* pi = (0.27, 0.24, 0.26, 0.23) */
    fmpq_set_si(p->pi+0, 27, 100);
    fmpq_set_si(p->pi+1, 24, 100);
    fmpq_set_si(p->pi+2, 26, 100);
    fmpq_set_si(p->pi+3, 23, 100);

    n = 10;
    m = 10;

    user_params_clear(p);

    return 0;
}
