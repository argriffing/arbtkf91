/*
 */

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
_get_log_gamma(arb_ptr * v, ulong n, const fmpq_t birth, const fmpq_t death,
               slong prec);
{
    /* n is the length of the vector v */
    ulong i;
    fmpq_t q;
    arb_t logq;

    if (n <= 0)
    {
        return;
    }

    fmpq_init(q);
    arb_init(logq);

    /* q = birth / death */
    fmpq_div(q, birth, death);

    /* v[0] = log(1 - q) */
    arb_set_fmpq(v+0, q, prec);
    arb_neg(v+0, v+0);
    arb_log1p(v+0, v+0, prec);

    /* v[i] = log(1 - q) + i*log(q) */
    for (i = 1; i < n; i++)
    {
        arb_set(v+i, v+0);
        arb_addmul_ui(v+i, logq, i, prec);
    }

    fmpq_clear(q);
    arb_clear(logq);
}

void
_get_dt(fmpq_t dt, const fmpq * pi, ulong n, const fmpq_t t)
{
    ulong i;

    /* dt = t / (1 - pi[0]^2 + pi[1]^2 + ... + pi[n-1]^2) */
    fmpq_one(dt);
    for (i = 0; i < n; i++)
    {
        fmpq_submul(delta, pi[i], pi[i]);
    }
    fmpq_div(dt, t, dt);
}

void
_get_logP(arb_mat_t logP, const fmpq *pi, const fmpq *qi, ulong n,
          const arb_t dt, slong prec)
{
    /* This function doesn't have to be fast. */
    ulong i, j;
    arb_t a, b, p, q;
    arb_t match, mismatch;

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

int
main(int argc, char *argc[])
{
    /* use hardcoded parameters but command-line sequences */
    slong i;

    fmpq_t birth, death, t;
    fmpq pi[4];
    fmpq qi[4];

    /* birth = 1 */
    fmpq_init(birth);
    fmpq_set_si(birth, 1, 1);

    /* death = 1/2 */
    fmpq_init(death);
    fmpq_set_si(death, 1, 2);

    /* t = 1/10 */
    fmpq_init(t);
    fmpq_set_si(t, 1, 10);

    /* pi = (0.27, 0.24, 0.26, 0.23) */
    for (i = 0; i < 4; i++)
    {
        fmpq_init(pi+i);
    }
    fmpq_set_si(pi+0, 27, 100);
    fmpq_set_si(pi+1, 24, 100);
    fmpq_set_si(pi+2, 26, 100);
    fmpq_set_si(pi+3, 23, 100);

    /* qi = 1 - pi */
    for (i = 0; i < 4; i++)
    {
        fmpq_init(qi+i);
    }

    n = 10;
    m = 10;

    /* initialize log gamma vector */
    v = 
    _get_log_gamma(arb_ptr * v, ulong n, const fmpq_t birth, const fmpq_t death,
                   slong prec);

    return 0;
}
