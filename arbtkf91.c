


typedef struct
{
    slong *data;
    ulong m;
    ulong n;
    ulong k;
    ulong l;
}
tableau_struct;

typedef tableau_struct tableau_t[1];

slong *
tableau_entry(tableau_t T, ulong x, ulong y, ulong z, ulong w)
{
    return T->data + T->l * (T->k * (T->n * x + y) + z) + w;
}

slong *
tableau_vec(tableau_t T, ulong x, ulong y, ulong z)
{
    return T->data + T->l * (T->k * (T->n * x + y) + z);
}

slong *tableau_M0_vec(tableau_t T, ulong x, ulong y) {
    return tableau_vec(T, x, y, 0);
}
slong *tableau_M1_vec(tableau_t T, ulong x, ulong y) {
    return tableau_vec(T, x, y, 1);
}
slong *tableau_M2_vec(tableau_t T, ulong x, ulong y) {
    return tableau_vec(T, x, y, 2);
}

void tableau_M0_add(tableau_t T, ulong x, ulong y, ulong w, slong value) {
    *tableau_entry(T, x, y, 0, w) += value;
}
void tableau_M1_add(tableau_t T, ulong x, ulong y, ulong w, slong value) {
    *tableau_entry(T, x, y, 1, w) += value;
}
void tableau_M2_add(tableau_t T, ulong x, ulong y, ulong w, slong value) {
    *tableau_entry(T, x, y, 2, w) += value;
}

void
tableau_init(tableau_t T, ulong m, ulong n, ulong k, ulong l,
        const int *A, const int *B)
{
    ulong x, y, z, w;
    ulong i, j;

    /* brutal error checking */
    if (k != 3) abort();
    if (l != 25) abort();

    /* allocate the table */
    T->data = (slong *) calloc(m*n*k*l, sizeof(long));
    T->m = m;
    T->n = n;
    T->k = k;
    T->l = l;
    
    /*
     * Initialize the entries of the tableau.
     * These are integer coefficients defining a linear combination
     * of complicated expressions.
     * Only the coefficients matter here, not the expression values.
     */

    /* M1(0,0) = log(1 - lam/mu) + log(1 - lam*b) */
    tableau_M1_add(T, 0, 0, IDX_LAMBDA_MU, 1);
    tableau_M1_add(T, 0, 0, IDX_LAMBDA_BETA, 1);

    /* M0(1,0) = log(1 - lam/mu) + log(1 - lam*b) + log(b) + log(pi_A_1) */
    tableau_M0_add(T, 1, 0, IDX_LAMBDA_MU, 1);
    tableau_M0_add(T, 1, 0, IDX_LAMBDA_BETA, 1);
    tableau_M0_add(T, 1, 0, IDX_BETA, 1);
    tableau_M0_add(T, 1, 0, IDX_PI + A[1], 1);

    /* M0(i,0) = M0(i-1,0) + 2*log(lam*b) + log(pi_A_i), i >= 2 */
    for (i = 2; i < m; i++)
    {
        _slong_vec_add_inplace(
                tableau_M0_vec(T, i, 0),
                tableau_M0_vec(T, i-1, 0), k);
        tableau_M0_add(T, i, 0, IDX_LAMDA, 2);
        tableau_M0_add(T, i, 0, IDX_BETA, 2);
        tableau_M0_add(T, i, 0, IDX_PI + A[i], 1);
    }

    /* M2(0,1) = log(1 - lam/mu) + log(1 - lam*b) + log(lam*b) + log(pi_B_1) */
    tableau_M2_add(T, 0, 1, IDX_LAMBDA_MU, 1);
    tableau_M2_add(T, 0, 1, IDX_LAMBDA_BETA, 1);
    tableau_M2_add(T, 0, 1, IDX_LAMBDA, 1);
    tableau_M2_add(T, 0, 1, IDX_BETA, 1);
    tableau_M2_add(T, 0, 1, IDX_PI + B[1], 1);

    /* M2(0,j) = M2(0,j-1) + log(lam*b) + log(pi_B_j), j >= 2 */
    for (j = 2; j < n; j++)
    {
        _slong_vec_add_inplace(
                tableau_M2_vec(T, 0, j),
                tableau_M2_vec(T, 0, j-1), k);
        tableau_M2_add(T, 0, j, IDX_LAMBDA, 1);
        tableau_M2_add(T, 0, j, IDX_BETA, 1);
        tableau_M2_add(T, 0, j, IDX_PI + B[j], 1);
    }

    /*
     * Initialize the rest of the tableau.
     * I think this corresponds to all entries such that min(i, j) >= 1.
     */
    for (i = 1; i < n+1; i++)
    {
        for (j = 1; j < m+1; j++)
        {
            /* M0(i,j) = log(lam*b*pi_A_i) */
            tableau_M0_add(T, i, j, IDX_LAMBDA, 1);
            tableau_M0_add(T, i, j, IDX_BETA, 1);
            tableau_M0_add(T, i, j, IDX_PI + A[i], 1);

            /* M1(i, j) += log(lam/mu * pi_A_i) */
            tableau_M1_add(T, i, j, IDX_LAMBDA, 1);
            tableau_M1_add(T, i, j, IDX_MU, -1);
            tableau_M1_add(T, i, j, IDX_PI + A[i], 1);

            /* M1(i, j) += log(1 - lam*b) + trans_Ai_Bj */
            tableau_M1_add(T, i, j, IDX_LAM_BETA, 1);
            tableau_M1_add(T, i, j, IDX_TRANS + A[i]*4 + B[j], 1);

            /* M2(i, j) += log(lam*b*pi_B_j) */
            tableau_M2_add(T, i, j, IDX_LAMBDA, 1);
            tableau_M2_add(T, i, j, IDX_BETA, 1);
            tableau_M2_add(T, i, j, IDX_PI + B[j], 1);
        }
    }

    for (x = 0; x < m; x++)
    {
        for (y = 0; y < n; y++)
        {
            tableau_entry(T, x, y, 0);
        }
    }
}

void
tableau_clear(tableau_t T)
{
    free(T->data);
    T->data = NULL;
    T->m = 0;
    T->n = 0;
    T->k = 0;
    T->l = 0;
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
        fmpq_submul(dt, pi+i, pi+i);
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


void
_fill_sequence_vector(slong *v, const char *str, size_t n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        switch(str[i])
        {
            case 'A' : v[i] = 0; break;
            case 'C' : v[i] = 1; break;
            case 'G' : v[i] = 2; break;
            case 'T' : v[i] = 3; break;
            default:
               abort();
        }
    }
}

int
main(int argc, char *argc[])
{
    user_params_t p;
    tableau_t T;

    const char strA[] = "ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA";
    const char strB[] = "AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA";

    slong *A;
    slong *B;

    size_t szA, szB;

    ulong m, n, k, l;

    slong i;


    /* initialize the parameter vector with hardcoded parameter values */

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


    /* initialize the tableau with hardcoded sequences */

    szA = strlen(strA);
    A = (slong *) malloc(szA * sizeof(slong+1));
    A[0] = -1;
    _fill_sequence_vector(A+1, strA, szA);

    szB = strlen(strB);
    B = (slong *) malloc(szA * sizeof(slong+1));
    B[0] = -1;
    _fill_sequence_vector(B+1, strB, szB);

    /* define the dimensions of the tableau */
    m = szA + 1;
    n = szB + 1;
    k = 3;
    l = 25;

    /* allocate the tableau */
    tableau_init(T, m, n, k, l);


    tableau_clear(T);
    user_params_clear(p);

    free(A);
    free(B);

    return 0;
}


int demo()
{
    fmpq_t a, b, t;
    fmpq_t at, bt;

    fmpq_init(a);
    fmpq_init(b);
    fmpq_init(t);

    fmpq_set_si(a, 1, 2);
    fmpq_set_si(b, 3, 4);
    fmpq_set_si(t, 5, 6);

    fmpq_init(at);
    fmpq_init(bt);

    fmpq_mul(at, a, t);
    fmpq_mul(bt, b, t);

    expr_t ax, bx;
    expr_t eat, ebt, aeat, bebt;
    expr_t num, den;
    expr_t res;

    expr_fmpq(ax, a);
    expr_fmpq(bx, b);
    expr_exp_fmpq(ebt, bt);
    expr_exp_fmpq(eat, at);
    expr_mul(bebt, bx, ebt);
    expr_mul(aeat, ax, eat);
    expr_sub(num, ebt, eat);
    expr_sub(den, bebt, aeat);
    expr_div(res, num, den);

    arb_t value;
    slong level;

    arb_init(value);

    for (level = 0; level < 8; level++)
    {
        flint_printf("evaluating ");
        expr_print(res);
        flint_printf(" at level %wd:\n", level);
        expr_eval(value, res, level);
        arb_print(value);
        flint_printf("\n\n");
    }

    arb_clear(value);
    fmpq_clear(a);
    fmpq_clear(b);
    fmpq_clear(t);
    expr_clear(ax);
    expr_clear(bx);
    expr_clear(eat);
    expr_clear(ebt);
    expr_clear(aeat);
    expr_clear(bebt);
    expr_clear(num);
    expr_clear(den);
    expr_clear(res);

    return 0;
}

int moredemo()
{
    /* unfinished */

    slong i, j;
    slong level;

    fmpq_t a, b, t;
    fmpq_t dt, negdt;
    fmpq pi[4];

    fmpq_init(a);
    fmpq_init(b);
    fmpq_init(t);
    fmpq_init(dt);
    fmpq_init(negdt);
    for (i = 0; i < 4; i++)
    {
        fmpq_init(pi+i);
    }

    fmpq_set_si(a, 1, 1);
    fmpq_set_si(b, 2, 1);
    fmpq_set_si(t, 1, 10);

    fmpq_set_si(pi+0, 27, 100);
    fmpq_set_si(pi+1, 24, 100);
    fmpq_set_si(pi+2, 26, 100);
    fmpq_set_si(pi+3, 23, 100);

    _get_dt(dt, pi, t);
    fmpq_neg(negdt, dt);

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            if (i == j)
            {
                ;
            }
            else
            {
                ;
            }
        }
    }

    fmpq_clear(a);
    fmpq_clear(b);
    fmpq_clear(t);
    fmpq_clear(dt);
    fmpq_clear(negdt);
    for (i = 0; i < 4; i++)
    {
        fmpq_clear(pi+i);
    }
}
