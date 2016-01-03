/*
 * Double precision tkf91 dynamic programming.
 */

#include <time.h>

#include "arb_mat.h"

#include "tkf91_dp.h"
#include "tkf91_dp_d.h"
#include "breadcrumbs.h"
#include "printutil.h"


typedef struct
{
    double m0;
    double m1;
    double m2;
} tnode_struct;
typedef tnode_struct tnode_t[1];
typedef tnode_struct * tnode_ptr;

typedef struct
{
    tnode_ptr data;
    slong r;
    slong c;
} tmat_struct;
typedef tmat_struct tmat_t[1];

static void tmat_init(tmat_t mat, slong nrows, slong ncols);
static void tmat_clear(tmat_t mat);
static void tmat_get_alignment(solution_t sol,
        const tmat_t mat, const slong *A, const slong *B);

static __inline__ slong
tmat_nrows(const tmat_t mat)
{
    return mat->r;
}

static __inline__ slong
tmat_ncols(const tmat_t mat)
{
    return mat->c;
}

static __inline__ tnode_ptr
tmat_entry(tmat_t mat, slong i, slong j)
{
    return mat->data + i * mat->c + j;
}

static __inline__ tnode_ptr
tmat_srcentry(const tmat_t mat, slong i, slong j)
{
    return mat->data + i * mat->c + j;
}

static __inline__ tnode_ptr
tmat_entry_top(tmat_t mat, slong i, slong j)
{
    return tmat_entry(mat, i-1, j);
}

static __inline__ tnode_ptr
tmat_entry_diag(tmat_t mat, slong i, slong j)
{
    return tmat_entry(mat, i-1, j-1);
}

static __inline__ tnode_ptr
tmat_entry_left(tmat_t mat, slong i, slong j)
{
    return tmat_entry(mat, i, j-1);
}

void
tmat_init(tmat_t mat, slong nrows, slong ncols)
{
    mat->data = flint_malloc(nrows * ncols * sizeof(tmat_struct));
    mat->r = nrows;
    mat->c = ncols;
}

void
tmat_clear(tmat_t mat)
{
    flint_free(mat->data);
}

void
tmat_get_alignment(
        solution_t sol,
        const tmat_t mat, const slong *A, const slong *B)
{
    slong i, j;
    char ACGT[4] = "ACGT";
    char tmp;
    slong len, nrows, ncols;
    double max3;
    tnode_ptr cell;
    char * sa;
    char * sb;

    sa = sol->A;
    sb = sol->B;

    nrows = tmat_nrows(mat);
    ncols = tmat_ncols(mat);
    i = nrows - 1;
    j = ncols - 1;
    len = 0;
    while (i > 0 || j > 0)
    {
        cell = tmat_srcentry(mat, i, j);
        max3 = fmax(cell->m0, fmax(cell->m1, cell->m2));
        if (cell->m0 == max3)
        {
            sa[len] = ACGT[A[i-1]];
            sb[len] = '-';
            i--;
        }
        else if (cell->m1 == max3)
        {
            sa[len] = ACGT[A[i-1]];
            sb[len] = ACGT[B[j-1]];
            i--;
            j--;
        }
        else if (cell->m2 == max3)
        {
            sa[len] = '-';
            sb[len] = ACGT[B[j-1]];
            j--;
        }
        else
        {
            flint_printf("lost the thread ");
            flint_printf("in the dynamic programing traceback\n");
            abort();
        }
        len++;
    }
    for (i = 0; i < len/2; i++)
    {
        j = len - 1 - i;
        tmp = sa[i]; sa[i] = sa[j]; sa[j] = tmp;
        tmp = sb[i]; sb[i] = sb[j]; sb[j] = tmp;
    }

    sol->len = len;
}

static __inline__ double
_arb_get_d(const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

/* helper function for converting the generator array to double precision */
/* m should be a column vector */
static __inline__ double
_doublify(slong i, const arb_mat_t m)
{
    return _arb_get_d(arb_mat_entry(m, i, 0));
}


void
tkf91_dynamic_programming_double_tmat(
        solution_t sol, const request_t req,
        const tkf91_generator_indices_t g,
        const arb_mat_t m,
        const slong *A, slong szA,
        const slong *B, slong szB);

void
tkf91_dynamic_programming_double_tmat(
        solution_t sol, const request_t req,
        const tkf91_generator_indices_t g,
        const arb_mat_t m,
        const slong *A, slong szA,
        const slong *B, slong szB)
{
    slong nrows, ncols;
    tmat_t tmat;
    double p0_max3, p1_max3, p2_max2;
    slong i, j;
    tnode_ptr cell, p0, p1, p2;
    slong nta, ntb;
    clock_t start;

    /* dynamic programming 'generators' as local variables */
    double m1_00;
    double m0_10;
    double m0_i0_incr[4];
    double m2_01;
    double m2_0j_incr[4];
    double c0_incr[4];
    double c1_incr[16];
    double c2_incr[4];

    /* values that are cached per row in the main loop */
    double c0_incr_nta;
    double * c1_incr_nta;

    /* start the clock */
    start = clock();

    /* init the dynamic programming 'generators' */
    m1_00 = _doublify(g->m1_00, m);
    m0_10 = _doublify(g->m0_10, m);
    m2_01 = _doublify(g->m2_01, m);
    for (i = 0; i < 4; i++)
    {
        m0_i0_incr[i] = _doublify(g->m0_i0_incr[i], m);
        m2_0j_incr[i] = _doublify(g->m2_0j_incr[i], m);
        c0_incr[i] = _doublify(g->c0_incr[i], m);
        for (j = 0; j < 4; j++)
        {
            c1_incr[i*4+j] = _doublify(g->c1_incr[i*4+j], m);
        }
        c2_incr[i] = _doublify(g->c2_incr[i], m);
    }

    nrows = szA + 1;
    ncols = szB + 1;

    tmat_init(tmat, nrows, ncols);

    /* corner */
    i = 0;
    j = 0;
    cell = tmat_entry(tmat, i, j);
    cell->m0 = -INFINITY;
    cell->m2 = -INFINITY;
    cell->m1 = m1_00;

    /* top edge */
    i = 0;
    for (j = 1; j < ncols; j++)
    {
        ntb = B[j - 1];
        cell = tmat_entry(tmat, i, j);
        cell->m0 = -INFINITY;
        cell->m1 = -INFINITY;
        if (j == 1)
        {
            cell->m2 = m2_01;
        }
        else
        {
            p2 = tmat_entry_left(tmat, i, j);
            p2_max2 = fmax(p2->m1, p2->m2);
            cell->m2 = p2_max2 + m2_0j_incr[ntb];
        }
    }

    /* left edge */
    j = 0;
    for (i = 1; i < nrows; i++)
    {
        nta = A[i - 1];
        cell = tmat_entry(tmat, i, j);
        cell->m1 = -INFINITY;
        cell->m2 = -INFINITY;
        if (i == 1)
        {
            cell->m0 = m0_10;
        }
        else
        {
            p0 = tmat_entry_top(tmat, i, j);
            p0_max3 = fmax(p0->m0, fmax(p0->m1, p0->m2));
            cell->m0 = p0_max3 + m0_i0_incr[nta];
        }
    }

    /* main loop */
    for (i = 1; i < nrows; i++)
    {
        nta = A[i - 1];

        /* precompute stuff for this row */
        c0_incr_nta = c0_incr[nta];
        c1_incr_nta = c1_incr + 4*nta;

        tnode_ptr prev_row = tmat->data + (i-1)*ncols;
        tnode_ptr curr_row = tmat->data + i*ncols;
        for (j = 1; j < ncols; j++)
        {
            ntb = B[j - 1];

            cell = curr_row + j;
            p0 = prev_row + j;
            p1 = prev_row + j-1;
            p2 = curr_row + j-1;

            p0_max3 = fmax(p0->m0, fmax(p0->m1, p0->m2));
            p1_max3 = fmax(p1->m0, fmax(p1->m1, p1->m2));
            p2_max2 = fmax(p2->m1, p2->m2);

            cell->m0 = p0_max3 + c0_incr_nta;
            cell->m1 = p1_max3 + c1_incr_nta[ntb];
            cell->m2 = p2_max2 + c2_incr[ntb];
        }
    }

    /* compute the log probability of the optimal alignment */
    double logp;
    cell = tmat_entry(tmat, nrows-1, ncols-1);
    logp = fmax(cell->m0, fmax(cell->m1, cell->m2));
    arb_set_d(sol->log_probability, logp);

    _fprint_elapsed(stderr, "forward dynamic programming", clock() - start);


    /* do the traceback if requested */
    if (req->trace)
    {
        start = clock();
        tmat_get_alignment(sol, tmat, A, B);
        _fprint_elapsed(stderr, "traceback", clock() - start);
    }

    start = clock();
    tmat_clear(tmat);
    _fprint_elapsed(stderr, "cleanup", clock() - start);
}

void
tkf91_dp_d(
        solution_t sol, const request_t req,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, size_t szA,
        const slong *B, size_t szB)
{
    slong level = 8;
    slong prec = 1 << level;

    arb_t x;
    arb_mat_t G;
    arb_mat_t expression_logs;
    arb_mat_t generator_logs;
    slong i;
    slong generator_count = fmpz_mat_nrows(mat);
    slong expression_count = fmpz_mat_ncols(mat);

    arb_init(x);

    /* initialize the arbitrary precision exponent matrix */
    arb_mat_init(G, generator_count, expression_count);
    arb_mat_set_fmpz_mat(G, mat);

    /* compute the expression logarithms */
    arb_mat_init(expression_logs, expression_count, 1);
    for (i = 0; i < expression_count; i++)
    {
        expr_eval(x, expressions_table[i], level);
        arb_log(arb_mat_entry(expression_logs, i, 0), x, prec);
    }

    /* compute the generator logarithms */
    arb_mat_init(generator_logs, generator_count, 1);
    arb_mat_mul(generator_logs, G, expression_logs, prec);

    tkf91_dynamic_programming_double_tmat(
            sol, req, g, generator_logs, A, szA, B, szB);

    arb_clear(x);
    arb_mat_clear(G);
    arb_mat_clear(expression_logs);
    arb_mat_clear(generator_logs);
}
