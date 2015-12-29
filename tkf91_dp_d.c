/*
 * Double precision tkf91 dynamic programming.
 */

#include <time.h>

#include "tkf91_dp_d.h"
#include "breadcrumbs.h"
#include "wavefront_double.h"
#include "printutil.h"

static __inline__
double max(double a, double b)
{
    return a > b ? a : b;
}


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
static void tmat_get_alignment(char **psa, char **psb,
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
tmat_get_alignment(char **psa, char **psb,
        const tmat_t mat, const slong *A, const slong *B)
{
    slong i, j;
    char ACGT[4] = "ACGT";
    char * sa;
    char * sb;
    char tmp;
    slong len, nrows, ncols;
    double max3;
    tnode_ptr cell;

    nrows = tmat_nrows(mat);
    ncols = tmat_ncols(mat);
    sa = calloc(nrows * ncols, sizeof(char));
    sb = calloc(nrows * ncols, sizeof(char));
    i = nrows - 1;
    j = ncols - 1;
    len = 0;
    while (i > 0 || j > 0)
    {
        cell = tmat_srcentry(mat, i, j);
        /* flint_printf("%lf %lf %lf\n", cell->m0, cell->m1, cell->m2); */
        max3 = max(cell->m0, max(cell->m1, cell->m2));
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
    *psa = sa;
    *psb = sb;
}






typedef struct
{
    double max2;
    double max3;
    breadcrumb_t crumb;
} dnode_struct;
typedef dnode_struct dnode_t[1];
typedef dnode_struct * dnode_ptr;

typedef struct
{
    dnode_ptr data;
    slong r;
    slong c;
} dmat_struct;
typedef dmat_struct dmat_t[1];

void dmat_init(dmat_t mat, slong nrows, slong ncols);
void dmat_clear(dmat_t mat);
void dmat_get_alignment(char **psa, char **psb,
        const dmat_t mat, const slong *A, const slong *B);

static __inline__ slong
dmat_nrows(const dmat_t mat)
{
    return mat->r;
}

static __inline__ slong
dmat_ncols(const dmat_t mat)
{
    return mat->c;
}

static __inline__ dnode_ptr
dmat_entry(dmat_t mat, slong i, slong j)
{
    /* return mat->data + (i % 2) * mat->c + j; */
    return mat->data + i * mat->c + j;
}

static __inline__ dnode_ptr
dmat_srcentry(const dmat_t mat, slong i, slong j)
{
    /* return mat->data + (i % 2) * mat->c + j; */
    return mat->data + i * mat->c + j;
}

static __inline__ dnode_ptr
dmat_entry_top(dmat_t mat, slong i, slong j)
{
    return dmat_entry(mat, i-1, j);
}

static __inline__ dnode_ptr
dmat_entry_diag(dmat_t mat, slong i, slong j)
{
    return dmat_entry(mat, i-1, j-1);
}

static __inline__ dnode_ptr
dmat_entry_left(dmat_t mat, slong i, slong j)
{
    return dmat_entry(mat, i, j-1);
}

void
dmat_init(dmat_t mat, slong nrows, slong ncols)
{
    /* mat->data = flint_malloc(2 * ncols * sizeof(dmat_struct)); */
    mat->data = flint_malloc(nrows * ncols * sizeof(dmat_struct));
    mat->r = nrows;
    mat->c = ncols;
}

void
dmat_clear(dmat_t mat)
{
    flint_free(mat->data);
}

void
dmat_get_alignment(char **psa, char **psb,
        const dmat_t mat, const slong *A, const slong *B)
{
    slong i, j;
    char ACGT[4] = "ACGT";
    char * sa;
    char * sb;
    char tmp;
    slong len, nrows, ncols, n;
    breadcrumb_t crumb;

    nrows = dmat_nrows(mat);
    ncols = dmat_ncols(mat);
    n = nrows * ncols;
    sa = calloc(n, sizeof(char));
    sb = calloc(n, sizeof(char));
    i = nrows - 1;
    j = ncols - 1;
    len = 0;
    while (i > 0 || j > 0)
    {
        crumb = dmat_srcentry(mat, i, j)->crumb;
        if (crumb & CRUMB_TOP)
        {
            sa[len] = ACGT[A[i-1]];
            sb[len] = '-';
            i--;
        }
        else if (crumb & CRUMB_DIAG)
        {
            sa[len] = ACGT[A[i-1]];
            sb[len] = ACGT[B[j-1]];
            i--;
            j--;
        }
        else if (crumb & CRUMB_LEFT)
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
    *psa = sa;
    *psb = sb;
}



typedef struct
{
    double m1_00;
    double m0_10;
    double m0_i0_incr[4];
    double m2_01;
    double m2_0j_incr[4];
    double c0_incr[4];
    double c1_incr[16];
    double c2_incr[4];
} tkf91_double_generators_struct;
typedef tkf91_double_generators_struct tkf91_double_generators_t[1];


void doublify_tkf91_generators(
        tkf91_double_generators_t h,
        tkf91_generator_indices_t g,
        arb_mat_t m);

void tkf91_dynamic_programming_double(
        const tkf91_generator_indices_t g,
        const arb_mat_t m,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB);


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


/* m should be a column vector of per-generator values */
void
doublify_tkf91_generators(
        tkf91_double_generators_t h,
        tkf91_generator_indices_t g,
        arb_mat_t m)
{
    slong i, j;
    h->m1_00 = _doublify(g->m1_00, m);
    h->m0_10 = _doublify(g->m0_10, m);
    h->m2_01 = _doublify(g->m2_01, m);
    for (i = 0; i < 4; i++)
    {
        h->m0_i0_incr[i] = _doublify(g->m0_i0_incr[i], m);
        h->m2_0j_incr[i] = _doublify(g->m2_0j_incr[i], m);
        h->c0_incr[i] = _doublify(g->c0_incr[i], m);
        for (j = 0; j < 4; j++)
        {
            h->c1_incr[i*4+j] = _doublify(g->c1_incr[i*4+j], m);
        }
        h->c2_incr[i] = _doublify(g->c2_incr[i], m);
    }
}

void
tkf91_dynamic_programming_double_tmat(
        const tkf91_generator_indices_t g,
        const arb_mat_t m,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB);

void
tkf91_dynamic_programming_double_tmat(
        const tkf91_generator_indices_t g,
        const arb_mat_t m,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB)
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
            p2_max2 = max(p2->m1, p2->m2);
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
            p0_max3 = max(p0->m0, max(p0->m1, p0->m2));
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

            /*
            cell = tmat_entry(tmat, i, j);
            p0 = tmat_entry_top(tmat, i, j);
            p1 = tmat_entry_diag(tmat, i, j);
            p2 = tmat_entry_left(tmat, i, j);
            */

            p0_max3 = max(p0->m0, max(p0->m1, p0->m2));
            p1_max3 = max(p1->m0, max(p1->m1, p1->m2));
            p2_max2 = max(p2->m1, p2->m2);

            cell->m0 = p0_max3 + c0_incr_nta;
            cell->m1 = p1_max3 + c1_incr_nta[ntb];
            cell->m2 = p2_max2 + c2_incr[ntb];

            /*
            cell->m0 = p0_max3 + c0_incr[nta];
            cell->m1 = p1_max3 + c1_incr[nta*4 + ntb];
            cell->m2 = p2_max2 + c2_incr[ntb];
            */
        }
    }

    /* report the score */
    double max3;
    cell = tmat_entry(tmat, nrows-1, ncols-1);
    max3 = max(cell->m0, max(cell->m1, cell->m2));
    flint_printf("score: %g\n", exp(max3));

    /* stop the clock */
    printf("pre-traceback dynamic programming ");
    _print_elapsed_time(clock() - start);


    /* do the traceback */
    if (trace_flag)
    {

        /* restart the clock for traceback */
        start = clock();

        char *sa, *sb;
        /* breadcrumb_mat_get_alignment(&sa, &sb, crumb_mat, A, B); */
        tmat_get_alignment(&sa, &sb, tmat, A, B);
        flint_printf("%s\n", sa);
        flint_printf("%s\n", sb);
        flint_printf("\n");
        free(sa);
        free(sb);

        printf("traceback ");
        _print_elapsed_time(clock() - start);
    }



    /* restart the clock for cleanup */
    start = clock();

    tmat_clear(tmat);

    printf("cleanup ");
    _print_elapsed_time(clock() - start);
}



void tkf91_dynamic_programming_double(
        const tkf91_generator_indices_t g,
        const arb_mat_t m,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB)
{
    slong nrows, ncols;
    breadcrumb_mat_t crumb_mat;
    /*
    breadcrumb_ptr pcrumb;
    breadcrumb_t crumb;
    */
    dmat_t dmat;
    double m0, m1, m2, max2, max3;
    slong i, j;
    dnode_ptr cell, p0, p1, p2;
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

    /* values that are cached per row */
    double m0_i0_incr_nta;
    double c0_incr_nta;
    double c1_incr_nta[4];

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

    /* redundantly init precomputed values to avoid compiler warning */
    nta = 0;
    m0_i0_incr_nta = m0_i0_incr[nta];
    c0_incr_nta = c0_incr[nta];
    for (ntb = 0; ntb < 4; ntb++)
    {
        c1_incr_nta[ntb] = c1_incr[4*nta + ntb];
    }

    nrows = szA + 1;
    ncols = szB + 1;

    if (trace_flag)
    {
        breadcrumb_mat_init(crumb_mat, nrows, ncols);
    }

    dmat_init(dmat, nrows, ncols);

    for (i = 0; i < nrows; i++)
    {
        /* precompute values that are constant along the row */
        if (i)
        {
            nta = A[i - 1];
            m0_i0_incr_nta = m0_i0_incr[nta];
            c0_incr_nta = c0_incr[nta];
            for (ntb = 0; ntb < 4; ntb++)
            {
                c1_incr_nta[ntb] = c1_incr[4*nta + ntb];
            }
        }

        for (j = 0; j < ncols; j++)
        {
            ntb = B[j - 1];

            if (i < 1 || j < 1)
            {
                m0 = -INFINITY;
                m1 = -INFINITY;
                m2 = -INFINITY;
                if (i == 0 && j == 0)
                {
                    m1 = m1_00;
                }
                else if (i == 1 && j == 0)
                {
                    m0 = m0_10;
                }
                else if (i == 0 && j == 1)
                {
                    m2 = m2_01;
                }
                else
                {
                    if (i == 0)
                    {
                        p2 = dmat_entry_left(dmat, i, j);
                        m2 = p2->max2 + m2_0j_incr[ntb];
                    }
                    else if (j == 0)
                    {
                        p0 = dmat_entry_top(dmat, i, j);
                        m0 = p0->max3 + m0_i0_incr_nta;
                    }
                }
            }
            else
            {
                p0 = dmat_entry_top(dmat, i, j);
                p1 = dmat_entry_diag(dmat, i, j);
                p2 = dmat_entry_left(dmat, i, j);
                m0 = p0->max3 + c0_incr_nta;
                m1 = p1->max3 + c1_incr_nta[ntb];
                m2 = p2->max2 + c2_incr[ntb];
            }
            max2 = max(m1, m2);
            max3 = max(m0, max2);

            cell = dmat_entry(dmat, i, j);
            cell->max2 = max2;
            cell->max3 = max3;

            /* fill the table for traceback */
            if (trace_flag)
            {
                cell->crumb = 0;
                if (m0 == max3) {
                    cell->crumb |= CRUMB_TOP;
                }
                else if (m1 == max3) {
                    cell->crumb |= CRUMB_DIAG;
                }
                else if (m2 == max3) {
                    cell->crumb |= CRUMB_LEFT;
                }
            }
        }
    }

    /* report the score */
    cell = dmat_entry(dmat, nrows-1, ncols-1);
    flint_printf("score: %g\n", exp(cell->max3));

    /* stop the clock */
    printf("pre-traceback dynamic programming ");
    _print_elapsed_time(clock() - start);


    /* do the traceback */
    if (trace_flag)
    {

        /* restart the clock for traceback */
        start = clock();

        char *sa, *sb;
        /* breadcrumb_mat_get_alignment(&sa, &sb, crumb_mat, A, B); */
        dmat_get_alignment(&sa, &sb, dmat, A, B);
        flint_printf("%s\n", sa);
        flint_printf("%s\n", sb);
        flint_printf("\n");
        free(sa);
        free(sb);

        printf("traceback ");
        _print_elapsed_time(clock() - start);
    }



    /* restart the clock for cleanup */
    start = clock();

    dmat_clear(dmat);
    if (trace_flag)
    {
        breadcrumb_mat_clear(crumb_mat);
    }

    printf("cleanup ");
    _print_elapsed_time(clock() - start);
}




void
tkf91_dp_d(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g,
        int trace_flag,
        slong *A, size_t szA,
        slong *B, size_t szB)
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

    /* tkf91_dynamic_programming_double( */
    tkf91_dynamic_programming_double_tmat(
            g, generator_logs, trace_flag, A, szA, B, szB);

    arb_clear(x);
    arb_mat_clear(G);
    arb_mat_clear(expression_logs);
    arb_mat_clear(generator_logs);
}
