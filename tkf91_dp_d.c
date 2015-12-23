/*
 * Double precision tkf91 dynamic programming.
 */

#include <time.h>

#include "tkf91_dp_d.h"
#include "breadcrumbs.h"
#include "wavefront_double.h"


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


double _doublify(slong i, arb_mat_t m);

void doublify_tkf91_generators(
        tkf91_double_generators_t h,
        tkf91_generator_indices_t g,
        arb_mat_t m);

void tkf91_dynamic_programming_double(tkf91_double_generators_t g,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB);


static __inline__
double max2(double a, double b)
{
    return a > b ? a : b;
}

static __inline__
double max3(double a, double b, double c)
{
    return max2(a, max2(b, c));
}


/* helper function for converting the generator array to double precision */
/* m should be a column vector */
double
_doublify(slong i, arb_mat_t m)
{
    arb_ptr p = arb_mat_entry(m, i, 0);
    return arf_get_d(arb_midref(p), ARF_RND_NEAR);
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



void tkf91_dynamic_programming_double(tkf91_double_generators_t g,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB)
{
    /*
     * Make the dynamic programming table.
     *
     * Inputs:
     *  - a structure whose member variables store generator indices
     *  - a list that maps expression indices to expression objects
     *  - a matrix G such that G_{ij} indicates the integer exponent of
     *    expression j in generator i
     *  - the two sequences to be aligned, as integer arrays
     *  - the two sequence lengths
     *
     * The dynamic programming can use a trick that
     * packs three diagonals into a table with two rows,
     * using an idea like the following diagram.
     * Notice that if you are iterating through rows of that diagram
     * using an algorithm that needs to track the most recent three rows,
     * the sparsity of the entries in that diagram allow you to use
     * a physical buffer of only two rows.
     *
     *    |    |  0 |  1 |  2 |  3 |  4 |  5
     *  --|----|----|----|----|----|----|---
     *  0 |    |    |    | 00 |    |    |          
     *  --|----|----|----|----|----|----|---
     *  1 |    |    | 10 |    | 01 |    |   
     *  --|----|----|----|----|----|----|---
     *  2 |    | 20 |    | 11 |    | 02 |   
     *  --|----|----|----|----|----|----|---
     *  3 |    |    | 21 |    | 12 |    | 03
     *  --|----|----|----|----|----|----|---
     *  4 |    |    |    | 22 |    | 13 |   
     *  --|----|----|----|----|----|----|---
     *  5 |    |    |    |    | 23 |    |   
     *
     */

    /*
     * Define the matrix to be used for the traceback.
     * The number of rows is one greater than the length of the first sequence,
     * and the number of columns is one greater than the length of the second
     * sequence.
     */
    slong nrows = szA + 1;
    slong ncols = szB + 1;
    breadcrumb_mat_t crumb_mat;

    if (trace_flag)
    {
        breadcrumb_mat_init(crumb_mat, nrows, ncols);
    }

    /* define the wavefront matrix */
    wave_mat_t wave;
    /* slong modulus = nrows + ncols - 1; */
    /* slong modulus = 3; */
    /* wave_mat_init(wave, nrows + ncols - 1, modulus); */
    wave_mat_init(wave, nrows + ncols - 1);

    /*
     * Let M_{ij} be the matrix created for traceback.
     * Let R_{kl} be the 'logical' wavefront matrix.
     * Then
     *  k = i + j
     *  l = nrows - 1 + j - i
     */


    /* iterate over anti-diagonal bands of the dynamic programming matrix */
    wave_value_ptr cell, p0, p1, p2;
    slong i, j, k, l;
    slong istart, jstart, lstart;
    for (k = 0; k < nrows + ncols - 1; k++)
    {
        if (k < nrows)
        {
            istart = k;
            jstart = 0;
        }
        else
        {
            istart = nrows - 1;
            jstart = k - (nrows - 1);
        }
        lstart = nrows - 1 + jstart - istart;

        /* iterate over entries of the diagonal */
        i = istart;
        j = jstart;
        l = lstart;
        slong nta, ntb;
        while (0 <= i && i < nrows && 0 <= j && j < ncols)
        {
            /* check some invariants */
            /*
            if (k != i + j)
            {
                flint_printf("wavefront indexing problem ");
                flint_printf("i=%wd j=%wd k=%wd\n", i, j, k);
                abort();
            }
            if (l != nrows - 1 + j - i)
            {
                flint_printf("wavefront indexing problem ");
                flint_printf("i=%wd j=%wd l=%wd\n", i, j, l);
                abort();
            }
            */

            cell = wave_mat_entry(wave, k, l);
            if (i < 1 || j < 1)
            {
                if (i == 0 && j == 0)
                {
                    cell->m0 = -INFINITY;
                    cell->m1 = g->m1_00;
                    cell->m2 = -INFINITY;
                }
                else if (i == 1 && j == 0)
                {
                    cell->m0 = g->m0_10;
                    cell->m1 = -INFINITY;
                    cell->m2 = -INFINITY;
                }
                else if (i == 0 && j == 1)
                {
                    cell->m0 = -INFINITY;
                    cell->m1 = -INFINITY;
                    cell->m2 = g->m2_01;
                }
                else
                {
                    if (i == 0)
                    {
                        ntb = B[j - 1];
                        p2 = wave_mat_entry_left(wave, k, l);
                        cell->m0 = -INFINITY;
                        cell->m1 = -INFINITY;
                        cell->m2 = p2->m2 + g->m2_0j_incr[ntb];
                    }
                    else if (j == 0)
                    {
                        nta = A[i - 1];
                        p0 = wave_mat_entry_top(wave, k, l);
                        cell->m0 = p0->m0 + g->m0_i0_incr[nta];
                        cell->m1 = -INFINITY;
                        cell->m2 = -INFINITY;
                    }
                }
            }
            else
            {
                nta = A[i - 1];
                ntb = B[j - 1];
                p0 = wave_mat_entry_top(wave, k, l);
                p1 = wave_mat_entry_diag(wave, k, l);
                p2 = wave_mat_entry_left(wave, k, l);
                cell->m0 = max3(p0->m0, p0->m1, p0->m2) + g->c0_incr[nta];
                cell->m1 = max3(p1->m0, p1->m1, p1->m2);
                cell->m1 += g->c1_incr[nta*4 + ntb];
                cell->m2 = max2(p2->m1, p2->m2) + g->c2_incr[ntb];
            }

            /* fill the table for traceback */
            if (trace_flag)
            {
                double best;
                best = max3(cell->m0, cell->m1, cell->m2);
                breadcrumb_ptr pcrumb = breadcrumb_mat_entry(crumb_mat, i, j);
                if (cell->m0 == best) {
                    *pcrumb |= CRUMB_TOP;
                }
                if (cell->m1 == best) {
                    *pcrumb |= CRUMB_DIAG;
                }
                if (cell->m2 == best) {
                    *pcrumb |= CRUMB_LEFT;
                }
            }

            /*
            flint_printf("%wd %wd %wd %wd ", i, j, k, l);
            flint_printf("%g %g %g\n",
                    exp(cell->m0),
                    exp(cell->m1),
                    exp(cell->m2));
            */

            l += 2;
            i--;
            j++;
        }
    }

    /* report the probability matrices */

    /*
    flint_printf("m0:\n");
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            k = i + j;
            l = nrows - 1 + j - i;
            cell = wave_mat_entry(wave, k, l);
            flint_printf("%.15lf ", exp(cell->m0));
        }
        flint_printf("\n");
    }
    flint_printf("\n");

    flint_printf("m1:\n");
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            k = i + j;
            l = nrows - 1 + j - i;
            cell = wave_mat_entry(wave, k, l);
            flint_printf("%.15lf ", exp(cell->m1));
        }
        flint_printf("\n");
    }
    flint_printf("\n");

    flint_printf("m2:\n");
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            k = i + j;
            l = nrows - 1 + j - i;
            cell = wave_mat_entry(wave, k, l);
            flint_printf("%.15lf ", exp(cell->m2));
        }
        flint_printf("\n");
    }
    flint_printf("\n");
    */

    /* report the score */
    i = nrows - 1;
    j = ncols - 1;
    k = i + j;
    l = nrows - 1 + j - i;
    cell = wave_mat_entry(wave, k, l);
    flint_printf("score: %g\n", exp(max3(cell->m0, cell->m1, cell->m2)));

    /* do the traceback */
    if (trace_flag)
    {
        char *sa, *sb;
        breadcrumb_mat_get_alignment(&sa, &sb, crumb_mat, A, B);
        flint_printf("%s\n", sa);
        flint_printf("%s\n", sb);
        flint_printf("\n");
        free(sa);
        free(sb);
    }

    wave_mat_clear(wave);
    if (trace_flag)
    {
        breadcrumb_mat_clear(crumb_mat);
    }
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
    tkf91_double_generators_t h;

    arb_init(x);

    /* initialize the arbitrary precision exponent matrix */
    arb_mat_init(G, generator_count, expression_count);
    arb_mat_set_fmpz_mat(G, mat);

    /* compute the expression logs */
    arb_mat_init(expression_logs, expression_count, 1);
    for (i = 0; i < expression_count; i++)
    {
        expr_eval(x, expressions_table[i], level);
        arb_log(arb_mat_entry(expression_logs, i, 0), x, prec);

        /*
        flint_printf("expression %wd : ", i);
        arb_printd(x, 15);
        flint_printf("\n");
        */
    }

    /* compute the generator logs */
    arb_mat_init(generator_logs, generator_count, 1);
    arb_mat_mul(generator_logs, G, expression_logs, prec);

    /* fill a structure with corresponding double precision values */
    doublify_tkf91_generators(h, g, generator_logs);

    {
        clock_t diff_b;
        clock_t start_b = clock();

        tkf91_dynamic_programming_double(h, trace_flag, A, szA, B, szB);

        diff_b = clock() - start_b;
        int msec_b = (diff_b * 1000) / CLOCKS_PER_SEC;
        printf("Internal dynamic programming time taken %d seconds %d milliseconds.\n",
                msec_b/1000, msec_b%1000);
    }

    arb_clear(x);
    arb_mat_clear(G);
    arb_mat_clear(expression_logs);
    arb_mat_clear(generator_logs);
}
