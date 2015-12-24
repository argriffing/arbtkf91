/*
 * tkf91 dynamic programming cell bounds
 */

#include <time.h>

#include "mag.h"

#include "tkf91_dp_bound.h"
#include "breadcrumbs.h"


typedef struct
{
    mag_t m1_00;
    mag_t m0_10;
    mag_struct m0_i0_incr[4];
    mag_t m2_01;
    mag_struct m2_0j_incr[4];
    mag_struct c0_incr[4];
    mag_struct c1_incr[16];
    mag_struct c2_incr[4];
} tkf91_values_struct;
typedef tkf91_values_struct tkf91_values_t[1];


void bound_tkf91_generators(
        tkf91_values_struct low,
        tkf91_values_struct high,
        tkf91_generator_indices_t g,
        arb_mat_t m);

void tkf91_dynamic_programming_bound(
        tkf91_values_t lb,
        tkf91_values_t ub,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB);

void tkf91_values_init(tkf91_values_t x, tkf91_generator_indices_t g, mag * m);



void
tkf91_values_init(
        tkf91_values_t h,
        tkf91_generator_indices_t g,
        mag * m)
{
    slong i, j;
    mag_init_set(h->m1_00, m+g->m1_00);
    mag_init_set(h->m0_10, m+g->m0_10);
    mag_init_set(h->m2_01, m+g->m2_01);
    for (i = 0; i < 4; i++)
    {
        mag_init_set(h->m0_i0_incr+i, m+g->m0_i0_incr[i]);
        mag_init_set(h->m2_0j_incr+i, m+g->m2_0j_incr[i]);
        mag_init_set(h->c0_incr+i, m+g->c0_incr[i]);
        for (j = 0; j < 4; j++)
        {
            mag_init_set(h->c1_incr+i*4+j, m+g->c1_incr[i*4+j]);
        }
        mag_init_set(h->c2_incr+i, m+g->c2_incr[i]);
    }
}



void tkf91_dynamic_programming_double(tkf91_double_generators_t g,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB)
{
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
    /* TODO do this optionally... */

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
_arb_mat_get_col(arb * v, const arb_mat_t mat, slong j)
{
    slong i, nr;
    nr = arb_mat_nrows(mat);
    for (i = 0; i < nr; i++)
    {
        arb_set(v+i, arb_mat_entry(mat, i, j));
    }
}


void
_bounds_init(tkf91_values_t lb, tkf91_values_t ub,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g)
{
    slong nr, nc, level, prec;
    arb_mat_t G, U, V;
    arb_t x;
    arb_ptr v;
    mag_ptr lb_arr, ub_arr;

    /* set the precision level used before converting to 30 bit magnitudes */
    level = 8;
    prec = 1 << level;

    /* count the generators and expressions respectively */
    nr = fmpz_mat_nrows(mat);
    nc = fmpz_mat_ncols(mat);

    /* initialize the arbitrary precision integer exponent matrix */
    arb_mat_init(G, nr, nc);
    arb_mat_set_fmpz_mat(G, mat);

    /* compute logs of expressions to the specified precision */
    arb_mat_init(U, nc, 1);
    for (i = 0; i < nc; i++)
    {
        expr_eval(x, expressions_table[i], level);
        arb_log(arb_mat_entry(U, i, 0), x, prec);
    }

    /* compute logs of generators */
    arb_mat_init(V, nr, 1);
    arb_mat_mul(V, G, expression_logs, prec);

    /* copy the column vector to a new vector and exponentiate its entries */
    v = _arb_vec_init(nr);
    _arb_mat_get_col(v, mat, 0);
    for (i = 0; i < nr; i++)
    {
        arb_exp(v+i, v+i, prec);
    }

    /* create the lb and ub arrays */
    lb_arr = _mag_vec_init(nr);
    ub_arr = _mag_vec_init(nr);
    for (i = 0; i < nr; i++)
    {
        arb_get_mag_lower(lb_arr, v+i);
        arb_get_mag(ub_arr, v+i);
    }

    /* initialize the lb and ub structures */
    tkf91_values_init(lb, g, lb_arr);
    tkf91_values_init(ub, g, ub_arr);

    _mag_vec_clear(lb_arr);
    _mag_vec_clear(ub_arr);
    _arb_vec_clear(v);
    arb_mat_clear(G);
    arb_mat_clear(U);
    arb_mat_clear(V);
    arb_clear(x);
}


void
tkf91_dp_bound(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g,
        int trace_flag,
        slong *A, size_t szA,
        slong *B, size_t szB)
{

    /*
     * Convert the first few args to generator lower and upper bounds.
     */
    tkf91_values_t lb;
    tkf91_values_t ub;
    _bounds_init(lb, ub, mat, expressions_table, g);

    tkf91_values_clear(lb);
    tkf91_values_clear(ub);

    /* */


    arb_t x;
    arb_mat_t expression_logs;
    arb_mat_t generator_logs;
    slong i;
    slong generator_count = fmpz_mat_nrows(mat);
    slong expression_count = fmpz_mat_ncols(mat);
    tkf91_double_generators_t h;

    arb_init(x);

    /* initialize the arbitrary precision exponent matrix */
    arb_mat_t G;
    arb_mat_init(G, generator_count, expression_count);
    arb_mat_set_fmpz_mat(G, mat);

    /* compute the expression logs */
    arb_mat_init(expression_logs, expression_count, 1);
    for (i = 0; i < expression_count; i++)
    {
        expr_eval(x, expressions_table[i], level);
        arb_log(arb_mat_entry(expression_logs, i, 0), x, prec);
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
