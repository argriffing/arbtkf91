#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "flint/flint.h"
#include "flint/fmpq.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

#include "arb.h"
#include "arb_mat.h"

#include "femtocas.h"
#include "tkf91_rationals.h"
#include "expressions.h"
#include "generators.h"
#include "rgenerators.h"
#include "tkf91_generators.h"
#include "tkf91_rgenerators.h"
#include "tkf91_generator_indices.h"
#include "wavefront_double.h"
#include "wavefront_hermite.h"
#include "breadcrumbs.h"

#define MAXSEQLEN 20000

typedef struct
{
    fmpq_t lambda;
    fmpq_t mu;
    fmpq_t tau;
    fmpq pi[4];
} user_params_struct;
typedef user_params_struct user_params_t[1];

void user_params_init(user_params_t p);
void user_params_clear(user_params_t p);
void user_params_print(const user_params_t p);

void
user_params_init(user_params_t p)
{
    slong i;
    fmpq_init(p->lambda);
    fmpq_init(p->mu);
    fmpq_init(p->tau);
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
    fmpq_clear(p->tau);
    for (i = 0; i < 4; i++)
    {
        fmpq_clear(p->pi+i);
    }
}

void
user_params_print(const user_params_t p)
{
    flint_printf("lambda: "); fmpq_print(p->lambda); flint_printf("\n");
    flint_printf("mu: "); fmpq_print(p->mu); flint_printf("\n");
    flint_printf("tau: "); fmpq_print(p->tau); flint_printf("\n");
    flint_printf("pa: "); fmpq_print(p->pi+0); flint_printf("\n");
    flint_printf("pc: "); fmpq_print(p->pi+1); flint_printf("\n");
    flint_printf("pg: "); fmpq_print(p->pi+2); flint_printf("\n");
    flint_printf("pt: "); fmpq_print(p->pi+3); flint_printf("\n");
}




/*
 * Hermitification of tkf91 generators.
 * The members of this structure will point to rows of an fmpz matrix.
 * If my ideas are right, the fmpz matrix will be the inverse of the
 * transform matrix associated with the hermite normal form transformation
 * of the generator matrix.
 * Each of the vectors will have length equal to the rank of the Hermite 
 * normal form, but the 'hermitification' procedures do not need to
 * know or care about this length.
 */
typedef struct
{
    fmpz * m1_00;
    fmpz * m0_10;
    fmpz * m0_i0_incr[4];
    fmpz * m2_01;
    fmpz * m2_0j_incr[4];
    fmpz * c0_incr[4];
    fmpz * c1_incr[16];
    fmpz * c2_incr[4];
} tkf91_hermite_generators_struct;
typedef tkf91_hermite_generators_struct tkf91_hermite_generators_t[1];

fmpz * _hermitify(slong i, fmpz_mat_t M);
void hermitify_tkf91_generators(
        tkf91_hermite_generators_t h,
        tkf91_generator_indices_t g,
        fmpz_mat_t M);

fmpz *
_hermitify(slong i, fmpz_mat_t M)
{
    return M->rows[i];
}

void
hermitify_tkf91_generators(
        tkf91_hermite_generators_t h,
        tkf91_generator_indices_t g,
        fmpz_mat_t M)
{
    slong i, j;
    h->m1_00 = _hermitify(g->m1_00, M);
    h->m0_10 = _hermitify(g->m0_10, M);
    h->m2_01 = _hermitify(g->m2_01, M);
    for (i = 0; i < 4; i++)
    {
        h->m0_i0_incr[i] = _hermitify(g->m0_i0_incr[i], M);
        h->m2_0j_incr[i] = _hermitify(g->m2_0j_incr[i], M);
        h->c0_incr[i] = _hermitify(g->c0_incr[i], M);
        for (j = 0; j < 4; j++)
        {
            h->c1_incr[i*4+j] = _hermitify(g->c1_incr[i*4+j], M);
        }
        h->c2_incr[i] = _hermitify(g->c2_incr[i], M);
    }
}



/*
 * "doublification" of tkf91 generators
 */

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





/*
 * This helper function converts a string to a list of indices.
 */

void _fill_sequence_vector(slong *v, const char *str, slong n);

void
_fill_sequence_vector(slong *v, const char *str, slong n)
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
                       {
                           flint_printf("unrecognized nucleotide\n");
                           abort();
                       }
        }
    }
}


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




void tkf91_dynamic_programming_double(tkf91_double_generators_t g,
        int trace_flag,
        slong *A, slong szA,
        slong *B, slong szB);

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
tkf91_double_precision(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g,
        int trace_flag,
        slong *A, size_t szA,
        slong *B, size_t szB);

void
tkf91_double_precision(
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






/*
 * Dynamic programming with some symbolic and numeric cleverness.
 * But undoubtedly much slower than the variant that uses double precision.
 * See that less sophisticated variant for more details.
 */

void
_arb_vec_dot_fmpz_vec(
        arb_t res, arb_srcptr vec1, const fmpz * vec2, slong len2, slong prec);

void
_arb_vec_dot_fmpz_vec(
        arb_t res, arb_srcptr vec1, const fmpz * vec2, slong len2, slong prec)
{
    slong i;
    arb_zero(res);
    for (i = 0; i < len2; i++)
    {
        arb_addmul_fmpz(res, vec1+i, vec2+i, prec);
    }
}


void
_hwave_element_add_vec(
        hwave_element_t e, const fmpz * vec1, const fmpz * vec2, arb_srcptr v,
        slong r, slong prec);

void
_hwave_element_set_vec(
        hwave_element_t e, const fmpz * vec, arb_srcptr v,
        slong r, slong prec);

void
_hwave_element_add_vec(
        hwave_element_t e, const fmpz * vec1, const fmpz * vec2, arb_srcptr v,
        slong r, slong prec)
{
    /*
    flint_printf("adding coefficients\n");
    flint_printf("  "); _fmpz_vec_print(vec1, r); flint_printf("\n");
    flint_printf("+ "); _fmpz_vec_print(vec2, r); flint_printf("\n");
    */

    _fmpz_vec_add(e->vec, vec1, vec2, r);
    _arb_vec_dot_fmpz_vec(e->value, v, e->vec, r, prec);
    e->status = HWAVE_STATUS_UNAMBIGUOUS;

    /*
    flint_printf("= "); _fmpz_vec_print(e->vec, r); flint_printf("\n\n");
    arb_t a, b;
    arb_init(a);
    arb_init(b);
    _arb_vec_dot_fmpz_vec(a, v, vec1, r, prec);
    _arb_vec_dot_fmpz_vec(b, v, vec2, r, prec);
    arb_printd(a, 15); flint_printf("\n");
    arb_printd(b, 15); flint_printf("\n");
    arb_clear(a);
    arb_clear(b);
    flint_printf("yielding function value ");
    arb_printd(e->value, 15);
    flint_printf("\n");
    */
}

void
_hwave_element_set_vec(
        hwave_element_t e, const fmpz * vec, arb_srcptr v,
        slong r, slong prec)
{
    /* Input:
     *  e : one of the three elements of an hwave cell
     *  vec : integer vector of r coefficients
     *  v : arb vector of r real balls
     *  r : the rank of the system
     *  prec : precision for arb
     */

    /* copy the integer coefficients */
    _fmpz_vec_set(e->vec, vec, r);

    /* update the dot product */
    _arb_vec_dot_fmpz_vec(e->value, v, e->vec, r, prec);

    /*
     * Mark that for this element the vector of integer coefficients
     * is known unambiguously.
     */
    e->status = HWAVE_STATUS_UNAMBIGUOUS;
}

hwave_element_ptr
_get_max3_checked(
        hwave_element_ptr e0,
        hwave_element_ptr e1,
        hwave_element_ptr e2,
        slong len);

hwave_element_ptr
_get_max3_checked(
        hwave_element_ptr e0,
        hwave_element_ptr e1,
        hwave_element_ptr e2,
        slong len)
{
    /*
     * FIXME
     * For now we abort if an ambiguity is detected numerically
     * but not symbolically.
     */
    /* find the element with the greatest midpoint */
    slong i;
    hwave_element_ptr x[] = {e0, e1, e2};

    /* for now require that all elements have unambiguous status */
    for (i = 0; i < 3; i++)
    {
        if (x[i]->status == HWAVE_STATUS_AMBIGUOUS)
        {
            flint_printf("ambiguous status -- ");
            flint_printf("handling this case is not yet implemented...\n");
            abort();
        }
    }

    /* determine the element with greatest midpoint */
    hwave_element_ptr best;
    best = e0;
    if (arf_cmp(arb_midref(best->value), arb_midref(e1->value)) < 0) {
        best = e1;
    }
    if (arf_cmp(arb_midref(best->value), arb_midref(e2->value)) < 0) {
        best = e2;
    }

    /*
     * If any element's value overlaps the value of the element
     * whose midpoint is greatest, then complain if the integer coefficient
     * vectors are not equal.
     */
    for (i = 0; i < 3; i++)
    {
        if (arb_overlaps(best->value, x[i]->value))
        {
            if (!_fmpz_vec_equal(best->vec, x[i]->vec, len))
            {
                flint_printf("an ambiguity was detected numericaly ");
                flint_printf("but not symbolically -- ");
                flint_printf("dealing with this situation ");
                flint_printf("is not yet implemented...\n");
                abort();
            }
        }
    }

    return best;
}

void tkf91_dynamic_programming_hermite(tkf91_hermite_generators_t g,
        arb_ptr v, slong rank, slong prec,
        slong *A, slong szA,
        slong *B, slong szB);

void tkf91_dynamic_programming_hermite(tkf91_hermite_generators_t g,
        arb_ptr v, slong rank, slong prec,
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
    breadcrumb_mat_init(crumb_mat, nrows, ncols);

    /* define the wavefront matrix */
    hwave_mat_t wave;
    slong modulus = 3;
    /* slong modulus = nrows + ncols - 1; */
    hwave_mat_init(wave, nrows + ncols - 1, modulus, rank);

    /* iterate over anti-diagonal bands of the dynamic programming matrix */
    hwave_element_ptr best;
    hwave_cell_ptr cell, p0, p1, p2;
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

            cell = hwave_mat_entry(wave, k, l);
            if (i == 0 && j == 0)
            {
                hwave_element_set_undefined(cell->m+0);
                _hwave_element_set_vec(cell->m+1, g->m1_00, v, rank, prec);
                hwave_element_set_undefined(cell->m+2);
            }
            else if (i == 1 && j == 0)
            {
                _hwave_element_set_vec(cell->m+0, g->m0_10, v, rank, prec);
                hwave_element_set_undefined(cell->m+1);
                hwave_element_set_undefined(cell->m+2);
            }
            else if (i == 0 && j == 1)
            {
                hwave_element_set_undefined(cell->m+0);
                hwave_element_set_undefined(cell->m+1);
                _hwave_element_set_vec(cell->m+2, g->m2_01, v, rank, prec);
            }
            else
            {
                if (i == 0)
                {
                    ntb = B[j - 1];
                    p2 = hwave_mat_entry_left(wave, k, l);
                    hwave_element_set_undefined(cell->m+0);
                    hwave_element_set_undefined(cell->m+1);
                    if (p2->m[2].status != HWAVE_STATUS_UNAMBIGUOUS)
                    {
                        flint_printf("found a not unambiguous status ");
                        flint_printf("while filling the first row\n");
                        abort();
                    }
                    _hwave_element_add_vec(
                            cell->m+2,
                            p2->m[2].vec,
                            g->m2_0j_incr[ntb],
                            v, rank, prec);
                }
                else if (j == 0)
                {
                    nta = A[i - 1];
                    p0 = hwave_mat_entry_top(wave, k, l);
                    if (p0->m[0].status != HWAVE_STATUS_UNAMBIGUOUS)
                    {
                        flint_printf("found a not unambiguous status ");
                        flint_printf("while filling the first column\n");
                        abort();
                    }
                    _hwave_element_add_vec(
                            cell->m+0,
                            p0->m[0].vec,
                            g->m0_i0_incr[nta],
                            v, rank, prec);
                    hwave_element_set_undefined(cell->m+1);
                    hwave_element_set_undefined(cell->m+2);
                }
                else
                {
                    nta = A[i - 1];
                    ntb = B[j - 1];
                    p0 = hwave_mat_entry_top(wave, k, l);
                    p1 = hwave_mat_entry_diag(wave, k, l);
                    p2 = hwave_mat_entry_left(wave, k, l);

                    best = _get_max3_checked(p0->m+0, p0->m+1, p0->m+2, rank);
                    if (best->status == HWAVE_STATUS_UNDEFINED)
                    {
                        hwave_element_set_undefined(cell->m+0);
                    }
                    else
                    {
                        _hwave_element_add_vec(
                                cell->m+0,
                                best->vec,
                                g->c0_incr[nta],
                                v, rank, prec);
                    }

                    best = _get_max3_checked(p1->m+0, p1->m+1, p1->m+2, rank);
                    if (best->status == HWAVE_STATUS_UNDEFINED)
                    {
                        hwave_element_set_undefined(cell->m+1);
                    }
                    else
                    {
                        _hwave_element_add_vec(
                                cell->m+1,
                                best->vec,
                                g->c1_incr[nta*4 + ntb],
                                v, rank, prec);
                    }

                    /* use max3 to check only max2 by duplicating one -- */
                    /* passing p2->m+1 twice is not a bug */
                    /* TODO make a max2 convenience function for this... */
                    best = _get_max3_checked(p2->m+1, p2->m+1, p2->m+2, rank);
                    if (best->status == HWAVE_STATUS_UNDEFINED)
                    {
                        hwave_element_set_undefined(cell->m+2);
                    }
                    else
                    {
                        _hwave_element_add_vec(
                                cell->m+2,
                                best->vec,
                                g->c2_incr[ntb],
                                v, rank, prec);
                    }
                }
            }

            /* fill the table for traceback */
            /* FIXME use the breadcrumb ambiguity bit as appropriate */
            breadcrumb_ptr pcrumb = breadcrumb_mat_entry(crumb_mat, i, j);
            best = _get_max3_checked(cell->m+0, cell->m+1, cell->m+2, rank);
            if (_fmpz_vec_equal(best->vec, cell->m[0].vec, rank))
            {
                *pcrumb |= CRUMB_TOP;
            }
            if (_fmpz_vec_equal(best->vec, cell->m[1].vec, rank))
            {
                *pcrumb |= CRUMB_DIAG;
            }
            if (_fmpz_vec_equal(best->vec, cell->m[2].vec, rank))
            {
                *pcrumb |= CRUMB_LEFT;
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
    slong c;
    arf_ptr mid;
    for (c = 0; c < 3; c++)
    {
        flint_printf("m%wd:\n", c);
        for (i = 0; i < nrows; i++)
        {
            for (j = 0; j < ncols; j++)
            {
                k = i + j;
                l = nrows - 1 + j - i;
                cell = hwave_mat_entry(wave, k, l);
                mid = arb_midref(cell->m[c].value);
                flint_printf("%.17lf ", exp(arf_get_d(mid, ARF_RND_NEAR)));
            }
            flint_printf("\n");
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
    cell = hwave_mat_entry(wave, k, l);
    best = _get_max3_checked(cell->m+0, cell->m+1, cell->m+2, rank);
    arb_t arbscore;
    arb_init(arbscore);
    _arb_vec_dot_fmpz_vec(arbscore, v, best->vec, rank, prec);
    arb_exp(arbscore, arbscore, prec);
    flint_printf("score: ");
    arb_printd(arbscore, 15);
    flint_printf("\n");
    arb_clear(arbscore);

    /* do the traceback */
    char *sa, *sb;
    breadcrumb_mat_get_alignment(&sa, &sb, crumb_mat, A, B);
    flint_printf("%s\n", sa);
    flint_printf("%s\n", sb);
    flint_printf("\n");
    free(sa);
    free(sb);

    /* clear the tables */
    hwave_mat_clear(wave);
    breadcrumb_mat_clear(crumb_mat);
}


void
tkf91_hermite(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g,
        slong *A, size_t szA,
        slong *B, size_t szB);

void
tkf91_hermite(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        tkf91_generator_indices_t g,
        slong *A, size_t szA,
        slong *B, size_t szB)
{
    /* Inputs:
     *   mat : the generator matrix -- mat_ij where i is a generator index
     *         and j is an expression index.
     *   expressions_table : map from expression index to expression object
     *   g : a struct with tkf91 generator indices
     */

    /*
     * For now, use an arbitrary precision.
     * In later versions of this program,
     * we could decide to adjust it dynamically if it is detected
     * to be insufficient.
     */
    slong level = 6;
    slong prec = 1 << level;

    arb_t x;
    arb_init(x);

    /* Compute a Hermite decomposition of the generator matrix. */
    /* U * mat = H */
    fmpz_mat_t H, U;
    fmpz_mat_init(H, fmpz_mat_nrows(mat), fmpz_mat_ncols(mat));
    fmpz_mat_init(U, fmpz_mat_nrows(mat), fmpz_mat_nrows(mat));
    fmpz_mat_hnf_transform(H, U, mat);

    /*
     * The number of nonzero rows of H should determine the rank of mat.
     * In H, rows with nonzero entries should precede rows without 
     * nonzero entries
     */
    slong i;
    slong rank = 0;
    for (i = 0; i < fmpz_mat_nrows(H); i++)
    {
        if (!fmpz_mat_is_zero_row(H, i))
        {
            if (i != rank)
            {
                flint_printf("expected each row in H ");
                flint_printf("containing a nonzero entry ");
                flint_printf("to precede each row containing only zeros\n");
                abort();
            }
            rank++;
        }
    }
    flint_printf("matrix rank ");
    flint_printf("as revealed by the Hermite normal form: %wd\n", rank);

    /*
     * Invert the transform matrix U.
     * This should be possible because U should be unimodular
     * and therefore nonsingular.
     * The inverse should have entries that are integers,
     * and its determinant should be 1.
     *
     * We only care about the first r columns of this inverse,
     * where r is the rank of the Hermite form H.
     */
    fmpz_mat_t V;
    int result;
    fmpz_t den;
    fmpz_mat_init(V, fmpz_mat_nrows(U), fmpz_mat_ncols(U));
    fmpz_init(den);
    result = fmpz_mat_inv(V, den, U);
    if (!result)
    {
        flint_printf("expected U to be nonsingular\n");
        abort();
    }
    if (!fmpz_is_one(den))
    {
        fmpz_t negden;
        fmpz_init(negden);
        fmpz_neg(negden, den);
        if (fmpz_is_one(negden))
        {
            fmpz_mat_neg(V, V);
        }
        else
        {
            flint_printf("expected U to be unimodular -- ");
            flint_printf("denominator of inverse of U: ");
            fmpz_print(den);
            flint_printf("\n");
            abort();
        }
        fmpz_clear(negden);
    }
    fmpz_clear(den);

    /*
     * 'hermitify' a structure with tkf91 generators,
     * by pointing the member variables to rows of V.
     */
    tkf91_hermite_generators_t h;
    hermitify_tkf91_generators(h, g, V);

    /*
     * Initialize an arbitrary precision matrix.
     * This will be part of the calculation H*log(y)
     * which gives the score of each quasi-generator.
     */
    arb_mat_t arbH;
    arb_mat_init(arbH, fmpz_mat_nrows(H), fmpz_mat_ncols(H));
    arb_mat_set_fmpz_mat(arbH, H);

    /*
     * Compute the expression logs.
     * This is like a log(y) column vector.
     */
    arb_mat_t expression_logs;
    arb_mat_init(expression_logs, fmpz_mat_ncols(mat), 1);
    for (i = 0; i < fmpz_mat_ncols(mat); i++)
    {
        expr_eval(x, expressions_table[i], level);
        arb_log(arb_mat_entry(expression_logs, i, 0), x, prec);

        /*
        flint_printf("expression %wd : ", i);
        arb_printd(x, 15);
        flint_printf("\n");
        */
    }

    /*
     * Compute the quasi-generator logs H*log(y).
     */
    arb_mat_t quasi_generator_logs;
    arb_mat_init(quasi_generator_logs, fmpz_mat_nrows(mat), 1);
    arb_mat_mul(quasi_generator_logs, arbH, expression_logs, prec);

    /*
     * Create an arb row vector by taking the first r elements
     * from the quasi_generator_logs column vector,
     * where r is the matrix rank of the generator matrix.
     */
    arb_ptr v;
    v = _arb_vec_init(rank);
    flint_printf("quasi generator logarithms:\n");
    for (i = 0; i < rank; i++)
    {
        arb_set(v+i, arb_mat_entry(quasi_generator_logs, i, 0));
        arb_printd(v+i, 15);
        flint_printf("\n");
    }
    flint_printf("\n");

    /*
     * Begin checking some invariants.
     * Require that the result of the calculation does not depend on the basis.
     *
     * Initialize an arbitrary precision matrix.
     * This will be part of the calculation G*log(y)
     * which gives the score of each generator.
     */
    arb_mat_t arbG;
    arb_mat_init(arbG, fmpz_mat_nrows(mat), fmpz_mat_ncols(mat));
    arb_mat_set_fmpz_mat(arbG, mat);

    /*
     * Compute the generator logs G*log(y).
     */
    arb_mat_t generator_logs;
    arb_mat_init(generator_logs, fmpz_mat_nrows(mat), 1);
    arb_mat_mul(generator_logs, arbG, expression_logs, prec);

    /*
    flint_printf("log probability for each generator:\n");
    arb_mat_printd(generator_logs, 15);
    flint_printf("\n");
    */

    /*
     * Compute dot products which should be equivalent.
     */
    /*
    flint_printf("a presumably equivalent calculation in a different basis\n");
    for (i = 0; i < fmpz_mat_nrows(mat); i++)
    {
        _arb_vec_dot_fmpz_vec(x, v, V->rows[i], rank, prec);

        flint_printf("first %wd columns of row %wd of V : ", rank, i);
        flint_printf("\n");
        _fmpz_vec_print(V->rows[i], rank);
        flint_printf("dot product : ");

        arb_printd(x, 15);
        flint_printf("\n");
    }
    flint_printf("\n");
    */

    /*
    flint_printf("V:\n"); fmpz_mat_print_pretty(V); flint_printf("\n");
    flint_printf("H:\n"); fmpz_mat_print_pretty(H); flint_printf("\n");
    */

    /*
     * Clear temporary variables.
     * Keep V because member variables of h point to its rows.
     * Keep v because it is used directly in the next stage.
     */
    arb_clear(x);
    fmpz_mat_clear(H);
    fmpz_mat_clear(U);
    arb_mat_clear(arbH);
    arb_mat_clear(arbG);
    arb_mat_clear(expression_logs);
    arb_mat_clear(generator_logs);
    arb_mat_clear(quasi_generator_logs);

    /* do the thing */
    tkf91_dynamic_programming_hermite(h, v, rank, prec, A, szA, B, szB);

    /* Clear the remaining variables. */
    fmpz_mat_clear(V);
    _arb_vec_clear(v, rank);
}



void run(const char *strA, const char *strB, const user_params_t params);

void
run(const char *strA, const char *strB, const user_params_t params)
{
    slong *A;
    slong *B;
    expr_ptr * expressions_table;

    size_t szA, szB;
    reg_t reg;
    tkf91_rationals_t r;
    tkf91_expressions_t p;
    tkf91_generator_indices_t g;
    fmpz_mat_t mat;

    szA = strlen(strA);
    A = flint_malloc(szA * sizeof(slong));
    _fill_sequence_vector(A, strA, szA);

    szB = strlen(strB);
    B = flint_malloc(szB * sizeof(slong));
    _fill_sequence_vector(B, strB, szB);

    int i;
    int runs = 1;


    for (i = 0; i < runs; i++)
    {
        clock_t diff;
        clock_t start = clock();

        reg_init(reg);
        tkf91_rationals_init(r,
                params->lambda, params->mu, params->tau, params->pi);
        tkf91_expressions_init(p, reg, r);

        rgen_reg_ptr rg = rgen_reg_new();
        tkf91_rgenerators_init(g, rg, r, p, A, szA, B, szB);
        rgen_reg_finalize(rg, reg);
        fmpz_mat_init(mat, rgen_reg_nrows(rg), rgen_reg_ncols(rg));
        rgen_reg_get_matrix(mat, rg);
        rgen_reg_clear(rg);

        expressions_table = reg_vec(reg);

        /* TODO use the trace flag provided by the user */
        int trace_flag = 1;

        /*
        tkf91_double_precision(mat, expressions_table, g, trace_flag,
                A, szA, B, szB);
        */

        tkf91_hermite(mat, expressions_table, g, A, szA, B, szB);

        fmpz_mat_clear(mat);
        flint_free(expressions_table);

        reg_clear(reg);
        tkf91_rationals_clear(r);
        tkf91_expressions_clear(p);

        diff = clock() - start;
        int msec = (diff * 1000) / CLOCKS_PER_SEC;
        printf("Time taken %d seconds %d milliseconds.\n",
                msec/1000, msec%1000);
    }

    flint_free(A);
    flint_free(B);
}


void bench(const user_params_t params, int trace_flag);

void
bench(const user_params_t params, int trace_flag)
{
    char strA[MAXSEQLEN];
    char strB[MAXSEQLEN];
    size_t szA, szB;

    slong *A;
    slong *B;
    expr_ptr * expressions_table;

    reg_t reg;
    tkf91_rationals_t r;
    tkf91_expressions_t p;
    tkf91_generator_indices_t g;
    fmpz_mat_t mat;

    /*
     * Read pairs of sequences from stdin.
     * Assume one sequence per line.
     */
    while (1)
    {

        printf("waiting for first sequence...\n");
        if (!fgets(strA, MAXSEQLEN, stdin))
        {
            break;
        }
        szA = strlen(strA);
        szA--; /* do not count the newline */
        if (!szA)
        {
            break;
        }
        printf("length of sequence A: %ld\n", szA);

        printf("waiting for second sequence...\n");
        if (!fgets(strB, MAXSEQLEN, stdin))
        {
            break;
        }
        szB = strlen(strB);
        szB--; /* do not count the newline */
        if (!szB)
        {
            break;
        }
        printf("length of sequence B: %ld\n", szB);


        A = flint_malloc(szA * sizeof(slong));
        _fill_sequence_vector(A, strA, szA);

        B = flint_malloc(szB * sizeof(slong));
        _fill_sequence_vector(B, strB, szB);


        clock_t diff;
        clock_t start = clock();

        reg_init(reg);
        tkf91_rationals_init(r,
                params->lambda, params->mu, params->tau, params->pi);
        tkf91_expressions_init(p, reg, r);

        rgen_reg_ptr rg = rgen_reg_new();
        tkf91_rgenerators_init(g, rg, r, p, A, szA, B, szB);
        rgen_reg_finalize(rg, reg);
        fmpz_mat_init(mat, rgen_reg_nrows(rg), rgen_reg_ncols(rg));
        rgen_reg_get_matrix(mat, rg);
        rgen_reg_clear(rg);

        expressions_table = reg_vec(reg);

        clock_t diff_b;
        clock_t start_b = clock();

        tkf91_double_precision(mat, expressions_table, g, trace_flag,
                A, szA, B, szB);
        /*
        tkf91_hermite(mat, expressions_table, g, A, szA, B, szB);
        */

        diff_b = clock() - start_b;
        int msec_b = (diff_b * 1000) / CLOCKS_PER_SEC;
        printf("Dynamic programming time taken %d seconds %d milliseconds.\n",
                msec_b/1000, msec_b%1000);


        fmpz_mat_clear(mat);
        flint_free(expressions_table);

        reg_clear(reg);
        tkf91_rationals_clear(r);
        tkf91_expressions_clear(p);

        diff = clock() - start;
        int msec = (diff * 1000) / CLOCKS_PER_SEC;
        printf("Total time taken %d seconds %d milliseconds.\n",
                msec/1000, msec%1000);

        flint_free(A);
        flint_free(B);
    }
}




int
main(int argc, char *argv[])
{
    int i;

    const char *Astr = NULL;
    const char *Bstr = NULL;

    user_params_t p;
    user_params_init(p);

    /* indicates benchmark mode where sequence pairs are read from stdin */
    int bench_flag = 0;
    int trace_flag = 0;

    slong lambda_num = 0;
    slong mu_num = 0;
    slong tau_num = 0;
    slong pi_num[4] = {0, 0, 0, 0};

    slong lambda_den = 1;
    slong mu_den = 1;
    slong tau_den = 1;
    slong pi_den[4] = {1, 1, 1, 1};

    for (i = 1; i < argc-1; i += 2)
    {
        if (strcmp(argv[i], "--bench") == 0) {
            flint_sscanf(argv[i + 1], "%d", &bench_flag);
        } else if (strcmp(argv[i], "--trace") == 0) {
            flint_sscanf(argv[i + 1], "%d", &trace_flag);
        }

        else if (strcmp(argv[i], "--sequence-1") == 0) {
            Astr = argv[i + 1];
        } else if (strcmp(argv[i], "--sequence-2") == 0) {
            Bstr = argv[i + 1];
        }
        
        else if (strcmp(argv[i], "--lambda-num") == 0) {
            flint_sscanf(argv[i + 1], "%ws", &lambda_num);
        } else if (strcmp(argv[i], "--mu-num") == 0) {
            flint_sscanf(argv[i + 1], "%ws", &mu_num);
        } else if (strcmp(argv[i], "--tau-num") == 0) {
            flint_sscanf(argv[i + 1], "%ws", &tau_num);
        } else if (strcmp(argv[i], "--pa-num") == 0) {
            flint_sscanf(argv[i + 1], "%ws", pi_num+0);
        } else if (strcmp(argv[i], "--pc-num") == 0) {
            flint_sscanf(argv[i + 1], "%ws", pi_num+1);
        } else if (strcmp(argv[i], "--pg-num") == 0) {
            flint_sscanf(argv[i + 1], "%ws", pi_num+2);
        } else if (strcmp(argv[i], "--pt-num") == 0) {
            flint_sscanf(argv[i + 1], "%ws", pi_num+3);
        }

        else if (strcmp(argv[i], "--lambda-den") == 0) {
            flint_sscanf(argv[i + 1], "%ws", &lambda_den);
        } else if (strcmp(argv[i], "--mu-den") == 0) {
            flint_sscanf(argv[i + 1], "%ws", &mu_den);
        } else if (strcmp(argv[i], "--tau-den") == 0) {
            flint_sscanf(argv[i + 1], "%ws", &tau_den);
        } else if (strcmp(argv[i], "--pa-den") == 0) {
            flint_sscanf(argv[i + 1], "%ws", pi_den+0);
        } else if (strcmp(argv[i], "--pc-den") == 0) {
            flint_sscanf(argv[i + 1], "%ws", pi_den+1);
        } else if (strcmp(argv[i], "--pg-den") == 0) {
            flint_sscanf(argv[i + 1], "%ws", pi_den+2);
        } else if (strcmp(argv[i], "--pt-den") == 0) {
            flint_sscanf(argv[i + 1], "%ws", pi_den+3);
        }
    }

    fmpq_set_si(p->lambda, lambda_num, lambda_den);
    fmpq_set_si(p->mu, mu_num, mu_den);
    fmpq_set_si(p->tau, tau_num, tau_den);
    for (i = 0; i < 4; i++)
    {
        fmpq_set_si(p->pi+i, pi_num[i], pi_den[i]);
    }

    flint_printf("user-provided parameter values:\n");
    user_params_print(p);
    flint_printf("\n");

    if (bench_flag)
    {
        bench(p, trace_flag);
    }
    else
    {
        run(Astr, Bstr, p);
    }

    user_params_clear(p);

    flint_cleanup();
    return 0;
}
