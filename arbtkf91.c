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
#include "expressions.h"
#include "generators.h"
#include "wavefront_double.h"

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
} named_double_generators_struct;
typedef named_double_generators_struct named_double_generators_t[1];

double _doublify(slong i, arb_mat_t m);
void doublify_named_generators(
        named_double_generators_t h,
        named_generators_t g,
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
doublify_named_generators(
        named_double_generators_t h,
        named_generators_t g,
        arb_mat_t m)
{
    slong i, j;
    h->m1_00 = _doublify(g->m1_00, m);
    h->m0_10 = _doublify(g->m0_10, m);
    h->m2_01 = _doublify(g->m2_01, m);
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            h->m0_i0_incr[i] = _doublify(g->m0_i0_incr[i], m);
            h->m2_0j_incr[i] = _doublify(g->m2_0j_incr[i], m);
            h->c0_incr[i] = _doublify(g->c0_incr[i], m);
            h->c1_incr[i*4+j] = _doublify(g->c1_incr[i*4+j], m);
            h->c2_incr[i] = _doublify(g->c2_incr[i], m);
        }
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
               abort();
        }
    }
}


double max2(double a, double b);
double max3(double a, double b, double c);

double
max2(double a, double b)
{
    return a > b ? a : b;
}

double
max3(double a, double b, double c)
{
    return max2(a, max2(b, c));
}




/* the matrix for the traceback stage of dynamic programming */

/*
 * Breadcrumbs for the traceback stage of dynamic programming may indicate
 * that multiple paths have identical scores.
 * To show that multiple traceback continuations from a partial solution
 * are equally valid, use a bitwise combination of these flags.
 */
#define CRUMB_TOP  0b00000001
#define CRUMB_DIAG 0b00000010
#define CRUMB_LEFT 0b00000100

typedef unsigned char breadcrumb_t;
typedef breadcrumb_t * breadcrumb_ptr;

typedef struct
{
    breadcrumb_t *data;
    slong nrows;
    slong ncols;
} breadcrumb_mat_struct;
typedef breadcrumb_mat_struct breadcrumb_mat_t[1];

void breadcrumb_mat_init(breadcrumb_mat_t mat, slong nrows, slong ncols);
void breadcrumb_mat_clear(breadcrumb_mat_t mat);
breadcrumb_ptr breadcrumb_mat_entry(breadcrumb_mat_t mat, slong i, slong j);
void breadcrumb_mat_get_alignment(char **psa, char **psb,
        breadcrumb_mat_t mat, const slong *A, const slong *B);

void
breadcrumb_mat_init(breadcrumb_mat_t mat, slong nrows, slong ncols)
{
    mat->data = calloc(nrows * ncols, sizeof(breadcrumb_t));
    mat->nrows = nrows;
    mat->ncols = ncols;
}

void
breadcrumb_mat_clear(breadcrumb_mat_t mat)
{
    free(mat->data);
}

breadcrumb_ptr
breadcrumb_mat_entry(breadcrumb_mat_t mat, slong i, slong j)
{
    if (i < 0 || i >= mat->nrows || j < 0 || j >= mat->ncols)
    {
        flint_printf("breadcrumb matrix indexing error\n");
        flint_printf("i=%wd j=%wd\n", i, j);
        flint_printf("nrows=%wd ncols=%wd\n", mat->nrows, mat->ncols);
        abort();
    }
    return mat->data + i * mat->ncols + j;
}

void
breadcrumb_mat_get_alignment(char **psa, char **psb,
        breadcrumb_mat_t mat, const slong *A, const slong *B)
{
    /* do the traceback */
    slong i, j;
    char ACGT[4] = "ACGT";
    slong n = mat->nrows + mat->ncols;
    char *sa = calloc(n, sizeof(char));
    char *sb = calloc(n, sizeof(char));
    slong len = 0;
    i = mat->nrows - 1;
    j = mat->ncols - 1;
    breadcrumb_t crumb;
    while (i > 0 || j > 0)
    {
        crumb = *breadcrumb_mat_entry(mat, i, j);
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
    char tmp;
    for (i = 0; i < len/2; i++)
    {
        j = len - 1 - i;
        tmp = sa[i]; sa[i] = sa[j]; sa[j] = tmp;
        tmp = sb[i]; sb[i] = sb[j]; sb[j] = tmp;
    }
    *psa = sa;
    *psb = sb;
}





void tkf91_dynamic_programming(named_double_generators_t g,
        slong *A, slong szA,
        slong *B, slong szB);

void tkf91_dynamic_programming(named_double_generators_t g,
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

    breadcrumb_mat_init(crumb_mat, nrows, ncols);

    /* define the wavefront matrix */
    wave_mat_t wave;
    slong modulus = nrows + ncols - 1;
    /* slong modulus = 3; */
    wave_mat_init(wave, nrows + ncols - 1, modulus);

    /*
     * Let M_{ij} be the matrix created for traceback.
     * Let R_{kl} be the 'logical' wavefront matrix.
     * Then
     *  k = i + j
     *  l = nrows - 1 + j - i
     */


    /* iterate over diagonals of the dynamic programming matrix */
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

            cell = wave_mat_entry(wave, k, l);
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
            }

            /* fill the table for traceback */
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
    char *sa, *sb;
    breadcrumb_mat_get_alignment(&sa, &sb, crumb_mat, A, B);
    flint_printf("%s\n", sa);
    flint_printf("%s\n", sb);
    flint_printf("\n");
    free(sa);
    free(sb);

    /* clear the tables */
    wave_mat_clear(wave);
    breadcrumb_mat_clear(crumb_mat);
}



void mess_with_generator_matrix(fmpz_mat_t A);

void
mess_with_generator_matrix(fmpz_mat_t A)
{
    fmpz_t den;
    fmpz_mat_t B, U, V, H, R;

    fmpz_init(den);

    fmpz_mat_init(B, fmpz_mat_ncols(A), fmpz_mat_nrows(A));
    fmpz_mat_init(U, fmpz_mat_nrows(B), fmpz_mat_nrows(B));
    fmpz_mat_init(V, fmpz_mat_nrows(B), fmpz_mat_nrows(B));
    fmpz_mat_init(H, fmpz_mat_nrows(B), fmpz_mat_ncols(B));
    fmpz_mat_init(R, fmpz_mat_nrows(A), fmpz_mat_nrows(A));

    fmpz_mat_transpose(B, A);

    /* hermite transform of B */
    fmpz_mat_hnf_transform(H, U, B);
    flint_printf("U * B = H\n");
    flint_printf("U:\n"); fmpz_mat_print_pretty(U);
    flint_printf("B:\n"); fmpz_mat_print_pretty(B);
    flint_printf("H:\n"); fmpz_mat_print_pretty(H);
    flint_printf("\n\n");

    /* inverse of U */
    fmpz_mat_inv(V, den, U);
    flint_printf("B = d*U^-1 * (1/d)H\n");
    flint_printf("B:\n"); fmpz_mat_print_pretty(B); flint_printf("\n");
    flint_printf("d: "); fmpz_print(den); flint_printf("\n");
    flint_printf("d * U^-1:\n"); fmpz_mat_print_pretty(V);
    flint_printf("H:\n"); fmpz_mat_print_pretty(H);
    flint_printf("\n\n");

    /* R */
    fmpz_mat_mul(R, A, V);
    flint_printf("A * (d*U^-1):\n");
    fmpz_mat_print_pretty(R);
    flint_printf("\n\n");

    fmpz_clear(den);

    fmpz_mat_clear(B);
    fmpz_mat_clear(U);
    fmpz_mat_clear(V);
    fmpz_mat_clear(H);
    fmpz_mat_clear(R);
}


void
tkf91_double_precision(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        named_generators_t g,
        slong *A, size_t szA,
        slong *B, size_t szB);

void
tkf91_double_precision(
        fmpz_mat_t mat, expr_ptr * expressions_table,
        named_generators_t g,
        slong *A, size_t szA,
        slong *B, size_t szB)
{
    slong level = 8;
    slong prec = 1 << level;

    arb_mat_t G;
    arb_mat_t expression_logs;
    arb_mat_t generator_logs;
    slong i;
    slong generator_count = fmpz_mat_nrows(mat);
    slong expression_count = fmpz_mat_ncols(mat);
    named_double_generators_t h;

    /* initialize the arbitrary precision exponent matrix */
    arb_mat_init(G, generator_count, expression_count);
    arb_mat_set_fmpz_mat(G, mat);

    /* compute the expression logs */
    arb_t x;
    arb_init(x);
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
    arb_clear(x);

    /* compute the generator logs */
    arb_mat_init(generator_logs, generator_count, 1);
    arb_mat_mul(generator_logs, G, expression_logs, prec);

    /* fill a structure with corresponding double precision values */
    doublify_named_generators(h, g, generator_logs);

    /* do the thing */
    tkf91_dynamic_programming(h, A, szA, B, szB);

    arb_mat_clear(G);
    arb_mat_clear(expression_logs);
    arb_mat_clear(generator_logs);
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
    tkf91_expressions_t p;
    generator_reg_t genreg;
    named_generators_t g;
    fmpz_mat_t mat;

    szA = strlen(strA);
    A = flint_malloc(szA * sizeof(slong));
    _fill_sequence_vector(A, strA, szA);

    szB = strlen(strB);
    B = flint_malloc(szB * sizeof(slong));
    _fill_sequence_vector(B, strB, szB);

    reg_init(reg);
    tkf91_expressions_init(p, reg,
            params->lambda, params->mu, params->tau, params->pi);
    generator_reg_init(genreg, reg->size);
    named_generators_init(g, genreg, p, A, szA, B, szB);

    fmpz_mat_init(mat,
            generator_reg_generators_len(genreg),
            generator_reg_expressions_len(genreg));
    generator_reg_get_matrix(mat, genreg);

    flint_printf("generator matrix:\n");
    fmpz_mat_print_pretty(mat);
    flint_printf("\n");

    /* report some transformation of the generator matrix */
    mess_with_generator_matrix(mat);

    expressions_table = reg_vec(reg);

    tkf91_double_precision(mat, expressions_table, g, A, szA, B, szB);

    reg_clear(reg);
    tkf91_expressions_clear(p);
    generator_reg_clear(genreg);
    named_generators_clear(g);
    fmpz_mat_clear(mat);

    flint_free(expressions_table);
    flint_free(A);
    flint_free(B);
}




int
main(int argc, char *argv[])
{
    int i;

    const char *Astr = NULL;
    const char *Bstr = NULL;

    user_params_t p;
    user_params_init(p);

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
        if (strcmp(argv[i], "--sequence-1") == 0) {
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

    run(Astr, Bstr, p);

    user_params_clear(p);

    flint_cleanup();
    return 0;
}
