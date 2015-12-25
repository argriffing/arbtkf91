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
#include "tkf91_dp_d.h"
#include "tkf91_dp_r.h"
#include "tkf91_dp_bound.h"


#define MAXSEQLEN 20000


typedef struct
{
    int trace_flag;
    fmpq_t lambda;
    fmpq_t mu;
    fmpq_t tau;
    fmpq pi[4];
} user_params_struct;
typedef user_params_struct user_params_t[1];

void user_params_init(user_params_t p);
void user_params_clear(user_params_t p);
void user_params_print(const user_params_t p);


/* function pointer typedef for the dynamic programming function */
typedef void (*tkf91_dp_fn)(
        fmpz_mat_t, expr_ptr *,
        tkf91_generator_indices_t,
        int,
        slong *, size_t,
        slong *, size_t);

void _fill_sequence_vector(slong *v, const char *str, slong n);
void bench(tkf91_dp_fn f, const user_params_t params);
void run(tkf91_dp_fn f, const user_params_t params,
        const char *strA, const char *strB);
void _run(tkf91_dp_fn f, const user_params_t p,
        slong *A, size_t szA, slong *B, size_t szB);


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


void
user_params_init(user_params_t p)
{
    slong i;
    p->trace_flag = 0;
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



void
_run(tkf91_dp_fn f, const user_params_t p,
        slong *A, size_t szA, slong *B, size_t szB)
{
    tkf91_rationals_t r;
    tkf91_expressions_t expressions;
    tkf91_generator_indices_t generators;
    fmpz_mat_t mat;
    expr_ptr * expressions_table;
    clock_t diff, start;
    int msec;

    /* timing */
    start = clock();

    /* expressions registry and (refining) generator registry */
    reg_t er;
    rgen_reg_ptr gr;

    reg_init(er);
    tkf91_rationals_init(r, p->lambda, p->mu, p->tau, p->pi);
    tkf91_expressions_init(expressions, er, r);

    gr = rgen_reg_new();
    tkf91_rgenerators_init(generators, gr, r, expressions, A, szA, B, szB);
    rgen_reg_finalize(gr, er);
    fmpz_mat_init(mat, rgen_reg_nrows(gr), rgen_reg_ncols(gr));
    rgen_reg_get_matrix(mat, gr);

    rgen_reg_clear(gr);
    tkf91_rationals_clear(r);

    expressions_table = reg_vec(er);

    f(mat, expressions_table, generators, p->trace_flag, A, szA, B, szB);

    /* timing */
    diff = clock() - start;
    msec = (diff * 1000) / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds.\n", msec/1000, msec%1000);

    fmpz_mat_clear(mat);
    flint_free(expressions_table);

    reg_clear(er);
    tkf91_expressions_clear(expressions);
}



void
run(tkf91_dp_fn f, const user_params_t params,
        const char *strA, const char *strB)
{
    slong *A;
    slong *B;
    size_t szA, szB;

    szA = strlen(strA);
    A = flint_malloc(szA * sizeof(slong));
    _fill_sequence_vector(A, strA, szA);

    szB = strlen(strB);
    B = flint_malloc(szB * sizeof(slong));
    _fill_sequence_vector(B, strB, szB);

    _run(f, params, A, szA, B, szB);

    flint_free(A);
    flint_free(B);
}


void
bench(tkf91_dp_fn f, const user_params_t params)
{
    char strA[MAXSEQLEN];
    char strB[MAXSEQLEN];
    size_t szA, szB;

    slong *A;
    slong *B;

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

        _run(f, params, A, szA, B, szB);

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

    slong lambda_num = 0;
    slong mu_num = 0;
    slong tau_num = 0;
    slong pi_num[4] = {0, 0, 0, 0};

    slong lambda_den = 1;
    slong mu_den = 1;
    slong tau_den = 1;
    slong pi_den[4] = {1, 1, 1, 1};

    tkf91_dp_fn f = NULL;

    for (i = 1; i < argc-1; i += 2)
    {
        if (strcmp(argv[i], "--bench") == 0) {
            flint_sscanf(argv[i + 1], "%d", &bench_flag);
        } else if (strcmp(argv[i], "--trace") == 0) {
            flint_sscanf(argv[i + 1], "%d", &(p->trace_flag));
        } else if (strcmp(argv[i], "--precision") == 0) {
            if (strcmp(argv[i + 1], "double") == 0) {
                f = &tkf91_dp_d;
            } else if (strcmp(argv[i + 1], "arbitrary") == 0) {
                f = &tkf91_dp_r;
            } else if (strcmp(argv[i + 1], "bound") == 0) {
                f = &tkf91_dp_bound;
            }
        }

        else if (strcmp(argv[i], "--sequence-1") == 0) {
            Astr = argv[i + 1];
        } else if (strcmp(argv[i], "--sequence-2") == 0) {
            Bstr = argv[i + 1];
        }
        
        else if (strcmp(argv[i], "--lambda-num") == 0) {
            flint_sscanf(argv[i + 1], "%wd", &lambda_num);
        } else if (strcmp(argv[i], "--mu-num") == 0) {
            flint_sscanf(argv[i + 1], "%wd", &mu_num);
        } else if (strcmp(argv[i], "--tau-num") == 0) {
            flint_sscanf(argv[i + 1], "%wd", &tau_num);
        } else if (strcmp(argv[i], "--pa-num") == 0) {
            flint_sscanf(argv[i + 1], "%wd", pi_num+0);
        } else if (strcmp(argv[i], "--pc-num") == 0) {
            flint_sscanf(argv[i + 1], "%wd", pi_num+1);
        } else if (strcmp(argv[i], "--pg-num") == 0) {
            flint_sscanf(argv[i + 1], "%wd", pi_num+2);
        } else if (strcmp(argv[i], "--pt-num") == 0) {
            flint_sscanf(argv[i + 1], "%wd", pi_num+3);
        }

        else if (strcmp(argv[i], "--lambda-den") == 0) {
            flint_sscanf(argv[i + 1], "%wd", &lambda_den);
        } else if (strcmp(argv[i], "--mu-den") == 0) {
            flint_sscanf(argv[i + 1], "%wd", &mu_den);
        } else if (strcmp(argv[i], "--tau-den") == 0) {
            flint_sscanf(argv[i + 1], "%wd", &tau_den);
        } else if (strcmp(argv[i], "--pa-den") == 0) {
            flint_sscanf(argv[i + 1], "%wd", pi_den+0);
        } else if (strcmp(argv[i], "--pc-den") == 0) {
            flint_sscanf(argv[i + 1], "%wd", pi_den+1);
        } else if (strcmp(argv[i], "--pg-den") == 0) {
            flint_sscanf(argv[i + 1], "%wd", pi_den+2);
        } else if (strcmp(argv[i], "--pt-den") == 0) {
            flint_sscanf(argv[i + 1], "%wd", pi_den+3);
        }
    }

    if (!f)
    {
        flint_printf("expected '--precision {double, arbitrary, bound}");
    }
    else
    {
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
            bench(f, p);
        }
        else
        {
            run(f, p, Astr, Bstr);
        }
    }

    user_params_clear(p);
    flint_cleanup();
    return 0;
}
