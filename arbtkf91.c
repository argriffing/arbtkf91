#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "flint/flint.h"
#include "flint/fmpq.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "femtocas.h"
#include "expressions.h"
#include "generators.h"


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



int run(const char *strA, const char *strB, const user_params_t params);

int
run(const char *strA, const char *strB, const user_params_t params)
{
    slong *A;
    slong *B;

    size_t szA, szB;

    szA = strlen(strA);
    A = (slong *) malloc(szA * sizeof(slong));
    _fill_sequence_vector(A, strA, szA);

    szB = strlen(strB);
    B = (slong *) malloc(szA * sizeof(slong));
    _fill_sequence_vector(B, strB, szB);

    reg_t reg;
    tkf91_expressions_t p;
    generator_reg_t genreg;
    named_generators_t named_generators;

    reg_init(reg);
    tkf91_expressions_init(p, reg,
            params->lambda, params->mu, params->tau, params->pi);
    generator_reg_init(genreg, reg->size);
    named_generators_init(named_generators, genreg, p, A, szA, B, szB);

    /* report the matrix of coefficients */
    {
        fmpz_mat_t mat;
        fmpz_mat_init(mat,
                generator_reg_generators_len(genreg),
                generator_reg_expressions_len(genreg));
        generator_reg_get_matrix(mat, genreg);
        fmpz_mat_print_pretty(mat);
        fmpz_mat_clear(mat);
        flint_printf("\n\n");
    }

    reg_clear(reg);
    tkf91_expressions_clear(p);
    generator_reg_clear(genreg);
    named_generators_clear(named_generators);

    free(A);
    free(B);

    return 0;
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

    int result = run(Astr, Bstr, p);

    user_params_clear(p);

    return result;
}
