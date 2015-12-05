/*
 * -I/usr/local/include/flint -lflint
 */

/*
--sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA
--sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA
--lambda 1 --mu 2 --tau 0.1 --pa 0.25 --pc 0.25 --pg 0.25 --pt 0.25
*/

/* this is kind of ridiculous... */
/* but parsing sucks so I want to be super explicit */
/*
--lambda-num 1 --lambda-den 1 --mu-num 2 --mu-den 1 --tau-num 1 --tau-den 10 --pa-num 1 --pa-den 4 --pc-num 1 --pc-den 4 --pg-num 1 --pg-den 4 --pt-num 1 --pt-den 4
*/


#include <stdio.h>
#include <stdlib.h>

#include "flint.h"
#include "fmpq.h"


typedef struct
{
    fmpq lambda[1];
    fmpq mu[1];
    fmpq tau[1];
    fmpq pi[4];
}
user_params_struct;

typedef user_params_struct user_params_t[1];

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
    slong i;
    const char dna[] = "";
    flint_printf("lambda: "); fmpq_print(p->lambda); flint_printf("\n");
    flint_printf("mu: "); fmpq_print(p->mu); flint_printf("\n");
    flint_printf("tau: "); fmpq_print(p->tau); flint_printf("\n");
    flint_printf("pa: "); fmpq_print(p->pi+0); flint_printf("\n");
    flint_printf("pc: "); fmpq_print(p->pi+1); flint_printf("\n");
    flint_printf("pg: "); fmpq_print(p->pi+2); flint_printf("\n");
    flint_printf("pt: "); fmpq_print(p->pi+3); flint_printf("\n");
}

int
main(int argc, char *argv[])
{
    int i;

    const char *s_A, *s_B;

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
            s_A = argv[i + 1];
        } else if (strcmp(argv[i], "--sequence-2") == 0) {
            s_B = argv[i + 1];
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

    user_params_print(p);

    user_params_clear(p);

    return 0;
}
