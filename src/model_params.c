#include <stdio.h>

#include "model_params.h"

void
model_params_init(model_params_t p)
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
model_params_clear(model_params_t p)
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


/* TODO this needs fmpq_fprint support...
void
_fmpq_fprint_named(FILE * file, const fmpq_t p, const char *name)
{
    flint_fprintf(file, "%s: ", name);
}
*/


void
model_params_print(const model_params_t p)
{
    flint_printf("lambda: "); fmpq_print(p->lambda); flint_printf("\n");
    flint_printf("mu: "); fmpq_print(p->mu); flint_printf("\n");
    flint_printf("tau: "); fmpq_print(p->tau); flint_printf("\n");
    flint_printf("pa: "); fmpq_print(p->pi+0); flint_printf("\n");
    flint_printf("pc: "); fmpq_print(p->pi+1); flint_printf("\n");
    flint_printf("pg: "); fmpq_print(p->pi+2); flint_printf("\n");
    flint_printf("pt: "); fmpq_print(p->pi+3); flint_printf("\n");
}
