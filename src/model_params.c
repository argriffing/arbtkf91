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




static int _assert_fmpq_positive(const fmpq_t x, const char *name);

int
_assert_fmpq_positive(const fmpq_t x, const char *name)
{
    if (fmpq_sgn(x) != 1 || !fmpq_is_canonical(x))
    {
        flint_fprintf(stderr, name);
        flint_fprintf(stderr, " must be positive and well defined\n");
        return -1;
    }
    return 0;
}


int
model_params_validate(const model_params_t p)
{
    int i, res, failflag;
    fmpq_t x;

    /* require all parameters be positive and well defined */
    res = _assert_fmpq_positive(p->lambda, "lambda"); if (res) return res;
    res = _assert_fmpq_positive(p->mu, "mu"); if (res) return res;
    res = _assert_fmpq_positive(p->tau, "tau"); if (res) return res;
    res = _assert_fmpq_positive(p->pi+0, "pa"); if (res) return res;
    res = _assert_fmpq_positive(p->pi+1, "pc"); if (res) return res;
    res = _assert_fmpq_positive(p->pi+2, "pg"); if (res) return res;
    res = _assert_fmpq_positive(p->pi+3, "pt"); if (res) return res;

    /* require lambda < mu */
    if (fmpq_cmp(p->lambda, p->mu) != -1)
    {
        flint_fprintf(stderr, "lambda must be less than mu\n");
        return -1;
    }

    /* require the probabilities to sum to 1 */
    fmpq_init(x);
    fmpq_zero(x);
    for (i = 0; i < 4; i++)
    {
        fmpq_add(x, x, p->pi+i);
    }
    failflag = !fmpq_is_one(x);
    fmpq_clear(x);
    if (failflag)
    {
        flint_fprintf(stderr, "probabilities must sum to 1\n");
        return -1;
    }

    return 0;
}
