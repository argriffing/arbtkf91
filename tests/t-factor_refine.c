#include <stdio.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_vec.h"
#include "factor_refine.h"

void fmpz_factor_randtest(fmpz_factor_t f, flint_rand_t state,
        slong num, mp_bitcnt_t bits);

void
fmpz_factor_randtest(fmpz_factor_t f, flint_rand_t state,
        slong num, mp_bitcnt_t bits)
{
    slong i;
    ulong n;
    int s;

    /* s is -1 or 1 or rarely 0 */
    s = 0;
    n = n_randint(state, 11);
    if (n)
    {
        s = (n % 2) ? 1 : -1;
    }

    _fmpz_factor_fit_length(f, num);
    for (i = 0; i < num; i++)
    {
        fmpz_randtest(f->p+i, state, bits);
        f->exp[i] = n_randint(state, 4);
    }

    f->num = num;
    f->sign = s;
}


int main(void)
{
    int iter, i;
    FLINT_TEST_INIT(state);

    flint_printf("factor_refine....");
    fflush(stdout);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpz_factor_t f, g;
        slong num;
        mp_bitcnt_t bits;

        bits = n_randint(state, 30) + 2;
        num = n_randint(state, 30);

        /* sample a factor structure that is not in canonical form */
        fmpz_factor_init(f);
        fmpz_factor_randtest(f, state, num, bits);

        /* compute the factor refinement */
        fmpz_factor_init(g);
        fmpz_factor_refine(g, f);

        /* each base must not be less than 2 */
        for (i = 0; i < g->num; i++)
        {
            if (fmpz_cmp_ui(g->p+i, 2) < 0)
            {
                flint_printf("FAIL: (base minimum)\n");
                fmpz_print(g->p+i); flint_printf("\n");
                abort();
            }
        }

        /* bases must be increasing */
        for (i = 0; i < g->num-1; i++)
        {
            if (fmpz_cmp(g->p+i, g->p+i+1) >= 0)
            {
                flint_printf("FAIL: (base sorting)\n");
                fmpz_print(g->p+i); flint_printf(" ");
                fmpz_print(g->p+i+1); flint_printf("\n");
                abort();
            }
        }

        /* each exponent must not be less than 1 */
        for (i = 0; i < g->num; i++)
        {
            if (g->exp[i] < 1)
            {
                flint_printf("FAIL: (exponent minimum)\n");
                flint_printf("%wd\n", g->exp[i]);
                abort();
            }
        }

        /* bases must be coprime */
        {
            slong u, v;
            fmpz_t x;
            fmpz_init(x);
            for (u = 0; u < g->num; u++)
            {
                for (v = 0; v < u; v++)
                {
                    fmpz_gcd(x, g->p+u, g->p+v);
                    if (!fmpz_is_one(x))
                    {
                        flint_printf("FAIL: (coprime bases)\n");
                        fmpz_print(g->p+u); flint_printf(" ");
                        fmpz_print(g->p+v); flint_printf("\n");
                        abort();
                    }
                }
            }
            fmpz_clear(x);
        }

        /* products must be equal when multiplied out */
        {
            fmpz_t x, y;

            fmpz_init(x);
            fmpz_init(y);

            fmpz_factor_expand(x, f);
            fmpz_factor_expand(y, g);

            if (!fmpz_equal(x, y))
            {
                flint_printf("FAIL: (products)\n");
                fmpz_factor_print(f); flint_printf(" : ");
                fmpz_print(x); flint_printf("\n");
                fmpz_factor_print(g); flint_printf(" : ");
                fmpz_print(y); flint_printf("\n");
                abort();
            }

            fmpz_clear(x);
            fmpz_clear(y);
        }

        fmpz_factor_clear(f);
        fmpz_factor_clear(g);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
