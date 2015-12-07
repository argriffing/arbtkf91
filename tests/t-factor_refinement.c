#include <stdio.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_vec.h"
#include "factor_refinement.h"

void _fmpz_vec_randtest_pos(fmpz * f, flint_rand_t state,
        slong len, mp_bitcnt_t bits);

void
_fmpz_vec_randtest_pos(fmpz * f, flint_rand_t state,
        slong len, mp_bitcnt_t bits)
{
    slong i;
    _fmpz_vec_randtest_unsigned(f, state, len, bits-1);
    for (i = 0; i < len; i++)
    {
        fmpz_add_ui(f+i, f+i, 1);
    }
}

int main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("factor_refinement....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        fmpz *x, *ybase, *yexp;
        slong xlen, ylen;

        mp_bitcnt_t bits;

        x = ybase = yexp = NULL;
        xlen = ylen = 0;

        bits = n_randint(state, 100) + 2;
        xlen = n_randint(state, 100) + 1;

        /* create the random input vector of positive integers */
        x = _fmpz_vec_init(xlen);
        _fmpz_vec_randtest_pos(x, state, xlen, bits);

        factor_refinement(&ybase, &yexp, &ylen, x, xlen);

        /* check that products are equal */
        {
            fmpz_t a, b, p;
            slong j, u;

            fmpz_init_set_ui(a, 1);
            for (j = 0; j < xlen; j++)
            {
                fmpz_mul(a, a, x+j);
            }

            fmpz_init_set_ui(b, 1);
            fmpz_init(p);
            for (j = 0; j < ylen; j++)
            {
                u = fmpz_get_ui(yexp+j);
                fmpz_pow_ui(p, ybase+j, u);
                fmpz_mul(b, b, p);
            }

            if (!fmpz_equal(a, b))
            {
                flint_printf("FAIL:\n");
                fmpz_print(a); flint_printf(" ");
                fmpz_print(b); flint_printf("\n");
                abort();
            }

            fmpz_clear(a);
            fmpz_clear(b);
            fmpz_clear(p);
        }

        /* check that elements of the base are pairwise coprime */
        {
            slong u, v;
            fmpz_t g;
            fmpz_init(g);
            for (u = 0; u < ylen; u++)
            {
                for (v = 0; v < u; v++)
                {
                    fmpz_gcd(g, ybase+u, ybase+v);
                    if (!fmpz_is_one(g))
                    {
                        flint_printf("FAIL:\n");
                        fmpz_print(ybase+u); flint_printf(" ");
                        fmpz_print(ybase+v); flint_printf("\n");
                        abort();
                    }
                }
            }
            fmpz_clear(g);
        }

        /* check that each input is a product of powers of outputs */
        {
            slong u, v;
            fmpz_t a;
            fmpz_init(a);
            for (u = 0; u < xlen; u++)
            {
                fmpz_set(a, x+u);
                for (v = 0; v < ylen && !fmpz_is_one(a); v++)
                {
                    while (fmpz_divisible(a, ybase+v))
                    {
                        fmpz_divexact(a, a, ybase+v);
                    }
                }
                if (!fmpz_is_one(a))
                {
                    flint_printf("FAIL:\n");
                    fmpz_print(ybase+u); flint_printf("\n");
                    abort();
                }
            }
            fmpz_clear(a);
        }

        _fmpz_vec_clear(x, xlen);
        _fmpz_vec_clear(ybase, ylen);
        _fmpz_vec_clear(yexp, ylen);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
