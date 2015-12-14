#include <stdio.h>
#include "flint/flint.h"
#include "flint/fmpq.h"
#include "arb.h"
#include "femtocas.h"

int main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("femtocas....");
    fflush(stdout);

    for (i = 0; i < 10; i++)
    {
        slong level;

        fmpq_t a, b;
        expr_t x, y;
        fmpq_init(a);
        fmpq_init(b);
        fmpq_set_si(a, -2, 3);
        fmpq_set_si(b, 5, 7);
        expr_fmpq(x, a);
        expr_fmpq(y, b);

        /* test addition */
        {
            fmpq_t desired;
            fmpq_init(desired);
            fmpq_add(desired, a, b);

            expr_t s;
            expr_add(s, x, y);

            for (level = 0; level < 10; level++)
            {
                arb_t actual;
                arb_init(actual);

                expr_eval(actual, s, level);
                if (!arb_contains_fmpq(actual, desired))
                {
                    flint_printf("FAIL: (add containment)\n");
                    abort();
                }

                arb_clear(actual);
            }

            fmpq_clear(desired);
            expr_clear(s);
        }

        /* test multiplication */
        {
            fmpq_t desired;
            fmpq_init(desired);
            fmpq_mul(desired, a, b);

            expr_t p;
            expr_mul(p, x, y);

            for (level = 0; level < 10; level++)
            {
                arb_t actual;
                arb_init(actual);

                expr_eval(actual, p, level);
                if (!arb_contains_fmpq(actual, desired))
                {
                    flint_printf("FAIL: (mul containment)\n");
                    abort();
                }

                arb_clear(actual);
            }

            fmpq_clear(desired);
            expr_clear(p);
        }

        fmpq_clear(a);
        fmpq_clear(b);
        expr_clear(x);
        expr_clear(y);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
