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

    /* test a more complicated expression */
    {
        fmpq_t a, b, t;
        fmpq_t at, bt;

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(t);

        fmpq_set_si(a, 1, 2);
        fmpq_set_si(b, 3, 4);
        fmpq_set_si(t, 5, 6);

        fmpq_init(at);
        fmpq_init(bt);

        fmpq_mul(at, a, t);
        fmpq_mul(bt, b, t);

        expr_t ax, bx;
        expr_t eat, ebt, aeat, bebt;
        expr_t num, den;
        expr_t res;

        expr_fmpq(ax, a);
        expr_fmpq(bx, b);
        expr_exp_fmpq(ebt, bt);
        expr_exp_fmpq(eat, at);
        expr_mul(bebt, bx, ebt);
        expr_mul(aeat, ax, eat);
        expr_sub(num, ebt, eat);
        expr_sub(den, bebt, aeat);
        expr_div(res, num, den);

        arb_t value;
        slong level;

        arb_init(value);

        for (level = 0; level < 8; level++)
        {
            /* flint_printf("evaluating "); */
            /* expr_print(res); */
            /* flint_printf(" at level %wd:\n", level); */
            expr_eval(value, res, level);
            /* arb_print(value); */
            /* flint_printf("\n\n"); */
        }

        arb_clear(value);
        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(t);
        expr_clear(ax);
        expr_clear(bx);
        expr_clear(eat);
        expr_clear(ebt);
        expr_clear(aeat);
        expr_clear(bebt);
        expr_clear(num);
        expr_clear(den);
        expr_clear(res);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
