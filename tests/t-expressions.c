#include "flint/flint.h"
#include "flint/fmpq.h"
#include "expressions.h"

int main()
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("expressions....");
    fflush(stdout);

    for (i = 0; i < 10; i++)
    {
        fmpq_t lambda, mu, tau;
        fmpq pi[4];
        reg_t reg;
        tkf91_expressions_t p;

        fmpq_init(lambda);
        fmpq_init(mu);
        fmpq_init(tau);

        /* lambda, mu, t = 1, 1/2, 1/10 */
        fmpq_set_si(lambda, 1, 1);
        fmpq_set_si(mu, 1, 2);
        fmpq_set_si(tau, 1, 10);

        /* pi = (0.27, 0.24, 0.26, 0.23) */
        fmpq_init(pi+0); fmpq_set_si(pi+0, 27, 100);
        fmpq_init(pi+1); fmpq_set_si(pi+1, 24, 100);
        fmpq_init(pi+2); fmpq_set_si(pi+2, 26, 100);
        fmpq_init(pi+3); fmpq_set_si(pi+3, 23, 100);

        reg_init(reg);
        tkf91_expressions_init(p, reg, lambda, mu, tau, pi);

        reg_clear(reg);
        tkf91_expressions_clear(p);
        fmpq_clear(lambda);
        fmpq_clear(mu);
        fmpq_clear(tau);
        fmpq_clear(pi+0);
        fmpq_clear(pi+1);
        fmpq_clear(pi+2);
        fmpq_clear(pi+3);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
