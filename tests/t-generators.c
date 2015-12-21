#include <string.h>

#include "flint/flint.h"
#include "flint/fmpq.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

#include "femtocas.h"
#include "expressions.h"
#include "generators.h"
#include "tkf91_generator_indices.h"
#include "tkf91_rationals.h"
#include "tkf91_generators.h"


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


int main()
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("generators....");
    fflush(stdout);


    for (i = 0; i < 10; i++)
    {
        const char strA[] = "ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA";
        const char strB[] = "AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA";

        slong *A;
        slong *B;

        size_t szA, szB;

        szA = strlen(strA);
        A = (slong *) malloc(szA * sizeof(slong));
        _fill_sequence_vector(A, strA, szA);

        szB = strlen(strB);
        B = (slong *) malloc(szA * sizeof(slong));
        _fill_sequence_vector(B, strB, szB);


        fmpq_t lambda, mu, tau;
        fmpq pi[4];
        reg_t reg;
        tkf91_expressions_t p;
        generator_reg_t genreg;
        tkf91_generator_indices_t tkf91_generators;
        tkf91_rationals_t r;

        fmpq_init(lambda);
        fmpq_init(mu);
        fmpq_init(tau);

        /* lambda, mu, tau = 1, 1/2, 1/10 */
        fmpq_set_si(lambda, 1, 1);
        fmpq_set_si(mu, 1, 2);
        fmpq_set_si(tau, 1, 10);

        /* pi = (0.27, 0.24, 0.26, 0.23) */
        fmpq_init(pi+0); fmpq_set_si(pi+0, 27, 100);
        fmpq_init(pi+1); fmpq_set_si(pi+1, 24, 100);
        fmpq_init(pi+2); fmpq_set_si(pi+2, 26, 100);
        fmpq_init(pi+3); fmpq_set_si(pi+3, 23, 100);

        reg_init(reg);
        tkf91_rationals_init(r, lambda, mu, tau, pi);
        tkf91_expressions_init(p, reg, r);
        generator_reg_init(genreg, reg->size);
        tkf91_generators_init(tkf91_generators, genreg, p, A, szA, B, szB);


        /* report the matrix of coefficients */
        {
            fmpz_mat_t mat;
            fmpz_mat_init(mat,
                    generator_reg_generators_len(genreg),
                    generator_reg_expressions_len(genreg));
            generator_reg_get_matrix(mat, genreg);
            /* fmpz_mat_print_pretty(mat); */
            /* flint_printf("\n\n"); */
            fmpz_mat_clear(mat);
        }

        reg_clear(reg);
        tkf91_rationals_clear(r);
        tkf91_expressions_clear(p);
        generator_reg_clear(genreg);
        tkf91_generators_clear(tkf91_generators);

        fmpq_clear(lambda);
        fmpq_clear(mu);
        fmpq_clear(tau);
        fmpq_clear(pi+0);
        fmpq_clear(pi+1);
        fmpq_clear(pi+2);
        fmpq_clear(pi+3);

        free(A);
        free(B);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

