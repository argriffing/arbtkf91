/*
 * How to use Flint2 to find the basis of a lattice.
 */

#include <stdio.h>

#include "flint.h"
#include "fmpz_mat.h"

#include "arb.h"

void
check_inf_times_zero()
{
    arb_t x, y;
    fmpz_t z;
    slong prec = 200;

    arb_init(x);
    arb_neg_inf(x);

    fmpz_init(z);
    fmpz_zero(z);

    arb_init(y);
    arb_mul_fmpz(y, x, z, prec);

    flint_printf("x: "); arb_print(x); flint_printf("\n");
    flint_printf("y: "); arb_print(y); flint_printf("\n");

    fmpz_clear(z);
    arb_clear(x);
}


int
main(int argc, const char *argv[])
{
    fmpz_mat_t H, U, A, V, B, R;
    slong i, j;
    slong r = 6;
    slong c = 5;
    fmpz_t den;

    slong v[] = {1, 2, 3, 4, 5};

    fmpz_init(den);

    fmpz_mat_init(A, r, c);

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < c; j++)
        {
            fmpz_set_si(fmpz_mat_entry(A, i, j), 2 * i * j);
        }
    }
    for (j = 0; j < c; j++)
    {
            fmpz_set_si(fmpz_mat_entry(A, 4, j), v[j]);
            fmpz_set_si(fmpz_mat_entry(A, 5, j), 3 * v[j]);
    }

    fmpz_mat_init(B, c, r);
    fmpz_mat_transpose(B, A);



    //fmpz_mat_init(B, r, c);
    //fmpz_mat_set(B, A);

    fmpz_mat_init(U, fmpz_mat_nrows(B), fmpz_mat_nrows(B));
    fmpz_mat_init(V, fmpz_mat_nrows(B), fmpz_mat_nrows(B));
    fmpz_mat_init(H, fmpz_mat_nrows(B), fmpz_mat_ncols(B));

    /* hermite transform of B */
    fmpz_mat_hnf_transform(H, U, B);
    flint_printf("U * B = H\n");
    flint_printf("U:\n"); fmpz_mat_print_pretty(U);
    flint_printf("B:\n"); fmpz_mat_print_pretty(B);
    flint_printf("H:\n"); fmpz_mat_print_pretty(H);
    flint_printf("\n\n");

    /* inverse of U */
    fmpz_mat_inv(V, den, U);
    flint_printf("B = d*U^-1 * (1/d)H\n");
    flint_printf("B:\n"); fmpz_mat_print_pretty(B); flint_printf("\n");
    flint_printf("d: "); fmpz_print(den); flint_printf("\n");
    flint_printf("d * U^-1:\n"); fmpz_mat_print_pretty(V);
    flint_printf("H:\n"); fmpz_mat_print_pretty(H);
    flint_printf("\n\n");

    /* R */
    fmpz_mat_init(R, r, c);
    fmpz_mat_mul(R, A, V);
    flint_printf("A * (d*U^-1):\n");
    fmpz_mat_print_pretty(R);
    flint_printf("\n\n");

    /* smith transform */
    /*
    fmpz_mat_snf(H, A);
    flint_printf("S * A * T = H\n");
    flint_printf("S: ???\n");
    flint_printf("A:\n"); fmpz_mat_print_pretty(A);
    flint_printf("T: ???\n");
    flint_printf("H:\n"); fmpz_mat_print_pretty(H);
    flint_printf("\n\n");
    */

    return 0;
}
