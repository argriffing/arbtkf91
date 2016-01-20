#include "flint/flint.h"
#include "flint/fmpz_mat.h"


void
_fmpz_mat_hnf_inverse_transform(
        fmpz_mat_t H, fmpz_mat_t V, slong * prank, const fmpz_mat_t A)
{
    /* U*A = H ; U^-1 = V */
    fmpz_mat_t U;
    fmpz_t den;
    slong i, rank;
    int result;

    if (fmpz_mat_nrows(A) != fmpz_mat_nrows(H) ||
        fmpz_mat_ncols(A) != fmpz_mat_ncols(H))
    {
        flint_printf("the matrix 'H' should have the same dimensions ");
        flint_printf("as the matrix 'A'\n");
        abort();
    }
    if (fmpz_mat_nrows(A) != fmpz_mat_nrows(V) ||
        fmpz_mat_nrows(A) != fmpz_mat_ncols(V))
    {
        flint_printf("the matrix 'V' should be square with dimensions (n, n) ");
        flint_printf("where n is the number of rows in the matrix 'A'\n");
        abort();
    }

    fmpz_mat_init(U, fmpz_mat_nrows(A), fmpz_mat_nrows(A));
    fmpz_init(den);

    fmpz_mat_hnf_transform(H, U, A);

    /*
     * The number of nonzero rows of H should determine the rank of mat.
     * In H, rows with nonzero entries should precede rows without 
     * nonzero entries
     */
    rank = 0;
    for (i = 0; i < fmpz_mat_nrows(H); i++)
    {
        if (!fmpz_mat_is_zero_row(H, i))
        {
            if (i != rank)
            {
                flint_printf("expected each row in H ");
                flint_printf("containing a nonzero entry ");
                flint_printf("to precede each row containing only zeros\n");
                abort();
            }
            rank++;
        }
    }

    /*
     * Invert the transform matrix U.
     * This should be possible because U should be unimodular
     * and therefore nonsingular.
     * The inverse should have entries that are integers,
     * and its determinant should be +1 or -1.
     *
     * We only care about the first r columns of this inverse,
     * where r is the rank of the Hermite form H.
     */
    result = fmpz_mat_inv(V, den, U);
    if (!result)
    {
        flint_printf("expected U to be nonsingular\n");
        abort();
    }
    if (!fmpz_is_pm1(den))
    {
        flint_printf("expected U to be unimodular -- ");
        flint_printf("denominator of inverse of U: ");
        fmpz_print(den);
        flint_printf("\n");
        abort();
    }
    if (!fmpz_is_one(den))
    {
        fmpz_mat_neg(V, V);
    }

    *prank = rank;

    fmpz_clear(den);
    fmpz_mat_clear(U);
}


int
main()
{
    slong rank = 0;
    slong nrows = 3;
    slong ncols = 3;

    fmpz_mat_t A, H, V;
    fmpz_mat_init(A, nrows, ncols);
    fmpz_mat_init(H, nrows, ncols);
    fmpz_mat_init(V, nrows, nrows); /* this is a square matrix */

    fmpz_set_si(fmpz_mat_entry(A, 0, 0), 1);
    fmpz_set_si(fmpz_mat_entry(A, 0, 1), -1);
    fmpz_set_si(fmpz_mat_entry(A, 0, 2), 1);

    fmpz_set_si(fmpz_mat_entry(A, 1, 0), 2);
    fmpz_set_si(fmpz_mat_entry(A, 1, 1), 0);
    fmpz_set_si(fmpz_mat_entry(A, 1, 2), -1);

    fmpz_set_si(fmpz_mat_entry(A, 2, 0), 3);
    fmpz_set_si(fmpz_mat_entry(A, 2, 1), -1);
    fmpz_set_si(fmpz_mat_entry(A, 2, 2), 0);

    _fmpz_mat_hnf_inverse_transform(H, V, &rank, A);

    flint_printf("rank: %wd\n", rank);
    fmpz_mat_print_pretty(A); flint_printf("\n");
    fmpz_mat_print_pretty(H); flint_printf("\n");
    fmpz_mat_print_pretty(V); flint_printf("\n");

    fmpz_mat_clear(A);
    fmpz_mat_clear(H);
    fmpz_mat_clear(V);

    return 0;
}
