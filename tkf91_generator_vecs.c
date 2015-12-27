#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

#include "tkf91_generator_indices.h"
#include "tkf91_generator_vecs.h"

fmpz * _vecify(slong i, const fmpz_mat_t M);

fmpz *
_vecify(slong i, const fmpz_mat_t M)
{
    return M->rows[i];
}

void
tkf91_generator_vecs_init(
        tkf91_generator_vecs_t h,
        const tkf91_generator_indices_t g,
        const fmpz_mat_t V,
        slong rank)
{
    slong i, j;

    /* initialize `h->M` to the first `rank` columns of `V` */
    fmpz_mat_struct * M = h->M;
    fmpz_mat_init(M, fmpz_mat_nrows(V), rank);
    for (i = 0; i < fmpz_mat_nrows(V); i++)
    {
        for (j = 0; j < rank; j++)
        {
            fmpz_set(fmpz_mat_entry(M, i, j), fmpz_mat_entry(V, i, j));
        }
    }

    h->m1_00 = _vecify(g->m1_00, M);
    h->m0_10 = _vecify(g->m0_10, M);
    h->m2_01 = _vecify(g->m2_01, M);
    for (i = 0; i < 4; i++)
    {
        h->m0_i0_incr[i] = _vecify(g->m0_i0_incr[i], M);
        h->m2_0j_incr[i] = _vecify(g->m2_0j_incr[i], M);
        h->c0_incr[i] = _vecify(g->c0_incr[i], M);
        for (j = 0; j < 4; j++)
        {
            h->c1_incr[i*4+j] = _vecify(g->c1_incr[i*4+j], M);
        }
        h->c2_incr[i] = _vecify(g->c2_incr[i], M);
    }
}


void tkf91_generator_vecs_clear(tkf91_generator_vecs_t h)
{
    fmpz_mat_clear(h->M);
}


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
