#include "wavefront_hermite.h"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_vec.h"



void
hwave_mat_init(hwave_mat_t mat, slong n, slong modulus, slong rank)
{
    /* 
     * The rank is the number of fmpz that we need
     * for each of m0, m1, m2 in each wavefront cell.
     */

    /*
     * The modulus should be at least 3.
     * If you want the probability matrices, then the modulus
     * should be n.
     */
    if (modulus < 3)
    {
        flint_printf("the wavefront modulus must be at least 3 ");
        flint_printf("because updating diagonal k+2 depends on the ");
        flint_printf("values of diagonals k and k+1\n");
        abort();
    }
    if (modulus > n)
    {
        flint_printf("it makes no sense to use a modulus larger than ");
        flint_printf("the size of the matrix\n");
        abort();
    }

    /* create a huge brick of initialized fmpz */
    mat->number_of_fmpz = modulus * n * 3 * rank;
    mat->all_the_fmpz = _fmpz_vec_init(mat->number_of_fmpz);

    mat->modulus = modulus;
    mat->n = n;
    mat->data = malloc(sizeof(hwave_value_struct) * modulus * n);

    /* in each wavefront cell, point m0, m1, m2 to a piece of the fmpz brick */
    fmpz * base;
    hwave_value_ptr cell;
    slong k, l;
    for (k = 0; k < modulus; k++)
    {
        for (l = 0; l < n; l++)
        {
            cell = hwave_mat_entry(mat, k, l);
            base = mat->all_the_fmpz + k * (n * 3 * rank) + l * (3 * rank);
            cell->m0_vec = base + 0 * rank;
            cell->m1_vec = base + 1 * rank;
            cell->m2_vec = base + 2 * rank;
        }
    }
}

void
hwave_mat_clear(hwave_mat_t mat)
{
    _fmpz_vec_clear(mat->all_the_fmpz, mat->number_of_fmpz);
    free(mat->data);
}

hwave_value_ptr
hwave_mat_entry(hwave_mat_t mat, slong k, slong l)
{
    slong idx = (k % mat->modulus)*mat->n + l;
    if (idx < 0 || idx > mat->n * mat->n)
    {
        flint_printf("error indexing the wavefront matrix\n");
        abort();
    }
    return mat->data + idx;
}

hwave_value_ptr
hwave_mat_entry_top(hwave_mat_t mat, slong k, slong l)
{
    return hwave_mat_entry(mat, k-1, l+1);
}

hwave_value_ptr
hwave_mat_entry_diag(hwave_mat_t mat, slong k, slong l)
{
    return hwave_mat_entry(mat, k-2, l);
}

hwave_value_ptr
hwave_mat_entry_left(hwave_mat_t mat, slong k, slong l)
{
    return hwave_mat_entry(mat, k-1, l-1);
}
