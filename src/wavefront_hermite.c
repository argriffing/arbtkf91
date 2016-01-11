#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_vec.h"

#include "arb.h"

#include "wavefront_hermite.h"


void
hwave_element_init(hwave_element_t p)
{
    arb_init(p->value);
}

void
hwave_element_clear(hwave_element_t p)
{
    arb_clear(p->value);
}

void
hwave_element_set_undefined(hwave_element_t p)
{
    arb_neg_inf(p->value);
    p->status = HWAVE_STATUS_UNDEFINED;
}

void
hwave_mat_init(hwave_mat_t mat, slong n, slong modulus, slong rank)
{
    /* 
     * The rank is the number of fmpz that we need
     * for each of m[0], m[1], m[2] in each wavefront cell.
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
    mat->data = malloc(sizeof(hwave_cell_struct) * modulus * n);

    /*
     * In each wavefront cell, point the coefficient vectors
     * of each of m[0], m[1], m[2] to a piece of the fmpz brick.
     * Also initialize the arb_t value in each element in each cell.
     */
    fmpz * base;
    hwave_cell_ptr cell;
    slong k, l, c;
    for (k = 0; k < modulus; k++)
    {
        for (l = 0; l < n; l++)
        {
            cell = hwave_mat_entry(mat, k, l);
            base = mat->all_the_fmpz + k * (n * 3 * rank) + l * (3 * rank);
            for (c = 0; c < 3; c++)
            {
                cell->m[c].vec = base + c * rank;
                hwave_element_init(cell->m + c);
            }
        }
    }
}

void
hwave_mat_clear(hwave_mat_t mat)
{
    /* Clear the arb_t value in each element in each cell. */
    slong k, l, c;
    hwave_cell_ptr cell;
    for (k = 0; k < mat->modulus; k++)
    {
        for (l = 0; l < mat->n; l++)
        {
            cell = hwave_mat_entry(mat, k, l);
            for (c = 0; c < 3; c++)
            {
                hwave_element_clear(cell->m + c);
            }
        }
    }

    /* Clear the huge matrix of integer coefficients. */
    _fmpz_vec_clear(mat->all_the_fmpz, mat->number_of_fmpz);

    /* Free the array of cells. */
    free(mat->data);
}

hwave_cell_ptr
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

hwave_cell_ptr
hwave_mat_entry_top(hwave_mat_t mat, slong k, slong l)
{
    return hwave_mat_entry(mat, k-1, l+1);
}

hwave_cell_ptr
hwave_mat_entry_diag(hwave_mat_t mat, slong k, slong l)
{
    return hwave_mat_entry(mat, k-2, l);
}

hwave_cell_ptr
hwave_mat_entry_left(hwave_mat_t mat, slong k, slong l)
{
    return hwave_mat_entry(mat, k-1, l-1);
}
