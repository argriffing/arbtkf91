#include <stdlib.h>

#include "flint/flint.h"

#include "wavefront_double.h"


void
/* wave_mat_init(wave_mat_t mat, slong n, slong modulus) */
wave_mat_init(wave_mat_t mat, slong n)
{
    /*
     * The modulus should be at least 3.
     * If you want the probability matrices, then the modulus
     * should be n.
     */
    int modulus = 3;
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
    /* mat->modulus = modulus; */
    mat->n = n;
    mat->data = malloc(sizeof(wave_value_struct) * modulus * n);
}

void
wave_mat_clear(wave_mat_t mat)
{
    free(mat->data);
}
