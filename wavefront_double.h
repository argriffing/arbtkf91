#ifndef WAVEFRONT_DOUBLE_H
#define WAVEFRONT_DOUBLE_H

#include "flint/flint.h"


/* double precision implementation of wavefront indexing */

/*
 * This placeholder uses only a double precision log likelihood value
 * for each entry in the table.
 */
typedef struct
{
    double m0;
    double m1;
    double m2;
} wave_value_struct;
typedef wave_value_struct * wave_value_ptr;


/* Assume that the access pattern needs only three physical rows. */
typedef struct
{
    wave_value_ptr data;
    slong n;
    /* slong modulus; */
} wave_mat_struct;
typedef wave_mat_struct wave_mat_t[1];



#ifdef __cplusplus
extern "C" {
#endif

/* void wave_mat_init(wave_mat_t mat, slong n, slong modulus); */
void wave_mat_init(wave_mat_t mat, slong n);
void wave_mat_clear(wave_mat_t mat);

static __inline__
wave_value_ptr wave_mat_entry(wave_mat_t mat, slong k, slong l)
{
    /*
    slong idx = (k % mat->modulus)*mat->n + l;
    if (idx < 0 || idx > mat->n * mat->n)
    {
        flint_printf("error indexing the wavefront matrix\n");
        abort();
    }
    return mat->data + idx;
    */
    return mat->data + (k % 3)*mat->n + l;
}

static __inline__
wave_value_ptr wave_mat_entry_top(wave_mat_t mat, slong k, slong l)
{
    return wave_mat_entry(mat, k-1, l+1);
}

static __inline__
wave_value_ptr wave_mat_entry_diag(wave_mat_t mat, slong k, slong l)
{
    return wave_mat_entry(mat, k-2, l);
}

static __inline__
wave_value_ptr wave_mat_entry_left(wave_mat_t mat, slong k, slong l)
{
    return wave_mat_entry(mat, k-1, l-1);
}



#ifdef __cplusplus
}
#endif

#endif
