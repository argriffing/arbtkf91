#ifndef WAVEFRONT_DOUBLE_H
#define WAVEFRONT_DOUBLE_H

#include "flint.h"

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
    slong modulus;
} wave_mat_struct;
typedef wave_mat_struct wave_mat_t[1];



#ifdef __cplusplus
extern "C" {
#endif

void wave_mat_init(wave_mat_t mat, slong n, slong modulus);
void wave_mat_clear(wave_mat_t mat);
wave_value_ptr wave_mat_entry(wave_mat_t mat, slong k, slong l);
wave_value_ptr wave_mat_entry_top(wave_mat_t mat, slong k, slong l);
wave_value_ptr wave_mat_entry_diag(wave_mat_t mat, slong k, slong l);
wave_value_ptr wave_mat_entry_left(wave_mat_t mat, slong k, slong l);


#ifdef __cplusplus
}
#endif

#endif
