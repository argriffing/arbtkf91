#ifndef WAVEFRONT_HERMITE_H
#define WAVEFRONT_HERMITE_H

#include "flint/flint.h"
#include "flint/fmpz.h"


typedef struct
{
    fmpz * m0_vec;
    fmpz * m1_vec;
    fmpz * m2_vec;
} hwave_value_struct;
typedef hwave_value_struct * hwave_value_ptr;


typedef struct
{
    hwave_value_ptr data;
    slong n;
    slong modulus;
    fmpz * all_the_fmpz;
    slong number_of_fmpz;
} hwave_mat_struct;
typedef hwave_mat_struct hwave_mat_t[1];



#ifdef __cplusplus
extern "C" {
#endif

void hwave_mat_init(hwave_mat_t mat, slong n, slong modulus, slong rank);
void hwave_mat_clear(hwave_mat_t mat);
hwave_value_ptr hwave_mat_entry(hwave_mat_t mat, slong k, slong l);
hwave_value_ptr hwave_mat_entry_top(hwave_mat_t mat, slong k, slong l);
hwave_value_ptr hwave_mat_entry_diag(hwave_mat_t mat, slong k, slong l);
hwave_value_ptr hwave_mat_entry_left(hwave_mat_t mat, slong k, slong l);

#ifdef __cplusplus
}
#endif

#endif
