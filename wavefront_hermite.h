#ifndef WAVEFRONT_HERMITE_H
#define WAVEFRONT_HERMITE_H

#include "flint/flint.h"
#include "flint/fmpz.h"

#include "arb.h"


/*
 * The top level of the data layout consists of a matrix of cells.
 * At the next deeper layer, each cell consists of three elements.
 * Each of the three elements has a status and a value and a pointer
 * to a vector of integer coefficients in some basis.
 */

/*
 * These statuses are mutually exclusive like an enum;
 * they are not binary flags.
 * They refer to the uninitialization, initialization, and loss
 * of the coefficient array associated to each cell element.
 * Loss of this array could occur due to numerical ambiguity.
 */
#define HWAVE_STATUS_UNDEFINED 0
#define HWAVE_STATUS_UNAMBIGUOUS 1
#define HWAVE_STATUS_AMBIGUOUS 2

typedef struct
{
    fmpz * vec;
    int status;
    arb_t value;
} hwave_element_struct;
typedef hwave_element_struct hwave_element_t[1];
typedef hwave_element_struct * hwave_element_ptr;

typedef struct
{
    hwave_element_struct m[3];
} hwave_cell_struct;
typedef hwave_cell_struct * hwave_cell_ptr;

typedef struct
{
    hwave_cell_ptr data;
    slong n;
    slong modulus;
    fmpz * all_the_fmpz;
    slong number_of_fmpz;
} hwave_mat_struct;
typedef hwave_mat_struct hwave_mat_t[1];



#ifdef __cplusplus
extern "C" {
#endif

void hwave_element_init(hwave_element_t);
void hwave_element_clear(hwave_element_t);
void hwave_element_set_undefined(hwave_element_t);

void hwave_mat_init(hwave_mat_t mat, slong n, slong modulus, slong rank);
void hwave_mat_clear(hwave_mat_t mat);
hwave_cell_ptr hwave_mat_entry(hwave_mat_t mat, slong k, slong l);
hwave_cell_ptr hwave_mat_entry_top(hwave_mat_t mat, slong k, slong l);
hwave_cell_ptr hwave_mat_entry_diag(hwave_mat_t mat, slong k, slong l);
hwave_cell_ptr hwave_mat_entry_left(hwave_mat_t mat, slong k, slong l);

#ifdef __cplusplus
}
#endif

#endif
