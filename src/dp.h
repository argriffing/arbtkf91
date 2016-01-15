#ifndef DP_H
#define DP_H

/*
 * This is meant to be a more sophisticated replacement of 'breadcrumbs'.
 *
 * The DP_* flags answer the following questions for each cell:
 * 1) Is max3 interesting?
 *    If so,
 * 1a) Is max3 a candidate to be on the optimal traceback path?
 * 1b) Which of {m0, m1, m2} are in the running for max3?
 * 2) Is max2 interesting?
 *    If so,
 * 2a) Which of {m1, m2} are in the running for max2?
 *
 * Other relevant per-cell questions that are deliberately not addressed
 * by these flags are:
 * 3) Which of {m0, m1, m2} are interesting? This can be deduced from
 *    the above flags.
 * 3a) For each interesting element, what is its associated data?
 *     (this data could consist of an arb_t value, or a (low, high)
 *      mag_t pair, or an integer vector of coefficients in the case
 *      of checking that numerical ties are also symbolic ties).
 *     This associated data is required only for the 'forward' pass,
 *     and only a couple 'wavefront' rows of the data needs to be
 *     active at a time, so this can be deferred to the forward algorithm
 *     implementation.
 *
 * In this notation a quantity is said to be 'interesting' if it needs
 * to be evaluated, as determined using forward and backward
 * dynamic programming passes and interval arithmetic.
 * The notation MAX3 represents the maximum value of {m0, m1, m2},
 * and the notation MAX2 represents the maximum value of {m1, m2},
 * where m0, m1, and m2 are three values tracked per cell of the dynamic
 * programming tableau.
 *
 */

#define DP_MAX3    0x01
#define DP_TRACE   0x02
#define DP_MAX3_M0 0x04
#define DP_MAX3_M1 0x08
#define DP_MAX3_M2 0x10
#define DP_MAX2    0x20
#define DP_MAX2_M1 0x40
#define DP_MAX2_M2 0x80


#include "flint/flint.h"


typedef unsigned char dp_t;
typedef dp_t * dp_ptr;

typedef struct
{
    dp_t *data;
    slong nrows;
    slong ncols;
} dp_mat_struct;
typedef dp_mat_struct dp_mat_t[1];
typedef dp_mat_struct * dp_mat_ptr;


#ifdef __cplusplus
extern "C" {
#endif

static __inline__ int
dp_m0_is_interesting(dp_t x)
{
    return ((x & DP_MAX3) && (x & DP_MAX3_M0));
}

static __inline__ int
dp_m1_is_interesting(dp_t x)
{
    return (((x & DP_MAX3) && (x & DP_MAX3_M1)) ||
            ((x & DP_MAX2) && (x & DP_MAX2_M1)));
}

static __inline__ int
dp_m2_is_interesting(dp_t x)
{
    return (((x & DP_MAX3) && (x & DP_MAX3_M2)) ||
            ((x & DP_MAX2) && (x & DP_MAX2_M2)));
}



static __inline__ slong
dp_mat_nrows(const dp_mat_t mat)
{
    return mat->nrows;
}

static __inline__ slong
dp_mat_ncols(const dp_mat_t mat)
{
    return mat->ncols;
}

static __inline__ dp_ptr
dp_mat_entry(dp_mat_t mat, slong i, slong j)
{
    return mat->data + i * mat->ncols + j;
}

void dp_mat_init(dp_mat_t mat, slong nrows, slong ncols);
void dp_mat_clear(dp_mat_t mat);
void dp_mat_set(dp_mat_t mat, const dp_mat_t src);
void dp_mat_get_alignment(char *sa, char *sb, slong *plen,
        dp_mat_t mat, const slong *A, const slong *B);
void dp_mat_fprint(FILE *stream, const dp_mat_t mat);
void dp_mat_check_alignment(
        int *p_is_optimal, int *p_is_canonical,
        dp_mat_t mat, const slong *A, const slong *B, slong len);
void dp_mat_backward(dp_mat_t mat);


#ifdef __cplusplus
}
#endif

#endif
